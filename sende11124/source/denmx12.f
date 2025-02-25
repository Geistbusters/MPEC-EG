c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine denmx12(mdon,mdoo,ibasis,iptrbit,ncsf,nfcn,            7d29s22
     $     vint,nctf,isymmrci,irel,ism,irefo,multh,ixw1,                7d28s22
     $     ixw2,nh0av,nroot,ihsdiag,nff0,iff0,nff1,iff1,nsing,ndoub,    3d9s21
     $     maxbx,norb,nff22,nfdat,vdnono,nec,ncsf2,ih0av,               7d28s22
     $     nbasisp,ncomp,lwrite,sr2,srh,ioooo,shift,bc,ibc,ionex,jmats, 7d10s23
     $     kmats,idoubo,iff,idarot,nder,natom,ngaus,ibdat,ihmat,iorb,   7d17s23
     $     noc,ipair,isym,iapair,ibstor,isstor,iptoh,idorel,ascale,     7d17s23
     $     i3x,i4x,i4xb,isend1,isend2,isend3,isend4,ih0copy,iunc,nff2,  9d5s23
     $     vdinout,mdoub,hss,iden1e,idenhvv,idenj,idenk,idenput,idenorg,10d16s23
     $     ncd,nbasdwsc,tovr,ibufvcvds,idvmt,potdws,idervcv,ipropsym)   10d30s24
      implicit real*8 (a-h,o-z)
      external second                                                   10d17s23
c                                                                       3d3s21
c     compute one and two particle density matrices in ao basis,        9d29s22
c     as well as orbital terms for gradients                            9d29s22
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
      integer*8 ihsdiag(mdoo+1,*),i18,i28,i48                           7d17s23
      character*10 dlabel
      include "common.hf"                                               7d28s22
      include "common.print"                                            6d14s24
      dimension ibasis(3,*),iptrbit(2,mdoo+1),ncsf(*),vint(nctf,*),        3d3s21
     $     irel(*),ism(*),irefo(*),multh(8,8),nh0av(*),iden(8),noc(*),  7d18s23
     $     nff1(mdoo+1,*),iff1(*),nff0(*),iff0(*),vdnono(*),ncsf2(4,*), 3d11s21
     $     nbasisp(*),iorbno(8),id4o(idbk),jmden(idbk),iptoh(8,8,8),    7d17s23
     $     kmden(idbk),id1x(idbk),id3x(idbk),ih0av(*),ioooo(*),iorbf(8),7d10s23
     $    ionex(*),jmats(*),kmats(*),iff(*),idoubo(*),iffp(8,3),iorb(8),7d18s23
     $    idora(3,8),itt(8),iunc(2),vdinout(*),idorth(4,8),nfdat(5,4,*),9d11s23
     $    idcont(8),hss(*),i3x(*),nxs(8),nbasdwsc(*),nff22(mdoo+1,2,*), 10d17s23
     $     nl(4),ichoice(11),idatta(7,3),data(3),iamatu(8),xjt(512),    5d20s24
     $     iaddr(36,2,6),naddr(36,3),itt4v(8,8),isou(8),isoub(8),       8d23s24
     $     isoua1(8),i3xb(512),ionexb(512),jmatt(512),kmatt(512),       6d12s24
     $     ionexc(512),kmatd(64),i3x3(512),ionexbt(512),d4v(4),         7d11s24
     $     jmtsd(512),kmtsd(512),ioood(512),ionxd(512),i3xxd(512),      7d11s24
     $     i4xbd(512),idervcv(2,8),ipropsym(*),iapair(3,*)              10d30s24
      include "common.store"
      include "common.basis"                                            4d29s24
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      common/drsigncm/drsign                                            8d20s24
      common/timerocm/tovrx,telapo(15)                                   4d26s18
      data loopx/200/
      ldebug=iprtr(32).ne.0
      if(lwrite)write(6,*)('Hi, my name is denmx12'),nroot,ibcoff
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
      ifder=ibcoff                                                      7d18s23
      ibcoff=ifder+nder*6                                               4d25s24
      call enough('denmx12.fder',bc,ibc)                                7d18s23
      do iz=ifder,ibcoff-1                                              7d18s23
       bc(iz)=0d0                                                       7d18s23
      end do                                                            7d18s23
      ibcoffo=ibcoff                                                    3d3s21
      igoal=1962320
      do isb=1,nsymb                                                    3d3s21
       iden(isb)=ibcoff                                                 3d3s21
       ibcoff=iden(isb)+nh0av(isb)*nh0av(isb)*nroot                     3d3s21
       iorbf(isb)=ibcoff                                                9d29s22
       ibcoff=iorbf(isb)+nh0av(isb)*nh0av(isb)*nroot                    9d29s22
      end do                                                            3d3s21
      do is=1,nsdlk                                                     7d28s22
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
      if(lwrite)write(6,*)('computing densities')                       6d14s24
      call second(time1)                                                10d17s23
      call hccsfd12(vint,nctf,ibasis,ncsf,iptrbit,nfcn,iden,id4o,nroot, 7d29s22
     $     mdon,nec,isymmrci,ixw1,ixw2,mdoo,nh0av,iorbf,bc,ibc,igoal)         11d10s22
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
      if(ndoub.gt.0)then                                                9d5s23
       if(iunc(2).eq.0)then                                             9d5s23
        nfdatd=ibcoff                                                   9d18s23
        ibcoff=nfdatd+4*(mdoo+1)                                        9d18s23
        call enough('denmx12.nfdatd',bc,ibc)                            9d18s23
        ipredorth=ibcoff                                                10d17s23
        call precder(idorth,idcont,nsymb,nfdat,nroot,mdon,mdoo,         9d11s23
     $       ncsf2,ibc(nff2),ibc(nfdatd),bc,ibc,nxs)                    10d16s23
        npredorth=ibcoff-ipredorth                                      10d17s23
        call second(time1)                                              10d17s23
        call hcdi12(ibc(nff2),nsymb,mdon,mdoo,multh,isymmrci,nvirt,     9d6s23
     $       ncsf,nec,ncsf2,nfdat,irel,ism,kmats,irefo,vdinout,ixw1,    9d6s23
     $       ixw2,norb,nff0,iff0,vint,nctf,nroot,.false.,srh,bc,ibc,    9d6s23
     $       kmden,vdnono,idorth,idcont,nsdlkk,isblkk,ibufvcvds,idvmt,  4d29s24
     $       nder,nff22)                                                4d29s24
        call second(time2)                                              10d17s23
        telap=time2-time1-tovr                                          10d17s23
        if(lwrite)write(6,*)('time for hcdi12: '),telap
c
        call dws_synca
        call ddi_iaccword0                                               8d11s22
        if(ldebug)then                                                  6d27s24
         ihds=ibcoff
         ibcoff=ihds+nroot*mdoub                                           11d25s20
         call enough('denmx12. 23',bc,ibc)
         do i=ihds,ibcoff-1                                                4d17s20
          bc(i)=0d0                                                        4d17s20
         end do                                                            4d17s20
         call hcdi(ibc(nff2),nsymb,mdon,mdoo,multh,isymmrci,nvirt,
     $      ncsf,nec,ncsf2,nfdat,irel,ism,kmats,irefo,idum,             9d5s23
     $      bc(ihds),ixw1,ixw2,norb,nff0,iff0,vint,                     9d5s23
     $      nctf,nroot,.false.,srh,bc,ibc)                              9d5s23
         jhds=ihds-1                                                    6d27s24
         dotdi=0d0
         do i=1,nroot*ndoub                                             6d27s24
          dotdi=dotdi+vdinout(i)*bc(jhds+i)                             6d27s24
         end do                                                         6d27s24
         write(6,*)('d*H*i '),dotdi,dotdi*2d0                                     6d27s24
         do i=ihds,ibcoff-1                                                4d17s20
          bc(i)=0d0                                                        4d17s20
         end do                                                            4d17s20
         nwiacc=0                                                         8d11s22
         write(6,*)('calling hcds '),ihds,maxbx,mdoub,ndoub
         call hcds(ihsdiag,nff1,iff1,nff22,nfdat,bc(ihds),vdnono,nsymb,  9d11s23
     $       mdon,mdoo,nec,multh,isymmrci,nvirt,ncsf,ncsf2,irel,ism,    9d11s23
     $       irefo,idum,ixw1,ixw2,norb,nroot,ih0av,nh0av,idum,i3x,idum, 9d11s23
     $       dum,mdoub,ndoub,.false.,maxbx,dum,dum,ionex,sr2,           9d12s23
     $       nwiacc,bc,ibc)                                             9d11s23
         call ddi_iaccword2(nwiacc)                                      9d12s23
         call dws_synca                                                 6d27s24
         call dws_gsumf(bc(ihds),nroot*mdoub)                           6d27s24
         dotds=0d0                                                      6d27s24
         do i=1,nroot*mdoub                                             6d27s24
          dotds=dotds+bc(jhds+i)*vdnono(i)                              6d27s24
         end do                                                         6d27s24
         write(6,*)('d*H*s '),dotds,2d0*dotds,2d0*(dotds+dotdi)
         do i=ihds,ibcoff-1                                                4d17s20
          bc(i)=0d0                                                        4d17s20
         end do                                                            4d17s20
         call hcddjko(nff22,nfdat,bc(ihds),vdnono,nsymb,mdon,mdoo,      6d27s24
     $      nec,multh,isymmrci,nvirt,ncsf,ncsf2,irel,ism,               9d28s23
     $      irefo,idum,ixw1,ixw2,norb,nroot,ih0av,nh0av,ioooo,          9d28s23
     $      jmats,kmats,ndoub,mdoub,0d0,dum,dum,sr2,srh,bc,ibc)         6d27s24
         call dws_synca                                                 6d27s24
         call dws_gsumf(bc(ihds),nroot*mdoub)                           6d27s24
         dotdd=0d0                                                      6d27s24
         do i=1,nroot*mdoub                                             6d27s24
          term=bc(jhds+i)*vdnono(i)                                     6d28s24
          dotdd=dotdd+term
         end do                                                         6d27s24
         write(6,*)('d*H*d'),dotdd,2d0*dotdd                                      6d27s24
         dottot=2d0*(dotdd+dotdi+dotds)                                 6d27s24
         write(6,*)('total dot '),dottot                                6d27s24
         ibcoff=ihds                                                    6d27s24
        end if
        igd=ibcoff                                                      9d13s23
        ibcoff=igd+mdoub*nroot                                          9d13s23
        call enough('denmx12.gd',bc,ibc)                                9d13s23
        do iz=igd,ibcoff-1                                              9d13s23
         bc(iz)=0d0                                                     9d13s23
        end do                                                          9d13s23
        call second(time1)                                              10d17s23
        igqqq=1
        call hcds12(ihsdiag,nff1,iff1,nff22,nfdat,bc(igd),vdnono,nsymb, 9d13s23
     $       mdon,mdoo,nec,multh,isymmrci,nvirt,ncsf,ncsf2,irel,ism,    9d11s23
     $       irefo,ixw1,ixw2,norb,nroot,ih0av,nh0av,i3x,mdoub,ndoub,    9d12s23
     $       .false.,maxbx,ionex,sr2,iden,id1x,id3x,idcont,idorth,      9d12s23
     $       vdinout,bc,ibc,nsdlk1,isblk1,ibc(nfdatd),ibufvcvds,idvmt,  4d29s24
     $       nder,igqqq)                                                      4d29s24
        call second(time2)
        telap=time2-time1-tovr
        if(lwrite)write(6,*)('time for hcds12'),telap
        call second(time1)
        call hcddjkd12(nff22,nfdat,bc(igd),vdnono,nsymb,mdon,mdoo,      9d28s23
     $      nec,multh,isymmrci,nvirt,ncsf,ncsf2,irel,ism,               9d28s23
     $      irefo,idum,ixw1,ixw2,norb,nroot,ih0av,nh0av,ioooo,          9d28s23
     $      jmats,kmats,ndoub,mdoub,dum,dum,dum,sr2,srh,bc,ibc,iden,    9d28s23
     $      id4o,jmden,kmden,idcont,idorth,vdinout,nsdlk,isblk,nsdlkk,  10d3s23
     $      isblkk)                                                     10d3s23
        call dws_synca                                                  10d17s23
        call dws_gsumf(bc(ipredorth),npredorth)                         10d17s23
        call second(time2)
        telap=time2-time1-tovr                                          10d17s23
        if(lwrite)write(6,*)('time for hcddjkd12 '),telap
        idoit=0
        jbufvcvds=ibufvcvds                                             10d17s23
        if(ldebug)then                                                  6d27s24
         dtest1=0d0                                                     6d27s24
         dtest2=0d0                                                     6d27s24
        end if                                                          6d27s24
        do isb=1,nsymb                                                  10d17s23
         iad1=idvmt+2*(isb-1)                                           10d17s23
         nwds=ibc(iad1+1)
         jbufv=1                                                        10d17s23
         call ilimts(1,nwds,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)   10d17s23
         nhere=ih+1-il                                                  10d17s23
         jbufvcvds=idervcv(1,isb)                                       8d7s24
         kbufvcvds=jbufvcvds
         lbufvcvds=jbufvcvds
         kbufv=jbufv                                                    10d17s23
         lbufv=jbufv                                                    10d17s23
         do l=1,4
          if(ldebug)then                                                6d27s24
           if(nfdat(2,l,isb).gt.0)then
            do i=0,nfdat(2,l,isb)*nfdat(3,l,isb)-1                      6d27s24
             iaddz=idorth(l,isb)+i
             dtest1=dtest1+bc(iaddz)*bc(nfdat(4,l,isb)+i)               6d27s24
            end do                                                      6d27s24
           end if
          end if                                                        6d27s24
          mbufvcvds=lbufvcvds                                           10d18s23
          mbufv=kbufv                                                   10d18s23
          do ider=1,nder                                                10d17s23
           ibot=mbufv                                                   10d18s23
           itop=mbufv+nfdat(2,l,isb)*nfdat(3,l,isb)-1                   10d18s23
           iad=ifder+ider-1+nder                                        4d25s24
           do i=max(ibot,il),min(itop,ih)                               10d17s23
            i1=i-il+jbufvcvds                                           10d17s23
            i2=i-mbufv+idorth(l,isb)                                    10d17s23
            term=bc(i1)*bc(i2)                                          10d17s23
            bc(iad)=bc(iad)+term                                        10d17s23
           end do                                                       10d17s23
           mbufvcvds=mbufvcvds+nfdat(2,l,isb)*nfdat(3,l,isb)
           mbufv=mbufv+nfdat(2,l,isb)*nfdat(3,l,isb)                    10d17s23
          end do                                                        10d17s23
          lbufvcvds=lbufvcvds+nfdat(2,l,isb)*nfdat(2,l,isb)*nder        10d17s23
          lbufv=lbufv+nfdat(2,l,isb)*nfdat(2,l,isb)*nder                10d17s23
          kbufv=kbufv+nfdat(2,l,isb)*nfdat(2,l,isb)*nder                10d17s23
         end do
         lvcv=nfdat(5,1,isb)                                            10d17s23
         lvcvd=idcont(isb)                                              10d17s23
         kbufv=lbufv                                                    10d17s23
         if(ldebug)then                                                 6d27s24
          lvcv=nfdat(5,1,isb)                                            10d17s23
          lvcvd=idcont(isb)                                              10d17s23
          ssum=0d0
          ifoff=0
          do nclo2=mdon,mdoo                                             10d17s23
           nclo2p=nclo2+1                                                10d17s23
           iarg2=nclo2p-mdon                                             10d17s23
           do if=0,nff22(nclo2p,1,isb)-1                                 10d17s23
            nspace=ibc(lvcv+1)
            do l=1,4
             nl(l)=ibc(lvcv+1+l)
             if(nl(l).gt.0)then                                            9d7s23
              iad1=lvcv+ibc(lvcv+5+l)                                      9d8s23
              iad2=iad1+nl(l)                                              9d7s23
              do i=0,nl(l)-1
               ii=ibc(iad1+i)
               do j=0,ncsf2(l,iarg2)-1
                ji=j+ncsf2(l,iarg2)*i
                term=bc(iad2+ji)*bc(lvcvd+ji)                              6d27s24
                dtest2=dtest2+term                                       6d27s24
               end do
              end do                                                    6d27s24
              lvcvd=lvcvd+ncsf2(l,iarg2)*nl(l)*nroot                     10d17s23
             end if                                                      10d17s23
            end do
            ifoff=ifoff+1
            lvcv=lvcv+nspace                                             10d17s23
           end do                                                        10d17s23
          end do                                                         10d17s23
         end if                                                         6d27s24
         do ider=1,nder                                                 10d17s23
          lvcv=nfdat(5,1,isb)                                            10d17s23
          lvcvd=idcont(isb)                                              10d17s23
          ssum=0d0
          ifoff=0
          do nclo2=mdon,mdoo                                             10d17s23
           nclo2p=nclo2+1                                                10d17s23
           iarg2=nclo2p-mdon                                             10d17s23
           do if=0,nff22(nclo2p,1,isb)-1                                 10d17s23
            nspace=ibc(lvcv+1)
            do l=1,4
             nl(l)=ibc(lvcv+1+l)
             if(nl(l).gt.0)then                                            9d7s23
              iad1=lvcv+ibc(lvcv+5+l)                                      9d8s23
              iad2=iad1+nl(l)                                              9d7s23
              ibot=kbufv                                                10d17s23
              itop=kbufv+ncsf2(l,iarg2)*nl(l)-1                         10d17s23
              iad=ifder+ider-1                                          10d17s23
              do i=max(ibot,il),min(itop,ih)                            10d17s23
               i1=i-il+jbufvcvds                                        10d17s23
               i2=i-kbufv+lvcvd                                         10d17s23
               term=bc(i1)*bc(i2)                                       10d17s23
               bc(iad)=bc(iad)+term
               ssum=ssum+term
              end do                                                    10d17s23
              kbufv=kbufv+ncsf2(l,iarg2)*nl(l)                          10d17s23
              lvcvd=lvcvd+ncsf2(l,iarg2)*nl(l)*nroot                     10d17s23
             end if                                                      10d17s23
            end do
            ifoff=ifoff+1
            lvcv=lvcv+nspace                                             10d17s23
           end do                                                        10d17s23
          end do                                                         10d17s23
         end do                                                         10d17s23
         jbufvcvds=jbufvcvds+nhere                                      10d17s23
        end do                                                          10d17s23
        if(ldebug)then                                                  6d27s24
         write(6,*)('dtest1 = '),dtest1
         write(6,*)('dtest2 = '),dtest2
        end if                                                          6d27s24
c
c     4v part
c
        ihdd=ibcoff
        ivdd=ihdd+ndoub*nroot
        ibcoff=ivdd+ndoub*nroot
        call enough('denmx12.ihdd',bc,ibc)
        do iz=ihdd,ibcoff-1
         bc(iz)=0d0
        end do
        jvdd=ivdd-1                                                     10d16s23
        do iz=1,ndoub*nroot                                             10d16s23
         bc(jvdd+iz)=vdinout(iz)                                        10d16s23
        end do                                                          10d16s23
        call reordergv(bc(ivdd),nroot,nfdat,nvirt,nsymb,multh,          10d16s23
     $       isymmrci,ndoub,1,bc,ibc)                                   10d16s23
        call hcdd4v(bc(ivdd),bc(ihdd),ncd,multh,iorb,nbasdwsc,natom,     10d16s23
     $      ngaus,ibdat,ipair,isym,iapair,ibstor,isstor,idorel,         11d2s20
     $      ascale,ndoub,nroot,nbasisp,ndoub,sr2,srh,bc,ibc)            10d16s23
        jhdd=ihdd-1
        if(ldebug)then                                                  6d20s24
         do iz=1,ndoub
          dot4v=dot4v+bc(jvdd+iz)*bc(jhdd+iz)                            10d16s23
         end do                                                          10d16s23
         write(6,*)('dot4v: '),dot4v                                     10d16s23
        end if                                                          6d20s24
       else                                                             9d5s23
        write(6,*)('we have no code for unmrci ders yet')               9d5s23
        write(6,*)('Sorry!')                                            9d5s23
        ibcoff=ibcoffo                                                  9d5s23
        return                                                          9d5s23
       end if                                                           9d5s23
      end if                                                            9d5s23
      if(ldebug)then                                                    6d14s24
       dot1e=0d0                                                         7d11s23
       do isb=1,nsymb                                                    7d11s23
        if(nh0av(isb).gt.0)then                                          7d11s23
         trace=0d0                                                       7d11s23
         do i1=0,nh0av(isb)-1
          do i2=0,i1-1
           i12=i1+nh0av(isb)*i2
           i21=i2+nh0av(isb)*i1
           dsum=bc(iden(isb)+i12)+bc(iden(isb)+i21)
           term=dsum*bc(ih0av(isb)+i12)
           trace=trace+term
          end do
          i11=i1*(nh0av(isb)+1)
          term=bc(iden(isb)+i11)*bc(ih0av(isb)+i11)
          trace=trace+term
         end do
         dot1e=dot1e+trace                                               7d11s23
         call dws_gsumf(dot1e,1)
         write(6,*)('trace for h0 of '),isb,(' is '),trace,              7d11s23
     $       (' dot1e so far '),dot1e                                   7d11s23
        end if                                                           7d11s23
       end do                                                            7d11s23
       dot2e=0d0                                                         7d10s23
       do is=1,nsdlk                                                     7d10s23
        if(isblk(1,is).eq.isblk(2,is))then                               7d28s22
         nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              7d28s22
         ncol=(irefo(isblk(3,is))*(irefo(isblk(3,is))+1))/2              7d28s22
        else                                                             7d28s22
         nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      7d28s22
         ncol=irefo(isblk(3,is))*irefo(isblk(4,is))                       7d28s22
        end if                                                           7d28s22
        if(min(nrow,ncol).gt.0)then                                      7d10s23
         trace=0d0                                                        7d10s23
         do ii=0,nrow*ncol-1                                              7d10s23
          trace=trace+bc(id4o(is)+ii)*bc(ioooo(is)+ii)                    7d10s23
         end do                                                           7d10s23
         dot2e=dot2e+trace                                                7d10s23
        end if                                                           7d10s23
       end do                                                            7d10s23
       if(max(ndoub,nsing).gt.0)then                                     10d16s23
        tt=0d0
        tt3=0d0                                                          10d16s23
        do is=1,nsdlk1                                                    7d28s22
         if(isblk1(1,is).eq.isblk1(2,is))then                             7d28s22
          nrow1=(irefo(isblk1(1,is))*(irefo(isblk1(1,is))+1))/2           7d28s22
          nrow3=(nvirt(isblk1(1,is))*(nvirt(isblk1(1,is))+1))/2          10d16s23
         else                                                             7d28s22
          nrow1=irefo(isblk1(1,is))*irefo(isblk1(2,is))                   7d28s22
          nrow3=nvirt(isblk1(1,is))*nvirt(isblk1(2,is))                  10d16s23
         end if                                                           7d28s22
         ncol=irefo(isblk1(3,is))*nvirt(isblk1(4,is))
         if(min(nrow1,ncol).gt.0)then                                     7d10s23
          trace=0d0                                                       7d10s23
          do ii=0,nrow1*ncol-1                                            7d10s23
           trace=trace+bc(id1x(is)+ii)*bc(ionex(is)+ii)                   7d10s23
          end do                                                          7d10s23
          dot2e=dot2e+trace                                             10d16s23
          tt=tt+trace                                                   10d16s23
         end if                                                          10d16s23
         call ilimts(irefo(isblk1(3,is)),nvirt(isblk1(4,is)),mynprocg,   10d16s23
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           10d16s23
         ncol=ih+1-il                                                    10d16s23
         if(min(nrow3,ncol).gt.0)then                                     7d10s23
          trace=0d0                                                       7d10s23
          do ii=0,nrow3*ncol-1                                            7d10s23
           trace=trace+bc(id3x(is)+ii)*bc(i3x(is)+ii)                    10d16s23
          end do                                                          7d10s23
          dot2e=dot2e+trace                                              10d16s23
          tt3=tt3+trace                                                  10d16s23
         end if                                                           7d10s23
        end do                                                            7d28s22
        ttt=tt
        call dws_gsumf(ttt,1)
        write(6,*)('total for 1x '),tt,ttt
        tt3t=tt3
        call dws_gsumf(tt3t,1)
        call dws_gsumf(tt3,1)                                              10d16s23
        write(6,*)('total for 3x '),tt3,tt3t
       end if
       if(max(nsing,ndoub).gt.0)then                                     7d11s23
        tt=0d0
        do is=1,nsdlk                                                     7d28s22
         if(isblk(1,is).eq.isblk(2,is))then                               7d28s22
          nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              7d28s22
         else                                                             7d28s22
          nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      7d28s22
         end if                                                           7d28s22
         call ilimts(nvirt(isblk(3,is)),nvirt(isblk(4,is)),mynprocg,      7d28s22
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            7d28s22
         ncol=ih+1-il                                                     7d28s22
         if(min(nrow,ncol).gt.0)then                                      7d10s23
          trace=0d0                                                       7d10s23
          do ii=0,nrow*ncol-1                                             7d10s23
           trace=trace+bc(jmden(is)+ii)*bc(jmats(is)+ii)                  7d10s23
          end do                                                          7d11s23
          dot2e=dot2e+trace                                               7d11s23
          tt=tt+trace
         end if                                                           7d11s23
        end do                                                            7d28s22
        ttt=tt
        call dws_gsumf(ttt,1)
        write(6,*)('total for J '),tt,ttt
        tt=0d0
        do is=1,nsdlkk                                                    7d28s22
         nrow=irefo(isblkk(1,is))*irefo(isblkk(2,is))                     7d28s22
         call ilimts(nvirt(isblkk(3,is)),nvirt(isblkk(4,is)),mynprocg,    7d28s22
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            7d28s22
         ncol=ih+1-il                                                     7d28s22
         if(min(nrow,ncol).gt.0)then                                      7d11s23
          trace=0d0                                                       7d11s23
          do ii=0,nrow*ncol-1                                             7d11s23
           trace=trace+bc(kmden(is)+ii)*bc(kmats(is)+ii)                  7d11s23
          end do                                                          7d11s23
          dot2e=dot2e+trace                                               7d11s23
          tt=tt+trace
         end if                                                           7d11s23
        end do                                                            7d28s22
        ttt=tt                                                           10d16s23
        call dws_gsumf(ttt,1)
        write(6,*)('total for K '),tt,ttt
       end if
       bc(ibcoff)=dot1e                                                  7d11s23
       bc(ibcoff+1)=dot2e                                                7d11s23
       call dws_gsumf(bc(ibcoff),2)                                      7d11s23
       dot1e=bc(ibcoff)                                                  7d11s23
       dot2e=bc(ibcoff+1)                                                7d11s23
       write(6,*)('full dot1e: '),dot1e                                  7d11s23
       write(6,*)('full dot2e: '),dot2e                                  7d11s23
       write(6,*)('shift: '),shift                                       7d11s23
       efromt=shift+dot1e+dot2e                                          7d11s23
       write(6,*)('energy from traces: '),efromt                         7d11s23
       write(6,*)('with 4v: '),efromt+dot4v                              10d16s23
      end if                                                            6d14s24
      call dws_gsumf(bc(ibcoffo),ndstot)                                 5d2s23
      ee=0d0                                                            7d17s23
      do isb=1,nsymb                                                    7d18s23
       do iz=0,nbasdws(isb)*nbasdws(isb)-1                              7d18s23
        bc(iff(isb)+iz)=0d0                                             7d18s23
       end do                                                           7d18s23
      end do                                                            7d18s23
      jdarot=idarot                                                     7d18s23
      do ider=0,nder-1                                                  7d18s23
       ipass=ibc(jdarot+2)
       jdarot=jdarot+4                                                  7d18s23
       do isb=1,nsymb                                                   7d18s23
        do ii=0,nbasdws(isb)*nbasdws(isb)-1                             7d18s23
         bc(iff(isb)+ii)=bc(iff(isb)+ii)+abs(bc(jdarot+ii))             7d18s23
        end do                                                          7d18s23
        jdarot=jdarot+nbasdws(isb)*nbasdws(isb)                         7d18s23
        if(ipass.ge.0)then                                              8d15s23
         jdarot=jdarot+nbasdws(isb)*nbasdws(isb)                        8d15s23
        end if                                                          8d15s23
       end do                                                           7d18s23
      end do                                                            7d18s23
      do isb=1,nsymb
       if(nh0av(isb).gt.0)then
        do i=0,nh0av(isb)-1                                             5d5s23
         do j=0,i-1                                                     5d5s23
          ji=iden(isb)+j+nh0av(isb)*i                                   5d5s23
          ij=iden(isb)+i+nh0av(isb)*j                                   5d5s23
          avg=0.5d0*(bc(ji)+bc(ij))                                     5d5s23
          bc(ji)=avg                                                    5d5s23
          bc(ij)=avg                                                    5d5s23
         end do                                                         5d5s23
        end do                                                          5d5s23
        if(ldebug)then                                                  6d20s24
         write(6,*)('for symmetry block '),isb,iden(isb)
         call prntm2(bc(iden(isb)),nh0av(isb),nh0av(isb),nh0av(isb))     5d2s23
        end if                                                          6d20s24
        iold=-1                                                         7d18s23
        nold=0                                                          7d18s23
        idora(1,isb)=0                                                  7d18s23
        idora(2,isb)=0                                                  7d18s23
        idora(3,isb)=0                                                  7d17s24
        do i=0,nbasdws(isb)-1                                             7d18s23
         iad=iff(isb)+nbasdws(isb)*i                                    7d18s23
         do j=i+1,nbasdws(isb)-1                                          7d18s23
          if(abs(bc(iad+j)).gt.1d-10)then                               7d18s23
           if(j.ne.iold)then                                            7d18s23
            nold=nold+1                                                 7d18s23
            if(nold.le.2)then                                           7d18s23
             idora(nold,isb)=j                                          7d18s23
            end if                                                      7d18s23
            iold=j                                                      7d18s23
           end if                                                       7d18s23
           go to 777                                                    7d18s23
          end if                                                        7d18s23
         end do                                                         7d18s23
  777    continue                                                       7d18s23
        end do                                                          7d18s23
        if(ldebug)then                                                  6d20s24
         trace=0d0                                                       7d17s23
         do i=0,nh0av(isb)*nh0av(isb)-1                                  7d17s23
          trace=trace+bc(iden(isb)+i)*bc(ih0av(isb)+i)                   7d17s23
         end do                                                          7d17s23
         write(6,*)('trace '),trace
         ee=ee+trace                                                     7d17s23
        end if                                                          6d20s24
       end if
      end do
      do isb=1,nsymb                                                    7d18s23
       if(idora(2,isb).eq.0.and.idora(1,isb).ne.0)then                  7d19s23
        idora(2,isb)=idora(1,isb)                                       7d19s23
        idora(1,isb)=0                                                  7d19s23
       end if                                                           7d19s23
       if(idora(2,isb).gt.0)idora(2,isb)=idora(2,isb)-idora(1,isb)      7d18s23
       idora(3,isb)=irefo(isb)-idora(2,isb)                             7d18s23
      end do                                                            7d18s23
      if(ldebug)write(6,*)('2-e densities at end of denmx12 '),bc(igoal)
      call postden(nsdlk,isblk,nsdlk1,isblk1,nsdlkk,isblkk,irefo,nvirt, 8d6s24
     $     id4o,ioooo,id1x,ionex,jmden,jmats,kmden,kmats,id3x,i3x,      8d6s24
     $     ldebug,nsing,ndoub,dot4v,shift,ee,igoal,bc,ibc)              8d6s24
      do isb=1,nsymb                                                    9d1s23
       itt(isb)=ibcoff                                                  9d1s23
       ibcoff=itt(isb)+nbasdws(isb)*nbasdws(isb)                        9d1s23
      end do                                                            9d1s23
      call enough('denmx12.itt',bc,ibc)                                 9d1s23
      do iz=itt(1),ibcoff-1                                             9d1s23
       bc(iz)=0d0                                                       9d1s23
      end do                                                            9d1s23
      if(ldebug)write(6,*)('4x '),bc(igoal)                             6d20s24
      call second(time1)                                                6d25s24
      ibcb4=ibcoff                                                      7d7s21
      do i1=1,nsymb                                                     7d22s14
       do i2=1,nsymb                                                    7d22s14
        do i3=1,nsymb                                                   7d22s14
         iptoh(i3,i2,i1)=0                                              7d22s14
        end do                                                          7d22s14
       end do                                                           7d22s14
      end do                                                            7d22s14
      do i=1,nsdlkh                                                     5d12s10
       iptoh(isblkh(1,i),isblkh(2,i),isblkh(3,i))=i                     5d12s10
      end do                                                            5d12s10
      idwsdeb=0
      if(ldebug)then                                                    6d20s24
       ifull4x=3                                                         8d18s23
      else                                                              6d20s24
       ifull4x=2                                                        6d20s24
      end if                                                            6d20s24
      if(min(nsing,ndoub).gt.0.and.ichoice(9).ne.0.and                  6d14s24
     $     .ichoice(11).ne.0)then                                       6d14s24
       isnd=ibcoff                                                      11d30s23
       ircv=isnd+mynprocg                                               11d30s23
       ipt=ircv+mynprocg                                                11d30s23
       ipf=ipt+mynprocg                                                 11d30s23
       nrowp=ipf+mynprocg                                               12d1s23
       ncolp=nrowp+mynprocg                                             12d1s23
       ibcoff=ncolp+mynprocg                                            12d1s23
       call enough('denmx12.snd',bc,ibc)                                11d30s23
       ipair=0                                                          4d17s24
       itthalf=ibcoff                                                   2d26s24
       jtthalf=itthalf                                                  2d26s24
       do isb=1,nsymb                                                   2d26s24
        jtthalf=jtthalf+irefo(isb)*nvirt(isb)                           4d17s24
       end do                                                           2d26s24
       ibcoff=jtthalf                                                   2d26s24
       call enough('denmx12.jtthalf',bc,ibc)                            2d26s24
       do iz=itthalf,ibcoff-1                                           2d26s24
        bc(iz)=0d0                                                      2d26s24
       end do                                                           2d26s24
       ibc(ibcoff)=1                                                    4d17s24
       ibc(ibcoff+1)=itthalf                                            2d26s24
       ibc(ibcoff+2)=id3x(1)                                            2d26s24
       ibufs=ibcoff                                                     2d26s24
       ibcoff=ibcoff+3
       do isb=1,nsymb                                                   4d4s24
        ibc(ibcoff)=irefo(isb)                                          4d4s24
        ibcoff=ibcoff+1                                                 4d4s24
       end do                                                           4d4s24
       ipair=1                                                          6d14s24
      else                                                              11d29s23
       ipair=0                                                          11d29s23
       ibufs=1                                                          11d29s23
      end if                                                            11d29s23
      if(ndoub.gt.1.and.ichoice(10).ne.0)then                           6d14s24
       if(idorel.eq.0)then                                              5d29s24
        idorelu=0                                                       5d29s24
       else                                                             5d29s24
        idorelu=3                                                       5d29s24
       end if                                                           5d29s24
       ivtmp=ibcoff                                                     6d5s24
       ibcoff=ivtmp+ndoub*nroot                                         6d5s24
       call enough('denmx12.vtmp',bc,ibc)                               6d5s24
       jvtmp=ivtmp-1                                                    6d5s24
       do i=1,ndoub*nroot                                               6d5s24
        bc(jvtmp+i)=vdinout(i)                                          6d5s24
       end do                                                           6d5s24
       call reordergv(bc(ivtmp),nroot,nfdat,nvirt,nsymb,multh,isymmrci, 6d5s24
     $      ndoub,1,bc,ibc)                                             6d5s24
       call vdmo2so(bc(ivtmp),iorb,multh,ncd,nbasdwsc,iaddr,naddr,nroot,6d5s24
     $     ndoub,nbasisp,idorelu,srh,1,bc,ibc)                          6d14s24
       do isbv2=1,nsymb
        do isbv1=1,isbv2
         itv12=((isbv2*(isbv2-1))/2)+isbv1
         nv1=nbasisp(isbv1)
         nv2=nbasisp(isbv2)
         nvvs=nv1*nv2
         nvvt=nvvs
         isub=1
         nvv=nvvs
         do ist=1,2
          if(naddr(itv12,ist).gt.0)then
           do ierri=1,ncomp*ncomp                                       6d11s24
            i1=iaddr(itv12,ist,ierri)
            i2=nvv
            i3=naddr(itv12,ist)
            if(ist.eq.1)then
            else
            end if
            itmpx=ibcoff                                                 6d11s24
            ibcoff=itmpx+nvv*naddr(itv12,ist)                            6d11s24
            call enough('denmx12.tmpxb',bc,ibc)                          6d11s24
            do ii=0,naddr(itv12,ist)-1
             do ivv=0,nvv-1
              iad=iaddr(itv12,ist,ierri)+ivv+nvv*ii
              jtmpx=itmpx+ii+naddr(itv12,ist)*ivv                        6d11s24
              bc(jtmpx)=bc(iad)                                          6d11s24
             end do
            end do
            do ii=0,naddr(itv12,ist)*nvv-1                              6d11s24
             bc(iaddr(itv12,ist,ierri)+ii)=bc(itmpx+ii)                 6d11s24
            end do                                                      6d11s24
            ibcoff=itmpx                                                6d11s24
           end do
          end if
          isub=-1
          nvv=nvvt
         end do
        end do
       end do
       jvtmp=ivtmp                                                      6d5s24
       xnan=-2d0
       do iz=0,ndoub-1
        bc(ivtmp+iz)=xnan
       end do
       do i=1,15                                                        8d23s24
        telapo(i)=0d0                                                   8d23s24
       end do                                                           8d23s24
       do isbv2=1,nsymb                                                 6d5s24
        do isbv1=1,isbv2                                                6d5s24
         itv12=((isbv2*(isbv2-1))/2)+isbv1                              6d5s24
         if(isbv1.eq.isbv2)then                                         6d5s24
          nvvs=(nvirt(isbv1)*(nvirt(isbv1)+1))/2                        6d5s24
          nvvt=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        6d5s24
          isw=0                                                         6d5s24
         else                                                           6d5s24
          nvvs=nvirt(isbv1)*nvirt(isbv2)                                6d5s24
          nvvt=nvvs                                                     6d5s24
          isw=1                                                         6d5s24
         end if                                                         6d5s24
         do ist=1,2                                                     6d5s24
          iaddr(itv12,ist,5)=jvtmp                                      6d5s24
          ioff=1                                                        6d5s24
          do jsb=1,nsymb                                                6d5s24
           jsbv12=multh(jsb,isymmrci)                                   6d5s24
           do jsbv1=1,nsymb                                             6d5s24
            jsbv2=multh(jsbv1,jsbv12)                                   6d5s24
            if(jsbv2.ge.jsbv1)then                                      6d5s24
             if(jsbv1.eq.jsbv2)then                                     6d5s24
              mvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                     6d5s24
              ioffvv=ioff                                               6d5s24
              ioff=ioff+nvirt(jsbv1)*nfdat(3,1,jsb)                     6d5s24
             else                                                       6d5s24
              mvv=nvirt(jsbv1)*nvirt(jsbv2)                             6d5s24
             end if                                                     6d5s24
             if(isbv1.eq.jsbv1.and.isbv2.eq.jsbv2)then                  6d5s24
              if(ist.eq.1.and.isw.eq.0)then                             6d5s24
               do i=0,nfdat(3,1,jsb)-1                                  6d5s24
                do jv=0,nvirt(jsbv1)-1                                  6d5s24
                 jtri=((jv*(jv+1))/2)+jv                                6d5s24
                 iad=jvtmp+jtri+nvvs*i                                  6d5s24
                 jad=ioffvv+jv+nvirt(isbv1)*i                           6d5s24
                 bc(iad)=vdinout(jad)*srh                               6d5s24
                end do                                                  6d5s24
               end do                                                   6d5s24
              end if                                                    6d5s24
              if(ist.eq.1)then                                          6d5s24
               joff=ioff                                                6d5s24
               do i=0,nfdat(3,1,jsb)-1                                  6d5s24
                jvv=0
                do jv2=0,nvirt(isbv2)-1                                 6d5s24
                 jv1top=jv2+isw*(nvirt(jsbv1)-jv2)-1                    6d5s24
                 do jv1=0,jv1top                                        6d5s24
                  jrec=jv1+nvirt(jsbv1)*jv2                             6d5s24
                  jtri=((jv2*(jv2+1))/2)+jv1                            6d5s24
                  jtri=jtri+isw*(jrec-jtri)                             6d5s24
                  iad=jvtmp+jtri+nvvs*i                                 6d5s24
                  bc(iad)=vdinout(joff)                                 6d5s24
                  jvv=jvv+1
                  joff=joff+1                                           6d5s24
                 end do                                                 6d5s24
                end do                                                  6d5s24
               end do                                                   6d5s24
               jvtmp=jvtmp+nvvs*nfdat(3,1,jsb)                          6d5s24
              else                                                      6d5s24
               joff=ioff+mvv*nfdat(3,1,jsb)                             6d5s24
               ii=0                                                     6d5s24
               do l=2,4                                                 6d5s24
                do i=0,nfdat(3,l,jsb)-1                                 6d5s24
                 jvv=0
                 do jv2=0,nvirt(isbv2)-1                                6d5s24
                  jv1top=jv2+isw*(nvirt(jsbv1)-jv2)-1                   6d5s24
                  do jv1=0,jv1top                                       6d5s24
                   jrec=jv1+nvirt(jsbv1)*jv2                            6d5s24
                   jtri=((jv2*(jv2-1))/2)+jv1                           6d5s24
                   jtri=jtri+isw*(jrec-jtri)                            6d5s24
                   iad=jvtmp+jtri+nvvt*ii                               6d5s24
                   bc(iad)=vdinout(joff)                                6d5s24
                   joff=joff+1                                          6d5s24
                   jvv=jvv+1
                  end do                                                6d5s24
                 end do                                                 6d5s24
                 ii=ii+1                                                6d5s24
                end do                                                  6d5s24
               end do                                                   6d5s24
               jvtmp=jvtmp+nvvt*ii                                      6d5s24
              end if                                                    6d5s24
             end if                                                     6d5s24
             do l=1,4                                                   6d5s24
              ioff=ioff+mvv*nfdat(3,l,jsb)                              6d5s24
             end do                                                     6d5s24
            end if                                                      6d5s24
           end do                                                       6d5s24
          end do                                                        6d5s24
         end do                                                         6d5s24
        end do                                                          6d5s24
       end do                                                           6d5s24
       do isbv2=1,nsymb                                                 6d5s24
        do isbv1=1,isbv2                                                6d5s24
         itv12=((isbv2*(isbv2-1))/2)+isbv1                              6d5s24
         if(isbv1.eq.isbv2)then                                         6d5s24
          nvvs=(nvirt(isbv1)*(nvirt(isbv1)+1))/2                        6d5s24
          nvvt=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        6d5s24
          isw=0                                                         6d5s24
         else                                                           6d5s24
          nvvs=nvirt(isbv1)*nvirt(isbv2)                                6d5s24
          nvvt=nvvs                                                     6d5s24
          isw=1                                                         6d5s24
         end if                                                         6d5s24
         nvv=nvvs                                                       6d5s24
         do ist=1,2                                                     6d5s24
          itmpx=ibcoff                                                  6d11s24
          ibcoff=itmpx+nvv*naddr(itv12,ist)                             6d11s24
          call enough('denmx12.tmpx',bc,ibc)                            6d11s24
          do i=0,naddr(itv12,ist)-1
           do ivv=0,nvv-1                                               6d5s24
            iad=iaddr(itv12,ist,5)+ivv+nvv*i                            6d5s24
            jtmpx=itmpx+i+naddr(itv12,ist)*ivv                          6d11s24
            bc(jtmpx)=bc(iad)                                           6d11s24
           end do                                                       6d5s24
          end do                                                        6d5s24
          iaddr(itv12,ist,6)=iaddr(itv12,ist,5)                         8d23s24
          iaddr(itv12,ist,5)=itmpx                                      8d23s24
          nvv=nvvt                                                      6d5s24
         end do
        end do                                                          6d5s24
       end do                                                           6d5s24
       ibc000=ibcoff                                                    8d27s24
       if(ldebug)then                                                   8d27s24
        do iperm=1,8                                                     8d23s24
         do isb=1,nsymb                                                   5d20s24
          itt4v(isb,iperm)=ibcoff                                               5d20s24
          ibcoff=ibcoff+nvirt(isb)*nbasisp(isb)*ncomp                     6d5s24
         end do
        end do                                                           5d20s24
        ittkeep1=1                                                      8d27s24
        ittkeep3=3                                                      8d27s24
       else                                                             8d27s24
        do isb=1,nsymb                                                  8d27s24
         do iperm=1,8                                                   8d27s24
          itt4v(isb,iperm)=-1                                           8d27s24
         end do                                                         8d27s24
         itt4v(isb,5)=ibcoff                                            8d27s24
         itt4v(isb,7)=itt4v(isb,5)+nvirt(isb)*nbasisp(isb)*ncomp        8d27s24
         ibcoff=itt4v(isb,7)+nvirt(isb)*nbasisp(isb)*ncomp              8d27s24
         ittkeep1=5                                                      8d27s24
         ittkeep3=7                                                      8d27s24
        end do                                                          8d27s24
       end if                                                           8d27s24
       call enough('denmx12.tt4v',bc,ibc)                               5d20s24
       do iz=ibc000,ibcoff-1                                            8d27s24
        bc(iz)=0d0                                                      5d20s24
       end do                                                           5d20s24
      else                                                              5d20s24
       itt4v(1,1)=0                                                       5d20s24
      end if                                                            5d20s24
      call paraeri(natom,ngaus,ibdat,idum,ihmat,iorb,noc,               4d24s20
     $             ipair,nhcolt,isym,iapair,ibstor,isstor,multh,iptoh,  11d29s12
     $     0,idwsdeb,idorel,ascale,nbasisp,ifull4x,ipair,ibufs,         4d5s24
     $     itt4v,iaddr,naddr,bc,ibc)                                    5d20s24
      if(ibufs.ne.1)then                                                4d22s24
       if(ldebug)write(6,*)('uu matrix ')                                         4d22s24
       jtthalf=itthalf                                                  4d22s24
       do isb=1,nsymb                                                   4d22s24
        if(min(irefo(isb),nvirt(isb)).gt.0)then                           4d22s24
         if(ldebug)then                                                 6d20s24
          write(6,*)('for symmetry '),isb,jtthalf                                4d22s24
          call prntm2(bc(jtthalf),nvirt(isb),irefo(isb),nvirt(isb))        4d22s24
          write(6,*)('copy to tt ')
         end if                                                         6d20s24
         jtth=jtthalf                                                   4d22s24
         do i=0,irefo(isb)-1                                            4d22s24
          ip=i+idoubo(isb)                                              4d22s24
          do j=0,nvirt(isb)-1                                           4d22s24
           jp=j+noc(isb)                                                4d22s24
           iad1=itt(isb)+jp+nbasdws(isb)*ip                             4d22s24
           iad2=itt(isb)+ip+nbasdws(isb)*jp                             4d22s24
           bc(iad2)=bc(iad2)+bc(jtth+j)                                 4d22s24
          end do                                                        4d22s24
          jtth=jtth+nvirt(isb)                                          4d22s24
         end do                                                         4d22s24
         jtthalf=jtthalf+nvirt(isb)*irefo(isb)                            4d22s24
        end if                                                          4d22s24
       end do                                                           4d22s24
      end if                                                            4d22s24
      if(ndoub.gt.1.and.ichoice(10).ne.0)then                           5d29s24
       if(ldebug)write(6,*)('in transform block ')                      6d20s24
       do isb=1,nsymb                                                   5d29s24
        if(nbasisp(isb).gt.0)then                                       5d29s24
         do iperm=1,4
          nb=nbasisp(isb)*ncomp                                          5d29s24
          if(ldebug)then                                                6d20s24
           write(6,*)('for perm = '),iperm
           write(6,*)
     $         ('tt4v in half primative basis for symmetry block '),
     $       isb,('address '),itt4v(isb,iperm)
           if(itt4v(isb,iperm).gt.0)then
            call prntm2(bc(itt4v(isb,iperm)),nvirt(isb),nb,nvirt(isb))           6d5s24
           end if
            ipermp=iperm+4
            write(6,*)('alternate: '),itt4v(isb,ipermp)
           if(itt4v(isb,ipermp).gt.0)then
            call prntm2(bc(itt4v(isb,ipermp)),nvirt(isb),nb,nvirt(isb))           6d5s24
           end if
           rmsd=0d0                                                     8d23s24
          end if                                                        6d20s24
          itt4vuse=-1                                                   8d27s24
          if(itt4v(isb,iperm).gt.0)then                                 8d27s24
           itt4vuse=itt4v(isb,iperm)                                    8d27s24
          else if(itt4v(isb,iperm+4).gt.0)then                          8d27s24
           itt4vuse=itt4v(isb,iperm+4)                                  8d27s24
          end if                                                        8d27s24
          if(itt4vuse.gt.0)then                                         8d27s24
           itmp=ibcoff                                                    5d29s24
           ibcoff=itmp+nvirt(isb)*nvirt(isb)                              6d5s24
           call enough                                                    5d29s24
           jorb=iorb(isb)+nb*noc(isb)                                     5d31s24
           call dgemm('n','n',nvirt(isb),nvirt(isb),nb,1d0,               6d6s24
     $        bc(itt4vuse),nvirt(isb),bc(jorb),nb,0d0,                  8d27s24
     $        bc(itmp),nvirt(isb))                                      6d6s24
           do i=0,nvirt(isb)*nvirt(isb)-1                                 6d5s24
            bc(itt4vuse+i)=bc(itmp+i)                                   8d27s24
           end do                                                         5d29s24
          end if                                                        8d27s24
          if(ldebug)then                                                6d20s24
           write(6,*)('in mo basis ')                                     5d29s24
           call prntm2(bc(itt4v(isb,iperm)),nvirt(isb),nvirt(isb),
     $         nvirt(isb))
          end if                                                        6d20s24
         end do
         if(ldebug)then                                                 8d27s24
          dif12=0d0                                                      6d11s24
          dif34=0d0                                                      6d11s24
          do i=0,nvirt(isb)*nvirt(isb)-1                                 6d11s24
           dif12=dif12+(bc(itt4v(isb,1)+i)-bc(itt4v(isb,2)+i))**2        7d16s24
           dif34=dif34+(bc(itt4v(isb,3)+i)-bc(itt4v(isb,4)+i))**2        7d16s24
          end do                                                         6d11s24
          dif12=sqrt(dif12)/dfloat(nvirt(isb))                           6d11s24
          dif34=sqrt(dif34)/dfloat(nvirt(isb))                           6d11s24
          write(6,*)('dif12: '),dif12                                    6d11s24
          write(6,*)('dif34: '),dif34                                    6d11s24
         end if                                                         8d27s24
         do i=0,nvirt(isb)*nvirt(isb)-1                                 6d11s24
          bc(itt4v(isb,ittkeep1)+i)=2d0*(bc(itt4v(isb,ittkeep1)+i)      8d27s24
     $          +bc(itt4v(isb,ittkeep3)+i))                             8d27s24
         end do                                                         6d11s24
         if(ldebug)then                                                 8d27s24
          write(6,*)('total: ')                                         8d27s24
          call prntm2(bc(itt4v(isb,ittkeep1)),nvirt(isb),nvirt(isb),    8d27s24
     $         nvirt(isb))                                              8d27s24
         end if                                                         6d20s24
         do i=0,nvirt(isb)-1                                            6d11s24
          ip=i+noc(isb)                                                 6d11s24
          do j=0,nvirt(isb)-1                                           6d11s24
           jp=j+noc(isb)                                                6d11s24
           iad1=itt4v(isb,ittkeep1)+j+nvirt(isb)*i                      8d27s24
           iad2=itt(isb)+jp+nbasdws(isb)*ip                             6d11s24
           bc(iad2)=bc(iad2)+bc(iad1)                                   6d11s24
          end do                                                        6d11s24
         end do                                                         6d11s24
        end if                                                          5d29s24
       end do                                                           5d29s24
      end if                                                            5d29s24
      call second(time2)                                                6d25s24
      telap=time2-time1                                                 6d25s24
      if(lwrite)write(6,*)('time for uu/vv and ints '),telap,telapo
      if(ldebug)then                                                    6d20s24
       call make4xdm(nsymb,nvirt,isblk4x,n4x,multh)
      else                                                              6d20s24
       n4x=0                                                            6d20s24
      end if                                                            6d20s24
      icol=ibcoff
      ibcoff=icol+mynprocg
      imsg=ibcoff                                                       6d3s10
      ibcoff=imsg+mynnode                                               6d3s10
      iter1=0
      idwsdeb=0
      call parajkfromh(ipair,ngaus,ibdat,nhcolt,ihmat,noc,icol,         11d29s12
     $     iorb,iapair,isstor,ibstor,multh,iptoh,imsg,jmats,kmats,      5d30s18
     $     ioooo,ionex,noc,iter1,idwsdeb,ncomp,nvirt,i3x,0,ih0,         7d18s23
     $     idoubo,shiftdup,nbasisp,idorel,ifull4x,i4x,i4xb,ibcb4,bc,ibc,8d8s23
     $     isend1,isend2,isend3,isend4,0,idum)                          6d14s24
      do isb=1,nsymb                                                    7d17s24
       idora(1,isb)=0                                                   7d17s24
       idora(2,isb)=0                                                   7d17s24
      end do                                                            7d17s24
      do isb=1,nsymb                                                    7d18s23
       iffp(isb,1)=ibcoff                                               7d18s23
       iffp(isb,2)=iffp(isb,1)+idora(1,isb)*idora(2,isb)                7d18s23
       iffp(isb,3)=iffp(isb,2)+idora(1,isb)*nvirt(isb)                  7d18s23
       ibcoff=iffp(isb,3)+idora(2,isb)*nvirt(isb)                       7d18s23
      end do                                                            7d18s23
      call enough('denmx12.iffp',bc,ibc)                                7d18s23
c
c     for h0 part ...
c
      call second(time1)                                                6d25s24
      if(mynowprog.eq.0)then                                            9d1s23
       jh0copy=ih0copy                                                   7d18s23
       do isb=1,nsymb                                                   9d1s23
       if(ichoice(1).ne.0)then                                          10d18s23
        do i=0,nbasdws(isb)-1
         iad1=itt(isb)+nbasdws(isb)*i
         iad2=jh0copy+nbasdws(isb)*i
         do j=0,idoubo(isb)-1
          bc(iad1+j)=4d0*bc(iad2+j)
         end do
        end do
         if(ldebug)then                                                 6d20s24
          write(6,*)('oo part of tt for sym '),isb,bc(igoal)
          call prntm2(bc(itt(isb)),nbasdws(isb),nbasdws(isb),
     $        nbasdws(isb))
         end if                                                         6d20s24
        end if
        if(ichoice(2).ne.0)then
         iad1=itt(isb)+idoubo(isb)
         iad2=jh0copy+idoubo(isb)
         if(nh0av(isb).gt.0)then                                        6d18s24
          call dgemm('n','n',nh0av(isb),nbasdws(isb),nh0av(isb),2d0,     10d19s23
     $        bc(iden(isb)),nh0av(isb),bc(iad2),nbasdws(isb),1d0,       10d19s23
     $        bc(iad1),nbasdws(isb),'denmx12.tt.h0*den')                10d19s23
          if(ldebug)then                                                6d20s24
           write(6,*)('aa part of tt for sym '),isb,bc(igoal)
           call prntm2(bc(itt(isb)),nbasdws(isb),nbasdws(isb),
     $        nbasdws(isb))
          end if                                                        6d20s24
         end if                                                         6d18s24
        end if
        nn=nbasisp(isb)*ncomp                                           10d18s23
        jh0copy=jh0copy+nn*nn                                           10d18s23
       end do                                                           9d1s23
      end if                                                            9d1s23
c
c     2e part of shift ...
c
      if(ichoice(3).ne.0)then
       if(ldebug)write(6,*)('2e part of shift ...'),loop,bc(igoal)
      call dorb2e(nsymb,idoubo,irefo,noc,nvirt,isblk,nsdlk,isblk1,      6d24s24
     $     nsdlk1,ionex,ioooo,bc,ibc,itt,nbasdws,igoal)                       6d20s24
       if(ldebug)write(6,*)('back from dorb2e '),loop,bc(igoal)
      end if
c
c     mixed double and active parts...
c
      if(ichoice(4).ne.0)then
       if(ldebug)write(6,*)('mixed da parts '),loop,bc(igoal)
       call dorbmix(nsymb,idoubo,irefo,noc,nvirt,nsdlk,isblk,           6d24s24
     $     nsdlk1,isblk1,iden,nh0av,ioooo,ionex,bc,ibc,itt,             6d17s24
     $      nbasdws,jmats,kmats,nsdlkk,isblkk,i3x)                      10d25s23
       if(ldebug)write(6,*)('back from dorbmin '),loop,bc(igoal)
      end if
c
c     id4o part ...
c
      if(ichoice(5).ne.0)then
       if(ldebug)write(6,*)('d4o '),loop,bc(igoal),loc(bc(igoal)),igoal
       call dorbd4o(nsymb,idoubo,irefo,noc,nvirt,nsdlk,isblk,           6d24s24
     $     nsdlk1,isblk1,id4o,ioooo,ionex,bc,ibc,itt,nbasdws,igoal)           6d14s24
       if(ldebug)write(6,*)('back '),loop
      end if
c
c     1dx part ...
c
      if(nsing.gt.0.and.ichoice(6).ne.0)then                            10d18s23
       if(ldebug)write(6,*)('1dx '),loop,bc(igoal)
       call dorbd1x(nsymb,idoubo,irefo,noc,nvirt,nsdlk,isblk,           6d24s24
     $     nsdlk1,isblk1,nsdlkk,isblkk,id1x,ioooo,ionex,jmats,kmats,    6d13s24
     $     bc,ibc,itt,nbasdws)                                          6d14s24
       if(ldebug)write(6,*)('back '),loop,bc(igoal)
      end if                                                            7d25s23
c
c     J and K parts
c
      if(max(nsing,ndoub).gt.0)then                                     7d31s23
       if(ichoice(7).ne.0)then                                          10d18s23
        if(ldebug)write(6,*)('dj '),loop,bc(igoal)
        call dorbdj(nsymb,idoubo,irefo,noc,nvirt,nsdlk,isblk,           6d24s24
     $     nsdlk1,isblk1,jmden,ionex,jmats,i3x,bc,ibc,itt,nbasdws)      6d14s24
        if(ldebug)write(6,*)('back '),loop
       end if                                                           10d18s23
       if(ichoice(8).ne.0)then
        if(ldebug)write(6,*)('dk '),loop,bc(igoal)
        call dorbdk(nsymb,idoubo,irefo,noc,nvirt,nsdlkk,isblkk,         6d24s24
     $     nsdlk1,isblk1,kmden,ionex,kmats,i3x,bc,ibc,itt,nbasdws)      6d14s24
        if(ldebug)write(6,*)('back '),loop
       end if                                                           10d18s23
      end if                                                            7d31s23
c                                                                       11d14s23
c     3x part                                                           11d14s23
c                                                                       11d14s23
      if(min(nsing,ndoub).gt.0.and.ichoice(9).ne.0)then                 11d14s23
       if(ldebug)write(6,*)('dorbd3x'),loop
       call dorbd3x(nsymb,idoubo,irefo,noc,nvirt,nsdlk,isblk,nsdlk1,    11d14s23
     $     isblk1,nsdlkk,isblkk,id3x,jmats,kmats,i3x,itt,nbasdws,bc,ibc)11d14s23
       if(ldebug)write(6,*)('back '),loop
      end if                                                            11d14s23
c
c     4v part
c
      if(ndoub.gt.0.and.ldebug)then                                     6d20s24
       write(6,*)('chc4v '),loop
       call chc4v(nsymb,multh,isymmrci,vdinout,nroot,nfdat,nvirt,sr2,    11d16s23
     $     bc,ibc,ndoub,dotx,i4xb,0)                                      11d16s23
       write(6,*)('dot from chc4v: '),dotx,('vs.'),dot4v,('diff'),
     $     dotx-dot4v,loop
       end if                                                           4d17s24
      if(ndoub.gt.0.and.ichoice(10).ne.0)then                           11d17s23
      if(ldebug)write(6,*)('d4x '),loop
       call dorbd4x(nsymb,idoubo,irefo,noc,nvirt,nsdlk1,isblk1,i3x,itt, 11d17s23
     $     nbasdws,isymmrci,multh,vdinout,nfdat,srh,ndoub,nroot,bc,ibc) 11d17s23
      end if                                                            11d17s23
      if(ndoub.gt.0)then                                                6d13s24
      call second(time2)                                                6d25s24
      telap=time2-time1                                                 6d25s24
      if(lwrite)write(6,*)('time for dorbs '),telap
       ivtmp=ibcoff                                                     6d5s24
       ibcoff=ivtmp+ndoub*nroot                                         6d5s24
       call enough('denmx12.vtmp',bc,ibc)                               6d5s24
       jvtmp=ivtmp-1                                                    6d5s24
       do i=1,ndoub*nroot                                               6d5s24
        bc(jvtmp+i)=vdinout(i)                                          6d5s24
       end do                                                           6d5s24
       call reordergv(bc(ivtmp),nroot,nfdat,nvirt,nsymb,multh,isymmrci, 6d5s24
     $      ndoub,1,bc,ibc)                                             6d5s24
       if(idorel.eq.0)then                                              5d29s24
        idorelu=0                                                       5d29s24
       else                                                             5d29s24
        idorelu=3                                                       5d29s24
       end if                                                           5d29s24
       call vdmo2so(bc(ivtmp),iorb,multh,ncd,nbasdwsc,iaddr,naddr,nroot,6d5s24
     $     ndoub,nbasisp,idorelu,srh,0,bc,ibc)                          6d14s24
       if(ldebug)then                                                   6d20s24
        write(6,*)('back from vdmo2so ')
        do isbv2=1,nsymb
         do isbv1=1,isbv2
          itv12=((isbv2*(isbv2-1))/2)+isbv1
          write(6,*)isbv1,isbv2,(naddr(itv12,iqq),iqq=1,3),
     $        ((iaddr(itv12,ist,ierri),ist=1,2),ierri=1,ncomp*ncomp)
         end do
        end do
       end if                                                           6d20s24
       do isbv2=1,nsymb                                                 6d13s24
        do isbv1=1,isbv2                                                6d13s24
         itv12=((isbv2*(isbv2-1))/2)+isbv1                              6d13s24
         nv1=nbasisp(isbv1)                                             6d13s24
         nv2=nbasisp(isbv2)                                             6d13s24
         do ierri=1,ncomp*ncomp                                         6d13s24
          if(isbv1.eq.isbv2.and.(ierri.eq.1.or.ierri.eq.4))then         6d13s24
           nvvs=(nv1*(nv1+1))/2                                         6d13s24
           nvvt=(nv1*(nv1-1))/2                                         6d13s24
           isw=0                                                        6d14s24
          else                                                          6d13s24
           nvvs=nv1*nv2                                                 6d13s24
           nvvt=nvvs                                                    6d13s24
           isw=1                                                        6d14s24
          end if                                                        6d13s24
          nvv=nvvs                                                      6d13s24
          isub=+1                                                       6d18s24
          do ist=1,2
           if(naddr(itv12,ist).gt.0)then
            i1=iaddr(itv12,ist,ierri)
            i2=nvv
            i3=naddr(itv12,ist)
            itmpx=ibcoff                                                 6d11s24
            ibcoff=itmpx+nvv*naddr(itv12,ist)                            6d11s24
            call enough('denmx12.tmpxb',bc,ibc)                          6d11s24
            if(ist.eq.1.and.isw.eq.0)then                               6d14s24
             do ii=0,naddr(itv12,ist)-1                                 6d14s24
              do iv=0,nbasisp(isbv1)-1                                  6d14s24
               ivv=((iv*(iv+1))/2)+iv                                   6d14s24
               iad=iaddr(itv12,ist,ierri)+ivv+nvv*ii                    6d14s24
               bc(iad)=bc(iad)*0.5d0                                    6d14s24
              end do                                                    6d14s24
             end do                                                     6d14s24
            end if                                                      6d14s24
            do ii=0,naddr(itv12,ist)-1
             do ivv=0,nvv-1
              iad=iaddr(itv12,ist,ierri)+ivv+nvv*ii
              jtmpx=itmpx+ii+naddr(itv12,ist)*ivv                        6d11s24
              bc(jtmpx)=bc(iad)                                          6d11s24
             end do
            end do
            do ii=0,naddr(itv12,ist)*nvv-1                              6d11s24
             bc(iaddr(itv12,ist,ierri)+ii)=bc(itmpx+ii)                 6d11s24
            end do                                                      6d11s24
            do iv2=0,nbasisp(isbv2)-1                                   6d18s24
             iv1top=iv2+1-ist                                           6d18s24
             iv1top=iv1top+isw*(nbasisp(isbv1)-1-iv1top)                6d18s24
             do iv1=0,iv1top                                            6d18s24
              irec=iv1+nbasisp(isbv1)*iv2                               6d18s24
              itri=((iv2*(iv2+isub))/2)+iv1                             6d18s24
              itri=itri+isw*(irec-itri)                                 6d18s24
              do ii=0,naddr(itv12,ist)-1                                6d18s24
               iad=iaddr(itv12,ist,ierri)+ii+naddr(itv12,ist)*itri      6d18s24
              end do                                                    6d18s24
             end do                                                     6d18s24
            end do                                                      6d18s24
            ibcoff=itmpx                                                6d11s24
           end if                                                       6d13s24
           nvv=nvvt                                                     6d13s24
           isub=-1                                                      6d18s24
          end do
         end do
        end do                                                          6d13s24
       end do                                                           6d13s24
      else                                                              6d13s24
       iaddr(1,1,1)=-1                                                  6d13s24
      end if                                                            6d13s24
      jdarot=idarot                                                     7d18s23
      ithit=0
      do ider=0,nder-1                                                  7d18s23
       ixyz=ibc(jdarot)                                                 7d13s23
       ia=ibc(jdarot+1)                                                 7d13s23
       ipass=ibc(jdarot+2)                                              7d13s23
       ia2=ibc(jdarot+3)                                                7d13s23
       if(ia2.gt.0)then                                                 7d11s24
        npass=2                                                         7d11s24
       else                                                             7d11s24
        npass=1                                                         7d11s24
       end if                                                           7d11s24
       if(ldebug)                                                       6d20s24
     $     write(6,*)('derivative code: '),ixyz,ia,ipass,ia2,jdarot,loop6d20s24
       jdarot=jdarot+4                                                  7d18s23
       jdarot0=jdarot
       if(mynprocg.eq.1.and.ldebug)then                                 6d20s24
        call playr(bc(jdarot),idora,iden,ih0copy,idoubo,irefo,noc,nvirt, 8d4s23
     $      nbasdws,nh0av,nbasisp,nsymb,isblk,nsdlk,isblk1,nsdlk1,      8d4s23
     $      ionex,ioooo,jmats,kmats,i3x,nsdlkk,isblkk,n4x,isblk4x,i4xb, 8d8s23
     $      multh,id4o,id1x,jmden,kmden,ipass,bc,ibc,itt,ichoice,       11d14s23
     $      id3x,isymmrci,vdinout,nfdat,sr2,ndoub,potdws,ncomp)         6d14s24
       end if                                                           4d18s24
       ibcxp=ibcoff                                                     6d11s24
       if(ipass.lt.0)then                                               4d25s24
        dnuc=0d0                                                        4d29s24
        if(ixyz.le.3)then                                               4d25s24
         npt=1                                                          4d29s24
         do i=1,7                                                       4d29s24
          idatta(i,1)=0                                                 4d29s24
         end do                                                         4d29s24
         idatta(ixyz,1)=1                                               4d29s24
         data(1)=-1d0                                                   4d29s24
         do ia=1,natom                                                  4d29s24
          dnuc=dnuc+atnum(1,ia)*xcart(ixyz,ia)                          4d29s24
         end do                                                         4d29s24
         write(dlabel,885)char(ichar('w')+ixyz)                         4d25s24
        else if(ixyz.eq.4)then                                          4d29s24
         dlabel='Q0 field  '
         npt=3                                                          4d29s24
         do j=1,3                                                       4d29s24
          do i=1,7                                                      4d29s24
           idatta(i,j)=0                                                4d29s24
          end do                                                        4d29s24
         end do                                                         4d29s24
         idatta(1,1)=2                                                  4d29s24
         data(1)=1d0/sqrt(6d0)                                          4d29s24
         idatta(2,2)=2                                                  4d29s24
         data(2)=1d0/sqrt(6d0)                                          4d29s24
         idatta(3,3)=2                                                  4d29s24
         data(3)=-2d0/sqrt(6d0)                                         4d29s24
         do ia=1,natom                                                  4d29s24
          dnuc=dnuc+atnum(1,ia)*(2d0*xcart(3,ia)*xcart(3,ia)            4d29s24
     $           -xcart(1,ia)*xcart(1,ia)-xcart(2,ia)*xcart(2,ia))      4d29s24
         end do                                                         4d29s24
         dnuc=dnuc/sqrt(6d0)                                            4d29s24
        else if(ixyz.eq.5)then                                          4d29s24
         dlabel='ReQ2 field'
         npt=2                                                          4d29s24
         do j=1,2                                                       4d29s24
          do i=1,7                                                      4d29s24
           idatta(i,j)=0                                                4d29s24
          end do                                                        4d29s24
         end do                                                         4d29s24
         idatta(1,1)=2                                                  4d29s24
         data(1)=-0.5d0                                                 4d29s24
         idatta(2,2)=2                                                  4d29s24
         data(2)=0.5d0                                                  3d31s23
         do ia=1,natom                                                  4d29s24
          dnuc=dnuc+atnum(1,ia)*(xcart(1,ia)*xcart(1,ia)                4d29s24
     $         -xcart(2,ia)*xcart(2,ia))                                4d29s24
         end do                                                         4d29s24
         dnuc=dnuc*0.5d0                                                4d29s24
        else if(ixyz.eq.6)then                                          4d29s24
         npt=1                                                          4d29s24
         dlabel='ImQ2 field'                                            4d25s24
         do i=1,7                                                       4d29s24
          idatta(i,1)=0                                                 4d29s24
         end do                                                         4d29s24
         idatta(1,1)=1                                                  4d29s24
         idatta(2,1)=1                                                  4d29s24
         data(1)=-1d0                                                   4d29s24
         do ia=1,natom                                                  4d29s24
          dnuc=dnuc+atnum(1,ia)*xcart(1,ia)*xcart(2,ia)                 4d29s24
         end do                                                         4d29s24
        else if(ixyz.eq.7)then                                          4d29s24
         dlabel='ReQ1 field'                                            4d25s24
         npt=1                                                          4d29s24
         do i=1,7                                                       4d29s24
          idatta(i,1)=0                                                 4d29s24
         end do                                                         4d29s24
         idatta(1,1)=1                                                  4d29s24
         idatta(3,1)=1                                                  4d29s24
         data(1)=1d0                                                    4d29s24
         do ia=1,natom                                                  4d29s24
          dnuc=dnuc+atnum(1,ia)*xcart(1,ia)*xcart(3,ia)                 4d29s24
         end do                                                         4d29s24
         dnuc=-dnuc                                                     4d29s24
        else                                                            4d29s24
         dlabel='ImQ1 field'                                            4d25s24
         npt=1                                                          4d29s24
         do i=1,7                                                       4d29s24
          idatta(i,1)=0                                                 4d29s24
         end do                                                         4d29s24
         idatta(2,1)=1                                                  4d29s24
         idatta(3,1)=1                                                  4d29s24
         data(1)=1d0                                                    4d29s24
         do ia=1,natom                                                  4d29s24
          dnuc=dnuc+atnum(1,ia)*xcart(2,ia)*xcart(3,ia)                 4d29s24
         end do                                                         4d29s24
         dnuc=-dnuc                                                     4d29s24
        end if                                                          4d29s24
        if(ldebug)write(6,*)('dnuc = '),dnuc                            6d20s24
        ifmat=ibcoff                                                    4d29s24
        do isb=1,nsymb                                                  4d29s24
         iamatu(isb)=ibcoff                                             4d29s24
         ibcoff=iamatu(isb)+nbasisp(isb)*nbasisp(isb)*ncomp*ncomp       4d29s24
        end do                                                          4d29s24
        nbb=ibcoff-ifmat                                                4d29s24
        call enough('denmx12.amatu',bc,ibc)                             4d29s24
        do iz=ifmat,ibcoff-1                                            4d29s24
         bc(iz)=0d0                                                     4d29s24
        end do                                                          4d29s24
        iptxx=1                                                         4d29s24
        call parap(natom,ngaus,ibdat,iamatu,isym,iapair,ibstor,         4d29s24
     $        isstor,idum,0,idorel,ascale,iptxx,npt,data,idatta,        8d30s22
     $        1,1,multh,nbb,nbasisp,nbasdws,iorb,dlabel,1,bc,ibc)       4d29s24
        dcore=0d0                                                       5d2s24
        drest=0d0                                                       4d29s24
        do isb=1,nsymb                                                  4d29s24
         if(ldebug)then                                                 6d20s24
          write(6,*)('matrix elements for isb = '),isb
          call prntm2(bc(iamatu(isb)),nbasdws(isb),nbasdws(isb),
     $        nbasdws(isb))
         end if                                                         6d20s24
         do i=0,idoubo(isb)-1                                           4d29s24
          iad=iamatu(isb)+i*(nbasdws(isb)+1)                            4d29s24
          dcore=dcore+bc(iad)*2d0                                       4d29s24
         end do                                                         4d29s24
         do i=0,nh0av(isb)-1                                            4d29s24
          ip=i+idoubo(isb)                                              4d29s24
          iad1=iamatu(isb)+idoubo(isb)+nbasdws(isb)*ip                  4d29s24
          iad2=iden(isb)+nh0av(isb)*i                                   4d29s24
          do j=0,nh0av(isb)-1                                           4d29s24
           drest=drest+bc(iad1+j)*bc(iad2+j)                            4d29s24
          end do                                                        4d29s24
         end do                                                         4d29s24
        end do                                                          4d29s24
        ddd=dnuc+dcore+drest                                            4d29s24
        if(ldebug)                                                      6d20s24
     $       write(6,*)('expectation value: nuc '),dnuc,('core '),dcore,6d20s24
     $       ('rest '),drest,('sum: '),ddd                              4d29s24
        iad=ifder+ider+nder*4                                           4d29s24
        bc(iad)=ddd/dfloat(mynprocg)                                    5d16s24
       else                                                             4d25s24
        derpotn=0d0                                                     6d11s24
        idersign=1                                                      6d11s24
        if(ia2.ne.0.and.ipass.eq.2)idersign=2                           6d11s24
        do ixpass=1,npass                                               7d11s24
         sig=1d0                                                        6d11s24
         if(ixpass.eq.1)then                                            6d11s24
          ja=ia                                                         6d11s24
         else                                                           6d11s24
          ja=ia2                                                        6d11s24
          if(idersign.eq.2)sig=-1d0                                     6d11s24
         end if                                                         6d11s24
         do iai=1,natom                                                 6d11s24
          if(iai.ne.ja)then                                             6d11s24
           dist=0d0                                                     6d11s24
           do jxyz=1,3                                                  6d11s24
            dist=dist+(xcart(jxyz,iai)-xcart(jxyz,ja))**2               6d11s24
           end do                                                       6d11s24
           dist=1d0/dist                                                6d11s24
           dists=sqrt(dist)                                             6d11s24
           vij=atnum(1,iai)*atnum(1,ja)*dists                           6d11s24
           dvij=(xcart(ixyz,iai)-xcart(ixyz,ja))*vij*dist               6d11s24
           derpotn=derpotn+drsign*dvij*sig                              8d20s24
          end if                                                        6d11s24
         end do                                                         6d11s24
        end do                                                          6d11s24
        if(ldebug)write(6,*)('derpotn = '),derpotn                                6d11s24
        nnp=0                                                           6d11s24
        nn=0                                                            6d11s24
        do isb=1,nsymb                                                  6d11s24
         nh=nbasdws(isb)                                                6d11s24
         nhp=nbasisp(isb)*ncomp                                         6d11s24
         isou(isb)=nnp                                                  6d11s24
         isoub(isb)=nn                                                  6d11s24
         isoua1(isb)=nnp                                                6d11s24
         nnp=nnp+nhp*nhp                                                6d11s24
         nn=nn+nh*nh                                                    6d11s24
        end do                                                          6d11s24
        nn1=nnp                                                         6d11s24
        ih0d=ibcoff                                                     6d11s24
        iovr=ih0d+nnp*2                                                 6d11s24
        ibcoff=iovr+nnp*2                                               6d11s24
        iovrdk=ibcoff                                                   6d11s24
        npropmat=iovrdk+nnp*2                                           6d11s24
        ibcoff=npropmat+nnp                                             6d11s24
        call enough('denmx12.h0d',bc,ibc)                               6d11s24
        iovrdd=ibcoff                                                   6d11s24
        idwsdeb=000                                                     6d13s24
        call dws_synca                                                  6d13s24
        call parah0grad(natom,ngaus,ibdat,nbasis,bc(ih0d),bc(iovr),     6d11s24
     $         bc(iovrdk),isym,iapair,ibstor,isstor,isou,nnp,           6d11s24
     $            idwsdeb,idorel,ascale,multh,ixyz,ia,1,idersign,       6d11s24
     $       nbasisp,.false.,bc(iovrdd),isoua1,nn1,bc,ibc)              6d11s24
        idwsdeb=0                                                       6d11s24
        dcore=0d0                                                       5d2s24
        drest=0d0                                                       4d29s24
        jh0d=ih0d                                                       6d11s24
        do isb=1,nsymb                                                  4d29s24
         if(nbasdws(isb).gt.0)then                                      6d11s24
          npdn=nbasisp(isb)*ncomp                                       6d11s24
          if(ldebug)then                                                6d20s24
           write(6,*)('for symmetry block '),isb
           write(6,*)('h0d'),jh0d                                        6d11s24
           call mpprnt2(bc(jh0d),npdn)                                   6d11s24
          end if                                                        6d20s24
          call square(bc(jh0d),npdn)                                    6d11s24
          if(ldebug)then                                                6d20s24
           write(6,*)('squared ')
           call prntm2(bc(jh0d),npdn,npdn,npdn)
          end if                                                        6d20s24
          itmph=ibcoff                                                  6d11s24
          ibcoff=itmph+npdn*nbasdws(isb)                                6d11s24
          call enough('denmx12.tmph',bc,ibc)                            6d11s24
          nrow=npdn                                                     6d11s24
          do ipx=1,2                                                    6d11s24
           call dgemm('n','n',nrow,nbasdws(isb),npdn,1d0,               6d11s24
     $          bc(jh0d),nrow,bc(iorb(isb)),npdn,0d0,                   6d11s24
     $          bc(itmph),nrow,'denmx12.tmph')                          6d11s24
           do i=0,nbasdws(isb)-1                                        6d11s24
            do j=0,nrow-1                                               6d11s24
             ji=itmph+j+nrow*i                                          6d11s24
             ij=jh0d+i+nbasdws(isb)*j                                   6d11s24
             bc(ij)=bc(ji)                                              6d11s24
            end do                                                      6d11s24
           end do                                                       6d11s24
           nrow=nbasdws(isb)                                            6d11s24
          end do                                                        6d11s24
          if(ldebug)then                                                6d20s24
           write(6,*)('in mo basis ')
           call prntm2(bc(jh0d),nbasdws(isb),nbasdws(isb),nbasdws(isb))  6d11s24
           write(6,*)('1e density: ')
           call prntm2(bc(iden(isb)),nh0av(isb),nh0av(isb),nh0av(isb))
          end if                                                        6d20s24
          do i=0,idoubo(isb)-1                                           4d29s24
           iad=jh0d+i*(nbasdws(isb)+1)                                  6d11s24
           dcore=dcore+bc(iad)*2d0                                       4d29s24
          end do                                                         4d29s24
          do i=0,nh0av(isb)-1                                            4d29s24
           ip=i+idoubo(isb)                                              4d29s24
           iad1=jh0d+idoubo(isb)+nbasdws(isb)*ip                        6d11s24
           iad2=iden(isb)+nh0av(isb)*i                                   4d29s24
           do j=0,nh0av(isb)-1                                           4d29s24
            drest=drest+bc(iad1+j)*bc(iad2+j)                            4d29s24
           end do                                                        4d29s24
          end do                                                         4d29s24
          jh0d=jh0d+npdn*npdn                                           6d11s24
         end if                                                         6d11s24
        end do                                                          6d11s24
        ddd=derpotn+dcore+drest                                         6d11s24
        if(ldebug)then                                                  6d20s24
         write(6,*)('expectation value: nuc '),derpotn,('1e core '),     6d11s24
     $       dcore,('1e rest '),drest,('sum: '),ddd                     6d11s24
        end if                                                          6d20s24
        iad=ifder+ider+nder*4                                           4d29s24
        ibcb4=ibcoff                                                    6d12s24
        do i=1,4                                                        6d13s24
         d4v(i)=0d0                                                     6d13s24
        end do                                                          6d13s24
        if(ldebug)then                                                  6d20s24
         ifull4x=3                                                       6d14s24
        else                                                            6d20s24
         ifull4x=2                                                      6d20s24
        end if                                                          6d20s24
        call second(time1)
        call paraeridd(natom,ngaus,ibdat,idum,ihmat,iorb,noc,           6d11s24
     $             ipair,nhcolt,isym,iapair,ibstor,isstor,multh,iptoh,  11d29s12
     $     0,idwsdeb,idorel,ascale,nbasisp,iaddr,naddr,ixyz,ia,         6d11s24
     $       idersign,d4v,ifull4x,bc,ibc)                               6d14s24
        call second(time2)
        telap=time2-time1
        if(lwrite)write(6,*)('time for paraeridd '),telap
        d4vsum=0d0                                                      6d14s24
        if(ndoub.gt.0)then                                              6d14s24
         if(ldebug)write(6,*)('what we have for d4v: '),d4v             6d20s24
         call dws_gsumf(d4v,4)                                           6d13s24
         d4vsum=d4v(1)+d4v(2)+d4v(3)+d4v(4)                              6d14s24
         if(ldebug)then                                                 6d20s24
          write(6,*)('global summed '),d4v                                6d13s24
          write(6,*)('summed '),d4vsum
         end if                                                         6d20s24
        end if                                                          6d14s24
        dshift=0d0                                                      6d12s24
        call parajkfromh(ipair,ngaus,ibdat,nhcolt,ihmat,noc,icol,         11d29s12
     $     iorb,iapair,isstor,ibstor,multh,iptoh,imsg,jmtsd,kmtsd,      7d11s24
     $     ioood,ionxd,irefo,iter1,idwsdeb,ncomp,nvirt,i3xxd,1,ih0d,    7d11s24
     $     idoubo,dshift,nbasisp,idorel,ifull4x,i4x,i4xbd,ibcb4,bc,ibc,  6d14s24
     $     isend1,isend2,isend3,isend4,0,idum)                          6d14s24
        if(ndoub.gt.0.and.ldebug)then                                   6d20s24
         dotq=0d0                                                        6d14s24
         call chc4v(nsymb,multh,isymmrci,vdinout,nroot,nfdat,nvirt,sr2,    11d16s23
     $     bc,ibc,ndoub,dotq,i4xbd,1)                                   7d11s24
         write(6,*)('dotq from chc4v with der ints: '),dotq
         call dws_gsumf(dotq,1)
         write(6,*)('global summed: '),dotq
         write(6,*)('diff: '),dotq-d4vsum
         write(6,*)('what we have for dshift: '),dshift
        end if                                                          6d14s24
        call trans2e(jmtsd,kmtsd,ionxd,i3xxd,i3xb,irefo,iunc(2),ionexb, 7d11s24
     $     jmatt,kmatt,ionexc,kmatd,i3x3,ionexbt,bc,ibc)                11d10s22
        jh0d=ih0d                                                       6d12s24
        dact=0d0                                                        6d12s24
        do isb=1,nsymb                                                  6d12s24
         if(nh0av(isb).gt.0)then                                        6d12s24
          if(ldebug)then                                                6d20s24
           write(6,*)('folded h0 for sym '),isb,jh0d                     6d12s24
           call prntm2(bc(jh0d),nh0av(isb),nh0av(isb),nh0av(isb))        6d12s24
          end if                                                        6d20s24
          do i=0,nh0av(isb)*nh0av(isb)-1                                6d12s24
           dact=dact+bc(jh0d+i)*bc(iden(isb)+i)                         6d12s24
          end do                                                        6d12s24
          jh0d=jh0d+nh0av(isb)*nh0av(isb)                               6d12s24
         end if                                                         6d12s24
        end do                                                          6d12s24
        ddd=derpotn+dshift+dact                                         6d12s24
        if(ldebug)then                                                  6d20s24
         write(6,*)('expectation value: nuc '),derpotn,('shift '),       6d12s24
     $       dshift,('1e '),dact,('sum: '),ddd                          6d12s24
        end if                                                          6d20s24
        iad=ifder+ider+nder*4                                           4d29s24
        bc(iad)=ddd/dfloat(mynprocg)                                    5d16s24
        dot2ea=0d0                                                       6d12s24
        dot2ed=0d0                                                      6d12s24
        do is=1,nsdlk                                                     7d10s23
         if(isblk(1,is).eq.isblk(2,is))then                               7d28s22
          nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              7d28s22
          ncol=(irefo(isblk(3,is))*(irefo(isblk(3,is))+1))/2              7d28s22
          isw=0                                                         6d12s24
         else                                                             7d28s22
          nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      7d28s22
          ncol=irefo(isblk(3,is))*irefo(isblk(4,is))                       7d28s22
          isw=1                                                         6d12s24
         end if                                                           7d28s22
         call ilimts(irefo(isblk(3,is)),irefo(isblk(4,is)),mynprocg,    6d12s24
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          6d12s24
         nhere=ih+1-il                                                  6d12s24
         if(min(nrow,ncol,nhere).gt.0)then                              6d12s24
          do i2=0,irefo(isblk(4,is))-1                                  6d12s24
           i1top=i2+isw*(irefo(isblk(3,is))-1-i2)                       6d12s24
           do i1=0,i1top                                                6d12s24
            icol=i1+1+irefo(isblk(3,is))*i2                             6d12s24
            if(icol.ge.il.and.icol.le.ih)then                           6d12s24
             icol0=icol
             icol=ioood(is)+nrow*(icol-il)                              7d11s24
             irec=i1+irefo(isblk(3,is))*i2                              6d12s24
             itri=((i2*(i2+1))/2)+i1                                    6d12s24
             itri=itri+isw*(irec-itri)                                  6d12s24
             idbase=id4o(is)+nrow*itri                                  6d12s24
             do i4=0,irefo(isblk(2,is))-1                               6d12s24
              i3top=i4+isw*(irefo(isblk(1,is))-1-i4)                    6d12s24
              jrec=irefo(isblk(1,is))*i4                                6d12s24
              jtri=((i4*(i4+1))/2)                                      6d12s24
              jtri=jtri+isw*(jrec-jtri)                                 6d12s24
              ktri=jtri+icol                                            6d12s24
              jtri=jtri+idbase                                          6d12s24
              do i3=0,i3top                                             6d12s24
               dot2ed=dot2ed+bc(jtri+i3)*bc(ktri+i3)                    6d12s24
              end do                                                    6d12s24
             end do                                                     6d12s24
            end if                                                      6d12s24
           end do                                                       6d12s24
          end do                                                        6d12s24
         end if                                                         6d12s24
         if(max(ndoub,nsing).gt.0.and.nrow.gt.0)then                    6d12s24
          call ilimts(nvirt(isblk(3,is)),nvirt(isblk(4,is)),mynprocg,   10d16s23
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           10d16s23
          nhere=ih+1-il                                                 6d12s24
          if(nhere.gt.0)then                                            6d12s24
           do ii=0,nrow*nhere-1                                         6d12s24
            dot2ed=dot2ed+bc(jmden(is)+ii)*bc(jmtsd(is)+ii)             7d11s24
           end do                                                       6d12s24
          end if                                                        6d12s24
         end if                                                         6d12s24
        end do                                                          6d12s24
        if(ldebug)then                                                  6d20s24
         write(6,*)('dot2ea after 4o,J: '),dot2ea,dot2ed                 6d12s24
         xx=dot2ed
         call dws_gsumf(xx,1)
         write(6,*)('gsummed '),xx
        end if                                                          6d20s24
        if(max(ndoub,nsing).gt.0)then                                   6d12s24
         do is=1,nsdlk1                                                    7d28s22
          if(isblk1(1,is).eq.isblk1(2,is))then                             7d28s22
           nrow1=(irefo(isblk1(1,is))*(irefo(isblk1(1,is))+1))/2           7d28s22
           nrow3=(nvirt(isblk1(1,is))*(nvirt(isblk1(1,is))+1))/2          10d16s23
          else                                                             7d28s22
           nrow1=irefo(isblk1(1,is))*irefo(isblk1(2,is))                   7d28s22
           nrow3=nvirt(isblk1(1,is))*nvirt(isblk1(2,is))                  10d16s23
          end if                                                           7d28s22
          ncol=irefo(isblk1(3,is))*nvirt(isblk1(4,is))
          call ilimts(irefo(isblk1(3,is)),nvirt(isblk1(4,is)),mynprocg,   10d16s23
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           10d16s23
          nhere=ih+1-il                                                 6d13s24
          if(min(nhere,nrow1,ncol).gt.0)then                            6d13s24
c
c     id1x was nvirt(4),nrow1,irefo(3) but dorbd1x transposed it to
c     nrow1,irefo(3),nvirt(4)
           iad=id1x(is)+nrow1*(il-1)                                    6d13s24
           do ii=0,nrow1*nhere-1                                        6d13s24
            dot2ed=dot2ed+bc(iad+ii)*bc(ionxd(is)+ii)                   7d11s24
           end do                                                          7d10s23
          end if                                                          10d16s23
          if(min(ndoub,nsing,nrow3,nhere).gt.0)then                     6d13s24
           do ii=0,nrow3*nhere-1                                        6d13s24
            dot2ed=dot2ed+bc(id3x(is)+ii)*bc(i3xxd(is)+ii)              7d11s24
           end do                                                          7d10s23
          end if                                                        6d12s24
         end do                                                            7d28s22
         if(ldebug)then                                                 6d20s24
          write(6,*)('dot2e after 1x,3x: '),dot2ea,dot2ed                6d12s24
          xx=dot2ed
          call dws_gsumf(xx,1)
          write(6,*)('gsummed '),xx
         end if                                                         6d20s24
         do is=1,nsdlkk                                                 6d12s24
          nrow=irefo(isblkk(1,is))*irefo(isblkk(2,is))                  6d12s24
          call ilimts(nvirt(isblkk(3,is)),nvirt(isblkk(4,is)),mynprocg, 6d13s24
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           10d16s23
          ncol=ih+1-il                                                    10d16s23
          if(min(nrow,ncol).gt.0)then                                   6d12s24
           do ii=0,nrow*ncol-1                                            7d10s23
            dot2ed=dot2ed+bc(kmden(is)+ii)*bc(kmtsd(is)+ii)             7d11s24
           end do                                                          7d10s23
          end if                                                        6d12s24
         end do                                                         6d12s24
         if(ldebug)write(6,*)('dot2e after K: '),dot2ed                 6d20s24
        end if                                                          6d12s24
        dot2ea=dot2ea+d4vsum                                            6d14s24
        call dws_gsumf(dot2ed,1)                                        6d12s24
        dot2e=dot2ea+dot2ed                                             6d12s24
        if(ldebug)write(6,*)('total 2e part: '),dot2e                   6d20s24
        iad=ifder+ider+nder*5                                           6d12s24
        bc(iad)=dot2e/dfloat(mynprocg)                                  6d12s24
        idwsdeb=0
       end if                                                           4d25s24
       ibcoff=ibcxp                                                     6d11s24
       ttranst=0d0
       ndaspace=4                                                       7d20s23
       trcsum=0d0                                                       10d22s23
       do isb=1,nsymb                                                   7d18s23
        if(ldebug)write(6,*)('for symmetry '),isb                       6d20s24
        ndaspace=ndaspace+nbasdws(isb)*nbasdws(isb)                     7d20s23
        trc=0d0                                                         10d18s23
        if(ldebug)then                                                  8d6s24
         write(6,*)('darot '),jdarot,jdarot-jdarot0+1
         call prntm2(bc(jdarot),nbasdws(isb),nbasdws(isb),nbasdws(isb))
         write(6,*)('tt '),itt(isb),bc(igoal),igoal-itt(isb)
         call prntm2(bc(itt(isb)),nbasdws(isb),nbasdws(isb),
     $        nbasdws(isb))
        end if                                                          8d6s24
        do i=0,nbasdws(isb)-1
         do j=0,nbasdws(isb)-1
          ji=itt(isb)+j+nbasdws(isb)*i
          ij=jdarot+i+nbasdws(isb)*j
          term=bc(ji)*bc(ij)
          orig=trc
          trc=trc+term
         end do
        end do
        call dws_gsumf(trc,1)
        trcsum=trcsum+trc                                               10d22s23
        if(ldebug)write(6,*)('trace with tt '),trc,trcsum               6d20s24
        jdarot=jdarot+nbasdws(isb)*nbasdws(isb)                         7d18s23
        if(ipass.gt.0)then                                              10d18s23
         if(ithit.le.nsymb.and.ldebug)then                              6d20s24
          write(6,*)('tt for transder: ')
          call prntm2(bc(itt(isb)),nbasdws(isb),nbasdws(isb),
     $         nbasdws(isb))
          write(6,*)('darot '),jdarot,jdarot-jdarot0+1
          call prntm2(bc(jdarot),nbasdws(isb),nbasdws(isb),nbasdws(isb))
          ithit=ithit+1
         end if                                                         10d18s23
         trct=0d0
         do i=0,nbasdws(isb)-1
          do j=0,nbasdws(isb)-1
           ji=itt(isb)+j+nbasdws(isb)*i
           ij=jdarot+i+nbasdws(isb)*j
           term=bc(ji)*bc(ij)
           orig=trct
           trct=trct+term
          end do
         end do
         call dws_gsumf(trct,1)
         if(ldebug)write(6,*)('trace for transder'),trct                6d20s24
         ttranst=ttranst+trct
         jdarot=jdarot+nbasdws(isb)*nbasdws(isb)                        10d18s23
        end if                                                          10d18s23
       end do                                                           7d18s23
       if(ipass.ge.0)then
        if(ldebug)then                                                  6d20s24
         write(6,*)('total tt: '),trcsum
         write(6,*)('total ttrans: '),ttranst
         write(6,*)('total orb+basis: '),trcsum+ttranst
        end if                                                          6d20s24
        iadr=ifder+ider+nder*3                                          4d25s24
        bc(iadr)=ttranst/dfloat(mynprocg)                               4d25s24
       end if
       iadr=ifder+ider+nder*2                                           4d25s24
       bc(iadr)=trcsum/dfloat(mynprocg)                                 4d25s24
      end do
      call dws_gsumf(bc(ifder),nder*6)                                  4d25s24
      if(lwrite)then                                                    6d20s24
       write(6,*)('>gsummed fder ...')
       write(6,889)                                                      4d25s24
  889  format(15x,'dcont',10x,'dcorth',9x,'dorb',11x,'doorth',9x,'1e',     4d25s24
     $     13x,'2e',18x,'total')                                                     4d25s24
       jdarot=idarot                                                     7d18s23
       do ix=1,nder                                                      4d25s24
        ixyz=ibc(jdarot)                                                 7d13s23
        ia=ibc(jdarot+1)                                                 7d13s23
        ipass=ibc(jdarot+2)                                              7d13s23
        ia2=ibc(jdarot+3)                                                7d13s23
        jdarot=jdarot+4                                                  7d18s23
        if(ipass.lt.0)then                                               4d25s24
         do isb=1,nsymb                                                  4d25s24
          jdarot=jdarot+nbasdws(isb)*nbasdws(isb)                         7d18s23
         end do                                                          4d25s24
         if(ixyz.le.3)then                                               8d3s23
          write(dlabel,885)char(ichar('w')+ixyz)                         4d25s24
  885     format('mu',a1,' field ')                                      4d25s24
         else if(ixyz.eq.4)then                                          8d3s23
          dlabel='Q0 field  '
         else if(ixyz.eq.5)then                                          8d3s23
          dlabel='ReQ2 field'
         else if(ixyz.eq.6)then                                          8d3s23
          dlabel='ImQ2 field'                                            4d25s24
         else if(ixyz.eq.7)then                                          8d3s23
          dlabel='ReQ1 field'                                            4d25s24
         else                                                            8d3s23
          dlabel='ImQ1 field'                                            4d25s24
         end if                                                          8d3s23
        else                                                             7d13s23
         do isb=1,nsymb                                                  4d25s24
          jdarot=jdarot+2*nbasdws(isb)*nbasdws(isb)                      4d25s24
         end do                                                          4d25s24
         if(ia2.eq.0)then                                                7d13s23
          write(dlabel,886)char(ichar('w')+ixyz),ia                      4d25s24
  886     format('d/d',a1,i2)                                           8d9s24
         else                                                            7d13s23
          if(ipass.eq.1)then                                             7d13s23
           write(dlabel,887)char(ichar('w')+ixyz),ia,('+'),ia2           4d25s24
  887      format('d/d',a1,i2,a1,i2)                                    8d9s24
          else                                                           7d13s23
           write(dlabel,887)char(ichar('w')+ixyz),ia,('-'),ia2           4d25s24
          end if                                                         7d13s23
         end if                                                          7d13s23
        end if                                                           7d13s23
        ixm=ix-1                                                         4d25s24
        sum=0d0                                                          4d25s24
        jfder=ifder+ixm                                                  4d25s24
        do j=0,5                                                         4d25s24
         sum=sum+bc(jfder)                                               4d25s24
         jfder=jfder+nder                                                4d25s24
        end do                                                           4d25s24
        jfder=ifder+ixm                                                  4d25s24
        write(6,888)dlabel,(bc(jfder+j*nder),j=0,5),sum                  4d25s24
  888   format(a10,x,6es15.7,4x,'>',es22.14)                              4d25s24
       end do                                                            4d25s24
      end if                                                            6d20s24
      return
      end
      subroutine fiddleda(idarot,fder,nder,ndaspace,bc,ibc)
      implicit real*8 (a-h,o-z)
      dimension fder(*)
      include "common.store"
      jdarot=idarot                                                     7d20s23
      do i=1,nder                                                       7d21s23
       ixyz=ibc(jdarot)                                                 7d13s23
       ia=ibc(jdarot+1)                                                 7d13s23
       ipass=ibc(jdarot+2)                                              7d13s23
       ia2=ibc(jdarot+3)                                                7d13s23
       jdarot=jdarot+ndaspace                                           7d20s23
       kdarot=jdarot                                                    7d21s23
       ffi=fder(i)                                                      7d21s23
       do k=i+1,nder                                                    7d21s23
        kxyz=ibc(kdarot)                                                 7d13s23
        ka=ibc(kdarot+1)                                                 7d13s23
        kpass=ibc(kdarot+2)                                              7d13s23
        ka2=ibc(kdarot+3)                                                7d13s23
        kdarot=kdarot+ndaspace                                          7d21s23
        ffk=fder(k)                                                     7d21s23
        if(ixyz.eq.kxyz.and.ipass.gt.0.and.                             7d21s23
     $        min(abs(ffi),abs(ffk)).gt.1d-5)then                       7d21s23
          if(abs(abs(ffi)-abs(ffk)).lt.1d-5)then                        7d21s23
           sum=ffi+ffk                                                  7d21s23
           dif=ffi-ffk                                                  7d21s23
           write(6,*)('sum,diff of '),ixyz,ia,ka,ffi,                   7d21s23
     $          ffk,(' is '),sum,dif                                    7d21s23
          end if                                                        7d21s23
        end if                                                          7d21s23
       end do                                                           7d20s23
       write(6,*)ixyz,ia,ipass,ia2,ffi                                  7d21s23
      end do                                                            7d20s23
      return
      end
      subroutine playr(darot,idora,iden,ih0copy,idoubo,irefo,noc,nvirt, 8d4s23
     $     nbasdws,nh0av,nbasisp,nsymb,isblk,nsdlk,isblk1,nsdlk1,       8d4s23
     $     ionex,ioooo,jmats,kmats,i3x,nsdlkk,isblkk,n4x,isblk4x,i4x,   8d8s23
     $     multha,id4o,id1x,jmden,kmden,ipass,bc,ibc,itt,ichoice,id3x,  11d17s23
     $     isymmrci,vdinout,nfdat,sr2,ndoub,potdws,ncomp)               6d14s24
      implicit real*8 (a-h,o-z)                                         8d4s23
      dimension darot(*),idora(3,*),iden(*),idoubo(*),irefo(*),noc(*),  8d4s23
     $     nvirt(*),nbasdws(*),idaprt(5,8),nh0av(*),trace(5,8),         8d4s23
     $     nbasisp(*),isblk(4,*),isblk1(4,*),ionex(*),ioooo(*),         8d4s23
     $     ms(28),multh(8,8),ms2(28),mms(8),mmo(8),jmats(*),kmats(*),   8d7s23
     $     i3x(*),isblk4x(5,*),i4x(*),j4o(512),isblkk(4,*),multha(8,8), 8d8s23
     $     nnqn(2,2,2,2),id4o(*),jonex(512),id1x(*),jmatd(512),jmden(*),8d14s23
     $     kmatd(512),kmden(*),itt(*),tder(8),ichoice(*),id3x(*),       11d14s23
     $     i3xd(512),i4xd(512,5),vdinout(*),nfdat(5,4,*),dcore2(5),     6d26s24
     $     dcore2k(5),eoad(8),deoa(5),de0x(5),de1x(5),dejj(5),dekk(5),  6d26s24
     $     iuseage(8,8),xjt(512),xjt2(512),dcore1(5),daa1(5),de3x(5),   5d10s24
     $     de4x(5),ie4v(8,5),sum4vt(5),sum4tv(5),iarot(8)               6d25s24
      data ms/1,5,1,5,1,2,3,6,7,5,5,1,2,3,6,7,1,5,2,3,1,4,5,8,1,6,7,5/  8d1s23
      data multh/1,2,3,4,5,6,7,8,2,1,4,3,6,5,8,7,3,4,1,2,7,8,5,6,
     $     4,3,2,1,8,7,6,5,5,6,7,8,1,2,3,4,6,5,8,7,2,1,4,3,
     $     7,8,5,6,3,4,1,2,8,7,6,5,4,3,2,1/
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      write(6,*)('Hi, my name is playr ')
      write(6,*)('isymmrci: '),isymmrci,loc(isymmrci)
      write(6,*)('ichoice: '),(ichoice(i),i=1,11)
      srh=1d0/sr2                                                       5d10s24
      do i=1,8
       do j=1,8
        iuseage(j,i)=0
       end do
      end do
      do i=1,512
       xjt(i)=0d0
       xjt2(i)=0d0
      end do
      ibcoffo=ibcoff
      npart=3                                                           10d19s23
      if(ipass.ge.0)npart=4                                             10d19s23
      if(nsymb.ne.1)then
       write(6,*)('warning!!! only nosym case has been debugged!')
       write(6,*)('J symmetries: ')
       do is=1,nsdlk
        write(6,177)(isblk(j,is),j=1,4)
  177   format(1x,4i1)
       end do
      else
       write(6,*)('setting up symmetry labels ...')
       do i=1,8
        mms(i)=0
       end do
       do i=1,28
        mms(ms(i))=mms(ms(i))+1
        ms2(i)=mms(ms(i))                                               8d7s23
        write(6,17)i,ms2(i),ms(i)
   17   format(2i3,'s',i1)                                               8d7s23
        if(i.eq.noc(1))then
         do j=1,8
          mmo(j)=mms(j)
         end do
        end if
       end do                                                           8d7s23
      end if
      ibcoffo=ibcoff                                                    8d4s23
      if(nsymb.eq.1)then
       igoal=2
       igoxl=3
       igoul=3
      else
       igoal=1
       igoxl=2
       igoul=2
      end if
      igoal=1
      igoxl=2
      igozl=2
      do isb=1,nsymb                                                    8d4s23
       tder(isb)=0d0                                                    10d18s23
       do ikind=1,npart                                                 10d19s23
        trace(ikind,isb)=0d0                                            8d4s23
       end do                                                           8d4s23
      end do                                                            8d4s23
      npartp=npart+1                                                    6d26s24
      do is=1,nsdlk                                                     8d7s23
       j4o(is)=ibcoff                                                   8d7s23
       if(isblk(1,is).eq.isblk(2,is))then                               8d7s23
        nrow=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                  8d7s23
        ncol=(noc(isblk(3,is))*(noc(isblk(3,is))+1))/2                  8d7s23
       else                                                             8d7s23
        nrow=noc(isblk(1,is))*noc(isblk(2,is))                          8d7s23
        ncol=noc(isblk(3,is))*noc(isblk(4,is))                          8d7s23
       end if                                                           8d7s23
       ibcoff=ibcoff+nrow*ncol*npartp                                   6d26s24
       ncol=nvirt(isblk(3,is))*nvirt(isblk(4,is))                       8d11s23
       jmatd(is)=ibcoff                                                 8d11s23
       ibcoff=ibcoff+nrow*ncol*npartp                                   6d26s24
      end do                                                            8d7s23
      do is=1,nsdlk1                                                    8d10s23
       jonex(is)=ibcoff                                                 8d10s23
       if(isblk1(1,is).eq.isblk1(2,is))then                             8d10s23
        nrow=(noc(isblk1(1,is))*(noc(isblk1(1,is))+1))/2                8d10s23
        nrow3=(nvirt(isblk1(1,is))*(nvirt(isblk1(1,is))+1))/2           11d14s23
       else                                                             8d10s23
        nrow=noc(isblk1(1,is))*noc(isblk1(2,is))                        8d10s23
        nrow3=nvirt(isblk1(1,is))*nvirt(isblk1(2,is))                   11d14s23
       end if                                                           8d10s23
       ncol=noc(isblk1(3,is))*nvirt(isblk1(4,is))                       8d10s23
       ibcoff=ibcoff+nrow*ncol*npartp                                   6d26s24
       ncol=irefo(isblk1(3,is))*nvirt(isblk1(4,is))                     11d14s23
       i3xd(is)=ibcoff                                                  11d14s23
       ibcoff=i3xd(is)+nrow3*ncol*npartp                                6d26s24
      end do                                                            8d10s23
      do is=1,nsdlkk
       kmatd(is)=ibcoff                                                 8d14s23
       nrow=noc(isblkk(1,is))*noc(isblkk(2,is))                         10d23s23
       ncol=nvirt(isblkk(3,is))*nvirt(isblkk(4,is))                     8d14s23
       ibcoff=ibcoff+nrow*ncol*npartp                                   6d26s24
      end do                                                            8d14s23
      do is=1,n4x                                                       11d17s23
       i4xd(is,1)=ibcoff                                                11d17s23
       if(isblk4x(1,is).eq.isblk4x(2,is))then                           11d17s23
        nab=(nvirt(isblk4x(1,is))*(nvirt(isblk4x(1,is))+1))/2           11d17s23
        ncd=(nvirt(isblk4x(3,is))*(nvirt(isblk4x(3,is))+1))/2           11d17s23
       else                                                             11d17s23
        nab=nvirt(isblk4x(1,is))*nvirt(isblk4x(2,is))                   11d17s23
        ncd=nvirt(isblk4x(3,is))*nvirt(isblk4x(4,is))                   11d17s23
       end if                                                           11d17s23
       itriab=((isblk4x(1,is)*(isblk4x(1,is)-1))/2)+isblk4x(2,is)       11d17s23
       itricd=((isblk4x(3,is)*(isblk4x(3,is)-1))/2)+isblk4x(4,is)       11d17s23
       if(itriab.eq.itricd)then                                         11d17s23
        nall=(nab*(nab+1))/2                                            11d17s23
       else                                                             11d17s23
        nall=nab*ncd                                                    11d17s23
       end if                                                           11d17s23
       ibcoff=i4xd(is,1)+nall                                           11d17s23
       write(6,*)('i4xd '),is,(isblk4x(j,is),j=1,4),nab,ncd,itriab,
     $      itricd,nall
       do j=2,npartp                                                    6d26s24
        i4xd(is,j)=ibcoff                                               11d17s23
        ibcoff=i4xd(is,j)+nall                                          11d17s23
       end do                                                           11d17s23
      end do                                                            11d17s23
      call enough('playr.j4o',bc,ibc)                                   8d7s23
      do iz=j4o(1),ibcoff-1                                             8d7s23
       bc(iz)=0d0                                                       8d7s23
      end do                                                            8d7s23
      ioff=1                                                            8d4s23
      jh0copy=ih0copy                                                   8d4s23
      tcannon=0d0                                                       6d26s24
      ax=0d0
      tx=0d0
      tcut=2d0
      ttop=4d0
      atop=0.017d0
c     acut=0.05 ok
c     acut=0.02 ok
c     acut=0.01 bad
      acut=0.016d0
      ass=1d0                                                           10d24s23
      ecore1=0d0                                                         5d2s24
      ecore2=0d0                                                         5d2s24
      ecore2k=0d0                                                       5d2s24
      e0x=0d0                                                           5d2s24
      e1x=0d0                                                           5d2s24
      e3x=0d0                                                           5d10s24
      e4x=0d0                                                           5d10s24
      eoa=0d0
      eaa1=0d0
      ejj=0d0
      ekk=0d0
      do kk=1,nsymb
       eoad(kk)=0d0
       do l=1,4
        ie4v(kk,l)=ibcoff                                                  5d15s24
        ibcoff=ie4v(kk,l)+nvirt(kk)*nbasdws(kk)                            5d15s24
       end do
      end do
      call enough('playr.e4v',bc,ibc)                                   5d15s24
      do iz=ie4v(1,1),ibcoff-1                                            5d15s24
       bc(iz)=0d0                                                       5d15s24
      end do                                                            5d15s24
      do kk=1,4
       daa1(kk)=0d0                                                     5d10s24
       dcore1(kk)=0d0                                                   5d9s24
       dcore2(kk)=0d0
       dcore2k(kk)=0d0
       deoa(kk)=0d0                                                     5d2s24
       de0x(kk)=0d0
       de1x(kk)=0d0                                                     5d2s24
       dejj(kk)=0d0
       dekk(kk)=0d0
       de3x(kk)=0d0
       de4x(kk)=0d0                                                     5d10s24
      end do
      do isb=1,nsymb                                                    8d4s23
       write(6,*)('for symmetry block '),isb                            8d4s23
       write(6,*)('darot is '),ioff                                          8d4s23
       iarot(isb)=ioff                                                  6d25s24
       call prntm2(darot(ioff),nbasdws(isb),nbasdws(isb),nbasdws(isb))  8d4s23
       write(6,*)('ax so far ...'),ax
       write(6,*)('integral scale factor to use: '),ass
       nn=nbasdws(isb)**2                                               8d4s23
       if(ipass.ge.0)then
        write(6,*)('transder is '),ioff+nn                                      10d18s23
        call prntm2(darot(ioff+nn),nbasdws(isb),nbasdws(isb),           10d18s23
     $       nbasdws(isb))                                              10d18s23
        write(6,*)('tx so far ... '),tx
        idaprt(4,isb)=ibcoff                                            10d19s23
        ibcoff=idaprt(4,isb)+nn                                         10d19s23
        call enough('playr.4',bc,ibc)                                   10d19s23
        do iz=0,nn-1                                                    10d19s23
         bc(idaprt(4,isb)+iz)=darot(ioff+nn+iz)                         10d19s23
        end do
       else
        idaprt(4,isb)=ibcoff                                            10d19s23
        ibcoff=idaprt(4,isb)+nn                                         10d19s23
        call enough('playr.4',bc,ibc)                                   10d19s23
        do iz=0,nn-1                                                    10d19s23
         bc(idaprt(4,isb)+iz)=0d0                                       6d24s24
        end do
       end if
       idaprt(1,isb)=ibcoff                                             8d4s23
       idaprt(2,isb)=idaprt(1,isb)+nn                                   8d4s23
       idaprt(3,isb)=idaprt(2,isb)+nn                                   8d4s23
       idaprt(5,isb)=idaprt(3,isb)+nn                                          8d4s23
       ibcoff=idaprt(5,isb)+nn                                          6d24s24
       call enough('playa.daprt',bc,ibc)                                8d4s23
       do iz=idaprt(1,isb),ibcoff-1                                     8d4s23
        bc(iz)=0d0                                                      8d4s23
       end do                                                           8d4s23
       do i=0,nn-1                                                      6d24s24
        bc(idaprt(5,isb)+i)=darot(ioff+i)                               6d24s24
       end do                                                           6d24s24
       do ia=0,idora(2,isb)-1                                           8d4s23
        iap=ia+idora(1,isb)                                             8d4s23
        do id=0,idora(1,isb)-1                                          8d4s23
         ida=id+iap*nbasdws(isb)                                         8d4s23
         iad=iap+id*nbasdws(isb)                                         8d4s23
         bc(idaprt(1,isb)+ida)=darot(ioff+ida)                          8d4s23
         bc(idaprt(1,isb)+iad)=darot(ioff+iad)                          8d4s23
        end do                                                          8d4s23
       end do                                                           8d4s23
       do iv=0,nvirt(isb)-1                                             8d4s23
        ivp=iv+noc(isb)                                                 8d4s23
        do id=0,idora(1,isb)-1                                          8d4s23
         idv=id+nbasdws(isb)*ivp
         ivd=ivp+nbasdws(isb)*id
         bc(idaprt(2,isb)+idv)=darot(ioff+idv)                          8d4s23
         bc(idaprt(2,isb)+ivd)=darot(ioff+ivd)                          8d4s23
        end do                                                          8d4s23
        do ia=0,idora(2,isb)-1                                          8d4s23
         iap=ia+idora(1,isb)                                            8d4s23
         iav=iap+nbasdws(isb)*ivp
         iva=ivp+nbasdws(isb)*iap
         bc(idaprt(3,isb)+iav)=darot(ioff+iav)                          8d4s23
         bc(idaprt(3,isb)+iva)=darot(ioff+iva)                          8d4s23
        end do                                                          8d4s23
       end do                                                           8d4s23
       write(6,*)('da part: ')
       call prntm2(bc(idaprt(1,isb)),nbasdws(isb),nbasdws(isb),
     $      nbasdws(isb))
       write(6,*)('dv part: ')
       call prntm2(bc(idaprt(2,isb)),nbasdws(isb),nbasdws(isb),
     $      nbasdws(isb))
       write(6,*)('av part: ')
       call prntm2(bc(idaprt(3,isb)),nbasdws(isb),nbasdws(isb),
     $      nbasdws(isb))
       itmp=ibcoff                                                      8d4s23
       ibcoff=itmp+nn                                                   8d4s23
       ecore1p=0d0                                                      5d10s24
       do i=0,idoubo(isb)-1                                             5d2s24
        ii=jh0copy+i*(nbasdws(isb)+1)                                   5d2s24
        ecore1=ecore1+2d0*bc(ii)
        ecore1p=ecore1p+2d0*bc(ii)                                      5d10s24
       end do
       write(6,*)('ecore1p '),ecore1p,ecore1
       do i=0,nh0av(isb)-1                                              5d2s24
        ip=i+idoubo(isb)                                                5d2s24
        do j=0,nh0av(isb)-1                                             5d2s24
         jp=j+idoubo(isb)                                               5d2s24
         iad1=jh0copy+jp+nbasdws(isb)*ip                                5d2s24
         iad2=iden(isb)+j+nh0av(isb)*i                                  5d2s24
         eaa1=eaa1+bc(iad1)*bc(iad2)                                    5d2s24
         eoad(isb)=eoad(isb)+bc(iad1)*bc(iad2)                          5d2s24
        end do                                                          5d2s24
       end do                                                           5d2s24
       if(nbasdws(isb).gt.0)then                                        7d11s24
        call dgemm('n','n',nbasdws(isb),nbasdws(isb),nbasdws(isb),       6d26s24
     $       1d0,bc(jh0copy),nbasdws(isb),darot(iarot(isb)),            6d26s24
     $       nbasdws(isb),0d0,bc(itmp),nbasdws(isb),'playr.h0')         8d4s23
        do i=0,nbasdws(isb)-1                                            6d26s24
         do j=0,i                                                        6d26s24
          ji=itmp+j+nbasdws(isb)*i                                       6d26s24
          ij=itmp+i+nbasdws(isb)*j                                       6d26s24
          sum=bc(ji)+bc(ij)                                              6d26s24
          bc(ji)=sum                                                     6d26s24
          bc(ij)=sum                                                     6d26s24
         end do                                                          6d26s24
        end do                                                           6d26s24
       end if                                                           7d11s24
       if(ichoice(1).ne.0)then                                          6d26s24
        do i=0,idoubo(isb)-1                                            8d4s23
         ii=itmp+i*(nbasdws(isb)+1)                                     8d4s23
         orig=tcannon
         tcannon=tcannon+2d0*bc(ii)                                     6d26s24
         if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannona '),orig,
     $        bc(ii),tcannon
        end do                                                          8d4s23
       end if                                                           6d26s24
       if(ichoice(2).ne.0)then                                          6d26s24
        do i=0,nh0av(isb)-1                                             6d26s24
         ip=i+idoubo(isb)                                               6d26s24
         iadd=iden(isb)+nh0av(isb)*i                                    6d26s24
         iadh=itmp+idoubo(isb)+nbasdws(isb)*ip                          6d26s24
         do j=0,nh0av(isb)-1                                            6d26s24
         orig=tcannon
          tcannon=tcannon+bc(iadd+j)*bc(iadh+j)                         6d26s24
         if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannonb '),orig,
     $        bc(iadd+j),bc(iadh+j),tcannon
         end do                                                         6d26s24
        end do                                                          6d26s24
       end if                                                           6d26s24
       do ikind=1,3                                                     8d4s23
        if(nbasdws(isb).gt.0)then                                       7d11s24
         call dgemm('n','n',nbasdws(isb),nbasdws(isb),nbasdws(isb),      8d4s23
     $       1d0,bc(jh0copy),nbasdws(isb),bc(idaprt(ikind,isb)),        8d4s23
     $       nbasdws(isb),0d0,bc(itmp),nbasdws(isb),'playr.h0')         8d4s23
         write(6,*)('for kind '),ikind
         do i=0,nbasdws(isb)-1                                           8d4s23
          do j=0,i                                                       8d4s23
           ji=itmp+j+nbasdws(isb)*i                                      8d4s23
           ij=itmp+i+nbasdws(isb)*j                                      8d4s23
           sum=bc(ji)+bc(ij)                                             8d4s23
           bc(ji)=sum                                                    8d4s23
           bc(ij)=sum                                                    8d4s23
          end do                                                         8d4s23
         end do                                                          8d4s23
        end if                                                          7d11s24
        do i=0,idoubo(isb)-1                                            5d10s24
         ii=itmp+i*(nbasdws(isb)+1)                                     5d10s24
         dcore1(ikind)=dcore1(ikind)+2d0*bc(ii)                         5d10s24
        end do                                                          5d10s24
        do i=0,nh0av(isb)-1                                              5d2s24
         ip=i+idoubo(isb)                                                5d2s24
         do j=0,nh0av(isb)-1                                             5d2s24
          jp=j+idoubo(isb)                                               5d2s24
          iad1=itmp+jp+nbasdws(isb)*ip                                  5d10s24
          iad2=iden(isb)+j+nh0av(isb)*i                                 5d10s24
          daa1(ikind)=daa1(ikind)+bc(iad1)*bc(iad2)                     5d10s24
         end do                                                          5d2s24
        end do                                                           5d2s24
        if(ichoice(1).ne.0)then
         do i=0,idoubo(isb)-1                                            8d4s23
          ii=itmp+i*(nbasdws(isb)+1)                                     8d4s23
          trace(ikind,isb)=trace(ikind,isb)+2d0*bc(ii)                   8d4s23
         end do                                                          8d4s23
         write(6,*)('trace after doub '),trace(ikind,isb)                8d4s23
         if(ikind.eq.3)then
          xx=trace(1,isb)+trace(2,isb)+trace(3,isb)
          write(6,*)('summed over kinds '),xx
         end if
        end if                                                          10d18s23
        if(ichoice(2).ne.0)then                                         10d18s23
         do i=0,nh0av(isb)-1                                             8d4s23
          ip=i+idoubo(isb)                                               8d4s23
          iadd=iden(isb)+nh0av(isb)*i                                    8d4s23
          iadh=itmp+idoubo(isb)+nbasdws(isb)*ip                          8d4s23
          do j=0,nh0av(isb)-1                                            8d4s23
           trace(ikind,isb)=trace(ikind,isb)+bc(iadd+j)*bc(iadh+j)       8d4s23
          end do                                                         8d4s23
         end do                                                          8d4s23
         write(6,*)('trace after rest '),trace(ikind,isb)                8d4s23
        end if
       end do                                                           8d4s23
       ioff=ioff+nn                                                     8d4s23
       if(ipass.ge.0)then                                               10d18s23
        if(nbasdws(isb).gt.0)then                                       7d11s24
         call dgemm('n','n',nbasdws(isb),nbasdws(isb),nbasdws(isb),      8d4s23
     $       1d0,bc(jh0copy),nbasdws(isb),darot(ioff),                  10d18s23
     $       nbasdws(isb),0d0,bc(itmp),nbasdws(isb),'playr.h0')         8d4s23
         write(6,*)('for transder ...')
         do i=0,nbasdws(isb)-1                                           8d4s23
          do j=0,i                                                       8d4s23
           ji=itmp+j+nbasdws(isb)*i                                      8d4s23
           ij=itmp+i+nbasdws(isb)*j                                      8d4s23
           sum=bc(ji)+bc(ij)                                             8d4s23
           bc(ji)=sum                                                    8d4s23
           bc(ij)=sum                                                    8d4s23
          end do                                                         8d4s23
         end do                                                          8d4s23
        end if                                                          7d11s24
        do i=0,idoubo(isb)-1                                            5d10s24
         ii=itmp+i*(nbasdws(isb)+1)                                     5d10s24
         dcore1(4)=dcore1(4)+2d0*bc(ii)                                 5d10s24
        end do                                                          5d10s24
        do i=0,nh0av(isb)-1                                             5d10s24
         ip=i+idoubo(isb)                                               5d10s24
         do j=0,nh0av(isb)-1                                            5d10s24
          jp=j+idoubo(isb)                                              5d10s24
          iad1=itmp+jp+nbasdws(isb)*ip                                  5d10s24
          iad2=iden(isb)+j+nh0av(isb)*i                                 5d10s24
          daa1(4)=daa1(4)+bc(iad1)*bc(iad2)                             5d10s24
         end do                                                         5d10s24
        end do                                                          5d10s24
        if(ichoice(1).ne.0)then                                         10d18s23
         do i=0,idoubo(isb)-1                                            8d4s23
          ii=itmp+i*(nbasdws(isb)+1)                                     8d4s23
          trace(4,isb)=trace(4,isb)+2d0*bc(ii)                                 10d18s23
         end do                                                          8d4s23
         write(6,*)('trace after doub '),trace(4,isb)                       10d18s23
        end if
        if(ichoice(2).ne.0)then                                         10d18s23
         do i=0,nh0av(isb)-1                                             8d4s23
          ip=i+idoubo(isb)                                               8d4s23
          iadd=iden(isb)+nh0av(isb)*i                                    8d4s23
          iadh=itmp+idoubo(isb)+nbasdws(isb)*ip                          8d4s23
          do j=0,nh0av(isb)-1                                            8d4s23
           trace(4,isb)=trace(4,isb)+bc(iadd+j)*bc(iadh+j)                     10d18s23
          end do                                                         8d4s23
         end do                                                          8d4s23
         write(6,*)('trace after rest '),trace(4,isb)                       10d18s23
        end if                                                          10d18s23
        ioff=ioff+nn                                                    10d18s23
       end if                                                           10d18s23
       jh0copy=jh0copy+(ncomp*nbasisp(isb))**2                          6d14s24
      end do                                                            8d4s23
      write(6,*)('potdws: '),potdws
      write(6,*)('ecore1: '),ecore1,ecore1+potdws
      if(nsymb.gt.1)then
       write(6,*)('summed over symmetry blocks: ')
       gtot=0d0
       do ikind=1,npart                                                 10d19s23
        write(6,*)('for kind '),ikind
        sum=trace(ikind,1)
        do i=2,nsymb
         sum=sum+trace(ikind,i)
        end do
        write(6,*)sum,gtot                                              10d19s23
        gtot=gtot+sum                                                   10d19s23
       end do
       write(6,*)('grand total: '),gtot
      end if
      val=0d0
      do isd=1,nsymb
       do isc=1,isd
        iscd=multha(isc,isd)                                            8d8s23
        icd=((isd*(isd-1))/2)+isc
        ncd=min(nbasdws(isc),nbasdws(isd))
        do isb=1,nsymb
         isbcd=multh(isb,iscd)                                          8d8s23
         do isa=1,isb
          iab=((isb*(isb-1))/2)+isa
          nab=min(nbasdws(isb),nbasdws(isa))
          if(icd.ge.iab.and.min(nab,ncd).gt.0.and.                      8d8s23
     $         multha(isa,isbcd).eq.1)then                              8d8s23
           nrow=nbasdws(isa)*nbasdws(isb)
           ncol=nbasdws(isc)*nbasdws(isd)
           itmp=ibcoff
           ibcoff=itmp+nrow*ncol
           call enough('playr.tmp4',bc,ibc)
           do iz=itmp,ibcoff-1
            bc(iz)=val
           end do
           sz4=0d0
           do i2=0,nvirt(isd)-1                                         1d11s24
            i2p=i2+noc(isd)                                             1d11s24
            do i1=0,nvirt(isc)-1                                        1d11s24
             i1p=i1+noc(isc)                                            1d11s24
             do i4=0,nvirt(isb)-1                                       1d11s24
              i4p=i4+noc(isb)                                           1d11s24
              do i3=0,nvirt(isa)-1                                      1d11s24
               i3p=i3+noc(isa)                                          1d11s24
               iad=itmp+i3p+nbasdws(isa)*(i4p+nbasdws(isb)*(i1p         1d11s24
     $                 +nbasdws(isc)*i2p))                              8d7s23
               if(ichoice(11).ne.0)then
                bc(iad)=get4int(i4x,isa,isb,isc,isd,                     1d11s24
     $               i3,i4,i1,i2,nvirt,idum,bc,ibc)                     1d11s24
                if(iad.eq.igozl)write(6,*)('gozl4x '),ibc(ibcoff)
                sz4=sz4+bc(iad)**2
                if(ibc(ibcoff).lt.0)then
                 write(6,*)('get4int returns zero!'),isa,isb,isc,isd,
     $               i1,i2,i3,i4
                 write(6,*)(ibc(ibcoff+iqy),iqy=0,8)
                 stop 'playr'
                end if
               else
                bc(iad)=0d0                                             6d14s24
               end if                                                   6d14s24
              end do                                                    1d11s24
             end do                                                     1d11s24
            end do                                                      1d11s24
           end do                                                       1d11s24
           sz4=sqrt(sz4/
     $          dfloat(nvirt(isa)*nvirt(isb)*nvirt(isc)*nvirt(isd)))
           write(6,*)('sz4 = '),sz4
           ihit=0                                                       8d7s23
           do is=1,nsdlk1                                               8d7s23
            if(isblk1(1,is).eq.isblk1(2,is))then                        8d7s23
             nrow1=(noc(isblk1(1,is))*(noc(isblk1(1,is))+1))/2          8d7s23
             nrow3=(nvirt(isblk1(1,is))*(nvirt(isblk1(1,is))+1))/2      8d7s23
             isw=0                                                      8d7s23
            else                                                        8d7s23
             nrow1=noc(isblk1(1,is))*noc(isblk1(2,is))                  8d7s23
             nrow3=nvirt(isblk1(1,is))*nvirt(isblk1(2,is))              8d7s23
             isw=1                                                      8d7s23
            end if                                                      8d7s23
            ncol3=noc(isblk1(3,is))*nvirt(isblk1(4,is))
            call ilimts(noc(isblk1(3,is)),nvirt(isblk1(4,is)),mynprocg, 8d7s23
     $           mynowprog,il,ih,i1s,i1e,i2s,i2e)                       8d7s23
            nhere=ih+1-il                                               8d7s23
            if(nhere.gt.0)then                                          8d7s23
             if(isblk1(1,is).eq.isa.and.isblk1(2,is).eq.isb.and.         8d7s23
     $         isblk1(3,is).eq.isd)then                                 8d7s23
              i10=i1s                                                   8d7s23
              i1n=noc(isblk1(3,is))                                     8d7s23
              iii=i3x(is)                                               8d7s23
              ii=ionex(is)                                              8d8s23
              do i2=i2s,i2e                                             8d7s23
               i2p=i2+noc(isblk1(4,is))-1                               8d7s23
               if(i2.eq.i2e)i1n=i1e                                     8d7s23
               do i1=i10,i1n                                            8d7s23
                i1m=i1-1                                                8d7s23
                do i4=0,nvirt(isblk1(2,is))-1                           8d7s23
                 i4p=i4+noc(isblk1(2,is))                               8d7s23
                 do i3=0,nvirt(isblk1(1,is))-1                          8d7s23
                  i3p=i3+noc(isblk1(1,is))                              8d7s23
                  irec=i3+nvirt(isblk1(1,is))*i4                        8d7s23
                  ix=max(i3,i4)                                         8d7s23
                  in=min(i3,i4)                                         8d7s23
                  itri=((ix*(ix+1))/2)+in                               8d7s23
                  itri=itri+isw*(irec-itri)+iii                         8d7s23
                  iad=itmp+i3p+nbasdws(isa)*(i4p+nbasdws(isb)*(i2p      8d7s23
     $                 +nbasdws(isc)*i1m))                              8d7s23
                  bc(iad)=bc(itri)                                      8d7s23
                  if(iad.eq.igozl)write(6,*)('gozlz '),bc(igozl)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                iii=iii+nrow3                                           8d7s23
                do i4=0,noc(isblk1(2,is))-1                             8d7s23
                 do i3=0,noc(isblk1(1,is))-1                            8d7s23
                  irec=i3+noc(isblk1(1,is))*i4                          8d7s23
                  ix=max(i3,i4)                                         8d7s23
                  in=min(i3,i4)                                         8d7s23
                  itri=((ix*(ix+1))/2)+in                               8d7s23
                  itri=itri+isw*(irec-itri)+ii                          8d7s23
                  iad=itmp+i3+nbasdws(isa)*(i4+nbasdws(isb)*(i2p        8d9s23
     $                 +nbasdws(isc)*i1m))                              8d7s23
                  bc(iad)=bc(itri)                                      8d7s23
                  if(iad.eq.igozl)write(6,*)('gozly '),bc(igozl),
     $                 ('onex of '),is,itri-ionex(is),
     $                 (isblk1(j,is),j=1,4)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                ii=ii+nrow1                                             8d7s23
               end do                                                   8d7s23
               i10=1                                                    8d7s23
              end do                                                    8d7s23
              ihit=ihit+1                                               8d7s23
             else if(isblk1(1,is).eq.isb.and.isblk1(2,is).eq.isa.and.    8d7s23
     $         isblk1(3,is).eq.isd)then                                 8d7s23
              i10=i1s                                                   8d7s23
              i1n=noc(isblk1(3,is))                                     8d7s23
              iii=i3x(is)                                               8d7s23
              ii=ionex(is)                                              8d8s23
              do i2=i2s,i2e                                             8d7s23
               i2p=i2+noc(isblk1(4,is))-1                               8d7s23
               if(i2.eq.i2e)i1n=i1e                                     8d7s23
               do i1=i10,i1n                                            8d7s23
                i1m=i1-1                                                8d7s23
                do i4=0,nvirt(isblk1(2,is))-1                           8d7s23
                 i4p=i4+noc(isblk1(2,is))                               8d7s23
                 do i3=0,nvirt(isblk1(1,is))-1                          8d7s23
                  i3p=i3+noc(isblk1(1,is))                              8d7s23
                  irec=i3+nvirt(isblk1(1,is))*i4+iii                    8d7s23
                  iad=itmp+i4p+nbasdws(isa)*(i3p+nbasdws(isb)*(i2p      8d7s23
     $                 +nbasdws(isc)*i1m))                              8d7s23
                  bc(iad)=bc(irec)                                      8d14s23
                  if(iad.eq.igozl)write(6,*)('gozlx '),bc(igozl)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                iii=iii+nrow3                                           8d7s23
                do i4=0,noc(isblk1(2,is))-1                             8d7s23
                 do i3=0,noc(isblk1(1,is))-1                            8d7s23
                  irec=i3+noc(isblk1(1,is))*i4+ii                       8d7s23
                  iad=itmp+i4+nbasdws(isa)*(i3+nbasdws(isb)*(i2p        8d7s23
     $                 +nbasdws(isc)*i1m))                              8d7s23
                  bc(iad)=bc(irec)                                      8d8s23
                  if(iad.eq.igozl)write(6,*)('gozlw '),bc(igozl)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                ii=ii+nrow1                                             8d7s23
               end do                                                   8d7s23
               i10=1                                                    8d7s23
              end do                                                    8d7s23
              ihit=ihit+1                                               8d7s23
             end if                                                      8d7s23
             if(isblk1(1,is).eq.isa.and.isblk1(2,is).eq.isb.and.         8d7s23
     $         isblk1(3,is).eq.isc)then                                 8d7s23
              i10=i1s                                                   8d7s23
              i1n=noc(isblk1(3,is))                                     8d7s23
              iii=i3x(is)                                               8d7s23
              ii=ionex(is)                                              8d8s23
              do i2=i2s,i2e                                             8d7s23
               i2p=i2+noc(isblk1(4,is))-1                               8d7s23
               if(i2.eq.i2e)i1n=i1e                                     8d7s23
               do i1=i10,i1n                                            8d7s23
                i1m=i1-1                                                8d7s23
                do i4=0,nvirt(isblk1(2,is))-1                           8d7s23
                 i4p=i4+noc(isblk1(2,is))                               8d7s23
                 do i3=0,nvirt(isblk1(1,is))-1                          8d7s23
                  i3p=i3+noc(isblk1(1,is))                              8d7s23
                  irec=i3+nvirt(isblk1(1,is))*i4                        8d7s23
                  ix=max(i3,i4)                                         8d7s23
                  in=min(i3,i4)                                         8d7s23
                  itri=((ix*(ix+1))/2)+in                               8d7s23
                  itri=itri+isw*(irec-itri)+iii                         8d7s23
                  iad=itmp+i3p+nbasdws(isa)*(i4p+nbasdws(isb)*(i1m      8d7s23
     $                 +nbasdws(isc)*i2p))                              8d7s23
                  bc(iad)=bc(itri)                                      8d7s23
                  if(iad.eq.igozl)write(6,*)('gozlv '),bc(igozl)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                iii=iii+nrow3                                           8d7s23
                do i4=0,noc(isblk1(2,is))-1                             8d7s23
                 do i3=0,noc(isblk1(1,is))-1                            8d7s23
                  irec=i3+noc(isblk1(1,is))*i4                          8d7s23
                  ix=max(i3,i4)                                         8d7s23
                  in=min(i3,i4)                                         8d7s23
                  itri=((ix*(ix+1))/2)+in                               8d7s23
                  itri=itri+isw*(irec-itri)+ii                          8d7s23
                  iad=itmp+i3+nbasdws(isa)*(i4+nbasdws(isb)*(i1m        8d7s23
     $                 +nbasdws(isc)*i2p))                              8d7s23
                  bc(iad)=bc(itri)                                      8d7s23
                  if(iad.eq.igozl)write(6,*)('gozlu '),bc(igozl)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                ii=ii+nrow1                                             8d7s23
               end do                                                   8d7s23
               i10=1                                                    8d7s23
              end do                                                    8d7s23
              ihit=ihit+1                                               8d7s23
             else if(isblk1(1,is).eq.isb.and.isblk1(2,is).eq.isa.and.    8d7s23
     $         isblk1(3,is).eq.isc)then                                 8d7s23
              i10=i1s                                                   8d7s23
              i1n=noc(isblk1(3,is))                                     8d7s23
              iii=i3x(is)                                               8d7s23
              ii=ionex(is)                                              8d8s23
              do i2=i2s,i2e                                             8d7s23
               i2p=i2+noc(isblk1(4,is))-1                               8d7s23
               if(i2.eq.i2e)i1n=i1e                                     8d7s23
               do i1=i10,i1n                                            8d7s23
                i1m=i1-1                                                8d7s23
                do i4=0,nvirt(isblk1(2,is))-1                           8d7s23
                 i4p=i4+noc(isblk1(2,is))                               8d7s23
                 do i3=0,nvirt(isblk1(1,is))-1                          8d7s23
                  i3p=i3+noc(isblk1(1,is))                              8d7s23
                  irec=i3+nvirt(isblk1(1,is))*i4+iii                    8d7s23
                  iad=itmp+i4p+nbasdws(isa)*(i3p+nbasdws(isb)*(i1m      8d7s23
     $                 +nbasdws(isc)*i2p))                              8d7s23
                  bc(iad)=bc(irec)                                      8d14s23
                  if(iad.eq.igozl)write(6,*)('gozlt '),bc(igozl),
     $                 (isblk1(j,is),j=1,4),isa,isb,isc,isd
                 end do                                                 8d7s23
                end do                                                  8d7s23
                iii=iii+nrow3                                           8d7s23
                do i4=0,noc(isblk1(2,is))-1                             8d7s23
                 do i3=0,noc(isblk1(1,is))-1                            8d7s23
                  irec=i3+noc(isblk1(1,is))*i4+ii                       8d7s23
                  iad=itmp+i4+nbasdws(isa)*(i3+nbasdws(isb)*(i1m        8d7s23
     $                 +nbasdws(isc)*i2p))                              8d7s23
                  bc(iad)=bc(irec)                                      8d7s23
                  if(iad.eq.igozl)write(6,*)('gozls '),bc(igozl)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                ii=ii+nrow1                                             8d7s23
               end do                                                   8d7s23
               i10=1                                                    8d7s23
              end do                                                    8d7s23
              ihit=ihit+1                                               8d7s23
             end if                                                      8d7s23
             if(isblk1(3,is).eq.isb)then                                8d7s23
              if(isblk1(1,is).eq.isc.and.isblk1(2,is).eq.isd)then       8d7s23
               i10=i1s                                                  8d7s23
               i1n=noc(isblk1(3,is))                                    8d7s23
               iii=i3x(is)                                              8d7s23
               ii=ionex(is)                                              8d8s23
               do i2=i2s,i2e                                            8d7s23
                i2p=i2+noc(isblk1(4,is))-1                              8d7s23
                if(i2.eq.i2e)i1n=i1e                                    8d7s23
                do i1=i10,i1n                                           8d7s23
                 i1m=i1-1                                               8d7s23
                 do i4=0,nvirt(isblk1(2,is))-1                          8d7s23
                  i4p=i4+noc(isblk1(2,is))                              8d7s23
                  do i3=0,nvirt(isblk1(1,is))-1                         8d7s23
                   i3p=i3+noc(isblk1(1,is))                             8d7s23
                   irec=i3+nvirt(isblk1(1,is))*i4                       8d7s23
                   ix=max(i3,i4)                                        8d7s23
                   in=min(i3,i4)                                        8d7s23
                   itri=((ix*(ix+1))/2)+in                              8d7s23
                   itri=itri+isw*(irec-itri)+iii                        8d7s23
                   iad=itmp+i2p+nbasdws(isa)*(i1m+nbasdws(isb)*(i3p     8d7s23
     $                  +nbasdws(isc)*i4p))                             8d8s23
                   bc(iad)=bc(itri)                                     8d7s23
                  if(iad.eq.igozl)write(6,*)('gozlr '),bc(igozl)
                  end do                                                8d7s23
                 end do                                                 8d7s23
                 iii=iii+nrow3                                          8d7s23
                 do i4=0,noc(isblk1(2,is))-1                            8d7s23
                  do i3=0,noc(isblk1(1,is))-1                           8d7s23
                   irec=i3+noc(isblk1(1,is))*i4                         8d7s23
                   ix=max(i3,i4)                                        8d7s23
                   in=min(i3,i4)                                        8d7s23
                   itri=((ix*(ix+1))/2)+in                              8d7s23
                   itri=itri+isw*(irec-itri)+ii                         8d7s23
                   iad=itmp+i2p+nbasdws(isa)*(i1m+nbasdws(isb)*(i3      8d7s23
     $                  +nbasdws(isc)*i4))                              8d8s23
                   bc(iad)=bc(itri)                                     8d7s23
                  if(iad.eq.igozl)write(6,*)('gozlq '),bc(igozl)
                  end do                                                8d7s23
                 end do                                                 8d7s23
                 ii=ii+nrow1                                            8d7s23
                end do                                                  8d7s23
                i10=1                                                   8d7s23
               end do                                                   8d7s23
               ihit=ihit+1                                               8d7s23
              else if(isblk1(1,is).eq.isd.and.isblk1(2,is).eq.isc)then  8d7s23
               i10=i1s                                                  8d7s23
               i1n=noc(isblk1(3,is))                                    8d7s23
               iii=i3x(is)                                              8d7s23
               ii=ionex(is)                                              8d8s23
               do i2=i2s,i2e                                            8d7s23
                i2p=i2+noc(isblk1(4,is))-1                              8d7s23
                if(i2.eq.i2e)i1n=i1e                                    8d7s23
                do i1=i10,i1n                                           8d7s23
                 i1m=i1-1                                               8d7s23
                 do i4=0,nvirt(isblk1(2,is))-1                          8d7s23
                  i4p=i4+noc(isblk1(2,is))                              8d7s23
                  do i3=0,nvirt(isblk1(1,is))-1                         8d7s23
                   i3p=i3+noc(isblk1(1,is))                             8d7s23
                   irec=iii+i3+nvirt(isblk1(1,is))*i4                       8d7s23
                   iad=itmp+i2p+nbasdws(isa)*(i1m+nbasdws(isb)*(i4p     8d7s23
     $                  +nbasdws(isc)*i3p))                             8d8s23
                   bc(iad)=bc(irec)                                     8d7s23
                  if(iad.eq.igozl)write(6,*)('gozlp '),bc(igozl)
                  end do                                                8d7s23
                 end do                                                 8d7s23
                 iii=iii+nrow3                                          8d7s23
                 do i4=0,noc(isblk1(2,is))-1                            8d7s23
                  do i3=0,noc(isblk1(1,is))-1                           8d7s23
                   irec=ii+i3+noc(isblk1(1,is))*i4                       8d7s23
                   iad=itmp+i2p+nbasdws(isa)*(i1m+nbasdws(isb)*(i4      8d7s23
     $                  +nbasdws(isc)*i3))                              8d8s23
                   bc(iad)=bc(irec)                                     8d7s23
                  if(iad.eq.igozl)write(6,*)('gozlo '),bc(igozl)
                  end do                                                8d7s23
                 end do                                                 8d7s23
                 ii=ii+nrow1                                            8d7s23
                end do                                                  8d7s23
                i10=1                                                   8d7s23
               end do                                                   8d7s23
               ihit=ihit+1                                               8d7s23
              end if                                                    8d7s23
             end if                                                     8d7s23
             if(isblk1(3,is).eq.isa)then                                8d7s23
              if(isblk1(1,is).eq.isc.and.isblk1(2,is).eq.isd)then       8d7s23
               i10=i1s                                                  8d7s23
               i1n=noc(isblk1(3,is))                                    8d7s23
               iii=i3x(is)                                              8d7s23
               ii=ionex(is)                                              8d8s23
               do i2=i2s,i2e                                            8d7s23
                i2p=i2+noc(isblk1(4,is))-1                              8d7s23
                if(i2.eq.i2e)i1n=i1e                                    8d7s23
                do i1=i10,i1n                                           8d7s23
                 i1m=i1-1                                               8d7s23
                 do i4=0,nvirt(isblk1(2,is))-1                          8d7s23
                  i4p=i4+noc(isblk1(2,is))                              8d7s23
                  do i3=0,nvirt(isblk1(1,is))-1                         8d7s23
                   i3p=i3+noc(isblk1(1,is))                             8d7s23
                   irec=i3+nvirt(isblk1(1,is))*i4                       8d7s23
                   ix=max(i3,i4)                                        8d7s23
                   in=min(i3,i4)                                        8d7s23
                   itri=((ix*(ix+1))/2)+in                              8d7s23
                   itri=itri+isw*(irec-itri)+iii                        8d7s23
                   iad=itmp+i1m+nbasdws(isa)*(i2p+nbasdws(isb)*(i3p     8d7s23
     $                  +nbasdws(isc)*i4p))                             8d8s23
                   bc(iad)=bc(itri)                                     8d7s23
                  if(iad.eq.igozl)write(6,*)('gozln '),bc(igozl),itri
                  end do                                                8d7s23
                 end do                                                 8d7s23
                 iii=iii+nrow3                                          8d7s23
                 do i4=0,noc(isblk1(2,is))-1                            8d7s23
                  do i3=0,noc(isblk1(1,is))-1                           8d7s23
                   irec=i3+noc(isblk1(1,is))*i4                         8d7s23
                   ix=max(i3,i4)                                        8d7s23
                   in=min(i3,i4)                                        8d7s23
                   itri=((ix*(ix+1))/2)+in                              8d7s23
                   itri=itri+isw*(irec-itri)+ii                         8d7s23
                   iad=itmp+i1m+nbasdws(isa)*(i2p+nbasdws(isb)*(i3      8d7s23
     $                  +nbasdws(isc)*i4))                              8d8s23
                   bc(iad)=bc(itri)                                     8d7s23
                  if(iad.eq.igozl)write(6,*)('gozlm '),bc(igozl)
                  end do                                                8d7s23
                 end do                                                 8d7s23
                 ii=ii+nrow1                                            8d7s23
                end do                                                  8d7s23
                i10=1                                                   8d7s23
               end do                                                   8d7s23
               ihit=ihit+1                                               8d7s23
              else if(isblk1(1,is).eq.isd.and.isblk1(2,is).eq.isc)then  8d7s23
               i10=i1s                                                  8d7s23
               i1n=noc(isblk1(3,is))                                    8d7s23
               iii=i3x(is)                                              8d7s23
               ii=ionex(is)                                              8d8s23
               do i2=i2s,i2e                                            8d7s23
                i2p=i2+noc(isblk1(4,is))-1                              8d7s23
                if(i2.eq.i2e)i1n=i1e                                    8d7s23
                do i1=i10,i1n                                           8d7s23
                 i1m=i1-1                                               8d7s23
                 do i4=0,nvirt(isblk1(2,is))-1                          8d7s23
                  i4p=i4+noc(isblk1(2,is))                              8d7s23
                  do i3=0,nvirt(isblk1(1,is))-1                         8d7s23
                   i3p=i3+noc(isblk1(1,is))                             8d7s23
                   irec=iii+i3+nvirt(isblk1(1,is))*i4                   8d7s23
                   iad=itmp+i1m+nbasdws(isa)*(i2p+nbasdws(isb)*(i4p     8d7s23
     $                  +nbasdws(isc)*i3p))                             8d8s23
                   bc(iad)=bc(irec)                                     8d7s23
                  if(iad.eq.igozl)write(6,*)('gozll '),bc(igozl),irec
                  end do                                                8d7s23
                 end do                                                 8d7s23
                 iii=iii+nrow3                                          8d7s23
                 do i4=0,noc(isblk1(2,is))-1                            8d7s23
                  do i3=0,noc(isblk1(1,is))-1                           8d7s23
                   irec=ii+i3+noc(isblk1(1,is))*i4                      8d7s23
                   iad=itmp+i1m+nbasdws(isa)*(i2p+nbasdws(isb)*(i4      8d7s23
     $                  +nbasdws(isc)*i3))                              8d8s23
                   bc(iad)=bc(irec)                                     8d7s23
                  if(iad.eq.igozl)write(6,*)('gozlk '),bc(igozl)
                  end do                                                8d7s23
                 end do                                                 8d7s23
                 ii=ii+nrow1                                            8d7s23
                end do                                                  8d7s23
                i10=1                                                   8d7s23
               end do                                                   8d7s23
               ihit=ihit+1                                               8d7s23
              end if                                                    8d7s23
             end if                                                     8d7s23
            end if                                                      8d7s23
           end do                                                       8d7s23
c               14 23
c     Kab^vu = (au|bv)=(ua|bv)=(au|vb)=(ua|vb)
           do is=1,nsdlkk                                               8d7s23
            call ilimts(nvirt(isblkk(3,is)),nvirt(isblkk(4,is)),        8d7s23
     $           mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)              8d7s23
            nhere=ih+1-il                                               8d7s23
            if(nhere.gt.0)then                                          8d7s23
             nrowk=noc(isblkk(1,is))*noc(isblkk(2,is))                  8d7s23
c     (au|bv)
             if(isblkk(1,is).eq.isa.and.isblkk(2,is).eq.isc.and.        8d7s23
     $          isblkk(3,is).eq.isd)then                                8d7s23
              i10=i1s                                                   8d7s23
              i1n=nvirt(isblkk(3,is))                                   8d7s23
              kk=kmats(is)                                              8d7s23
              do i2=i2s,i2e                                             8d7s23
               i2p=i2+noc(isblkk(4,is))-1                               8d7s23
               if(i2.eq.i2e)i1n=i1e                                     8d7s23
               do i1=i10,i1n                                            8d7s23
                i1p=i1+noc(isblkk(3,is))-1                              8d7s23
                do i4=0,noc(isblkk(2,is))-1                             8d7s23
                 do i3=0,noc(isblkk(1,is))-1                            8d7s23
                  irec=kk+i3+noc(isblkk(1,is))*i4                       8d10s23
                  iad=itmp+i3+nbasdws(isa)*(i2p+nbasdws(isb)*(i4        8d7s23
     $                 +nbasdws(isc)*i1p))                              8d7s23
                  bc(iad)=bc(irec)                                      8d7s23
                  if(iad.eq.igozl)write(6,*)('gozlj '),bc(igozl)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                kk=kk+nrowk                                             8d7s23
               end do                                                   8d7s23
               i10=1
              end do                                                    8d7s23
             end if                                                     8d8s23
c     (ua|bv)
             if(isblkk(1,is).eq.isb.and.isblkk(2,is).eq.isc.and.        8d8s23
     $          isblkk(3,is).eq.isd)then                                8d7s23
              i10=i1s                                                   8d7s23
              i1n=nvirt(isblkk(3,is))                                   8d7s23
              kk=kmats(is)                                              8d7s23
              do i2=i2s,i2e                                             8d7s23
               i2p=i2+noc(isblkk(4,is))-1                               8d7s23
               if(i2.eq.i2e)i1n=i1e                                     8d7s23
               do i1=i10,i1n                                            8d7s23
                i1p=i1+noc(isblkk(3,is))-1                              8d7s23
                do i4=0,noc(isblkk(2,is))-1                             8d7s23
                 do i3=0,noc(isblkk(1,is))-1                            8d7s23
                  irec=kk+i3+noc(isblkk(1,is))*i4                       8d10s23
                  iad=itmp+i2p+nbasdws(isa)*(i3+nbasdws(isb)*(i4        8d7s23
     $                 +nbasdws(isc)*i1p))                              8d7s23
                  bc(iad)=bc(irec)                                      8d7s23
                  if(iad.eq.igozl)write(6,*)('gozli '),bc(igozl)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                kk=kk+nrowk                                             8d7s23
               end do                                                   8d7s23
              end do                                                    8d7s23
             end if                                                     8d7s23
c     (au|vb)
             if(isblkk(1,is).eq.isa.and.isblkk(2,is).eq.isd.and.        8d8s23
     $          isblkk(3,is).eq.isc)then                                8d7s23
              i10=i1s                                                   8d7s23
              i1n=nvirt(isblkk(3,is))                                   8d7s23
              kk=kmats(is)                                              8d7s23
              do i2=i2s,i2e                                             8d7s23
               i2p=i2+noc(isblkk(4,is))-1                               8d7s23
               if(i2.eq.i2e)i1n=i1e                                     8d7s23
               do i1=i10,i1n                                            8d7s23
                i1p=i1+noc(isblkk(3,is))-1                              8d7s23
                do i4=0,noc(isblkk(2,is))-1                             8d7s23
                 do i3=0,noc(isblkk(1,is))-1                            8d7s23
                  irec=kk+i3+noc(isblkk(1,is))*i4                       8d10s23
                  iad=itmp+i3+nbasdws(isa)*(i2p+nbasdws(isb)*(i1p       8d8s23
     $                 +nbasdws(isc)*i4))                               8d8s23
                  bc(iad)=bc(irec)                                      8d7s23
                  if(iad.eq.igozl)write(6,*)('gozlh '),bc(igozl)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                kk=kk+nrowk                                             8d7s23
               end do                                                   8d7s23
              end do                                                    8d7s23
             end if                                                     8d7s23
c     (ua|vb)
             if(isblkk(1,is).eq.isb.and.isblkk(2,is).eq.isd.and.        8d8s23
     $          isblkk(3,is).eq.isc)then                                8d7s23
              i10=i1s                                                   8d7s23
              i1n=nvirt(isblkk(3,is))                                   8d7s23
              kk=kmats(is)                                              8d7s23
              do i2=i2s,i2e                                             8d7s23
               i2p=i2+noc(isblkk(4,is))-1                               8d7s23
               if(i2.eq.i2e)i1n=i1e                                     8d7s23
               do i1=i10,i1n                                            8d7s23
                i1p=i1+noc(isblkk(3,is))-1                              8d7s23
                do i4=0,noc(isblkk(2,is))-1                             8d7s23
                 do i3=0,noc(isblkk(1,is))-1                            8d7s23
                  irec=kk+i3+noc(isblkk(1,is))*i4                       8d10s23
                  iad=itmp+i2p+nbasdws(isa)*(i3+nbasdws(isb)*(i1p       8d8s23
     $                 +nbasdws(isc)*i4))                               8d8s23
                  bc(iad)=bc(irec)                                      8d7s23
                  if(iad.eq.igozl)write(6,*)('gozlg '),bc(igozl)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                kk=kk+nrowk                                             8d7s23
               end do                                                   8d7s23
              end do                                                    8d7s23
             end if                                                     8d7s23
            end if                                                      8d7s23
           end do                                                       8d7s23
           do is=1,nsdlk                                                8d7s23
            call ilimts(nvirt(isblk(3,is)),nvirt(isblk(4,is)),mynprocg, 8d7s23
     $           mynowprog,il,ih,i1s,i1e,i2s,i2e)                       8d7s23
            nhere=ih+1-il                                               8d7s23
            if(nhere.gt.0)then                                          8d7s23
             if(isblk(1,is).eq.isblk(2,is))then                         8d7s23
              nrowj=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2            8d7s23
              isw=0                                                     8d7s23
             else                                                       8d7s23
              nrowj=noc(isblk(1,is))*noc(isblk(2,is))                    8d7s23
              isw=1                                                     8d7s23
             end if                                                     8d7s23
             if(isblk(1,is).eq.isa.and.isblk(2,is).eq.isb.and.          8d11s23
     $            isblk(3,is).eq.isc)then                               8d7s23
              jj=jmats(is)                                              8d7s23
              i10=i1s                                                   8d7s23
              i1n=nvirt(isblk(3,is))                                    8d7s23
              do i2=i2s,i2e                                             8d7s23
               i2p=i2+noc(isblk(4,is))-1                                8d7s23
               if(i2.eq.i2e)i1n=i1e                                     8d7s23
               do i1=i10,i1n                                            8d7s23
                i1p=i1+noc(isblk(3,is))-1                               8d7s23
                do i4=0,noc(isblk(2,is))-1                              8d7s23
                 do i3=0,noc(isblk(1,is))-1                             8d7s23
                  irec=i3+noc(isblk(1,is))*i4                           8d7s23
                  ix=max(i3,i4)                                         8d7s23
                  in=min(i3,i4)                                         8d7s23
                  itri=((ix*(ix+1))/2)+in                               8d7s23
                  itri=itri+isw*(irec-itri)+jj                          8d7s23
                  iad=itmp+i3+nbasdws(isa)*(i4+nbasdws(isb)*(i1p        8d7s23
     $                 +nbasdws(isc)*i2p))                              8d7s23
                  bc(iad)=bc(itri)                                      8d7s23
                  if(iad.eq.igozl)write(6,*)('gozlf '),bc(igozl),
     $                 (isblk(j,is),j=1,4),itri-jmats(is)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                jj=jj+nrowj                                              8d7s23
               end do                                                   8d7s23
               i10=1                                                    8d7s23
              end do                                                    8d7s23
             else if(isblk(2,is).eq.isa.and.isblk(1,is).eq.isb.and.     8d7s23
     $            isblk(3,is).eq.isc)then                               8d7s23
              jj=jmats(is)                                              8d7s23
              i10=i1s                                                   8d7s23
              i1n=nvirt(isblk(3,is))                                    8d7s23
              do i2=i2s,i2e                                             8d7s23
               i2p=i2+noc(isblk(4,is))-1                                8d7s23
               if(i2.eq.i2e)i1n=i1e                                     8d7s23
               do i1=i10,i1n                                            8d7s23
                i1p=i1+noc(isblk(3,is))-1                               8d7s23
                do i4=0,noc(isblk(2,is))-1                              8d7s23
                 do i3=0,noc(isblk(1,is))-1                             8d7s23
                  irec=jj+i3+noc(isblk(1,is))*i4                        8d7s23
                  iad=itmp+i4+nbasdws(isa)*(i3+nbasdws(isb)*(i1p        8d7s23
     $                 +nbasdws(isc)*i2p))                              8d7s23
                  bc(iad)=bc(irec)                                      8d7s23
                  if(iad.eq.igozl)write(6,*)('gozle '),bc(igozl)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                jj=jj+nrowj                                              8d7s23
               end do                                                   8d7s23
               i10=1                                                    8d7s23
              end do                                                    8d7s23
             end if                                                     8d7s23
             if(isblk(1,is).eq.isc.and.isblk(2,is).eq.isd.and.          8d11s23
     $            isblk(3,is).eq.isa)then                               8d7s23
              jj=jmats(is)                                              8d7s23
              i10=i1s                                                   8d7s23
              i1n=nvirt(isblk(3,is))                                    8d7s23
              do i2=i2s,i2e                                             8d7s23
               i2p=i2+noc(isblk(4,is))-1                                8d7s23
               if(i2.eq.i2e)i1n=i1e                                     8d7s23
               do i1=i10,i1n                                            8d7s23
                i1p=i1+noc(isblk(3,is))-1                               8d7s23
                do i4=0,noc(isblk(2,is))-1                              8d7s23
                 do i3=0,noc(isblk(1,is))-1                             8d7s23
                  irec=i3+noc(isblk(1,is))*i4                           8d7s23
                  ix=max(i3,i4)                                         8d7s23
                  in=min(i3,i4)                                         8d7s23
                  itri=((ix*(ix+1))/2)+in                               8d7s23
                  itri=itri+isw*(irec-itri)+jj                          8d7s23
                  iad=itmp+i1p+nbasdws(isa)*(i2p+nbasdws(isb)*(i3       8d9s23
     $                 +nbasdws(isc)*i4))                               8d7s23
                  bc(iad)=bc(itri)                                      8d7s23
                  if(iad.eq.igozl)write(6,*)('gozld '),bc(igozl)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                jj=jj+nrowj                                              8d7s23
               end do                                                   8d7s23
               i10=1                                                    8d7s23
              end do                                                    8d7s23
             else if(isblk(2,is).eq.isc.and.isblk(1,is).eq.isd.and.     8d7s23
     $            isblk(3,is).eq.isa)then                               8d7s23
              jj=jmats(is)                                              8d7s23
              i10=i1s                                                   8d7s23
              i1n=nvirt(isblk(3,is))                                    8d7s23
              do i2=i2s,i2e                                             8d7s23
               i2p=i2+noc(isblk(4,is))-1                                8d7s23
               if(i2.eq.i2e)i1n=i1e                                     8d7s23
               do i1=i10,i1n                                            8d7s23
                i1p=i1+noc(isblk(3,is))-1                               8d7s23
                do i4=0,noc(isblk(2,is))-1                              8d7s23
                 do i3=0,noc(isblk(1,is))-1                             8d7s23
                  irec=jj+i3+noc(isblk(1,is))*i4                        8d7s23
                  iad=itmp+i1p+nbasdws(isa)*(i2p+nbasdws(isb)*(i4       8d9s23
     $                 +nbasdws(isc)*i3))                               8d7s23
                  bc(iad)=bc(irec)                                      8d7s23
                  if(iad.eq.igozl)write(6,*)('gozlc '),bc(igozl)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                jj=jj+nrowj                                              8d7s23
               end do                                                   8d7s23
               i10=1                                                    8d7s23
              end do                                                    8d7s23
             end if                                                     8d7s23
            end if                                                      8d7s23
           end do                                                       8d7s23
           do is=1,nsdlk                                                8d7s23
            call ilimts(noc(isblk(3,is)),noc(isblk(4,is)),mynprocg,     8d7s23
     $           mynowprog,il,ih,i1s,i1e,i2s,i2e)                       8d7s23
            nhere=ih+1-il                                               8d7s23
            if(nhere.gt.0)then                                          8d7s23
             if(isblk(1,is).eq.isblk(2,is))then                         8d7s23
              nrow4=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2            8d7s23
              isw=0                                                     8d7s23
             else                                                       8d7s23
              nrow4=noc(isblk(1,is))*noc(isblk(2,is))                    8d7s23
              isw=1                                                     8d7s23
             end if
             if(isblk(1,is).eq.isa.and.isblk(2,is).eq.isb.and.          8d7s23
     $          isblk(3,is).eq.isc)then                                 8d7s23
              ii=ioooo(is)                                              8d7s23
              i10=i1s                                                   8d7s23
              i1n=noc(isblk(3,is))                                      8d7s23
              do i2=i2s,i2e                                             8d7s23
               i2m=i2-1                                                 8d7s23
               if(i2.eq.i2e)i1n=i1e                                     8d7s23
               do i1=i10,i1n                                            8d7s23
                i1m=i1-1                                                8d7s23
                do i4=0,noc(isblk(2,is))-1                              8d7s23
                 do i3=0,noc(isblk(1,is))-1                             8d7s23
                  irec=i3+noc(isblk(1,is))*i4                           8d7s23
                  ix=max(i3,i4)                                         8d7s23
                  in=min(i3,i4)                                         8d7s23
                  itri=((ix*(ix+1))/2)+in                               8d7s23
                  itri=itri+isw*(irec-itri)+ii                          8d7s23
                  iad=itmp+i3+nbasdws(isa)*(i4+nbasdws(isb)*(i1m        8d7s23
     $                 +nbasdws(isc)*i2m))                              8d7s23
                  bc(iad)=bc(itri)                                      8d7s23
                  if(iad.eq.igozl)write(6,*)('gozlb '),bc(igozl)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                ii=ii+nrow4                                              8d7s23
               end do                                                   8d7s23
               i10=1                                                    8d7s23
              end do                                                    8d7s23
             else if(isblk(2,is).eq.isa.and.isblk(1,is).eq.isb.and.     8d7s23
     $          isblk(3,is).eq.isc)then                                 8d7s23
              ii=ioooo(is)                                              8d7s23
              i10=i1s                                                   8d7s23
              i1n=noc(isblk(3,is))                                      8d7s23
              do i2=i2s,i2e                                             8d7s23
               i2m=i2-1                                                 8d7s23
               if(i2.eq.i2e)i1n=i1e                                     8d7s23
               do i1=i10,i1n                                            8d7s23
                i1m=i1-1                                                8d7s23
                do i4=0,noc(isblk(2,is))-1                              8d7s23
                 do i3=0,noc(isblk(1,is))-1                             8d7s23
                  irec=ii+i3+noc(isblk(1,is))*i4                        8d7s23
                  iad=itmp+i4+nbasdws(isa)*(i3+nbasdws(isb)*(i1m        8d7s23
     $                 +nbasdws(isc)*i2m))                              8d7s23
                  bc(iad)=bc(irec)                                      8d7s23
                  if(iad.eq.igozl)write(6,*)('gozla '),bc(igozl)
                 end do                                                 8d7s23
                end do                                                  8d7s23
                ii=ii+nrow4                                              8d7s23
               end do                                                   8d7s23
               i10=1                                                    8d7s23
              end do                                                    8d7s23
             end if                                                     8d7s23
            end if                                                      8d7s23
           end do                                                       8d7s23
           do i4=1,2                                                    8d9s23
            do i3=1,2                                                   8d9s23
             do i2=1,2                                                  8d9s23
              do i1=1,2                                                 8d9s23
               nnqn(i1,i2,i3,i4)=0                                      8d9s23
              end do                                                    8d9s23
             end do                                                     8d9s23
            end do                                                      8d9s23
           end do                                                       8d9s23
           do i4=0,nbasdws(isd)-1
            do i3=0,nbasdws(isc)-1
             do i2=0,nbasdws(isb)-1
              do i1=0,nbasdws(isa)-1
               iad=itmp+i1+nbasdws(isa)*(i2+nbasdws(isb)*(i3
     $              +nbasdws(isc)*i4))
               if(bc(iad).ne.bc(iad))then
                j1=1
                if(i1.ge.noc(isa))j1=2
                j2=1
                if(i2.ge.noc(isb))j2=2
                j3=1
                if(i3.ge.noc(isc))j3=2
                j4=1
                if(i4.ge.noc(isd))j4=2
                nnqn(j1,j2,j3,j4)=nnqn(j1,j2,j3,j4)+1
               end if
              end do
             end do
            end do
           end do
           do i4=1,2                                                    8d9s23
            do i3=1,2                                                   8d9s23
             do i2=1,2                                                  8d9s23
              do i1=1,2                                                 8d9s23
               if(nnqn(i1,i2,i3,i4).gt.0)then                           8d9s23
               end if                                                   8d9s23
              end do                                                    8d9s23
             end do                                                     8d9s23
            end do                                                      8d9s23
           end do                                                       8d9s23
           call e0xd(bc(itmp),isa,isb,isc,isd,nbasdws,idoubo,irefo,     5d2s24
     $         isblk,nsdlk,e0x,id4o,bc,ibc,iuseage)                             5d2s24
           call e1xd(bc(itmp),isa,isb,isc,isd,nbasdws,idoubo,irefo,     5d2s24
     $         nvirt,noc,isblk1,nsdlk1,e1x,id1x,bc,ibc,iuseage(1,2))                 5d2s24
           call ejjd(bc(itmp),isa,isb,isc,isd,nbasdws,idoubo,irefo,     5d2s24
     $         nvirt,noc,isblk,nsdlk,ejj,jmden,bc,ibc,iuseage(1,3),0,   5d6s24
     $          xjt)
           call ekkd(bc(itmp),isa,isb,isc,isd,nbasdws,idoubo,irefo,     5d2s24
     $        nvirt,noc,isblkk,nsdlkk,ekk,kmden,bc,ibc,iuseage(1,4))    5d6s24
           call e3xd(bc(itmp),isa,isb,isc,isd,nbasdws,idoubo,irefo,     5d2s24
     $         nvirt,noc,isblk1,nsdlk1,e3x,id3x,bc,ibc,iuseage(1,5))                 5d2s24
           if(ndoub.gt.0)then                                           6d12s24
            call e4vd(bc(itmp),isa,isb,isc,isd,nbasdws,nvirt,noc,e4x,    5d10s24
     $       vdinout,nfdat,multh,nsymb,isymmrci,bc,ibc,iuseage(1,6),srh)5d10s24
            call e4tvd(bc(itmp),isa,isb,isc,isd,nbasdws,nvirt,noc,ie4v,  5d15s24
     $       vdinout,nfdat,multh,nsymb,isymmrci,bc,ibc,iuseage(1,6),srh,6d2s24
     $          itmp)
           end if                                                       6d12s24
           itmpx=ibcoff                                                 5d2s24
           itmpy=itmpx+nrow*ncol                                        5d2s24
           itmpz=itmpy+nrow*ncol                                        5d2s24
           ibcoff=itmpz+nrow*ncol                                       5d2s24
           call enough('playr.xy',bc,ibc)                               5d2s24
           do ikind=1,npartp                                            6d26s24
            nrowx=nrow*nbasdws(isc)
            if(ikind.eq.npartp)then                                     6d26s24
             call dgemm('n','n',nrowx,nbasdws(isd),nbasdws(isd),1d0,     5d2s24
     $          bc(itmp),nrowx,darot(iarot(isd)),nbasdws(isd),0d0,      6d26s24
     $          bc(itmpz),nrowx,'playr.d')                              5d2s24
            else                                                        6d26s24
             call dgemm('n','n',nrowx,nbasdws(isd),nbasdws(isd),1d0,     5d2s24
     $          bc(itmp),nrowx,bc(idaprt(ikind,isd)),nbasdws(isd),0d0,  5d2s24
     $          bc(itmpz),nrowx,'playr.d')                              5d2s24
            end if                                                      6d26s24
            do ic=0,nbasdws(isc)-1
             do id=0,nbasdws(isd)-1
              iad1=itmp+nrow*(ic+nbasdws(isc)*id)
              iad2=itmpx+nrow*(id+nbasdws(isd)*ic)
              do jab=0,nrow-1
               bc(iad2+jab)=bc(iad1+jab)
              end do
             end do
            end do
            nrowx=nrow*nbasdws(isd)
            if(ikind.eq.npartp)then                                     6d26s24
             call dgemm('n','n',nrowx,nbasdws(isc),nbasdws(isc),1d0,     5d2s24
     $          bc(itmpx),nrowx,darot(iarot(isc)),nbasdws(isc),0d0,     6d26s24
     $          bc(itmpy),nrowx,'playr.c')                              5d2s24
            else
             call dgemm('n','n',nrowx,nbasdws(isc),nbasdws(isc),1d0,     5d2s24
     $          bc(itmpx),nrowx,bc(idaprt(ikind,isc)),nbasdws(isc),0d0,  5d2s24
     $          bc(itmpy),nrowx,'playr.c')                              5d2s24
            end if                                                      6d26s24
            do ic=0,nbasdws(isc)-1                                      5d2s24
             do id=0,nbasdws(isd)-1                                     5d9s24
              iad1=itmpz+nrow*(ic+nbasdws(isc)*id)                      5d2s24
              iad2=itmpy+nrow*(id+nbasdws(isd)*ic)                      5d2s24
              do jab=0,nrow-1                                           5d2s24
               bc(iad1+jab)=bc(iad1+jab)+bc(iad2+jab)                   5d2s24
              end do                                                    5d2s24
             end do                                                     5d2s24
            end do                                                      5d2s24
            do ib=0,nbasdws(isb)-1                                      5d2s24
             do ia=0,nbasdws(isa)-1                                     5d2s24
              do jcd=0,ncol-1                                           5d2s24
               iad1=itmp+ia+nbasdws(isa)*(ib+nbasdws(isb)*jcd)          5d2s24
               iad2=itmpx+jcd+ncol*(ia+nbasdws(isa)*ib)                 5d2s24
               bc(iad2)=bc(iad1)                                        5d2s24
              end do                                                    5d2s24
             end do                                                     5d2s24
            end do                                                      5d2s24
            nrowx=ncol*nbasdws(isa)
            if(ikind.eq.npartp)then                                     6d26s24
             call dgemm('n','n',nrowx,nbasdws(isb),nbasdws(isb),1d0,     5d2s24
     $          bc(itmpx),nrowx,darot(iarot(isb)),nbasdws(isb),0d0,     6d26s24
     $          bc(itmpy),nrowx,'playr.b')                              5d2s24
            else                                                        6d26s24
             call dgemm('n','n',nrowx,nbasdws(isb),nbasdws(isb),1d0,     5d2s24
     $          bc(itmpx),nrowx,bc(idaprt(ikind,isb)),nbasdws(isb),0d0,  5d2s24
     $          bc(itmpy),nrowx,'playr.b')                              5d2s24
            end if                                                      6d26s24
            do ib=0,nbasdws(isb)-1                                      5d2s24
             do ia=0,nbasdws(isa)-1                                     5d2s24
              do jcd=0,ncol-1                                           5d2s24
               iad1=itmpz+ia+nbasdws(isa)*(ib+nbasdws(isb)*jcd)          5d2s24
               iad2=itmpy+jcd+ncol*(ia+nbasdws(isa)*ib)                 5d2s24
               bc(iad1)=bc(iad2)+bc(iad1)                                        5d2s24
               iad1=itmp+ia+nbasdws(isa)*(ib+nbasdws(isb)*jcd)
               iad2=itmpx+jcd+ncol*(ib+nbasdws(isb)*ia)
               bc(iad2)=bc(iad1)
              end do                                                    5d2s24
             end do                                                     5d2s24
            end do                                                      5d2s24
            nrowx=ncol*nbasdws(isb)
            if(ikind.eq.npartp)then                                     6d26s24
             call dgemm('n','n',nrowx,nbasdws(isa),nbasdws(isa),1d0,     5d2s24
     $          bc(itmpx),nrowx,darot(iarot(isa)),nbasdws(isa),0d0,     6d26s24
     $          bc(itmpy),nrowx,'playr.a')                              5d2s24
            else                                                        6d26s24
            call dgemm('n','n',nrowx,nbasdws(isa),nbasdws(isa),1d0,     5d2s24
     $          bc(itmpx),nrowx,bc(idaprt(ikind,isa)),nbasdws(isa),0d0,  5d2s24
     $          bc(itmpy),nrowx,'playr.a')                              5d2s24
            end if                                                      6d26s24
            do ib=0,nbasdws(isb)-1                                      5d2s24
             do ia=0,nbasdws(isa)-1                                     5d2s24
              do jcd=0,ncol-1                                           5d2s24
               iad1=itmpz+ia+nbasdws(isa)*(ib+nbasdws(isb)*jcd)          5d2s24
               iad2=itmpy+jcd+ncol*(ib+nbasdws(isb)*ia)                 5d2s24
               bc(iad1)=bc(iad2)+bc(iad1)                                        5d2s24
              end do                                                    5d2s24
             end do                                                     5d2s24
            end do                                                      5d2s24
            call e0xd(bc(itmpz),isa,isb,isc,isd,nbasdws,idoubo,irefo,     5d2s24
     $         isblk,nsdlk,de0x(ikind),id4o,bc,ibc,iuseage)                             5d2s24
            call e1xd(bc(itmpz),isa,isb,isc,isd,nbasdws,idoubo,irefo,     5d2s24
     $         nvirt,noc,isblk1,nsdlk1,de1x(ikind),id1x,bc,ibc,         5d6s24
     $           iuseage(1,2))                                          5d6s24
             itcode=0
            call ejjd(bc(itmpz),isa,isb,isc,isd,nbasdws,idoubo,irefo,     5d2s24
     $      nvirt,noc,isblk,nsdlk,dejj(ikind),jmden,bc,ibc,iuseage(1,3),5d6s24
     $           itcode,xjt)
            call ekkd(bc(itmpz),isa,isb,isc,isd,nbasdws,idoubo,irefo,     5d2s24
     $    nvirt,noc,isblkk,nsdlkk,dekk(ikind),kmden,bc,ibc,iuseage(1,4))    5d6s24
            call e3xd(bc(itmpz),isa,isb,isc,isd,nbasdws,idoubo,irefo,     5d2s24
     $     nvirt,noc,isblk1,nsdlk1,de3x(ikind),id3x,bc,ibc,iuseage(1,5))5d10s24
            if(ndoub.gt.0)then                                          6d12s24
             call e4vd(bc(itmpz),isa,isb,isc,isd,nbasdws,nvirt,noc,
     $           de4x(ikind),
     $       vdinout,nfdat,multh,nsymb,isymmrci,bc,ibc,iuseage(1,6),srh)5d10s24
            end if                                                      6d12s24
            if(isa.eq.isb.and.min(idoubo(isa),idoubo(isc)).gt.0)then     5d2s24
             ff=2d0                                                      5d10s24
             if(isa.ne.isc)ff=4d0                                        5d10s24
             do jcd=0,idoubo(isc)-1                                      5d2s24
              do jab=0,idoubo(isa)-1                                      5d2s24
               iad=itmpz+jab+nbasdws(isa)*(jab+nbasdws(isa)*(jcd          5d2s24
     $             +nbasdws(isc)*jcd))                                  5d2s24
               dcore2(ikind)=dcore2(ikind)+ff*bc(iad)                   5d10s24
              end do                                                     5d2s24
             end do                                                      5d2s24
            end if                                                       5d2s24
            if(isa.eq.isc.and.min(idoubo(isa),idoubo(isb)).gt.0)then     5d2s24
             ff=1d0                                                      5d2s24
             if(isa.ne.isb)ff=2d0                                       5d10s24
             do ibd=0,idoubo(isb)-1                                      5d2s24
              do iac=0,idoubo(isa)-1                                     5d2s24
               iad=itmpz+iac+nbasdws(isa)*(ibd+nbasdws(isb)*(iac          5d2s24
     $             +nbasdws(isc)*ibd))                                  5d2s24
               dcore2k(ikind)=dcore2k(ikind)-ff*bc(iad)                                  5d2s24
              end do                                                     5d2s24
             end do                                                      5d2s24
            end if                                                       5d2s24
            if(isa.eq.isc)then                                           5d2s24
             if(min(idoubo(isa),nh0av(isd)).gt.0)then
              do iac=0,idoubo(isa)-1                                      5d2s24
               do id=0,nh0av(isd)-1
                idp=id+idoubo(isd)
                do ib=0,nh0av(isb)-1
                 ibp=ib+idoubo(isb)
                 iad1=iden(isb)+ib+nh0av(isb)*id
                 iad2=itmpz+iac+nbasdws(isa)*(ibp+nbasdws(isb)*(iac        5d2s24
     $               +nbasdws(isc)*idp))
                 deoa(ikind)=deoa(ikind)-bc(iad1)*bc(iad2)
                end do
               end do
              end do
             end if
             if(isa.ne.isb)then                                          5d2s24
              if(min(idoubo(isd),nh0av(isc)).gt.0)then
               do ibd=0,idoubo(isd)-1                                     5d2s24
                do ic=0,nh0av(isc)-1                                      5d2s24
                 icp=ic+idoubo(isc)                                       5d2s24
                 do ia=0,nh0av(isa)-1                                     5d2s24
                  iap=ia+idoubo(isa)                                      5d2s24
                  iad1=iden(isa)+ia+nh0av(isa)*ic
                  iad2=itmpz+iap+nbasdws(isa)*(ibd+nbasdws(isb)*(icp       5d2s24
     $                +nbasdws(isc)*ibd))                                5d2s24
                  deoa(ikind)=deoa(ikind)-bc(iad1)*bc(iad2)                               5d2s24
                 end do                                                   5d2s24
                end do                                                    5d2s24
               end do                                                     5d2s24
              end if
             end if                                                      5d2s24
            end if                                                       5d2s24
            if(isa.eq.isb)then
             if(min(idoubo(isa),nh0av(isd)).gt.0)then
              do id=0,nh0av(isd)-1
               idp=id+idoubo(isd)
               do ic=0,nh0av(isc)-1
                icp=ic+idoubo(isc)
                iad1=iden(isd)+ic+nh0av(isc)*id
                do jab=0,idoubo(isa)-1
                 iad2=itmpz+jab+nbasdws(isa)*(jab+nbasdws(isa)*(icp
     $               +nbasdws(isc)*idp))
                 deoa(ikind)=deoa(ikind)+2d0*bc(iad2)*bc(iad1)
                end do
               end do
              end do
             end if
             if(isc.ne.isa)then
              if(min(idoubo(isc),nh0av(isa)).gt.0)then
               do jcd=0,idoubo(isc)-1
                do ib=0,nh0av(isb)-1
                 ibp=ib+idoubo(isb)
                 do ia=0,nh0av(isa)-1
                  iap=ia+idoubo(isa)
                  iad1=iden(isa)+ia+nh0av(isa)*ib
                  iad2=itmpz+iap+nbasdws(isa)*(ibp+nbasdws(isb)*(jcd
     $                +nbasdws(isc)*jcd))
                  deoa(ikind)=deoa(ikind)+2d0*bc(iad2)*bc(iad1)
                 end do
                end do
               end do
              end if
             end if
            end if
           end do
           ibcoff=itmpx
           if(isa.eq.isc)then                                           5d2s24
            if(min(idoubo(isa),nh0av(isd)).gt.0)then
             write(6,*)('eoaK '),isa,isd,isa,isb,isc,isd
             do iac=0,idoubo(isa)-1                                      5d2s24
              do id=0,nh0av(isd)-1
               idp=id+idoubo(isd)
               do ib=0,nh0av(isb)-1
                ibp=ib+idoubo(isb)
                iad1=iden(isb)+ib+nh0av(isb)*id
                iad2=itmp+iac+nbasdws(isa)*(ibp+nbasdws(isb)*(iac        5d2s24
     $               +nbasdws(isc)*idp))
                eoa=eoa-bc(iad1)*bc(iad2)
                eoad(isb)=eoad(isb)-bc(iad1)*bc(iad2)
               end do
              end do
             end do
            end if
            if(isa.ne.isb)then                                          5d2s24
             if(min(idoubo(isd),nh0av(isc)).gt.0)then
              write(6,*)('eoaK '),isd,isa,isa,isb,isc,isd
              do ibd=0,idoubo(isd)-1                                     5d2s24
               do ic=0,nh0av(isc)-1                                      5d2s24
                icp=ic+idoubo(isc)                                       5d2s24
                do ia=0,nh0av(isa)-1                                     5d2s24
                 iap=ia+idoubo(isa)                                      5d2s24
                 iad1=iden(isa)+ia+nh0av(isa)*ic
                 iad2=itmp+iap+nbasdws(isa)*(ibd+nbasdws(isb)*(icp       5d2s24
     $                +nbasdws(isc)*ibd))                                5d2s24
                 eoa=eoa-bc(iad1)*bc(iad2)                               5d2s24
                 eoad(isa)=eoad(isa)-bc(iad1)*bc(iad2)                               5d2s24
                end do                                                   5d2s24
               end do                                                    5d2s24
              end do                                                     5d2s24
             end if
            end if                                                      5d2s24
           end if                                                       5d2s24
           if(isa.eq.isb)then
            if(min(idoubo(isa),nh0av(isd)).gt.0)then
             write(6,*)('eoaJ '),isa,isd,isa,isb,isc,isd
             do id=0,nh0av(isd)-1
              idp=id+idoubo(isd)
              do ic=0,nh0av(isc)-1
               icp=ic+idoubo(isc)
               iad1=iden(isd)+ic+nh0av(isc)*id
               do jab=0,idoubo(isa)-1
                iad2=itmp+jab+nbasdws(isa)*(jab+nbasdws(isa)*(icp
     $               +nbasdws(isc)*idp))
                eoa=eoa+2d0*bc(iad2)*bc(iad1)
                eoad(isd)=eoad(isd)+2d0*bc(iad2)*bc(iad1)
               end do
              end do
             end do
            end if
            if(isc.ne.isa)then
             if(min(idoubo(isc),nh0av(isa)).gt.0)then
              write(6,*)('eoaJ '),isc,isa,isa,isb,isc,isd
              do jcd=0,idoubo(isc)-1
               do ib=0,nh0av(isb)-1
                ibp=ib+idoubo(isb)
                do ia=0,nh0av(isa)-1
                 iap=ia+idoubo(isa)
                 iad1=iden(isa)+ia+nh0av(isa)*ib
                 iad2=itmp+iap+nbasdws(isa)*(ibp+nbasdws(isb)*(jcd
     $                +nbasdws(isc)*jcd))
                 eoa=eoa+2d0*bc(iad2)*bc(iad1)
                 eoad(isa)=eoad(isa)+2d0*bc(iad2)*bc(iad1)
                end do
               end do
              end do
             end if
            end if
           end if
           if(isa.eq.isb.and.min(idoubo(isa),idoubo(isc)).gt.0)then     5d2s24
            ff=2d0                                                      5d10s24
            if(isa.ne.isc)ff=4d0                                        5d10s24
            do jcd=0,idoubo(isc)-1                                      5d2s24
             do jab=0,idoubo(isa)-1                                      5d2s24
              iad=itmp+jab+nbasdws(isa)*(jab+nbasdws(isa)*(jcd          5d2s24
     $             +nbasdws(isc)*jcd))                                  5d2s24
              ecore2=ecore2+ff*bc(iad)
             end do                                                     5d2s24
            end do                                                      5d2s24
           end if                                                       5d2s24
           if(isa.eq.isc.and.min(idoubo(isa),idoubo(isb)).gt.0)then     5d2s24
            ff=1d0                                                      5d2s24
            if(isa.ne.isb)ff=2d0                                        5d2s24
            do ibd=0,idoubo(isb)-1                                      5d2s24
             do iac=0,idoubo(isa)-1                                     5d2s24
              iad=itmp+iac+nbasdws(isa)*(ibd+nbasdws(isb)*(iac          5d2s24
     $             +nbasdws(isc)*ibd))                                  5d2s24
              ecore2k=ecore2k-ff*bc(iad)                                  5d2s24
             end do                                                     5d2s24
            end do                                                      5d2s24
           end if                                                       5d2s24
           itmpt=ibcoff                                                 8d7s23
           ibcoff=itmpt+nrow*ncol                                       8d7s23
           call enough('playr.tmpt.4x',bc,ibc)                          8d7s23
           do ikind=1,npartp                                            6d26s24
            write(6,*)('for ikind= '),ikind                             8d7s23
            do iperm=1,4
             write(6,*)('for iperm: '),iperm
             igotit=0
             if(iperm.eq.1)then                                         8d8s23
              nmul=nrow*nbasdws(isc)                                    8d9s23
              write(6,*)('transform 4th index ...'),isa,isb,isc,isd                        8d7s23
              if(ikind.eq.npartp)then                                   6d26s24
               call dgemm('n','n',nmul,nbasdws(isd),nbasdws(isd),1d0,       8d7s23
     $          bc(itmp),nmul,darot(iarot(isd)),nbasdws(isd),0d0,       6d26s24
     $           bc(itmpt),nmul,'playr.tmpt4x')                         8d7s23
              else                                                      6d26s24
               call dgemm('n','n',nmul,nbasdws(isd),nbasdws(isd),1d0,       8d7s23
     $          bc(itmp),nmul,bc(idaprt(ikind,isd)),nbasdws(isd),0d0,   8d7s23
     $           bc(itmpt),nmul,'playr.tmpt4x')                         8d7s23
              end if
              if(itmpt.le.igoxl.and.itmpt+nmul*nbasdws(isd).ge.igoxl)
     $             then
               write(6,*)('igoxla')
               idelta=igoxl-itmpt
               ii2=idelta/nmul
               ii1=idelta-nmul*ii2
               ij4=ii1/nrow
               ij3=ii1-nrow*ij4
               ij2=ij3/nbasdws(isa)
               ij1=ij3-nbasdws(isa)*ij2
               write(6,*)ii1,nrow,ij4
               write(6,*)ij3,nbasdws(isa),ij2,ij1
               ij3=ij1
               ij2=ij2+1
               ij3=ij3+1
               ij4=ij4+1
               sum=0d0
               do i=0,nbasdws(isd)-1
                iad1=itmp+ii1+nmul*i
                if(ikind.eq.npartp)then                                 6d26s24
                 iad2=iarot(isd)+i+nbasdws(isd)*ii2                     6d26s24
                 value=darot(iad2)                                      6d26s24
                else                                                    6d26s24
                 iad2=idaprt(ikind,isd)+i+nbasdws(isd)*ii2
                 value=bc(iad2)                                         6d26s24
                end if                                                  6d26s24
                term=bc(iad1)*value                                     6d26s24
                sum=sum+term
                if(abs(term).gt.1d-10)then
                 write(6,*)bc(iad1),value,sum,
     $               iad1,iad2,ii1,i,ii2,nmul,nrow
                end if
               end do
              end if
              jsa=isa
              jsb=isb
              jsc=isc
              jsd=isd
              igotit=1
             else if(iperm.eq.2)then
              if(isc.ne.isd)then                                        8d8s23
               write(6,*)('transform 3rd index ')                       8d8s23
               itnew=ibcoff                                             8d8s23
               ibcoff=itnew+nrow*ncol                                   8d8s23
               call enough('playr.tnew',bc,ibc)                         8d8s23
               do i4=0,nbasdws(isd)-1                                   8d8s23
                do i3=0,nbasdws(isc)-1                                  8d8s23
                 i34=itmp+nrow*(i3+nbasdws(isc)*i4)                     8d8s23
                 i43=itnew+nrow*(i4+nbasdws(isd)*i3)                    8d8s23
                 do i12=0,nrow-1                                        8d8s23
                  bc(i43+i12)=bc(i34+i12)                               8d8s23
                 end do                                                 8d8s23
                end do                                                  8d8s23
               end do                                                   8d8s23
               nmul=nrow*nbasdws(isd)                                   8d8s23
               if(ikind.eq.npartp)then                                  6d26s24
                call dgemm('n','n',nmul,nbasdws(isc),nbasdws(isc),1d0,       8d7s23
     $          bc(itnew),nmul,darot(iarot(isc)),nbasdws(isc),0d0,      6d26s24
     $           bc(itmpt),nmul,'playr.tmpt4x')                         8d7s23
               else                                                     6d26s24
                call dgemm('n','n',nmul,nbasdws(isc),nbasdws(isc),1d0,       8d7s23
     $          bc(itnew),nmul,bc(idaprt(ikind,isc)),nbasdws(isc),0d0,   8d7s23
     $           bc(itmpt),nmul,'playr.tmpt4x')                         8d7s23
               end if                                                   6d26s24
              if(ikind.eq.1)then
               do i=0,nbasdws(isc)-1
                do j=0,nbasdws(isc)-1
                 do k=0,nmul-1
                  iad1=itnew+k+nmul*j
                  iad2=idaprt(ikind,isc)+j+nbasdws(isc)*i
                  iad3=itmpt+k+nmul*i
                  term=bc(iad1)*bc(iad2)
                 end do
                end do
               end do
              end if
              if(itmpt.le.igoxl.and.itmpt+nmul*nbasdws(isc).gt.igoxl)
     $             then
               write(6,*)('we multiplied ')
               call prntm2(bc(itnew),nbasdws(isa)*nbasdws(isb),
     $              nbasdws(isc)*nbasdws(isd),nbasdws(isa)*nbasdws(isb))
               write(6,*)('by ')
               if(ikind.eq.npartp)then                                  6d26s24
                call prntm2(darot(iarot(isc)),nbasdws(isc),             6d26s24
     $              nbasdws(isc),nbasdws(isc))
               else                                                     6d26s24
                call prntm2(bc(idaprt(ikind,isc)),nbasdws(isc),
     $              nbasdws(isc),nbasdws(isc))
               end if                                                   6d26s24
               write(6,*)('to yield ')
               call prntm2(bc(itmpt),nbasdws(isa)*nbasdws(isb),
     $              nbasdws(isc)*nbasdws(isd),nbasdws(isa)*nbasdws(isb))
               do ic=0,nbasdws(isc)-1
                do id=0,nbasdws(isd)-1
                 do ib=0,nbasdws(isb)-1
                  do ia=0,nbasdws(isa)-1
                   iad=itmpt+ia+nbasdws(isa)*(ib+nbasdws(isb)*(id
     $                  +nbasdws(isd)*ic))
                   if(abs(bc(iad)).gt.1d-10)write(6,*)ia,ia-noc(isa),
     $                  ib,ib-noc(isb),id,id-noc(isd),ic,ic-noc(isc),
     $                  bc(iad),iad,iad-itmpt
                  end do
                 end do
                end do
               end do
               write(6,*)('igoxlb')
               idelta=igoxl-itmpt
               ii2=idelta/nmul
               ii1=idelta-nmul*ii2
               sum=0d0
               do i=0,nbasdws(isc)-1
                iad1=itnew+ii1+nmul*i
                if(ikind.eq.npartp)then                                 6d26s24
                 iad2=iarot(isc)+i+nbasdws(isc)*ii2                     6d26s24
                 value=darot(iad2)                                      6d26s24
                else                                                    6d26s24
                 iad2=idaprt(ikind,isc)+i+nbasdws(isc)*ii2
                 value=bc(iad2)                                         6d26s24
                end if                                                  6d26s24
                term=bc(iad1)*value                                     6d26s24
                sum=sum+term
                if(abs(term).gt.-1d-10)write(6,*)bc(iad1),value,sum,
     $               iad1,iad2,ii1,i,ii2,nmul,nrow
               end do
              end if
               ibcoff=itnew                                             8d8s23
               jsa=isa
               jsb=isb
               jsc=isd
               jsd=isc
               igotit=1
              end if                                                    8d8s23
             else if(iperm.eq.3)then
              if(isb.ne.isc.and.isb.ne.isd)then                         8d8s23
               write(6,*)('transform 2nd index ')                       8d8s23
               itnew=ibcoff                                             8d8s23
               ibcoff=itnew+nrow*ncol                                   8d8s23
               call enough('playr.tnew',bc,ibc)                         8d8s23
               do i34=0,ncol-1                                          8d8s23
                do i2=0,nbasdws(isb)-1                                  8d8s23
                 do i1=0,nbasdws(isa)-1                                 8d8s23
                  iado=itmp+i1+nbasdws(isa)*(i2+nbasdws(isb)*i34)       8d8s23
                  iadn=itnew+i34+ncol*(i1+nbasdws(isa)*i2)              8d8s23
                  bc(iadn)=bc(iado)                                     8d8s23
                 end do                                                 8d8s23
                end do                                                  8d8s23
               end do                                                   8d8s23
               nmul=ncol*nbasdws(isa)
               if(ikind.eq.npartp)then                                  6d26s24
                call dgemm('n','n',nmul,nbasdws(isb),nbasdws(isb),1d0,       8d7s23
     $          bc(itnew),nmul,darot(iarot(isb)),nbasdws(isb),0d0,      6d26s24
     $           bc(itmpt),nmul,'playr.tmpt4x')                         8d7s23
               else
                call dgemm('n','n',nmul,nbasdws(isb),nbasdws(isb),1d0,       8d7s23
     $          bc(itnew),nmul,bc(idaprt(ikind,isb)),nbasdws(isb),0d0,   8d7s23
     $           bc(itmpt),nmul,'playr.tmpt4x')                         8d7s23
               end if                                                   6d26s24
              if(ikind.eq.1)then
               do i=0,nbasdws(isb)-1
                do j=0,nbasdws(isb)-1
                 do k=0,nmul-1
                  iad1=itnew+k+nmul*j
                  iad2=idaprt(ikind,isb)+j+nbasdws(isb)*i
                  iad3=itmpt+k+nmul*i
                  term=bc(iad1)*bc(iad2)
                 end do
                end do
               end do
              end if
              if(itmpt.le.igoxl.and.itmpt+nmul*nbasdws(isb).gt.igoxl)
     $             then
               write(6,*)('igoxlc')
               idelta=igoxl-itmpt
               ii2=idelta/nmul
               ii1=idelta-nmul*ii2
               sum=0d0
               do i=0,nbasdws(isb)-1
                iad1=itnew+ii1+nmul*i
                if(ikind.eq.npartp)then                                 6d26s24
                 iad2=iarot(isb)+i+nbasdws(isb)*ii2                     6d26s24
                 value=darot(iad2)                                      6d26s24
                else                                                    6d26s24
                 iad2=idaprt(ikind,isb)+i+nbasdws(isb)*ii2
                 value=bc(iad2)                                         6d26s24
                end if                                                  6d26s24
                term=bc(iad1)*value                                     6d26s24
                sum=sum+term
                if(abs(term).gt.1d-10)write(6,*)bc(iad1),value,sum,
     $               iad1,iad2,ii1,i,ii2,nmul,nrow
               end do
              end if
               ibcoff=itnew                                             8d8s23
               jsa=isc                                                  8d8s23
               jsb=isd
               jsc=isa
               jsd=isb
               igotit=1
              end if                                                    8d8s23
             else if(iperm.eq.4)then
              if(isa.ne.isb.and.isa.ne.isc.and.isa.ne.isd)then          8d8s23
               write(6,*)('transform 1st index ')                       8d8s23
               itnew=ibcoff                                             8d8s23
               ibcoff=itnew+nrow*ncol                                   8d8s23
               call enough('playr.tnew',bc,ibc)                         8d8s23
               do i34=0,ncol-1                                          8d8s23
                do i2=0,nbasdws(isb)-1                                  8d8s23
                 do i1=0,nbasdws(isa)-1                                 8d8s23
                  iado=itmp+i1+nbasdws(isa)*(i2+nbasdws(isb)*i34)       8d8s23
                  iadn=itnew+i34+ncol*(i2+nbasdws(isb)*i1)              8d8s23
                  bc(iadn)=bc(iado)                                     8d8s23
                  if(iadn.eq.2019039)write(6,*)('getting '),bc(iado),
     $                 ('from '),iado,i1,i2,i34
                 end do                                                 8d8s23
                end do                                                  8d8s23
               end do                                                   8d8s23
               nmul=ncol*nbasdws(isb)
               if(ikind.eq.npartp)then                                  6d26s24
                call dgemm('n','n',nmul,nbasdws(isa),nbasdws(isa),1d0,       8d7s23
     $          bc(itnew),nmul,darot(iarot(isa)),nbasdws(isa),0d0,      6d26s24
     $           bc(itmpt),nmul,'playr.tmpt4x')                         8d7s23
               else                                                     6d26s24
                call dgemm('n','n',nmul,nbasdws(isa),nbasdws(isa),1d0,       8d7s23
     $          bc(itnew),nmul,bc(idaprt(ikind,isa)),nbasdws(isa),0d0,   8d7s23
     $           bc(itmpt),nmul,'playr.tmpt4x')                         8d7s23
               end if                                                   6d26s24
              if(ikind.eq.1)then
               do i=0,nbasdws(isa)-1
                do j=0,nbasdws(isa)-1
                 do k=0,nmul-1
                  iad1=itnew+k+nmul*j
                  iad2=idaprt(ikind,isa)+j+nbasdws(isa)*i
                  iad3=itmpt+k+nmul*i
                  term=bc(iad1)*bc(iad2)
                 end do
                end do
               end do
              end if
              if(itmpt.le.igoxl.and.itmpt+nmul*nbasdws(isa).gt.igoxl)
     $             then
               write(6,*)('igoxld')
               idelta=igoxl-itmpt
               ii2=idelta/nmul
               ii1=idelta-nmul*ii2
               sum=0d0
               do i=0,nbasdws(isa)-1
                iad1=itnew+ii1+nmul*i
                if(ikind.eq.npartp)then                                 6d26s24
                 iad2=iarot(isa)+i+nbasdws(isa)*ii2                     6d26s24
                 value=darot(iad2)                                      6d26s24
                else                                                    6d26s24
                 iad2=idaprt(ikind,isa)+i+nbasdws(isa)*ii2
                 value=bc(iad2)                                         6d26s24
                end if                                                  6d26s24
                term=bc(iad1)*value                                     6d26s24
                sum=sum+term
                if(abs(term).gt.1d-10)write(6,*)bc(iad1),value,sum,     6d26s24
     $               iad1,iad2,ii1,i,ii2,nmul,nrow
               end do
              end if
               ibcoff=itnew                                             8d8s23
               jsa=isc                                                  8d8s23
               jsb=isd
               jsc=isb
               jsd=isa
               igotit=1
              end if                                                    8d8s23
             end if
             if(igotit.ne.0)then                                        8d8s23
              write(6,*)('to yield ')
              nrowu=nbasdws(jsa)*nbasdws(jsb)                           8d8s23
              ncolu=nbasdws(jsc)*nbasdws(jsd)                           8d8s23
              sz=0d0
              do i=0,nrowu*ncolu-1
               sz=sz+bc(itmpt+i)**2
              end do
              sz=sqrt(sz/dfloat(nrowu*ncolu))
              write(6,*)('sz: '),sz,nrowu*ncolu
              nc=noc(jsc)*nbasdws(jsd)
              do is1=1,nsdlk1                                           8d10s23
               if(isblk1(1,is1).eq.isblk1(2,is1))then                   8d10s23
                nrow1=(noc(isblk1(1,is1))*(noc(isblk1(1,is1))+1))/2     10d23s23
                nrow3=(nvirt(isblk1(1,is1))*(nvirt(isblk1(1,is1))+1))/2 8d10s23
                isw=0                                                   8d10s23
               else                                                     8d10s23
                nrow1=noc(isblk1(1,is1))*noc(isblk1(2,is1))             10d23s23
                nrow3=nvirt(isblk1(1,is1))*nvirt(isblk1(2,is1))         8d10s23
                isw=1                                                   8d10s23
               end if                                                   8d10s23
               ncolx=noc(isblk1(3,is1))*nvirt(isblk1(4,is1))            10d23s23
               ncolx3=irefo(isblk1(3,is1))*nvirt(isblk1(4,is1))         11d14s23
               j1x=jonex(is1)+nrow1*ncolx*(ikind-1)                     8d10s23
               j3xd=i3xd(is1)+nrow3*ncolx3*(ikind-1)                    11d14s23
               if(jsd.eq.isblk1(4,is1))then                             8d10s23
                if(jsa.eq.isblk1(1,is1).and.jsb.eq.isblk1(2,is1))then   8d10s23
                 do i4=0,nvirt(isblk1(4,is1))-1                         8d10s23
                  i4p=i4+noc(isblk1(4,is1))                             8d10s23
                  do i3=0,noc(isblk1(3,is1))-1                          10d23s23
                   icol=j1x+nrow1*(i3+noc(isblk1(3,is1))*i4)            10d23s23
                   do i2=0,noc(isblk1(2,is1))-1                         10d23s23
                    i1top=i2+isw*(noc(isblk1(1,is1))-1-i2)              10d23s23
                    do i1=0,i1top                                       8d10s23
                     irec=i1+noc(isblk1(1,is1))*i2                      10d23s23
                     itri=((i2*(i2+1))/2)+i1                            8d10s23
                     irec=itri+isw*(irec-itri)+icol                     8d10s23
                     iad=itmpt+i1+nbasdws(jsa)*(i2+nbasdws(jsb)*(i3     10d23s23
     $                    +nbasdws(jsc)*i4p))                           8d10s23
                     orig=bc(igoal)
                     bc(irec)=bc(irec)+bc(iad)                          8d10s23
                     if(abs(orig-bc(igoal)).gt.1d-12)write(6,*)
     $                    ('goala'),orig,bc(iad)*ass,bc(igoal),iad
                    end do                                              8d10s23
                   end do                                               8d10s23
                  end do                                                8d10s23
                 end do                                                 8d10s23
                 do i4=0,nvirt(isblk1(4,is1))-1                         8d10s23
                  i4p=i4+noc(isblk1(4,is1))                             8d10s23
                  do i3=0,irefo(isblk1(3,is1))-1                        11d14s23
                   i3p=i3+idoubo(isblk1(3,is1))                         11d14s23
                   icol=j3xd+nrow3*(i3+irefo(isblk1(3,is1))*i4)          11d14s23
                   do i2=0,nvirt(isblk1(2,is1))-1                       11d14s23
                    i2p=i2+noc(isblk1(2,is1))                           11d14s23
                    i1top=i2+isw*(nvirt(isblk1(1,is1))-1-i2)            11d14s23
                    do i1=0,i1top                                       8d10s23
                     i1p=i1+noc(isblk1(1,is1))                          11d14s23
                     irec=i1+nvirt(isblk1(1,is1))*i2                    11d14s23
                     itri=((i2*(i2+1))/2)+i1                            8d10s23
                     irec=itri+isw*(irec-itri)+icol                     8d10s23
                     iad=itmpt+i1p+nbasdws(jsa)*(i2p+nbasdws(jsb)*(i3p  11d14s23
     $                    +nbasdws(jsc)*i4p))                           8d10s23
                     orig=bc(igoal)
                     bc(irec)=bc(irec)+bc(iad)                          8d10s23
                     if(abs(orig-bc(igoal)).gt.1d-12)write(6,*)
     $                    ('goal3xa'),orig,bc(iad)*ass,bc(igoal),iad
                    end do                                              8d10s23
                   end do                                               8d10s23
                  end do                                                8d10s23
                 end do                                                 8d10s23
                else if(jsa.eq.isblk1(2,is1).and.                       8d10s23
     $               jsb.eq.isblk1(1,is1))then                          8d10s23
                 do i4=0,nvirt(isblk1(4,is1))-1                         8d10s23
                  i4p=i4+noc(isblk1(4,is1))                             8d10s23
                  do i3=0,noc(isblk1(3,is1))-1                          10d23s23
                   icol=j1x+nrow1*(i3+noc(isblk1(3,is1))*i4)            10d23s23
                   do i2=0,noc(isblk1(2,is1))-1                         10d23s23
                    i1top=i2+isw*(noc(isblk1(1,is1))-1-i2)              10d23s23
                    do i1=0,i1top                                       8d10s23
                     irec=i1+noc(isblk1(1,is1))*i2+icol                 10d24s23
                     iad=itmpt+i2+nbasdws(jsa)*(i1+nbasdws(jsb)*(i3     10d23s23
     $                    +nbasdws(jsc)*i4p))                           8d10s23
                     orig=bc(igoal)
                     bc(irec)=bc(irec)+bc(iad)                          8d10s23
                     if(abs(orig-bc(igoal)).gt.1d-12)write(6,*)
     $                    ('goalb'),orig,bc(iad)*ass,bc(igoal),iad
                    end do                                              8d10s23
                   end do                                               8d10s23
                  end do                                                8d10s23
                 end do                                                 8d10s23
                 do i4=0,nvirt(isblk1(4,is1))-1                         8d10s23
                  i4p=i4+noc(isblk1(4,is1))                             8d10s23
                  do i3=0,irefo(isblk1(3,is1))-1                        11d14s23
                   i3p=i3+idoubo(isblk1(3,is1))                         11d14s23
                   icol=j3xd+nrow3*(i3+irefo(isblk1(3,is1))*i4)          11d14s23
                   do i2=0,nvirt(isblk1(2,is1))-1                       11d14s23
                    i2p=i2+noc(isblk1(2,is1))                           11d14s23
                    i1top=i2+isw*(nvirt(isblk1(1,is1))-1-i2)            11d14s23
                    do i1=0,i1top                                       8d10s23
                     i1p=i1+noc(isblk1(1,is1))                          11d14s23
                     irec=i1+nvirt(isblk1(1,is1))*i2+icol               11d14s23
                     iad=itmpt+i2p+nbasdws(jsa)*(i1p+nbasdws(jsb)*(i3p  11d14s23
     $                    +nbasdws(jsc)*i4p))                           8d10s23
                     orig=bc(igoal)
                     bc(irec)=bc(irec)+bc(iad)                          8d10s23
                     if(abs(orig-bc(igoal)).gt.1d-12)write(6,*)
     $                    ('goal3b'),orig,bc(iad)*ass,bc(igoal),iad
                    end do                                              8d10s23
                   end do                                               8d10s23
                  end do                                                8d10s23
                 end do                                                 8d10s23
                end if                                                  8d10s23
               end if                                                   8d10s23
               if(jsd.eq.isblk1(3,is1))then                             8d10s23
                if(jsa.eq.isblk1(1,is1).and.jsb.eq.isblk1(2,is1))then   8d10s23
                 do i4=0,nvirt(isblk1(4,is1))-1                         8d10s23
                  i4p=i4+noc(isblk1(4,is1))                             8d10s23
                  do i3=0,noc(isblk1(3,is1))-1                          10d23s23
                   icol0=i3+noc(isblk1(3,is1))*i4                       10d30s23
                   icol=j1x+nrow1*icol0                                 10d30s23
                   do i2=0,noc(isblk1(2,is1))-1                         10d23s23
                    i1top=i2+isw*(noc(isblk1(1,is1))-1-i2)              10d23s23
                    do i1=0,i1top                                       8d10s23
                     irec=i1+noc(isblk1(1,is1))*i2                      10d23s23
                     itri=((i2*(i2+1))/2)+i1                            8d10s23
                     itri0=itri
                     irec=itri+isw*(irec-itri)+icol                     8d10s23
                     iad=itmpt+i1+nbasdws(jsa)*(i2+nbasdws(jsb)*(i4p    10d23s23
     $                    +nbasdws(jsc)*i3))                            10d23s23
                     orig=bc(igoal)
                     bc(irec)=bc(irec)+bc(iad)                          8d10s23
                     if(abs(orig-bc(igoal)).gt.1d-12)write(6,*)
     $                    ('goalc'),orig,bc(iad)*ass,bc(igoal),iad,
     $                    ikind,i1,i2,i3,i4,itri0,icol0,nrow1,irec-j1x,
     $                    j1x,irec-jonex(is1),is1
                    end do                                              8d10s23
                   end do                                               8d10s23
                  end do                                                8d10s23
                 end do                                                 8d10s23
                 do i4=0,nvirt(isblk1(4,is1))-1                         8d10s23
                  i4p=i4+noc(isblk1(4,is1))                             8d10s23
                  do i3=0,irefo(isblk1(3,is1))-1                        11d14s23
                   i3p=i3+idoubo(isblk1(3,is1))                         11d14s23
                   icol0=i3+irefo(isblk1(3,is1))*i4                     11d14s23
                   icol=j3xd+nrow3*icol0                                 10d30s23
                   do i2=0,nvirt(isblk1(2,is1))-1                       11d14s23
                    i2p=i2+noc(isblk1(2,is1))                           11d14s23
                    i1top=i2+isw*(nvirt(isblk1(1,is1))-1-i2)            11d14s23
                    do i1=0,i1top                                       8d10s23
                     i1p=i1+noc(isblk1(1,is1))                          11d14s23
                     irec=i1+nvirt(isblk1(1,is1))*i2                    11d14s23
                     itri=((i2*(i2+1))/2)+i1                            8d10s23
                     itri0=itri
                     irec=itri+isw*(irec-itri)+icol                     8d10s23
                     iad=itmpt+i1p+nbasdws(jsa)*(i2p+nbasdws(jsb)*(i4p  11d14s23
     $                    +nbasdws(jsc)*i3p))                           11d14s23
                     orig=bc(igoal)
                     bc(irec)=bc(irec)+bc(iad)                          8d10s23
                     if(abs(orig-bc(igoal)).gt.1d-12)write(6,*)
     $                    ('goal3c'),orig,bc(iad)*ass,bc(igoal),iad,
     $                    ikind,i1,i2,i3,i4,itri0,icol0,nrow1,irec-j1x,
     $                    j1x,irec-jonex(is1),is1
                    end do                                              8d10s23
                   end do                                               8d10s23
                  end do                                                8d10s23
                 end do                                                 8d10s23
                else if(jsa.eq.isblk1(2,is1).and.                       8d10s23
     $               jsb.eq.isblk1(1,is1))then                          8d10s23
                 do i4=0,nvirt(isblk1(4,is1))-1                         8d10s23
                  i4p=i4+noc(isblk1(4,is1))                             8d10s23
                  do i3=0,noc(isblk1(3,is1))-1                          10d23s23
                   icol=j1x+nrow1*(i3+noc(isblk1(3,is1))*i4)            10d23s23
                   do i2=0,noc(isblk1(2,is1))-1                         10d23s23
                    i1top=i2+isw*(noc(isblk1(1,is1))-1-i2)              10d23s23
                    do i1=0,i1top                                       8d10s23
                     irec=i1+noc(isblk1(1,is1))*i2+icol                 10d23s23
                     iad=itmpt+i2+nbasdws(jsa)*(i1+nbasdws(jsb)*(i4p    10d23s23
     $                    +nbasdws(jsc)*i3))                            10d23s23
                     orig=bc(igoal)
                     bc(irec)=bc(irec)+bc(iad)                          8d10s23
                     if(abs(orig-bc(igoal)).gt.1d-12)write(6,*)
     $                    ('goald'),orig,bc(iad)*ass,bc(igoal),iad
                    end do                                              8d10s23
                   end do                                               8d10s23
                  end do                                                8d10s23
                 end do                                                 8d10s23
                 do i4=0,nvirt(isblk1(4,is1))-1                         8d10s23
                  i4p=i4+noc(isblk1(4,is1))                             8d10s23
                  do i3=0,irefo(isblk1(3,is1))-1                        11d14s23
                   i3p=i3+idoubo(isblk1(3,is1))                         11d14s23
                   icol=j3xd+nrow3*(i3+irefo(isblk1(3,is1))*i4)         11d15s23
                   do i2=0,nvirt(isblk1(2,is1))-1                       11d14s23
                    i2p=i2+noc(isblk1(2,is1))                           11d14s23
                    i1top=i2+isw*(nvirt(isblk1(1,is1))-1-i2)            11d14s23
                    do i1=0,i1top                                       8d10s23
                     i1p=i1+noc(isblk1(1,is1))                          11d14s23
                     irec=i1+nvirt(isblk1(1,is1))*i2+icol               11d14s23
                     iad=itmpt+i2p+nbasdws(jsa)*(i1p+nbasdws(jsb)*(i4p  11d14s23
     $                    +nbasdws(jsc)*i3p))                           11d14s23
                     orig=bc(igoal)
                     bc(irec)=bc(irec)+bc(iad)                          8d10s23
                     if(abs(orig-bc(igoal)).gt.1d-12)write(6,*)
     $                    ('goal3d'),orig,bc(iad)*ass,bc(igoal),iad
                    end do                                              8d10s23
                   end do                                               8d10s23
                  end do                                                8d10s23
                 end do                                                 8d10s23
                end if                                                  8d10s23
               end if                                                   8d10s23
               if(jsd.eq.isblk1(2,is1))then                             8d10s23
                if(jsa.eq.isblk1(3,is1).and.jsb.eq.isblk1(4,is1))then   8d10s23
                 do i4=0,nvirt(isblk1(4,is1))-1                         8d10s23
                  i4p=i4+noc(isblk1(4,is1))                             8d10s23
                  do i3=0,noc(isblk1(3,is1))-1                          10d23s23
                   icol=j1x+nrow1*(i3+noc(isblk1(3,is1))*i4)            10d23s23
                   do i2=0,noc(isblk1(2,is1))-1                         10d23s23
                    i1top=i2+isw*(noc(isblk1(1,is1))-1-i2)              10d23s23
                    do i1=0,i1top                                       8d10s23
                     irec=i1+noc(isblk1(1,is1))*i2                      10d23s23
                     itri=((i2*(i2+1))/2)+i1                            8d10s23
                     irec=itri+isw*(irec-itri)+icol                     8d10s23
                     iad=itmpt+i3+nbasdws(jsa)*(i4p+nbasdws(jsb)*(i1    10d23s23
     $                    +nbasdws(jsc)*i2))                            10d23s23
                     orig=bc(igoal)
                     bc(irec)=bc(irec)+bc(iad)                          8d10s23
                     if(abs(orig-bc(igoal)).gt.1d-12)write(6,*)
     $                    ('goale'),orig,bc(iad)*ass,bc(igoal),iad,
     $                    i1,i2,i3,i4,irec-j1x,nrow1,isw
                    end do                                              8d10s23
                   end do                                               8d10s23
                  end do                                                8d10s23
                 end do                                                 8d10s23
                 do i4=0,nvirt(isblk1(4,is1))-1                         8d10s23
                  i4p=i4+noc(isblk1(4,is1))                             8d10s23
                  do i3=0,irefo(isblk1(3,is1))-1                          10d23s23
                   i3p=i3+idoubo(isblk1(3,is1))
                   icol=j3xd+nrow3*(i3+irefo(isblk1(3,is1))*i4)          11d14s23
                   do i2=0,nvirt(isblk1(2,is1))-1                       11d14s23
                    i2p=i2+noc(isblk1(2,is1))                           11d14s23
                    i1top=i2+isw*(nvirt(isblk1(1,is1))-1-i2)            11d14s23
                    do i1=0,i1top                                       8d10s23
                     i1p=i1+noc(isblk1(1,is1))                          11d14s23
                     irec=i1+nvirt(isblk1(1,is1))*i2                    11d14s23
                     itri=((i2*(i2+1))/2)+i1                            8d10s23
                     irec=itri+isw*(irec-itri)+icol                     8d10s23
                     iad=itmpt+i3p+nbasdws(jsa)*(i4p+nbasdws(jsb)*(i1p  11d14s23
     $                    +nbasdws(jsc)*i2p))                           11d14s23
                     orig=bc(igoal)
                     bc(irec)=bc(irec)+bc(iad)                          8d10s23
                     if(abs(orig-bc(igoal)).gt.1d-12)write(6,*)
     $                    ('goal3e'),orig,bc(iad)*ass,bc(igoal),iad,
     $                    i1,i2,i3,i4,irec-j1x,nrow1,isw
                    end do                                              8d10s23
                   end do                                               8d10s23
                  end do                                                8d10s23
                 end do                                                 8d10s23
                else if(jsa.eq.isblk1(4,is1).and.                       8d10s23
     $               jsb.eq.isblk1(3,is1))then                          8d10s23
                 do i4=0,nvirt(isblk1(4,is1))-1                         8d10s23
                  i4p=i4+noc(isblk1(4,is1))                             8d10s23
                  do i3=0,noc(isblk1(3,is1))-1                          10d23s23
                   icol=j1x+nrow1*(i3+noc(isblk1(3,is1))*i4)            10d23s23
                   do i2=0,noc(isblk1(2,is1))-1                         10d23s23
                    i1top=i2+isw*(noc(isblk1(1,is1))-1-i2)              10d23s23
                    do i1=0,i1top                                       8d10s23
                     irec=i1+noc(isblk1(1,is1))*i2                      10d23s23
                     irec0=irec
                     itri=((i2*(i2+1))/2)+i1                            8d10s23
                     irec=itri+isw*(irec-itri)+icol                     8d10s23
                     iad=itmpt+i4p+nbasdws(jsa)*(i3+nbasdws(jsb)*(i1    10d23s23
     $                    +nbasdws(jsc)*i2))                            10d23s23
                     orig=bc(igoal)
                     bc(irec)=bc(irec)+bc(iad)                          8d10s23
                     if(abs(orig-bc(igoal)).gt.1d-12)write(6,*)
     $                    ('goalf'),orig,bc(iad)*ass,bc(igoal),iad,
     $                    i1,i2,i3,i4,irec0
                    end do                                              8d10s23
                   end do                                               8d10s23
                  end do                                                8d10s23
                 end do                                                 8d10s23
                 do i4=0,nvirt(isblk1(4,is1))-1                         8d10s23
                  i4p=i4+noc(isblk1(4,is1))                             8d10s23
                  do i3=0,irefo(isblk1(3,is1))-1                          10d23s23
                   i3p=i3+idoubo(isblk1(3,is1))                         11d14s23
                   icol=j3xd+nrow3*(i3+irefo(isblk1(3,is1))*i4)          11d14s23
                   do i2=0,nvirt(isblk1(2,is1))-1                       11d14s23
                    i2p=i2+noc(isblk1(2,is1))                           11d14s23
                    i1top=i2+isw*(nvirt(isblk1(1,is1))-1-i2)            11d14s23
                    do i1=0,i1top                                       8d10s23
                     i1p=i1+noc(isblk1(1,is1))                          11d14s23
                     irec=i1+nvirt(isblk1(1,is1))*i2                    11d14s23
                     irec0=irec
                     itri=((i2*(i2+1))/2)+i1                            8d10s23
                     irec=itri+isw*(irec-itri)+icol                     8d10s23
                     iad=itmpt+i4p+nbasdws(jsa)*(i3p+nbasdws(jsb)*(i1p  11d14s23
     $                    +nbasdws(jsc)*i2p))                           11d14s23
                     orig=bc(igoal)
                     bc(irec)=bc(irec)+bc(iad)                          8d10s23
                     if(abs(orig-bc(igoal)).gt.1d-12)write(6,*)
     $                    ('goal3f'),orig,bc(iad)*ass,bc(igoal),iad,
     $                    i1,i2,i3,i4,irec0
                    end do                                              8d10s23
                   end do                                               8d10s23
                  end do                                                8d10s23
                 end do                                                 8d10s23
                end if                                                  8d10s23
               end if                                                   8d10s23
               if(jsd.eq.isblk1(1,is1))then                             8d10s23
                if(jsa.eq.isblk1(3,is1).and.jsb.eq.isblk1(4,is1))then   8d10s23
                 do i4=0,nvirt(isblk1(4,is1))-1                         8d10s23
                  i4p=i4+noc(isblk1(4,is1))                             8d10s23
                  do i3=0,noc(isblk1(3,is1))-1                          10d23s23
                   icol=j1x+nrow1*(i3+noc(isblk1(3,is1))*i4)            10d23s23
                   do i2=0,noc(isblk1(2,is1))-1                         10d23s23
                    i1top=i2+isw*(noc(isblk1(1,is1))-1-i2)              10d23s23
                    do i1=0,i1top                                       8d10s23
                     irec=i1+noc(isblk1(1,is1))*i2                      10d23s23
                     itri=((i2*(i2+1))/2)+i1                            8d10s23
                     irec=itri+isw*(irec-itri)+icol                     8d10s23
                     iad=itmpt+i3+nbasdws(jsa)*(i4p+nbasdws(jsb)*(i2    10d23s23
     $                    +nbasdws(jsc)*i1))                            10d23s23
                     orig=bc(igoal)
                     bc(irec)=bc(irec)+bc(iad)                          8d10s23
                     if(abs(orig-bc(igoal)).gt.1d-12)write(6,*)
     $                    ('goalg'),orig,bc(iad)*ass,bc(igoal),iad
                    end do                                              8d10s23
                   end do                                               8d10s23
                  end do                                                8d10s23
                 end do                                                 8d10s23
                 do i4=0,nvirt(isblk1(4,is1))-1                         8d10s23
                  i4p=i4+noc(isblk1(4,is1))                             8d10s23
                  do i3=0,irefo(isblk1(3,is1))-1                        11d14s23
                   i3p=i3+idoubo(isblk1(3,is1))                         11d14s23
                   icol=j3xd+nrow3*(i3+irefo(isblk1(3,is1))*i4)          11d14s23
                   do i2=0,nvirt(isblk1(2,is1))-1                       11d14s23
                    i2p=i2+noc(isblk1(2,is1))                           11d14s23
                    i1top=i2+isw*(nvirt(isblk1(1,is1))-1-i2)            11d14s23
                    do i1=0,i1top                                       8d10s23
                     i1p=i1+noc(isblk1(1,is1))                          11d14s23
                     irec=i1+nvirt(isblk1(1,is1))*i2                    11d14s23
                     itri=((i2*(i2+1))/2)+i1                            8d10s23
                     irec=itri+isw*(irec-itri)+icol                     8d10s23
                     iad=itmpt+i3p+nbasdws(jsa)*(i4p+nbasdws(jsb)*(i2p  11d14s23
     $                    +nbasdws(jsc)*i1p))                           11d14s23
                     orig=bc(igoal)
                     bc(irec)=bc(irec)+bc(iad)                          8d10s23
                     if(abs(orig-bc(igoal)).gt.1d-12)write(6,*)
     $                    ('goal3g'),orig,bc(iad)*ass,bc(igoal),iad
                    end do                                              8d10s23
                   end do                                               8d10s23
                  end do                                                8d10s23
                 end do                                                 8d10s23
                else if(jsa.eq.isblk1(4,is1).and.                       8d10s23
     $               jsb.eq.isblk1(3,is1))then                          8d10s23
                 do i4=0,nvirt(isblk1(4,is1))-1                         8d10s23
                  i4p=i4+noc(isblk1(4,is1))                             8d10s23
                  do i3=0,noc(isblk1(3,is1))-1                          10d23s23
                   icol=j1x+nrow1*(i3+noc(isblk1(3,is1))*i4)            10d23s23
                   do i2=0,noc(isblk1(2,is1))-1                         10d23s23
                    i1top=i2+isw*(noc(isblk1(1,is1))-1-i2)              10d23s23
                    do i1=0,i1top                                       8d10s23
                     irec=i1+noc(isblk1(1,is1))*i2                      10d23s23
                     itri=((i2*(i2+1))/2)+i1                            8d10s23
                     irec=itri+isw*(irec-itri)+icol                     8d10s23
                     iad=itmpt+i4p+nbasdws(jsa)*(i3+nbasdws(jsb)*(i2    10d23s23
     $                    +nbasdws(jsc)*i1))                            10d23s23
                     orig=bc(igoal)
                     bc(irec)=bc(irec)+bc(iad)                          8d10s23
                     if(abs(orig-bc(igoal)).gt.1d-12)write(6,*)
     $                    ('goalh'),orig,bc(iad)*ass,bc(igoal),iad
                    end do                                              8d10s23
                   end do                                               8d10s23
                  end do                                                8d10s23
                 end do                                                 8d10s23
                 do i4=0,nvirt(isblk1(4,is1))-1                         8d10s23
                  i4p=i4+noc(isblk1(4,is1))                             8d10s23
                  do i3=0,irefo(isblk1(3,is1))-1                          10d23s23
                   i3p=i3+idoubo(isblk1(3,is1))
                   icol=j3xd+nrow3*(i3+irefo(isblk1(3,is1))*i4)            10d23s23
                   do i2=0,nvirt(isblk1(2,is1))-1                         10d23s23
                    i2p=i2+noc(isblk1(2,is1))
                    i1top=i2+isw*(nvirt(isblk1(1,is1))-1-i2)              10d23s23
                    do i1=0,i1top                                       8d10s23
                     i1p=i1+noc(isblk1(1,is1))                          11d14s23
                     irec=i1+nvirt(isblk1(1,is1))*i2                      10d23s23
                     itri=((i2*(i2+1))/2)+i1                            8d10s23
                     irec=itri+isw*(irec-itri)+icol                     8d10s23
                     iad=itmpt+i4p+nbasdws(jsa)*(i3p+nbasdws(jsb)*(i2p    10d23s23
     $                    +nbasdws(jsc)*i1p))                            10d23s23
                     orig=bc(igoal)
                     bc(irec)=bc(irec)+bc(iad)                          8d10s23
                     if(abs(orig-bc(igoal)).gt.1d-12)write(6,*)
     $                    ('goal3h'),orig,bc(iad)*ass,bc(igoal),iad
                    end do                                              8d10s23
                   end do                                               8d10s23
                  end do                                                8d10s23
                 end do                                                 8d10s23
                end if                                                  8d10s23
               end if                                                   8d10s23
              end do                                                    8d10s23
              do is=1,nsdlkk                                               8d14s23
               nrowk=noc(isblkk(1,is))*noc(isblkk(2,is))                10d23s23
               ncolk=nvirt(isblkk(3,is))*nvirt(isblkk(4,is))            8d14s23
               if(jsd.eq.isblkk(4,is))then                              8d14s23
c                     ab cd
c     Kab^uv=(av|ub)=(ub|av)
c                     32 14
                if(isblkk(3,is).eq.jsa.and.isblkk(2,is).eq.jsb)then     8d14s23
                 kk=kmatd(is)+nrowk*ncolk*(ikind-1)                     8d14s23
                 do i4=0,nvirt(isblkk(4,is))-1                          8d14s23
                  i4p=i4+noc(isblkk(4,is))                              8d14s23
                  do i3=0,nvirt(isblkk(3,is))-1                         8d14s23
                   i3p=i3+noc(isblkk(3,is))                             8d14s23
                   jcol=i3+nvirt(isblkk(3,is))*i4                       8d14s23
                   icol=kk+nrowk*jcol                                   8d14s23
                   do i2=0,noc(isblkk(2,is))-1                          10d23s23
                    do i1=0,noc(isblkk(1,is))-1                         10d23s23
                     jrow=i1+noc(isblkk(1,is))*i2                       10d23s23
                     irow=icol+jrow                                     8d14s23
                     iad=itmpt+i3p+nbasdws(jsa)*(i2+nbasdws(jsb)        10d23s23
     $                    *(i1+nbasdws(jsc)*i4p))                       10d23s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                          8d14s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                    ('kgoala '),orig,bc(iad),bc(igoal),iad,
     $                    i1,i2,i3,i4,jrow,nrowk,jcol,(isblkk(j,is),
     $                    j=1,4),kmatd(is),kk-kmatd(is)
                    end do                                              8d14s23
                   end do                                               8d14s23
                  end do                                                8d14s23
                 end do                                                 8d14s23
c                     ab cd
c     Kab^uv=(av|ub)=(bu|av)
c                     23 14
                else if(isblkk(2,is).eq.jsa.and.isblkk(3,is).eq.jsb)then8d14s23
                 kk=kmatd(is)+nrowk*ncolk*(ikind-1)                     8d14s23
                 do i4=0,nvirt(isblkk(4,is))-1                          8d14s23
                  i4p=i4+noc(isblkk(4,is))                              8d14s23
                  do i3=0,nvirt(isblkk(3,is))-1                         8d14s23
                   i3p=i3+noc(isblkk(3,is))                             8d14s23
                   jcol=i3+nvirt(isblkk(3,is))*i4                       8d14s23
                   icol=kk+nrowk*jcol                                   8d14s23
                   do i2=0,noc(isblkk(2,is))-1                          10d23s23
                    do i1=0,noc(isblkk(1,is))-1                         10d23s23
                     jrow=i1+noc(isblkk(1,is))*i2                       10d23s23
                     irow=icol+jrow                                     8d14s23
                     iad=itmpt+i2+nbasdws(jsa)*(i3p+nbasdws(jsb)        10d23s23
     $                    *(i1+nbasdws(jsc)*i4p))                       10d23s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                          8d14s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                    ('kgoalb '),orig,bc(iad),bc(igoal),iad,
     $                    i1,i2,i3,i4,jrow,nrowk,jcol,(isblkk(j,is),
     $                    j=1,4),kmatd(is),kk-kmatd(is)
                    end do                                              8d14s23
                   end do                                               8d14s23
                  end do                                                8d14s23
                 end do                                                 8d14s23
                end if                                                  8d14s23
               end if                                                   8d14s23
               if(jsd.eq.isblkk(3,is))then                              8d14s23
c                     ab cd
c     Kab^uv=(av|ub)=(av|bu)
c                     14 23
                if(isblkk(1,is).eq.jsa.and.isblkk(4,is).eq.jsb)then     8d14s23
                 kk=kmatd(is)+nrowk*ncolk*(ikind-1)                     8d14s23
                 do i4=0,nvirt(isblkk(4,is))-1                          8d14s23
                  i4p=i4+noc(isblkk(4,is))                              8d14s23
                  do i3=0,nvirt(isblkk(3,is))-1                         8d14s23
                   i3p=i3+noc(isblkk(3,is))                             8d14s23
                   jcol=i3+nvirt(isblkk(3,is))*i4                       8d14s23
                   icol=kk+nrowk*jcol                                   8d14s23
                   do i2=0,noc(isblkk(2,is))-1                          10d23s23
                    do i1=0,noc(isblkk(1,is))-1                         10d23s23
                     jrow=i1+noc(isblkk(1,is))*i2                       10d23s23
                     irow=icol+jrow                                     8d14s23
                     iad=itmpt+i1+nbasdws(jsa)*(i4p+nbasdws(jsb)        10d23s23
     $                    *(i2+nbasdws(jsc)*i3p))                       10d23s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                          8d14s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                    ('kgoalc '),orig,bc(iad),bc(igoal),iad,
     $                    i1,i2,i3,i4,jrow,nrowk,jcol,(isblkk(j,is),
     $                    j=1,4),kmatd(is),kk-kmatd(is)
                    end do                                              8d14s23
                   end do                                               8d14s23
                  end do                                                8d14s23
                 end do                                                 8d14s23
c                     ab cd
c     Kab^uv=(av|ub)=(va|bu)
c                     41 23
                else if(isblkk(1,is).eq.jsb.and.isblkk(4,is).eq.jsa)then8d14s23
                 kk=kmatd(is)+nrowk*ncolk*(ikind-1)                     8d14s23
                 do i4=0,nvirt(isblkk(4,is))-1                          8d14s23
                  i4p=i4+noc(isblkk(4,is))                              8d14s23
                  do i3=0,nvirt(isblkk(3,is))-1                         8d14s23
                   i3p=i3+noc(isblkk(3,is))                             8d14s23
                   jcol=i3+nvirt(isblkk(3,is))*i4                       8d14s23
                   icol=kk+nrowk*jcol                                   8d14s23
                   do i2=0,noc(isblkk(2,is))-1                          10d23s23
                    do i1=0,noc(isblkk(1,is))-1                         10d23s23
                     jrow=i1+noc(isblkk(1,is))*i2                       10d23s23
                     irow=icol+jrow                                     8d14s23
                     iad=itmpt+i4p+nbasdws(jsa)*(i1+nbasdws(jsb)        10d23s23
     $                    *(i2+nbasdws(jsc)*i3p))                       10d23s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                          8d14s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                    ('kgoald '),orig,bc(iad),bc(igoal),iad,
     $                    i1,i2,i3,i4,jrow,nrowk,jcol,(isblkk(j,is),
     $                    j=1,4),kmatd(is),kk-kmatd(is)
                    end do                                              8d14s23
                   end do                                               8d14s23
                  end do                                                8d14s23
                 end do                                                 8d14s23
                end if                                                  8d14s23
               end if                                                   8d14s23
               if(jsd.eq.isblkk(2,is))then                              8d14s23
c             ab cd
c     Kab^uv=(av|ub)
c             14 32
                if(isblkk(1,is).eq.jsa.and.isblkk(4,is).eq.jsb)then     8d14s23
                 kk=kmatd(is)+nrowk*ncolk*(ikind-1)                     8d14s23
                 do i4=0,nvirt(isblkk(4,is))-1                          8d14s23
                  i4p=i4+noc(isblkk(4,is))                              8d14s23
                  do i3=0,nvirt(isblkk(3,is))-1                         8d14s23
                   i3p=i3+noc(isblkk(3,is))                             8d14s23
                   jcol=i3+nvirt(isblkk(3,is))*i4                       8d14s23
                   icol=kk+nrowk*jcol                                   8d14s23
                   do i2=0,noc(isblkk(2,is))-1                          10d23s23
                    do i1=0,noc(isblkk(1,is))-1                         10d23s23
                     jrow=i1+noc(isblkk(1,is))*i2                       10d23s23
                     irow=icol+jrow                                     8d14s23
                     iad=itmpt+i1+nbasdws(jsa)*(i4p+nbasdws(jsb)        10d23s23
     $                    *(i3p+nbasdws(jsc)*i2))                       10d23s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                          8d14s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                    ('kgoale '),orig,bc(iad),bc(igoal),
     $                    bc(igoal)*ass,iad,
     $                    i1,i2,i3,i4,jrow,nrowk,jcol,(isblkk(j,is),
     $                    j=1,4),kmatd(is),kk-kmatd(is)
                    end do                                              8d14s23
                   end do                                               8d14s23
                  end do                                                8d14s23
                 end do                                                 8d14s23
c             ab cd
c     Kab^uv=(va|ub)
c             41 32
                else if(isblkk(1,is).eq.jsb.and.isblkk(4,is).eq.jsa)then8d14s23
                 kk=kmatd(is)+nrowk*ncolk*(ikind-1)                     8d14s23
                 do i4=0,nvirt(isblkk(4,is))-1                          8d14s23
                  i4p=i4+noc(isblkk(4,is))                              8d14s23
                  do i3=0,nvirt(isblkk(3,is))-1                         8d14s23
                   i3p=i3+noc(isblkk(3,is))                             8d14s23
                   jcol=i3+nvirt(isblkk(3,is))*i4                       8d14s23
                   icol=kk+nrowk*jcol                                   8d14s23
                   do i2=0,noc(isblkk(2,is))-1                          10d23s23
                    do i1=0,noc(isblkk(1,is))-1                         10d23s23
                     jrow=i1+noc(isblkk(1,is))*i2                       10d23s23
                     irow=icol+jrow                                     8d14s23
                     iad=itmpt+i4p+nbasdws(jsa)*(i1+nbasdws(jsb)        10d23s23
     $                    *(i3p+nbasdws(jsc)*i2))                       10d23s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                          8d14s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                    ('kgoalf '),orig,bc(iad),bc(igoal),iad,
     $                    i1,i2,i3,i4,jrow,nrowk,jcol,(isblkk(j,is),
     $                    j=1,4),kmatd(is),kk-kmatd(is)
                    end do                                              8d14s23
                   end do                                               8d14s23
                  end do                                                8d14s23
                 end do                                                 8d14s23
                end if                                                  8d14s23
               end if                                                   8d14s23
               if(jsd.eq.isblkk(1,is))then                              8d14s23
c                             ab cd
c     Kab^uv=(av|ub)=(ub|av)=(ub|va)
c                             32 41
                if(isblkk(3,is).eq.jsa.and.isblkk(2,is).eq.jsb)then     8d14s23
                 kk=kmatd(is)+nrowk*ncolk*(ikind-1)                     8d14s23
                 do i4=0,nvirt(isblkk(4,is))-1                          8d14s23
                  i4p=i4+noc(isblkk(4,is))                              8d14s23
                  do i3=0,nvirt(isblkk(3,is))-1                         8d14s23
                   i3p=i3+noc(isblkk(3,is))                             8d14s23
                   jcol=i3+nvirt(isblkk(3,is))*i4                       8d14s23
                   icol=kk+nrowk*jcol                                   8d14s23
                   do i2=0,noc(isblkk(2,is))-1                          10d23s23
                    do i1=0,noc(isblkk(1,is))-1                         10d23s23
                     jrow=i1+noc(isblkk(1,is))*i2                       10d23s23
                     irow=icol+jrow
                     iad=itmpt+i3p+nbasdws(jsa)*(i2+nbasdws(jsb)        10d23s23
     $                    *(i4p+nbasdws(jsc)*i1))                       10d23s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                          8d14s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                    ('kgoalg '),orig,bc(iad),bc(igoal),iad,
     $                    i1,i2,i3,i4,jrow,nrowk,jcol,(isblkk(j,is),
     $                    j=1,4),kmatd(is),kk-kmatd(is)
                    end do                                              8d14s23
                   end do                                               8d14s23
                  end do                                                8d14s23
                 end do                                                 8d14s23
c                             ab cd
c     Kab^uv=(av|ub)=(ub|av)=(bu|va)
c                             23 41
                else if(isblkk(2,is).eq.jsa.and.isblkk(3,is).eq.jsb)then8d14s23
                 kk=kmatd(is)+nrowk*ncolk*(ikind-1)                     8d14s23
                 do i4=0,nvirt(isblkk(4,is))-1                          8d14s23
                  i4p=i4+noc(isblkk(4,is))                              8d14s23
                  do i3=0,nvirt(isblkk(3,is))-1                         8d14s23
                   i3p=i3+noc(isblkk(3,is))                             8d14s23
                   jcol=i3+nvirt(isblkk(3,is))*i4                       8d14s23
                   icol=kk+nrowk*jcol                                   8d14s23
                   do i2=0,noc(isblkk(2,is))-1                          10d23s23
                    do i1=0,noc(isblkk(1,is))-1                         10d23s23
                     jrow=i1+noc(isblkk(1,is))*i2                       10d23s23
                     irow=icol+jrow                                     8d14s23
                     iad=itmpt+i2+nbasdws(jsa)*(i3p+nbasdws(jsb)        10d23s23
     $                    *(i4p+nbasdws(jsc)*i1))                       10d23s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                          8d14s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                    ('kgoalh '),orig,bc(iad),bc(igoal),iad,
     $                    i1,i2,i3,i4,jrow,nrowk,jcol,(isblkk(j,is),
     $                    j=1,4),kmatd(is),kk-kmatd(is)
                    end do                                              8d14s23
                   end do                                               8d14s23
                  end do                                                8d14s23
                 end do                                                 8d14s23
                end if                                                  8d14s23
               end if                                                   8d14s23
              end do                                                       8d14s23
              do is4=1,n4x                                                 11d17s23
               if(isblk4x(1,is4).eq.isblk4x(2,is4))then                 11d21s23
                mab=(nvirt(isblk4x(1,is4))*(nvirt(isblk4x(1,is4))+1))/2 11d17s23
                mcd=(nvirt(isblk4x(3,is4))*(nvirt(isblk4x(3,is4))+1))/2 11d17s23
                iswab=0                                                 11d17s23
               else                                                     11d17s23
                mab=nvirt(isblk4x(1,is4))*nvirt(isblk4x(2,is4))         11d17s23
                mcd=nvirt(isblk4x(3,is4))*nvirt(isblk4x(4,is4))         11d17s23
                iswab=1                                                 11d17s23
               end if                                                           11d17s23
               itriab=((isblk4x(1,is4)*(isblk4x(1,is4)-1))/2)           11d17s23
     $              +isblk4x(2,is4)                                     11d17s23
               itricd=((isblk4x(3,is4)*(isblk4x(3,is4)-1))/2)           11d17s23
     $              +isblk4x(4,is4)                                     11d17s23
               if(itriab.eq.itricd)then                                 11d17s23
                iswac=0                                                 11d17s23
                nall=(mab*(mab+1))/2                                    11d17s23
               else                                                     11d17s23
                iswac=1                                                 11d17s23
                nall=mab*mcd                                            11d17s23
               end if                                                   11d17s23
                jj=i4xd(is4,ikind)                                      11d21s23
               if(jsd.eq.isblk4x(4,is4))then                            11d17s23
                if(isblk4x(1,is4).eq.jsa.and.isblk4x(2,is4).eq.jsb)then 11d17s23
                 i34=0                                                  11d17s23
                 do i4=0,nvirt(isblk4x(4,is4))-1                        11d17s23
                  i3top=i4+iswab*(nvirt(isblk4x(3,is4))-1-i4)           11d17s23
                  i4p=i4+noc(isblk4x(4,is4))                            11d17s23
                  do i3=0,i3top                                         11d17s23
                   i3p=i3+noc(isblk4x(3,is4))                           11d17s23
                   i12top=i34+iswac*(mab-1-i34)                         11d21s23
                   irec=i3+nvirt(isblk4x(3,is4))*i4                     11d17s23
                   itri=((i4*(i4+1))/2)+i3                              11d17s23
                   icol=itri+iswab*(irec-itri)                          11d17s23
                   i12=0                                                11d17s23
                   do i2=0,nvirt(isblk4x(2,is4))-1                      11d17s23
                    i1top=i2+iswab*(nvirt(isblk4x(1,is4))-1-i2)         11d17s23
                    i2p=i2+noc(isblk4x(2,is4))                          11d17s23
                    do i1=0,i1top                                       11d17s23
                     if(i12.le.i12top)then                              11d17s23
                      irec=i1+nvirt(isblk4x(1,is4))*i2                  11d17s23
                      itri=((i2*(i2+1))/2)+i1                           11d17s23
                      irow=itri+iswab*(irec-itri)                       11d17s23
                      irec=irow+mab*icol                                11d17s23
                      itri=((icol*(icol+1))/2)+irow                     11d17s23
                      itri=itri+iswac*(irec-itri)+jj                    11d17s23
                      i1p=i1+noc(isblk4x(1,is4))                          11d17s23
                      iad=itmpt+i1p+nbasdws(jsa)*(i2p+nbasdws(jsb)*(i3p 11d17s23
     $                    +nbasdws(jsc)*i4p))                           11d17s23
                orig=bc(igoal)
                      bc(itri)=bc(itri)+bc(iad)                         11d17s23
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('goal4ih '),orig,bc(iad),bc(igoal),
     $                     bc(iad)*ass
                     end if                                             11d17s23
                     i12=i12+1                                          11d17s23
                    end do                                              11d17s23
                   end do                                               11d17s23
                   i34=i34+1                                            11d17s23
                  end do                                                11d17s23
                 end do                                                 11d17s23
                else if(isblk4x(2,is4).eq.jsa.and.                      11d17s23
     $               isblk4x(1,is4).eq.jsb)then                         11d17s23
                 i34=0                                                  11d17s23
                 do i4=0,nvirt(isblk4x(4,is4))-1                        11d17s23
                  i3top=i4+iswab*(nvirt(isblk4x(3,is4))-1-i4)           11d17s23
                  i4p=i4+noc(isblk4x(4,is4))                            11d17s23
                  do i3=0,i3top                                         11d17s23
                   i3p=i3+noc(isblk4x(3,is4))                           11d17s23
                   i12top=i34+iswac*(mab-1-i34)                         11d21s23
                   irec=i3+nvirt(isblk4x(3,is4))*i4                     11d17s23
                   itri=((i4*(i4+1))/2)+i3                              11d17s23
                   icol=itri+iswab*(irec-itri)                          11d17s23
                   i12=0                                                11d17s23
                   do i2=0,nvirt(isblk4x(2,is4))-1                      11d17s23
                    i1top=i2+iswab*(nvirt(isblk4x(1,is4))-1-i2)         11d17s23
                    i2p=i2+noc(isblk4x(2,is4))                          11d17s23
                    do i1=0,i1top                                       11d17s23
                     if(i12.le.i12top)then                              11d17s23
                      irec=i1+nvirt(isblk4x(1,is4))*i2                  11d17s23
                      itri=((i2*(i2+1))/2)+i1                           11d17s23
                      irow=itri+iswab*(irec-itri)                       11d17s23
                      irec=irow+mab*icol                                11d17s23
                      itri=((icol*(icol+1))/2)+irow                     11d17s23
                      itri=itri+iswac*(irec-itri)+jj                    11d17s23
                      i1p=i1+noc(isblk4x(1,is4))                          11d17s23
                      iad=itmpt+i2p+nbasdws(jsa)*(i1p+nbasdws(jsb)*(i3p 11d17s23
     $                    +nbasdws(jsc)*i4p))                           11d17s23
                orig=bc(igoal)
                      bc(itri)=bc(itri)+bc(iad)                         11d17s23
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('goal4ig '),orig,bc(iad),bc(igoal),
     $                     bc(iad)*ass
                     end if                                             11d17s23
                     i12=i12+1                                          11d17s23
                    end do                                              11d17s23
                   end do                                               11d17s23
                   i34=i34+1                                            11d17s23
                  end do                                                11d17s23
                 end do                                                 11d17s23
                end if                                                  11d17s23
               end if                                                   11d17s23
               if(jsd.eq.isblk4x(3,is4))then                            11d17s23
                if(isblk4x(1,is4).eq.jsa.and.isblk4x(2,is4).eq.jsb)then 11d17s23
                 i34=0                                                  11d17s23
                 do i4=0,nvirt(isblk4x(4,is4))-1                        11d17s23
                  i3top=i4+iswab*(nvirt(isblk4x(3,is4))-1-i4)           11d17s23
                  i4p=i4+noc(isblk4x(4,is4))                            11d17s23
                  do i3=0,i3top                                         11d17s23
                   i3p=i3+noc(isblk4x(3,is4))                           11d17s23
                   i12top=i34+iswac*(mab-1-i34)                         11d21s23
                   irec=i3+nvirt(isblk4x(3,is4))*i4                     11d17s23
                   itri=((i4*(i4+1))/2)+i3                              11d17s23
                   icol=itri+iswab*(irec-itri)                          11d17s23
                   icol0=icol
                   i12=0                                                11d17s23
                   do i2=0,nvirt(isblk4x(2,is4))-1                      11d17s23
                    i1top=i2+iswab*(nvirt(isblk4x(1,is4))-1-i2)         11d17s23
                    i2p=i2+noc(isblk4x(2,is4))                          11d17s23
                    do i1=0,i1top                                       11d17s23
                     if(i12.le.i12top)then                              11d17s23
                      irec=i1+nvirt(isblk4x(1,is4))*i2                  11d17s23
                      itri=((i2*(i2+1))/2)+i1                           11d17s23
                      irow=itri+iswab*(irec-itri)                       11d17s23
                      irow0=irow
                      irec=irow+mab*icol                                11d17s23
                      irec0=irec
                      itri=((icol*(icol+1))/2)+irow                     11d17s23
                      itri0=itri
                      itri=itri+iswac*(irec-itri)+jj                    11d17s23
                      i1p=i1+noc(isblk4x(1,is4))                          11d17s23
                      iad=itmpt+i1p+nbasdws(jsa)*(i2p+nbasdws(jsb)*(i4p 11d17s23
     $                    +nbasdws(jsc)*i3p))                           11d17s23
                orig=bc(igoal)
                      bc(itri)=bc(itri)+bc(iad)                         11d17s23
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('goal4if '),orig,bc(iad),bc(igoal),
     $                     bc(iad)*ass,ikind,iad,(itri-i4xd(is4,iiiii),
     $                     iiiii=1,3),(i4xd(is4,iiiii),iiiii=1,3),
     $                     itri,is4,irow0,icol0,irec0,itri0,iswac
                     end if                                             11d17s23
                     i12=i12+1                                          11d17s23
                    end do                                              11d17s23
                   end do                                               11d17s23
                   i34=i34+1                                            11d17s23
                  end do                                                11d17s23
                 end do                                                 11d17s23
                else if(isblk4x(2,is4).eq.jsa.and.                      11d17s23
     $               isblk4x(1,is4).eq.jsb)then                         11d17s23
                 i34=0                                                  11d17s23
                 do i4=0,nvirt(isblk4x(4,is4))-1                        11d17s23
                  i3top=i4+iswab*(nvirt(isblk4x(3,is4))-1-i4)           11d17s23
                  i4p=i4+noc(isblk4x(4,is4))                            11d17s23
                  do i3=0,i3top                                         11d17s23
                   i3p=i3+noc(isblk4x(3,is4))                           11d17s23
                   i12top=i34+iswac*(mab-1-i34)                         11d21s23
                   irec=i3+nvirt(isblk4x(3,is4))*i4                     11d17s23
                   itri=((i4*(i4+1))/2)+i3                              11d17s23
                   icol=itri+iswab*(irec-itri)                          11d17s23
                   i12=0                                                11d17s23
                   do i2=0,nvirt(isblk4x(2,is4))-1                      11d17s23
                    i1top=i2+iswab*(nvirt(isblk4x(1,is4))-1-i2)         11d17s23
                    i2p=i2+noc(isblk4x(2,is4))                          11d17s23
                    do i1=0,i1top                                       11d17s23
                     if(i12.le.i12top)then                              11d17s23
                      irec=i1+nvirt(isblk4x(1,is4))*i2                  11d17s23
                      itri=((i2*(i2+1))/2)+i1                           11d17s23
                      irow=itri+iswab*(irec-itri)                       11d17s23
                      irec=irow+mab*icol                                11d17s23
                      itri=((icol*(icol+1))/2)+irow                     11d17s23
                      itri=itri+iswac*(irec-itri)+jj                    11d17s23
                      i1p=i1+noc(isblk4x(1,is4))                          11d17s23
                      iad=itmpt+i2p+nbasdws(jsa)*(i1p+nbasdws(jsb)*(i4p 11d17s23
     $                    +nbasdws(jsc)*i3p))                           11d17s23
                orig=bc(igoal)
                      bc(itri)=bc(itri)+bc(iad)                         11d17s23
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('goal4ie '),orig,bc(iad),bc(igoal),
     $                     bc(iad)*ass
                     end if                                             11d17s23
                     i12=i12+1                                          11d17s23
                    end do                                              11d17s23
                   end do                                               11d17s23
                   i34=i34+1                                            11d17s23
                  end do                                                11d17s23
                 end do                                                 11d17s23
                end if                                                  11d17s23
               end if                                                   11d17s23
               if(jsd.eq.isblk4x(2,is4))then                            11d17s23
                if(isblk4x(3,is4).eq.jsa.and.isblk4x(4,is4).eq.jsb)then 11d17s23
                 i34=0                                                  11d17s23
                 do i4=0,nvirt(isblk4x(4,is4))-1                        11d17s23
                  i3top=i4+iswab*(nvirt(isblk4x(3,is4))-1-i4)           11d17s23
                  i4p=i4+noc(isblk4x(4,is4))                            11d17s23
                  do i3=0,i3top                                         11d17s23
                   i3p=i3+noc(isblk4x(3,is4))                           11d17s23
                   i12top=i34+iswac*(mab-1-i34)                         11d21s23
                   irec=i3+nvirt(isblk4x(3,is4))*i4                     11d17s23
                   itri=((i4*(i4+1))/2)+i3                              11d17s23
                   icol=itri+iswab*(irec-itri)                          11d17s23
                   i12=0                                                11d17s23
                   do i2=0,nvirt(isblk4x(2,is4))-1                      11d17s23
                    i1top=i2+iswab*(nvirt(isblk4x(1,is4))-1-i2)         11d17s23
                    i2p=i2+noc(isblk4x(2,is4))                          11d17s23
                    do i1=0,i1top                                       11d17s23
                     if(i12.le.i12top)then                              11d17s23
                      irec=i1+nvirt(isblk4x(1,is4))*i2                  11d17s23
                      itri=((i2*(i2+1))/2)+i1                           11d17s23
                      irow=itri+iswab*(irec-itri)                       11d17s23
                      irec=irow+mab*icol                                11d17s23
                      itri=((icol*(icol+1))/2)+irow                     11d17s23
                      itri=itri+iswac*(irec-itri)+jj                    11d17s23
                      i1p=i1+noc(isblk4x(1,is4))                          11d17s23
                      iad=itmpt+i3p+nbasdws(jsa)*(i4p+nbasdws(jsb)*(i1p 11d17s23
     $                    +nbasdws(jsc)*i2p))                           11d17s23
                orig=bc(igoal)
                      bc(itri)=bc(itri)+bc(iad)                         11d17s23
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('goal4ic '),orig,bc(iad),bc(igoal),
     $                     bc(iad)*ass
                     end if                                             11d17s23
                     i12=i12+1                                          11d17s23
                    end do                                              11d17s23
                   end do                                               11d17s23
                   i34=i34+1                                            11d17s23
                  end do                                                11d17s23
                 end do                                                 11d17s23
                else if(isblk4x(4,is4).eq.jsa.and.                      11d17s23
     $               isblk4x(3,is4).eq.jsb)then                         11d17s23
                 i34=0                                                  11d17s23
                 do i4=0,nvirt(isblk4x(4,is4))-1                        11d17s23
                  i3top=i4+iswab*(nvirt(isblk4x(3,is4))-1-i4)           11d17s23
                  i4p=i4+noc(isblk4x(4,is4))                            11d17s23
                  do i3=0,i3top                                         11d17s23
                   i3p=i3+noc(isblk4x(3,is4))                           11d17s23
                   i12top=i34+iswac*(mab-1-i34)                         11d21s23
                   irec=i3+nvirt(isblk4x(3,is4))*i4                     11d17s23
                   itri=((i4*(i4+1))/2)+i3                              11d17s23
                   icol=itri+iswab*(irec-itri)                          11d17s23
                   i12=0                                                11d17s23
                   do i2=0,nvirt(isblk4x(2,is4))-1                      11d17s23
                    i1top=i2+iswab*(nvirt(isblk4x(1,is4))-1-i2)         11d17s23
                    i2p=i2+noc(isblk4x(2,is4))                          11d17s23
                    do i1=0,i1top                                       11d17s23
                     if(i12.le.i12top)then                              11d17s23
                      irec=i1+nvirt(isblk4x(1,is4))*i2                  11d17s23
                      itri=((i2*(i2+1))/2)+i1                           11d17s23
                      irow=itri+iswab*(irec-itri)                       11d17s23
                      irec=irow+mab*icol                                11d17s23
                      itri=((icol*(icol+1))/2)+irow                     11d17s23
                      itri=itri+iswac*(irec-itri)+jj                    11d17s23
                      i1p=i1+noc(isblk4x(1,is4))                          11d17s23
                      iad=itmpt+i4p+nbasdws(jsa)*(i3p+nbasdws(jsb)*(i1p 11d17s23
     $                    +nbasdws(jsc)*i2p))                           11d17s23
                orig=bc(igoal)
                      bc(itri)=bc(itri)+bc(iad)                         11d17s23
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('goal4idd '),orig,bc(iad),bc(igoal),
     $                     bc(iad)*ass
                     end if                                             11d17s23
                     i12=i12+1                                          11d17s23
                    end do                                              11d17s23
                   end do                                               11d17s23
                   i34=i34+1                                            11d17s23
                  end do                                                11d17s23
                 end do                                                 11d17s23
                end if                                                  11d17s23
               end if                                                   11d17s23
               if(jsd.eq.isblk4x(1,is4))then                            11d17s23
                if(isblk4x(3,is4).eq.jsa.and.isblk4x(4,is4).eq.jsb)then 11d17s23
                 i34=0                                                  11d17s23
                 do i4=0,nvirt(isblk4x(4,is4))-1                        11d17s23
                  i3top=i4+iswab*(nvirt(isblk4x(3,is4))-1-i4)           11d17s23
                  i4p=i4+noc(isblk4x(4,is4))                            11d17s23
                  do i3=0,i3top                                         11d17s23
                   i3p=i3+noc(isblk4x(3,is4))                           11d17s23
                   i12top=i34+iswac*(mab-1-i34)                         11d21s23
                   irec=i3+nvirt(isblk4x(3,is4))*i4                     11d17s23
                   itri=((i4*(i4+1))/2)+i3                              11d17s23
                   icol=itri+iswab*(irec-itri)                          11d17s23
                   i12=0                                                11d17s23
                   do i2=0,nvirt(isblk4x(2,is4))-1                      11d17s23
                    i1top=i2+iswab*(nvirt(isblk4x(1,is4))-1-i2)         11d17s23
                    i2p=i2+noc(isblk4x(2,is4))                          11d17s23
                    do i1=0,i1top                                       11d17s23
                     if(i12.le.i12top)then                              11d17s23
                      irec=i1+nvirt(isblk4x(1,is4))*i2                  11d17s23
                      itri=((i2*(i2+1))/2)+i1                           11d17s23
                      irow=itri+iswab*(irec-itri)                       11d17s23
                      irec=irow+mab*icol                                11d17s23
                      itri=((icol*(icol+1))/2)+irow                     11d17s23
                      itri=itri+iswac*(irec-itri)+jj                    11d17s23
                      i1p=i1+noc(isblk4x(1,is4))                          11d17s23
                      iad=itmpt+i3p+nbasdws(jsa)*(i4p+nbasdws(jsb)*(i2p 11d21s23
     $                    +nbasdws(jsc)*i1p))                           11d21s23
                      orig=bc(igoal)
                      bc(itri)=bc(itri)+bc(iad)                         11d17s23
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('goal4ia '),orig,bc(iad),bc(igoal),
     $                     bc(iad)*ass
                     end if                                             11d17s23
                     i12=i12+1                                          11d17s23
                    end do                                              11d17s23
                   end do                                               11d17s23
                   i34=i34+1                                            11d17s23
                  end do                                                11d17s23
                 end do                                                 11d17s23
                else if(isblk4x(4,is4).eq.jsa.and.                      11d17s23
     $               isblk4x(3,is4).eq.jsb)then                         11d17s23
                 i34=0                                                  11d17s23
                 do i4=0,nvirt(isblk4x(4,is4))-1                        11d17s23
                  i3top=i4+iswab*(nvirt(isblk4x(3,is4))-1-i4)           11d17s23
                  i4p=i4+noc(isblk4x(4,is4))                            11d17s23
                  do i3=0,i3top                                         11d17s23
                   i3p=i3+noc(isblk4x(3,is4))                           11d17s23
                   i12top=i34+iswac*(mab-1-i34)                         11d21s23
                   irec=i3+nvirt(isblk4x(3,is4))*i4                     11d17s23
                   itri=((i4*(i4+1))/2)+i3                              11d17s23
                   icol=itri+iswab*(irec-itri)                          11d17s23
                   i12=0                                                11d17s23
                   do i2=0,nvirt(isblk4x(2,is4))-1                      11d17s23
                    i1top=i2+iswab*(nvirt(isblk4x(1,is4))-1-i2)         11d17s23
                    i2p=i2+noc(isblk4x(2,is4))                          11d17s23
                    do i1=0,i1top                                       11d17s23
                     if(i12.le.i12top)then                              11d17s23
                      irec=i1+nvirt(isblk4x(1,is4))*i2                  11d17s23
                      itri=((i2*(i2+1))/2)+i1                           11d17s23
                      irow=itri+iswab*(irec-itri)                       11d17s23
                      irec=irow+mab*icol                                11d17s23
                      itri=((icol*(icol+1))/2)+irow                     11d17s23
                      itri=itri+iswac*(irec-itri)+jj                    11d17s23
                      i1p=i1+noc(isblk4x(1,is4))                          11d17s23
                      iad=itmpt+i4p+nbasdws(jsa)*(i3p+nbasdws(jsb)*(i2p 11d17s23
     $                    +nbasdws(jsc)*i1p))                           11d17s23
                orig=bc(igoal)
                      bc(itri)=bc(itri)+bc(iad)                         11d17s23
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('goal4ib '),orig,bc(iad),bc(igoal),
     $                     bc(iad)*ass,itri-i4xd(is4,ikind),itri-jj
                     end if                                             11d17s23
                     i12=i12+1                                          11d17s23
                    end do                                              11d17s23
                   end do                                               11d17s23
                   i34=i34+1                                            11d17s23
                  end do                                                11d17s23
                 end do                                                 11d17s23
                end if                                                  11d17s23
               end if                                                   11d17s23
              end do                                                       11d17s23
              do is=1,nsdlk                                               8d7s23
               if(isblk(1,is).eq.isblk(2,is))then                         8d8s23
                nrowj=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2           8d7s23
                ncolj=(noc(isblk(3,is))*(noc(isblk(3,is))+1))/2           8d7s23
                nrowjd=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2    8d11s23
                isw=0                                                     8d7s23
               else                                                       8d7s23
                nrowj=noc(isblk(1,is))*noc(isblk(2,is))                   8d7s23
                ncolj=noc(isblk(3,is))*noc(isblk(4,is))                   8d7s23
                nrowjd=irefo(isblk(1,is))*irefo(isblk(2,is))            8d11s23
                isw=1                                                     8d7s23
               end if                                                     8d7s23
               ncoljd=nvirt(isblk(3,is))*nvirt(isblk(4,is))             8d11s23
               if(jsd.eq.isblk(4,is))then                                 8d7s23
                if(isblk(1,is).eq.jsa.and.isblk(2,is).eq.jsb)then         8d7s23
                 jj=j4o(is)+nrowj*ncolj*(ikind-1)                         8d7s23
                 do i4=0,noc(isblk(4,is))-1                               8d7s23
                  i3top=i4+isw*(noc(isblk(3,is))-1-i4)                    8d7s23
                  do i3=0,i3top                                           8d7s23
                   irec=i3+noc(isblk(3,is))*i4                            8d7s23
                   itri=((i4*(i4+1))/2)+i3                                8d7s23
                   icol=jj+nrowj*(itri+isw*(irec-itri))                   8d7s23
                   do i2=0,noc(isblk(2,is))-1                             8d7s23
                    i1top=i2+isw*(noc(isblk(1,is))-1-i2)                  8d7s23
                    do i1=0,i1top                                         8d7s23
                     irec=i1+noc(isblk(1,is))*i2                          8d7s23
                     itri=((i2*(i2+1))/2)+i1                              8d7s23
                     irow=icol+itri+isw*(irec-itri)                       8d7s23
                     iad=itmpt+i1+nbasdws(jsa)*(i2+nbasdws(jsb)*(i3       8d7s23
     $                    +nbasdws(jsc)*i4))                              8d7s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                            8d7s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                   ('goal_a'),bc(iad),iad
                    end do                                                8d7s23
                   end do                                                 8d7s23
                  end do                                                  8d7s23
                 end do                                                   8d7s23
                 jj=jmatd(is)+nrowj*ncoljd*(ikind-1)                    10d23s23
                 do i4=0,nvirt(isblk(4,is))-1                           8d11s23
                  i4p=i4+noc(isblk(4,is))                               8d11s23
                  do i3=0,nvirt(isblk(3,is))-1                          8d11s23
                   i3p=i3+noc(isblk(3,is))                              8d11s23
                   irec=i3+nvirt(isblk(3,is))*i4                        8d11s23
                   icol=jj+nrowj*irec                                   10d23s23
                   do i2=0,noc(isblk(2,is))-1                           10d23s23
                    i1top=i2+isw*(noc(isblk(1,is))-1-i2)                10d23s23
                    do i1=0,i1top                                         8d7s23
                     irec=i1+noc(isblk(1,is))*i2                        10d23s23
                     itri=((i2*(i2+1))/2)+i1                              8d7s23
                     irow=icol+itri+isw*(irec-itri)                       8d7s23
                     iad=itmpt+i1+nbasdws(jsa)*(i2+nbasdws(jsb)*(i3p    10d23s23
     $                    +nbasdws(jsc)*i4p))                           8d11s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                            8d7s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                   ('goal Ja'),orig,bc(iad)*ass,bc(igoal),iad
                    end do                                                8d7s23
                   end do                                                 8d7s23
                  end do                                                  8d7s23
                 end do                                                   8d7s23
                else if(isblk(2,is).eq.jsa.and.isblk(1,is).eq.jsb)then   8d8s23
                 jj=j4o(is)+nrowj*ncolj*(ikind-1)                         8d7s23
                 do i4=0,noc(isblk(4,is))-1                               8d7s23
                  do i3=0,noc(isblk(3,is))-1                              8d7s23
                   icol=jj+nrowj*(i3+noc(isblk(3,is))*i4)                 8d7s23
                   jtri=i3+noc(isblk(3,is))*i4
                   do i2=0,noc(isblk(2,is))-1                             8d7s23
                    do i1=0,noc(isblk(1,is))-1                            8d7s23
                     irow=icol+i1+noc(isblk(1,is))*i2                     8d7s23
                     ktri=i1+noc(isblk(1,is))*i2
                     iad=itmpt+i2+nbasdws(jsa)*(i1+nbasdws(jsb)*(i3       8d7s23
     $                    +nbasdws(jsc)*i4))                              8d7s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                            8d7s23
                     if(abs(orig-bc(igoal)).gt.1d-10)then
                      write(6,*)('4ogoal')
                     end if
                    end do                                                8d7s23
                   end do                                                 8d7s23
                  end do                                                  8d7s23
                 end do                                                   8d7s23
                 jj=jmatd(is)+nrowj*ncoljd*(ikind-1)                    10d23s23
                 do i4=0,nvirt(isblk(4,is))-1                           8d11s23
                  i4p=i4+noc(isblk(4,is))                               8d11s23
                  do i3=0,nvirt(isblk(3,is))-1                          8d11s23
                   i3p=i3+noc(isblk(3,is))                              8d11s23
                   icol=jj+nrowj*(i3+nvirt(isblk(3,is))*i4)             10d23s23
                   jtri=i3+nvirt(isblk(3,is))*i4
                   do i2=0,noc(isblk(2,is))-1                           10d23s23
                    do i1=0,noc(isblk(1,is))-1                          10d23s23
                     irow=icol+i1+noc(isblk(1,is))*i2                   10d23s23
                     iad=itmpt+i2+nbasdws(jsa)*(i1+nbasdws(jsb)*(i3p    10d23s23
     $                    +nbasdws(jsc)*i4p))                              8d7s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                            8d7s23
                     if(abs(orig-bc(igoal)).gt.1d-10)then
                      write(6,*)('jgoal'),orig,bc(iad)*ass,bc(igoal)
                     end if
                    end do                                                8d7s23
                   end do                                                 8d7s23
                  end do                                                  8d7s23
                 end do                                                   8d7s23
                end if                                                    8d7s23
               end if
               if(jsd.eq.isblk(3,is))then                                 8d7s23
                if(isblk(1,is).eq.jsa.and.isblk(2,is).eq.jsb)then         8d7s23
                 jj=j4o(is)+nrowj*ncolj*(ikind-1)                         8d7s23
                 do i4=0,noc(isblk(4,is))-1                               8d7s23
                  i3top=i4+isw*(noc(isblk(3,is))-1-i4)                    8d7s23
                  do i3=0,i3top                                           8d7s23
                   irec=i3+noc(isblk(3,is))*i4                            8d7s23
                   itri=((i4*(i4+1))/2)+i3                                8d7s23
                   jtri=itri+isw*(irec-itri)
                   icol=jj+nrowj*(itri+isw*(irec-itri))                   8d7s23
                   do i2=0,noc(isblk(2,is))-1                             8d7s23
                    i1top=i2+isw*(noc(isblk(1,is))-1-i2)                  8d7s23
                    do i1=0,i1top                                         8d7s23
                     irec=i1+noc(isblk(1,is))*i2                          8d7s23
                     itri=((i2*(i2+1))/2)+i1                              8d7s23
                     ktri=itri+isw*(irec-itri)
                     irow=icol+itri+isw*(irec-itri)                       8d7s23
                     iad=itmpt+i1+nbasdws(jsa)*(i2+nbasdws(jsb)*(i4       8d7s23
     $                    +nbasdws(jsc)*i3))                              8d7s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                            8d7s23
                     if(abs(orig-bc(igoal)).gt.1d-10)then
                      write(6,*)('4ogoalb'),orig,bc(iad)*ass,bc(irow)
                     end if
                    end do                                                8d7s23
                   end do                                                 8d7s23
                  end do                                                  8d7s23
                 end do                                                   8d7s23
                 jj=jmatd(is)+nrowj*ncoljd*(ikind-1)                    10d23s23
                 do i4=0,nvirt(isblk(4,is))-1                           8d11s23
                  i4p=i4+noc(isblk(4,is))                               8d11s23
                  do i3=0,nvirt(isblk(3,is))-1                          8d11s23
                   i3p=i3+noc(isblk(3,is))                              8d11s23
                   irec=i3+nvirt(isblk(3,is))*i4                        8d11s23
                   jtri=irec
                   icol=jj+nrowj*irec                                   10d23s23
                   do i2=0,noc(isblk(2,is))-1                           10d23s23
                    i1top=i2+isw*(noc(isblk(1,is))-1-i2)                10d23s23
                    do i1=0,i1top                                         8d7s23
                     irec=i1+noc(isblk(1,is))*i2                        10d23s23
                     itri=((i2*(i2+1))/2)+i1                              8d7s23
                     ktri=itri+isw*(irec-itri)
                     irow=icol+itri+isw*(irec-itri)                       8d7s23
                     iad=itmpt+i1+nbasdws(jsa)*(i2+nbasdws(jsb)*(i4p    10d23s23
     $                    +nbasdws(jsc)*i3p))                           8d11s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                            8d7s23
                     if(abs(orig-bc(igoal)).gt.1d-10)then
                      write(6,*)('jgoalb'),orig,bc(iad)*ass,bc(igoal),
     $                     bc(igoal)*ass
                     end if
                    end do                                                8d7s23
                   end do                                                 8d7s23
                  end do                                                  8d7s23
                 end do                                                   8d7s23
                else if(isblk(2,is).eq.jsa.and.isblk(1,is).eq.jsb)then    8d7s23
                 jj=j4o(is)+nrowj*ncolj*(ikind-1)                         8d7s23
                 do i4=0,noc(isblk(4,is))-1                               8d7s23
                  do i3=0,noc(isblk(3,is))-1                              8d7s23
                   icol=jj+nrowj*(i3+noc(isblk(3,is))*i4)                 8d7s23
                   do i2=0,noc(isblk(2,is))-1                             8d7s23
                    do i1=0,noc(isblk(1,is))-1                            8d7s23
                     irow=icol+i1+noc(isblk(1,is))*i2                     8d7s23
                     iad=itmpt+i2+nbasdws(jsa)*(i1+nbasdws(jsb)*(i4       8d7s23
     $                    +nbasdws(jsc)*i3))                              8d7s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                            8d7s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                    ('4ogoal d'),orig,bc(iad)*ass,bc(igoal),iad
                    end do                                                8d7s23
                   end do                                                 8d7s23
                  end do                                                  8d7s23
                 end do                                                   8d7s23
                 jj=jmatd(is)+nrowj*ncoljd*(ikind-1)                    10d23s23
                 do i4=0,nvirt(isblk(4,is))-1                           8d11s23
                  i4p=i4+noc(isblk(4,is))                               8d11s23
                  do i3=0,nvirt(isblk(3,is))-1                          8d11s23
                   i3p=i3+noc(isblk(3,is))                              8d11s23
                   icol=jj+nrowj*(i3+nvirt(isblk(3,is))*i4)             10d23s23
                   do i2=0,noc(isblk(2,is))-1                           10d23s23
                    do i1=0,noc(isblk(1,is))-1                          10d23s23
                     irow=icol+i1+noc(isblk(1,is))*i2                   10d23s23
                     iad=itmpt+i2+nbasdws(jsa)*(i1+nbasdws(jsb)*(i4p    10d23s23
     $                    +nbasdws(jsc)*i3p))                           8d11s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                            8d7s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                    ('jgoal d'),orig,bc(iad)*ass,bc(irow),iad,
     $                    jsa,isblk(2,is),jsb,isblk(1,is),jsc,
     $                    isblk(4,is),jsd,isblk(3,is),i2,i1,i4,i4p,i3,
     $                    i3p,i2+nbasdws(jsa)*i1,i4p+nbasdws(jsc)*i3p,
     $                    nbasdws(jsa),nbasdws(jsb),nbasdws(jsc),
     $                    nbasdws(jsd)
                    end do                                                8d7s23
                   end do                                                 8d7s23
                  end do                                                  8d7s23
                 end do                                                   8d7s23
                end if                                                    8d7s23
               end if
               if(isblk(2,is).eq.jsd)then                                 8d7s23
                if(isblk(3,is).eq.jsa.and.isblk(4,is).eq.jsb)then         8d7s23
                 jj=j4o(is)+nrowj*ncolj*(ikind-1)                         8d7s23
                 do i4=0,noc(isblk(4,is))-1                               8d7s23
                  i3top=i4+isw*(noc(isblk(3,is))-1-i4)                    8d7s23
                  do i3=0,i3top                                           8d7s23
                   irec=i3+noc(isblk(3,is))*i4                            8d7s23
                   itri=((i4*(i4+1))/2)+i3                                8d7s23
                   icol=jj+nrowj*(itri+isw*(irec-itri))                   8d7s23
                   do i2=0,noc(isblk(2,is))-1                             8d7s23
                    i1top=i2+isw*(noc(isblk(1,is))-1-i2)                  8d7s23
                    do i1=0,i1top                                         8d7s23
                     irec=i1+noc(isblk(1,is))*i2                          8d7s23
                     itri=((i2*(i2+1))/2)+i1                              8d7s23
                     irow=icol+itri+isw*(irec-itri)                       8d7s23
                     iad=itmpt+i3+nbasdws(jsa)*(i4+nbasdws(jsb)*(i1       8d7s23
     $                    +nbasdws(jsc)*i2))                              8d7s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                            8d7s23
                     if(abs(orig-bc(igoal)).gt.1d-10)then
                      write(6,*)('4ogoal e'),orig,bc(iad)*ass,bc(irow),
     $                     iad
                     end if
                    end do                                                8d7s23
                   end do                                                 8d7s23
                  end do                                                  8d7s23
                 end do                                                   8d7s23
                 jj=jmatd(is)+nrowj*ncoljd*(ikind-1)                    10d23s23
                 do i4=0,nvirt(isblk(4,is))-1                           8d11s23
                  i4p=i4+noc(isblk(4,is))                               8d11s23
                  do i3=0,nvirt(isblk(3,is))-1                          8d11s23
                   i3p=i3+noc(isblk(3,is))                              8d11s23
                   irec=i3+nvirt(isblk(3,is))*i4                        8d11s23
                   icol=jj+nrowj*irec                                   10d23s23
                   do i2=0,noc(isblk(2,is))-1                           10d23s23
                    i1top=i2+isw*(noc(isblk(1,is))-1-i2)                10d23s23
                    do i1=0,i1top                                         8d7s23
                     irec=i1+noc(isblk(1,is))*i2                        10d23s23
                     itri=((i2*(i2+1))/2)+i1                              8d7s23
                     irow=icol+itri+isw*(irec-itri)                       8d7s23
                     iad=itmpt+i3p+nbasdws(jsa)*(i4p+nbasdws(jsb)*(i1   10d23s23
     $                    +nbasdws(jsc)*i2))                            10d23s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                            8d7s23
                     if(abs(orig-bc(igoal)).gt.1d-10)then
                      write(6,*)('jgoal e'),orig,bc(iad)*ass,bc(igoal)
                     end if
                    end do                                                8d7s23
                   end do                                                 8d7s23
                  end do                                                  8d7s23
                 end do                                                   8d7s23
                else if(isblk(4,is).eq.jsa.and.isblk(3,is).eq.jsb)then   8d8s23
                 jj=j4o(is)+nrowj*ncolj*(ikind-1)                         8d7s23
                 do i4=0,noc(isblk(4,is))-1                               8d7s23
                  do i3=0,noc(isblk(3,is))-1                              8d7s23
                   icol=jj+nrowj*(i3+noc(isblk(3,is))*i4)                 8d7s23
                   do i2=0,noc(isblk(2,is))-1                             8d7s23
                    do i1=0,noc(isblk(1,is))-1                            8d7s23
                     irow=icol+i1+noc(isblk(1,is))*i2                     8d7s23
                     iad=itmpt+i4+nbasdws(jsa)*(i3+nbasdws(jsb)*(i1       8d7s23
     $                    +nbasdws(jsc)*i2))                              8d7s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                            8d7s23
                     if(abs(orig-bc(igoal)).gt.1d-10)then
                      write(6,*)('goal F'),orig,bc(iad),bc(igoal)
                     end if
                    end do                                                8d7s23
                   end do                                                 8d7s23
                  end do                                                  8d7s23
                 end do                                                   8d7s23
                 jj=jmatd(is)+nrowj*ncoljd*(ikind-1)                    10d23s23
                 do i4=0,nvirt(isblk(4,is))-1                           8d11s23
                  i4p=i4+noc(isblk(4,is))                               8d11s23
                  do i3=0,nvirt(isblk(3,is))-1                          8d11s23
                   i3p=i3+noc(isblk(3,is))                              8d11s23
                   icol=jj+nrowj*(i3+nvirt(isblk(3,is))*i4)             10d23s23
                   do i2=0,noc(isblk(2,is))-1                           10d23s23
                    do i1=0,noc(isblk(1,is))-1                          10d23s23
                     irow=icol+i1+noc(isblk(1,is))*i2                   10d23s23
                     iad=itmpt+i4p+nbasdws(jsa)*(i3p+nbasdws(jsb)*(i1   10d23s23
     $                    +nbasdws(jsc)*i2))                            10d23s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                            8d7s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                    ('goalf2'),orig,bc(iad)*ass,bc(irow)
                    end do                                                8d7s23
                   end do                                                 8d7s23
                  end do                                                  8d7s23
                 end do                                                   8d7s23
                end if                                                    8d7s23
               end if                                                     8d7s23
               if(isblk(1,is).eq.jsd)then                                 8d7s23
                if(isblk(3,is).eq.jsa.and.isblk(4,is).eq.jsb)then         8d7s23
                 jj=j4o(is)+nrowj*ncolj*(ikind-1)                         8d7s23
                 do i4=0,noc(isblk(4,is))-1                               8d7s23
                  i3top=i4+isw*(noc(isblk(3,is))-1-i4)                    8d7s23
                  do i3=0,i3top                                           8d7s23
                   irec=i3+noc(isblk(3,is))*i4                            8d7s23
                   itri=((i4*(i4+1))/2)+i3                                8d7s23
                   jtri=itri+isw*(irec-itri)
                   icol=jj+nrowj*(itri+isw*(irec-itri))                   8d7s23
                   do i2=0,noc(isblk(2,is))-1                             8d7s23
                    i1top=i2+isw*(noc(isblk(1,is))-1-i2)                  8d7s23
                    do i1=0,i1top                                         8d7s23
                     irec=i1+noc(isblk(1,is))*i2                          8d7s23
                     itri=((i2*(i2+1))/2)+i1                              8d7s23
                     ktri=itri+isw*(irec-itri)
                     irow=icol+itri+isw*(irec-itri)                       8d7s23
                     iad=itmpt+i3+nbasdws(jsa)*(i4+nbasdws(jsb)*(i2       8d7s23
     $                    +nbasdws(jsc)*i1))                              8d7s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                            8d7s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                    ('goalg2'),orig,bc(iad)*ass,bc(igoal)
                    end do                                                8d7s23
                   end do                                                 8d7s23
                  end do                                                  8d7s23
                 end do                                                   8d7s23
                 jj=jmatd(is)+nrowj*ncoljd*(ikind-1)                    10d23s23
                 do i4=0,nvirt(isblk(4,is))-1                           8d11s23
                  i4p=i4+noc(isblk(4,is))                               8d11s23
                  do i3=0,nvirt(isblk(3,is))-1                          8d11s23
                   i3p=i3+noc(isblk(3,is))                              8d11s23
                   irec=i3+nvirt(isblk(3,is))*i4                        8d11s23
                   jtri=irec
                   icol=jj+nrowj*irec                                   10d23s23
                   do i2=0,noc(isblk(2,is))-1                           10d23s23
                    i1top=i2+isw*(noc(isblk(1,is))-1-i2)                10d23s23
                    do i1=0,i1top                                         8d7s23
                     irec=i1+noc(isblk(1,is))*i2                        10d23s23
                     itri=((i2*(i2+1))/2)+i1                              8d7s23
                     ktri=itri+isw*(irec-itri)
                     irow=icol+itri+isw*(irec-itri)                       8d7s23
                     iad=itmpt+i3p+nbasdws(jsa)*(i4p+nbasdws(jsb)*(i2   10d23s23
     $                    +nbasdws(jsc)*i1))                            10d23s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                            8d7s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                    ('goalg2'),orig,bc(iad)*ass,bc(igoal)
                    end do                                                8d7s23
                   end do                                                 8d7s23
                  end do                                                  8d7s23
                 end do                                                   8d7s23
                else if(isblk(4,is).eq.jsa.and.isblk(3,is).eq.jsb)then    8d7s23
                 jj=j4o(is)+nrowj*ncolj*(ikind-1)                         8d7s23
                 do i4=0,noc(isblk(4,is))-1                               8d7s23
                  do i3=0,noc(isblk(3,is))-1                              8d7s23
                   icol=jj+nrowj*(i3+noc(isblk(3,is))*i4)                 8d7s23
                   do i2=0,noc(isblk(2,is))-1                             8d7s23
                    do i1=0,noc(isblk(1,is))-1                            8d7s23
                     irow=icol+i1+noc(isblk(1,is))*i2                     8d7s23
                     iad=itmpt+i4+nbasdws(jsa)*(i3+nbasdws(jsb)*(i2       8d7s23
     $                    +nbasdws(jsc)*i1))                              8d7s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                            8d7s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                    ('goal i'),bc(iad),iad
                    end do                                                8d7s23
                   end do                                                 8d7s23
                  end do                                                  8d7s23
                 end do                                                   8d7s23
                 jj=jmatd(is)+nrowj*ncoljd*(ikind-1)                    10d23s23
                 do i4=0,nvirt(isblk(4,is))-1                           8d11s23
                  i4p=i4+noc(isblk(4,is))                               8d11s23
                  do i3=0,nvirt(isblk(3,is))-1                          8d11s23
                   i3p=i3+noc(isblk(3,is))                              8d11s23
                   icol=jj+nrowj*(i3+nvirt(isblk(3,is))*i4)             10d23s23
                   do i2=0,noc(isblk(2,is))-1                           10d23s23
                    do i1=0,noc(isblk(1,is))-1                          10d23s23
                     irow=icol+i1+noc(isblk(1,is))*i2                   10d23s23
                     iad=itmpt+i4p+nbasdws(jsa)*(i3p+nbasdws(jsb)*(i2   10d23s23
     $                    +nbasdws(jsc)*i1))                            10d23s23
                     orig=bc(igoal)
                     bc(irow)=bc(irow)+bc(iad)                            8d7s23
                     if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                    ('goal i'),orig,bc(iad)*ass,bc(igoal),iad
                    end do                                                8d7s23
                   end do                                                 8d7s23
                  end do                                                  8d7s23
                 end do                                                   8d7s23
                end if                                                    8d7s23
               end if                                                     8d7s23
              end do                                                      8d7s23
             end if                                                     8d8s23
            end do                                                      8d8s23
           end do                                                       8d7s23
          end if                                                        8d7s23
         end do
        end do
       end do
      end do
      write(6,*)('ecore2, ecore2k '),ecore2,ecore2k
      write(6,*)('sum(shift) '),potdws+ecore1+ecore2+ecore2k
      write(6,*)('eoa: '),eoa,('eaa1 '),eaa1,eoa+eaa1
      write(6,*)('eoad: '),(eoad(isb),isb=1,nsymb)
      write(6,*)('e0x'),e0x,(iuseage(j,1),j=1,8)
      write(6,*)('e1x'),e1x,(iuseage(j,2),j=1,8)
      write(6,*)('ejj'),ejj,(iuseage(j,3),j=1,8)
      write(6,*)('ekk'),ekk,(iuseage(j,4),j=1,8)
      write(6,*)('e3x'),e3x,(iuseage(j,5),j=1,8)
      write(6,*)('e4x'),e4x,(iuseage(j,6),j=1,8)
      sumt=0d0                                                          5d10s24
      do kk=1,npart
       sum=0d0
       if(ichoice(1).ne.0)then
        write(6,*)('dcore1 for part '),kk,dcore1(kk)                     5d10s24
        sum=sum+dcore1(kk)
       end if
       if(ichoice(2).ne.0)then
        write(6,*)('daa1 for part '),kk,daa1(kk)
        sum=sum+daa1(kk)
       end if
       if(ichoice(3).ne.0)then
        write(6,*)('dcore2 for part '),kk,dcore2(kk)+dcore2k(kk),
     $      dcore2(kk),dcore2k(kk)
        sum=sum+dcore2(kk)+dcore2k(kk)
       end if
       write(6,*)('dshift for part '),kk,dcore1(kk)+dcore2(kk)
     $      +dcore2k(kk)
       if(ichoice(4).ne.0)then
        write(6,*)('deoa for part '),kk,deoa(kk)
        sum=sum+deoa(kk)
       end if
       if(ichoice(5).ne.0)then
        write(6,*)('de0x for part '),kk,de0x(kk)
        sum=sum+de0x(kk)
       end if
       if(ichoice(6).ne.0)then
        write(6,*)('de1x for part '),kk,de1x(kk)
        sum=sum+de1x(kk)
       end if
       if(ichoice(7).ne.0)then
        write(6,*)('dejj for part '),kk,dejj(kk)
        sum=sum+dejj(kk)
       end if
       if(ichoice(8).ne.0)then
        write(6,*)('dekk for part '),kk,dekk(kk)
        sum=sum+dekk(kk)
       end if
       if(ichoice(9).ne.0)then                                          5d10s24
        write(6,*)('de3x for part '),kk,de3x(kk)                        5d10s24
        sum=sum+de3x(kk)                                                5d10s24
       end if
       if(ichoice(10).ne.0)then                                         5d13s24
        write(6,*)('de4x for part '),kk,de4x(kk)                        5d10s24
        sum=sum+de4x(kk)                                                5d10s24
       end if
       write(6,*)('sum '),sum
       sumt=sumt+sum                                                    5d10s24
      end do
      write(6,*)('grand total: '),sumt                                  5d10s24
      do ipart=1,npart                                                  5d15s24
       sum4vt(ipart)=0d0                                                       5d15s24
       sum4tv(ipart)=0d0                                                       5d15s24
      end do                                                            5d15s24
      sumxxxx=0d0
       do isb=1,nsymb                                                   5d15s24
        if(nvirt(isb).gt.0)then                                         5d15s24
         write(6,*)('for symmetry block '),isb
         do l=1,4
          iad=ie4v(isb,l)+nvirt(isb)*noc(isb)
          write(6,*)('e4x for l = '),l,ie4v(isb,l),iad
          call prntm2(bc(iad),nvirt(isb),nvirt(isb),
     $         nvirt(isb))
         end do
         write(6,*)('2+4')
         do i=0,nvirt(isb)*nbasdws(isb)-1
          bc(ie4v(isb,2)+i)=bc(ie4v(isb,2)+i)+bc(ie4v(isb,4)+i)
          bc(ie4v(isb,4)+i)=0d0
         end do
         iad=ie4v(isb,2)+nvirt(isb)*noc(isb)
         call prntm2(bc(iad),nvirt(isb),nvirt(isb),
     $         nvirt(isb))
         do i=0,nvirt(isb)*nbasdws(isb)-1
          bc(ie4v(isb,1)+i)=bc(ie4v(isb,1)+i)+bc(ie4v(isb,2)+i)
     $         +bc(ie4v(isb,3)+i)
          bc(ie4v(isb,2)+i)=0d0
          bc(ie4v(isb,3)+i)=0d0
         end do
         write(6,*)('1+2+3+4')
         iad=ie4v(isb,1)+nvirt(isb)*noc(isb)
         call prntm2(bc(iad),nvirt(isb),nvirt(isb),
     $         nvirt(isb))
         do ipart=1,npart
          do i=0,nbasdws(isb)-1                                          5d15s24
           do j=0,nvirt(isb)-1                                           5d15s24
            jp=j+noc(isb)                                                5d15s24
            iad1=ie4v(isb,1)+j+nvirt(isb)*i
            iad2=ie4v(isb,2)+j+nvirt(isb)*i
            iad3=ie4v(isb,3)+j+nvirt(isb)*i
            iad4=ie4v(isb,4)+j+nvirt(isb)*i
            add=bc(iad1)+bc(iad2)+bc(iad3)+bc(iad4)
            ij=idaprt(ipart,isb)+i+nbasdws(isb)*jp                         5d15s24
            ji=idaprt(ipart,isb)+jp+nbasdws(isb)*i                      5d15s24
            sum4tv(ipart)=sum4tv(ipart)+add*bc(ij)
           end do
          end do
          do i=0,nvirt(isb)-1
           ip=i+noc(isb)
           do j=0,nvirt(isb)-1
            jp=j+noc(isb)
            iad1=ie4v(isb,1)+j+nvirt(isb)*ip
            ij=idaprt(ipart,isb)+ip+nbasdws(isb)*jp
            sum4vt(ipart)=sum4vt(ipart)+bc(iad1)*bc(ij)
           end do
          end do
         end do
        end if                                                          5d15s24
       end do                                                           5d15s24
       write(6,*)('sum4tv*1: '),(1d0*sum4tv(ipart),ipart=1,npart)
       write(6,*)('sum4vt: '),(sum4vt(ipart),ipart=1,npart)
       write(6,*)('sum4vt-vt: '),(sum4tv(ipart)-sum4vt(ipart),
     $      ipart=1,npart)
      if(ichoice(3).ne.0)then                                           10d18s23
      write(6,*)('2e part of shift: ')
      do is=1,nsdlk                                                     8d7s23
       if(isblk(1,is).eq.isblk(2,is))then                               8d7s23
        if(min(idoubo(isblk(1,is)),idoubo(isblk(3,is))).gt.0)then       8d8s23
         nrow=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                  8d7s23
         ncol=(noc(isblk(3,is))*(noc(isblk(3,is))+1))/2                  8d7s23
         jj=j4o(is)                                                     8d8s23
         do ikind=1,npartp                                              6d26s24
          sum=0d0                                                       6d26s24
          write(6,*)('sum for J '),(isblk(i14,is),i14=1,4)
          do i34=1,idoubo(isblk(3,is))                                  8d8s23
           icol=((i34*(i34+1))/2)-1                                     8d8s23
           do i12=1,idoubo(isblk(1,is))                                 8d8s23
            irow=((i12*(i12+1))/2)-1                                    8d8s23
            iad=jj+irow+nrow*icol                                       8d8s23
            sum=sum+2d0*bc(iad)                                         6d26s24
            if(abs(bc(iad)).gt.1d-10)write(6,*)('J '),bc(iad)*ass,ikind,
     $           (isblk(j,is),j=1,4),iad,i12,i34,trace(ikind,1),bc(iad)
           end do                                                       8d8s23
          end do                                                        8d8s23
          if(ikind.eq.npartp)then                                       6d26s24
         orig=tcannon
           tcannon=tcannon+sum                                          6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannonc '),
     $          orig,sum,tcannon
          else                                                          6d26s24
           trace(ikind,1)=trace(ikind,1)+sum                            6d26s24
          end if                                                        6d26s24
          jj=jj+nrow*ncol                                               8d8s23
         end do                                                         8d8s23
        end if                                                          8d8s23
       end if                                                           8d8s23
       if(isblk(1,is).eq.isblk(4,is))then                               8d8s23
        if(min(idoubo(isblk(1,is)),idoubo(isblk(2,is))).gt.0)then       8d8s23
         if(isblk(1,is).eq.isblk(2,is))then                             8d8s23
          nrow=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                  8d7s23
          ncol=(noc(isblk(3,is))*(noc(isblk(3,is))+1))/2                  8d7s23
          isw=0                                                         8d8s23
         else                                                           8d8s23
          nrow=noc(isblk(1,is))*noc(isblk(2,is))                        8d8s23
          ncol=noc(isblk(3,is))*noc(isblk(4,is))                        8d8s23
          isw=1                                                         8d8s23
         end if                                                         8d8s23
         jj=j4o(is)                                                     8d8s23
         ff=-1d0                                                        8d8s23
         if(isblk(1,is).ne.isblk(2,is))ff=-2d0                          8d8s23
         do ikind=1,npartp                                              6d26s24
          sum=0d0                                                       6d26s24
          write(6,*)('sum for K '),(isblk(i14,is),i14=1,4)
          do i14=0,idoubo(isblk(1,is))-1                                8d8s23
           do i23=0,idoubo(isblk(2,is))-1                               8d8s23
            icol=i23+noc(isblk(3,is))*i14                               8d8s23
            irow=i14+noc(isblk(1,is))*i23                               8d8s23
            ix=max(i14,i23)                                             8d8s23
            in=min(i14,i23)                                             8d8s23
            itri=((ix*(ix+1))/2)+in                                     8d8s23
            icol=itri+isw*(icol-itri)                                   8d8s23
            irow=itri+isw*(irow-itri)                                   8d8s23
            iad=jj+irow+nrow*icol                                       8d8s23
            sum=sum+ff*bc(iad)                                          6d26s24
            if(abs(bc(iad)).gt.1d-10)write(6,*)('K '),ff,bc(iad)*ass,
     $           ikind,
     $           (isblk(j,is),j=1,4),iad,i14,i23,bc(igoal),iad-igoal,
     $           trace(ikind,1),bc(iad)
           end do                                                       8d8s23
          end do
          if(ikind.eq.npartp)then                                       6d26s24
         orig=tcannon
           tcannon=tcannon+sum                                          6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannond '),
     $          orig,sum,tcannon
          else                                                          6d26s24
           trace(ikind,1)=trace(ikind,1)+sum                            6d26s24
          end if                                                        6d26s24
          jj=jj+nrow*ncol                                               8d8s23
         end do                                                         8d8s23
        end if                                                          8d8s23
       end if                                                           8d7s23
      end do                                                            8d7s23
      end if                                                            10d18s23
      if(ichoice(4).ne.0)then                                           10d18s23
      write(6,*)('mixed da part ')
      do isb=1,nsymb                                                    8d8s23
       if(irefo(isb).ge.0)then                                          8d8s23
        write(6,*)('for symmetry block '),isb
        do is=1,nsdlk1                                                  10d22s23
         if(isblk1(3,is).eq.isblk1(4,is).and.isblk1(3,is).eq.isb)then   7d20s23
          jsb=isblk1(1,is)                                              10d22s23
          if(idoubo(jsb).gt.0)then                                      10d22s23
           nr=(noc(jsb)*(noc(jsb)+1))/2                                 10d22s23
           ii=jonex(is)                                                 10d22s23
           write(6,*)('block 0 '),(isblk1(j,is),j=1,4),trace(2,1)
           do ikind=1,npartp                                            6d26s24
            sum=0d0                                                     6d26s24
            do i2=0,nvirt(isb)-1                                        10d22s23
             i2p=i2+irefo(isb)                                          10d22s23
             iadd=iden(isb)+nh0av(isb)*i2p                              10d22s23
             do i1=0,irefo(isb)-1                                       10d23s23
              i1p=i1+idoubo(isb)                                        10d23s23
              iii=ii+nr*(i1p+noc(isb)*i2)                               10d23s23
              ff=4d0*bc(iadd+i1)                                         10d22s23
              do id=1,idoubo(jsb)                                        10d22s23
               icol=iii+((id*(id+1))/2)-1                               10d23s23
               sum=sum+ff*bc(icol)                                      6d26s24
              end do                                                     10d22s23
             end do                                                      10d22s23
            end do                                                      10d22s23
            if(ikind.eq.npartp)then                                     6d26s24
         orig=tcannon
             tcannon=tcannon+sum                                        6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannone '),
     $          orig,sum,tcannon
            else                                                        6d26s24
             trace(ikind,1)=trace(ikind,1)+sum                          6d26s24
            end if                                                      6d26s24
            ii=ii+nr*noc(isb)*nvirt(isb)                                10d23s23
           end do                                                       10d22s23
          end if                                                        10d22s23
         end if                                                         10d22s23
         if(isblk1(4,is).eq.isb)then                                    10d22s23
          if(isblk1(4,is).eq.isblk1(1,is))then                          10d22s23
           jsb=isblk1(2,is)                                             10d22s23
           if(idoubo(jsb).gt.0)then                                     10d22s23
            write(6,*)('block 1: '),(isblk1(j,is),j=1,4)
            if(isblk1(1,is).eq.isblk1(2,is))then                         10d22s23
             nr=(noc(isb)*(noc(isb)+1))/2                                10d22s23
             isw=0                                                       10d22s23
            else                                                         10d22s23
             nr=noc(isb)*noc(jsb)                                        10d22s23
             isw=1                                                       10d22s23
            end if                                                       10d22s23
            ii=jonex(is)                                                 10d22s23
            do ikind=1,npartp                                           6d26s24
             sum=0d0                                                    6d26s24
             do iv=0,nvirt(isb)-1                                       10d22s23
              ivp=iv+irefo(isb)                                         10d22s23
              iadd=iden(isb)+nh0av(isb)*ivp                             10d22s23
              do id=0,idoubo(jsb)-1                                     10d22s23
               iii=ii+nr*(id+noc(jsb)*iv)                               10d22s23
               do in=0,irefo(isb)-1                                     10d22s23
                inp=in+idoubo(isb)                                      10d22s23
                irec=inp+noc(isb)*id                                    10d22s23
                itri=((inp*(inp+1))/2)+id                               10d22s23
                itri=itri+isw*(irec-itri)
                sum=sum-2d0*bc(iadd+in)*bc(itri+iii)                    6d26s24
               end do                                                   10d22s23
              end do                                                    10d22s23
             end do                                                     10d22s23
             ii=ii+nvirt(isb)*noc(jsb)*nr                               10d22s23
             if(ikind.eq.npartp)then                                    6d26s24
         orig=tcannon
              tcannon=tcannon+sum                                       6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannonf '),
     $          orig,sum,tcannon
             else                                                       6d26s24
              trace(ikind,1)=trace(ikind,1)+sum                         6d26s24
             end if                                                     6d26s24
            end do                                                      10d22s23
           end if                                                       10d22s23
          else if(isblk1(4,is).eq.isblk1(2,is))then                     10d22s23
           jsb=isblk1(1,is)                                             10d22s23
           if(idoubo(jsb).gt.0)then                                     10d23s23
            write(6,*)('block 2: '),(isblk1(j,is),j=1,4)
            nr=noc(isb)*noc(jsb)                                         10d22s23
            ii=jonex(is)                                                 10d22s23
            do ikind=1,npartp                                           6d26s24
             sum=0d0                                                    6d26s24
             do iv=0,nvirt(isb)-1                                        10d22s23
              ivp=iv+irefo(isb)                                          10d22s23
              iadd=iden(isb)+nh0av(isb)*ivp                              10d22s23
              do id=0,idoubo(jsb)-1                                      10d22s23
               iii=ii+nr*(id+noc(jsb)*iv)                                10d22s23
               do in=0,irefo(isb)-1                                      10d22s23
                inp=in+idoubo(isb)                                       10d22s23
                irec=id+noc(jsb)*inp
                sum=sum-2d0*bc(iadd+in)*bc(irec+iii)                    6d26s24
               end do                                                    10d22s23
              end do                                                     10d22s23
             end do                                                      10d22s23
             ii=ii+nr*noc(jsb)*nvirt(isb)                                10d22s23
             if(ikind.eq.npartp)then                                    6d26s24
         orig=tcannon
              tcannon=tcannon+sum                                       6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannong '),
     $          orig,sum,tcannon
             else                                                       6d26s24
              trace(ikind,1)=trace(ikind,1)+sum                         6d26s24
             end if                                                     6d26s24
            end do                                                       10d22s23
           end if
          end if                                                        10d23s23
         end if                                                         10d22s23
        end do                                                          7d20s23
        do is=1,nsdlk
         if(isblk(1,is).eq.isblk(2,is).and.isblk(3,is).eq.isb)then       8d8s23
          jsb=isblk(1,is)                                                8d8s23
          write(6,*)('consider J for isb,jsb = '),isb,jsb
          if(idoubo(jsb).gt.0)then                                      8d8s23
           write(6,*)('J from '),(isblk(j,is),j=1,4),trace(2,1)
           nr=(noc(jsb)*(noc(jsb)+1))/2                                 8d8s23
           nc=(noc(isb)*(noc(isb)+1))/2
           jj=j4o(is)                                                   8d8s23
           do ikind=1,npartp                                            6d26s24
            sum=0d0                                                     6d26s24
            write(6,*)('for kind = '),ikind
            do i4=0,irefo(isb)-1                                         8d8s23
             i4p=i4+idoubo(isb)                                          8d8s23
             do i3=0,irefo(isb)-1                                        8d8s23
              i3p=i3+idoubo(isb)                                         8d8s23
              iadd=iden(isb)+i3+nh0av(isb)*i4                           8d8s23
              ix=max(i3p,i4p)                                           8d8s23
              in=min(i3p,i4p)                                           8d8s23
              icol=jj+nr*(((ix*(ix+1))/2)+in)-1
              ff=2d0*bc(iadd)
              do i12=1,idoubo(jsb)                                      8d8s23
               iad=icol+((i12*(i12+1))/2)
               sum=sum+ff*bc(iad)                                       6d26s24
              end do
             end do
            end do
            if(ikind.eq.npartp)then                                     6d26s24
         orig=tcannon
             tcannon=tcannon+sum                                        6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannonh '),
     $          orig,sum,tcannon
            else                                                        6d26s24
             trace(ikind,1)=trace(ikind,1)+sum                          6d26s24
            end if                                                      6d26s24
            jj=jj+nr*nc
           end do
           jjj=jmatd(is)                                                10d22s23
           do ikind=1,npartp                                            6d26s24
            sum=0d0                                                     6d26s24
            sz=0d0
            do i=0,nr*nvirt(isb)*nvirt(isb)-1
             sz=sz+bc(jjj+i)**2
            end do
            sz=sqrt(sz/dfloat(nr*nvirt(isb)*nvirt(isb)))
            write(6,*)('for kind = '),ikind,sz
            do i2=0,nvirt(isb)-1                                        10d22s23
             i2p=i2+irefo(isb)                                          10d22s23
             iadd=iden(isb)+irefo(isb)+nh0av(isb)*i2p                   10d22s23
             do i1=0,nvirt(isb)-1                                       10d22s23
              ff=2d0*bc(iadd+i1)                                        7d20s23
              do id=1,idoubo(jsb)                                       7d20s23
               icol=jjj+((id*(id+1))/2)-1                               10d22s23
               sum=sum+ff*bc(icol)                                      6d26s24
              end do                                                    7d20s23
              jjj=jjj+nr                                                10d22s23
             end do                                                     7d20s23
            end do                                                      7d20s23
            if(ikind.eq.npartp)then                                     6d26s24
         orig=tcannon
             tcannon=tcannon+sum                                        6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannoni '),
     $          orig,sum,tcannon
            else                                                        6d26s24
             trace(ikind,1)=trace(ikind,1)+sum                          6d26s24
            end if                                                      6d26s24
           end do
          end if                                                        8d8s23
         end if                                                         8d8s23
         if(isblk(2,is).eq.isb.and.isblk(3,is).eq.isb)then              8d8s23
          jsb=isblk(1,is)                                               8d8s23
          if(idoubo(jsb).gt.0)then                                      8d8s23
           write(6,*)('K from '),(isblk(j,is),j=1,4)
           if(isb.eq.jsb)then                                           8d8s23
            nr=(noc(isb)*(noc(isb)+1))/2                                8d8s23
            isw=0                                                       8d8s23
           else                                                         8d8s23
            nr=noc(isb)*noc(jsb)                                        8d8s23
            isw=1                                                       8d8s23
           end if                                                       8d8s23
           jj=j4o(is)                                                   8d8s23
           do ikind=1,npartp                                            6d26s24
            write(6,*)('for ikind = '),ikind
            call prntm2(bc(jj),nr,nr,nr)
            sum=0d0                                                     6d26s24
            do i14=0,idoubo(jsb)-1                                       8d8s23
             do i3=0,irefo(isb)-1                                        8d8s23
              i3p=i3+idoubo(isb)                                         8d8s23
              irecc=i3p+noc(isb)*i14                                      8d8s23
              itri=((i3p*(i3p+1))/2)+i14                                 8d8s23
              icol=jj+nr*(itri+isw*(irecc-itri))
              do i2=0,irefo(isb)-1                                       8d8s23
               i2p=i2+idoubo(isb)                                        8d8s23
               iadd=iden(isb)+i2+nh0av(isb)*i3                           8d8s23
               irec=i14+noc(jsb)*i2p                                    8d8s23
               itri=((i2p*(i2p+1))/2)+i14                               8d8s23
               iad=icol+itri+isw*(irec-itri)                            8d8s23
               sum=sum-bc(iad)*bc(iadd)                                 6d26s24
              end do                                                    8d8s23
             end do                                                     8d8s23
            end do                                                      8d8s23
            if(ikind.eq.npartp)then                                     6d26s24
         orig=tcannon
             tcannon=tcannon+sum                                        6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannonj '),
     $          orig,sum,tcannon
            else                                                        6d26s24
             trace(ikind,1)=trace(ikind,1)+sum                          6d26s24
            end if                                                      6d26s24
            jj=jj+nr*nr                                                 8d8s23
           end do                                                       8d8s23
          end if                                                        8d8s23
         else if(isblk(1,is).ne.isblk(2,is).and.isblk(1,is).eq.isb.and. 8d8s23
     $         isblk(4,is).eq.isb)then                                  8d8s23
          jsb=isblk(2,is)                                               8d8s23
          if(idoubo(jsb).gt.0)then                                      8d8s23
           write(6,*)('K from '),(isblk(j,is),j=1,4)
           if(isb.eq.jsb)then                                           8d8s23
            nr=(noc(isb)*(noc(isb)+1))/2                                8d8s23
            isw=0                                                       8d8s23
           else                                                         8d8s23
            nr=noc(isb)*noc(jsb)                                        8d8s23
            isw=1                                                       8d8s23
           end if                                                       8d8s23
           jj=j4o(is)                                                   8d8s23
           do ikind=1,npartp                                            6d26s24
            sum=0d0                                                     6d26s24
            write(6,*)('for ikind = '),ikind
            call prntm2(bc(jj),nr,nr,nr)
            do i4=0,irefo(isb)-1                                        8d8s23
             i4p=i4+idoubo(isb)                                         8d8s23
             do i23=0,idoubo(jsb)-1                                     8d8s23
              irecc=i23+noc(jsb)*i4p                                    8d8s23
              itri=((i4p*(i4p+1))/2)+i23                                8d8s23
              icol=jj+nr*(itri+isw*(irecc-itri))                        8d8s23
              do i1=0,irefo(isb)-1                                      8d8s23
               i1p=i1+idoubo(isb)                                       8d8s23
               iadd=iden(isb)+i1+nh0av(isb)*i4                          8d8s23
               irec=i1p+noc(isb)*i23                                    8d8s23
               itri=((i1p*(i1p+1))/2)+i23                               8d8s23
               iad=icol+itri+isw*(irec-itri)                            8d8s23
               sum=sum-bc(iadd)*bc(iad)                                 6d26s24
              end do                                                    8d8s23
             end do                                                     8d8s23
            end do                                                      8d8s23
            if(ikind.eq.npartp)then                                     6d26s24
         orig=tcannon
             tcannon=tcannon+sum                                        6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannonk '),
     $          orig,sum,tcannon
            else                                                        6d26s24
             trace(ikind,1)=trace(ikind,1)+sum                          6d26s24
            end if                                                      6d26s24
            jj=jj+nr*nr                                                 8d8s23
           end do                                                       8d8s23
          end if                                                        8d8s23
         end if                                                          8d8s23
        end do
       end if
       do is=1,nsdlkk                                                    10d22s23
        if(isblkk(3,is).eq.isblkk(4,is).and.isblkk(3,is).eq.isb)then    10d22s23
         jsb=isblkk(1,is)                                               10d22s23
         if(idoubo(jsb).gt.0)then                                       10d22s23
          write(6,*)('looking at K of '),(isblkk(j,is),j=1,4),npart
          nr=noc(jsb)*noc(jsb)                                          10d23s23
          kk=kmatd(is)                                                  10d22s23
          do ikind=1,npartp                                             6d26s24
           sum=0d0                                                      6d26s24
           sz=0d0
           do i=0,nr*nvirt(isb)*nvirt(isb)-1
            sz=sz+bc(kk+i)**2
           end do
           sz=sqrt(sz/dfloat(nr*nvirt(isb)*nvirt(isb)))
           write(6,*)('for kind '),ikind,('size is '),sz
           do iv=0,nvirt(isb)-1                                         10d22s23
            ivp=iv+irefo(isb)                                           10d22s23
            iadd=iden(isb)+irefo(isb)+nh0av(isb)*ivp                    10d22s23
            do jv=0,nvirt(isb)-1                                        10d22s23
             do id=0,idoubo(jsb)-1                                      10d22s23
              irec=id+noc(jsb)*id+kk                                    10d22s23
              sum=sum-bc(iadd+jv)*bc(irec)                              6d26s24
             end do                                                     10d22s23
             kk=kk+nr                                                   10d22s23
            end do                                                      10d22s23
           end do                                                       10d22s23
           if(ikind.eq.npartp)then                                      6d26s24
         orig=tcannon
            tcannon=tcannon+sum                                         6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannonl '),
     $          orig,sum,tcannon
           else                                                         6d26s24
            trace(ikind,1)=trace(ikind,1)+sum                           6d26s24
           end if                                                       6d26s24
          end do                                                        10d22s23
         end if                                                         10d22s23
        end if                                                          10d22s23
       end do                                                           10d22s23
      end do                                                            8d8s23
      end if
      if(ichoice(5).ne.0)then
      write(6,*)('d4o part')
      do is=1,nsdlk                                                     8d9s23
       if(min(irefo(isblk(1,is)),irefo(isblk(2,is)),irefo(isblk(3,is)), 8d9s23
     $      irefo(isblk(4,is))).gt.0)then                               8d9s23
        if(isblk(1,is).eq.isblk(2,is))then                              8d9s23
         nrowd=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2            8d9s23
         nrowi=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                8d9s23
         ncol=(noc(isblk(3,is))*(noc(isblk(3,is))+1))/2                 8d9s23
         jj=j4o(is)                                                     8d9s23
         do ikind=1,npartp                                              6d26s24
          sum=0d0                                                       6d26s24
          do i4=0,irefo(isblk(4,is))-1                                   8d9s23
           i4p=i4+idoubo(isblk(4,is))                                    8d9s23
           do i3=0,i4                                                    8d9s23
            i3p=i3+idoubo(isblk(4,is))                                   8d9s23
            icold=id4o(is)+nrowd*(((i4*(i4+1))/2)+i3)                   8d9s23
            icoli=jj+nrowi*(((i4p*(i4p+1))/2)+i3p)                      8d9s23
            jtri=((i4p*(i4p+1))/2)+i3p
            do i2=0,irefo(isblk(2,is))-1                                8d9s23
             i2p=i2+idoubo(isblk(2,is))                                 8d9s23
             iadd=icold+((i2*(i2+1))/2)                                 8d9s23
             iadi=icoli+((i2p*(i2p+1))/2)+idoubo(isblk(1,is))           8d9s23
             do i1=0,i2                                                 8d9s23
              sum=sum+bc(iadd+i1)*bc(iadi+i1)                           6d26s24
             end do                                                     8d9s23
            end do                                                      8d9s23
           end do                                                       8d9s23
          end do                                                        8d9s23
          if(ikind.eq.npartp)then                                       6d26s24
         orig=tcannon
           tcannon=tcannon+sum                                          6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannonm '),
     $          orig,sum,tcannon
          else                                                          6d26s24
           trace(ikind,1)=trace(ikind,1)+sum                            6d26s24
          end if                                                        6d26s24
          jj=jj+nrowi*ncol                                              8d9s23
         end do                                                         8d9s23
        else                                                            8d9s23
         nrowd=irefo(isblk(1,is))*irefo(isblk(2,is))                    8d9s23
         nrowi=noc(isblk(1,is))*noc(isblk(2,is))                        8d9s23
         ncol=noc(isblk(3,is))*noc(isblk(4,is))                         8d9s23
         jj=j4o(is)                                                     8d9s23
         do ikind=1,npartp                                              6d26s24
          sum=0d0                                                       6d26s24
          do i4=0,irefo(isblk(4,is))-1                                   8d9s23
           i4p=i4+idoubo(isblk(4,is))                                    8d9s23
           do i3=0,irefo(isblk(3,is))-1                                 8d9s23
            i3p=i3+idoubo(isblk(3,is))                                   8d9s23
            icold=id4o(is)+nrowd*(i3+irefo(isblk(3,is))*i4)             8d9s23
            icoli=jj+nrowi*(i3p+noc(isblk(3,is))*i4p)                   8d9s23
            do i2=0,irefo(isblk(2,is))-1                                8d9s23
             i2p=i2+idoubo(isblk(2,is))                                 8d9s23
             iadd=icold+irefo(isblk(1,is))*i2                           8d9s23
             iadi=icoli+noc(isblk(1,is))*i2p+idoubo(isblk(1,is))        8d9s23
             do i1=0,irefo(isblk(1,is))-1                               8d9s23
              orig=trace(3,1)
              sum=sum+bc(iadd+i1)*bc(iadi+i1)                           6d26s24
             end do                                                     8d9s23
            end do                                                      8d9s23
           end do                                                       8d9s23
          end do                                                        8d9s23
          if(ikind.eq.npartp)then                                       6d26s24
         orig=tcannon
           tcannon=tcannon+sum                                          6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannonn '),
     $          orig,sum,tcannon
          else                                                          6d26s24
           trace(ikind,1)=trace(ikind,1)+sum                            6d26s24
          end if                                                        6d26s24
          jj=jj+nrowi*ncol                                              8d9s23
         end do                                                         8d9s23
        end if                                                          8d9s23
       end if                                                           8d9s23
      end do                                                            8d9s23
      end if
      if(ichoice(6).ne.0)then
      write(6,*)('donex part '),(trace(j,1),j=1,npart)                                         8d10s
      zum=0d0
      do is1=1,nsdlk1                                                   8d10s23
       if(isblk1(1,is1).eq.isblk1(2,is1))then                           8d10s23
        nrow1=(irefo(isblk1(1,is1))*(irefo(isblk1(1,is1))+1))/2         8d10s23
        nrowi=(noc(isblk1(1,is1))*(noc(isblk1(1,is1))+1))/2             10d23s23
        nrow3=(nvirt(isblk1(1,is1))*(nvirt(isblk1(1,is1))+1))/2         8d10s23
        isw=0                                                           8d10s23
       else                                                             8d10s23
        nrow1=irefo(isblk1(1,is1))*irefo(isblk1(2,is1))                 8d10s23
        nrowi=noc(isblk1(1,is1))*noc(isblk1(2,is1))                     10d23s23
        nrow3=nvirt(isblk1(1,is1))*nvirt(isblk1(2,is1))                 8d10s23
        isw=1                                                           8d10s23
       end if                                                           8d10s23
       ncol=noc(isblk1(3,is1))*nvirt(isblk1(4,is1))                     10d30s23
       ntot1=nrowi*ncol                                                 8d10s23
       j1x=jonex(is1)                                                   8d10s23
       write(6,*)('jonex for '),(isblk1(j,is1),j=1,4)
       call prntm2(bc(j1x),nrowi,ncol,nrowi)
       do ikind=1,npartp                                                6d26s24
        jj=j1x                                                          8d10s23
        kk=id1x(is1)                                                    8d10s23
        sum=0d0
        do iv=0,nvirt(isblk1(4,is1))-1
         ivp=iv+noc(isblk1(4,is1))                                      10d30s23
         do ia=0,irefo(isblk1(3,is1))-1
          iap=ia+idoubo(isblk1(3,is1))                                  10d30s23
          kcol=ia+irefo(isblk1(3,is1))*iv
          icol0=iap+noc(isblk1(3,is1))*iv                               10d30s23
          jjj=jj+nrowi*icol0                                            10d30s23
          kkk=id1x(is1)+nrow1*kcol                                      10d30s23
          do i4=0,irefo(isblk1(2,is1))-1
           i4p=i4+idoubo(isblk1(2,is1))                                 10d30s23
           i3top=i4+isw*(irefo(isblk1(1,is1))-1-i4)                     8d10s23
           do i3=0,i3top                                                8d10s23
            krec=i3+irefo(isblk1(1,is1))*i4
            ktri=((i4*(i4+1))/2)+i3
            ktri=ktri+isw*(krec-ktri)+kkk                               10d30s23
            i3p=i3+idoubo(isblk1(1,is1))                                10d30s23
            irec=i3p+noc(isblk1(1,is1))*i4p                             10d23s23
            itri=((i4p*(i4p+1))/2)+i3p                                  10d23s23
            itri0=itri
            itri=itri+isw*(irec-itri)+jjj                               10d23s23
            term=bc(itri)*bc(kk+i3)                                     10d30s23
            orig=sum
            sum=sum+term                                                8d10s23
            if(ikind.eq.4)zum=zum+term*ass
           end do
           kk=kk+i3top+1
          end do
         end do
        end do
        jj=jj+nrowi*noc(isblk1(3,is1))*nvirt(isblk1(4,is1))             10d23s23
        if(ikind.eq.npartp)then                                         6d26s24
         orig=tcannon
         tcannon=tcannon+sum                                            6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannono '),
     $          orig,sum,tcannon
        else                                                            6d26s24
         trace(ikind,1)=trace(ikind,1)+sum                               8d11s23
        end if                                                          6d26s24
        j1x=j1x+ntot1                                                   8d10s23
       end do                                                           8d10s23
      end do                                                            8d10s23
      end if
      if(ichoice(7).ne.0)then
                write(6,*)('what is under goal '),bc(igoal)
      write(6,*)('dj part')                                                 8d11s23
      do is=1,nsdlk                                                     8d11s23
       zum=0d0
       if(isblk(1,is).eq.isblk(2,is))then                               8d11s23
        nrow=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                  10d23s23
        nrowr=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2                  10d23s23
        isw=0                                                           8d11s23
       else                                                             8d11s23
        nrow=noc(isblk(1,is))*noc(isblk(2,is))                          10d23s23
        nrowr=irefo(isblk(1,is))*irefo(isblk(2,is))                          10d23s23
        isw=1                                                           8d11s23
       end if
       ncol=nvirt(isblk(3,is))*nvirt(isblk(4,is))                       8d11s23
       ntot=nrow*ncol                                                   8d11s23
       jj=jmatd(is)                                                     8d11s23
       write(6,*)('for sym '),(isblk(j,is),j=1,4)
       write(6,*)('density ')
       call prntm2(bc(jmden(is)),ncol,nrowr,ncol)
       do ikind=1,npartp                                                6d26s24
        sum=0d0                                                         8d11s23
        jjj=jj                                                          8d11s23
        if(ikind.eq.3)then
         write(6,*)('transformed integrals ')
         nvv=nvirt(isblk(3,is))*nvirt(isblk(4,is))
         call prntm2(bc(jjj),nrow,nvv,nrow)
        end if
        do i4=0,nvirt(isblk(4,is))-1                                    8d11s23
         do i3=0,nvirt(isblk(3,is))-1                                   8d11s23
          kkk=jmden(is)+i3+nvirt(isblk(3,is))*i4                        8d11s23
          do i2=0,irefo(isblk(2,is))-1                                  8d11s23
           i2p=i2+idoubo(isblk(2,is))                                   10d23s23
           i1top=i2+isw*(irefo(isblk(1,is))-1-i2)                       8d11s23
           do i1=0,i1top                                                8d11s23
            i1p=i1+idoubo(isblk(1,is))                                  10d23s23
            ix=max(i1p,i2p)                                             10d23s23
            in=min(i1p,i2p)                                             10d23s23
            itri=((ix*(ix+1))/2)+in                                     10d23s23
            irec=i1p+noc(isblk(1,is))*i2p                               11d1s23
            itri=itri+isw*(irec-itri)+jjj
            term=bc(itri)*bc(kkk)                                       10d23s23
            sum=sum+term                                                8d11s23
            xjt2(is)=xjt2(is)+term
            kkk=kkk+ncol
           end do                                                       8d11s23
          end do                                                        8d11s23
          jjj=jjj+nrow                                                  10d23s23
         end do                                                         8d11s23
        end do                                                          8d11s23
        if(ikind.eq.npartp)then                                         6d26s24
         orig=tcannon
         tcannon=tcannon+sum                                            6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannonp '),
     $          orig,sum,tcannon
        else                                                            6d26s24
         trace(ikind,1)=trace(ikind,1)+sum                               8d11s23
        end if                                                          6d26s24
        jj=jj+ntot                                                      8d11s23
       end do                                                           8d11s23
      end do                                                            8d11s23
      end if
      if(ichoice(8).ne.0)then                                           10d18s23
      write(6,*)('dk part')                                                 8d11s23
      do is=1,nsdlkk                                                     8d11s23
       nrow=noc(isblkk(1,is))*noc(isblkk(2,is))                         10d23s23
       nrowr=irefo(isblkk(1,is))*irefo(isblkk(2,is))                    11d6s23
       ncol=nvirt(isblkk(3,is))*nvirt(isblkk(4,is))                     8d14s23
       ntot=nrow*ncol                                                   8d11s23
       kk=kmatd(is)                                                     8d11s23
       write(6,*)('for syms '),(isblkk(j,is),j=1,4)
       write(6,*)('density ')
       call prntm2(bc(kmden(is)),ncol,nrowr,ncol)
       do ikind=1,npartp                                                6d26s24
        if(ikind.eq.1)then
         write(6,*)('integrals ')
         call prntm2(bc(kk),nrow,ncol,nrow)
        end if
        sum=0d0                                                         8d11s23
        kkk=kk                                                          8d11s23
        do i4=0,nvirt(isblkk(4,is))-1                                    8d11s23
         i4p=i4+noc(isblkk(4,is))+1
         do i3=0,nvirt(isblkk(3,is))-1                                   8d11s23
          i3p=i3+noc(isblkk(3,is))+1
          irow=i3+nvirt(isblkk(3,is))*i4+1
          jjj=kmden(is)+i3+nvirt(isblkk(3,is))*i4                        8d11s23
          icol=0
          do i2=0,irefo(isblkk(2,is))-1                                  8d11s23
           i2p=i2+idoubo(isblkk(2,is))                                  11d7s23
           kkkk=kkk+idoubo(isblkk(1,is))+noc(isblkk(1,is))*i2p          11d7s23
           do i1=0,irefo(isblkk(1,is))-1                                8d14s23
            icol=icol+1
            term=bc(kkkk+i1)*bc(jjj)                                    10d23s23
            sum=sum+term                                                8d11s23
            jjj=jjj+ncol
           end do                                                       8d11s23
          end do                                                        8d11s23
          kkk=kkk+nrow                                                  11d6s23
         end do                                                         8d11s23
        end do                                                          8d11s23
        if(ikind.eq.npartp)then                                         6d26s24
         orig=tcannon
         tcannon=tcannon+sum                                            6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannonq '),
     $          orig,sum,tcannon
        else                                                            6d26s24
         trace(ikind,1)=trace(ikind,1)+sum                               8d11s23
        end if                                                          6d26s24
        kk=kk+ntot                                                      8d14s23
       end do                                                           8d11s23
      end do                                                            8d11s23
      end if
      if(ichoice(9).ne.0)then                                           11d14s23
       write(6,*)('d3x part ')                                          11d14s23
       zumt=0d0
       do is=1,nsdlk1                                                    7d28s22
        if(isblk1(1,is).eq.isblk1(2,is))then                             7d28s22
         nrow3=(nvirt(isblk1(1,is))*(nvirt(isblk1(1,is))+1))/2           7d28s22
        else                                                             7d28s22
         nrow3=nvirt(isblk1(1,is))*nvirt(isblk1(2,is))                   7d28s22
        end if                                                           7d28s22
        call ilimts(irefo(isblk1(3,is)),nvirt(isblk1(4,is)),mynprocg,    7d28s22
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            7d28s22
        ncol=ih+1-il                                                     7d28s22
        if(min(nrow3,ncol).gt.0)then                                     10d18s23
         write(6,*)('for symmetries '),(isblk1(j,is),j=1,4)
         j3xd=i3xd(is)                                                  11d14s23
         do ikind=1,npartp                                              6d26s24
          write(6,*)('for kind '),ikind
          zum=0d0
          sum=0d0                                                       6d26s24
          do i=0,nrow3*ncol-1                                             10d18s23
           term=bc(id3x(is)+i)*bc(j3xd+i)                               11d14s23
           sum=sum+term                                                 6d26s24
           zum=zum+term*ass
           zumt=zumt+term*ass
           if(abs(term).gt.1d-10.and.ikind.eq.2)
     $          write(6,*)bc(id3x(is)+i),
     $          bc(j3xd+i)*ass,zum,zumt,j3xd+i
          end do                                                        11d14s23
          if(abs(zum).gt.1d-10.and.ikind.eq.2)write(6,*)('>> '),
     $         (isblk1(j,is),j=1,4),zum,zumt,trace(ikind,1)
          j3xd=j3xd+nrow3*ncol                                          11d14s23
          if(ikind.eq.npartp)then                                       6d26s24
         orig=tcannon
           tcannon=tcannon+sum                                          6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannonr '),
     $          orig,sum,tcannon
          else                                                          6d26s24
           trace(ikind,1)=trace(ikind,1)+sum                            6d26s24
          end if                                                        6d26s24
         end do
        end if                                                           10d18s23
       end do                                                            7d28s22
       if(ichoice(11).ne.0)then                                         1d9s24
        trx=0d0
        trx2=0d0
        ioff=1
        do isb=1,nsymb                                                  1d9s24
         nn=nbasdws(isb)**2
         if(min(nvirt(isb),irefo(isb)).gt.0)then                        1d10s24
          write(6,*)('for symmetry '),isb
          iuu=ibcoff                                                     1d9s24
          ibcoff=iuu+nvirt(isb)*irefo(isb)                               1d9s24
          call enough('iuu.denmx12',bc,ibc)                              1d9s24
          do iz=iuu,ibcoff-1                                             1d9s24
           bc(iz)=0d0                                                    1d9s24
          end do                                                         1d9s24
          do is3=1,nsdlk1                                                1d9s24
           if(isblk1(3,is3).eq.isb)then                                  1d9s24
            write(6,*)('3x density symmetry '),(isblk1(j,is3),j=1,4)
            if(isblk1(1,is3).eq.isblk1(2,is3))then                       1d9s24
             isw=0                                                       1d9s24
             nrow=(nvirt(isblk1(1,is3))*(nvirt(isblk1(1,is3))+1))/2      1d9s24
            else                                                         1d9s24
             isw=1                                                       1d9s24
             nrow=nvirt(isblk1(1,is3))*nvirt(isblk1(2,is3))              1d9s24
            end if                                                       1d9s24
            look=2
            zum=0d0
            do i2=0,nvirt(isblk1(4,is3))-1                               1d10s24
             do i1=0,irefo(isblk1(3,is3))-1                              1d10s24
              idad=id3x(is3)+nrow*(i1+irefo(isblk1(3,is3))*i2)           1d10s24
              do iv=0,nvirt(isb)-1                                       1d10s24
               juu=iuu+iv+nvirt(isb)*i1                                 1d11s24
               do i4=0,nvirt(isblk1(2,is3))-1                            1d10s24
                i3top=i4+isw*(nvirt(isblk1(1,is3))-1-i4)                 1d10s24
                do i3=0,i3top                                            1d10s24
                 itri=((i4*(i4+1))/2)+i3                                 1d10s24
                 irec=i3+nvirt(isblk1(1,is3))*i4                         1d10s24
                 itri=itri+isw*(irec-itri)+idad                          1d10s24
                 xinta=get4int(i4x,isblk1(1,is3),isblk1(2,is3),isb,     1d10s24
     $               isblk1(4,is3),i3,i4,iv,i2,nvirt,idum,bc,ibc)       1d10s24
                 orig=bc(look)
                 bc(juu)=bc(juu)+xinta*bc(itri)                          1d10s24
                 if(abs(orig-bc(look)).gt.1d-14)then
                  zum=zum+xinta*bc(itri)
                  write(6,*)
     $                orig,bc(itri),xinta,bc(look),zum,i3,i4,iv,i2
                 end if
                end do                                                   1d10s24
               end do                                                    1d10s24
              end do                                                     1d10s24
             end do                                                      1d10s24
            end do                                                       1d10s24
           end if                                                        1d9s24
          end do                                                         1d9s24
          write(6,*)('what we have for uu matrix: '),iuu
          call prntm2(bc(iuu),nvirt(isb),irefo(isb),nvirt(isb))
          tcur=0d0
          do i2=0,irefo(isb)-1
           i2p=i2+idoubo(isb)
           do i1=0,nvirt(isb)-1
            i1p=i1+noc(isb)
            iad1=iuu+i1+nvirt(isb)*i2
            iad2=ioff+i1p+nbasdws(isb)*i2p
            if(abs(darot(iad2)).gt.1d-10)write(6,*)i1,i2,darot(iad2),
     $           bc(iad1),iad1
            tcur=tcur+darot(iad2)*bc(iad1)
           end do
          end do
          trx=trx+tcur
          write(6,*)('trace with darot: '),tcur,trx
          if(ipass.ge.0)then
           iip=ioff+nn
           tcur=0d0
           do i2=0,irefo(isb)-1
            i2p=i2+idoubo(isb)
            do i1=0,nvirt(isb)-1
             i1p=i1+noc(isb)
             iad1=iuu+i1+nvirt(isb)*i2
             iad2=iip+i1p+nbasdws(isb)*i2p
             if(abs(darot(iad2)).gt.1d-10)write(6,*)i1,i2,darot(iad2),
     $           bc(iad1),iad1
             tcur=tcur+darot(iad2)*bc(iad1)
            end do
           end do
           trx2=trx2+tcur
           write(6,*)('trace with transder: '),tcur,trx2
          end if
          ibcoff=iuu                                                     1d9s24
         end if
         ioff=ioff+nn
         if(ipass.ge.0)ioff=ioff+nn
        end do                                                          1d9s24
       end if
      end if                                                            11d14s23
      if(ichoice(10).ne.0.and.ndoub.gt.0)then                           4d17s24
       write(6,*)('d4x part '),bc(igoal)                                          11d17s23
       do ikind=1,npartp                                                6d26s24
        write(6,*)('for kind '),ikind,bc(igoal)
        if(ndoub.gt.0)then                                              6d12s24
         write(6,*)('calling chc4v for kind '),ikind
         call chc4v(nsymb,multha,isymmrci,vdinout,1,nfdat,nvirt,sr2,      11d17s23
     $      bc,ibc,ndoub,dotx,i4xd(1,ikind),1)                            11d17s23
         write(6,*)('4v part: '),dotx                                    11d17s23
         if(ikind.eq.npartp)then                                        6d26s24
         orig=tcannon
          tcannon=tcannon+dotx                                          6d26s24
           if(abs(orig-tcannon).gt.1d-12)write(6,*)('tcannons '),
     $          orig,sum,tcannon
         else                                                           6d26s24
          trace(ikind,1)=trace(ikind,1)+dotx                              11d17s23
         end if                                                         6d26s24
        end if                                                          6d12s24
       end do                                                           11d17s23
      end if                                                            11d17s23
      write(6,*)('final traces ')
      call prntm2(trace,npart,nsymb,4)
      if(nsymb.gt.1)then
       write(6,*)('summed over symmetry blocks: ')
       gsum=0d0
       do ikind=1,npart
        sum=trace(ikind,1)
        do i=2,nsymb
         sum=sum+trace(ikind,i)
        end do
        write(6,*)sum
        gsum=gsum+sum
       end do
      else
       gsum=trace(1,1)+trace(1,2)+trace(3,1)
      end if
      write(6,*)('grand total: '),gsum
      write(6,*)('tcannon: '),tcannon
      if(npart.eq.4)write(6,*)('w/o tt part: '),gsum-sum
      ibcoff=ibcoffo                                                    8d4s23
      return                                                            8d4s23
      end                                                               8d4s23
          subroutine e0xd(tmp,isa,isb,isc,isd,nbasdws,idoubo,irefo,     5d2s24
     $         isblk,nsdlk,e0x,id4o,bc,ibc,iuseage)                     5d6s24
          implicit real*8 (a-h,o-z)
          include "common.store"
          dimension tmp(*),nbasdws(*),idoubo(*),irefo(*),isblk(4,*),    5d2s24
     $         id4o(*),iuseage(8)                                       5d6s24
c     isb gt isa
c     isd gt isc
c     icd ge iab
          iab=((isb*(isb-1))/2)+isa                                     5d2s24
          icd=((isd*(isd-1))/2)+isc                                     5d2s24
          do is=1,nsdlk
           ix=max(isblk(1,is),isblk(2,is))                              5d2s24
           in=min(isblk(1,is),isblk(2,is))                              5d2s24
           i12=((ix*(ix-1))/2)+in                                       5d2s24
           ix=max(isblk(3,is),isblk(4,is))                              5d2s24
           in=min(isblk(3,is),isblk(4,is))                              5d2s24
           i34=((ix*(ix-1))/2)+in                                       5d2s24
           ix=max(i12,i34)
           in=min(i12,i34)
           if(ix.eq.icd.and.in.eq.iab)then
            if(isblk(1,is).eq.isblk(2,is))then
             nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2         5d2s24
             ncol=(irefo(isblk(3,is))*(irefo(isblk(4,is))+1))/2         5d2s24
             isw=0
            else
             nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                 5d2s24
             nrow=irefo(isblk(3,is))*irefo(isblk(4,is))                 5d2s24
             isw=1
            end if
            if(i12.eq.iab)then                                          5d2s24
             if(isblk(1,is).eq.isa)then                                 5d2s24
              if(isblk(3,is).eq.isc)then                                5d2s24
               iuseage(1)=iuseage(1)+1                                  5d6s24
               do i4=0,irefo(isblk(4,is))-1                             5d2s24
                i3top=i4+isw*(irefo(isblk(3,is))-1-i4)                  5d2s24
                i4p=i4+idoubo(isblk(4,is))                              5d2s24
                do i3=0,i3top                                           5d2s24
                 i3p=i3+idoubo(isblk(3,is))
                 irec34=i3+irefo(isblk(3,is))*i4                        5d2s24
                 itri34=((i4*(i4+1))/2)+i3                              5d2s24
                 itri34=itri34+isw*(irec34-itri34)                      5d2s24
                 do i2=0,irefo(isblk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk(2,is))
                  i1top=i2+isw*(irefo(isblk(1,is))-1-i2)                5d2s24
                  do i1=0,i1top                                         5d2s24
                   i1p=i1+idoubo(isblk(1,is))
                   irec12=i1+irefo(isblk(1,is))*i2
                   itri12=((i2*(i2+1))/2)+i1                            5d2s24
                   itri12=itri12+isw*(irec12-itri12)
                   iad1=id4o(is)+itri12+nrow*itri34                     5d2s24
                   iad2=1+i1p+nbasdws(isa)*(i2p+nbasdws(isb)*(i3p        5d2s24
     $                 +nbasdws(isc)*i4p))                               5d2s24
                   e0x=e0x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(2)=iuseage(2)+1                                  5d6s24
               do i4=0,irefo(isblk(4,is))-1                             5d2s24
                i4p=i4+idoubo(isblk(4,is))                              5d2s24
                do i3=0,irefo(isblk(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk(3,is))
                 do i2=0,irefo(isblk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk(2,is))
                  do i1=0,irefo(isblk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk(1,is))
                   iad1=id4o(is)+i1+irefo(isblk(1,is))*(i2              5d2s24
     $                  +irefo(isblk(2,is))*(i3+irefo(isblk(3,is))*i4)) 5d2s24
                   iad2=1+i1p+nbasdws(isa)*(i2p+nbasdws(isb)*(i4p        5d2s24
     $                 +nbasdws(isc)*i3p))                               5d2s24
                   e0x=e0x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             else                                                       5d2s24
              if(isblk(3,is).eq.isc)then                                5d2s24
               iuseage(3)=iuseage(3)+1                                  5d6s24
               do i4=0,irefo(isblk(4,is))-1                             5d2s24
                i4p=i4+idoubo(isblk(4,is))                              5d2s24
                do i3=0,irefo(isblk(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk(3,is))
                 do i2=0,irefo(isblk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk(2,is))
                  do i1=0,irefo(isblk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk(1,is))
                   iad1=id4o(is)+i1+irefo(isblk(1,is))*(i2              5d2s24
     $                  +irefo(isblk(2,is))*(i3+irefo(isblk(3,is))*i4)) 5d2s24
                   iad2=1+i2p+nbasdws(isa)*(i1p+nbasdws(isb)*(i3p        5d2s24
     $                 +nbasdws(isc)*i4p))                               5d2s24
                   e0x=e0x+bc(iad1)*tmp(iad2)                           5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(4)=iuseage(4)+1                                  5d6s24
               do i4=0,irefo(isblk(4,is))-1                             5d2s24
                i4p=i4+idoubo(isblk(4,is))                              5d2s24
                do i3=0,irefo(isblk(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk(3,is))
                 do i2=0,irefo(isblk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk(2,is))
                  do i1=0,irefo(isblk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk(1,is))
                   iad1=id4o(is)+i1+irefo(isblk(1,is))*(i2              5d2s24
     $                  +irefo(isblk(2,is))*(i3+irefo(isblk(3,is))*i4)) 5d2s24
                   iad2=1+i2p+nbasdws(isa)*(i1p+nbasdws(isb)*(i4p        5d2s24
     $                 +nbasdws(isc)*i3p))                               5d2s24
                   e0x=e0x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             end if                                                     5d2s24
            else                                                        5d2s24
             if(isblk(3,is).eq.isa)then                                 5d2s24
              if(isblk(1,is).eq.isc)then                                5d2s24
               iuseage(5)=iuseage(5)+1                                  5d6s24
               do i4=0,irefo(isblk(4,is))-1                             5d2s24
                i3top=i4+isw*(irefo(isblk(3,is))-1-i4)                  5d2s24
                i4p=i4+idoubo(isblk(4,is))                              5d2s24
                do i3=0,i3top                                           5d2s24
                 i3p=i3+idoubo(isblk(3,is))
                 irec34=i3+irefo(isblk(3,is))*i4                        5d2s24
                 itri34=((i4*(i4+1))/2)+i3                              5d2s24
                 itri34=itri34+isw*(irec34-itri34)                      5d2s24
                 do i2=0,irefo(isblk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk(2,is))
                  i1top=i2+isw*(irefo(isblk(1,is))-1-i2)                5d2s24
                  do i1=0,i1top                                         5d2s24
                   i1p=i1+idoubo(isblk(1,is))
                   irec12=i1+irefo(isblk(1,is))*i2
                   itri12=((i2*(i2+1))/2)+i1                            5d2s24
                   itri12=itri12+isw*(irec12-itri12)
                   iad1=id4o(is)+itri12+nrow*itri34                     5d2s24
                   iad2=1+i3p+nbasdws(isa)*(i4p+nbasdws(isb)*(i1p        5d2s24
     $                 +nbasdws(isc)*i2p))                               5d2s24
                   e0x=e0x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(6)=iuseage(6)+1                                  5d6s24
               do i4=0,irefo(isblk(4,is))-1                             5d2s24
                i4p=i4+idoubo(isblk(4,is))                              5d2s24
                do i3=0,irefo(isblk(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk(3,is))
                 do i2=0,irefo(isblk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk(2,is))
                  do i1=0,irefo(isblk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk(1,is))
                   iad1=id4o(is)+i1+irefo(isblk(1,is))*(i2              5d2s24
     $                  +irefo(isblk(2,is))*(i3+irefo(isblk(3,is))*i4)) 5d2s24
                   iad2=1+i3p+nbasdws(isa)*(i4p+nbasdws(isb)*(i2p        5d2s24
     $                 +nbasdws(isc)*i1p))                               5d2s24
                   e0x=e0x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             else                                                       5d2s24
              if(isblk(1,is).eq.isc)then                                5d2s24
               iuseage(7)=iuseage(7)+1                                  5d6s24
               do i4=0,irefo(isblk(4,is))-1                             5d2s24
                i4p=i4+idoubo(isblk(4,is))                              5d2s24
                do i3=0,irefo(isblk(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk(3,is))
                 do i2=0,irefo(isblk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk(2,is))
                  do i1=0,irefo(isblk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk(1,is))
                   iad1=id4o(is)+i1+irefo(isblk(1,is))*(i2              5d2s24
     $                  +irefo(isblk(2,is))*(i3+irefo(isblk(3,is))*i4)) 5d2s24
                   iad2=1+i4p+nbasdws(isa)*(i3p+nbasdws(isb)*(i1p        5d2s24
     $                 +nbasdws(isc)*i2p))                               5d2s24
                   e0x=e0x+bc(iad1)*tmp(iad2)                           5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(8)=iuseage(8)+1                                  5d6s24
               do i4=0,irefo(isblk(4,is))-1                             5d2s24
                i4p=i4+idoubo(isblk(4,is))                              5d2s24
                do i3=0,irefo(isblk(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk(3,is))
                 do i2=0,irefo(isblk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk(2,is))
                  do i1=0,irefo(isblk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk(1,is))
                   iad1=id4o(is)+i1+irefo(isblk(1,is))*(i2              5d2s24
     $                  +irefo(isblk(2,is))*(i3+irefo(isblk(3,is))*i4)) 5d2s24
                   iad2=1+i4p+nbasdws(isa)*(i3p+nbasdws(isb)*(i2p        5d2s24
     $                 +nbasdws(isc)*i1p))                               5d2s24
                   e0x=e0x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             end if                                                     5d2s24
            end if                                                      5d2s24
           end if
          end do                                                        5d2s24
          return                                                        5d2s24
          end                                                           5d2s24
          subroutine e1xd(tmp,isa,isb,isc,isd,nbasdws,idoubo,irefo,     5d2s24
     $         nvirt,noc,isblk1,nsdlk1,e1x,id1x,bc,ibc,iuseage)         5d6s24
          implicit real*8 (a-h,o-z)
          include "common.store"
          dimension tmp(*),nbasdws(*),idoubo(*),irefo(*),isblk1(4,*),    5d2s24
     $         id1x(*),nvirt(*),noc(*),iuseage(*)                       5d6s24
c     isb gt isa
c     isd gt isc
c     icd ge iab
          iab=((isb*(isb-1))/2)+isa                                     5d2s24
          icd=((isd*(isd-1))/2)+isc                                     5d2s24
          do is=1,nsdlk1
           ix=max(isblk1(1,is),isblk1(2,is))                              5d2s24
           in=min(isblk1(1,is),isblk1(2,is))                              5d2s24
           i12=((ix*(ix-1))/2)+in                                       5d2s24
           ix=max(isblk1(3,is),isblk1(4,is))                              5d2s24
           in=min(isblk1(3,is),isblk1(4,is))                              5d2s24
           i34=((ix*(ix-1))/2)+in                                       5d2s24
           ix=max(i12,i34)
           in=min(i12,i34)
           if(ix.eq.icd.and.in.eq.iab)then
            if(isblk1(1,is).eq.isblk1(2,is))then
             nrow=(irefo(isblk1(1,is))*(irefo(isblk1(1,is))+1))/2         5d2s24
             isw=0
            else
             nrow=irefo(isblk1(1,is))*irefo(isblk1(2,is))                 5d2s24
             isw=1
            end if
            if(i12.eq.iab)then                                          5d2s24
             if(isblk1(1,is).eq.isa)then                                 5d2s24
              if(isblk1(3,is).eq.isc)then                                5d2s24
               iuseage(1)=iuseage(1)+1                                  5d6s24
               do i4=0,nvirt(isblk1(4,is))-1                            5d2s24
                i4p=i4+noc(isblk1(4,is))                                5d2s24
                do i3=0,irefo(isblk1(3,is))-1                           5d2s24
                 i3p=i3+idoubo(isblk1(3,is))
                 do i2=0,irefo(isblk1(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk1(2,is))
                  i1top=i2+isw*(irefo(isblk1(1,is))-1-i2)                5d2s24
                  do i1=0,i1top                                         5d2s24
                   i1p=i1+idoubo(isblk1(1,is))
                   irec12=i1+irefo(isblk1(1,is))*i2
                   itri12=((i2*(i2+1))/2)+i1                            5d2s24
                   itri12=itri12+isw*(irec12-itri12)
                   iad1=id1x(is)+itri12+nrow*(i3+irefo(isblk1(3,is))*i4)5d2s24
                   iad2=1+i1p+nbasdws(isa)*(i2p+nbasdws(isb)*(i3p        5d2s24
     $                 +nbasdws(isc)*i4p))                               5d2s24
                   e1x=e1x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(2)=iuseage(2)+1                                  5d6s24
               do i4=0,nvirt(isblk1(4,is))-1                            5d2s24
                i4p=i4+noc(isblk1(4,is))                                5d2s24
                do i3=0,irefo(isblk1(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk1(3,is))
                 do i2=0,irefo(isblk1(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk1(2,is))
                  do i1=0,irefo(isblk1(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk1(1,is))
                   iad1=id1x(is)+i1+irefo(isblk1(1,is))*(i2              5d2s24
     $                 +irefo(isblk1(2,is))*(i3+irefo(isblk1(3,is))*i4)) 5d2s24
                   iad2=1+i1p+nbasdws(isa)*(i2p+nbasdws(isb)*(i4p        5d2s24
     $                 +nbasdws(isc)*i3p))                               5d2s24
                   e1x=e1x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             else                                                       5d2s24
              if(isblk1(3,is).eq.isc)then                                5d2s24
               iuseage(3)=iuseage(3)+1                                  5d6s24
               do i4=0,nvirt(isblk1(4,is))-1                            5d2s24
                i4p=i4+noc(isblk1(4,is))                                5d2s24
                do i3=0,irefo(isblk1(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk1(3,is))
                 do i2=0,irefo(isblk1(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk1(2,is))
                  do i1=0,irefo(isblk1(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk1(1,is))
                   iad1=id1x(is)+i1+irefo(isblk1(1,is))*(i2              5d2s24
     $                 +irefo(isblk1(2,is))*(i3+irefo(isblk1(3,is))*i4)) 5d2s24
                   iad2=1+i2p+nbasdws(isa)*(i1p+nbasdws(isb)*(i3p        5d2s24
     $                 +nbasdws(isc)*i4p))                               5d2s24
                   e1x=e1x+bc(iad1)*tmp(iad2)                           5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(4)=iuseage(4)+1                                  5d6s24
               do i4=0,nvirt(isblk1(4,is))-1                             5d2s24
                i4p=i4+noc(isblk1(4,is))                                5d2s24
                do i3=0,irefo(isblk1(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk1(3,is))
                 do i2=0,irefo(isblk1(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk1(2,is))
                  do i1=0,irefo(isblk1(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk1(1,is))
                   iad1=id1x(is)+i1+irefo(isblk1(1,is))*(i2              5d2s24
     $                 +irefo(isblk1(2,is))*(i3+irefo(isblk1(3,is))*i4)) 5d2s24
                   iad2=1+i2p+nbasdws(isa)*(i1p+nbasdws(isb)*(i4p        5d2s24
     $                 +nbasdws(isc)*i3p))                               5d2s24
                   e1x=e1x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             end if                                                     5d2s24
            else                                                        5d2s24
             if(isblk1(3,is).eq.isa)then                                 5d2s24
              if(isblk1(1,is).eq.isc)then                                5d2s24
               iuseage(5)=iuseage(5)+1                                  5d6s24
               do i4=0,nvirt(isblk1(4,is))-1                             5d2s24
                i4p=i4+noc(isblk1(4,is))                                5d2s24
                do i3=0,irefo(isblk1(3,is))-1                           5d2s24
                 i3p=i3+idoubo(isblk1(3,is))
                 do i2=0,irefo(isblk1(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk1(2,is))
                  i1top=i2+isw*(irefo(isblk1(1,is))-1-i2)                5d2s24
                  do i1=0,i1top                                         5d2s24
                   i1p=i1+idoubo(isblk1(1,is))
                   irec12=i1+irefo(isblk1(1,is))*i2
                   itri12=((i2*(i2+1))/2)+i1                            5d2s24
                   itri12=itri12+isw*(irec12-itri12)
                   iad1=id1x(is)+itri12+nrow*(i3+irefo(isblk1(3,is))*i4)5d2s24
                   iad2=1+i3p+nbasdws(isa)*(i4p+nbasdws(isb)*(i1p        5d2s24
     $                 +nbasdws(isc)*i2p))                               5d2s24
                   e1x=e1x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(6)=iuseage(6)+1                                  5d6s24
               do i4=0,nvirt(isblk1(4,is))-1                             5d2s24
                i4p=i4+noc(isblk1(4,is))                                5d2s24
                do i3=0,irefo(isblk1(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk1(3,is))
                 do i2=0,irefo(isblk1(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk1(2,is))
                  do i1=0,irefo(isblk1(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk1(1,is))
                   iad1=id1x(is)+i1+irefo(isblk1(1,is))*(i2              5d2s24
     $                 +irefo(isblk1(2,is))*(i3+irefo(isblk1(3,is))*i4)) 5d2s24
                   iad2=1+i3p+nbasdws(isa)*(i4p+nbasdws(isb)*(i2p        5d2s24
     $                 +nbasdws(isc)*i1p))                               5d2s24
                   e1x=e1x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             else                                                       5d2s24
              if(isblk1(1,is).eq.isc)then                                5d2s24
               iuseage(7)=iuseage(7)+1                                  5d6s24
               do i4=0,nvirt(isblk1(4,is))-1                            5d2s24
                i4p=i4+noc(isblk1(4,is))                                5d2s24
                do i3=0,irefo(isblk1(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk1(3,is))
                 do i2=0,irefo(isblk1(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk1(2,is))
                  do i1=0,irefo(isblk1(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk1(1,is))
                   iad1=id1x(is)+i1+irefo(isblk1(1,is))*(i2              5d2s24
     $                 +irefo(isblk1(2,is))*(i3+irefo(isblk1(3,is))*i4)) 5d2s24
                   iad2=1+i4p+nbasdws(isa)*(i3p+nbasdws(isb)*(i1p        5d2s24
     $                 +nbasdws(isc)*i2p))                               5d2s24
                   e1x=e1x+bc(iad1)*tmp(iad2)                           5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(8)=iuseage(8)+1                                  5d6s24
               do i4=0,nvirt(isblk1(4,is))-1                             5d2s24
                i4p=i4+noc(isblk1(4,is))                                5d2s24
                do i3=0,irefo(isblk1(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk1(3,is))
                 do i2=0,irefo(isblk1(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk1(2,is))
                  do i1=0,irefo(isblk1(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk1(1,is))
                   iad1=id1x(is)+i1+irefo(isblk1(1,is))*(i2              5d2s24
     $                 +irefo(isblk1(2,is))*(i3+irefo(isblk1(3,is))*i4)) 5d2s24
                   iad2=1+i4p+nbasdws(isa)*(i3p+nbasdws(isb)*(i2p        5d2s24
     $                 +nbasdws(isc)*i1p))                               5d2s24
                   e1x=e1x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             end if                                                     5d2s24
            end if                                                      5d2s24
           end if
          end do                                                        5d2s24
          return                                                        5d2s24
          end                                                           5d2s24
          subroutine e3xd(tmp,isa,isb,isc,isd,nbasdws,idoubo,irefo,     5d2s24
     $         nvirt,noc,isblk1,nsdlk1,e1x,id3x,bc,ibc,iuseage)         5d6s24
          implicit real*8 (a-h,o-z)
          include "common.store"
          dimension tmp(*),nbasdws(*),idoubo(*),irefo(*),isblk1(4,*),    5d2s24
     $         id3x(*),nvirt(*),noc(*),iuseage(*)                       5d6s24
c     isb gt isa
c     isd gt isc
c     icd ge iab
          iab=((isb*(isb-1))/2)+isa                                     5d2s24
          icd=((isd*(isd-1))/2)+isc                                     5d2s24
          do is=1,nsdlk1
           ix=max(isblk1(1,is),isblk1(2,is))                              5d2s24
           in=min(isblk1(1,is),isblk1(2,is))                              5d2s24
           i12=((ix*(ix-1))/2)+in                                       5d2s24
           ix=max(isblk1(3,is),isblk1(4,is))                              5d2s24
           in=min(isblk1(3,is),isblk1(4,is))                              5d2s24
           i34=((ix*(ix-1))/2)+in                                       5d2s24
           ix=max(i12,i34)
           in=min(i12,i34)
           if(ix.eq.icd.and.in.eq.iab)then
            if(isblk1(1,is).eq.isblk1(2,is))then
             nrow=(nvirt(isblk1(1,is))*(nvirt(isblk1(1,is))+1))/2       5d10s24
             isw=0
            else
             nrow=nvirt(isblk1(1,is))*nvirt(isblk1(2,is))               5d10s24
             isw=1
            end if
            if(i12.eq.iab)then                                          5d2s24
             if(isblk1(1,is).eq.isa)then                                 5d2s24
              if(isblk1(3,is).eq.isc)then                                5d2s24
               iuseage(1)=iuseage(1)+1                                  5d6s24
               do i4=0,nvirt(isblk1(4,is))-1                            5d2s24
                i4p=i4+noc(isblk1(4,is))                                5d2s24
                do i3=0,irefo(isblk1(3,is))-1                           5d2s24
                 i3p=i3+idoubo(isblk1(3,is))
                 do i2=0,nvirt(isblk1(2,is))-1                           5d2s24
                  i2p=i2+noc(isblk1(2,is))                              5d10s24
                  i1top=i2+isw*(nvirt(isblk1(1,is))-1-i2)               5d10s24
                  do i1=0,i1top                                         5d2s24
                   i1p=i1+noc(isblk1(1,is))                             5d10s24
                   irec12=i1+nvirt(isblk1(1,is))*i2                     5d10s24
                   itri12=((i2*(i2+1))/2)+i1                            5d2s24
                   itri12=itri12+isw*(irec12-itri12)
                   iad1=id3x(is)+itri12+nrow*(i3+irefo(isblk1(3,is))*i4)5d2s24
                   iad2=1+i1p+nbasdws(isa)*(i2p+nbasdws(isb)*(i3p        5d2s24
     $                 +nbasdws(isc)*i4p))                               5d2s24
                   e1x=e1x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(2)=iuseage(2)+1                                  5d6s24
               do i4=0,nvirt(isblk1(4,is))-1                            5d2s24
                i4p=i4+noc(isblk1(4,is))                                5d2s24
                do i3=0,irefo(isblk1(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk1(3,is))
                 do i2=0,nvirt(isblk1(2,is))-1                          5d10s24
                  i2p=i2+noc(isblk1(2,is))                              5d10s24
                  do i1=0,nvirt(isblk1(1,is))-1                         5d10s24
                   i1p=i1+noc(isblk1(1,is))                             5d10s24
                   iad1=id3x(is)+i1+nvirt(isblk1(1,is))*(i2             5d10s24
     $                 +nvirt(isblk1(2,is))*(i3+irefo(isblk1(3,is))*i4))5d10s24
                   iad2=1+i1p+nbasdws(isa)*(i2p+nbasdws(isb)*(i4p        5d2s24
     $                 +nbasdws(isc)*i3p))                               5d2s24
                   e1x=e1x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             else                                                       5d2s24
              if(isblk1(3,is).eq.isc)then                                5d2s24
               iuseage(3)=iuseage(3)+1                                  5d6s24
               do i4=0,nvirt(isblk1(4,is))-1                            5d2s24
                i4p=i4+noc(isblk1(4,is))                                5d2s24
                do i3=0,irefo(isblk1(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk1(3,is))
                 do i2=0,nvirt(isblk1(2,is))-1                          5d10s24
                  i2p=i2+noc(isblk1(2,is))                              5d10s24
                  do i1=0,nvirt(isblk1(1,is))-1                         5d10s24
                   i1p=i1+noc(isblk1(1,is))                             5d10s24
                   iad1=id3x(is)+i1+nvirt(isblk1(1,is))*(i2             5d10s24
     $                 +nvirt(isblk1(2,is))*(i3+irefo(isblk1(3,is))*i4))5d10s24
                   iad2=1+i2p+nbasdws(isa)*(i1p+nbasdws(isb)*(i3p        5d2s24
     $                 +nbasdws(isc)*i4p))                               5d2s24
                   e1x=e1x+bc(iad1)*tmp(iad2)                           5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(4)=iuseage(4)+1                                  5d6s24
               do i4=0,nvirt(isblk1(4,is))-1                             5d2s24
                i4p=i4+noc(isblk1(4,is))                                5d2s24
                do i3=0,irefo(isblk1(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk1(3,is))
                 do i2=0,nvirt(isblk1(2,is))-1                           5d2s24
                  i2p=i2+noc(isblk1(2,is))                              5d10s24
                  do i1=0,nvirt(isblk1(1,is))-1                         5d10s24
                   i1p=i1+noc(isblk1(1,is))                             5d10s24
                   iad1=id3x(is)+i1+nvirt(isblk1(1,is))*(i2             5d10s24
     $                 +nvirt(isblk1(2,is))*(i3+irefo(isblk1(3,is))*i4))5d10s24
                   iad2=1+i2p+nbasdws(isa)*(i1p+nbasdws(isb)*(i4p        5d2s24
     $                 +nbasdws(isc)*i3p))                               5d2s24
                   e1x=e1x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             end if                                                     5d2s24
            else                                                        5d2s24
             if(isblk1(3,is).eq.isa)then                                 5d2s24
              if(isblk1(1,is).eq.isc)then                                5d2s24
               iuseage(5)=iuseage(5)+1                                  5d6s24
               do i4=0,nvirt(isblk1(4,is))-1                             5d2s24
                i4p=i4+noc(isblk1(4,is))                                5d2s24
                do i3=0,irefo(isblk1(3,is))-1                           5d2s24
                 i3p=i3+idoubo(isblk1(3,is))
                 do i2=0,nvirt(isblk1(2,is))-1                          5d10s24
                  i2p=i2+noc(isblk1(2,is))                              5d10s24
                  i1top=i2+isw*(nvirt(isblk1(1,is))-1-i2)               5d10s24
                  do i1=0,i1top                                         5d2s24
                   i1p=i1+noc(isblk1(1,is))                             5d10s24
                   irec12=i1+nvirt(isblk1(1,is))*i2                     5d10s24
                   itri12=((i2*(i2+1))/2)+i1                            5d2s24
                   itri12=itri12+isw*(irec12-itri12)
                   iad1=id3x(is)+itri12+nrow*(i3+irefo(isblk1(3,is))*i4)5d2s24
                   iad2=1+i3p+nbasdws(isa)*(i4p+nbasdws(isb)*(i1p        5d2s24
     $                 +nbasdws(isc)*i2p))                               5d2s24
                   e1x=e1x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(6)=iuseage(6)+1                                  5d6s24
               do i4=0,nvirt(isblk1(4,is))-1                             5d2s24
                i4p=i4+noc(isblk1(4,is))                                5d2s24
                do i3=0,irefo(isblk1(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk1(3,is))
                 do i2=0,nvirt(isblk1(2,is))-1                          5d10s24
                  i2p=i2+noc(isblk1(2,is))                              5d10s24
                  do i1=0,nvirt(isblk1(1,is))-1                         5d10s24
                   i1p=i1+noc(isblk1(1,is))                             5d10s24
                   iad1=id3x(is)+i1+nvirt(isblk1(1,is))*(i2             5d10s24
     $                 +nvirt(isblk1(2,is))*(i3+irefo(isblk1(3,is))*i4))5d10s24
                   iad2=1+i3p+nbasdws(isa)*(i4p+nbasdws(isb)*(i2p        5d2s24
     $                 +nbasdws(isc)*i1p))                               5d2s24
                   e1x=e1x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             else                                                       5d2s24
              if(isblk1(1,is).eq.isc)then                                5d2s24
               iuseage(7)=iuseage(7)+1                                  5d6s24
               do i4=0,nvirt(isblk1(4,is))-1                            5d2s24
                i4p=i4+noc(isblk1(4,is))                                5d2s24
                do i3=0,irefo(isblk1(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk1(3,is))
                 do i2=0,nvirt(isblk1(2,is))-1                          5d10s24
                  i2p=i2+noc(isblk1(2,is))                              5d10s24
                  do i1=0,nvirt(isblk1(1,is))-1                         5d10s24
                   i1p=i1+noc(isblk1(1,is))                             5d10s24
                   iad1=id3x(is)+i1+nvirt(isblk1(1,is))*(i2             5d10s24
     $                 +nvirt(isblk1(2,is))*(i3+irefo(isblk1(3,is))*i4))5d10s24
                   iad2=1+i4p+nbasdws(isa)*(i3p+nbasdws(isb)*(i1p        5d2s24
     $                 +nbasdws(isc)*i2p))                               5d2s24
                   e1x=e1x+bc(iad1)*tmp(iad2)                           5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(8)=iuseage(8)+1                                  5d6s24
               do i4=0,nvirt(isblk1(4,is))-1                             5d2s24
                i4p=i4+noc(isblk1(4,is))                                5d2s24
                do i3=0,irefo(isblk1(3,is))-1                            5d2s24
                 i3p=i3+idoubo(isblk1(3,is))
                 do i2=0,nvirt(isblk1(2,is))-1                          5d10s24
                  i2p=i2+noc(isblk1(2,is))                              5d10s24
                  do i1=0,nvirt(isblk1(1,is))-1                         5d10s24
                   i1p=i1+noc(isblk1(1,is))                             5d10s24
                   iad1=id3x(is)+i1+nvirt(isblk1(1,is))*(i2             5d10s24
     $                 +nvirt(isblk1(2,is))*(i3+irefo(isblk1(3,is))*i4))5d10s24
                   iad2=1+i4p+nbasdws(isa)*(i3p+nbasdws(isb)*(i2p        5d2s24
     $                 +nbasdws(isc)*i1p))                               5d2s24
                   e1x=e1x+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             end if                                                     5d2s24
            end if                                                      5d2s24
           end if
          end do                                                        5d2s24
          return                                                        5d2s24
          end                                                           5d2s24
          subroutine ejjd(tmp,isa,isb,isc,isd,nbasdws,idoubo,irefo,     5d2s24
     $         nvirt,noc,isblk,nsdlk,ejj,idjj,bc,ibc,iuseage,itcode,    5d6s24
     $     xjt)                                                         5d6s24
          implicit real*8 (a-h,o-z)
          include "common.store"
          dimension tmp(*),nbasdws(*),idoubo(*),irefo(*),isblk(4,*),    5d2s24
     $         idjj(*),nvirt(*),noc(*),iuseage(*),xjt(*)                5d6s24
c     isb gt isa
c     isd gt isc
c     icd ge iab
          iab=((isb*(isb-1))/2)+isa                                     5d2s24
          icd=((isd*(isd-1))/2)+isc                                     5d2s24
          do is=1,nsdlk
           ix=max(isblk(1,is),isblk(2,is))                              5d2s24
           in=min(isblk(1,is),isblk(2,is))                              5d2s24
           i12=((ix*(ix-1))/2)+in                                       5d2s24
           ix=max(isblk(3,is),isblk(4,is))                              5d2s24
           in=min(isblk(3,is),isblk(4,is))                              5d2s24
           i34=((ix*(ix-1))/2)+in                                       5d2s24
           ix=max(i12,i34)
           in=min(i12,i34)
           if(ix.eq.icd.and.in.eq.iab)then
            if(isblk(1,is).eq.isblk(2,is))then
             nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2         5d2s24
             isw=0
            else
             nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                 5d2s24
             isw=1
            end if
            if(i12.eq.iab)then                                          5d2s24
             if(isblk(1,is).eq.isa)then                                 5d2s24
              if(isblk(3,is).eq.isc)then                                5d2s24
               iuseage(1)=iuseage(1)+1                                  5d6s24
               do i4=0,nvirt(isblk(4,is))-1                             5d2s24
                i4p=i4+noc(isblk(4,is))                                 5d2s24
                do i3=0,nvirt(isblk(3,is))-1                            5d2s24
                 i3p=i3+noc(isblk(3,is))
                 do i2=0,irefo(isblk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk(2,is))
                  i1top=i2+isw*(irefo(isblk(1,is))-1-i2)                5d2s24
                  do i1=0,i1top                                         5d2s24
                   i1p=i1+idoubo(isblk(1,is))
                   irec12=i1+irefo(isblk(1,is))*i2
                   itri12=((i2*(i2+1))/2)+i1                            5d2s24
                   itri12=itri12+isw*(irec12-itri12)
c 1234
                   iad1=idjj(is)+i3+nvirt(isblk(3,is))*(i4
     $                  +nvirt(isblk(4,is))*itri12)                     5d6s24
                   iad2=1+i1p+nbasdws(isa)*(i2p+nbasdws(isb)*(i3p        5d2s24
     $                 +nbasdws(isc)*i4p))                               5d2s24
                   ejj=ejj+bc(iad1)*tmp(iad2)                            5d2s24
                   if(itcode.ne.0)then
                    xjt(is)=xjt(is)+bc(iad1)*tmp(iad2)
                    if(is.eq.169)write(6,*)('xjt1 '),bc(iad1),tmp(iad2)
                   end if
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(2)=iuseage(2)+1                                  5d6s24
               do i4=0,nvirt(isblk(4,is))-1                             5d2s24
                i4p=i4+noc(isblk(4,is))                                 5d2s24
                do i3=0,nvirt(isblk(3,is))-1                            5d2s24
                 i3p=i3+noc(isblk(3,is))                                5d2s24
                 do i2=0,irefo(isblk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk(2,is))
                  do i1=0,irefo(isblk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk(1,is))
                   iad1=idjj(is)+i3+nvirt(isblk(3,is))*(i4              5d6s24
     $                  +nvirt(isblk(4,is))*(i1+irefo(isblk(1,is))*i2)) 5d6s24
c     1243
                   iad2=1+i1p+nbasdws(isa)*(i2p+nbasdws(isb)*(i4p        5d2s24
     $                 +nbasdws(isc)*i3p))                               5d2s24
                   ejj=ejj+bc(iad1)*tmp(iad2)                            5d2s24
                   if(itcode.ne.0)then
                    xjt(is)=xjt(is)+bc(iad1)*tmp(iad2)
                    if(is.eq.169)write(6,*)('xjt2 '),bc(iad1),tmp(iad2)
                   end if
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             else                                                       5d2s24
              if(isblk(3,is).eq.isc)then                                5d2s24
               iuseage(3)=iuseage(3)+1                                  5d6s24
               do i4=0,nvirt(isblk(4,is))-1                             5d2s24
                i4p=i4+noc(isblk(4,is))                                 5d2s24
                do i3=0,nvirt(isblk(3,is))-1                            5d2s24
                 i3p=i3+noc(isblk(3,is))
                 do i2=0,irefo(isblk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk(2,is))
                  do i1=0,irefo(isblk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk(1,is))
                   iad1=idjj(is)+i3+nvirt(isblk(3,is))*(i4+             5d6s24
     $                  nvirt(isblk(4,is))*(i1+irefo(isblk(1,is))*i2))  5d6s24
c     2134
                   iad2=1+i2p+nbasdws(isa)*(i1p+nbasdws(isb)*(i3p        5d2s24
     $                 +nbasdws(isc)*i4p))                               5d2s24
                   ejj=ejj+bc(iad1)*tmp(iad2)                           5d2s24
                   if(itcode.ne.0)then
                    xjt(is)=xjt(is)+bc(iad1)*tmp(iad2)
                    if(is.eq.169)write(6,*)('xjt3 '),bc(iad1),tmp(iad2),
     $                   iad2,(isblk(j,is),j=1,4),isa,isb,isc,isd,
     $                   i1,i1p,i2,i2p,i3,i3p,i4,i4p
                   end if
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(4)=iuseage(4)+1                                  5d6s24
               do i4=0,nvirt(isblk(4,is))-1                             5d2s24
                i4p=i4+noc(isblk(4,is))                                 5d2s24
                do i3=0,nvirt(isblk(3,is))-1                            5d2s24
                 i3p=i3+noc(isblk(3,is))
                 do i2=0,irefo(isblk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk(2,is))
                  do i1=0,irefo(isblk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk(1,is))
                   iad1=idjj(is)+i3+nvirt(isblk(3,is))*(i4              5d6s24
     $                  +nvirt(isblk(4,is))*(i1+irefo(isblk(1,is))*i2)) 5d6s24
c     2143
                   iad2=1+i2p+nbasdws(isa)*(i1p+nbasdws(isb)*(i4p        5d2s24
     $                 +nbasdws(isc)*i3p))                               5d2s24
                   ejj=ejj+bc(iad1)*tmp(iad2)                            5d2s24
                   if(itcode.ne.0)then
                    xjt(is)=xjt(is)+bc(iad1)*tmp(iad2)
                    if(is.eq.169)write(6,*)('xjt4 '),bc(iad1),tmp(iad2)
                   end if
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             end if                                                     5d2s24
            else                                                        5d2s24
             if(isblk(3,is).eq.isa)then                                 5d2s24
              if(isblk(1,is).eq.isc)then                                5d2s24
               iuseage(5)=iuseage(5)+1                                  5d6s24
               do i4=0,nvirt(isblk(4,is))-1                             5d2s24
                i4p=i4+noc(isblk(4,is))                                 5d2s24
                do i3=0,nvirt(isblk(3,is))-1                            5d2s24
                 i3p=i3+noc(isblk(3,is))
                 do i2=0,irefo(isblk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk(2,is))
                  i1top=i2+isw*(irefo(isblk(1,is))-1-i2)                5d2s24
                  do i1=0,i1top                                         5d2s24
                   i1p=i1+idoubo(isblk(1,is))
                   irec12=i1+irefo(isblk(1,is))*i2
                   itri12=((i2*(i2+1))/2)+i1                            5d2s24
                   itri12=itri12+isw*(irec12-itri12)
                   iad1=idjj(is)+i3+nvirt(isblk(3,is))*(i4              5d6s24
     $                  +nvirt(isblk(4,is))*itri12)                     5d6s24
c     3412
                   iad2=1+i3p+nbasdws(isa)*(i4p+nbasdws(isb)*(i1p        5d2s24
     $                 +nbasdws(isc)*i2p))                               5d2s24
                   ejj=ejj+bc(iad1)*tmp(iad2)                            5d2s24
                   if(itcode.ne.0)then
                    xjt(is)=xjt(is)+bc(iad1)*tmp(iad2)
                    if(is.eq.169)write(6,*)('xjt5 '),bc(iad1),tmp(iad2)
                   end if
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(6)=iuseage(6)+1                                  5d6s24
               do i4=0,nvirt(isblk(4,is))-1                             5d2s24
                i4p=i4+noc(isblk(4,is))                                 5d2s24
                do i3=0,nvirt(isblk(3,is))-1                            5d2s24
                 i3p=i3+noc(isblk(3,is))
                 do i2=0,irefo(isblk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk(2,is))
                  do i1=0,irefo(isblk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk(1,is))
                   iad1=idjj(is)+i3+nvirt(isblk(3,is))*(i4              5d6s24
     $                  +nvirt(isblk(4,is))*(i1+irefo(isblk(1,is))*i2)) 5d6s24
c     3421
                   iad2=1+i3p+nbasdws(isa)*(i4p+nbasdws(isb)*(i2p        5d2s24
     $                 +nbasdws(isc)*i1p))                               5d2s24
                   ejj=ejj+bc(iad1)*tmp(iad2)                            5d2s24
                   if(itcode.ne.0)then
                    xjt(is)=xjt(is)+bc(iad1)*tmp(iad2)
                    if(is.eq.169)write(6,*)('xjt6 '),bc(iad1),tmp(iad2)
                   end if
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             else                                                       5d2s24
              if(isblk(1,is).eq.isc)then                                5d2s24
               iuseage(7)=iuseage(7)+1                                  5d6s24
               do i4=0,nvirt(isblk(4,is))-1                             5d2s24
                i4p=i4+noc(isblk(4,is))                                 5d2s24
                do i3=0,nvirt(isblk(3,is))-1                            5d2s24
                 i3p=i3+noc(isblk(3,is))
                 do i2=0,irefo(isblk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk(2,is))
                  do i1=0,irefo(isblk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk(1,is))
                   iad1=idjj(is)+i3+nvirt(isblk(3,is))*(i4              5d6s24
     $                  +nvirt(isblk(4,is))*(i1+irefo(isblk(1,is))*i2)) 5d6s24
c     4312
                   iad2=1+i4p+nbasdws(isa)*(i3p+nbasdws(isb)*(i1p        5d2s24
     $                 +nbasdws(isc)*i2p))                               5d2s24
                   ejj=ejj+bc(iad1)*tmp(iad2)                           5d2s24
                   if(itcode.ne.0)then
                    xjt(is)=xjt(is)+bc(iad1)*tmp(iad2)
                    if(is.eq.169)write(6,*)('xjt7 '),bc(iad1),tmp(iad2)
                   end if
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(8)=iuseage(8)+1                                  5d6s24
               do i4=0,nvirt(isblk(4,is))-1                             5d2s24
                i4p=i4+noc(isblk(4,is))                                 5d2s24
                do i3=0,nvirt(isblk(3,is))-1                            5d2s24
                 i3p=i3+noc(isblk(3,is))
                 do i2=0,irefo(isblk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblk(2,is))
                  do i1=0,irefo(isblk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblk(1,is))
                   iad1=idjj(is)+i3+nvirt(isblk(3,is))*(i4              5d6s24
     $                  +nvirt(isblk(4,is))*(i1+irefo(isblk(1,is))*i2)) 5d6s24
c     4321
                   iad2=1+i4p+nbasdws(isa)*(i3p+nbasdws(isb)*(i2p        5d2s24
     $                 +nbasdws(isc)*i1p))                               5d2s24
                   ejj=ejj+bc(iad1)*tmp(iad2)                            5d2s24
                   if(itcode.ne.0)then
                    xjt(is)=xjt(is)+bc(iad1)*tmp(iad2)
                    if(is.eq.169)write(6,*)('xjt8 '),bc(iad1),tmp(iad2)
                   end if
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             end if                                                     5d2s24
            end if                                                      5d2s24
           end if
          end do                                                        5d2s24
          return                                                        5d2s24
          end                                                           5d2s24
          subroutine ekkd(tmp,isa,isb,isc,isd,nbasdws,idoubo,irefo,     5d2s24
     $         nvirt,noc,isblkk,nsdlkk,ekk,idkk,bc,ibc,iuseage)         5d6s24
          implicit real*8 (a-h,o-z)
          include "common.store"
          dimension tmp(*),nbasdws(*),idoubo(*),irefo(*),isblkk(4,*),    5d2s24
     $         idkk(*),nvirt(*),noc(*),iuseage(*)
c     isb gt isa
c     isd gt isc
c     icd ge iab
          iab=((isb*(isb-1))/2)+isa                                     5d2s24
          icd=((isd*(isd-1))/2)+isc                                     5d2s24
          do is=1,nsdlkk
c
c             i14=iab
c     Kab^uv=(av|ub) 1=a,4=b
c            (av|bu)
c            (va|ub)
c            (va|bu)
c            i24=iab
c            (ub|av)
c            (ub|va)
c            (bu|av)
c            (bu|va)
           ix=max(isblkk(1,is),isblkk(4,is))                              5d2s24
           in=min(isblkk(1,is),isblkk(4,is))                              5d2s24
           i14=((ix*(ix-1))/2)+in                                       5d2s24
           ix=max(isblkk(2,is),isblkk(3,is))                              5d2s24
           in=min(isblkk(2,is),isblkk(3,is))                              5d2s24
           i23=((ix*(ix-1))/2)+in                                       5d2s24
           ix=max(i14,i23)
           in=min(i14,i23)
           if(ix.eq.icd.and.in.eq.iab)then
            if(i14.eq.iab)then                                          5d2s24
             if(isblkk(1,is).eq.isa)then                                 5d2s24
              if(isblkk(2,is).eq.isd)then                                5d2s24
               iuseage(1)=iuseage(1)+1                                  5d6s24
               do i4=0,nvirt(isblkk(4,is))-1                             5d2s24
                i4p=i4+noc(isblkk(4,is))                                 5d2s24
                do i3=0,nvirt(isblkk(3,is))-1                            5d2s24
                 i3p=i3+noc(isblkk(3,is))
                 do i2=0,irefo(isblkk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblkk(2,is))
                  do i1=0,irefo(isblkk(1,is))-1
                   i1p=i1+idoubo(isblkk(1,is))
c 1324
                   iad1=idkk(is)+i3+nvirt(isblkk(3,is))*(i4             5d6s24
     $                 +nvirt(isblkk(4,is))*(i1+irefo(isblkk(1,is))*i2))5d6s24
                   iad2=1+i1p+nbasdws(isa)*(i4p+nbasdws(isb)*(i3p        5d2s24
     $                 +nbasdws(isc)*i2p))                               5d2s24
                   ekk=ekk+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(2)=iuseage(2)+1                                  5d6s24
               do i4=0,nvirt(isblkk(4,is))-1                             5d2s24
                i4p=i4+noc(isblkk(4,is))                                 5d2s24
                do i3=0,nvirt(isblkk(3,is))-1                            5d2s24
                 i3p=i3+noc(isblkk(3,is))                                5d2s24
                 do i2=0,irefo(isblkk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblkk(2,is))
                  do i1=0,irefo(isblkk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblkk(1,is))
                   iad1=idkk(is)+i3+nvirt(isblkk(3,is))*(i4             5d6s24
     $                 +nvirt(isblkk(4,is))*(i1+irefo(isblkk(1,is))*i2))5d6s24
c     1423
                   iad2=1+i1p+nbasdws(isa)*(i4p+nbasdws(isb)*(i2p        5d2s24
     $                 +nbasdws(isc)*i3p))                               5d2s24
                   ekk=ekk+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             else                                                       5d2s24
              if(isblkk(2,is).eq.isd)then                                5d2s24
               iuseage(3)=iuseage(3)+1                                  5d6s24
               do i4=0,nvirt(isblkk(4,is))-1                             5d2s24
                i4p=i4+noc(isblkk(4,is))                                 5d2s24
                do i3=0,nvirt(isblkk(3,is))-1                            5d2s24
                 i3p=i3+noc(isblkk(3,is))
                 do i2=0,irefo(isblkk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblkk(2,is))
                  do i1=0,irefo(isblkk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblkk(1,is))
                   iad1=idkk(is)+i3+nvirt(isblkk(3,is))*(i4             5d6s24
     $                 +nvirt(isblkk(4,is))*(i1+irefo(isblkk(1,is))*i2))5d6s24
c     4132
                   iad2=1+i4p+nbasdws(isa)*(i1p+nbasdws(isb)*(i3p        5d2s24
     $                 +nbasdws(isc)*i2p))                               5d2s24
                   ekk=ekk+bc(iad1)*tmp(iad2)                           5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(4)=iuseage(4)+1                                  5d6s24
               do i4=0,nvirt(isblkk(4,is))-1                             5d2s24
                i4p=i4+noc(isblkk(4,is))                                 5d2s24
                do i3=0,nvirt(isblkk(3,is))-1                            5d2s24
                 i3p=i3+noc(isblkk(3,is))
                 do i2=0,irefo(isblkk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblkk(2,is))
                  do i1=0,irefo(isblkk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblkk(1,is))
                   iad1=idkk(is)+i3+nvirt(isblkk(3,is))*(i4             5d6s24
     $                 +nvirt(isblkk(4,is))*(i1+irefo(isblkk(1,is))*i2))5d6s24
c     4123
                   iad2=1+i4p+nbasdws(isa)*(i1p+nbasdws(isb)*(i2p        5d2s24
     $                 +nbasdws(isc)*i3p))                               5d2s24
                   ekk=ekk+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             end if                                                     5d2s24
            else                                                        5d2s24
             if(isblkk(2,is).eq.isb)then                                 5d2s24
              if(isblkk(1,is).eq.isc)then                                5d2s24
               iuseage(5)=iuseage(5)+1                                  5d6s24
               do i4=0,nvirt(isblkk(4,is))-1                             5d2s24
                i4p=i4+noc(isblkk(4,is))                                 5d2s24
                do i3=0,nvirt(isblkk(3,is))-1                            5d2s24
                 i3p=i3+noc(isblkk(3,is))
                 do i2=0,irefo(isblkk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblkk(2,is))
                  do i1=0,irefo(isblkk(1,is))-1
                   i1p=i1+idoubo(isblkk(1,is))
                   irec12=i1+irefo(isblkk(1,is))*i2
                   iad1=idkk(is)+i3+nvirt(isblkk(3,is))*(i4             5d6s24
     $                 +nvirt(isblkk(4,is))*(i1+irefo(isblkk(1,is))*i2))5d6s24
c     3214
                   iad2=1+i3p+nbasdws(isa)*(i2p+nbasdws(isb)*(i1p        5d2s24
     $                 +nbasdws(isc)*i4p))                               5d2s24
                   ekk=ekk+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(6)=iuseage(6)+1                                  5d6s24
               do i4=0,nvirt(isblkk(4,is))-1                             5d2s24
                i4p=i4+noc(isblkk(4,is))                                 5d2s24
                do i3=0,nvirt(isblkk(3,is))-1                            5d2s24
                 i3p=i3+noc(isblkk(3,is))
                 do i2=0,irefo(isblkk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblkk(2,is))
                  do i1=0,irefo(isblkk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblkk(1,is))
                   iad1=idkk(is)+i3+nvirt(isblkk(3,is))*(i4             5d6s24
     $                 +nvirt(isblkk(4,is))*(i1+irefo(isblkk(1,is))*i2))5d6s24
c     3241
                   iad2=1+i3p+nbasdws(isa)*(i2p+nbasdws(isb)*(i4p        5d2s24
     $                 +nbasdws(isc)*i1p))                               5d2s24
                   ekk=ekk+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             else                                                       5d2s24
              if(isblkk(1,is).eq.isc)then                                5d2s24
               iuseage(7)=iuseage(7)+1                                  5d6s24
               do i4=0,nvirt(isblkk(4,is))-1                             5d2s24
                i4p=i4+noc(isblkk(4,is))                                 5d2s24
                do i3=0,nvirt(isblkk(3,is))-1                            5d2s24
                 i3p=i3+noc(isblkk(3,is))
                 do i2=0,irefo(isblkk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblkk(2,is))
                  do i1=0,irefo(isblkk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblkk(1,is))
                   iad1=idkk(is)+i3+nvirt(isblkk(3,is))*(i4             5d6s24
     $                 +nvirt(isblkk(4,is))*(i1+irefo(isblkk(1,is))*i2))5d6s24
c     2314
                   iad2=1+i2p+nbasdws(isa)*(i3p+nbasdws(isb)*(i1p        5d2s24
     $                 +nbasdws(isc)*i4p))                               5d2s24
                   ekk=ekk+bc(iad1)*tmp(iad2)                           5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              else                                                      5d2s24
               iuseage(8)=iuseage(8)+1                                  5d6s24
               do i4=0,nvirt(isblkk(4,is))-1                             5d2s24
                i4p=i4+noc(isblkk(4,is))                                 5d2s24
                do i3=0,nvirt(isblkk(3,is))-1                            5d2s24
                 i3p=i3+noc(isblkk(3,is))
                 do i2=0,irefo(isblkk(2,is))-1                           5d2s24
                  i2p=i2+idoubo(isblkk(2,is))
                  do i1=0,irefo(isblkk(1,is))-1                          5d2s24
                   i1p=i1+idoubo(isblkk(1,is))
                   iad1=idkk(is)+i3+nvirt(isblkk(3,is))*(i4             5d6s24
     $                 +nvirt(isblkk(4,is))*(i1+irefo(isblkk(1,is))*i2))5d6s24
c     2341
                   iad2=1+i2p+nbasdws(isa)*(i3p+nbasdws(isb)*(i4p        5d2s24
     $                 +nbasdws(isc)*i1p))                               5d2s24
                   ekk=ekk+bc(iad1)*tmp(iad2)                            5d2s24
                  end do                                                5d2s24
                 end do                                                 5d2s24
                end do                                                  5d2s24
               end do                                                   5d2s24
              end if
             end if                                                     5d2s24
            end if                                                      5d2s24
           end if
          end do                                                        5d2s24
          return                                                        5d2s24
          end                                                           5d2s24
          subroutine e4vd(tmp,lsa,lsb,lsc,lsd,nbasdws,nvirt,noc,e4v,    5d10s24
     $     vd,nfdat,multh,nsymb,isymmrci,bc,ibc,iuseage,srh)            5d10s24
          implicit real*8 (a-h,o-z)                                     5d10s24
          include "common.store"                                        5d10s24
          dimension tmp(*),nbasdws(*),nvirt(*),noc(*),vd(*),            5d10s24
     $         nfdat(5,4,*),multh(8,8),iuseage(*),iperm(4),iiii(4)      5d10s24
c     lsb gt lsa
c     lsd gt lsc
c     lcd ge lab
          lab=((lsb*(lsb-1))/2)+lsa                                     5d2s24
          lcd=((lsd*(lsd-1))/2)+lsc                                     5d2s24
          ioff=1                                                        5d10s24
          do isb=1,nsymb                                                5d10s24
           isbv12=multh(isb,isymmrci)                                   5d10s24
           do isbv1=1,nsymb                                             5d10s24
            isbv2=multh(isbv1,isbv12)                                   5d10s24
            if(isbv2.ge.isbv1)then                                      5d10s24
             if(isbv2.eq.isbv1)then                                     5d10s24
              nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                     5d10s24
              isw=0                                                     5d10s24
              ioffvv=ioff                                               5d10s24
              ioff=ioff+nvirt(isbv1)*nfdat(3,1,isb)                     5d10s24
             else                                                       5d10s24
              nvv=nvirt(isbv1)*nvirt(isbv2)                             5d10s24
              isw=1                                                     5d10s24
             end if                                                     5d10s24
             joff=1                                                     5d10s24
             do jsb=1,nsymb                                             5d10s24
              jsbv12=multh(jsb,isymmrci)                                5d10s24
              do jsbv1=1,nsymb                                          5d10s24
               jsbv2=multh(jsbv1,jsbv12)                                5d10s24
               if(jsbv2.ge.jsbv1)then                                   5d10s24
                if(jsbv2.eq.jsbv1)then                                  5d10s24
                 mvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                  5d10s24
                 jsw=0                                                  5d10s24
                 joffvv=joff                                            5d10s24
                 joff=joff+nvirt(jsbv1)*nfdat(3,1,jsb)                  5d10s24
                else                                                    5d10s24
                 mvv=nvirt(jsbv1)*nvirt(jsbv2)                          5d10s24
                 jsw=1                                                  5d10s24
                end if                                                  5d10s24
                if(isb.eq.jsb)then                                      5d10s24
c     (isbv1 jsbv1|isbv2 jsbv2)
c     tfact*(isbv1 jsbv2|isbv2 jsbv1)
                 ix=max(isbv1,jsbv1)                                    5d10s24
                 in=min(isbv1,jsbv1)                                    5d10s24
                 i12=((ix*(ix-1))/2)+in                                       5d2s24
                 ix=max(isbv2,jsbv2)                                    5d10s24
                 in=min(isbv2,jsbv2)                                    5d10s24
                 i34=((ix*(ix-1))/2)+in                                       5d2s24
                 ix=max(isbv1,jsbv2)                                    5d10s24
                 in=min(isbv1,jsbv2)                                    5d10s24
                 i12r=((ix*(ix-1))/2)+in                                       5d2s24
                 ix=max(isbv2,jsbv1)                                    5d10s24
                 in=min(isbv2,jsbv1)                                    5d10s24
                 i34r=((ix*(ix-1))/2)+in                                       5d2s24
                 ix=max(i12,i34)
                 in=min(i12,i34)
                 if(ix.eq.lcd.and.in.eq.lab)then                        5d10s24
                  if(i12.eq.lab)then                                    5d10s24
                   if(isbv1.eq.lsa)then                                 5d10s24
                    if(isbv2.eq.lsc)then
                     iperm(1)=1                                         5d10s24
                     iperm(2)=2                                         5d10s24
                     iperm(3)=3                                         5d10s24
                     iperm(4)=4                                         5d10s24
                     iuseage(1)=iuseage(1)+1                            5d10s24
                    else                                                5d10s24
                     iperm(1)=1                                         5d10s24
                     iperm(2)=2                                         5d10s24
                     iperm(3)=4                                         5d10s24
                     iperm(4)=3                                         5d10s24
                     iuseage(2)=iuseage(2)+1                            5d10s24
                    end if                                              5d10s24
                   else                                                 5d10s24
                    if(isbv2.eq.lsc)then                                5d10s24
                     iperm(1)=2                                         5d10s24
                     iperm(2)=1                                         5d10s24
                     iperm(3)=3                                         5d10s24
                     iperm(4)=4                                         5d10s24
                     iuseage(3)=iuseage(3)+1                            5d10s24
                    else                                                5d10s24
                     iperm(1)=2                                         5d10s24
                     iperm(2)=1                                         5d10s24
                     iperm(3)=4                                         5d10s24
                     iperm(4)=3                                         5d10s24
                     iuseage(4)=iuseage(4)+1                            5d10s24
                    end if                                              5d10s24
                   end if                                               5d10s24
                  else                                                  5d10s24
                   if(isbv2.eq.lsa)then                                 5d10s24
                    if(isbv1.eq.lsc)then                                5d10s24
                     iperm(1)=3                                         5d10s24
                     iperm(2)=4                                         5d10s24
                     iperm(3)=1                                         5d10s24
                     iperm(4)=2                                         5d10s24
                     iuseage(5)=iuseage(5)+1                            5d10s24
                    else                                                5d10s24
                     iperm(1)=3                                         5d10s24
                     iperm(2)=4                                         5d10s24
                     iperm(3)=2                                         5d10s24
                     iperm(4)=1                                         5d10s24
                     iuseage(6)=iuseage(6)+1                            5d10s24
                    end if                                              5d10s24
                   else                                                 5d10s24
                    if(isbv1.eq.lsc)then                                5d10s24
                     iperm(1)=4                                         5d10s24
                     iperm(2)=3                                         5d10s24
                     iperm(3)=1                                         5d10s24
                     iperm(4)=2                                         5d10s24
                     iuseage(7)=iuseage(7)+1                            5d10s24
                    else                                                5d10s24
                     iperm(1)=4                                         5d10s24
                     iperm(2)=3                                         5d10s24
                     iperm(3)=2                                         5d10s24
                     iperm(4)=1                                         5d10s24
                     iuseage(8)=iuseage(8)+1                            5d10s24
                    end if                                              5d10s24
                   end if                                               5d10s24
                  end if                                                5d10s24
                  if(isbv12.eq.1.and.jsbv12.eq.1)then                    5d10s24
                   do iv=0,nvirt(isbv1)-1                               5d10s24
                    iiii(1)=iv+noc(isbv1)                               5d10s24
                    iiii(3)=iv+noc(isbv1)                               5d10s24
                    do jv=0,nvirt(jsbv1)-1                              5d10s24
                     iiii(2)=jv+noc(jsbv1)                              5d10s24
                     iiii(4)=jv+noc(jsbv2)                              5d10s24
                     iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                     sum=0d0                                            5d10s24
                     iiad=ioffvv+iv
                     jjad=joffvv+jv
                     do i=0,nfdat(3,1,isb)-1                            5d10s24
                      sum=sum+vd(iiad)*vd(jjad)                         5d10s24
                      iiad=iiad+nvirt(isbv1)                            5d10s24
                      jjad=jjad+nvirt(jsbv1)                            5d10s24
                     end do                                             5d10s24
                     e4v=e4v+tmp(iad)*sum*0.5d0                         5d10s24
                    end do                                              5d10s24
                   end do                                               5d10s24
                  end if                                                5d10s24
                  iioff=ioff                                            5d10s24
                  jjoff=joff                                            5d10s24
                  if(jsbv12.eq.1)then                                   5d10s24
                   do jv=0,nvirt(jsbv1)-1                               5d10s24
                    iiii(4)=jv+noc(jsbv1)                               5d10s24
                    iiii(2)=iiii(4)                                     5d10s24
                    do iv2=0,nvirt(isbv2)-1                             5d10s24
                     iiii(3)=iv2+noc(isbv2)                             5d10s24
                     iv1top=iv2-1+isw*(nvirt(isbv1)-iv2)                5d10s24
                     do iv1=0,iv1top                                    5d10s24
                      iiii(1)=iv1+noc(isbv1)                            5d10s24
                      itri=((iv2*(iv2-1))/2)+iv1                        5d10s24
                      irec=iv1+nvirt(isbv1)*iv2                         5d10s24
                      itri=itri+isw*(irec-itri)                         5d10s24
                      iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                      sum=0d0                                            5d10s24
                      iiad=ioff+itri                                    5d10s24
                      jjad=joffvv+jv                                    5d10s24
                      do i=0,nfdat(3,1,isb)-1                            5d10s24
                       sum=sum+vd(iiad)*vd(jjad)                         5d10s24
                       iiad=iiad+nvv                                    5d10s24
                       jjad=jjad+nvirt(jsbv1)                           5d10s24
                      end do                                            5d10s24
                      e4v=e4v+sum*srh*tmp(iad)                          5d10s24
                     end do                                             5d10s24
                    end do                                              5d10s24
                   end do                                               5d10s24
                  end if                                                5d10s24
                  if(isbv12.eq.1)then                                   5d10s24
                   do jv2=0,nvirt(jsbv2)-1                              5d10s24
                    jv1top=jv2-1+jsw*(nvirt(jsbv1)-jv2)                 5d10s24
                    iiii(4)=jv2+noc(jsbv2)                              5d10s24
                    do jv1=0,jv1top                                     5d10s24
                     iiii(2)=jv1+noc(jsbv1)                             5d10s24
                     jtri=((jv2*(jv2-1))/2)+jv1                         5d10s24
                     jrec=jv1+nvirt(jsbv1)*jv2                          5d10s24
                     jtri=jtri+jsw*(jrec-jtri)                          5d10s24
                     do iv=0,nvirt(isbv1)-1                             5d10s24
                      iiii(1)=iv+noc(isbv1)                             5d10s24
                      iiii(3)=iiii(1)                                   5d10s24
                      iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                      sum=0d0                                            5d10s24
                      iiad=ioffvv+iv                                    5d10s24
                      jjad=jjoff+jtri                                   5d10s24
                      do i=0,nfdat(3,1,isb)-1                            5d10s24
                       sum=sum+vd(iiad)*vd(jjad)                         5d10s24
                       iiad=iiad+nvirt(isbv1)                            5d10s24
                       jjad=jjad+mvv                                    5d10s24
                      end do                                            5d10s24
                      e4v=e4v+sum*srh*tmp(iad)                          5d10s24
                     end do                                             5d10s24
                    end do                                              5d10s24
                   end do                                               5d10s24
                  end if                                                5d10s24
                  do l=1,4                                              5d10s24
                   do jv2=0,nvirt(jsbv2)-1                              5d10s24
                    jv1top=jv2-1+jsw*(nvirt(jsbv1)-jv2)                 5d10s24
                    iiii(4)=jv2+noc(jsbv2)                              5d10s24
                    do jv1=0,jv1top                                     5d10s24
                     iiii(2)=jv1+noc(jsbv1)                             5d10s24
                     jtri=((jv2*(jv2-1))/2)+jv1                         5d10s24
                     jrec=jv1+nvirt(jsbv1)*jv2                          5d10s24
                     jtri=jtri+jsw*(jrec-jtri)                          5d10s24
                     do iv2=0,nvirt(isbv2)-1                            5d10s24
                      iiii(3)=iv2+noc(isbv2)                            5d10s24
                      iv1top=iv2-1+isw*(nvirt(isbv1)-iv2)               5d10s24
                      do iv1=0,iv1top                                   5d10s24
                       iiii(1)=iv1+noc(isbv1)                           5d10s24
                       itri=((iv2*(iv2-1))/2)+iv1                       5d10s24
                       irec=iv1+nvirt(isbv1)*iv2                        5d10s24
                       itri=itri+isw*(irec-itri)                        5d10s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       sum=0d0                                            5d10s24
                       iiad=iioff+itri                                  5d10s24
                       jjad=jjoff+jtri                                  5d10s24
                       do i=0,nfdat(3,l,isb)-1                          5d10s24
                        sum=sum+vd(iiad)*vd(jjad)                         5d10s24
                        iiad=iiad+nvv                                   5d10s24
                        jjad=jjad+mvv                                   5d10s24
                       end do                                             5d10s24
                       e4v=e4v+tmp(iad)*sum                             5d10s24
                      end do
                     end do                                             5d10s24
                    end do                                              5d10s24
                   end do                                               5d10s24
                   iioff=iioff+nvv*nfdat(3,l,isb)                       5d10s24
                   jjoff=jjoff+mvv*nfdat(3,l,isb)                       5d10s24
                  end do                                                5d10s24
                 end if                                                 5d10s24
                 ix=max(i12r,i34r)
                 in=min(i12r,i34r)
                 if(ix.eq.lcd.and.in.eq.lab)then                        5d10s24
                  if(i12r.eq.lab)then                                   5d13s24
                   if(isbv1.eq.lsa)then                                 5d10s24
                    if(isbv2.eq.lsc)then
                     iperm(1)=1                                         5d10s24
                     iperm(2)=2                                         5d10s24
                     iperm(3)=3                                         5d10s24
                     iperm(4)=4                                         5d10s24
                     iuseage(1)=iuseage(1)+1                            5d10s24
                    else                                                5d10s24
                     iperm(1)=1                                         5d10s24
                     iperm(2)=2                                         5d10s24
                     iperm(3)=4                                         5d10s24
                     iperm(4)=3                                         5d10s24
                     iuseage(2)=iuseage(2)+1                            5d10s24
                    end if                                              5d10s24
                   else                                                 5d10s24
                    if(isbv2.eq.lsc)then                                5d10s24
                     iperm(1)=2                                         5d10s24
                     iperm(2)=1                                         5d10s24
                     iperm(3)=3                                         5d10s24
                     iperm(4)=4                                         5d10s24
                     iuseage(3)=iuseage(3)+1                            5d10s24
                    else                                                5d10s24
                     iperm(1)=2                                         5d10s24
                     iperm(2)=1                                         5d10s24
                     iperm(3)=4                                         5d10s24
                     iperm(4)=3                                         5d10s24
                     iuseage(4)=iuseage(4)+1                            5d10s24
                    end if                                              5d10s24
                   end if                                               5d10s24
                  else                                                  5d10s24
                   if(isbv2.eq.lsa)then                                 5d10s24
                    if(isbv1.eq.lsc)then                                5d10s24
                     iperm(1)=3                                         5d10s24
                     iperm(2)=4                                         5d10s24
                     iperm(3)=1                                         5d10s24
                     iperm(4)=2                                         5d10s24
                     iuseage(5)=iuseage(5)+1                            5d10s24
                    else                                                5d10s24
                     iperm(1)=3                                         5d10s24
                     iperm(2)=4                                         5d10s24
                     iperm(3)=2                                         5d10s24
                     iperm(4)=1                                         5d10s24
                     iuseage(6)=iuseage(6)+1                            5d10s24
                    end if                                              5d10s24
                   else                                                 5d10s24
                    if(isbv1.eq.lsc)then                                5d10s24
                     iperm(1)=4                                         5d10s24
                     iperm(2)=3                                         5d10s24
                     iperm(3)=1                                         5d10s24
                     iperm(4)=2                                         5d10s24
                     iuseage(7)=iuseage(7)+1                            5d10s24
                    else                                                5d10s24
                     iperm(1)=4                                         5d10s24
                     iperm(2)=3                                         5d10s24
                     iperm(3)=2                                         5d10s24
                     iperm(4)=1                                         5d10s24
                     iuseage(8)=iuseage(8)+1                            5d10s24
                    end if                                              5d10s24
                   end if                                               5d10s24
                  end if                                                5d10s24
                  if(isbv12.eq.1.and.jsbv12.eq.1)then                    5d10s24
                   do iv=0,nvirt(isbv1)-1                               5d10s24
                    iiii(1)=iv+noc(isbv1)                               5d10s24
                    iiii(3)=iv+noc(isbv1)                               5d10s24
                    do jv=0,nvirt(jsbv1)-1                              5d10s24
                     iiii(2)=jv+noc(jsbv1)                              5d10s24
                     iiii(4)=jv+noc(jsbv2)                              5d10s24
                     iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                     sum=0d0                                            5d10s24
                     iiad=ioffvv+iv
                     jjad=joffvv+jv
                     do i=0,nfdat(3,1,isb)-1                            5d10s24
                      sum=sum+vd(iiad)*vd(jjad)                         5d10s24
                      iiad=iiad+nvirt(isbv1)                            5d10s24
                      jjad=jjad+nvirt(jsbv1)                            5d10s24
                     end do                                             5d10s24
                     e4v=e4v+tmp(iad)*sum*0.5d0                         5d10s24
                    end do                                              5d10s24
                   end do                                               5d10s24
                  end if                                                5d10s24
                  iioff=ioff                                            5d10s24
                  jjoff=joff                                            5d10s24
                  if(jsbv12.eq.1)then                                   5d10s24
                   do jv=0,nvirt(jsbv1)-1                               5d10s24
                    iiii(4)=jv+noc(jsbv1)                               5d10s24
                    iiii(2)=iiii(4)                                     5d10s24
                    do iv2=0,nvirt(isbv2)-1                             5d10s24
                     iiii(3)=iv2+noc(isbv2)                             5d10s24
                     iv1top=iv2-1+isw*(nvirt(isbv1)-iv2)                5d10s24
                     do iv1=0,iv1top                                    5d10s24
                      iiii(1)=iv1+noc(isbv1)                            5d10s24
                      itri=((iv2*(iv2-1))/2)+iv1                        5d10s24
                      irec=iv1+nvirt(isbv1)*iv2                         5d10s24
                      itri=itri+isw*(irec-itri)                         5d10s24
                      iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                      sum=0d0                                            5d10s24
                      iiad=ioff+itri                                    5d10s24
                      jjad=joffvv+jv                                    5d10s24
                      do i=0,nfdat(3,1,isb)-1                            5d10s24
                       sum=sum+vd(iiad)*vd(jjad)                         5d10s24
                       iiad=iiad+nvv                                    5d10s24
                       jjad=jjad+nvirt(jsbv1)                           5d10s24
                      end do                                            5d10s24
                      e4v=e4v+sum*srh*tmp(iad)                          5d10s24
                     end do                                             5d10s24
                    end do                                              5d10s24
                   end do                                               5d10s24
                  end if                                                5d10s24
                  if(isbv12.eq.1)then                                   5d10s24
                   do jv2=0,nvirt(jsbv2)-1                              5d10s24
                    jv1top=jv2-1+jsw*(nvirt(jsbv1)-jv2)                 5d10s24
                    iiii(2)=jv2+noc(jsbv2)                              5d10s24
                    do jv1=0,jv1top                                     5d10s24
                     iiii(4)=jv1+noc(jsbv1)                             5d10s24
                     jtri=((jv2*(jv2-1))/2)+jv1                         5d10s24
                     jrec=jv1+nvirt(jsbv1)*jv2                          5d10s24
                     jtri=jtri+jsw*(jrec-jtri)                          5d10s24
                     do iv=0,nvirt(isbv1)-1                             5d10s24
                      iiii(1)=iv+noc(isbv1)                             5d10s24
                      iiii(3)=iiii(1)                                   5d10s24
                      iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                      sum=0d0                                            5d10s24
                      iiad=ioffvv+iv                                    5d10s24
                      jjad=jjoff+jtri                                   5d10s24
                      do i=0,nfdat(3,1,isb)-1                            5d10s24
                       sum=sum+vd(iiad)*vd(jjad)                         5d10s24
                       iiad=iiad+nvirt(isbv1)                            5d10s24
                       jjad=jjad+mvv                                    5d10s24
                      end do                                            5d10s24
                      e4v=e4v+sum*srh*tmp(iad)                          5d10s24
                     end do                                             5d10s24
                    end do                                              5d10s24
                   end do                                               5d10s24
                  end if                                                5d10s24
                  tf=1d0                                                5d10s24
                  do l=1,4                                              5d10s24
                   do jv2=0,nvirt(jsbv2)-1                              5d10s24
                    jv1top=jv2-1+jsw*(nvirt(jsbv1)-jv2)                 5d10s24
                    iiii(2)=jv2+noc(jsbv2)                              5d10s24
                    do jv1=0,jv1top                                     5d10s24
                     iiii(4)=jv1+noc(jsbv1)                             5d10s24
                     jtri=((jv2*(jv2-1))/2)+jv1                         5d10s24
                     jrec=jv1+nvirt(jsbv1)*jv2                          5d10s24
                     jtri=jtri+jsw*(jrec-jtri)                          5d10s24
                     do iv2=0,nvirt(isbv2)-1                            5d10s24
                      iiii(3)=iv2+noc(isbv2)                            5d10s24
                      iv1top=iv2-1+isw*(nvirt(isbv1)-iv2)               5d10s24
                      do iv1=0,iv1top                                   5d10s24
                       iiii(1)=iv1+noc(isbv1)                           5d10s24
                       itri=((iv2*(iv2-1))/2)+iv1                       5d10s24
                       irec=iv1+nvirt(isbv1)*iv2                        5d10s24
                       itri=itri+isw*(irec-itri)                        5d10s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       sum=0d0                                            5d10s24
                       iiad=iioff+itri                                  5d10s24
                       jjad=jjoff+jtri                                  5d10s24
                       do i=0,nfdat(3,l,isb)-1                          5d10s24
                        sum=sum+vd(iiad)*vd(jjad)                         5d10s24
                        iiad=iiad+nvv                                   5d10s24
                        jjad=jjad+mvv                                   5d10s24
                       end do                                             5d10s24
                       e4v=e4v+tmp(iad)*sum*tf                           5d10s24
                      end do
                     end do                                             5d10s24
                    end do                                              5d10s24
                   end do                                               5d10s24
                   iioff=iioff+nvv*nfdat(3,l,isb)                       5d10s24
                   jjoff=jjoff+mvv*nfdat(3,l,isb)                       5d10s24
                   tf=-1d0                                              5d10s24
                  end do                                                5d10s24
                 end if                                                 5d10s24
                end if
                do l=1,4                                                5d10s24
                 joff=joff+nfdat(3,l,jsb)*mvv                           5d10s24
                end do                                                  5d10s24
               end if                                                   5d10s24
              end do                                                    5d10s24
             end do                                                     5d10s24
             do l=1,4                                                   5d10s24
              ioff=ioff+nvv*nfdat(3,l,isb)                              5d10s24
             end do                                                     5d10s24
            end if                                                      5d10s24
           end do                                                       5d10s24
          end do                                                        5d10s24
          return                                                        5d2s24
          end                                                           5d2s24
          subroutine e4tvd(tmp,lsa,lsb,lsc,lsd,nbasdws,nvirt,noc,ie4v,    5d10s24
     $     vd,nfdat,multh,nsymb,isymmrci,bc,ibc,iuseage,srh,itmp)       6d2s24
          implicit real*8 (a-h,o-z)                                     5d10s24
          include "common.store"                                        5d10s24
          logical ldoj,ldok
          dimension tmp(*),nbasdws(*),nvirt(*),noc(*),vd(*),            5d10s24
     $         nfdat(5,4,*),multh(8,8),iuseage(*),iperm(4),iiii(4),     5d15s24
     $         ie4v(8,*)                                                  5d15s24
c     lsb gt lsa
c     lsd gt lsc
c     lcd ge lab
          ldoj=.true.
          ldok=.true.
          igoal=2
          lab=((lsb*(lsb-1))/2)+lsa                                     5d2s24
          lcd=((lsd*(lsd-1))/2)+lsc                                     5d2s24
          ioff=1                                                        5d10s24
          do isb=1,nsymb                                                5d10s24
           isbv12=multh(isb,isymmrci)                                   5d10s24
           do isbv1=1,nsymb                                             5d10s24
            isbv2=multh(isbv1,isbv12)                                   5d10s24
            if(isbv2.ge.isbv1)then                                      5d10s24
             if(isbv2.eq.isbv1)then                                     5d10s24
              nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                     5d10s24
              isw=0                                                     5d10s24
              ioffvv=ioff                                               5d10s24
              ioff=ioff+nvirt(isbv1)*nfdat(3,1,isb)                     5d10s24
             else                                                       5d10s24
              nvv=nvirt(isbv1)*nvirt(isbv2)                             5d10s24
              isw=1                                                     5d10s24
             end if                                                     5d10s24
             joff=1                                                     5d10s24
             do jsb=1,nsymb                                             5d10s24
              jsbv12=multh(jsb,isymmrci)                                5d10s24
              do jsbv1=1,nsymb                                          5d10s24
               jsbv2=multh(jsbv1,jsbv12)                                5d10s24
               if(jsbv2.ge.jsbv1)then                                   5d10s24
                if(jsbv2.eq.jsbv1)then                                  5d10s24
                 mvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                  5d10s24
                 jsw=0                                                  5d10s24
                 joffvv=joff                                            5d10s24
                 joff=joff+nvirt(jsbv1)*nfdat(3,1,jsb)                  5d10s24
                else                                                    5d10s24
                 mvv=nvirt(jsbv1)*nvirt(jsbv2)                          5d10s24
                 jsw=1                                                  5d10s24
                end if                                                  5d10s24
                if(isb.eq.jsb)then                                      5d10s24
c     (isbv1 jsbv1|isbv2 jsbv2)
c     tfact*(isbv1 jsbv2|isbv2 jsbv1)
                 ix=max(isbv1,jsbv1)                                    5d10s24
                 in=min(isbv1,jsbv1)                                    5d10s24
                 i12=((ix*(ix-1))/2)+in                                       5d2s24
                 ix=max(isbv2,jsbv2)                                    5d10s24
                 in=min(isbv2,jsbv2)                                    5d10s24
                 i34=((ix*(ix-1))/2)+in                                       5d2s24
                 ix=max(isbv1,jsbv2)                                    5d10s24
                 in=min(isbv1,jsbv2)                                    5d10s24
                 i12r=((ix*(ix-1))/2)+in                                       5d2s24
                 ix=max(isbv2,jsbv1)                                    5d10s24
                 in=min(isbv2,jsbv1)                                    5d10s24
                 i34r=((ix*(ix-1))/2)+in                                       5d2s24
                 ix=max(i12,i34)
                 in=min(i12,i34)
                 if(ldoj)then                                           6d11s24
                 if(ix.eq.lcd.and.in.eq.lab)then                        5d10s24
                  if(i12.eq.lab)then                                    5d10s24
                   if(isbv1.eq.lsa)then                                 5d10s24
                    if(isbv2.eq.lsc)then
                     iperm(1)=1                                         5d10s24
                     iperm(2)=2                                         5d10s24
                     iperm(3)=3                                         5d10s24
                     iperm(4)=4                                         5d10s24
                     iuseage(1)=iuseage(1)+1                            5d10s24
                    else                                                5d10s24
                     iperm(1)=1                                         5d10s24
                     iperm(2)=2                                         5d10s24
                     iperm(3)=4                                         5d10s24
                     iperm(4)=3                                         5d10s24
                     iuseage(2)=iuseage(2)+1                            5d10s24
                    end if                                              5d10s24
                   else                                                 5d10s24
                    if(isbv2.eq.lsc)then                                5d10s24
                     iperm(1)=2                                         5d10s24
                     iperm(2)=1                                         5d10s24
                     iperm(3)=3                                         5d10s24
                     iperm(4)=4                                         5d10s24
                     iuseage(3)=iuseage(3)+1                            5d10s24
                    else                                                5d10s24
                     iperm(1)=2                                         5d10s24
                     iperm(2)=1                                         5d10s24
                     iperm(3)=4                                         5d10s24
                     iperm(4)=3                                         5d10s24
                     iuseage(4)=iuseage(4)+1                            5d10s24
                    end if                                              5d10s24
                   end if                                               5d10s24
                  else                                                  5d10s24
                   if(isbv2.eq.lsa)then                                 5d10s24
                    if(isbv1.eq.lsc)then                                5d10s24
                     iperm(1)=3                                         5d10s24
                     iperm(2)=4                                         5d10s24
                     iperm(3)=1                                         5d10s24
                     iperm(4)=2                                         5d10s24
                     iuseage(5)=iuseage(5)+1                            5d10s24
                    else                                                5d10s24
                     iperm(1)=3                                         5d10s24
                     iperm(2)=4                                         5d10s24
                     iperm(3)=2                                         5d10s24
                     iperm(4)=1                                         5d10s24
                     iuseage(6)=iuseage(6)+1                            5d10s24
                    end if                                              5d10s24
                   else                                                 5d10s24
                    if(isbv1.eq.lsc)then                                5d10s24
                     iperm(1)=4                                         5d10s24
                     iperm(2)=3                                         5d10s24
                     iperm(3)=1                                         5d10s24
                     iperm(4)=2                                         5d10s24
                     iuseage(7)=iuseage(7)+1                            5d10s24
                    else                                                5d10s24
                     iperm(1)=4                                         5d10s24
                     iperm(2)=3                                         5d10s24
                     iperm(3)=2                                         5d10s24
                     iperm(4)=1                                         5d10s24
                     iuseage(8)=iuseage(8)+1                            5d10s24
                    end if                                              5d10s24
                   end if                                               5d10s24
                  end if                                                5d10s24
                  if(isbv12.eq.1.and.jsbv12.eq.1)then                    5d10s24
                   do iv=0,nvirt(isbv1)-1                               5d10s24
                    do jv=0,nvirt(jsbv1)-1                              5d10s24
                     sum=0d0                                            5d10s24
                     iiad=ioffvv+iv
                     jjad=joffvv+jv
                     ihitj=0
                     do i=0,nfdat(3,1,isb)-1                            5d10s24
                      sum=sum+vd(iiad)*vd(jjad)                         5d10s24
                      iiad=iiad+nvirt(isbv1)                            5d10s24
                      jjad=jjad+nvirt(jsbv1)                            5d10s24
                     end do                                             5d10s24
                     iiii(1)=iv+noc(isbv1)                               5d10s24
                     iiii(3)=iv+noc(isbv2)                              5d15s24
                     iiii(2)=jv+noc(jsbv1)                              5d10s24
                     iiii(4)=jv+noc(jsbv2)                              5d10s24
                     do it=0,nbasdws(isbv1)-1                             5d15s24
                      iiii(1)=it                                        5d15s24
                      iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                      jad=ie4v(isbv1,1)+iv+nvirt(isbv1)*it                5d15s24
                      orig=bc(igoal)
                      bc(jad)=bc(jad)+tmp(iad)*sum*0.5d0                5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalssj1 '),orig,tmp(iad),sum,bc(igoal),
     $                     lsa,lsb,lsc,lsd,iperm
                     end do                                             5d15s24
                     iiii(1)=iv+noc(isbv1)                               5d10s24
                     iiii(3)=iv+noc(isbv2)                               5d10s24
                     iiii(2)=jv+noc(jsbv1)                              5d10s24
                     iiii(4)=jv+noc(jsbv2)                              5d10s24
                     do it=0,nbasdws(isbv2)-1                             5d15s24
                      iiii(3)=it                                        5d15s24
                      iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                      jad=ie4v(isbv2,3)+iv+nvirt(isbv2)*it                5d15s24
                      orig=bc(igoal)
                      bc(jad)=bc(jad)+tmp(iad)*sum*0.5d0                5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalssj3 '),orig,tmp(iad),sum,bc(igoal),
     $                     lsa,lsb,lsc,lsd,iperm
                     end do                                             5d15s24
                     iiii(1)=iv+noc(isbv1)                               5d10s24
                     iiii(3)=iv+noc(isbv2)                               5d10s24
                     iiii(2)=jv+noc(jsbv1)                              5d10s24
                     iiii(4)=jv+noc(jsbv2)                              5d10s24
                     do it=0,nbasdws(jsbv1)-1                             5d15s24
                      iiii(2)=it                                        5d15s24
                      iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                      jad=ie4v(jsbv1,2)+jv+nvirt(jsbv1)*it                5d15s24
                      orig=bc(igoal)
                      bc(jad)=bc(jad)+tmp(iad)*sum*0.5d0                5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalssj2 '),orig,tmp(iad),sum,bc(igoal),
     $                     lsa,lsb,lsc,lsd,iperm
                     end do                                             5d15s24
                     iiii(1)=iv+noc(isbv1)                               5d10s24
                     iiii(3)=iv+noc(isbv2)                               5d10s24
                     iiii(2)=jv+noc(jsbv1)                              5d10s24
                     iiii(4)=jv+noc(jsbv2)                              5d10s24
                     do it=0,nbasdws(jsbv2)-1                             5d15s24
                      iiii(4)=it                                        5d15s24
                      iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                      jad=ie4v(jsbv2,4)+jv+nvirt(jsbv2)*it                5d15s24
                      orig=bc(igoal)
                      bc(jad)=bc(jad)+tmp(iad)*sum*0.5d0                5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalssj4 '),orig,tmp(iad),sum,bc(igoal),
     $                     lsa,lsb,lsc,lsd,iperm
                     end do                                             5d15s24
                    end do                                              5d10s24
                   end do                                               5d10s24
                  end if                                                5d10s24
                  iioff=ioff                                            5d10s24
                  jjoff=joff                                            5d10s24
                  if(jsbv12.eq.1)then                                   5d10s24
                   do jv=0,nvirt(jsbv1)-1                               5d10s24
                    do iv2=0,nvirt(isbv2)-1                             5d10s24
                     iv1top=iv2-1+isw*(nvirt(isbv1)-iv2)                5d10s24
                     do iv1=0,iv1top                                    5d10s24
                      itri=((iv2*(iv2-1))/2)+iv1                        5d10s24
                      irec=iv1+nvirt(isbv1)*iv2                         5d10s24
                      itri=itri+isw*(irec-itri)                         5d10s24
                      sum=0d0                                            5d10s24
                      iiad=ioff+itri                                    5d10s24
                      jjad=joffvv+jv                                    5d10s24
                      do i=0,nfdat(3,1,isb)-1                            5d10s24
                       sum=sum+vd(iiad)*vd(jjad)                         5d10s24
                       iiad=iiad+nvv                                    5d10s24
                       jjad=jjad+nvirt(jsbv1)                           5d10s24
                      end do                                            5d10s24
                      iiii(1)=iv1+noc(isbv1)                            5d10s24
                      iiii(3)=iv2+noc(isbv2)                             5d10s24
                      iiii(2)=jv+noc(jsbv1)                               5d10s24
                      iiii(4)=jv+noc(jsbv2)                               5d10s24
                      do it=0,nbasdws(isbv1)-1                            5d15s24
                       iiii(1)=it                                       5d15s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                      orig=bc(igoal)
                       jad=ie4v(isbv1,1)+iv1+nvirt(isbv1)*it              5d15s24
                       bc(jad)=bc(jad)+sum*srh*tmp(iad)                 5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalsnj1 '),orig,tmp(iad),sum,bc(igoal),
     $                      lsa,lsb,lsc,lsd,iperm
                      end do                                            5d15s24
                      iiii(1)=iv1+noc(isbv1)                            5d10s24
                      iiii(3)=iv2+noc(isbv2)                             5d10s24
                      iiii(2)=jv+noc(jsbv1)                               5d10s24
                      iiii(4)=jv+noc(jsbv2)                               5d10s24
                      do it=0,nbasdws(isbv2)-1                            5d15s24
                       iiii(3)=it                                       5d15s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       jad=ie4v(isbv2,3)+iv2+nvirt(isbv2)*it              5d15s24
                      orig=bc(igoal)
                       bc(jad)=bc(jad)+sum*srh*tmp(iad)                 5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalsnj3 '),orig,tmp(iad),sum,bc(igoal),
     $                      lsa,lsb,lsc,lsd,iperm
                      end do                                            5d15s24
                      iiii(1)=iv1+noc(isbv1)                            5d10s24
                      iiii(3)=iv2+noc(isbv2)                             5d10s24
                      iiii(2)=jv+noc(jsbv1)                               5d10s24
                      iiii(4)=jv+noc(jsbv2)                               5d10s24
                      do it=0,nbasdws(jsbv1)-1                            5d15s24
                       iiii(2)=it                                       5d15s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       jad=ie4v(jsbv1,2)+jv+nvirt(jsbv1)*it              5d15s24
                      orig=bc(igoal)
                       bc(jad)=bc(jad)+sum*srh*tmp(iad)                 5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalsnj2 '),orig,tmp(iad),sum,bc(igoal),
     $                      lsa,lsb,lsc,lsd,iperm
                      end do                                            5d15s24
                      iiii(1)=iv1+noc(isbv1)                            5d10s24
                      iiii(3)=iv2+noc(isbv2)                             5d10s24
                      iiii(2)=jv+noc(jsbv1)                               5d10s24
                      iiii(4)=jv+noc(jsbv2)                               5d10s24
                      do it=0,nbasdws(jsbv2)-1                            5d15s24
                       iiii(4)=it                                       5d15s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       jad=ie4v(jsbv2,4)+jv+nvirt(jsbv2)*it              5d15s24
                      orig=bc(igoal)
                       bc(jad)=bc(jad)+sum*srh*tmp(iad)                 5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalsnj4 '),orig,tmp(iad),sum,bc(igoal),
     $                      lsa,lsb,lsc,lsd,iperm
                      end do                                            5d15s24
                     end do                                             5d10s24
                    end do                                              5d10s24
                   end do                                               5d10s24
                  end if                                                5d10s24
                  if(isbv12.eq.1)then                                   5d10s24
                   do jv2=0,nvirt(jsbv2)-1                              5d10s24
                    jv1top=jv2-1+jsw*(nvirt(jsbv1)-jv2)                 5d10s24
                    do jv1=0,jv1top                                     5d10s24
                     jtri=((jv2*(jv2-1))/2)+jv1                         5d10s24
                     jrec=jv1+nvirt(jsbv1)*jv2                          5d10s24
                     jtri=jtri+jsw*(jrec-jtri)                          5d10s24
                     do iv=0,nvirt(isbv1)-1                             5d10s24
                      sum=0d0                                            5d10s24
                      iiad=ioffvv+iv                                    5d10s24
                      jjad=jjoff+jtri                                   5d10s24
                      do i=0,nfdat(3,1,isb)-1                            5d10s24
                       sum=sum+vd(iiad)*vd(jjad)                         5d10s24
                       iiad=iiad+nvirt(isbv1)                            5d10s24
                       jjad=jjad+mvv                                    5d10s24
                      end do                                            5d10s24
                      iiii(1)=iv+noc(isbv1)                             5d10s24
                      iiii(3)=iv+noc(isbv2)                             5d10s24
                      iiii(2)=jv1+noc(jsbv1)                             5d10s24
                      iiii(4)=jv2+noc(jsbv2)                              5d10s24
                      do it=0,nbasdws(isbv1)-1                            5d15s24
                       iiii(1)=it                                       5d15s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       jad=ie4v(isbv1,1)+iv+nvirt(isbv1)*it               5d15s24
                      orig=bc(igoal)
                       bc(jad)=bc(jad)+sum*srh*tmp(iad)                 5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalnsj1 '),orig,tmp(iad),sum,bc(igoal),
     $                      lsa,lsb,lsc,lsd,iperm
                      end do                                            5d15s24
                      iiii(1)=iv+noc(isbv1)                             5d10s24
                      iiii(3)=iv+noc(isbv2)                             5d10s24
                      iiii(2)=jv1+noc(jsbv1)                             5d10s24
                      iiii(4)=jv2+noc(jsbv2)                              5d10s24
                      do it=0,nbasdws(isbv2)-1                            5d15s24
                       iiii(3)=it                                       5d15s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       jad=ie4v(isbv2,3)+iv+nvirt(isbv2)*it               5d15s24
                      orig=bc(igoal)
                       bc(jad)=bc(jad)+sum*srh*tmp(iad)                 5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalnsj3 '),orig,tmp(iad),sum,bc(igoal),
     $                      lsa,lsb,lsc,lsd,iperm
                      end do                                            5d15s24
                      iiii(1)=iv+noc(isbv1)                             5d10s24
                      iiii(3)=iv+noc(isbv2)                             5d10s24
                      iiii(2)=jv1+noc(jsbv1)                             5d10s24
                      iiii(4)=jv2+noc(jsbv2)                              5d10s24
                      do it=0,nbasdws(jsbv1)-1                            5d15s24
                       iiii(2)=it                                       5d15s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       jad=ie4v(jsbv1,2)+jv1+nvirt(jsbv1)*it               5d15s24
                      orig=bc(igoal)
                       bc(jad)=bc(jad)+sum*srh*tmp(iad)                 5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalnsj2 '),orig,tmp(iad),sum,bc(igoal),
     $                      lsa,lsb,lsc,lsd,iperm
                      end do                                            5d15s24
                      iiii(1)=iv+noc(isbv1)                             5d10s24
                      iiii(3)=iv+noc(isbv2)                             5d10s24
                      iiii(2)=jv1+noc(jsbv1)                             5d10s24
                      iiii(4)=jv2+noc(jsbv2)                              5d10s24
                      do it=0,nbasdws(jsbv2)-1                            5d15s24
                       iiii(4)=it                                       5d15s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       jad=ie4v(jsbv2,4)+jv2+nvirt(jsbv2)*it            5d15s24
                      orig=bc(igoal)
                       bc(jad)=bc(jad)+sum*srh*tmp(iad)                 5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalnsj4 '),orig,tmp(iad),sum,bc(igoal),
     $                      lsa,lsb,lsc,lsd,iperm
                      end do                                            5d15s24
                     end do                                             5d10s24
                    end do                                              5d10s24
                   end do                                               5d10s24
                  end if                                                5d10s24
                  do l=1,4                                              5d10s24
                   do jv2=0,nvirt(jsbv2)-1                              5d10s24
                    jv1top=jv2-1+jsw*(nvirt(jsbv1)-jv2)                 5d10s24
                    do jv1=0,jv1top                                     5d10s24
                     jtri=((jv2*(jv2-1))/2)+jv1                         5d10s24
                     jrec=jv1+nvirt(jsbv1)*jv2                          5d10s24
                     jtri=jtri+jsw*(jrec-jtri)                          5d10s24
                     do iv2=0,nvirt(isbv2)-1                            5d10s24
                      iv1top=iv2-1+isw*(nvirt(isbv1)-iv2)               5d10s24
                      do iv1=0,iv1top                                   5d10s24
                       itri=((iv2*(iv2-1))/2)+iv1                       5d10s24
                       irec=iv1+nvirt(isbv1)*iv2                        5d10s24
                       itri=itri+isw*(irec-itri)                        5d10s24
                       sum=0d0                                            5d10s24
                       iiad=iioff+itri                                  5d10s24
                       jjad=jjoff+jtri                                  5d10s24
                       do i=0,nfdat(3,l,isb)-1                          5d10s24
                        orig=sum
                        sum=sum+vd(iiad)*vd(jjad)                         5d10s24
                        iiad=iiad+nvv                                   5d10s24
                        jjad=jjad+mvv                                   5d10s24
                       end do                                             5d10s24
                       iiii(1)=iv1+noc(isbv1)                           5d10s24
                       iiii(3)=iv2+noc(isbv2)                            5d10s24
                       iiii(2)=jv1+noc(jsbv1)                             5d10s24
                       iiii(4)=jv2+noc(jsbv2)                              5d10s24
                       do it=0,nbasdws(isbv1)-1                           5d15s24
                        iiii(1)=it                                      5d15s24
                        iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                        jad=ie4v(isbv1,1)+iv1+nvirt(isbv1)*it             5d15s24
                      orig=bc(igoal)
                        bc(jad)=bc(jad)+tmp(iad)*sum                    5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalnnj1 '),orig,tmp(iad),sum,bc(igoal),
     $                       lsa,lsb,lsc,lsd,iperm,iad+itmp-1,iv1,it,
     $                       iv1,jv1,iv2,jv2,iiii
                       end do                                           5d15s24
                       iiii(1)=iv1+noc(isbv1)                           5d10s24
                       iiii(3)=iv2+noc(isbv2)                            5d10s24
                       iiii(2)=jv1+noc(jsbv1)                             5d10s24
                       iiii(4)=jv2+noc(jsbv2)                              5d10s24
                       do it=0,nbasdws(isbv2)-1                           5d15s24
                        iiii(3)=it                                      5d15s24
                        iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                        jad=ie4v(isbv2,3)+iv2+nvirt(isbv2)*it             5d15s24
                      orig=bc(igoal)
                        bc(jad)=bc(jad)+tmp(iad)*sum                    5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalnnj3 '),orig,tmp(iad),sum,bc(igoal),
     $                       lsa,lsb,lsc,lsd,iperm
                       end do                                           5d15s24
                       iiii(1)=iv1+noc(isbv1)                           5d10s24
                       iiii(3)=iv2+noc(isbv2)                            5d10s24
                       iiii(2)=jv1+noc(jsbv1)                             5d10s24
                       iiii(4)=jv2+noc(jsbv2)                              5d10s24
                       do it=0,nbasdws(jsbv1)-1                           5d15s24
                        iiii(2)=it                                      5d15s24
                        iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                        jad=ie4v(jsbv1,2)+jv1+nvirt(jsbv1)*it             5d15s24
                      orig=bc(igoal)
                        bc(jad)=bc(jad)+tmp(iad)*sum                    5d15s24
                        if(abs(orig-bc(igoal)).gt.1d-10)then
                         write(6,*)
     $                     ('ggoalnnj2 '),orig,tmp(iad),sum,bc(igoal),
     $                       lsa,lsb,lsc,lsd,iperm
                       end if
                       end do                                           5d15s24
                       iiii(1)=iv1+noc(isbv1)                           5d10s24
                       iiii(3)=iv2+noc(isbv2)                            5d10s24
                       iiii(2)=jv1+noc(jsbv1)                             5d10s24
                       iiii(4)=jv2+noc(jsbv2)                              5d10s24
                       do it=0,nbasdws(jsbv2)-1                           5d15s24
                        iiii(4)=it                                      5d15s24
                        iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                        jad=ie4v(jsbv2,4)+jv2+nvirt(jsbv2)*it             5d15s24
                      orig=bc(igoal)
                        bc(jad)=bc(jad)+tmp(iad)*sum                    5d15s24
                        if(abs(orig-bc(igoal)).gt.1d-10)then
                         write(6,*)
     $                     ('ggoalnnj4 '),orig,tmp(iad),sum,bc(igoal),
     $                       lsa,lsb,lsc,lsd,iperm
                       end if
                       end do                                           5d15s24
                      end do
                     end do                                             5d10s24
                    end do                                              5d10s24
                   end do                                               5d10s24
                   iioff=iioff+nvv*nfdat(3,l,isb)                       5d10s24
                   jjoff=jjoff+mvv*nfdat(3,l,isb)                       5d10s24
                  end do                                                5d10s24
                 end if                                                 5d10s24
                 end if                                                 6d11s24
                 if(ldok)then                                           6d11s24
                 ix=max(i12r,i34r)
                 in=min(i12r,i34r)
                 if(ix.eq.lcd.and.in.eq.lab)then                        5d10s24
                  if(i12r.eq.lab)then                                   5d13s24
                   if(isbv1.eq.lsa)then                                 5d10s24
                    if(isbv2.eq.lsc)then
                     iperm(1)=1                                         5d10s24
                     iperm(2)=2                                         5d10s24
                     iperm(3)=3                                         5d10s24
                     iperm(4)=4                                         5d10s24
                     iuseage(1)=iuseage(1)+1                            5d10s24
                    else                                                5d10s24
                     iperm(1)=1                                         5d10s24
                     iperm(2)=2                                         5d10s24
                     iperm(3)=4                                         5d10s24
                     iperm(4)=3                                         5d10s24
                     iuseage(2)=iuseage(2)+1                            5d10s24
                    end if                                              5d10s24
                   else                                                 5d10s24
                    if(isbv2.eq.lsc)then                                5d10s24
                     iperm(1)=2                                         5d10s24
                     iperm(2)=1                                         5d10s24
                     iperm(3)=3                                         5d10s24
                     iperm(4)=4                                         5d10s24
                     iuseage(3)=iuseage(3)+1                            5d10s24
                    else                                                5d10s24
                     iperm(1)=2                                         5d10s24
                     iperm(2)=1                                         5d10s24
                     iperm(3)=4                                         5d10s24
                     iperm(4)=3                                         5d10s24
                     iuseage(4)=iuseage(4)+1                            5d10s24
                    end if                                              5d10s24
                   end if                                               5d10s24
                  else                                                  5d10s24
                   if(isbv2.eq.lsa)then                                 5d10s24
                    if(isbv1.eq.lsc)then                                5d10s24
                     iperm(1)=3                                         5d10s24
                     iperm(2)=4                                         5d10s24
                     iperm(3)=1                                         5d10s24
                     iperm(4)=2                                         5d10s24
                     iuseage(5)=iuseage(5)+1                            5d10s24
                    else                                                5d10s24
                     iperm(1)=3                                         5d10s24
                     iperm(2)=4                                         5d10s24
                     iperm(3)=2                                         5d10s24
                     iperm(4)=1                                         5d10s24
                     iuseage(6)=iuseage(6)+1                            5d10s24
                    end if                                              5d10s24
                   else                                                 5d10s24
                    if(isbv1.eq.lsc)then                                5d10s24
                     iperm(1)=4                                         5d10s24
                     iperm(2)=3                                         5d10s24
                     iperm(3)=1                                         5d10s24
                     iperm(4)=2                                         5d10s24
                     iuseage(7)=iuseage(7)+1                            5d10s24
                    else                                                5d10s24
                     iperm(1)=4                                         5d10s24
                     iperm(2)=3                                         5d10s24
                     iperm(3)=2                                         5d10s24
                     iperm(4)=1                                         5d10s24
                     iuseage(8)=iuseage(8)+1                            5d10s24
                    end if                                              5d10s24
                   end if                                               5d10s24
                  end if                                                5d10s24
                  if(isbv12.eq.1.and.jsbv12.eq.1)then                    5d10s24
                   do iv=0,nvirt(isbv1)-1                               5d10s24
                    do jv=0,nvirt(jsbv1)-1                              5d10s24
                     sum=0d0                                            5d10s24
                     iiad=ioffvv+iv
                     jjad=joffvv+jv
                     do i=0,nfdat(3,1,isb)-1                            5d10s24
                      sum=sum+vd(iiad)*vd(jjad)                         5d10s24
                      iiad=iiad+nvirt(isbv1)                            5d10s24
                      jjad=jjad+nvirt(jsbv1)                            5d10s24
                     end do                                             5d10s24
                     iiii(1)=iv+noc(isbv1)                               5d10s24
                     iiii(3)=iv+noc(isbv2)                               5d10s24
                     iiii(2)=jv+noc(jsbv2)                              5d10s24
                     iiii(4)=jv+noc(jsbv1)                              5d10s24
                     do it=0,nbasdws(isbv1)-1                             5d15s24
                      iiii(1)=it                                        5d15s24
                      iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                      jad=ie4v(isbv1,1)+iv+nvirt(isbv1)*it                5d15s24
                      orig=bc(igoal)
                      bc(jad)=bc(jad)+tmp(iad)*sum*0.5d0                5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalssk1 '),orig,tmp(iad),sum,bc(igoal),
     $                     lsa,lsb,lsc,lsd,iperm
                     end do                                             5d15s24
                     iiii(1)=iv+noc(isbv1)                               5d10s24
                     iiii(3)=iv+noc(isbv2)                               5d10s24
                     iiii(2)=jv+noc(jsbv2)                              5d10s24
                     iiii(4)=jv+noc(jsbv1)                              5d10s24
                     do it=0,nbasdws(isbv2)-1                             5d15s24
                      iiii(3)=it                                        5d15s24
                      iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                      jad=ie4v(isbv2,3)+iv+nvirt(isbv2)*it                5d15s24
                      orig=bc(igoal)
                      bc(jad)=bc(jad)+tmp(iad)*sum*0.5d0                5d15s24
                     end do                                             5d15s24
                     iiii(1)=iv+noc(isbv1)                               5d10s24
                     iiii(3)=iv+noc(isbv2)                               5d10s24
                     iiii(2)=jv+noc(jsbv2)                              5d10s24
                     iiii(4)=jv+noc(jsbv1)                              5d10s24
                     do it=0,nbasdws(jsbv2)-1                             5d15s24
                      iiii(2)=it                                        5d15s24
                      iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                      jad=ie4v(jsbv2,2)+jv+nvirt(jsbv2)*it                5d15s24
                      orig=bc(igoal)
                      bc(jad)=bc(jad)+tmp(iad)*sum*0.5d0                5d15s24
                     end do                                             5d15s24
                     iiii(1)=iv+noc(isbv1)                               5d10s24
                     iiii(3)=iv+noc(isbv2)                               5d10s24
                     iiii(2)=jv+noc(jsbv2)                              5d10s24
                     iiii(4)=jv+noc(jsbv1)                              5d10s24
                     do it=0,nbasdws(jsbv1)-1                             5d15s24
                      iiii(4)=it                                        5d15s24
                      iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                      jad=ie4v(jsbv1,4)+jv+nvirt(jsbv1)*it                5d15s24
                      orig=bc(igoal)
                      bc(jad)=bc(jad)+tmp(iad)*sum*0.5d0                5d15s24
                     end do                                             5d15s24
                    end do                                              5d10s24
                   end do                                               5d10s24
                  end if                                                5d10s24
                  iioff=ioff                                            5d10s24
                  jjoff=joff                                            5d10s24
                  if(jsbv12.eq.1)then                                   5d10s24
                   do jv=0,nvirt(jsbv1)-1                               5d10s24
                    do iv2=0,nvirt(isbv2)-1                             5d10s24
                     iv1top=iv2-1+isw*(nvirt(isbv1)-iv2)                5d10s24
                     do iv1=0,iv1top                                    5d10s24
                      itri=((iv2*(iv2-1))/2)+iv1                        5d10s24
                      irec=iv1+nvirt(isbv1)*iv2                         5d10s24
                      itri=itri+isw*(irec-itri)                         5d10s24
                      sum=0d0                                            5d10s24
                      iiad=ioff+itri                                    6d11s24
                      jjad=joffvv+jv                                    5d10s24
                      do i=0,nfdat(3,1,isb)-1                            5d10s24
                       sum=sum+vd(iiad)*vd(jjad)                         5d10s24
                       iiad=iiad+nvv                                    5d10s24
                       jjad=jjad+nvirt(jsbv1)                           5d10s24
                      end do                                            5d10s24
                      iiii(1)=iv1+noc(isbv1)                            5d10s24
                      iiii(3)=iv2+noc(isbv2)                             5d10s24
                      iiii(4)=jv+noc(jsbv1)                               5d10s24
                      iiii(2)=jv+noc(jsbv2)                               5d10s24
                      do it=0,nbasdws(isbv1)-1                            5d15s24
                       iiii(1)=it                                       5d15s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       jad=ie4v(isbv1,1)+iv1+nvirt(isbv1)*it              5d15s24
                      orig=bc(igoal)
                       bc(jad)=bc(jad)+sum*srh*tmp(iad)                 5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalsnk1 '),orig,tmp(iad),sum,bc(igoal),
     $                      lsa,lsb,lsc,lsd,iperm
                      end do                                            5d15s24
                      iiii(1)=iv1+noc(isbv1)                            5d10s24
                      iiii(3)=iv2+noc(isbv2)                             5d10s24
                      iiii(4)=jv+noc(jsbv1)                               5d10s24
                      iiii(2)=jv+noc(jsbv2)                               5d10s24
                      do it=0,nbasdws(isbv2)-1                            5d15s24
                       iiii(3)=it                                       5d15s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       jad=ie4v(isbv2,3)+iv2+nvirt(isbv2)*it              5d15s24
                      orig=bc(igoal)
                       bc(jad)=bc(jad)+sum*srh*tmp(iad)                 5d15s24
                      end do                                            5d15s24
                      iiii(1)=iv1+noc(isbv1)                            5d10s24
                      iiii(3)=iv2+noc(isbv2)                             5d10s24
                      iiii(4)=jv+noc(jsbv1)                               5d10s24
                      iiii(2)=jv+noc(jsbv2)                               5d10s24
                      do it=0,nbasdws(jsbv2)-1                            5d15s24
                       iiii(2)=it                                       5d15s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       jad=ie4v(jsbv2,2)+jv+nvirt(jsbv2)*it              5d15s24
                      orig=bc(igoal)
                       bc(jad)=bc(jad)+sum*srh*tmp(iad)                 5d15s24
                      end do                                            5d15s24
                      iiii(1)=iv1+noc(isbv1)                            5d10s24
                      iiii(3)=iv2+noc(isbv2)                             5d10s24
                      iiii(4)=jv+noc(jsbv1)                               5d10s24
                      iiii(2)=jv+noc(jsbv2)                               5d10s24
                      do it=0,nbasdws(jsbv1)-1                            5d15s24
                       iiii(4)=it                                       5d15s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       jad=ie4v(jsbv1,4)+jv+nvirt(jsbv1)*it              5d15s24
                      orig=bc(igoal)
                       bc(jad)=bc(jad)+sum*srh*tmp(iad)                 5d15s24
                      end do                                            5d15s24
                     end do                                             5d10s24
                    end do                                              5d10s24
                   end do                                               5d10s24
                  end if                                                5d10s24
                  if(isbv12.eq.1)then                                   5d10s24
                   do jv2=0,nvirt(jsbv2)-1                              5d10s24
                    jv1top=jv2-1+jsw*(nvirt(jsbv1)-jv2)                 5d10s24
                    do jv1=0,jv1top                                     5d10s24
                     jtri=((jv2*(jv2-1))/2)+jv1                         5d10s24
                     jrec=jv1+nvirt(jsbv1)*jv2                          5d10s24
                     jtri=jtri+jsw*(jrec-jtri)                          5d10s24
                     do iv=0,nvirt(isbv1)-1                             5d10s24
                      sum=0d0                                            5d10s24
                      iiad=ioffvv+iv                                    5d10s24
                      jjad=jjoff+jtri                                   5d10s24
                      do i=0,nfdat(3,1,isb)-1                            5d10s24
                       sum=sum+vd(iiad)*vd(jjad)                         5d10s24
                       iiad=iiad+nvirt(isbv1)                            5d10s24
                       jjad=jjad+mvv                                    5d10s24
                      end do                                            5d10s24
                      iiii(1)=iv+noc(isbv1)                             5d10s24
                      iiii(3)=iv+noc(isbv2)                             5d10s24
                      iiii(2)=jv2+noc(jsbv2)                              5d10s24
                      iiii(4)=jv1+noc(jsbv1)                             5d10s24
                      do it=0,nbasdws(isbv1)-1                            5d15s24
                       iiii(1)=it                                       5d15s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       jad=ie4v(isbv1,1)+iv+nvirt(isbv1)*it               5d15s24
                      orig=bc(igoal)
                       bc(jad)=bc(jad)+sum*srh*tmp(iad)                 5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                     ('ggoalnsk1 '),orig,tmp(iad),sum,bc(igoal),
     $                      lsa,lsb,lsc,lsd,iperm
                      end do                                            5d15s24
                      iiii(1)=iv+noc(isbv1)                             5d10s24
                      iiii(3)=iv+noc(isbv2)                             5d10s24
                      iiii(2)=jv2+noc(jsbv2)                              5d10s24
                      iiii(4)=jv1+noc(jsbv1)                             5d10s24
                      do it=0,nbasdws(isbv2)-1                            5d15s24
                       iiii(3)=it                                       5d15s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       jad=ie4v(isbv2,3)+iv+nvirt(isbv2)*it               5d15s24
                      orig=bc(igoal)
                       bc(jad)=bc(jad)+sum*srh*tmp(iad)                 5d15s24
                      end do                                            5d15s24
                      iiii(1)=iv+noc(isbv1)                             5d10s24
                      iiii(3)=iv+noc(isbv2)                             5d10s24
                      iiii(2)=jv2+noc(jsbv2)                              5d10s24
                      iiii(4)=jv1+noc(jsbv1)                             5d10s24
                      do it=0,nbasdws(jsbv2)-1                            5d15s24
                       iiii(2)=it                                       5d15s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       jad=ie4v(jsbv2,2)+jv2+nvirt(jsbv2)*it               5d15s24
                      orig=bc(igoal)
                       bc(jad)=bc(jad)+sum*srh*tmp(iad)                 5d15s24
                      end do                                            5d15s24
                      iiii(1)=iv+noc(isbv1)                             5d10s24
                      iiii(3)=iv+noc(isbv2)                             5d10s24
                      iiii(2)=jv2+noc(jsbv2)                              5d10s24
                      iiii(4)=jv1+noc(jsbv1)                             5d10s24
                      do it=0,nbasdws(jsbv1)-1                            5d15s24
                       iiii(4)=it                                       5d15s24
                       iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                       jad=ie4v(jsbv1,4)+jv1+nvirt(jsbv1)*it               5d15s24
                      orig=bc(igoal)
                       bc(jad)=bc(jad)+sum*srh*tmp(iad)                 5d15s24
                      end do                                            5d15s24
                     end do                                             5d10s24
                    end do                                              5d10s24
                   end do                                               5d10s24
                  end if                                                5d10s24
                  tf=1d0                                                5d10s24
                  do l=1,4                                              5d10s24
                   do jv2=0,nvirt(jsbv2)-1                              5d10s24
                    jv1top=jv2-1+jsw*(nvirt(jsbv1)-jv2)                 5d10s24
                    do jv1=0,jv1top                                     5d10s24
                     jtri=((jv2*(jv2-1))/2)+jv1                         5d10s24
                     jrec=jv1+nvirt(jsbv1)*jv2                          5d10s24
                     jtri=jtri+jsw*(jrec-jtri)                          5d10s24
                     do iv2=0,nvirt(isbv2)-1                            5d10s24
                      iv1top=iv2-1+isw*(nvirt(isbv1)-iv2)               5d10s24
                      do iv1=0,iv1top                                   5d10s24
                       itri=((iv2*(iv2-1))/2)+iv1                       5d10s24
                       irec=iv1+nvirt(isbv1)*iv2                        5d10s24
                       itri=itri+isw*(irec-itri)                        5d10s24
                       sum=0d0                                            5d10s24
                       iiad=iioff+itri                                  5d10s24
                       jjad=jjoff+jtri                                  5d10s24
                       do i=0,nfdat(3,l,isb)-1                          5d10s24
                        orig=sum
                        sum=sum+vd(iiad)*vd(jjad)                         5d10s24
                        iiad=iiad+nvv                                   5d10s24
                        jjad=jjad+mvv                                   5d10s24
                       end do                                             5d10s24
                       iiii(1)=iv1+noc(isbv1)                           5d10s24
                       iiii(3)=iv2+noc(isbv2)                            5d10s24
                       iiii(4)=jv1+noc(jsbv1)                             5d10s24
                       iiii(2)=jv2+noc(jsbv2)                              5d10s24
                       do it=0,nbasdws(isbv1)-1                         5d15s24
                        iiii(1)=it                                      5d15s24
                        iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                        jad=ie4v(isbv1,1)+iv1+nvirt(isbv1)*it             5d15s24
                      orig=bc(igoal)
                        bc(jad)=bc(jad)+tmp(iad)*sum*tf                 5d15s24
                      if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)
     $                    ('ggoalnnk1 '),orig,tmp(iad),sum*tf,bc(igoal),
     $                       lsa,lsb,lsc,lsd,iperm,iv1,jv2,iv2,jv1,iiii,
     $                       tf
                       end do                                           5d15s24
                       iiii(1)=iv1+noc(isbv1)                           5d10s24
                       iiii(3)=iv2+noc(isbv2)                            5d10s24
                       iiii(4)=jv1+noc(jsbv1)                             5d10s24
                       iiii(2)=jv2+noc(jsbv2)                              5d10s24
                       do it=0,nbasdws(isbv2)-1                         5d15s24
                        iiii(3)=it                                      5d15s24
                        iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                        jad=ie4v(isbv2,3)+iv2+nvirt(isbv2)*it             5d15s24
                      orig=bc(igoal)
                        bc(jad)=bc(jad)+tmp(iad)*sum*tf                 5d15s24
                       end do                                           5d15s24
                       iiii(1)=iv1+noc(isbv1)                           5d10s24
                       iiii(3)=iv2+noc(isbv2)                            5d10s24
                       iiii(4)=jv1+noc(jsbv1)                             5d10s24
                       iiii(2)=jv2+noc(jsbv2)                              5d10s24
                       do it=0,nbasdws(jsbv1)-1                         5d15s24
                        iiii(4)=it                                      5d15s24
                        iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                        jad=ie4v(jsbv1,4)+jv1+nvirt(jsbv1)*it             5d15s24
                      orig=bc(igoal)
                        bc(jad)=bc(jad)+tmp(iad)*sum*tf                 5d15s24
                       end do                                           5d15s24
                       iiii(1)=iv1+noc(isbv1)                           5d10s24
                       iiii(3)=iv2+noc(isbv2)                            5d10s24
                       iiii(4)=jv1+noc(jsbv1)                             5d10s24
                       iiii(2)=jv2+noc(jsbv2)                              5d10s24
                       do it=0,nbasdws(jsbv2)-1                         5d15s24
                        iiii(2)=it                                      5d15s24
                        iad=1+iiii(iperm(1))+nbasdws(lsa)*                 5d10s24
     $                    (iiii(iperm(2))+nbasdws(lsb)*                 5d10s24
     $                    (iiii(iperm(3))+nbasdws(lsc)*iiii(iperm(4)))) 5d10s24
                        jad=ie4v(jsbv2,2)+jv2+nvirt(jsbv2)*it             5d15s24
                      orig=bc(igoal)
                        bc(jad)=bc(jad)+tmp(iad)*sum*tf                 5d15s24
                       end do                                           5d15s24
                      end do
                     end do                                             5d10s24
                    end do                                              5d10s24
                   end do                                               5d10s24
                   iioff=iioff+nvv*nfdat(3,l,isb)                       5d10s24
                   jjoff=jjoff+mvv*nfdat(3,l,isb)                       5d10s24
                   tf=-1d0                                              5d10s24
                  end do                                                5d10s24
                 end if                                                 5d10s24
                 end if                                                 6d11s24
                end if
                do l=1,4                                                5d10s24
                 joff=joff+nfdat(3,l,jsb)*mvv                           5d10s24
                end do                                                  5d10s24
               end if                                                   5d10s24
              end do                                                    5d10s24
             end do                                                     5d10s24
             do l=1,4                                                   5d10s24
              ioff=ioff+nvv*nfdat(3,l,isb)                              5d10s24
             end do                                                     5d10s24
            end if                                                      5d10s24
           end do                                                       5d10s24
          end do                                                        5d10s24
          return                                                        5d2s24
          end                                                           5d2s24
