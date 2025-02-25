c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hc1cd(ivecb,ivectb,ivecab,iveck,ivectk,ivecak,ilc,ihc, 5d27s22
     $     ilct,ihct,nroot,nherec,                                      5d27s22
     $                nherect,iaorb,iborb,nalpha,nbeta,                 3d27s17
     $                m12sym,norb,nsymb,idata1,idata2,idatb1,idatb2,    3d10s17
     $             nsbeta,m1c,idatac,idatbc,iden1,jdenpt,wgt,l2e,bc,ibc,3d15s23
     $     mdenoff,igoal)                                                     3d15s23
      implicit real*8 (a-h,o-z)
      integer*2 idwstest2(4)                                            12d13s06
      integer*8 idwstest8                                               12d13s06
      logical l2e                                                       6d2s22
      equivalence (idwstest8,idwstest2)                                 12d13s06
c
c     parallel formation of density matrix from h*c
c
      dimension ivecb(8),ivectb(8),iden1(8),jdenpt(*),                    3d27s17
     $     nherec(8),nsbeta(8),idata1(8),idata2(8),wgt(nroot),          3d27s17
     $     idatb1(8),idatb2(8),ilct(8),ihct(8),nherect(8),              3d10s17
     $     ilc(8),ihc(8),m12sym(8,2,2),ivecab(8),idatac(36),idatbc(36),  3d15s17
     $     m1c(36,2),iveck(*),ivectk(*),ivecak(*)                       5d27s22
      include "common.cas"
      include "common.store"
      ism=ibcoff                                                        8d29s06
      irelo=ism+norb                                                    8d29s06
      ibcoff=irelo+norb                                                 8d29s06
      call enough('hc1cd.  1',bc,ibc)
      jsm=ism-1                                                         8d29s06
      jrelo=irelo-1                                                     8d29s06
      do i=1,nsymb
       do j=1,iacto(i)                                                  8d29s06
        ibc(jsm+j)=i                                                    8d29s06
        ibc(jrelo+j)=j                                                  8d29s06
       end do                                                           8d29s06
       jsm=jsm+iacto(i)                                                 8d29s06
       jrelo=jrelo+iacto(i)                                             8d29s06
      end do                                                            8d29s06
c
c     beta and alpha single and double differences
c
      do isb=1,nsymb                                                    8d29s06
       kadd=ilc(isb)-1+ioffdeta(isb)                                    12d20s06
       jsb=nsbeta(isb)                                                  3d14s17
       numbx=nbdet(jsb)                                                 3d14s17
       call dosingd(ibc(idatb1(jsb)),m12sym(jsb,1,2),bc(ivecb(isb)),    5d27s22
     $     bc(iveck(isb)),numbx,nroot,nherec(isb),iborb,nbeta,          5d27s22
     $       iaorb,nalpha,kadd,ibc(ism),ibc(irelo),iacto,               12d20s06
     $       ioffdetb(jsb),iden1,jdenpt,wgt,l2e,bc,ibc,mdenoff,igoal)         3d15s23
       if(l2e)call dodoubd(ibc(idatb2(jsb)),m12sym(jsb,2,2),            6d2s22
     $      bc(ivecb(isb)),bc(iveck(isb)),numbx,nroot,nherec(isb),iacto,6d2s22
     $      ibc(ism),ibc(irelo),jdenpt,wgt,bc,ibc)                      11d14s22
       kadd=ilct(isb)-1+ioffdetb(isb)                                    12d20s06
       numbx=nadet(nsbeta(isb))                                         3d10s17
       call dosingd(ibc(idata1(jsb)),m12sym(jsb,1,1),bc(ivectb(isb)),     3d10s17
     $       bc(ivectk(isb)),numbx,nroot,nherect(isb),iaorb,nalpha,     5d27s22
     $       iborb,nbeta,kadd,ibc(ism),ibc(irelo),iacto,                3d14s17
     $       ioffdeta(nsbeta(isb)),iden1,jdenpt,wgt,l2e,bc,ibc,mdenoff, 3d16s23
     $      igoal)                                                      3d16s23
       if(l2e)call dodoubd(ibc(idata2(jsb)),m12sym(jsb,2,1),            6d2s22
     $      bc(ivectb(isb)),bc(ivectk(isb)),numbx,nroot,nherect(isb),   6d2s22
     $      iacto,ibc(ism),ibc(irelo),jdenpt,wgt,bc,ibc)                11d14s22
      end do
c
c     alpha-beta terms
c
      if(l2e)then                                                       6d2s22
        do isb=1,nsymb                                                    3d10s17
        jsb=nsbeta(isb)                                                  3d14s17
        numbx=nbdet(jsb)                                                 3d14s17
        call doab1d(ibc(idatb1(jsb)),m12sym(jsb,1,2),ibc(idata1(isb)),    3d15s17
     $      m12sym(isb,1,1),bc(ivecab(isb)),bc(ivecak(isb)),nadet(isb), 5d27s22
     $      numbx,nroot,nherec(isb),ilc(isb),ihc(isb),ibc(irelo),       5d27s22
     $      ibc(ism),iacto,jdenpt,wgt,bc,ibc)                           11d14s22
        do isk=1,isb-1                                                   3d15s17
         iada=((isb*(isb-1))/2)+isk                                      3d15s17
         jsk=nsbeta(isk)                                                 3d15s17
         ix=max(jsb,jsk)                                                 3d15s17
         in=min(jsb,jsk)                                                 3d15s17
         iadb=((ix*(ix-1))/2)+in                                         3d15s17
         if(m1c(iada,1)*m1c(iadb,2).gt.0)then                            3d15s17
          if(ix.eq.jsb)then
           iu1=1
           iu2=2
          else
           iu1=2
           iu2=1
          end if                                                         3d16s17
          call doab2d(ibc(idatac(iada)),m1c(iada,1),ibc(idatbc(iadb)),   3d27s17
     $        m1c(iadb,2),nadet(isb),                                   3d16s17
     $        nadet(isk),nbdet(jsb),nbdet(jsk),nroot,                   3d27s17
     $        nherec(isb),ilc(isb),ihc(isb),nherec(isk),                3d27s17
     $        ilc(isk),ihc(isk),bc(ivecab(isb)),bc(ivecak(isk)),        5d27s22
     $        bc(ivecab(jsb)),bc(ivecak(jsk)),ibc(irelo),ibc(ism),      5d27s22
     $        iacto,iu1,iu2,jdenpt,wgt,bc,ibc)                          12d2s22
         end if                                                          3d15s17
        end do                                                           3d15s17
       end do                                                            3d10s17
      end if                                                            6d2s22
      ibcoff=ism                                                        3d10s17
      return
      end
