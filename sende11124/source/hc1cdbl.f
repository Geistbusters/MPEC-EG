c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hc1cdbl(ivec,ivect,iveca,ilc,ihc,ilct,ihct,nherec,
     $                nherect,iaorb,iborb,nalpha,nbeta,                 3d27s17
     $                m12sym,norb,nsymb,idata1,idata2,idatb1,idatb2,    3d10s17
     $             nsbeta,m1c,idatac,idatbc,jdenpt,bc,ibc,nroot,mdenoff,3d15s23
     $     iden1,igoal)                                                       3d15s23
      implicit real*8 (a-h,o-z)
      integer*2 idwstest2(4)                                            12d13s06
      integer*8 idwstest8                                               12d13s06
      equivalence (idwstest8,idwstest2)                                 12d13s06
c
c     parallel formation of bi-linear density matrix from h*c
c
      dimension ivec(*),ivect(*),jdenpt(*),                             7d21s22
     $     nherec(8),nsbeta(8),idata1(8),idata2(8),                     7d21s22
     $     idatb1(8),idatb2(8),ilct(8),ihct(8),nherect(8),              3d10s17
     $     ilc(8),ihc(8),m12sym(8,2,2),iveca(8),idatac(36),idatbc(36),  3d15s17
     $     m1c(36,2),iden1(*)                                           3d15s23
      include "common.cas"
      include "common.store"
      ism=ibcoff                                                        8d29s06
      irelo=ism+norb                                                    8d29s06
      ibcoff=irelo+norb                                                 8d29s06
      call enough('hc1cdbl.  1',bc,ibc)
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
       call dosingdbl(ibc(idatb1(jsb)),m12sym(jsb,1,2),bc(ivec(isb)),   7d21s22
     $       numbx,nherec(isb),iborb,nbeta,                             7d21s22
     $       iaorb,nalpha,kadd,ibc(ism),ibc(irelo),iacto,               12d20s06
     $       ioffdetb(jsb),jdenpt,bc,ibc,nroot,mdenoff,iden1,igoal)           3d15s23
       call dodoubdbl(ibc(idatb2(jsb)),m12sym(jsb,2,2),bc(ivec(isb)),   7d21s22
     $      numbx,nherec(isb),iacto,ibc(ism),ibc(irelo),jdenpt,bc,ibc,  3d14s23
     $      nroot,mdenoff,igoal)                                        3d15s23
       kadd=ilct(isb)-1+ioffdetb(isb)                                    12d20s06
       numbx=nadet(nsbeta(isb))                                         3d10s17
       call dosingdbl(ibc(idata1(jsb)),m12sym(jsb,1,1),bc(ivect(isb)),  7d21s22
     $       numbx,nherect(isb),iaorb,nalpha,                           7d21s22
     $       iborb,nbeta,kadd,ibc(ism),ibc(irelo),iacto,                3d14s17
     $       ioffdeta(nsbeta(isb)),jdenpt,bc,ibc,nroot,mdenoff,iden1,   3d15s23
     $      igoal)                                                      3d15s23
       call dodoubdbl(ibc(idata2(jsb)),m12sym(jsb,2,1),bc(ivect(isb)),  7d21s22
     $      numbx,nherect(isb),iacto,ibc(ism),ibc(irelo),jdenpt,bc,ibc, 3d14s23
     $      nroot,mdenoff,igoal)                                        3d15s23
      end do
c
c     alpha-beta terms
c
      do isb=1,nsymb                                                    3d10s17
       jsb=nsbeta(isb)                                                  3d14s17
       numbx=nbdet(jsb)                                                 3d14s17
       call doab1dbl(ibc(idatb1(jsb)),m12sym(jsb,1,2),ibc(idata1(isb)), 7d21s22
     $      m12sym(isb,1,1),bc(iveca(isb)),nadet(isb),numbx,nherec(isb),7d21s22
     $      ilc(isb),ihc(isb),ibc(irelo),ibc(ism),iacto,jdenpt,bc,ibc,  3d15s23
     $      nroot,mdenoff,igoal)                                        3d15s23
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
         ff=1d0
          call doab2dbl(ibc(idatac(iada)),m1c(iada,1),ibc(idatbc(iadb)),7d22s22
     $        m1c(iadb,2),nadet(isb),nadet(isk),nbdet(jsb),nbdet(jsk),  7d22s22
     $        nherec(isb),ilc(isb),ihc(isb),nherec(isk),                3d27s17
     $        ilc(isk),ihc(isk),bc(iveca(isb)),bc(iveca(isk)),          7d22s22
     $        ibc(irelo),ibc(ism),iacto,iu1,iu2,jdenpt,bc,ibc,nroot,    3d15s23
     $        mdenoff,igoal,ff)                                         3d27s23
        end if                                                          3d15s17
       end do                                                           3d15s17
      end do                                                            3d10s17
      ibcoff=ism                                                        3d10s17
      return
      end
