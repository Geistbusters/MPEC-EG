c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hc1cnona1(ivec,ivect,iveca,ig,igt,nroot,ih0,nherec,    6d10s22
     $          ilc,ihc,nherect,ilct,ihct,iaorb,iborb,nalpha,nbeta,i2e, 6d13s22
     $                nsblkder,isblkder,norb,nsymb,ipuse,multh,         6d13s22
     $     m12symnon,idata1non,idata2non,idatb1non,idatb2non,ipxder,    6d13s22
     $     nsbeta,m12sym,idata1,idatb1,m1c,idatac,idatbc,bc,ibc,        3d27s23
     $     igoal)
      implicit real*8 (a-h,o-z)
      integer*2 idwstest2(4)                                            12d13s06
      integer*8 idwstest8                                               12d13s06
      logical lbail                                                     5d10s18
      equivalence (idwstest8,idwstest2)                                 12d13s06
      external second                                                   4d9s18
c
c     parallel formation of h*c
c     this version is for non-totally symmetric interactions.           6d10s22
c
      dimension ig(8),igt(8),ivec(8),ivect(8),ih0(8),multh(8,8),        6d10s22
     $     i2e(*),isblkder(4,*),nherec(8),nsbeta(8),idata1(8),idata2(8),6d10s22
     $     idatb1(8),idatb2(8),ilct(8),ihct(8),nherect(8),              3d10s17
     $     ilc(8),ihc(8),m12sym(8,2,2),iveca(8),idatac(36),idatbc(36),  3d15s17
     $     m1c(36,2),mst(8,8,2),mnz(3,64,36,2),nz(36,2),
     $     m12symnon(8,2,2),idata1non(*),idata2non(*),idatb1non(*),     6d13s22
     $     idatb2non(*),ipxder(4,8,8,8)                                 6d28s22
      include "common.cas"
      include "common.store"
      data icall/0/
      save
      ism=ibcoff                                                        8d29s06
      irelo=ism+norb                                                    8d29s06
      ibcoff=irelo+norb                                                 8d29s06
      call enough('hc1cnona1.  1',bc,ibc)
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
      do isb=1,nsymb                                                    6d10s22
       jsb=multh(ipuse,nsbeta(isb))                                     6d10s22
       do k=0,nherec(isb)*nroot*nbdet(jsb)-1                            6d10s22
        bc(ig(isb)+k)=0d0                                               3d10s17
       end do                                                           3d10s17
       do k=0,nherect(isb)*nroot*nadet(jsb)-1                           6d10s22
        bc(igt(isb)+k)=0d0                                              3d10s17
       end do                                                           3d10s17
      end do
c
c     beta and alpha single and double differences
c
c
c          ig             ivec
c          bra            ket
c   alpha   beta     alpha  beta
c     1 ipuse*nsbeta   1   nsbeta,  alpha distributed, nherec
c          igt            ivect
c          bra            ket
c   beta   alpha     beta   alpha
c     1 ipuse*nsbeta   1   nsbeta,  beta distributed, nherect
c
      do isb=1,nsymb                                                    6d13s22
       kadd=ilc(isb)-1+ioffdeta(isb)                                    12d20s06
       jsbk=nsbeta(isb)                                                  6d13s22
       jsbb=multh(jsbk,ipuse)                                           6d13s22
       call dosingnona1(ibc(idatb1non(jsbb)),m12symnon(jsbb,1,2),       6d13s22
     $      bc(ivec(isb)),bc(ig(isb)),nroot,nherec(isb),ih0,iborb,      6d13s22
     $      nbeta,i2e,iaorb,nalpha,kadd,ibc(ism),ibc(irelo),iacto,      6d13s22
     $       ioffdetb(jsbb),ipxder,bc,ibc,igoal)                              11d14s22
       call dodoubnona1(ibc(idatb2non(jsbb)),m12symnon(jsbb,2,2),             6d13s22
     $      bc(ivec(isb)),bc(ig(isb)),nroot,nherec(isb),i2e,iacto,      6d13s22
     $      ibc(ism),ibc(irelo),ipxder,bc,ibc,igoal)                          11d14s22
       kadd=ilct(isb)-1+ioffdetb(isb)                                    12d20s06
       call dosingnona1(ibc(idata1non(jsbb)),m12symnon(jsbb,1,1),       6d13s22
     $      bc(ivect(isb)),bc(igt(isb)),nroot,nherect(isb),ih0,iaorb,   6d13s22
     $      nalpha,i2e,iborb,nbeta,kadd,ibc(ism),ibc(irelo),iacto,      6d13s22
     $       ioffdeta(jsbb),ipxder,bc,ibc,igoal)                              11d14s22
       call dodoubnona1(ibc(idata2non(jsbb)),m12symnon(jsbb,2,1),       6d13s22
     $      bc(ivect(isb)),bc(igt(isb)),nroot,nherect(isb),i2e,iacto,   6d13s22
     $      ibc(ism),ibc(irelo),ipxder,bc,ibc,igoal)                          11d14s22
      end do
c
c     alpha-beta terms
c
c          ig             ivec
c          bra            ket
c   alpha   beta     alpha  beta
c     1 ipuse*nsbeta   1   nsbeta,  alpha distributed, nherec
c          igt            ivect
c          bra            ket
c   beta   alpha     beta   alpha
c     1 ipuse*nsbeta   1   nsbeta,  beta distributed, nherect
c                         iveca
c                    alpha  beta
c                      1   nsbeta, beta full.
c     m12symnon singles have sym(bra)*sym(ket)=ipuse
c     these go with m12sym singles which have sym(bra)*sym(ket)=1.
c     we also have m1c singles with other symmetries whose products are
c     ipuse...
c
c    bdet adet
c     jsk=isk*isinfo(1,i)
c     jsb=isb*isinfo(1,i)*ipuse
c     if isb*isk=isbk, then isb*isinfo(1,i)*isk*isinfo(1,i)=isbk, so
c     jsb*jsk=isbk*ipuse.
c
      do isk=1,nsymb                                                    7d25s22
       jsk=nsbeta(isk)                                                  7d25s22
       do isb=1,nsymb                                                   7d25s22
        jsb=multh(nsbeta(isb),ipuse)                                    7d25s22
        if(isb.ge.isk)then                                               7d25s22
         iu1=1                                                           7d25s22
         iu2=2                                                           7d25s22
         ii=((isb*(isb-1))/2)+isk                                       7d25s22
        else                                                             7d25s22
         iu1=2                                                           7d25s22
         iu2=1                                                           7d25s22
         ii=((isk*(isk-1))/2)+isb                                       7d25s22
        end if                                                           7d25s22
        if(jsb.ge.jsk)then                                               7d25s22
         ju1=1                                                           7d25s22
         ju2=2                                                           7d25s22
         jj=((jsb*(jsb-1))/2)+jsk                                       7d25s22
        else                                                             7d25s22
         ju1=2                                                           7d25s22
         ju2=1                                                           7d25s22
         jj=((jsk*(jsk-1))/2)+jsb                                       7d25s22
        end if                                                           7d25s22
        if(isb.eq.isk)then                                              7d25s22
         if(min(m12sym(isb,1,1),m1c(jj,2)).gt.0)then                    7d25s22
c     ig is a distributed, beta full
          call doab1non(ibc(idata1(isb)),m12sym(isb,1,1),               7d25s22
     $         ibc(idatbc(jj)),m1c(jj,2),bc(ig(isb)),nherec(isb),       7d25s22
     $         ilc(isb),ihc(isb),bc(iveca(isk)),nadet(isk),i2e,ibc(ism),7d25s22
     $         ibc(irelo),iacto,iu1,iu2,ju1,ju2,ipxder,bc,ibc,igoal,    3d27s23
     $         nroot)                                                   3d27s23
          call doab1non(ibc(idata1(isb)),m12sym(isb,1,1),               7d25s22
     $         ibc(idatbc(jj)),m1c(jj,2),bc(ig(isb)),nherec(isb),       7d25s22
     $         ilc(isb),ihc(isb),bc(iveca(isk)),nadet(isk),i2e,ibc(ism),7d25s22
     $         ibc(irelo),iacto,iu2,iu1,ju1,ju2,ipxder,bc,ibc,igoal,    3d27s23
     $         nroot)                                                   3d27s23
         end if
        else if(jsb.eq.jsk)then                                         7d25s22
         if(min(m1c(ii,1),m12sym(jsb,1,2)).gt.0)then                    7d25s22
          call doab1non(ibc(idatac(ii)),m1c(ii,1),ibc(idatb1(jsb)),     7d25s22
     $        m12sym(jsb,1,2),bc(ig(isb)),nherec(isb),ilc(isb),         7d25s22
     $        ihc(isb),bc(iveca(isk)),nadet(isk),i2e,ibc(ism),          7d25s22
     $        ibc(irelo),iacto,iu1,iu2,ju1,ju2,ipxder,bc,ibc,igoal,     3d27s23
     $        nroot)                                                    3d27s23
          call doab1non(ibc(idatac(ii)),m1c(ii,1),ibc(idatb1(jsb)),     7d25s22
     $        m12sym(jsb,1,2),bc(ig(isb)),nherec(isb),ilc(isb),         7d25s22
     $        ihc(isb),bc(iveca(isk)),nadet(isk),i2e,ibc(ism),          7d25s22
     $        ibc(irelo),iacto,iu1,iu2,ju2,ju1,ipxder,bc,ibc,igoal,     3d27s23
     $        nroot)                                                    3d27s23
         end if
        else                                                            7d25s22
         if(min(m1c(ii,1),m1c(jj,2)).gt.0)then                          7d25s22
          call doab1non(ibc(idatac(ii)),m1c(ii,1),ibc(idatbc(jj)),      7d25s22
     $        m1c(jj,2),bc(ig(isb)),nherec(isb),ilc(isb),ihc(isb),      7d25s22
     $        bc(iveca(isk)),nadet(isk),i2e,ibc(ism),ibc(irelo),iacto,  7d25s22
     $        iu1,iu2,ju1,ju2,ipxder,bc,ibc,igoal,nroot)                4d19s23
         end if                                                         7d25s22
        end if                                                          7d25s22
       end do                                                           7d25s22
      end do                                                            7d25s22
      ibcoff=ism                                                        3d10s17
      return
      end
