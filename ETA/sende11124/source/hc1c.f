c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hc1c(ivec,ivect,iveca,ig,igt,ilc,ihc,ilct,ihct,        3d10s17
     $                nroot,hdig,ih0,nherec,                            3d10s17
     $                nherect,iaorb,iborb,nalpha,nbeta,i2e,             3d10s17
     $                m12sym,norb,nsymb,idata1,idata2,idatb1,idatb2,    3d10s17
     $                nsbeta,m1c,idatac,idatbc,mst,mnz,nz,lbail,bc,ibc) 11d10s22
      implicit real*8 (a-h,o-z)
      integer*2 idwstest2(4)                                            12d13s06
      integer*8 idwstest8                                               12d13s06
      logical lbail                                                     5d10s18
      equivalence (idwstest8,idwstest2)                                 12d13s06
      external second                                                   4d9s18
c
c     parallel formation of h*c
c
      dimension ig(8),igt(8),ivec(8),ivect(8),hdig(1),ih0(8),           3d10s17
     $     i2e(1),nherec(8),nsbeta(8),idata1(8),idata2(8),              3d10s17
     $     idatb1(8),idatb2(8),ilct(8),ihct(8),nherect(8),              3d10s17
     $     ilc(8),ihc(8),m12sym(8,2,2),iveca(8),idatac(36),idatbc(36),  3d15s17
     $     m1c(36,2),mst(8,8,2),mnz(3,64,36,2),nz(36,2)                 4d30s18
      include "common.cas"
      include "common.store"
      data icall/0/
      save
      ism=ibcoff                                                        8d29s06
      irelo=ism+norb                                                    8d29s06
      ibcoff=irelo+norb                                                 8d29s06
      call enough('hc1c.  1',bc,ibc)
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
c     diagononals
c
      if(nalpha.eq.nbeta.and.nsbeta(1).eq.-1)then                        12d13s06
       write(6,*)('zeroing g rather than diag part ')
       do isb=1,nsymb                                                   8d29s06
        numbx=nbdet(nsbeta(isb))                                        8d29s06
        if(nherec(isb).gt.0)then                                        12d14s06
         do k=0,nherec(isb)*nroot*numbx-1                                3d10s17
          bc(ig(isb)+k)=0d0                                              3d10s17
         end do                                                          3d10s17
        end if                                                          12d14s06
       end do
      else
       joff=1                                                           3d14s17
       ncount=0
       do isb=1,nsymb                                                   8d29s06
        numbx=nbdet(nsbeta(isb))                                        8d29s06
        jcvec=ivec(isb)                                                 3d10s17
        jg=ig(isb)                                                      3d14s17
        do k=0,numbx-1                                                  3d14s17
         do i=0,nherec(isb)-1                                           3d14s17
          ioff=joff+k+numbx*i                                           3d14s17
          ncount=ncount+1
          do j=0,nroot-1                                                3d14s17
           bc(jg)=bc(jcvec)*hdig(ioff)                                  3d14s17
           jcvec=jcvec+1                                                3d14s17
           jg=jg+1                                                      3d14s17
          end do
         end do
        end do
        joff=joff+numbx*nherec(isb)                                     8d29s06
       end do                                                           8d29s06
      end if
      do isb=1,nsymb
       do k=0,nroot*nherect(isb)*nadet(nsbeta(isb))-1                   3d10s17
        bc(igt(isb)+k)=0d0                                              3d10s17
       end do                                                           3d10s17
      end do                                                            3d10s17
c
c     beta and alpha single and double differences
c
      do isb=1,nsymb                                                    8d29s06
       kadd=ilc(isb)-1+ioffdeta(isb)                                    12d20s06
       jsb=nsbeta(isb)                                                  3d14s17
       numbx=nbdet(jsb)                                                 3d14s17
       call dosing(ibc(idatb1(jsb)),m12sym(jsb,1,2),bc(ivec(isb)),      3d10s17
     $     bc(ig(isb)),numbx,nroot,nherec(isb),ih0,iborb,nbeta,i2e,     3d10s17
     $       iaorb,nalpha,kadd,ibc(ism),ibc(irelo),iacto,               12d20s06
     $       ioffdetb(jsb),bc,ibc)                                      11d14s22
       call dodoub(ibc(idatb2(jsb)),m12sym(jsb,2,2),bc(ivec(isb)),      3d10s17
     $      bc(ig(isb)),numbx,nroot,nherec(isb),i2e,iacto,ibc(ism),     3d10s17
     $      ibc(irelo),bc,ibc)                                          11d14s22
       kadd=ilct(isb)-1+ioffdetb(isb)                                    12d20s06
       numbx=nadet(nsbeta(isb))                                         3d10s17
       call dosing(ibc(idata1(jsb)),m12sym(jsb,1,1),bc(ivect(isb)),     3d10s17
     $   bc(igt(isb)),numbx,nroot,nherect(isb),ih0,iaorb,nalpha,i2e,    3d14s17
     $       iborb,nbeta,kadd,ibc(ism),ibc(irelo),iacto,                3d14s17
     $       ioffdeta(nsbeta(isb)),bc,ibc)                              11d14s22
       call dodoub(ibc(idata2(jsb)),m12sym(jsb,2,1),bc(ivect(isb)),     3d14s17
     $      bc(igt(isb)),numbx,nroot,nherect(isb),i2e,iacto,ibc(ism),   3d10s17
     $      ibc(irelo),bc,ibc)                                          11d14s22
      end do
c
c     alpha-beta terms
c
      do isb=1,nsymb                                                    3d10s17
       jsb=nsbeta(isb)                                                  3d14s17
       numbx=nbdet(jsb)                                                 3d14s17
       if(nroot.gt.1)then                                               4d17s18
        call doab1(ibc(idatb1(jsb)),m12sym(jsb,1,2),ibc(idata1(isb)),    3d15s17
     $      m12sym(isb,1,1),bc(iveca(isb)),bc(ig(isb)),nadet(isb),numbx,3d10s17
     $      nroot,nherec(isb),ilc(isb),ihc(isb),i2e,ibc(irelo),         3d15s17
     $      ibc(ism),iacto,mst(1,isb,1),mst(1,jsb,2),nsymb,bc,ibc)      11d14s22
       else                                                             4d17s18
        call doab11(ibc(idatb1(jsb)),m12sym(jsb,1,2),ibc(idata1(isb)),  4d17s18
     $      m12sym(isb,1,1),bc(iveca(isb)),bc(ig(isb)),nadet(isb),numbx,3d10s17
     $      nherec(isb),ilc(isb),ihc(isb),i2e,ibc(irelo),               4d17s18
     $      ibc(ism),iacto,mst(1,isb,1),mst(1,jsb,2),nsymb,isb,lbail,   11d14s22
     $       bc,ibc)                                                    11d14s22
       end if                                                           4d17s18
       do isk=1,isb-1                                                   3d15s17
        iada=((isb*(isb-1))/2)+isk                                      3d15s17
        jsk=nsbeta(isk)                                                 3d15s17
        ix=max(jsb,jsk)                                                 3d15s17
        in=min(jsb,jsk)                                                 3d15s17
        iadb=((ix*(ix-1))/2)+in                                         3d15s17
        if(min(m1c(iada,1),m1c(iadb,2)).gt.0)then                       9d12s24
         if(ix.eq.jsb)then
          iu1=1
          iu2=2
         else
          iu1=2
          iu2=1
         end if                                                         3d16s17
         if(nroot.gt.1)then                                             4d17s18
          call doab2(ibc(idatac(iada)),m1c(iada,1),ibc(idatbc(iadb)),    3d16s17
     $        m1c(iadb,2),nadet(isb),                                   3d16s17
     $        nadet(isk),nbdet(jsb),nbdet(jsk),nroot,bc(ig(isb)),       3d16s17
     $        nherec(isb),ilc(isb),ihc(isb),bc(ig(isk)),nherec(isk),    3d16s17
     $        ilc(isk),ihc(isk),bc(iveca(isb)),bc(iveca(isk)),          3d16s17
     $        bc(iveca(jsb)),bc(iveca(jsk)),i2e,ibc(irelo),ibc(ism),    3d16s17
     $        iacto,iu1,iu2,mnz(1,1,iada,1),mnz(1,1,iadb,2),            4d30s18
     $         nz(iada,1),nz(iadb,2),bc,ibc)                            11d14s22
         else                                                           4d17s18
          call doab21(ibc(idatac(iada)),m1c(iada,1),ibc(idatbc(iadb)),  4d17s18
     $        m1c(iadb,2),nadet(isb),                                   4d17s18
     $        nadet(isk),nbdet(jsb),nbdet(jsk),bc(ig(isb)),             4d17s18
     $        nherec(isb),ilc(isb),ihc(isb),bc(ig(isk)),nherec(isk),    4d17s18
     $        ilc(isk),ihc(isk),bc(iveca(isb)),bc(iveca(isk)),          4d17s18
     $        bc(iveca(jsb)),bc(iveca(jsk)),i2e,ibc(irelo),ibc(ism),    4d17s18
     $        iacto,iu1,iu2,mnz(1,1,iada,1),mnz(1,1,iadb,2),            4d30s18
     $         nz(iada,1),nz(iadb,2),bc,ibc,norb)                            11d14s22
         end if                                                         4d17s18
        end if                                                          3d15s17
       end do                                                           3d15s17
      end do                                                            3d10s17
      ibcoff=ism                                                        3d10s17
      return
      end
