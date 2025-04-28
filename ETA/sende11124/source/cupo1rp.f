c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cupo1rp(idif,ndim,vecb,ndetb,ncsfb,iasb,nopenb,         8d30s21
     $     veck,ndetk,ncsfk,iask,nopenk,iout,nnot,n1den,ntype1,bc,ibc)  11d14s22
      implicit real*8 (a-h,o-z)                                         7d9s19
      real*16 fl                                                        7d11s19
c
c
      logical lb4,ldebug,log1,log2                                      9d9s19
      parameter (id8x=5 000)
      integer*1 iasb(nopenb,*),iask(nopenk,*),idif(*),iorb(64,2,2),     11d7s19
     $     ipackc1(8),id8d(id8x),i2mat(64),ipack14(4)                   8d31s21
      integer*1 imap(64,2),idet(64,2),idetu(64,2),ideta(64,2),iodif(4,2)1d17s23
      integer*2 ipack2(4)                                               7d9s19
      integer*4 ipack48(2),ipack4                                       8d31s21
      parameter (idos=60,idx=6,idu=64)                                  11d30s19
      integer*8 ipower2(idos)                                           11d29s19
      integer*8 ipack,ipackc,ipack8                                     8d31s21
      equivalence (ipack,ipack2),(ipackc,ipackc1),(ipack14,ipack4),     8d31s21
     $     (ipack8,ipack48)                                             8d31s21
      dimension vecb(ndetb,*),veck(ndetk,*),inot(4),                    8d30s21
     $     iou(2,2),igofor(64,2),ngofor(2),idatad(64),idatao(64),       11d30s19
     $     iordk(idos,2),idlast(idos),inpart(idx,idu,2),nblk(2),        11d30s19
     $     ispec(2,idu,2),iwhere(2,idos),ioffc(2,64)                    12d1s19
      include "common.store"                                            7d9s19
      include "common.print"                                            2d14s20
      COMMON/FACT16/FL(922),NCALL                                       5d27s19
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      save
      data icall/0/
      icall=icall+1                                                     7d10s19
      ntype1=0                                                          9d2s21
      n1den=0                                                           9d2s21
      nmat=0                                                            9d2s21
      if(icall.eq.1)then
       mxnmidm=0                                                        12d1s19
       ipower2(1)=2                                                      5d28s19
       do i=2,idos                                                       5d28s19
        im=i-1                                                           5d28s19
        ipower2(i)=ipower2(im)*2                                         5d28s19
       end do                                                            5d28s19
      end if                                                            11d29s19
      ibcoffo=ibcoff                                                    7d9s19
       ldebug=.false.                                                   2d14s20
       nmat=0                                                           9d1s21
      idelta=0                                                          9d6s21
      if(ndetk.gt.0)then                                                9d6s21
       do j=1,nopenk                                                    9d6s21
        idelta=idelta+iask(j,1)                                         9d6s21
        if(ldebug)write(6,*)('for ket, add '),iask(j,1),('to delta '),
     $       idelta
       end do                                                           9d6s21
      end if                                                            9d6s21
      if(ndetb.gt.0)then                                                9d6s21
       do j=1,nopenb                                                    9d6s21
        idelta=idelta-iasb(j,1)                                         9d6s21
        if(ldebug)write(6,*)('for bra, subtract '),iasb(j,1),
     $       ('from delta '),
     $       idelta
       end do                                                           9d6s21
      end if                                                            9d6s21
      if(ldebug)then                                                    12d14s19
       write(6,*)('cupo1r for call '),icall
       write(6,*)('iasb: ')
       do i=1,ndetb
        write(6,*)i,(iasb(j,i),j=1,nopenb)
       end do
       write(6,*)('vectors')
       call prntm2(vecb,ndetb,ncsfb,ndetb)                              8d30s21
       write(6,*)('iask: ')
       do i=1,ndetk
        write(6,*)i,(iask(j,i),j=1,nopenk)
       end do
       write(6,*)('vectors')
       call prntm2(veck,ndetk,ncsfk,ndetk)                              8d30s21
      end if
      nnot=0                                                            7d9s19
      ndoub=0                                                           8d5s19
      ntop=0                                                            4d6s20
      if(ldebug)write(6,*)('ndim in cupo2 '),ndim
      do i=1,ndim                                                       7d10s19
       if(idif(i).gt.0)then                                             7d9s19
        if(idif(i).ne.3)then                                            7d9s19
         nnot=nnot+1                                                    7d9s19
         if(idif(i).eq.4.or.idif(i).eq.6)ndoub=ndoub+1                  8d5s19
         inot(min(4,nnot))=i                                            7d9s19
        end if                                                          7d9s19
        ntop=i                                                          7d9s19
       end if                                                           7d9s19
      end do                                                            7d9s19
      if(ntop.eq.0)then                                                 4d6s20
       iout=ibcoff                                                      4d6s20
       nnot=1                                                           4d6s20
       return                                                           4d6s20
      end if                                                            4d6s20
      if(ldebug)write(6,*)('dif code: '),(idif(i),i=1,ntop)             12d14s19
c
c     new version to use dets ordered as in iasb,k
c
      npb=0                                                             9d9s19
      ndb=0                                                             9d9s19
      npk=0                                                             9d9s19
      ndk=0                                                             9d9s19
      n3s=0                                                             8d31s21
      nb=1                                                              8d31s21
      nk=1                                                              8d31s21
      do i=1,ndim                                                       9d9s19
       if(idif(i).eq.1.or.idif(i).eq.7)then                             9d9s19
        npb=npb+1                                                       9d9s19
        ndb=ndb+1                                                       9d9s19
        imap(npb,1)=ndb                                                 9d9s19
        idet(ndb,1)=i                                                   9d9s19
        ideta(ndb,1)=i                                                  9d9s19
        if(idif(i).eq.1)then                                            8d31s21
         iou(nb,1)=i                                                    8d31s21
         nb=min(2,nb+1)                                                 8d31s21
        else                                                            8d31s21
         iou(nk,2)=i                                                    8d31s21
         nk=min(2,nk+1)                                                 8d31s21
        end if                                                          8d31s21
       else if(idif(i).eq.2.or.idif(i).eq.5)then                        9d9s19
        npk=npk+1                                                       9d9s19
        ndk=ndk+1                                                       9d9s19
        imap(npk,2)=ndk                                                 9d9s19
        idet(ndk,2)=i                                                   9d9s19
        ideta(ndk,2)=i                                                  9d9s19
        if(idif(i).eq.2)then                                            8d31s21
         iou(nk,2)=i                                                    8d31s21
         nk=min(2,nk+1)                                                 8d31s21
        else                                                            8d31s21
         iou(nb,1)=i                                                    8d31s21
         nb=min(2,nb+1)                                                 8d31s21
        end if                                                          8d31s21
       else if(idif(i).eq.3)then                                        9d9s19
        n3s=n3s+1                                                       8d31s21
        idatad(n3s)=i                                                   8d31s21
        npb=npb+1                                                       9d9s19
        ndb=ndb+1                                                       9d9s19
        imap(npb,1)=ndb                                                 9d9s19
        idet(ndb,1)=i                                                   9d9s19
        ideta(ndb,1)=i                                                  9d9s19
        npk=npk+1                                                       9d9s19
        ndk=ndk+1                                                       9d9s19
        imap(npk,2)=ndk                                                 9d9s19
        idet(ndk,2)=i                                                   9d9s19
        ideta(ndk,2)=i                                                  9d9s19
       end if                                                           9d9s19
       if(idif(i).eq.4.or.idif(i).eq.5)then                             9d9s19
        ndb=ndb+1                                                       9d9s19
        idet(ndb,1)=-i                                                  9d9s19
        ideta(ndb,1)=i                                                  9d9s19
        ndb=ndb+1                                                       9d9s19
        idet(ndb,1)=+i                                                  9d9s19
        ideta(ndb,1)=i                                                  9d9s19
        if(idif(i).eq.4)then                                            8d31s21
         iou(nb,1)=i                                                    8d31s21
         iou(nb+1,1)=i                                                  8d31s21
         nb=3
        end if                                                          8d31s21
       else if(idif(i).eq.6.or.idif(i).eq.7)then                        9d9s19
        ndk=ndk+1                                                       9d9s19
        idet(ndk,2)=-i                                                  9d9s19
        ideta(ndk,2)=i                                                  9d9s19
        ndk=ndk+1                                                       9d9s19
        idet(ndk,2)=+i                                                  9d9s19
        ideta(ndk,2)=i                                                  9d9s19
        if(idif(i).eq.6)then                                            8d31s21
         iou(nk,2)=i                                                    8d31s21
         iou(nk+1,2)=i                                                  8d31s21
         nk=3                                                           8d31s21
        end if                                                          8d31s21
       end if                                                           9d9s19
      end do                                                            9d9s19
      ndbp=ndb+1                                                        9d9s19
      if(ldebug)then                                                    9d9s19
       write(6,*)('bra det: '),(idet(i,1),i=1,ndb)                      9d9s19
       write(6,*)('map: '),(imap(i,1),i=1,npb)
       write(6,*)('ket det: '),(idet(i,2),i=1,ndk)                      9d9s19
       write(6,*)('map: '),(imap(i,2),i=1,npk)
      end if                                                            9d9s19
      if(ndb.ne.ndk)then                                                9d9s19
       write(6,*)('ndb = '),ndb,(' ne ndk = '),ndk                      9d9s19
       write(6,*)('for call '),icall
       call dws_sync                                                    9d9s19
       call dws_finalize                                                9d9s19
       stop 'cupo1r'                                                     9d9s19
      end if                                                            9d9s19
      if(npb.ne.nopenb)then
       write(6,*)('npb = '),npb,(' ne nopenb = '),nopenb
       write(6,*)('for call '),icall
       call dws_sync
       call dws_finalize
       stop 'cupo1r'
      end if
      if(npk.ne.nopenk)then
       write(6,*)('npk = '),npk,(' ne nopenk = '),nopenk
       write(6,*)('for call '),icall
       call dws_sync
       call dws_finalize
       stop 'cupo1r'
      end if
c
      nhit=0                                                            7d10s19
c
      if(ldebug)write(6,*)('what we have for nnot: '),nnot                        8d31s21
      if(nnot.le.2)then                                                 7d9s19
       iout=ibcoff                                                      7d9s19
       idelta=iabs(idelta)                                              9d9s21
       n1x=1                                                            8d31s21
       if(nnot.eq.0)then
        if(idelta.eq.0)then                                             8d31s21
         if(ldebug)write(6,*)('nnot=0 and idelta=0 ')
         nnot=1                                                         9d1s21
        else if(idelta.eq.2)then                                        8d31s21
         if(ldebug)write(6,*)('idelta is 2, so reset nnot to 2')                  8d31s21
         nnot=2                                                         8d31s21
         n1x=n3s                                                        8d31s21
        else                                                            8d31s21
         if(ldebug)write(6,*)('idelta is gt 2, so reset nnot to 4')               8d31s21
         nnot=4                                                         8d31s21
        end if                                                          8d31s21
       else if(nnot.eq.2)then                                           8d31s21
        if(ldebug)write(6,*)('nb,nk = '),nb,nk
        if(min(nb,nk).gt.2)then                                         8d31s21
         if(ldebug)write(6,*)('turn this into double ')                           8d31s21
         nnot=4                                                         8d31s21
        end if                                                          8d31s21
       end if                                                           8d31s21
       if(nnot.gt.2)then                                                9d10s21
        nnot=0                                                          9d9s21
        return                                                          9d9s21
       end if                                                           9d9s21
       ntype1=0                                                         8d31s21
       n1den=0                                                          9d1s21
       if(nnot.eq.2)then                                                7d11s19
c     (bk)+(bk|cc)-(bc|ck)
        if(ldebug)write(6,*)('bk: '),iou(1,1),iou(1,2)                            8d31s21
        if(idelta.eq.0)then                                             8d31s21
         n1den=2*n1x                                                    8d31s21
        else                                                            8d31s21
         n1den=1*n1x                                                    8d31s21
        end if                                                          8d31s21
        nmat=n1den                                                      9d9s21
        itype=ibcoff                                                    8d31s21
        imat=itype+nmat                                                 8d31s21
        ibcoff=imat+nmat*ncsfb*ncsfk                                    8d31s21
        if(ldebug)write(6,*)('3s: '),(idatad(i),i=1,n3s)
       else if(nnot.eq.1)then                                           9d1s21
c     (cc|CC)-(cC|Cc)
        m1den=0
        itype=ibcoff                                                    7d11s19
        n1den=2*n3s                                                     9d1s21
        ntype1=0                                                        9d1s21
        nmat=n1den                                                      9d9s21
        imat=itype+nmat                                                 7d11s19
        ibcoff=imat+ncsfb*ncsfk*nmat                                    8d31s21
       end if                                                           7d9s19
       if(ldebug)write(6,*)('n3s, n1den,nmat '),n3s,n1den,nmat
       call enough('cupo1rp.  1',bc,ibc)
       do i=0,nmat-1                                                    8d31s21
        ibc(itype+i)=0                                                  8d31s21
       end do                                                           8d31s21
       do i=0,ncsfb*ncsfk*nmat-1                                        8d31s21
        bc(imat+i)=0d0                                                  12d1s19
       end do                                                           12d1s19
       ibtmpa=ibcoff                                                    7d9s19
       ibcoff=ibtmpa+nmat*max(1,ncsfb)                                  8d31s21
       call enough('cupo1rp.  2',bc,ibc)
   10    format(i6,5x,32i2)                                             7d9s19
   11   format(11x,32i2)                                                6d18s19
c
       do ikdet=1,max(1,ndetk)                                          11d4s19
        do i=0,nmat*max(ncsfb,1)-1                                      8d31s21
         ibc(ibtmpa+i)=0                                                 7d9s19
        end do                                                           7d9s19
c
c     new code
c
        do i=1,ndk                                                      9d9s19
         idetu(i,2)=idet(i,2)                                           9d9s19
        end do                                                          9d9s19
        do i=1,npk                                                      9d9s19
         idetu(imap(i,2),2)=idetu(imap(i,2),2)*iask(i,ikdet)            9d9s19
        end do                                                          9d9s19
c
c
c     code = 1 means o bra
c     code = 2 means o ket
c     code = 3 means o bra and ket
c     code = 4 means c bra
c     code = 5 means c bra o ket
c     code = 6 means c ket                                              6d24s19
c     code = 7 means c ket o bra
c
        do ibdet=1,max(ndetb,1)                                         11d3s19
         if(ldebug)write(6,*)('for bra det '),ibdet,('ket det '),ikdet
         do i=1,ndb                                                      9d9s19
          idetu(i,1)=idet(i,1)                                           9d9s19
         end do                                                          9d9s19
         if(ndetb.gt.0)then                                             11d3s19
          do i=1,npb                                                      9d9s19
           idetu(imap(i,1),1)=idetu(imap(i,1),1)*iasb(i,ibdet)            9d9s19
          end do                                                          9d9s19
         end if                                                         11d3s19
         iqsuma=-1                                                      9d9s19
         iqsumb=-1                                                      9d9s19
         jb=1                                                            9d9s19
         jk=1                                                            9d9s19
         idifo=0                                                         9d9s19
         idifb=0                                                        9d9s19
         idifk=0                                                        9d9s19
  256    continue                                                        9d9s19
          if(idetu(jb,1).eq.idetu(jk,2))then
           jb=jb+1
           jk=jk+1
          else if(ideta(jb,1).eq.ideta(jk,2))then                       9d9s19
           if(jb.lt.ndb.and.ideta(jb,1).eq.ideta(jb+1,1))then           9d9s19
            idifo=idifo+1                                               9d9s19
            if(idifo.gt.4)go to 257                                     9d9s19
            idifb=idifb+1                                                9d9s19
            iodif(idifb,1)=jb                                            9d9s19
            jb=jb+1                                                       9d9s19
           else
            idifo=idifo+1                                               9d9s19
            if(idifo.gt.4)go to 257                                     9d9s19
            idifk=idifk+1                                                9d9s19
            iodif(idifk,2)=jk                                            9d9s19
            jk=jk+1                                                       9d9s19
           end if                                                       9d9s19
          else if(ideta(jb,1).gt.ideta(jk,2))then                        9d9s19
           idifo=idifo+1                                                 9d9s19
           if(idifo.gt.4)go to 257                                       9d9s19
           idifk=idifk+1                                                9d9s19
           iodif(idifk,2)=jk                                            9d9s19
           jk=jk+1                                                       9d9s19
          else if(ideta(jb,1).lt.ideta(jk,2))then                        9d9s19
           idifo=idifo+1                                                 9d9s19
           if(idifo.gt.4)go to 257                                       9d9s19
           idifb=idifb+1                                                9d9s19
           iodif(idifb,1)=jb                                            9d9s19
           jb=jb+1                                                       9d9s19
          end if                                                         9d9s19
         if(min(jb,jk).le.ndb.and.max(jb,jk).le.ndb)go to 256
         do itry=1,ndbp-min(jb,jk)                                      9d9s19
          if(jb.le.ndb)then                                              9d9s19
           idifo=idifo+1                                                 9d9s19
           if(idifo.gt.4)go to 257                                       9d9s19
           idifb=idifb+1                                                 9d9s19
           iodif(idifb,1)=jb                                             9d9s19
           jb=jb+1
          end if                                                         9d9s19
          if(jk.le.ndk)then                                              9d9s19
           idifo=idifo+1                                                 9d9s19
           if(idifo.gt.4)go to 257                                       9d9s19
           idifk=idifk+1                                                 9d9s19
           iodif(idifk,2)=jk                                             9d9s19
           jk=jk+1
          end if                                                         9d9s19
         end do
         if(ldebug)then
          write(6,*)('at bottom with idifo = '),idifo                     9d9s19
          write(6,*)(idetu(jb,1),jb=1,ndb)                               9d9s19
          write(6,*)(idetu(jk,2),jk=1,ndk)                               9d9s19
          write(6,*)(iodif(jb,1),jb=1,idifb)
          write(6,*)(iodif(jk,2),jk=1,idifk)
         end if
         go to 258                                                       9d9s19
  257    continue                                                        9d9s19
         if(ldebug)write(6,*)('no coupling ')                                      9d9s19
         idifo=10                                                       8d31s21
  258    continue                                                        9d9s19
         if(idifo.ge.0.and.idifo.le.2)then                              9d9s21
          if(idifo.eq.0)then                                            8d31s21
           do ic=1,n3s                                                   8d31s21
            iduse=idetu(ic,1)                                           9d5s21
            ipackc1(1)=iduse                                            9d6s21
            ipackc1(2)=ipackc1(1)                                       8d31s21
            do i=3,8                                                     8d31s21
             ipackc1(i)=0                                                8d31s21
            end do                                                       8d31s21
            itype1=-1                                                    8d31s21
            do i=0,ntype1-1                                              8d31s21
             if(ipackc.eq.ibc(itype+i))itype1=i                          8d31s21
            end do                                                       8d31s21
            if(itype1.lt.0)then                                          8d31s21
             ibc(itype+ntype1)=ipackc                                    8d31s21
             itype1=ntype1                                               8d31s21
             ntype1=ntype1+1                                             8d31s21
             if(ntype1.gt.n1den)then                                     8d31s21
              write(6,*)('ntype1 exceeds n1den!!! '),ntype1,n1den,icall,
     $            ic
              do i=0,ntype1-1
               ipackc=ibc(itype+i)
               write(6,*)i+1,ipackc1
              end do
              stop 'cupo1r'
             end if                                                      8d31s21
            end if                                                       8d31s21
            jbtmpa=ibtmpa+ncsfb*itype1-1                                 8d31s21
            do ib=1+mynowprog,ncsfb,mynprocg                            3d23s22
             bc(jbtmpa+ib)=bc(jbtmpa+ib)+vecb(ibdet,ib)                 8d31s21
            end do                                                       8d31s21
           end do                                                       8d31s21
          else if(idifo.eq.2)then                                       8d31s21
           if(ldebug)write(6,*)('single difference '),
     $          (idetu(iodif(1,jx),jx),jx=1,2)
           ff=1d0                                                       8d31s21
           if(iodif(1,1).eq.iodif(1,2))then                             8d31s21
            if(ldebug)write(6,*)('there is no phase difference ')
           else                                                         8d31s21
            idd=iodif(1,1)-iodif(1,2)                                   8d31s21
            idd=iabs(idd)                                               8d31s21
            if(mod(idd,2).ne.0)ff=-1d0                                  8d31s21
            if(ldebug)write(6,*)('phase difference! '),ff
           end if
           ipackc1(1)=idetu(iodif(1,1),1)
           ipackc1(2)=idetu(iodif(1,2),2)
           ff1=ff                                                       9d5s21
           do i=3,8                                                     8d31s21
            ipackc1(i)=0                                                8d31s21
           end do                                                       8d31s21
           itype1=-1                                                    8d31s21
           do i=0,ntype1-1                                              8d31s21
            if(ipackc.eq.ibc(itype+i))itype1=i                          8d31s21
           end do                                                       8d31s21
           if(itype1.lt.0)then                                          8d31s21
            ibc(itype+ntype1)=ipackc                                    8d31s21
            itype1=ntype1                                               8d31s21
            ntype1=ntype1+1                                             8d31s21
            if(ntype1.gt.n1den)then                                     8d31s21
             write(6,*)('ntype1 exceeds n1den!!! '),ntype1,n1den,icall
              do i=0,ntype1-1
               ipackc=ibc(itype+i)
               write(6,*)i+1,ipackc1
              end do
             stop 'cupo1r'
            end if                                                      8d31s21
           end if                                                       8d31s21
           if(ldebug)write(6,*)('itype1 = '),itype1,ntype1
           jbtmpa=ibtmpa+ncsfb*itype1-1                                 8d31s21
           do ib=1+mynowprog,ncsfb,mynprocg                             3d23s22
            bc(jbtmpa+ib)=bc(jbtmpa+ib)+vecb(ibdet,ib)*ff1              9d5s21
           end do                                                       8d31s21
          end if                                                        8d31s21
         end if                                                         8d31s21
 1676    continue                                                       12d1s19
        end do                                                          7d9s19
        do i=0,nmat-1                                                   8d31s21
         jbtmpa=ibtmpa+ncsfb*i                                          8d31s21
         do ik=1,ncsfk                                                  9d1s21
          jmat=imat+ncsfb*(ik-1+ncsfk*i)                                9d1s21
          do ib=mynowprog,ncsfb-1,mynprocg                              3d23s22
           bc(jmat+ib)=bc(jmat+ib)+bc(jbtmpa+ib)*veck(ikdet,ik)         9d1s21
          end do                                                        9d1s21
         end do                                                         9d1s21
        end do                                                          9d1s21
       end do                                                           9d1s21
      else                                                              7d9s19
       nnot=0                                                           7d9s19
       nmat=0                                                           8d31s21
      end if                                                            7d9s19
      if(ntype1.eq.0)then                                                8d31s21
       nnot=0                                                           8d31s21
       nmat=0                                                           8d31s21
      end if                                                            8d31s21
      if(nnot.gt.0)then                                                 11d29s19
       ibcoff=ibtmpa                                                    8d31s21
      else                                                              11d29s19
       ibcoff=ibcoffo                                                   11d29s19
      end if                                                            11d29s19
      nnsum=ntype1*ncsfb*ncsfk                                          3d23s22
      call dws_gsumf(bc(imat),nnsum)                                    3d23s22
      if(nnot.gt.0.and.ldebug)then                                      9d1s21
       do i=0,ntype1-1                                                    8d31s21
        ipackc=ibc(itype+i)                                             8d31s21
        write(6,*)('for '),i,('integral type '),ipackc1                        8d31s21
        imatx=imat+ncsfb*ncsfk*i                                        8d31s21
        call prntm2(bc(imatx),ncsfb,ncsfk,ncsfb)                        8d31s21
       end do                                                           8d31s21
 1066  continue                                                         9d1s21
       write(6,*)('ibcoff vs. imatx '),ibcoff,imat+ncsfb*ncsfk*nmat     8d31s21
      end if                                                            8d31s21
      return                                                            7d9s19
      end                                                               7d9s19
