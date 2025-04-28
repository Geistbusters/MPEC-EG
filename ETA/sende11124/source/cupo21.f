c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cupo21(idif,ndim,ilb,ipb,vecb,ndetb,ncsfb,iasb,nopenb, 4d7s20
     $     ilk,ipk,veck,ndetk,ncsfk,iask,nopenk,iout,nnot,bc,ibc)       11d14s22
      implicit real*8 (a-h,o-z)                                         7d9s19
      real*16 fl                                                        7d11s19
c
c     this is like cupo2, except only return one particle density.      4d7s20
c
      logical lb4,ldebug,log1,log2                                      9d9s19
      parameter (id8x=1 500 000)
      integer*1 iasb(nopenb,*),iask(nopenk,*),idif(*),iorb(64,2,2),     11d7s19
     $     ipackc1(8),id8d(id8x),i2mat(64)                              12d1s19
      integer*1 imap(64,2),idet(64,2),idetu(64,2),ideta(64,2),iodif(2,2)9d9s19
      integer*2 ipack2(4)                                               7d9s19
      parameter (idos=60,idx=6,idu=64)                                  11d30s19
      integer*8 ipower2(idos)                                           11d29s19
      integer*8 ipack,ipackc                                            11d7s19
      equivalence (ipack,ipack2),(ipackc,ipackc1)                       11d7s19
      dimension vecb(*),veck(*),ilb(2,*),ipb(*),ilk(2,*),ipk(*),inot(4),11d27s19
     $     iou(2,2),igofor(64,2),ngofor(2),idatad(64),idatao(64),       11d30s19
     $     iordk(idos,2),idlast(idos),inpart(idx,idu,2),nblk(2),        11d30s19
     $     ispec(2,idu,2),iwhere(2,idos),ioffc(2,64),ialphabit(2),      5d7s20
     $     ibetabit(2),ialphabitu(2),ibetabitu(2),nab(2)                5d7s20
      include "common.store"                                            7d9s19
      include "common.print"                                            2d14s20
      COMMON/FACT16/FL(922),NCALL                                       5d27s19
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      save
      data icall/0/
      icall=icall+1                                                     7d10s19
      if(icall.eq.1)then
       mxnmidm=0                                                        12d1s19
       ipower2(1)=2                                                      5d28s19
       do i=2,idos                                                       5d28s19
        im=i-1                                                           5d28s19
        ipower2(i)=ipower2(im)*2                                         5d28s19
       end do                                                            5d28s19
      end if                                                            11d29s19
      ibcoffo=ibcoff                                                    7d9s19
      if(iprtr(15).eq.0)then                                            2d14s20
       ldebug=.false.                                                   2d14s20
      else                                                              2d14s20
       ldebug=.true.                                                    2d14s20
      end if                                                            2d14s20
      if(ldebug)then                                                    12d14s19
       write(6,*)('iasb: ')
       do i=1,ndetb
        write(6,*)i,(iasb(j,i),j=1,nopenb)
       end do
       write(6,*)('packed vectors ')
       do i=1,ndetb
        write(6,*)ilb(1,i),ilb(2,i)
        write(6,*)(ipb(k),k=ilb(1,i),ilb(2,i))
        write(6,*)(vecb(k),k=ilb(1,i),ilb(2,i))
       end do
      end if                                                            12d14s19
      if(ldebug)then                                                    12d14s19
       write(6,*)('iask: ')
       do i=1,ndetk
        write(6,*)i,(iask(j,i),j=1,nopenk)
       end do
       write(6,*)('packed vectors ')
       do i=1,ndetk
        write(6,*)ilk(1,i),ilk(2,i)
        write(6,*)(ipk(k),k=ilk(1,i),ilk(2,i))
        write(6,*)(veck(k),k=ilk(1,i),ilk(2,i))
       end do
      end if
      nnot=0                                                            7d9s19
      ndoub=0                                                           8d5s19
      ntop=0                                                            4d6s20
      if(ldebug)write(6,*)('ndim in cupo21 '),ndim
c     code = 1 means o bra
c     code = 2 means o ket
c     code = 3 means o bra and ket
c     code = 4 means c bra
c     code = 5 means c bra o ket
c     code = 6 means c ket                                              6d24s19
c     code = 7 means c ket o bra
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
      if(ntop.eq.0.or.nnot.eq.0)then                                    4d7s20
       iout=ibcoff                                                      4d6s20
       nnot=1                                                           4d6s20
       return                                                           4d6s20
      end if                                                            4d6s20
      if(ldebug)write(6,*)('dif code: '),(idif(i),i=1,ntop)             12d14s19
      nnot0=nnot                                                        11d30s19
c
c     new version to use dets ordered as in iasb,k
c
      npb=0                                                             9d9s19
      ndb=0                                                             9d9s19
      npk=0                                                             9d9s19
      ndk=0                                                             9d9s19
      ialphabit(1)=0                                                    5d7s20
      ialphabit(2)=0                                                    5d7s20
      ibetabit(1)=0                                                     5d7s20
      ibetabit(2)=0                                                     5d7s20
      do i=1,ndim                                                       9d9s19
       if(idif(i).eq.1.or.idif(i).eq.7)then                             9d9s19
        npb=npb+1                                                       9d9s19
        ndb=ndb+1                                                       9d9s19
        imap(npb,1)=ndb                                                 9d9s19
        idet(ndb,1)=i                                                   9d9s19
        ideta(ndb,1)=i                                                  9d9s19
       else if(idif(i).eq.2.or.idif(i).eq.5)then                        9d9s19
        npk=npk+1                                                       9d9s19
        ndk=ndk+1                                                       9d9s19
        imap(npk,2)=ndk                                                 9d9s19
        idet(ndk,2)=i                                                   9d9s19
        ideta(ndk,2)=i                                                  9d9s19
       else if(idif(i).eq.3)then                                        9d9s19
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
        ialphabit(1)=ibset(ialphabit(1),i)
        ibetabit(1)=ibset(ibetabit(1),i)
        ndb=ndb+1                                                       9d9s19
        idet(ndb,1)=-i                                                  9d9s19
        ideta(ndb,1)=i                                                  9d9s19
        ndb=ndb+1                                                       9d9s19
        idet(ndb,1)=+i                                                  9d9s19
        ideta(ndb,1)=i                                                  9d9s19
       else if(idif(i).eq.6.or.idif(i).eq.7)then                        9d9s19
        ialphabit(2)=ibset(ialphabit(2),i)
        ibetabit(2)=ibset(ibetabit(2),i)
        ndk=ndk+1                                                       9d9s19
        idet(ndk,2)=-i                                                  9d9s19
        ideta(ndk,2)=i                                                  9d9s19
        ndk=ndk+1                                                       9d9s19
        idet(ndk,2)=+i                                                  9d9s19
        ideta(ndk,2)=i                                                  9d9s19
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
       stop 'cupo21'                                                     9d9s19
      end if                                                            9d9s19
      if(npb.ne.nopenb)then
       write(6,*)('npb = '),npb,(' ne nopenb = '),nopenb
       write(6,*)('for call '),icall
       call dws_sync
       call dws_finalize
       stop 'cupo21'
      end if
      if(npk.ne.nopenk)then
       write(6,*)('npk = '),npk,(' ne nopenk = '),nopenk
       write(6,*)('for call '),icall
       call dws_sync
       call dws_finalize
       stop 'cupo21'
      end if
c
      nhit=0                                                            7d10s19
      if(nnot.eq.2)then                                                 4d10s20
       iout=ibcoff                                                      7d9s19
       imat=ibcoff                                                      4d7s20
       nmat=1                                                           4d7s20
       itmp1=imat+nmat*ncsfb*ncsfk                                      4d7s20
       ntmp1=nmat*ncsfb*max(1,ndetk)                                    4d7s20
       ibcoff=itmp1+ntmp1                                               7d9s19
       call enough('cupo21.  1',bc,ibc)
       do i=0,ncsfb*ncsfk*nmat-1                                        12d1s19
        bc(imat+i)=0d0                                                  12d1s19
       end do                                                           12d1s19
       do i=0,ntmp1-1                                                   7d9s19
        bc(itmp1+i)=0d0                                                 7d9s19
       end do                                                           7d9s19
   10    format(i6,5x,32i2)                                             7d9s19
   11   format(11x,32i2)                                                6d18s19
c
c     lets sort dets to make it easy to suss out couplings
c     entries with non-3 dif codes will be most significant.
c
        if(min(ndetb,ndetk).gt.3)then                                   4d17s20
         if(ldebug)then                                                 12d1s19
          write(6,*)('sorting dets '),(inot(i),i=1,nnot)
          write(6,*)('ndb vs ndk: '),ndb,ndk
          write(6,*)('dif code: '),(idif(i),i=1,ntop)                       7d9s19
          write(6,*)('nnot = '),nnot                                       8d5s19
         nmidm=ntop-nnot0                                               11d30s19
         nleft=ndb-nmidm                                                11d30s19
         nleftp=nleft+1                                                 11d30s19
         nmid=nmidm+1                                                   11d30s19
          write(6,*)('nmidm, '),nmidm,ntop,nnot0,nleft
         end if                                                         12d1s19
         nmidm=ntop-nnot0                                               11d30s19
         nleft=ndb-nmidm                                                11d30s19
         nleftp=nleft+1                                                 11d30s19
         nmid=nmidm+1                                                   11d30s19
         if(nmidm.gt.mxnmidm)then
          mxnmidm=nmidm                                                 12d1s19
         end if                                                         12d1s19
         if(nleft.gt.idx)then                                           11d30s19
          stop 'nleft:idx'                                              11d30s19
         end if                                                         11d30s19
c
c     iwhere will point in id8d where the extra strings are
c     it is address by no. of beta electrons(+1)
c
         do i=1,nmid                                                    11d30s19
          iwhere(1,i)=0                                                 11d30s19
         end do                                                         11d30s19
         jd8d=1                                                         11d30s19
         nbphase=0                                                      12d1s19
         do ibk=1,2
          if(ibk.eq.1)then                                              11d30s19
           ndet=ndetb
          else                                                          11d30s19
           ndet=ndetk                                                   11d30s19
          end if                                                        11d30s19
          isorta=ibcoff                                                  11d29s19
          itmpa=isorta+ndet                                             11d29s19
          ibcoff=itmpa+ndet                                              11d29s19
          call enough('cupo21.  2',bc,ibc)
          if(ibk.eq.1)then
           nduse=ndb
           npuse=npb
           isortb=isorta                                                11d30s19
          else
           nduse=ndk
           npuse=npk                                                    11d30s19
           isortk=isorta                                                11d30s19
          end if
c
c     reorder so that "extras" occur before "principles".
c     keep track of phase as well. to phase, imagine starting at the
c     beginning of the det and moving extras to the end, one at a time.
c     if the first is an extra, the phase will be (-1)**(nduse-1).
c     same for the second, etc, until we hit a principle. then the
c     phase will be (-1)**(nduse-2). Then we move the principles. Each
c     move will pick up the phase factor (-1)**(nduse-1).
c
          jorda=1                                                        11d30s19
          nb=0                                                          12d1s19
          nphase=nduse-1                                                12d1s19
          do i=1,nduse                                                     11d29s19
           itry=idet(i,ibk)                                                11d29s19
           itry=iabs(itry)                                               11d29s19
           do j=1,nnot                                                   11d29s19
            if(itry.eq.inot(j))then                                      11d29s19
             nphase=nphase-1                                            12d1s19
             go to 1216                                                  11d29s19
            end if                                                       11d29s19
           end do                                                        11d29s19
           iordk(jorda,ibk)=i                                                11d30s19
           jorda=jorda+1                                                 11d29s19
           nb=nb+nphase                                                 12d1s19
 1216      continue                                                      11d29s19
          end do                                                         11d29s19
          nphase=nduse-1                                                12d1s19
          do i=1,nduse                                                     11d29s19
           itry=idet(i,ibk)                                                11d29s19
           itry=iabs(itry)                                               11d29s19
           do j=1,nnot                                                   11d29s19
            if(itry.eq.inot(j))then                                      11d29s19
             iordk(jorda,ibk)=i                                              11d30s19
             jorda=jorda+1                                               11d29s19
             nb=nb+nphase                                               12d1s19
             go to 1226                                                  11d29s19
            end if                                                       11d29s19
           end do                                                        11d29s19
 1226      continue                                                      11d29s19
          end do                                                         11d29s19
          nbphase=nbphase+mod(nb,2)                                     12d1s19
          do jdet=1,ndet                                                11d30s19
           do i=1,nduse                                                 11d30s19
            idetu(i,ibk)=idet(i,ibk)                                    11d30s19
           end do                                                          9d9s19
           if(ibk.eq.1)then                                             11d30s19
            do i=1,npuse                                                 11d30s19
             idetu(imap(i,ibk),ibk)=idetu(imap(i,ibk),ibk)*iasb(i,jdet) 11d30s19
            end do                                                          9d9s19
           else                                                         11d30s19
            do i=1,npuse                                                 11d30s19
             idetu(imap(i,ibk),ibk)=idetu(imap(i,ibk),ibk)*iask(i,jdet)            9d9s19
            end do                                                          9d9s19
           end if                                                       11d30s19
           jtmpa=itmpa+jdet-1                                           11d29s19
           ibc(jtmpa)=0                                                  11d29s19
           do i=1,nduse                                                 11d30s19
            if(idetu(iordk(i,ibk),ibk).lt.0)                            11d30s19
     $           ibc(jtmpa)=ibc(jtmpa)+ipower2(i)                       11d30s19
           end do                                                        11d29s19
          end do                                                         11d29s19
          call idsortdws(ibc(itmpa),ibc(isorta),ndet)                   1d18s23
          ibcoff=itmpa                                                  11d30s19
          jsorta=isorta-1                                                11d29s19
          nblk(ibk)=1                                                        11d30s19
          do jdeto=1,ndet                                               11d30s19
           jdet=ibc(jsorta+jdeto)                                       11d30s19
           do i=1,nduse                                                 11d30s19
            idetu(i,ibk)=idet(i,ibk)                                    11d30s19
           end do                                                          9d9s19
           if(ibk.eq.1)then                                             11d30s19
            do i=1,npuse                                                11d30s19
             idetu(imap(i,ibk),ibk)=idetu(imap(i,ibk),ibk)*iasb(i,jdet) 11d30s19
            end do                                                          9d9s19
           else                                                         11d30s19
            do i=1,npuse                                                11d30s19
             idetu(imap(i,ibk),ibk)=idetu(imap(i,ibk),ibk)*iask(i,jdet) 11d30s19
            end do                                                          9d9s19
           end if                                                       11d30s19
           if(jdeto.eq.1)then                                           11d30s19
            ii=1
            do i=nmid,nduse                                             11d30s19
             idlast(ii)=idetu(iordk(i,ibk),ibk)                         11d30s19
             inpart(ii,nblk(ibk),ibk)=idlast(ii)                             11d30s19
             ii=ii+1
            end do
            ispec(2,nblk(ibk),ibk)=jdeto                                  11d30s19
            nneg=0                                                      11d30s19
            do i=1,nmidm                                                11d30s19
             if(idetu(iordk(i,ibk),ibk).lt.0)nneg=nneg+1                11d30s19
            end do                                                      11d30s19
            ispec(1,nblk(ibk),ibk)=nneg                                 11d30s19
            nnegp=nneg+1                                                11d30s19
            if(iwhere(1,nnegp).eq.0)then                                11d30s19
             log1=.true.                                                11d30s19
             iwhere(1,nnegp)=jd8d                                       11d30s19
             iwhere(2,nnegp)=1                                          11d30s19
             if(jd8d+nmid.gt.id8x)then
              write(6,*)('id8x crashing in cupo21 ')
              write(6,*)('call no. '),icall
              write(6,*)jd8d,nmid,id8x
              stop 'id8x'
             end if
             do i=1,nmidm                                               11d30s19
              id8d(jd8d)=idetu(iordk(i,ibk),ibk)                        11d30s19
              jd8d=jd8d+1                                               11d30s19
             end do                                                     11d30s19
            else                                                        11d30s19
             log1=.false.                                               11d30s19
            end if                                                      11d30s19
            nblk(ibk)=nblk(ibk)+1                                                 11d30s19
           else                                                          11d30s19
            iddf=0
            ii=1
            do i=nmid,nduse                                             11d30s19
             itry=idlast(ii)-idetu(iordk(i,ibk),ibk)                    11d30s19
             iddf=iddf+iabs(itry)
             ii=ii+1
            end do
            if(iddf.ne.0)then
             ii=1
             if(nblk(ibk).gt.idu)stop 'nblk:idu'                             11d30s19
             do i=nmid,nduse                                            11d30s19
              idlast(ii)=idetu(iordk(i,ibk),ibk)                        11d30s19
              inpart(ii,nblk(ibk),ibk)=idlast(ii)                             11d30s19
              ii=ii+1
             end do
             ispec(2,nblk(ibk),ibk)=jdeto                               11d30s19
             nneg=0                                                      11d30s19
             do i=1,nmidm                                               11d30s19
              if(idetu(iordk(i,ibk),ibk).lt.0)nneg=nneg+1                11d30s19
             end do                                                      11d30s19
             ispec(1,nblk(ibk),ibk)=nneg                                 11d30s19
             nnegp=nneg+1                                               11d30s19
             if(iwhere(1,nnegp).eq.0)then                                11d30s19
              log1=.true.                                                11d30s19
              iwhere(1,nnegp)=jd8d                                       11d30s19
              iwhere(2,nnegp)=1                                          11d30s19
              if(jd8d+nmid.gt.id8x)then
               write(6,*)('id8x crashing in cupo21 2nd time')
               write(6,*)('call no. '),icall
               write(6,*)jd8d,nmid,id8x
               stop 'id8x'
              end if
              do i=1,nmidm                                              11d30s19
               id8d(jd8d)=idetu(iordk(i,ibk),ibk)                        11d30s19
               jd8d=jd8d+1                                               11d30s19
              end do                                                     11d30s19
             else                                                        11d30s19
              log1=.false.                                               11d30s19
             end if                                                      11d30s19
             nblk(ibk)=nblk(ibk)+1                                                11d30s19
            else if(log1)then                                           11d30s19
             iwhere(2,nnegp)=iwhere(2,nnegp)+1                          11d30s19
              if(jd8d+nmid.gt.id8x)then
               write(6,*)('id8x crashing in cupo21 third time')
               write(6,*)('call no. '),icall
               write(6,*)jd8d,nmid,id8x
               stop 'id8x'
              end if
             do i=1,nmidm                                               11d30s19
              id8d(jd8d)=idetu(iordk(i,ibk),ibk)                        11d30s19
              jd8d=jd8d+1                                               11d30s19
             end do                                                     11d30s19
            end if
           end if
          end do                                                         11d29s19
          nblk(ibk)=nblk(ibk)-1                                                    11d30s19
         end do                                                         11d30s19
         do ikb=1,nblk(2)                                               11d30s19
          do i=1,nleft                                                  11d30s19
           idetu(i,2)=inpart(i,ikb,2)                                   11d30s19
           ideta(i,2)=iabs(inpart(i,ikb,2))                             11d30s19
          end do                                                        11d30s19
          do ibb=1,nblk(1)
           if(iabs(ispec(1,ikb,2)-ispec(1,ibb,1)).gt.1)go to 1258       4d17s20
           do i=1,nleft                                                  11d30s19
            idetu(i,1)=inpart(i,ibb,1)                                   11d30s19
            ideta(i,1)=iabs(inpart(i,ibb,1))                             11d30s19
           end do                                                        11d30s19
           iqsuma=-1                                                      9d9s19
           iqsumb=-1                                                      9d9s19
           jb=1                                                            9d9s19
           jk=1                                                            9d9s19
           idifo=0                                                         9d9s19
           idifb=0                                                        9d9s19
           idifk=0                                                        9d9s19
 1256      continue                                                        9d9s19
           if(idetu(jb,1).eq.idetu(jk,2))then
            jb=jb+1
            jk=jk+1
           else if(ideta(jb,1).eq.ideta(jk,2))then                       9d9s19
            if(jb.lt.nleft.and.ideta(jb,1).eq.ideta(jb+1,1))then           9d9s19
             idifo=idifo+1                                               9d9s19
             if(idifo.gt.2)go to 1257                                     9d9s19
             idifb=idifb+1                                                9d9s19
             iodif(idifb,1)=jb                                            9d9s19
             jb=jb+1                                                       9d9s19
            else
             idifo=idifo+1                                               9d9s19
             if(idifo.gt.2)go to 1257                                     9d9s19
             idifk=idifk+1                                                9d9s19
             iodif(idifk,2)=jk                                            9d9s19
             jk=jk+1                                                       9d9s19
            end if                                                       9d9s19
           else if(ideta(jb,1).gt.ideta(jk,2))then                        9d9s19
            idifo=idifo+1                                                 9d9s19
            if(idifo.gt.2)go to 1257                                       9d9s19
            idifk=idifk+1                                                9d9s19
            iodif(idifk,2)=jk                                            9d9s19
            jk=jk+1                                                       9d9s19
           else if(ideta(jb,1).lt.ideta(jk,2))then                        9d9s19
            idifo=idifo+1                                                 9d9s19
            if(idifo.gt.2)go to 1257                                       9d9s19
            idifb=idifb+1                                                9d9s19
            iodif(idifb,1)=jb                                            9d9s19
            jb=jb+1                                                       9d9s19
           end if                                                         9d9s19
           if(min(jb,jk).le.nleft.and.max(jb,jk).le.nleft)go to 1256
           do itry=1,nleftp-min(jb,jk)                                  11d30s19
            if(jb.le.nleft)then                                              9d9s19
            idifo=idifo+1                                                 9d9s19
            if(idifo.gt.2)go to 1257                                       9d9s19
            idifb=idifb+1                                                 9d9s19
            iodif(idifb,1)=jb                                             9d9s19
            jb=jb+1
           end if                                                         9d9s19
           if(jk.le.nleft)then                                              9d9s19
            idifo=idifo+1                                                 9d9s19
            if(idifo.gt.2)go to 1257                                       9d9s19
            idifk=idifk+1                                                 9d9s19
            iodif(idifk,2)=jk                                             9d9s19
            jk=jk+1
           end if                                                         9d9s19
          end do
          mb=iodif(1,1)-iodif(1,2)                                      9d9s19
          nbp=iabs(mb)                                                  9d9s19
          if(idetu(iodif(1,1),1).gt.0)then
           iqsuma=2                                                     9d9s19
           iqsumb=0                                                     9d9s19
          else                                                          9d9s19
           iqsuma=0                                                     9d9s19
           iqsumb=2                                                     9d9s19
          end if                                                        9d9s19
          iou(1,1)=ideta(iodif(1,1),1)                                  9d9s19
          iou(1,2)=ideta(iodif(1,2),2)                                  9d9s19
          iqsum=iqsuma+iqsumb                                           11d30s19
          isuma=iqsuma                                                  12d1s19
          isumb=iqsumb                                                  12d1s19
          if(ispec(1,ikb,2).eq.ispec(1,ibb,1))then                      11d30s19
           if(isuma.eq.2)then                                           7d9s19
            iu=1                                                        7d9s19
            iup=+1                                                      9d9s19
           else                                                         7d9s19
            iu=2                                                        7d9s19
            iup=-1                                                      9d9s19
           end if                                                       7d9s19
           nb=nbp+nbphase                                               12d1s19
           if(mod(nb,2).eq.0)then                                       7d9s19
            ff=1d0                                                      7d9s19
           else                                                         7d9s19
            ff=-1d0                                                     7d9s19
           end if                                                       7d9s19
           nnegp=ispec(1,ikb,2)+1                                       12d1s19
           idetks=ispec(2,ikb,2)-1                                      12d1s19
           idetbs=ispec(2,ibb,1)-1                                      12d1s19
           do isk=1,iwhere(2,nnegp)                                     12d1s19
            ikdet=ibc(isortk+idetks)                                    12d1s19
            ibdet=ibc(isortb+idetbs)                                    12d1s19
            do l=ilk(1,ikdet),ilk(2,ikdet)                              11d29s19
             ktmp=imat-1+ncsfb*(ipk(l)-1)                               12d1s19
             gg=ff*veck(l)                                              11d29s19
             do j=ilb(1,ibdet),ilb(2,ibdet)                             11d29s19
              iad=ktmp+ipb(j)                                           11d29s19
              bc(iad)=bc(iad)+gg*vecb(j)                                11d29s19
             end do                                                     11d29s19
            end do                                                      11d29s19
            idetks=idetks+1                                             12d1s19
            idetbs=idetbs+1                                             12d1s19
           end do                                                       12d1s19
          end if                                                        11d30s19
          go to 1258                                                       9d9s19
 1257     continue                                                        9d9s19
 1258     continue                                                        9d9s19
          end do
         end do
          ibcoff=itmp1                                                  12d1s19
          return                                                        12d1s19
        end if                                                          11d29s19
c
       do ikdet=1,max(1,ndetk)                                          11d4s19
c
c     new code
c
        do i=1,ndk                                                      9d9s19
         idetu(i,2)=idet(i,2)                                           9d9s19
        end do                                                          9d9s19
        ialphabitu(2)=ialphabit(2)                                      5d7s20
        ibetabitu(2)=ibetabit(2)                                        5d7s20
        do i=1,npk                                                      9d9s19
         idetu(imap(i,2),2)=idetu(imap(i,2),2)*iask(i,ikdet)            9d9s19
         if(iask(i,ikdet).gt.0)then                                     5d7s20
          ialphabitu(2)=ibset(ialphabitu(2),ideta(imap(i,2),2))                          5d7s20
         else                                                           5d7s20
          ibetabitu(2)=ibset(ibetabitu(2),ideta(imap(i,2),2))                            5d7s20
         end if                                                         5d7s20
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
         ialphabitu(1)=ialphabit(1)                                     5d7s20
         ibetabitu(1)=ibetabit(1)                                       5d7s20
         if(ndetb.gt.0)then                                             11d3s19
          do i=1,npb                                                      9d9s19
           if(iasb(i,ibdet).gt.0)then
            ialphabitu(1)=ibset(ialphabitu(1),ideta(imap(i,1),1))                        5d7s20
           else                                                         5d7s20
            ibetabitu(1)=ibset(ibetabitu(1),ideta(imap(i,1),1))                          5d7s20
           end if                                                       5d7s20
          end do                                                          9d9s19
         end if                                                         11d3s19
         nbb=0                                                          5d7s20
         iatest=ieor(ialphabitu(1),ialphabitu(2))                       5d7s20
         iadiff=popcnt(iatest)                                          5d7s20
         if(iadiff.le.2)then                                            5d7s20
          ibtest=ieor(ibetabitu(1),ibetabitu(2))                         5d7s20
          ibdiff=popcnt(ibtest)                                          5d7s20
          if(iadiff+ibdiff.gt.2)then                                    5d7s20
           iadiff=-1                                                    5d7s20
           ibdiff=-1                                                    5d7s20
          else                                                          5d7s20
           if(iadiff.eq.2)then                                          5d7s20
            nbb=0                                                        5d7s20
            iq=0                                                        5d7s20
            do i=1,ndim                                                 5d7s20
             if(btest(iatest,i))then                                    5d7s20
              iq=iq+1                                                   5d7s20
              nab(iq)=i                                                 5d7s20
              if(iq.eq.2)go to 3000                                     5d7s20
             end if                                                     5d7s20
            end do                                                      5d7s20
 3000       continue                                                    5d7s20
            if(btest(ialphabitu(1),nab(1)))then                         5d7s20
             iu=1                                                       5d7s20
            else
             iu=2                                                       5d7s20
            end if                                                      5d7s20
            do i=nab(1)+1,nab(2)-1
             if(btest(ialphabitu(iu),i).or.btest(ibetabitu(iu),i))
     $            nbb=nbb+1
            end do
            if(btest(ialphabitu(1),nab(2)).and.
     $           btest(ibetabitu(1),nab(2)))nbb=nbb+1
            if(btest(ialphabitu(2),nab(2)).and.
     $           btest(ibetabitu(2),nab(2)))nbb=nbb-1
            nbb=iabs(nbb)
           else if(ibdiff.eq.2)then                                     5d7s20
            nbb=0                                                        5d7s20
            iq=0                                                        5d7s20
            do i=1,ndim                                                 5d7s20
             if(btest(ibtest,i))then                                    5d7s20
              iq=iq+1                                                   5d7s20
              nab(iq)=i                                                 5d7s20
              if(iq.eq.2)go to 3001                                     5d7s20
             end if                                                     5d7s20
            end do                                                      5d7s20
 3001       continue                                                    5d7s20
            if(btest(ibetabitu(1),nab(1)))then                          5d7s20
             iu=1                                                       5d7s20
            else                                                        5d7s20
             iu=2                                                       5d7s20
            end if                                                      5d7s20
            if(btest(ialphabitu(iu),nab(1)))nbb=nbb+1
            do i=nab(1)+1,nab(2)-1
             if(btest(ialphabitu(iu),i).or.btest(ibetabitu(iu),i))
     $            nbb=nbb+1
            end do
           end if                                                       5d7s20
          end if                                                        5d7s20
         else                                                           5d7s20
          iadiff=-1                                                     5d7s20
          ibdiff=-1                                                     5d7s20
         end if                                                         5d7s20
         isuma=iadiff
         isumb=ibdiff
         nbp=nbb
         isum=isuma+isumb                                               7d9s19
         if(isum.eq.2)then                                              4d10s20
          if(ldebug)write(6,*)('for bra '),ibdet,(' kdet '),ikdet
           nb=nbp                                                       9d9s19
           if(ldebug)write(6,*)('nb for phase '),nb
           if(mod(nb,2).eq.0)then                                       7d9s19
            ff=1d0                                                      7d9s19
           else                                                         7d9s19
            ff=-1d0                                                     7d9s19
           end if                                                       7d9s19
           if(ldebug)write(6,*)('overall phase '),ff,ndetb,ndetk                              7d9s19
           jtmp=itmp1-1+ncsfb*(ikdet-1)                                 7d9s19
           nhit=nhit+1                                                  7d10s19
           if(ndetb.gt.0)then                                           11d3s19
            if(ndetk.eq.0)then                                          11d29s19
             do j=ilb(1,ibdet),ilb(2,ibdet)                              11d29s19
              iad=jtmp+ipb(j)                                            11d29s19
              bc(iad)=bc(iad)+ff*vecb(j)                                 11d29s19
             end do                                                      11d29s19
            else                                                        11d29s19
             do l=ilk(1,ikdet),ilk(2,ikdet)                             11d29s19
              ktmp=itmp1-1+ncsfb*(ipk(l)-1)                             11d29s19
              gg=ff*veck(l)                                             11d29s19
              do j=ilb(1,ibdet),ilb(2,ibdet)                            11d29s19
               iad=ktmp+ipb(j)                                          11d29s19
               orig=bc(iad)
               bc(iad)=bc(iad)+gg*vecb(j)                               11d29s19
              end do                                                    11d29s19
             end do                                                     11d29s19
            end if                                                      11d29s19
           else                                                         11d3s19
            if(ndetk.eq.0)then                                          11d29s19
             bc(jtmp+1)=bc(jtmp+1)+ff                                    11d4s19
            else                                                        11d29s19
             do l=ilk(1,ikdet),ilk(2,ikdet)                             11d29s19
              ktmp=itmp1+ncsfb*(ipk(l)-1)                               11d29s19
              bc(ktmp)=bc(ktmp)+ff*veck(l)                              11d29s19
             end do                                                     11d29s19
            end if                                                      11d29s19
           end if                                                       11d3s19
          end if                                                        4d17s20
        end do                                                          7d9s19
       end do                                                           7d9s19
       jtmp=itmp1                                                       7d8s19
       jmat=imat                                                        7d8s19
       do i=1,nmat                                                      7d8s19
        do ik0=0,max(1,ncsfk)-1                                         11d29s19
         ltmp=jtmp+ncsfb*ik0                                            11d29s19
         lmat=jmat+ncsfb*ik0                                            11d29s19
         do ib0=0,max(1,ncsfb)-1                                        11d29s19
          bc(lmat+ib0)=bc(ltmp+ib0)                                     11d29s19
         end do                                                         11d29s19
        end do                                                          11d29s19
        jtmp=jtmp+ncsfb*max(1,ndetk)                                    11d4s19
        jmat=jmat+ncsfb*ncsfk                                           7d9s19
       end do                                                           7d8s19
       ibcoff=itmp1                                                     7d9s19
       if(nnot.eq.0)nnot=1                                              7d11s19
      else                                                              7d9s19
       nnot=0                                                           7d9s19
      end if                                                            7d9s19
      if(nhit.eq.0)nnot=0                                               7d10s19
      if(nnot.gt.0)then                                                 11d29s19
       ibcoff=jmat                                                      11d29s19
      else                                                              11d29s19
       ibcoff=ibcoffo                                                   11d29s19
      end if                                                            11d29s19
      return                                                            7d9s19
      end                                                               7d9s19
