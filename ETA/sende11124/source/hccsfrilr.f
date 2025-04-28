c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hccsfrilr(gl,ndl,vecr,ndr,ibasisl,ibasisr,ncsf,nfcnl,  3d23s21
     $     nfcnr,ih0a,i2e,nrootz,mdon,nec,multh,nkeep,ikeep,ibasisk,    3d24s21
     $     mdoo,iptrbl,iptrbr,ixw1,ixw2,iptkb,bc,ibc)                   11d10s22
      implicit real*8 (a-h,o-z)                                         7d11s19
      external second                                                   8d1s19
      integer*1 nab(2)                                                  3d23s21
      dimension vecr(ndr,nrootz),gl(ndl,nrootz),ih0a(*),i2e(*),         3d23s21
     $     ibasisl(3,*),ibasisr(3,*),ncsf(*),iptrbr(2,*),iptrbl(2,*),   3d23s21
     $     multh(8,8),igkv(8,8),iekv(8,8),ibasisk(3,*),iptkb(2,*),      3d24s21
     $     ikeep(*)                                                     3d24s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      include "common.hf"                                               3d31s20
      include "common.store"                                            7d11s19
      include "common.mrci"                                             6d19s19
      include "common.print"                                            1d13s20
      nstat=0                                                           4d17s20
      ibcoffo=ibcoff                                                    4d8s20
      ipass=0
      teacup=0d0                                                        8d1s19
      do i=1,nrootz                                                     7d11s19
       do j=1,ndl                                                       3d25s21
        gl(j,i)=0d0                                                     7d11s19
       end do                                                           7d11s19
      end do                                                            7d11s19
      innot1=0                                                          4d8s20
      krps=1                                                             4d8s20
      do kr=1,min(nfcnr,mynowprog)                                       10d20s20
       nclor=ibasisr(1,kr)                                               10d15s20
       nclorp=nclor+1                                                     7d11s19
       iargr=nclorp-mdon                                                  7d11s19
       krps=krps+ncsf(iargr)                                               4d8s20
      end do                                                            4d8s20
      ipss=1                                                            4d16s20
      ifs=1                                                             4d16s20
      jpss=1                                                            4d16s20
      jfs=1                                                             4d16s20
      do kr=1+mynowprog,nfcnr,mynprocg                                   10d15s20
       nclor=ibasisr(1,kr)                                               10d15s20
       nclorp=nclor+1                                                     7d11s19
       iargr=nclorp-mdon                                                  7d11s19
       nopenr=nec-2*nclor                                                 7d11s19
       kkcr=iptrbr(1,nclorp)+ibasisr(2,kr)-1                               10d15s20
       kkor=iptrbr(2,nclorp)+ibasisr(3,kr)-1                               10d15s20
       nrow=ncsf(iargr)*nrootz                                           4d7s20
       jps=jpss                                                         4d16s20
       jfsn=nfcnl+1                                                     3d23s21
       nhit=0                                                           3d24s21
       do jl=jfs,nfcnl                                                  3d23s21
        nclol=ibasisl(1,jl)                                              7d11s19
        if(nclol.gt.nclorp)go to 12                                      4d16s20
        nclopl=nclol+1                                                  6d11s19
        iargl=nclopl-mdon                                                4d7s20
        if(iabs(nclol-nclor).le.1)then                                  4d7s20
         if(jfsn.gt.jl)then                                             4d16s20
          jfsn=jl                                                       4d16s20
          jpss=jps                                                      4d16s20
         end if                                                         4d16s20
         nopenl=nec-2*nclol
         jjcl=iptrbl(1,nclopl)+ibasisl(2,jl)-1                              5d7s20
         jjol=iptrbl(2,nclopl)+ibasisl(3,jl)-1                              5d7s20
         call gandc(ibc(jjcl),ibc(jjol),ibc(kkcr),ibc(kkor),nopenl,
     $        nopenr,iargl,iargr,ncsf,norb,ixw1,ixw2,nnot,nab,iwpb,iwpk,
     $        ncsfmid,bc,ibc)                                           11d14s22
         if(nnot.eq.1)then                                              4d7s20
          nhit=nhit+1                                                   3d24s21
          do i=1,norb                                                   3d24s21
           if(btest(ibc(jjcl),i))then                                   3d23s21
            isi=ism(i)                                                  4d8s20
            igi=irel(i)-1                                               4d8s20
            ih0u=ih0a(isi)+igi*(irefo(isi)+1)                            4d8s20
            sum=2d0*bc(ih0u)                                             4d8s20
            do ir=1,nrootz                                               4d8s20
             do j=0,ncsf(iargl)-1                                         4d8s20
              gl(jps+j,ir)=gl(jps+j,ir)+sum*vecr(krps+j,ir)             3d24s21
             end do                                                      4d8s20
            end do                                                       4d8s20
           end if                                                       3d23s21
           if(btest(ibc(jjol),i))then                                   3d23s21
            isi=ism(i)                                                  4d8s20
            igi=irel(i)-1                                               4d8s20
            ih0u=ih0a(isi)+igi*(irefo(isi)+1)                            4d8s20
            sum=bc(ih0u)                                                 4d8s20
            do ir=1,nrootz                                               4d8s20
             do j=0,ncsf(iargl)-1                                         4d8s20
              gl(jps+j,ir)=gl(jps+j,ir)+sum*vecr(krps+j,ir)             3d24s21
             end do                                                      4d8s20
            end do                                                       4d8s20
           end if                                                       3d23s21
          end do                                                        4d7s20
         else if(nnot.eq.2)then                                         4d7s20
          nhit=1                                                        3d24s21
          isi=ism(nab(1))                                               4d8s20
          igi=irel(nab(1))-1                                            4d8s20
          igj=irel(nab(2))-1                                            4d8s20
          ih0u=ih0a(isi)+igi+irefo(isi)*igj                             4d8s20
          sum=bc(ih0u)                                                  4d8s20
          call xtimesn(ncsf(iargl),ncsf(iargr),ncsfmid,nrootz,iwpb,iwpk,  12d4s20
     $         vecr(krps,1),ndr,gl(jps,1),ndl,sum,1d0,bc,ibc)           11d10s22
         end if
        end if                                                          4d7s20
        jps=jps+ncsf(iargl)                                              4d7s20
       end do                                                           4d7s20
   12  continue                                                         4d16s20
       jfs=jfsn                                                         4d16s20
       do kri=kr,min(nfcnr,kr+mynprocg-1)                               3d23s21
        nclor=ibasisr(1,kri)                                             10d15s20
        nclorp=nclor+1                                                  4d8s20
        iargr=nclorp-mdon                                                4d8s20
        krps=krps+ncsf(iargr)                                           3d24s21
       end do                                                           4d8s20
      end do                                                            4d7s20
c
      ipss=1                                                            4d16s20
      ifs=1                                                             4d16s20
      jpss=1                                                            4d16s20
      jfs=1                                                             4d16s20
      do kf=1+mynowprog,nkeep,mynprocg                                  3d24s21
       kff=ikeep(kf)                                                    3d24s21
       nclok=ibasisk(1,kff)                                             3d24s21
       nclokp=nclok+1                                                     7d11s19
       karg=nclokp-mdon                                                  7d11s19
       nopenk=nec-2*nclok                                                 7d11s19
       kkc=iptkb(1,nclokp)+ibasisk(2,kff)-1                             3d24s21
       kko=iptkb(2,nclokp)+ibasisk(3,kff)-1                             3d24s21
c
c     gamma^ki_kl times vec.
c                                                                       4d7s20
       nrow=ncsf(karg)*nrootz                                           4d7s20
       do isb=1,nsymb                                                   4d17s20
        do isa=1,isb                                                    4d17s20
         igkv(isa,isb)=ibcoff                                           4d17s20
         if(isa.eq.isb)then                                             4d17s20
          nrr=(irefo(isa)*(irefo(isa)+1))/2                             4d17s20
         else                                                           4d17s20
          nrr=irefo(isa)*irefo(isb)                                     4d17s20
         end if                                                         4d17s20
         ibcoff=igkv(isa,isb)+nrow*nrr                                  4d17s20
        end do                                                          4d17s20
       end do                                                           4d17s20
       call enough('hccsfrilr.  1',bc,ibc)
       do i=igkv(1,1),ibcoff-1                                               4d7s20
        bc(i)=0d0                                                       4d7s20
       end do                                                           4d7s20
       ips=ipss                                                         4d16s20
       ifsn=nfcnr+1                                                      4d16s20
       nhit=0                                                           3d24s21
       do ifr=ifs,nfcnr                                                     4d7s20
        nclor=ibasisr(1,ifr)                                               3d11s20
        if(nclor.gt.nclokp)go to 1
        nclorp=nclor+1                                                     7d11s19
        iargr=nclorp-mdon                                                  7d11s19
        if(iabs(nclor-nclok).le.1)then                                   4d7s20
         if(ifsn.gt.ifr)then                                             4d16s20
          ifsn=ifr                                                       4d16s20
          ipss=ips                                                      4d16s20
         end if                                                         4d16s20
         nopenr=nec-2*nclor                                             3d23s21
         iic=iptrbr(1,nclorp)+ibasisr(2,ifr)-1                              5d7s20
         iio=iptrbr(2,nclorp)+ibasisr(3,ifr)-1                              5d7s20
         call gandc(ibc(kkc),ibc(kko),ibc(iic),ibc(iio),nopenk,nopenr,   12d4s20
     $        karg,iargr,ncsf,norb,ixw1,ixw2,nnot,nab,iwpb,iwpk,ncsfmid,11d14s22
     $        bc,ibc)                                                   11d14s22
         if(nnot.eq.1)then                                              4d7s20
          nhit=nhit+1                                                   3d24s21
          innot1=1                                                      4d8s20
          do i=1,norb                                                   3d24s21
           if(btest(ibc(kkc),i))then                                    3d23s21
            isa=ism(i)                                                  4d17s20
            iga=irel(i)                                                 4d17s20
            icol=((iga*(iga+1))/2)-1                                       4d7s20
            jgkv=igkv(isa,isa)+ncsf(karg)*nrootz*icol                    4d17s20
            do ir=1,nrootz                                               4d7s20
             do j=0,ncsf(karg)-1                                          4d7s20
              bc(jgkv+j)=bc(jgkv+j)+vecr(ips+j,ir)                       4d7s20
             end do                                                      4d7s20
             jgkv=jgkv+ncsf(karg)                                        4d7s20
            end do                                                       4d7s20
           end if                                                       3d23s21
           if(btest(ibc(kko),i))then                                    3d23s21
            isa=ism(i)                                                  4d17s20
            iga=irel(i)                                                 4d17s20
            icol=((iga*(iga+1))/2)-1                                       4d7s20
            jgkv=igkv(isa,isa)+ncsf(karg)*nrootz*icol                    4d17s20
            do ir=1,nrootz                                               4d7s20
             do j=0,ncsf(karg)-1                                          4d7s20
              bc(jgkv+j)=bc(jgkv+j)+vecr(ips+j,ir)*0.5d0                 4d7s20
             end do                                                      4d7s20
             jgkv=jgkv+ncsf(karg)                                        4d7s20
            end do                                                       4d7s20
           end if                                                       3d23s21
          end do                                                        4d7s20
         else if(nnot.eq.2)then                                         4d7s20
          nhit=nhit+1                                                   3d24s21
          if(nab(1).gt.nab(2))then                                      4d7s20
           ione=nab(2)                                                  4d17s20
           itwo=nab(1)                                                  4d17s20
          else                                                          4d7s20
           ione=nab(1)                                                  4d17s20
           itwo=nab(2)                                                  4d17s20
          end if                                                        4d7s20
          isa=ism(ione)                                                 4d17s20
          iga=irel(ione)-1                                              4d17s20
          isb=ism(itwo)                                                 4d17s20
          igb=irel(itwo)-1                                              4d17s20
          if(isa.eq.isb)then                                            4d17s20
           icol=((igb*(igb+1))/2)+iga                                   4d17s20
          else                                                          4d17s20
           icol=iga+irefo(isa)*igb                                      4d17s20
          end if                                                        4d17s20
          jgkv=igkv(isa,isb)+ncsf(karg)*nrootz*icol                     4d17s20
          call xtimesn(ncsf(karg),ncsf(iargr),ncsfmid,nrootz,iwpb,iwpk, 12d4s20
     $          vecr(ips,1),ndr,bc(jgkv),ncsf(karg),0.5d0,1d0,bc,ibc)   11d10s22
         end if
        end if                                                          4d7s20
        ips=ips+ncsf(iargr)                                              4d7s20
       end do                                                           4d7s20
    1  continue                                                         4d16s20
       if(ifsn.le.nfcnr)ifs=ifsn                                        3d23s21
c
c     mult by (ij|kl)
c
       do isb=1,nsymb                                                   4d17s20
        do isa=1,isb                                                    4d17s20
         iekv(isa,isb)=ibcoff                                           4d17s20
         if(isa.eq.isb)then                                             4d17s20
          nrr=(irefo(isa)*(irefo(isa)+1))/2                             4d17s20
         else                                                           4d17s20
          nrr=irefo(isa)*irefo(isb)                                     4d17s20
         end if                                                         4d17s20
         ibcoff=iekv(isa,isb)+nrr*nrow                                  4d17s20
        end do                                                          4d17s20
       end do                                                           4d17s20
       call enough('hccsfrilr.  2',bc,ibc)
       do i=iekv(1,1),ibcoff-1                                          4d17s20
        bc(i)=0d0                                                       4d7s20
       end do                                                           4d7s20
       do isl=1,nsymb                                                   4d17s20
        do isk=1,isl                                                    4d17s20
         if(isk.eq.isl)then                                             4d17s20
          nrrkl=(irefo(isk)*(irefo(isk)+1))/2                           4d17s20
         else                                                           4d17s20
          nrrkl=irefo(isk)*irefo(isl)                                   4d17s20
         end if                                                         4d17s20
         if(nrrkl.gt.0)then                                             4d17s20
          islk=multh(isl,isk)                                           4d17s20
          do isj=1,nsymb                                                 4d17s20
           isi=multh(isj,islk)                                          4d17s20
           if(isi.le.isj)then                                           4d17s20
            if(isi.eq.isj)then                                          4d17s20
             nrrij=(irefo(isj)*(irefo(isj)+1))/2                        4d17s20
            else                                                        4d17s20
             nrrij=irefo(isi)*irefo(isj)                                4d17s20
            end if                                                      4d17s20
            if(nrrij.gt.0)then                                          4d17s20
             i2eu=ifind2(isi,isj,isk,isl,icase)                           4d17s20
             if(icase.eq.1)then                                         4d17s20
              do kl=0,nrrkl-1                                            4d17s20
               jgkv=igkv(isk,isl)+nrow*kl                                4d17s20
               do ij=0,nrrij-1                                           4d17s20
                jekv=iekv(isi,isj)+nrow*ij                               4d17s20
                xint=bc(i2e(i2eu)+ij+nrrij*kl)                          4d17s20
                do irow=0,nrow-1                                         4d17s20
                 bc(jekv+irow)=bc(jekv+irow)+bc(jgkv+irow)*xint          4d17s20
                end do                                                   4d17s20
               end do                                                    4d17s20
              end do                                                     4d17s20
             else if(icase.eq.4)then                                    9d2s20
              do kl=0,nrrkl-1                                           9d2s20
               jgkv=igkv(isk,isl)+nrow*kl                               9d2s20
               do j=0,irefo(isj)-1                                      9d2s20
                do i=0,irefo(isi)-1                                     9d2s20
                 ij=i+irefo(isi)*j                                      9d2s20
                 ji=j+irefo(isj)*i                                      9d2s20
                 jekv=iekv(isi,isj)+nrow*ij                             9d2s20
                 xint=bc(i2e(i2eu)+ji+nrrij*kl)                         9d2s20
                 do irow=0,nrow-1                                       9d2s20
                  bc(jekv+irow)=bc(jekv+irow)+bc(jgkv+irow)*xint        9d2s20
                 end do                                                 9d2s20
                end do                                                  9d2s20
               end do                                                   9d2s20
              end do                                                    9d2s20
             else if(icase.eq.2)then                                    4d17s20
              do l=0,irefo(isl)-1                                       4d17s20
               do k=0,irefo(isk)-1                                      4d17s20
                kl=k+irefo(isk)*l                                       4d17s20
                jgkv=igkv(isk,isl)+nrow*kl                                4d17s20
                lk=l+irefo(isl)*k                                       4d17s20
                do j=0,irefo(isj)-1                                     4d17s20
                 do i=0,irefo(isi)-1                                    4d17s20
                  ij=i+irefo(isi)*j                                     4d17s20
                  ji=j+irefo(isj)*i                                     4d17s20
                  jekv=iekv(isi,isj)+nrow*ij                               4d17s20
                  xint=bc(i2e(i2eu)+ji+nrrij*lk)                            4d17s20
                  do irow=0,nrow-1                                         4d17s20
                   bc(jekv+irow)=bc(jekv+irow)+bc(jgkv+irow)*xint          4d17s20
                  end do                                                4d17s20
                 end do                                                 4d17s20
                end do                                                   4d17s20
               end do                                                    4d17s20
              end do                                                     4d17s20
             else                                                       4d17s20
              write(6,*)('icase for '),isi,isj,isk,isl,(' is '),icase
              call dws_sync
              call dws_finalize
              stop
             end if                                                     4d17s20
            end if                                                      4d17s20
           end if                                                       4d17s20
          end do                                                         4d17s20
         end if                                                         4d17s20
        end do                                                          4d17s20
       end do
c
c     multiply by density into gx
c
       jps=jpss                                                         4d16s20
       jfsn=nfcnl+1                                                     3d24s21
       nhit=0
       do jfl=jfs,nfcnl                                                   4d16s20
        nclol=ibasisl(1,jfl)                                            3d23s21
        if(nclol.gt.nclokp)go to 2                                      4d16s20
        nclopl=nclol+1                                                  6d11s19
        iargl=nclopl-mdon                                                4d7s20
        if(iabs(nclol-nclok).le.1)then                                  4d7s20
         if(jfsn.gt.jfl)then                                             4d16s20
          jfsn=jfl                                                       4d16s20
          jpss=jps                                                      4d16s20
         end if                                                         4d16s20
         nopenl=nec-2*nclol
         jjc=iptrbl(1,nclopl)+ibasisl(2,jfl)-1                              5d7s20
         jjo=iptrbl(2,nclopl)+ibasisl(3,jfl)-1                          3d23s21
         call gandc(ibc(jjc),ibc(jjo),ibc(kkc),ibc(kko),nopenl,nopenk,  3d23s21
     $        iargl,karg,ncsf,norb,ixw1,ixw2,nnot,nab,iwpb,iwpk,ncsfmid,11d14s22
     $        bc,ibc)                                                   11d14s22
         if(nnot.eq.1)then                                              4d7s20
          nhit=1
          do i=1,norb                                                   3d24s21
           if(btest(ibc(jjc),i))then                                    3d23s21
            isa=ism(i)                                                  4d17s20
            iga=irel(i)                                                 4d17s20
            icol=((iga*(iga+1))/2)-1                                     4d17s20
            jekv=iekv(isa,isa)+ncsf(iargl)*nrootz*icol                  3d25s21
            do ir=1,nrootz                                               4d7s20
             do j=0,ncsf(iargl)-1                                          4d7s20
              gl(jps+j,ir)=gl(jps+j,ir)+2d0*bc(jekv+j)                   4d8s20
             end do                                                      4d7s20
             jekv=jekv+ncsf(iargl)                                        4d7s20
            end do                                                       4d7s20
           end if                                                       3d23s21
           if(btest(ibc(jjo),i))then                                    3d23s21
            isa=ism(i)                                                  4d17s20
            iga=irel(i)                                                 4d17s20
            icol=((iga*(iga+1))/2)-1                                     4d17s20
            jekv=iekv(isa,isa)+ncsf(iargl)*nrootz*icol                  3d25s21
            do ir=1,nrootz                                               4d7s20
             do j=0,ncsf(iargl)-1                                          4d7s20
              gl(jps+j,ir)=gl(jps+j,ir)+bc(jekv+j)                       4d8s20
             end do                                                      4d7s20
             jekv=jekv+ncsf(iargl)                                        4d7s20
            end do                                                       4d7s20
           end if                                                       3d23s21
          end do                                                        4d7s20
         else if(nnot.eq.2)then                                         4d7s20
          nhit=1
          if(nab(1).gt.nab(2))then                                      4d7s20
           ione=nab(2)                                                  4d17s20
           itwo=nab(1)                                                  4d17s20
          else                                                          4d7s20
           ione=nab(1)                                                  4d17s20
           itwo=nab(2)                                                  4d17s20
          end if                                                        4d7s20
          isa=ism(ione)                                                 4d17s20
          iga=irel(ione)-1                                              4d17s20
          isb=ism(itwo)                                                 4d17s20
          igb=irel(itwo)-1                                              4d17s20
          if(isa.eq.isb)then                                            4d17s20
           icol=((igb*(igb+1))/2)+iga                                   4d17s20
          else                                                          4d17s20
           icol=iga+irefo(isa)*igb                                      4d17s20
          end if                                                        4d17s20
          jekv=iekv(isa,isb)+ncsf(karg)*nrootz*icol                     4d17s20
          call xtimesn(ncsf(iargl),ncsf(karg),ncsfmid,nrootz,iwpb,iwpk, 12d4s20
     $          bc(jekv),ncsf(karg),gl(jps,1),ndl,1d0,1d0,bc,ibc)       11d10s22
         end if
        end if                                                          4d7s20
        jps=jps+ncsf(iargl)                                              4d7s20
       end do                                                           4d7s20
    2  continue                                                         4d16s20
       if(jfsn.le.nfcnl)jfs=jfsn                                        3d23s21
       ibcoff=igkv(1,1)                                                 4d17s20
      end do                                                            4d7s20
      nwds=ndl*nrootz                                                   3d23s21
      call dws_gsumf(gl,nwds)
      ibcoff=ibcoffo                                                    4d8s20
      return
      end                                                               7d11s19
