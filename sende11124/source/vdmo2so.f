c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine vdmo2so(vecin,iorb,multh,ncd,nbasispc,iaddr,nfcn,      12d13s18
     $     nrootu,ntot,nbasisp,idorel,srh,itype,bc,ibc)                 6d7s24
      implicit real*8 (a-h,o-z)
c
c     itype=0: normal operation
c     itype ne 0: going for ders, so don't fold vectors in ao basis     6d7s24
c
      include "common.store"
      include "common.hf"
      include "common.mrci"
      logical lprint
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      dimension vecin(*),iorb(8),multh(8,8),ncd(3,8,2),nbasispc(8),     11d19s18
     $     nfcn(36,3),iaddr(36,2,4),nbasisp(*)                          11d10s20
c
      if(idorel.eq.0.or.idorel.eq.1)then                                10d28s20
       nerimat=1                                                        10d28s20
       nv2=1                                                            2d10s21
      else if(idorel.eq.2.or.idorel.eq.-2.or.idorel.eq.-4)then          10d28s20
       nerimat=3                                                        2d16s21
       nv2=2                                                            2d10s21
      else                                                              10d28s20
       nerimat=4                                                        2d10s21
       nv2=2                                                            2d10s21
      end if                                                            10d28s20
      nn=(nsymb*(nsymb+1))/2                                            11d21s18
      xnorm=0d0
      do ist=1,2
       do i=1,nn                                                         11d21s18
        nfcn(i,ist)=0                                                   11d10s20
        iaddr(i,ist,1)=0                                                  11d10s20
        iaddr(i,ist,2)=0                                                  11d10s20
        iaddr(i,ist,3)=0                                                  11d10s20
        iaddr(i,ist,4)=0                                                  11d10s20
       end do                                                            11d21s18
      end do                                                            11d10s20
      xnan=-2d0
      iaddstart=ibcoff                                                  12d14s18
      do isb=1,nsymb                                                    11d21s18
       nn=nvirt(isb)+irefo(isb)+idoubo(isb)                             11d29s18
       isbv12=multh(isb,isymmrci)                                       11d21s18
       do isbv1=1,nsymb                                                 11d21s18
        isbv2=multh(isbv1,isbv12)                                       11d21s18
        if(isbv1.le.isbv2)then                                          11d21s18
         ii=((isbv2*(isbv2-1))/2)+isbv1                                 11d21s18
         if(isbv1.eq.isbv2.and.itype.eq.0)then                          6d7s24
          nvvs=(nbasisp(isbv1)*(nbasisp(isbv1)+1))/2                    10d31s20
          nvvt=(nbasisp(isbv1)*(nbasisp(isbv1)-1))/2                    10d31s20
          nvvx=nbasisp(isbv1)*nbasisp(isbv1)                            11d5s20
         else                                                           11d23s18
          nvvs=nbasisp(isbv1)*nbasisp(isbv2)                            10d31s20
          nvvt=nvvs                                                     11d23s18
          nvvx=nvvs                                                     11d10s20
         end if                                                         11d23s18
         nfcn(ii,3)=isb                                                 11d10s20
         if(nvvx.gt.0)then                                              6d18s24
          nfcn(ii,1)=nrootu*ncd(1,isb,1)                                6d18s24
          nfcn(ii,2)=nrootu*(ncd(1,isb,2)+ncd(2,isb,2)+ncd(3,isb,2))    6d18s24
         end if                                                         6d18s24
         do in=1,nerimat                                                11d10s20
          iaddr(ii,1,in)=ibcoff                                         11d10s20
          if(in.eq.1.or.in.eq.4)then                                    2d16s21
           ibcoff=iaddr(ii,1,in)+nfcn(ii,1)*nvvs                        11d10s20
          else                                                          11d10s20
           ibcoff=iaddr(ii,1,in)+nfcn(ii,1)*nvvx                        11d10s20
          end if                                                        11d10s20
          iaddr(ii,2,in)=ibcoff                                         11d10s20
          if(in.eq.1.or.in.eq.4)then                                    2d16s21
           ibcoff=iaddr(ii,2,in)+nfcn(ii,2)*nvvt                        11d10s20
          else                                                          11d10s20
           ibcoff=iaddr(ii,2,in)+nfcn(ii,2)*nvvx                        11d10s20
          end if                                                        11d10s20
         end do                                                         11d10s20
         do i=iaddr(ii,1,1),ibcoff-1                                    11d10s20
          bc(i)=xnan                                                    11d10s20
         end do                                                         12d14s18
        end if                                                          11d21s18
       end do                                                           11d21s18
      end do                                                            11d21s18
      naddr=ibcoff-iaddstart                                            12d14s18
      ibcoffo=ibcoff                                                    11d21s18
      ioff=1
      iout=ibcoff                                                       11d15s18
      jout=iout                                                         11d15s18
      do isb=1,nsymb
       isbv12=multh(isb,isymmrci)
       do isbv1=1,nsymb                                                 11d15s18
        isbv2=multh(isbv1,isbv12)                                       11d15s18
        if(isbv1.le.isbv2)then                                          11d15s18
         iiii=((isbv2*(isbv2-1))/2)+isbv1                               11d21s18
         if(isbv1.eq.isbv2)then
          nvvs=(nvirt(isbv1)*(nvirt(isbv1)+1))/2                        11d15s18
          nvvt=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        11d15s18
         else                                                           11d15s18
          nvvs=nvirt(isbv1)*nvirt(isbv2)                                11d15s18
          nvvt=nvvs                                                     11d15s18
         end if                                                         11d15s18
         do ipass=1,2                                                   11d19s18
          kstor=0                                                       12d7s18
          if(ipass.eq.1)then                                            11d19s18
           isx=1                                                        11d19s18
          else                                                          11d19s18
           isx=3                                                        11d19s18
          end if                                                        11d19s18
          do is=1,isx                                                   11d19s18
           nrr=ncd(is,isb,ipass)*nrootu                                 12d13s18
           nwds=nrr*nbasisp(isbv1)*nbasisp(isbv2)*nv2                   2d10s21
           ibcoff=jout+nwds                                               11d15s18
           call enough('vdmo2so.  1',bc,ibc)
           do i=0,nwds-1
            bc(jout+i)=0d0
           end do                                                         11d15s18
           kout=jout                                                    2d10s21
           nrow=nbasisp(isbv2)*nrr                                      2d10s21
           ncol=nvirt(isbv1)
           do jv2=1,nv2                                                 2d10s21
            iv12=0                                                       9d6s19
            do iv2=0,nvirt(isbv2)-1                                        11d15s18
             if(isbv1.eq.isbv2)then                                        11d15s18
              if(ipass.eq.1)then                                         11d19s18
               ivtop=iv2                                                    11d15s18
              else                                                       11d19s18
               ivtop=iv2-1                                               11d19s18
              end if                                                     11d19s18
             else                                                          11d15s18
              ivtop=nvirt(isbv1)-1                                         11d15s18
             end if                                                      11d28s18
             jorb=iorb(isbv2)+nbasispc(isbv2)*(idoubo(isbv2)            11d2s20
     $            +irefo(isbv2)+iv2)                                    11d2s20
             if(jv2.ge.2)jorb=jorb+nbasisp(isbv2)                       2d10s21
             do iv1=0,ivtop                                                11d15s18
              if(mod(iv12,mynprocg).eq.mynowprog)then                    9d6s19
               fact=1d0                                                   11d26s18
               if(isbv1.eq.isbv2)then                                       11d15s18
                if(ipass.eq.1)then                                        11d19s18
                 i12=((iv2*(iv2+1))/2)+iv1                                   11d15s18
                 if(iv1.eq.iv2)fact=srh                                   11d26s18
                else                                                      11d19s18
                 i12=((iv2*(iv2-1))/2)+iv1                                   11d15s18
                end if                                                    11d19s18
               else                                                         11d15s18
                i12=iv1+nvirt(isbv1)*iv2                                    11d15s18
               end if                                                       11d15s18
               do iil=0,ncd(is,isb,ipass)-1                               4d9s19
                do l=0,nrootu-1                                           4d9s19
                 ii=l+nrootu*iil                                          4d9s19
                 iad=ioff+iil+ncd(is,isb,ipass)*i12+ntot*l                4d10s19
                 val=vecin(iad)*fact                                       4d9s19
c
c     iv2 ii iv1
c
                 jad=kout+nbasisp(isbv2)*(ii+nrr*iv1)                   11d2s20
                 do id=0,nbasisp(isbv2)-1                               11d2s20
                  bc(jad+id)=bc(jad+id)+val*bc(jorb+id)                     11d15s18
                 end do
                end do                                                    4d9s19
               end do                                                       11d15s18
              end if                                                     9d6s19
              iv12=iv12+1                                                9d6s19
             end do                                                        11d15s18
            end do                                                         11d15s18
            itmp=ibcoff                                                 2d10s21
            ibcoff=itmp+nrow*ncol                                       2d10s21
            call enough('vdmo2so.  2',bc,ibc)
c     ncol=nvirt(isbv1)
c     nrow: nbasisp(isbv2),root,ncd,
            do i=0,ncol-1                                                  11d15s18
             do j=0,nrow-1                                                 11d15s18
              ji=kout+j+nrow*i                                             11d15s18
              ij=itmp+i+ncol*j                                             11d15s18
              bc(ij)=bc(ji)                                                11d15s18
             end do                                                        11d15s18
            end do                                                         11d15s18
            do i=0,nrow*ncol-1                                          2d10s21
             bc(kout+i)=bc(itmp+i)                                      2d10s21
            end do                                                      2d10s21
            ibcoff=itmp                                                 2d10s21
            kout=kout+nrow*ncol                                         2d10s21
           end do                                                       2d10s21
           do in=0,nerimat-1                                            2d10s21
            inp=in+1                                                    2d10s21
            if(nrow*ncol.gt.0)then                                       11d28s18
c     this is no longer correct since || : gsumf is required
c
c     iv1 iv2 ii
c
             itmp=ibcoff                                                2d10s21
             ibcoff=itmp+nrow*ncol                                      2d10s21
             call enough('vdmo2so.  3',bc,ibc)
             jorb=iorb(isbv1)+nbasispc(isbv1)*(idoubo(isbv1)            11d2s20
     $           +irefo(isbv1))                                         11d29s18
             if(nbasispc(isbv1).gt.0)then                               2d10s21
              jorbp=jorb+nbasisp(isbv1)                                 2d10s21
              kout=jout+nrow*ncol                                       2d10s21
              if(in.eq.0)then                                           2d10s21
               call dgemm('n','n',nbasisp(isbv1),nrow,ncol,1d0,          11d2s20
     $        bc(jorb),nbasispc(isbv1),bc(jout),ncol,0d0,               11d15s18
     $        bc(itmp),nbasisp(isbv1),                                  11d2s20
     d' vdmo2so.  1')
              else if(in.eq.1)then                                      2d10s21
               call dgemm('n','n',nbasisp(isbv1),nrow,ncol,1d0,          11d2s20
     $        bc(jorb),nbasispc(isbv1),bc(kout),ncol,0d0,               11d15s18
     $        bc(itmp),nbasisp(isbv1),                                  11d2s20
     d' vdmo2so.  2')
              else if(in.eq.2)then                                      2d16s21
                call dgemm('n','n',nbasisp(isbv1),nrow,ncol,1d0,          11d2s20
     $        bc(jorbp),nbasispc(isbv1),bc(jout),ncol,0d0,              2d16s21
     $        bc(itmp),nbasisp(isbv1),                                  11d2s20
     d' vdmo2so.  3')
              else                                                      2d10s21
               call dgemm('n','n',nbasisp(isbv1),nrow,ncol,1d0,          11d2s20
     $        bc(jorbp),nbasispc(isbv1),bc(kout),ncol,0d0,              2d10s21
     $        bc(itmp),nbasisp(isbv1),                                  11d2s20
     d' vdmo2so.  4')
              end if                                                    2d10s21
             end if                                                     2d10s21
c     this is no longer correct since || : gsumf is required
             if(isbv1.eq.isbv2.and.(in.eq.0.or.in.eq.3)                 6d7s24
     $            .and.itype.eq.0)then                                  6d7s24
              if(ipass.eq.1)then                                          11d23s18
               mvv=(nbasisp(isbv1)*(nbasisp(isbv1)+1))/2                11d2s20
               do i3=0,nbasisp(isbv1)-1                                 11d2s20
                do i2=0,i3                                              11d10s20
                 i23=((i3*(i3+1))/2)+i2                                   11d23s18
                 do i1=0,nrr-1                                            11d23s18
                  iad=itmp+i2+nbasisp(isbv1)*(i3+nbasisp(isbv1)*i1)     2d10s21
                  iadt=itmp+i3+nbasisp(isbv1)*(i2+nbasisp(isbv1)*i1)    2d10s21
                  iad2=iaddr(iiii,1,inp)+i23+mvv*i1                     11d10s20
                  bc(iad2)=bc(iad)+bc(iadt)                               11d23s18
                 end do                                                   11d23s18
                end do                                                    11d23s18
               end do                                                     11d23s18
              else                                                        11d23s18
               mvv=(nbasisp(isbv1)*(nbasisp(isbv1)-1))/2                11d2s20
               do i3=0,nbasisp(isbv1)-1                                 11d2s20
                do i2=0,i3-1                                            11d10s20
                 i23=((i3*(i3-1))/2)+i2                                   11d23s18
                 do i1=0,nrr-1                                            11d23s18
                  iad=itmp+i2+nbasisp(isbv1)*(i3+nbasisp(isbv1)*i1)     2d10s21
                  iadt=itmp+i3+nbasisp(isbv1)*(i2+nbasisp(isbv1)*i1)    2d10s21
                  iad2=iaddr(iiii,2,inp)+i23+mvv*(i1+kstor)             11d10s20
                  bc(iad2)=bc(iad)-bc(iadt)                               11d23s18
                 end do                                                   11d23s18
                end do                                                    11d23s18
               end do                                                     11d23s18
              end if
             else                                                         11d23s18
              do i3=0,nbasisp(isbv2)-1                                  11d2s20
               do i2=0,nbasisp(isbv1)-1                                 11d2s20
                do i1=0,nrr-1                                              11d21s18
                 iad1=itmp+i2+nbasisp(isbv1)*(i3+nbasisp(isbv2)*i1)     2d10s21
                 iad2=iaddr(iiii,ipass,inp)+                            11d10s20
     $               i2+nbasisp(isbv1)*(i3+nbasisp(isbv2)*(kstor+i1))   11d10s20
                 bc(iad2)=bc(iad1)                                         11d21s18
                end do
               end do
              end do                                                       11d21s18
             end if                                                       11d23s18
            end if                                                       11d28s18
           end do                                                       11d2s20
           kstor=kstor+nrr                                              11d21s18
           if(ipass.eq.1)then                                           11d19s18
            nvv=nvvs                                                    11d19s18
           else                                                         11d19s18
            nvv=nvvt                                                    11d19s18
           end if                                                       11d19s18
           ioff=ioff+nvv*ncd(is,isb,ipass)                              4d9s19
          end do                                                        11d19s18
c     this is no longer correct since || : gsumf is required
         end do                                                         11d19s18
        end if                                                          11d15s18
       end do                                                           11d15s18
      end do
      call dws_gsumf(bc(iaddstart),naddr)                               12d14s18
      ibcoff=ibcoffo                                                    11d21s18
      return
      end
