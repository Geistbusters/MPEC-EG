c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine vdmo2sobk4(idv4,iorb,multh,nfdat,isymw,ntype,nbasispc,
     $     iaddr,nfcn,nrootu,nbasisp,srh,idorel,bc,ibc)                 11d14s22
      implicit real*8 (a-h,o-z)
c
      include "common.store"
      include "common.hf"
      include "common.mrci"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      dimension idv4(2,4,8),iorb(8),multh(8,8),nfdat(5,4,*),            12s17s21
     $     nfcn(36,2),iaddr(36,16),nbasisp(*),nbasispc(*),noff(8),      1d25s22
     $     jdv4(4,8),idv4l(4)                                           2d17s22
c
      do isb=1,nsymb
       noff(isb)=idoubo(isb)+irefo(isb)                                 12d20s20
      end do
      if(idorel.eq.2)then                                               1d31s22
       nerimat=ntype*2                                                  1d31s22
      else                                                              1d31s22
       nerimat=ntype*4
      end if                                                            1d31s22
      nn=(nsymb*(nsymb+1))/2                                            11d21s18
      xnorm=0d0
      do i=1,nn                                                         11d21s18
       nfcn(i,1)=0                                                      12d18s21
       nfcn(i,2)=0                                                      12d18s21
      end do                                                            12d18s21
      xnan=-2d0
      iaddstart=ibcoff                                                  12d14s18
      do isb=1,nsymb                                                    11d21s18
       nfsum=nfdat(2,1,isb)                                             1d26s22
       do l=2,4                                                         1d26s22
        nfsum=nfsum+nfdat(2,l,isb)                                      1d26s22
       end do                                                           1d26s22
       nfsum=nfsum*nrootu                                               1d26s22
       isbv12=multh(isb,isymw)                                          12d20s20
       do isbv1=1,nsymb                                                 11d21s18
        isbv2=multh(isbv1,isbv12)                                       11d21s18
        if(isbv1.le.isbv2)then                                          11d21s18
         ii=((isbv2*(isbv2-1))/2)+isbv1                                 11d21s18
         nvvx=nbasisp(isbv1)*nbasisp(isbv2)                             12d18s21
         if(isbv1.eq.isbv2)then                                         12d20s20
          nvisv=nvirt(isbv1)                                            12d20s20
          nvnotv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                      12d20s20
         else                                                           12d20s20
          nvisv=0                                                       12d20s20
          nvnotv=nvirt(isbv1)*nvirt(isbv2)                              12d20s20
         end if                                                         12d20s20
         nfcn(ii,2)=isb                                                 11d10s20
         nfcn(ii,1)=0                                                   12d18s21
         if(max(nvisv,nvnotv).gt.0)nfcn(ii,1)=nfsum                     1d26s22
         do in=1,nerimat                                                11d10s20
          iaddr(ii,in)=ibcoff                                           12d18s21
          ibcoff=iaddr(ii,in)+nfcn(ii,1)*nvvx                           12d18s21
         end do                                                         11d10s20
         call enough('vdmo2sobk4.  1',bc,ibc)
         do i=iaddr(ii,1),ibcoff-1                                      12d18s21
          bc(i)=0d0                                                     12d20s20
         end do                                                         12d14s18
        end if                                                          11d21s18
       end do                                                           11d21s18
      end do                                                            11d21s18
      naddr=ibcoff-iaddstart                                            12d14s18
      ibcoffo=ibcoff                                                    11d21s18
      do lb=1,4                                                         1d25s22
       ioff=idv4(2,lb,isymw)                                            1d25s22
       do isbv=1,nsymb                                                   1d25s22
        jdv4(lb,isbv)=ioff                                              1d25s22
        ioff=ioff+nvirt(isbv)*nrootu*nfdat(2,lb,isymw)*ntype            1d25s22
       end do                                                           1d25s22
      end do                                                            1d25s22
      do isb=1,nsymb
       isbv12=multh(isb,isymw)                                          12d20s20
       do l=1,4                                                         2d17s22
        idv4l(l)=idv4(1,l,isb)                                          2d17s22
       end do                                                           2d17s22
       do isbv1=1,nsymb                                                 11d15s18
        isbv2=multh(isbv1,isbv12)                                       11d15s18
        if(isbv1.le.isbv2)then                                          11d15s18
         iiii=((isbv2*(isbv2-1))/2)+isbv1                               11d21s18
         idell=nbasisp(isbv1)*nbasisp(isbv2)*nfcn(iiii,1)
c     dimension of idv4 is nrowx,ncol=nvv&root,nfdatb&ntype
         if(isbv1.eq.isbv2)then
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         12d19s20
          nrowx=nvirt(isbv1)*nrootu                                     12d19s20
          ioff=0                                                        1d25s22
          do lb=1,4                                                     1d25s22
           ncol=nfdat(2,lb,isb)*ntype                                   1d25s22
           if(min(nrowx,ncol).gt.0)then                                  12d20s20
            ncolx=ncol*nrootu                                            12d20s20
            nn=nrootu*nfdat(2,lb,isb)                                      12d20s20
            call ilimts(ncolx,nvirt(isbv1),mynprocg,mynowprog,il,ih,i1s,1d26s22
     $           i1e,i2s,i2e)                                           1d26s22
            if(i2s.eq.i2e)then                                          1d26s22
             n1d=i1e+1-i1s                                              1d26s22
             i1sm=i1s-1                                                    12d20s20
            else                                                        1d26s22
             n1d=ncolx                                                  1d26s22
             i1sm=0                                                     2d14s22
            end if                                                      1d26s22
            nn2=i2e+1-i2s                                               1d26s22
            itmp=ibcoff                                                   12d20s20
            ndel=nn2*nbasisp(isbv2)*n1d                                 1d26s22
            ibcoff=itmp+ndel*2                                            12d20s20
            call enough('vdmo2sobk4.  2',bc,ibc)
            do iz=itmp,ibcoff-1                                           12d20s20
             bc(iz)=0d0                                                   12d20s20
            end do                                                        12d20s20
            i10=i1s                                                     1d26s22
            i1n=ncolx                                                   1d26s22
            do i2=i2s,i2e                                               1d26s22
             if(i2.eq.i2e)i1n=i1e                                       1d26s22
             iv=i2-1                                                    1d26s22
             ivm=i2-i2s                                                 1d26s22
             jorb1=iorb(isbv1)+nbasispc(isbv1)*(noff(isbv1)+iv)          12d20s20
             jorb2=jorb1+nbasisp(isbv1)                                   12d20s20
             do i1=i10,i1n                                              1d26s22
              i1m=i1-1                                                    12d20s20
              iadf=jdv4(lb,isbv1)+iv+nvirt(isbv1)*i1m                   1d25s22
              xx=bc(iadf)
              iad1=itmp+nbasisp(isbv2)*(ivm+nn2*(i1m-i1sm))             2d15s22
              iad2=iad1+ndel                                              12d20s20
              do id=0,nbasisp(isbv2)-1                                    12d20s20
               bc(iad1+id)=bc(iad1+id)+xx*bc(jorb1+id)                    12d20s20
               bc(iad2+id)=bc(iad2+id)+xx*bc(jorb2+id)                    12d20s20
              end do                                                      12d20s20
             end do                                                       12d20s20
             i10=1                                                      2d15s22
            end do                                                      1d26s22
            do i2=i2s,i2e                                               1d26s22
             iv=i2-1                                                    1d26s22
             ivm=i2-i2s                                                 1d26s22
             jorb1=iorb(isbv1)+nbasispc(isbv1)*(noff(isbv1)+iv)         1d26s22
             jorb2=jorb1+nbasisp(isbv1)                                 1d26s22
             do kd=0,1                                                    12d20s20
              do icol=0,n1d-1                                             12d20s20
               ifrm=itmp+nbasisp(isbv2)*(ivm+nn2*(icol+n1d*kd))         1d26s22
               idecode=icol+i1sm                                          12d20s20
               it=idecode/nn                                              12d20s20
               itp=it+1                                                   12d20s20
               lcol=idecode-nn*it+ioff                                  1d25s22
               if(idorel.eq.2)then                                       1d31s22
                in=itp+ntype*kd                                         1d31s22
                if(kd.eq.0)then                                         1d31s22
                 jorb=jorb2                                             1d31s22
                else                                                    1d31s22
                 jorb=jorb1                                             1d31s22
                end if                                                  1d31s22
                do id=0,nbasisp(isbv2)-1                                   12d20s20
                 iad=iaddr(iiii,in)+nbasisp(isbv1)*(id                   12d20s20
     $             +nbasisp(isbv2)*lcol)                                12d20s20
                 do ic=0,nbasisp(isbv1)-1                                  12d20s20
                  bc(iad+ic)=bc(iad+ic)+bc(ifrm+id)*bc(jorb+ic)         1d31s22
                 end do                                                 1d31s22
                end do
               else                                                     1d31s22
                in1=itp+ntype*2*kd                                         12d20s20
                in2=in1+ntype                                              12d20s20
                do id=0,nbasisp(isbv2)-1                                   12d20s20
                 iad1=iaddr(iiii,in1)+nbasisp(isbv1)*(id                   12d20s20
     $             +nbasisp(isbv2)*lcol)                                12d20s20
                 iad2=iaddr(iiii,in2)+nbasisp(isbv1)*(id                   12d20s20
     $               +nbasisp(isbv2)*lcol)                                12d20s20
                 do ic=0,nbasisp(isbv1)-1                                  12d20s20
                  bc(iad1+ic)=bc(iad1+ic)+bc(ifrm+id)*bc(jorb1+ic)         12d20s20
                  bc(iad2+ic)=bc(iad2+ic)+bc(ifrm+id)*bc(jorb2+ic)         12d20s20
                 end do                                                    12d20s20
                end do                                                     12d20s20
               end if                                                   1d31s22
              end do                                                      12d20s20
             end do                                                       12d20s20
            end do                                                        12d20s20
            ibcoff=itmp                                                   12d20s20
            ioff=ioff+nn                                                1d27s22
           end if                                                        12d20s20
          end do                                                        1d25s22
          isw=0                                                         12d19s20
         else                                                           11d15s18
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 12d19s20
          isw=1                                                         12d19s20
         end if                                                         11d15s18
         if(nvv.gt.0)then                                               12d20s20
          nrowx=nvv*nrootu
          ioff=0                                                         12d19s20
          do l=1,4                                                       12d19s20
           if(nfdat(2,l,isb).gt.0)then                                   12d19s20
            ncol=nfdat(2,l,isb)*ntype
            nn=nrootu*nfdat(2,l,isb)                                     12d20s20
            ncolx=ncol*nrootu                                            12d19s20
            call ilimts(ncolx,nvv,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,  12d19s20
     $          i2e)                                                    12d19s20
c     root,nfdat,ntype
            if(i2s.eq.i2e)then                                           12d20s20
             n1d=i1e+1-i1s                                               12d20s20
             i1sm=i1s-1                                                 2d14s22
            else                                                         12d20s20
             n1d=ncolx                                                   12d20s20
             i1sm=0                                                     2d14s22
            end if                                                       12d20s20
            itmp=ibcoff                                                  12d20s20
            ndel=n1d*nvirt(isbv1)*nbasisp(isbv2)                         12d20s20
            ntmp=itmp+ndel*2                                             12d20s20
            ibcoff=ntmp+nvirt(isbv1)                                     12d20s20
            call enough('vdmo2sobk4.  3',bc,ibc)
            do iz=itmp,ntmp-1                                            12d20s20
             bc(iz)=0d0                                                  12d20s20
            end do                                                       12d20s20
            do iz=0,nvirt(isbv1)-1                                       12d20s20
             ibc(ntmp+iz)=0d0                                            12d20s20
            end do                                                       12d20s20
            i10=i1s                                                      12d19s20
            i1n=ncolx                                                    12d19s20
            do i2=i2s,i2e                                                12d19s20
             if(i2.eq.i2e)i1n=i1e                                        12d19s20
             i2m=i2-1
c
c     if isbv1=isbv2, i2=((iv2*(iv2-1))/2)+iv1+1,
c     2*(i2-iv1-1)=iv2**2-iv2=iv2**2-iv2+0.25-0.25=(iv2-0.5)**2-0.25.
c     2*(i2-1)+0.25 ge 2*(i2-iv1-1)+0.25=(iv2-0.5)**2
c     sqrt(2*(i2-1)+0.25).ge.iv2-0.5
c     if isbv1.ne.isbv2, i2=iv1+nvirt(isbv1)*iv2+1
c
             iv2r=i2m/nvirt(isbv1)                                       12d19s20
             iv1r=i2m-nvirt(isbv1)*iv2r                                  12d19s20
             xx=sqrt(dfloat(2*i2m)+0.25d0)+0.500001d0                    12d19s20
             iv2t=int(xx)                                                12d19s20
             iv1t=i2m-((iv2t*(iv2t-1))/2)                                12d19s20
             iv2=iv2t+isw*(iv2r-iv2t)                                    12d19s20
             iv1=iv1t+isw*(iv1r-iv1t)                                    12d19s20
             ibc(ntmp+iv1)=1                                             12d20s20
             jorb21=iorb(isbv2)+nbasispc(isbv2)*(noff(isbv2)+iv2)       12d20s20
             jorb22=jorb21+nbasisp(isbv2)                                12d19s20
             do i1=i10,i1n                                               12d19s20
              i1m=i1-1                                                   12d19s20
              iad=idv4l(l)+i2m+nvv*i1m                                  2d17s22
              iadt=itmp+nbasisp(isbv2)*(iv1+nvirt(isbv1)*(i1m-i1sm))    2d15s22
              iadt2=iadt+ndel                                            12d20s20
              do id=0,nbasisp(isbv2)-1                                   12d19s20
               bc(iadt+id)=bc(iadt+id)+bc(iad)*bc(jorb21+id)             12d20s20
               bc(iadt2+id)=bc(iadt2+id)+bc(iad)*bc(jorb22+id)           12d20s20
              end do                                                     12d19s20
             end do                                                      12d19s20
             i10=1                                                       12d20s20
            end do                                                       12d19s20
            do iv1=0,nvirt(isbv1)-1                                      12d20s20
             if(ibc(ntmp+iv1).ne.0)then                                  12d20s20
              jorb11=iorb(isbv1)+nbasispc(isbv1)*(noff(isbv1)+iv1)      12d20s20
              jorb12=jorb11+nbasisp(isbv1)                               12d20s20
              do kd=0,1                                                  12d20s20
               do icol=0,n1d-1                                           12d20s20
                ifrm=itmp+nbasisp(isbv2)*(iv1+nvirt(isbv1)*(icol        12d20s20
     $               +n1d*kd))                                          12d20s20
c     icol is r,nf,spin,
                idecode=icol+i1sm
                it=idecode/nn                                            12d20s20
                itp=it+1
                idecode=idecode-nn*it                                    12d20s20
                lcol=idecode+ioff                                        12d20s20
                if(idorel.eq.2)then
                 in=itp+ntype*kd                                         1d31s22
                 if(kd.eq.0)then                                         1d31s22
                  jorb=jorb12                                             1d31s22
                 else                                                    1d31s22
                  jorb=jorb11                                           1d31s22
                 end if                                                  1d31s22
                 do id=0,nbasisp(isbv2)-1                                  12d20s20
                  iad=iaddr(iiii,in)+nbasisp(isbv1)*(id                 1d31s22
     $               +nbasisp(isbv2)*lcol)                              12d20s20
                  do ic=0,nbasisp(isbv1)-1                                 12d20s20
                   bc(iad+ic)=bc(iad+ic)+bc(ifrm+id)*bc(jorb+ic)        1d31s22
                  end do                                                1d31s22
                 end do                                                 1d31s22
                else                                                    1d31s22
                 in1=itp+ntype*2*kd                                       12d20s20
                 in2=in1+ntype                                            12d20s20
                 do id=0,nbasisp(isbv2)-1                                  12d20s20
                  iad1=iaddr(iiii,in1)+nbasisp(isbv1)*(id                 12d20s20
     $               +nbasisp(isbv2)*lcol)                              12d20s20
                  iad2=iaddr(iiii,in2)+nbasisp(isbv1)*(id                 12d20s20
     $               +nbasisp(isbv2)*lcol)                              12d20s20
                  do ic=0,nbasisp(isbv1)-1                                 12d20s20
                   bc(iad1+ic)=bc(iad1+ic)+bc(ifrm+id)*bc(jorb11+ic)      12d20s20
                   bc(iad2+ic)=bc(iad2+ic)+bc(ifrm+id)*bc(jorb12+ic)      12d20s20
                  end do                                                  12d20s20
                 end do                                                   12d20s20
                end if                                                  1d31s22
               end do                                                    12d20s20
              end do                                                     12d20s20
             end if                                                      12d20s20
            end do                                                       12d20s20
            ibcoff=itmp                                                  12d20s20
            ioff=ioff+nn                                                 12d20s20
            idv4l(l)=idv4l(l)+nvv*nfdat(2,l,isb)*nrootu*ntype           2d17s22
           end if                                                        12d19s20
          end do                                                         12d19s20
         end if                                                         12d20s20
        end if                                                          11d15s18
       end do                                                           11d15s18
      end do
      call dws_gsumf(bc(iaddstart),naddr)                               12d14s18
      ibcoff=ibcoffo                                                    11d21s18
      return
      end
