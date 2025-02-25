c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine trans3xden(idorel,ngaus,ibdat,iapair,ipair,id3x,idoubo,11d28s23
     $    irefo,nvirt,nbasisp,nsymb,isnd,ircv,ipt,ipf,iorb,multh,ibufs, 11d28s23
     $     isstor,nrowp,ncolp,bc,ibc)                                   12d1s23
      implicit real*8 (a-h,o-z)
c
c     transform 3 external density to ao basis (except for the internal part)
c     to compute 4 external contribution to orbital response derivative.
c
      include "common.store"
      integer*8 isstor(*)                                               11d29s23
      dimension iapair(3,*),id3x(*),idoubo(*),nvirt(*),nbasisp(*),      11d28s23
     $     irefo(*),isnd(*),ircv(*),iorb(*),multh(8,8),ipt(*),ipf(*),   11d28s23
     $     itsc(8),nrowp(*),ncolp(*),nrowps(8,8),ncolps(8,8)            12d1s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      common/kmfind/invk1(2,8,8,8,2)                                    11d28s23
      data loopx/100000000/
      loop=0
      write(6,*)('Hi, my name is trans3xden ')
      write(6,*)('ibcoff '),ibcoff,ibdat,ngaus,idorel
      ncomp=1                                                           11d27s23
      if(idorel.ne.0)ncomp=2                                            11d27s23
      nn=(ngaus*(ngaus+1))/2                                            2d19s10
      nn2=nn*2                                                          2d19s10
      ipair=ibcoff                                                      2d19s10
      ibcoff=ipair+nn*2                                                 3d4s10
      i12=ibcoff                                                        2d19s10
      ibcoff=i12+nn*4                                                   2d19s10
      call enough('trans3xden.ipair',bc,ibc)
      j12=i12                                                           2d19s10
      jbdat=ibdat-1                                                     2d19s10
      ngaus3=ngaus*3                                                    5d12s10
      ngaus7=ngaus*7                                                    5d3s10
      do i1=1,ngaus                                                     2d19s10
       do i2=1,i1                                                       2d19s10
        ibc(j12)=-((ibc(jbdat+i1)+ibc(jbdat+i2)+2)**2)                  2d19s10
        if(iapair(1,ibc(jbdat+i1+ngaus7)).gt.0)ibc(j12)=ibc(j12)*2      5d3s10
        if(iapair(1,ibc(jbdat+i2+ngaus7)).gt.0)ibc(j12)=ibc(j12)*2      5d3s10
        ibc(j12+nn)=i1                                                  2d19s10
        ibc(j12+nn2)=i2                                                 2d19s10
        j12=j12+1                                                       2d19s10
       end do                                                           2d19s10
      end do                                                            2d19s10
      is12=i12+nn*3                                                     2d19s10
      nn8=nn                                                            2d19s10
      call idsortdws(ibc(i12),ibc(is12),nn8)                            1d18s23
      ii1=i12+nn                                                        2d19s10
      ii2=i12+nn2                                                       2d19s10
      jpair=ipair                                                       2d19s10
      do i=1,nn                                                         2d19s10
       j=ibc(is12+i-1)-1                                                2d19s10
       ibc(jpair)=ibc(ii1+j)                                            2d19s10
       jpair=jpair+1                                                    2d19s10
       ibc(jpair)=ibc(ii2+j)                                            2d19s10
       jpair=jpair+1                                                    2d19s10
      end do                                                            2d19s10
      do i=1,mynprocg                                                   12d1s23
       nrowp(i)=0                                                       12d1s23
       ncolp(i)=0                                                       12d1s23
      end do                                                            12d1s23
      do isb=1,nsymb                                                    12d1s23
       do isa=1,nsymb                                                   12d1s23
        nrowps(isa,isb)=0                                               12d1s23
        ncolps(isa,isb)=0                                               12d5s23
       end do                                                           12d1s23
      end do                                                            12d1s23
      nbnb=nsymb*nsymb                                                  11d27s23
      if(idorel.eq.0.or.idorel.eq.1)then                                11d28s23
       inc=1                                                            11d28s23
      else                                                              11d28s23
       inc=2                                                            11d28s23
      end if                                                            11d28s23
      do i=1,nn                                                         11d30s23
       im=i-1                                                           11d27s23
       ip=mod(im,mynprocg)                                              11d27s23
       ipp=ip+1                                                         12d1s23
       jpair=ipair+2*(i-1)                                              2d19s10
       la=ibc(jbdat+ibc(jpair))                                         3d18s10
       na=(2*la+1)                                                      11d28s23
       lb=ibc(jbdat+ibc(jpair+1))                                       3d18s10
       write(6,*)('address? '),i,jpair,ibc(jpair),ibc(jpair+1)
       nb=(2*lb+1)                                                      11d28s23
       naa=na                                                           5d12s10
       nba=nb                                                           5d12s10
       if(iapair(1,ibc(jbdat+ibc(jpair)+ngaus7)).gt.0)then              5d12s10
        naa=naa*2                                                       5d12s10
       end if                                                           5d12s10
       if(iapair(1,ibc(jbdat+ibc(jpair+1)+ngaus7)).gt.0)then            5d12s10
        nba=nba*2                                                       5d12s10
       end if                                                           5d12s10
       do ib=1,nba                                                      11d28s23
        ibb=ib+ibc(jbdat+ibc(jpair+1)+ngaus3)                           6d4s10
        isb=isstor(ibb)                                                 5d12s10
        do ia=1,naa                                                     11d28s23
         iaa=ia+ibc(jbdat+ibc(jpair)+ngaus3)                            6d4s10
         isa=isstor(iaa)                                                5d12s10
         if(ibb.le.iaa)then                                             1d30s12
          if(ip.eq.mynowprog)then                                       11d28s23
           iad=nrowps(isa,isb)                                          12d1s23
           nrowps(isa,isb)=nrowps(isa,isb)+inc                          12d1s23
          end if                                                        11d28s23
          nrowp(ipp)=nrowp(ipp)+inc                                     12d1s23
         end if                                                         11d27s23
        end do                                                          11d27s23
       end do                                                           11d27s23
      end do                                                            11d27s23
      write(6,*)('nrowp: '),(nrowp(k),k=1,mynprocg)
      mcolps=0                                                          12d5s23
      do ip=0,mynprocg-1                                                11d28s23
       ipp=ip+1                                                         12d1s23
       do isd=1,nsymb                                                    11d28s23
        do isc=1,nsymb                                                   11d28s23
         call ilimts(irefo(isc),nvirt(isd),mynprocg,ip,il,ih,i1s,i1e,   11d28s23
     $        i2s,i2e)                                                  11d28s23
         nhere=ih+1-il                                                  11d28s23
         ncolp(ipp)=ncolp(ipp)+nhere                                    12d1s23
         write(6,*)ip,isc,isd,nhere,ncolp(ipp)
         if(ip.eq.mynowprog)then                                        12d5s23
          ncolps(isc,isd)=ncolps(isc,isd)+nhere                         12d5s23
          mcolps=mcolps+nhere                                           12d5s23
         end if                                                         12d5s23
        end do                                                          11d28s23
       end do                                                           11d28s23
      end do                                                            11d28s23
      write(6,*)('ncolp: '),(ncolp(k),k=1,mynprocg)
      write(6,*)('send/rcv')
      nst=0                                                             11d28s23
      nrt=0                                                             11d28s23
      do ip=0,mynprocg-1                                                11d28s23
       ipp=ip+1
       isnd(ipp)=ncolp(mynowprog+1)*nrowp(ipp)                          12d1s23
       ircv(ipp)=ncolp(ipp)*nrowp(mynowprog+1)                          12d1s23
       write(6,*)ip,isnd(ipp),ncolp(mynowprog+1),nrowp(ipp),            12d1s23
     $      ircv(ipp),ncolp(ipp),nrowp(mynowprog+1)                     12d1s23
       nst=nst+isnd(ipp)                                                11d28s23
       nrt=nrt+ircv(ipp)                                                11d28s23
      end do                                                            11d27s23
      write(6,*)('size of send buffer: '),nst                           11d28s23
      write(6,*)('size of receive buffer: '),nrt                        11d28s23
      nrtrn=0                                                           11d28s23
      do isa=1,nsymb                                                    11d28s23
       do isb=1,isa                                                     11d28s23
        isab=multh(isa,isb)                                             11d28s23
        do isc=1,nsymb                                                  11d28s23
         isd=multh(isc,isab)                                            11d28s23
         nrtrn=nrtrn+nrowps(isa,isb)*irefo(isc)*nbasisp(isd)*ncomp      12d1s23
        end do                                                          11d28s23
       end do                                                           11d28s23
      end do                                                            11d28s23
      write(6,*)('no. words to return: '),nrtrn,nst,nrt                         11d28s23
      ibufs=ibcoff                                                      11d28s23
      ibufr=ibufs+max(nrtrn,nst)                                        11d28s23
      ibcoff=ibufr+nrt                                                  11d28s23
      call enough('trans3xden.bufs',bc,ibc)                             12d1s23
      jbufs=ibufs                                                       11d28s23
      jbufr=0                                                           11d28s23
      do ipp=1,mynprocg                                                 11d28s23
       ipt(ipp)=jbufs                                                   11d28s23
       jbufs=jbufs+isnd(ipp)                                            11d28s23
       ipf(ipp)=jbufr                                                   11d28s23
       jbufr=jbufr+ircv(ipp)                                            11d28s23
       write(6,*)('ipt of '),ipp,(' is '),ipt(ipp)
      end do                                                            11d28s23
      do isa=1,nsymb                                                    11d28s23
       nbba=nbasisp(isa)*ncomp                                          12d5s23
       do isb=1,isa                                                     11d28s23
        isab=multh(isa,isb)                                             11d28s23
        nbb=nbasisp(isb)*ncomp                                          11d28s23
        itmp1=ibcoff                                                    11d28s23
        itmp2=itmp1+nbba*nbb*mcolps                                     12d5s23
        ibcoff=itmp2+nbba*nbb*mcolps                                    12d5s23
        call enough('trans3xder.tmp1',bc,ibc)                           11d28s23
        jtmp1=itmp1                                                     11d28s23
        nn=0                                                            12d1s23
        do isc=1,nsymb                                                  11d28s23
         isd=multh(isab,isc)                                            11d28s23
         call ilimts(irefo(isc),nvirt(isd),mynprocg,mynowprog,il,ih,    11d28s23
     $        i1s,i1e,i2s,i2e)                                          11d28s23
         nhere=ih+1-il                                                  11d28s23
         nn=nn+nhere                                                    12d1s23
         if(nhere.gt.0)then                                             11d28s23
          i2eu=invk1(1,isa,isb,isc,2)                                   11d30s23
          if(isa.eq.isb)then                                            11d28s23
           irow=0                                                       11d27s23
           nrow=(nvirt(isa)*(nvirt(isa)+1))/2                           11d28s23
           do i2=0,nvirt(isa)-1                                         11d28s23
            ff=0.5d0                                                    12d1s23
            do i1=0,i2                                                  11d27s23
             if(i1.eq.i2)ff=1d0                                         11d27s23
             do icol=0,nhere-1                                          11d27s23
              iad2=jtmp1+i2+nvirt(isa)*(i1+nvirt(isa)*icol)             11d28s23
              iad1=jtmp1+i1+nvirt(isa)*(i2+nvirt(isa)*icol)             11d28s23
              iad3=id3x(i2eu)+irow+nrow*icol                            11d28s23
              bc(iad1)=bc(iad3)*ff                                      11d27s23
              bc(iad2)=bc(iad3)*ff                                      11d27s23
              write(6,*)isc,isd,icol,bc(iad3),bc(iad1),iad1-itmp1
             end do                                                     11d27s23
             irow=irow+1                                                11d27s23
            end do                                                      11d27s23
           end do                                                       11d27s23
          else                                                          11d27s23
           irow=0                                                       11d27s23
           nrow=nvirt(isa)*nvirt(isb)                                   11d28s23
           do i2=0,nvirt(isb)-1                                         11d28s23
            do i1=0,nvirt(isa)-1                                        11d28s23
             do icol=0,nhere-1                                          11d27s23
              iad2=jtmp1+i2+nvirt(isb)*(i1+nvirt(isa)*icol)             11d28s23
              iad3=id3x(i2eu)+irow+nrow*icol                            11d28s23
              bc(iad2)=bc(iad3)                                         11d27s23
             end do                                                     11d27s23
             irow=irow+1                                                11d27s23
            end do                                                      11d27s23
           end do                                                       11d27s23
          end if                                                        11d28s23
          jtmp1=jtmp1+nvirt(isb)*nvirt(isa)*nhere                       11d28s23
         end if                                                         11d28s23
        end do                                                          11d28s23
        jorb=iorb(isb)+nbb*(idoubo(isb)+irefo(isb))                     11d28s23
        nx=nvirt(isa)*nn                                                12d1s23
        if(nn.gt.0)then                                                 12d5s23
         call prntm2(bc(itmp1),nvirt(isa)*nvirt(isb),nn,                 12d1s23
     $       nvirt(isa)*nvirt(isb))
         call dgemm('n','n',nbb,nx,nvirt(isb),1d0,bc(jorb),nbb,         12d5s23
     $        bc(itmp1),nvirt(isb),0d0,bc(itmp2),nbb,'trans3xden.1')    12d5s23
         if(idorel.eq.0.or.idorel.eq.1)then                              11d28s23
          do i=0,nx-1                                                     11d28s23
           do j=0,nbb-1                                                   11d28s23
            ji=itmp2+j+nbb*i                                              11d28s23
            ij=itmp1+i+nx*j                                               11d28s23
            bc(ij)=bc(ji)                                                 11d28s23
           end do                                                         11d28s23
          end do                                                          11d28s23
         else                                                            11d28s23
          do i=0,nx-1                                                    11d28s23
           do j=0,nbasisp(isb)-1                                         11d28s23
            ji=itmp2+j+nbb*i                                             11d28s23
            ij=itmp1+i+nx*j                                              11d28s23
            bc(ij)=bc(ji)                                                11d28s23
            ji=ji+nbasisp(isb)                                           11d28s23
            ij=ij+nx*nbasisp(isb)                                        11d28s23
            bc(ij)=bc(ji)                                                11d28s23
           end do                                                        11d28s23
          end do                                                         11d28s23
         end if                                                          11d28s23
c
c     tmp1 is now nvirt(isa),ncolp,nbb.
c     if rel, primitive ints are
c     in=1 is llll, 2 is llss, 3 is ssll, 4 is ssss
c     i.e. ab must either be ll or ss.
c
         jorb=iorb(isa)+nbba*(idoubo(isa)+irefo(isa))                   12d5s23
         nx=nbasisp(isb)*nn                                              12d1s23
         call dgemm('n','n',nbasisp(isa),nx,nvirt(isa),1d0,bc(jorb),    12d5s23
     $        nbba,bc(itmp1),nvirt(isa),0d0,bc(itmp2),nbasisp(isa),     12d5s23
     $       'trans3xden.2')                                            11d28s23
         nhere=nn                                                        12d1s23
         if(.not.(idorel.eq.0.or.idorel.eq.1))then                       11d28s23
          jorb=jorb+nbasisp(isa)                                         11d28s23
          jtmp1=jtmp1+nvirt(isa)*nx                                      11d28s23
          jtmp2=itmp2+nbasisp(isa)*nx                                    11d28s23
          call dgemm('n','n',nbasisp(isa),nx,nvirt(isa),1d0,bc(jorb),   12d5s23
     $         nbba,bc(jtmp1),nvirt(isa),0d0,bc(jtmp2),nbasisp(isa),    12d5s23
     $       'trans3xden.3')                                            11d28s23
          nhere=nhere*2                                                  11d28s23
         end if                                                          11d28s23
         do ip=0,mynprocg-1                                              11d28s23
          ipp=ip+1                                                       11d28s23
          do i=1+ip,nn,mynprocg                                          12d5s23
           jpair=ipair+2*(i-1)                                           11d28s23
           la=ibc(jbdat+ibc(jpair))                                      11d28s23
           na=(2*la+1)                                                   11d28s23
           lb=ibc(jbdat+ibc(jpair+1))                                    11d28s23
           nb=(2*lb+1)                                                   11d28s23
           naa=na                                                        11d28s23
           nba=nb                                                        11d28s23
           if(iapair(1,ibc(jbdat+ibc(jpair)+ngaus7)).gt.0)then           11d28s23
            naa=naa*2                                                    11d28s23
           end if                                                        11d28s23
           if(iapair(1,ibc(jbdat+ibc(jpair+1)+ngaus7)).gt.0)then         11d28s23
            nba=nba*2                                                    11d28s23
           end if                                                        11d28s23
           do ib=1,nba                                                   11d28s23
            ibb=ib+ibc(jbdat+ibc(jpair+1)+ngaus3)                        11d28s23
            iisb=isstor(ibb)                                             11d28s23
            do ia=1,naa                                                  11d28s23
             iaa=ia+ibc(jbdat+ibc(jpair)+ngaus3)                         11d28s23
             iisa=isstor(iaa)                                            11d28s23
             if(ibb.le.iaa)then                                          11d28s23
              if(isa.eq.iisa.and.isb.eq.iisb)then                        11d28s23
               iad=itmp1+iaa-1+nbasisp(isa)*nhere*(ibb-1)                11d28s23
               do ii=0,nhere-1                                           11d28s23
                bc(ipt(ipp))=bc(iad)                                     11d28s23
                ipt(ipp)=ipt(ipp)+1                                      11d28s23
                iad=iad+nbasisp(isa)                                     11d28s23
               end do                                                    11d28s23
              end if                                                     11d28s23
             end if                                                      11d28s23
            end do                                                       11d28s23
           end do                                                        11d28s23
          end do                                                         11d28s23
         end do                                                          11d28s23
        end if                                                          12d5s23
        write(6,*)('done storing in ibufs ')
        ibcoff=itmp1                                                    11d28s23
       end do                                                           11d28s23
      end do                                                            11d28s23
      jbufs=0                                                           11d28s23
      do ip=1,mynprocg                                                  11d28s23
       ipt(ip)=jbufs                                                    11d28s23
       jbufs=jbufs+isnd(ip)                                             11d28s23
      end do                                                            11d28s23
      call dws_sync                                                     12d5s23
      call dws_all2allvb(bc(ibufs),isnd,ipt,bc(ibufr),ircv,ipf)         11d30s23
      jbufr=ibufr                                                       11d28s23
      do ip=1,mynprocg                                                  11d28s23
       ipf(ip)=jbufr                                                    11d28s23
       jbufr=jbufr+ircv(ip)                                             11d28s23
      end do                                                            11d28s23
      jbufs=ibufs                                                       11d28s23
c
c     under ipf we have
c     do isa=1,nsymb
c      do isb=1,isa                                                     11d28s23
c       do ij=1,nrowps(isa,isb)
c        do isc=1,nsymb
c         do il,ih
      do isa=1,nsymb                                                    11d28s23
       do isb=1,isa                                                     11d28s23
        isab=multh(isa,isb)                                             11d28s23
        if(nrowps(isa,isb).gt.0)then                                    12d1s23
         do isc=1,nsymb                                                 11d28s23
          isd=multh(isab,isc)                                           11d28s23
          itsc(isc)=ibcoff                                              11d28s23
          ibcoff=itsc(isc)+nrowps(isa,isb)*irefo(isc)*nvirt(isd)        12d1s23
         end do                                                         11d28s23
         call enough('trans3xden.tsc',bc,ibc)                           11d28s23
         do ij=0,nrowps(isa,isb)-1                                      12d1s23
          do isc=1,nsymb                                                  11d28s23
           isd=multh(isab,isc)                                            11d28s23
           if(min(irefo(isc),nvirt(isd)).gt.0)then                        11d28s23
            do ip=0,mynprocg-1                                            11d28s23
             ipp=ip+1                                                   11d28s23
             call ilimts(irefo(isc),nvirt(isd),mynprocg,ip,il,ih,i1s,   11d28s23
     $            i1e,i2s,i2e)                                          11d28s23
             i10=i1s                                                      11d28s23
             i1n=irefo(isc)                                               11d28s23
             do i2=i2s,i2e                                                11d28s23
              if(i2.eq.i2e)i1n=i1e                                        11d28s23
              do i1=i10,i1n                                               11d28s23
               iad=itsc(isc)+i2-1+nvirt(isd)*(i1-1+irefo(isc)*ij)       11d28s23
               bc(iad)=bc(ipf(ipp))                                     11d28s23
               ipf(ipp)=ipf(ipp)+1                                      11d28s23
              end do                                                      11d28s23
              i10=1                                                       11d28s23
             end do                                                       11d28s23
            end do                                                        11d28s23
           end if                                                         11d28s23
          end do                                                        11d30s23
         end do                                                          11d28s23
         do isc=1,nsymb                                                 11d28s23
          isd=multh(isc,isab)                                           11d28s23
          if(min(irefo(isc),nvirt(isd)).gt.0)then                       11d28s23
           nx=irefo(isc)*nrowps(isa,isb)                                12d1s23
           nbd=nbasisp(isd)*ncomp                                       11d28s23
           jorb=iorb(isd)+nbd*(irefo(isd)+idoubo(isd))                  11d28s23
           call dgemm('n','n',nbd,nx,nvirt(isd),1d0,bc(jorb),nbd,       11d28s23
     $          bc(itsc(isc)),nvirt(isd),0d0,bc(jbufs),nbd,             11d28s23
     $          'trans3xden.bufs')                                      11d28s23
           jbufs=jbufs+nbd*nx                                           11d28s23
          end if                                                        11d28s23
         end do                                                         11d28s23
        end if                                                          11d28s23
       end do                                                           11d28s23
      end do                                                            11d28s23
      ibcoff=jbufs                                                      11d28s23
      write(6,*)('words returned '),ibcoff-ibufs
      return
      end
