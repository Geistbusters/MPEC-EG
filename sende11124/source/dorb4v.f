c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dorb4v(ieraw,nsymbx,multh,nbasisp,ncsz0,nerimat,ncomp,  2d29s24
     $     iorb,i3xdao,itthalf,i3xd,irefo8,ipassa,npassa,ipassb,npassb, 4d5s24
     $     ja4,nszaa0,jb4,nszba0,isstor,ibstor,iptra,iptrb,nsza0,nszb0, 4d2s24
     $     srh,iflag,nocc,nvirtx,nduse,nbasdwsx,iapair,isnd,ircv,ipr,   4d15s24
     $     ipt,bc,ibc)                                                  4d15s24
      implicit real*8 (a-h,o-z)
      integer*8 isstor(*),ibstor(*),irefo8(*)                           4d5s24
      character*4 stype
      dimension ieraw(8,*),multh(8,*),nbasisp(*),iorb(*),naa(8),        4d5s24
     $     nbb(8),iptra(8),iptrb(8),iblock(8,8,8),iblockx(8,8),         4d10s24
     $    i4v(8,8,8),nvirtx(*),nduse(*),itmpab(8,8),nbasdwsx(*),        4d10s24
     $     noc(8),i3xdd(512),iuun(8),i3xddd(8,8,8),nan(8),nbn(8),       4d11s24
     $     i3xddp(8,8,8),j3xddp(8,8,8),nupdate(8,8,8),iapair(3,*),      4d15s24
     $     isnd(*),ircv(*),ipr(*),ipt(*),iuum(8),k3xddp(8,8,8)          4d16s24
      include "common.store"                                            2d29s24
      include "common.hf"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      4d15s24
     $     mynnode                                                      4d15s24
      data ifirst/0/
      data loopx/10000/
      save                                                              3d4s24
      loop=0
      if(ifirst.eq.0)then                                               3d19s24
       do isb=1,nsymb                                                   4d4s24
        noc(isb)=nbasdws(isb)-nvirt(isb)                                4d4s24
       end do                                                           4d4s24
       itop=ibcoff                                                      3d19s24
       do ipass=1,3
c
c     pass 1, count
c     pass 2, stuff into send buffer
c     pass 3, unpack from recv buffer
c
        if(ipass.eq.1)then                                              4d15s24
         do ip=1,mynprocg                                                4d15s24
          isnd(ip)=0d0                                                  4d15s24
          ircv(ip)=0d0                                                  4d15s24
         end do                                                          4d15s24
        end if                                                          4d15s24
        ipair=ipassa                                                     4d11s24
        nn=npassa
        ngaus=ipassb
        jbdat=nocc                                                       4d12s24
        if(ipass.eq.3)then                                              4d16s24
         i0=1+mynowprog                                                 4d16s24
         nskip=mynprocg                                                 4d17s24
        else                                                            4d16s24
         i0=1                                                           4d17s24
         nskip=1                                                        4d17s24
        end if                                                          4d16s24
        do i=i0,nn,nskip                                                4d17s24
         im=i-1                                                         4d15s24
         ip=mod(im,mynprocg)                                            4d15s24
         ipp=ip+1                                                       4d15s24
         jpair=ipair+2*(i-1)                                             4d11s24
         ja=jbdat+ibc(jpair)                                              2d19s10
         ja2=ja+ngaus                                                     2d19s10
         ja3=ja2+ngaus                                                    2d19s10
         jja4=ja3+ngaus                                                    2d19s10
         ja5=jja4+ngaus                                                    2d19s10
         ja6=ja5+ngaus                                                    2d19s10
         ja7=ja6+ngaus                                                    2d19s10
         ja8=ja7+ngaus                                                    5d11s10
         msza=2*ibc(ja)+1                                                 5d11s10
         msza0=msza
         msza=msza*ncomp                                                  8d24s15
         mszaa=msza                                                       5d12s10
         mszaa0=msza0                                                     5d4s20
         if(iapair(1,ibc(ja8)).gt.0)then                                  5d11s10
          mszaa=mszaa*2                                                   5d12s10
          mszaa0=mszaa0*2                                                 5d4s20
         end if                                                           5d11s10
         jb=jbdat+ibc(jpair+1)                                            2d19s10
         jb2=jb+ngaus                                                     2d19s10
         jb3=jb2+ngaus                                                    2d19s10
         jjb4=jb3+ngaus                                                    2d19s10
         jb5=jjb4+ngaus                                                    2d19s10
         jb6=jb5+ngaus                                                    2d19s10
         jb7=jb6+ngaus                                                    2d19s10
         jb8=jb7+ngaus                                                    5d11s10
         mszb=2*ibc(jb)+1                                                 5d11s10
         mszb0=mszb
         mszb=mszb*ncomp                                                  8d24s15
         mszba=mszb                                                       5d12s10
         mszba0=mszb0                                                     5d4s20
         if(iapair(1,ibc(jb8)).gt.0)then                                  5d11s10
          mszba=mszba*2                                                   5d12s10
          mszba0=mszba0*2                                                 5d4s20
         end if                                                           5d11s10
         do iss=1,nsymb                                                   4d10s24
          nan(iss)=0                                                      4d10s24
          nbn(iss)=0                                                      4d10s24
         end do                                                           4d9s24
         ipta=ibcoff                                                      4d9s24
         iptb=ipta+mszaa0                                                 4d9s24
         ibcoff=iptb+mszba0                                               4d9s24
         call enough('dorb.pta',bc,ibc)                                  4d11s24
         jpta=ipta-1                                                      4d9s24
         jptb=iptb-1                                                      4d9s24
         do ia=1,mszaa0                                                   4d9s24
          iaa=ia+ibc(jja4)                                                 4d9s24
          isa=isstor(iaa)                                                 4d9s24
          ibc(jpta+ia)=nan(isa)                                           4d9s24
          nan(isa)=nan(isa)+1                                             4d9s24
         end do                                                           4d9s24
         do ib=1,mszba0                                                   4d9s24
          ibb=ib+ibc(jjb4)                                                 4d9s24
          isb=isstor(ibb)                                                 4d9s24
          ibc(jptb+ib)=nbn(isb)                                           4d9s24
          nbn(isb)=nbn(isb)+1                                             4d9s24
         end do                                                           4d9s24
         do isa=1,nsymb                                                   4d10s24
          iptra(isa)=ibcoff                                               4d10s24
          iptrb(isa)=iptra(isa)+nan(isa)                                  4d10s24
          ibcoff=iptrb(isa)+nbn(isa)                                      4d10s24
         end do                                                           4d10s24
         do ia=1,mszaa0                                                   4d10s24
          iaa=ia+ibc(jja4)                                                 4d10s24
          isa=isstor(iaa)                                                 4d10s24
          ibc(iptra(isa)+ibc(jpta+ia))=ibstor(iaa)-1                      4d10s24
         end do                                                           4d10s24
         do ib=1,mszba0                                                   4d10s24
          ibb=ib+ibc(jjb4)                                                 4d10s24
          isb=isstor(ibb)                                                 4d10s24
          ibc(iptrb(isb)+ibc(jptb+ib))=ibstor(ibb)-1                     4d11s24
         end do                                                          4d11s24
         do is2=1,nsymb                                                  4d11s24
          do is1=1,nsymb                                                 4d11s24
           is12=multh(is1,is2)                                           4d12s24
           if(is1.eq.is2)then                                               4d10s24
            isw=0                                                           4d10s24
            nrow=(nbasisp(is1)*(nbasisp(is1)+1))/2                          4d10s24
           else                                                             4d10s24
            isw=1                                                           4d10s24
            nrow=nbasisp(is1)*nbasisp(is2)                                  4d10s24
           end if                                                           4d10s24
           ipt22=ibcoff                                                 4d16s24
           jpt22=ipt22                                                  4d16s24
           if(is1.ge.is2)then                                            4d11s24
            it1=is1                                                      4d12s24
            it2=is2                                                      4d12s24
            do i2=0,nbn(is2)-1                                           4d11s24
             i22=ibc(iptrb(is2)+i2)                                      4d10s24
             do i1=0,nan(is1)-1                                          4d11s24
              i11=ibc(iptra(is1)+i1)                                     4d10s24
              if(i22.le.i11.or.isw.eq.1)then                             4d11s24
               itri=((i11*(i11+1))/2)+i22                                4d11s24
               irec=i11+nbasisp(is1)*i22                                 4d11s24
               itri=itri+isw*(irec-itri)                                 4d11s24
               ibc(jpt22)=itri                                           4d11s24
               jpt22=jpt22+1                                             4d11s24
              end if                                                     4d11s24
             end do                                                      4d11s24
            end do
           else if(jja4.ne.jjb4)then                                       4d11s24
            it1=is2                                                      4d12s24
            it2=is1                                                      4d12s24
            do i2=0,nbn(is2)-1                                           4d11s24
             i22=ibc(iptrb(is2)+i2)                                      4d11s24
             do i1=0,nan(is1)-1                                          4d11s24
              i11=ibc(iptra(is1)+i1)                                     4d16s24
              irec=i22+nbasisp(is2)*i11                                  4d11s24
              ibc(jpt22)=irec                                            4d11s24
              jpt22=jpt22+1                                              4d11s24
             end do                                                      4d11s24
            end do                                                       4d11s24
           end if                                                        4d11s24
           npt22=jpt22-ipt22                                             4d11s24
           if(npt22.gt.0)then                                            4d11s24
            if(ipass.eq.1)then                                          4d15s24
             do is3=1,nsymb                                             4d15s24
              is4=multh(is12,is3)                                       4d15s24
              ncol=irefo8(is3)*nvirt(is4)                               4d15s24
              if(is1.ge.is2.or.jja4.ne.jjb4)then                        4d16s24
               npc=ncol*npt22*ncomp                                     4d16s24
               call ilimts(1,ncol,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,  4d15s24
     $             i2e)                                                 4d15s24
               nhere=ih+1-il                                             4d15s24
               if(nhere.gt.0)then                                        4d15s24
                isnd(ipp)=isnd(ipp)+nhere*ncomp*npt22                    4d15s24
               end if                                                    4d15s24
               if(ip.eq.mynowprog)then                                  4d16s24
                j3xddp(is1,is2,is3)=j3xddp(is1,is2,is3)+npc              4d16s24
                do ipppz=0,mynprocg-1                                   4d17s24
                 ipppzp=ipppz+1                                         4d16s24
                 call ilimts(1,ncol,mynprocg,ipppz,il,ih,i1s,i1e,i2s,   4d16s24
     $             i2e)                                                 4d15s24
                 nhere=ih+1-il                                             4d15s24
                 if(nhere.gt.0)then                                        4d15s24
                  ncp=nhere*ncomp*npt22                                   4d15s24
                  ircv(ipppzp)=ircv(ipppzp)+ncp                         4d16s24
                 end if                                                   4d15s24
                end do                                                  4d16s24
               end if                                                   4d16s24
              end if                                                    4d16s24
             end do                                                     4d15s24
            else if(ipass.eq.2)then                                     4d15s24
             do is3=1,nsymb                                             4d15s24
              is4=multh(is12,is3)                                       4d15s24
              ncol=irefo8(is3)*nvirt(is4)                               4d15s24
              call ilimts(1,ncol,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,  4d15s24
     $             i2e)                                                 4d15s24
              nhere=ih+1-il                                             4d15s24
              if(nhere.gt.0)then                                        4d15s24
               if(is1.ge.is2)then                                         4d12s24
                do icol=0,nhere*ncomp-1                                 4d15s24
                 do irow=0,npt22-1                                         4d11s24
                  irowx=ibc(ipt22+irow)                                    4d11s24
                  iad1=i3xddd(is1,is2,is3)+irowx+nrow*icol                 4d11s24
                  bc(ipt(ipp))=bc(iad1)                                 4d15s24
                  ipt(ipp)=ipt(ipp)+1                                   4d15s24
                 end do                                                    4d11s24
                end do                                                     4d11s24
               else if(jja4.ne.jjb4)then                                4d15s24
                do icol=0,nhere*ncomp-1                                 4d15s24
                 do irow=0,npt22-1                                         4d11s24
                  irowx=ibc(ipt22+irow)                                    4d11s24
                  iad1=i3xddd(is2,is1,is3)+irowx+nrow*icol                 4d11s24
                  bc(ipt(ipp))=bc(iad1)                                 4d15s24
                  ipt(ipp)=ipt(ipp)+1                                   4d15s24
                 end do                                                    4d11s24
                end do                                                     4d11s24
               end if                                                      4d11s24
              end if                                                    4d15s24
             end do                                                     4d15s24
            else if(ipass.eq.3)then                                     4d15s24
             do is3=1,nsymb                                               4d12s24
              nb=irefo8(is3)                                              4d11s24
              is4=multh(is12,is3)                                        4d12s24
              ncol=nb*nvirt(is4)                                        4d15s24
              ncolc=ncol*ncomp                                          4d15s24
              do ippz=0,mynprocg-1                                       4d17s24
               ippzp=ippz+1                                              4d17s24
               call ilimts(1,ncol,mynprocg,ippz,il,ih,i1s,i1e,i2s,i2e)   4d15s24
               nhere=ih+1-il                                            4d15s24
               if(nhere.gt.0)then                                       4d15s24
                il=il-1                                                 4d15s24
                if(is1.ge.is2)then                                         4d12s24
                 do ic=0,ncomp-1                                        4d15s24
                  do icol=0,nhere-1                                     4d15s24
                   icolp=icol+il+ncol*ic                                4d15s24
                   do irow=0,npt22-1                                    4d15s24
                    iad2=j3xddp(is1,is2,is3)+icolp+ncolc*irow           4d15s24
                    bc(iad2)=bc(ipr(ippzp))                              4d15s24
 1888               format(4i1)
 1889               format(a8,x,a4,es15.7,4i5,3i8)
                    ipr(ippzp)=ipr(ippzp)+1                               4d15s24
                   end do                                                    4d11s24
                  end do                                                     4d11s24
                 end do                                                 4d15s24
                else if(jja4.ne.jjb4)then                               4d16s24
                 do ic=0,ncomp-1                                        4d15s24
                  do icol=0,nhere-1                                     4d15s24
                   icolp=icol+il+ncol*ic                                4d15s24
                   do irow=0,npt22-1                                         4d11s24
                    iad2=j3xddp(is1,is2,is3)+icolp+ncolc*irow           4d16s24
                    bc(iad2)=bc(ipr(ippzp))
                    ipr(ippzp)=ipr(ippzp)+1                               4d15s24
                   end do                                                    4d11s24
                  end do                                                     4d11s24
                 end do                                                 4d15s24
                end if                                                  4d15s24
               end if                                                      4d11s24
              end do                                                       4d11s24
              j3xddp(is1,is2,is3)=j3xddp(is1,is2,is3)+ncolc*npt22       4d15s24
             end do                                                     4d17s24
            end if                                                      4d15s24
           end if                                                        4d11s24
          end do                                                         4d11s24
         end do                                                          4d11s24
         ibcoff=ipta                                                     4d11s24
        end do                                                           4d11s24
        if(ipass.eq.1)then                                              4d15s24
         nsendt=0                                                       4d15s24
         nrecvt=0                                                       4d15s24
         do ip=0,mynprocg-1                                             4d15s24
          ipp=ip+1                                                      4d15s24
          nsendt=nsendt+isnd(ipp)                                       4d15s24
          nrecvt=nrecvt+ircv(ipp)                                       4d15s24
         end do                                                         4d15s24
         bc(ibcoff)=dfloat(nsendt)                                      4d15s24
         bc(ibcoff+1)=dfloat(nrecvt)                                    4d15s24
         call dws_gsumf(bc(ibcoff),2)                                   4d15s24
         nst=nint(bc(ibcoff))                                           4d15s24
         nrt=nint(bc(ibcoff+1))                                         4d15s24
         if(nst.ne.nrt)then                                             4d15s24
          write(6,*)('total to send: '),nst                              4d15s24
          write(6,*)('total to receive: '),nrt                           4d15s24
          write(6,*)('not the same!!!')                                 4d17s24
          call dws_synca                                                4d15s24
          call dws_finalize                                             4d15s24
          stop 'dorb4v'                                                 4d15s24
         end if                                                         4d15s24
         if(isnd(mynowprog+1).ne.ircv(mynowprog+1))then
          write(6,*)('i''m not sending the same number of receiving '),
     $         ('from myself!!'),isnd(mynowprog+1),ircv(mynowprog+1)
          stop 'dorb4v'
         end if
         isbuff=ibcoff                                                  4d15s24
         irbuff=isbuff+max(nsendt,nrecvt)                               4d15s24
         ibcoff=irbuff+nrecvt                                           4d15s24
         call enough('dorb4v.srbuff',bc,ibc)                            4d15s24
         ibcoff=irbuff                                                  4d17s24
         ioff=i3xd
         do is=1,nsdlk1
          mm=nbasisp(isblk1(1,is))*nbasisp(isblk1(2,is))                6d25s24
          if(isblk1(1,is).eq.isblk1(2,is))then
           nrow=(nvirt(isblk1(1,is))*(nvirt(isblk1(1,is))+1))/2
           mrow=(nbasisp(isblk1(1,is))*(nbasisp(isblk1(2,is))+1))/2
          else
           nrow=nvirt(isblk1(1,is))*nvirt(isblk1(2,is))
           mrow=nbasisp(isblk1(1,is))*nbasisp(isblk1(2,is))
          end if
          ncol=irefo8(isblk1(3,is))*nvirt(isblk1(4,is))
          call ilimts(1,ncol,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)    4d15s24
          nhere=ih+1-il                                                   4d15s24
          i3xdd(is)=-1                                                    4d5s24
          i3xddd(isblk1(1,is),isblk1(2,is),isblk1(3,is))=-1               4d12s24
          if(min(nhere,nrow).gt.0)then                                    4d15s24
           mcol=nhere*ncomp                                               4d15s24
           i3xdd(is)=ibcoff                                               4d5s24
           i3xddd(isblk1(1,is),isblk1(2,is),isblk1(3,is))=ibcoff          4d10s24
           ibcoff=i3xdd(is)+mm*mcol                                     6d25s24
           is1=isblk1(1,is)
           is2=isblk1(2,is)
           is3=isblk1(3,is)
           is4=isblk1(4,is)
           nb=nhere                                                       4d15s24
           call abtrans3(bc(i3xdd(is)),bc(ioff),is1,is2,is3,is4,nb,       4d8s24
     $        noc,nbasisp,nvirt,ncomp,iorb,nrow,bc,ibc)                 4d8s24
           ibcoff=i3xdd(is)+mrow*mcol
           ioff=ioff+nrow*nhere                                           4d15s24
          end if
         end do
         ibcoff=max(ibcoff,irbuff+nrecvt)                               4d17s24
         ipt(1)=isbuff                                                  4d15s24
         ipr(1)=irbuff                                                  4d15s24
         do ip=1,mynprocg-1                                             4d15s24
          ipp=ip+1                                                      4d15s24
          ipt(ipp)=ipt(ip)+isnd(ip)                                     4d15s24
          ipr(ipp)=ipr(ip)+ircv(ip)                                     4d15s24
         end do                                                         4d15s24
        else if(ipass.eq.2)then                                         4d15s24
         ipt(1)=0                                                       4d15s24
         ipr(1)=0                                                       4d15s24
         do ip=1,mynprocg-1                                             4d15s24
          ipp=ip+1                                                      4d15s24
          ipt(ipp)=ipt(ip)+isnd(ip)                                     4d15s24
          ipr(ipp)=ipr(ip)+ircv(ip)                                     4d15s24
         end do                                                         4d15s24
         call dws_all2allvb(bc(isbuff),isnd,ipt,bc(irbuff),ircv,ipr)    4d15s24
         ipr(1)=irbuff                                                  4d15s24
         do ip=1,mynprocg-1                                             4d15s24
          ipp=ip+1                                                      4d15s24
          ipr(ipp)=ipr(ip)+ircv(ip)                                     4d15s24
         end do                                                         4d15s24
         jsbuff=isbuff                                                  4d16s24
         nscheck=0                                                      4d17s24
         do is3=1,nsymb                                                 4d15s24
          do is2=1,nsymb                                                4d15s24
           do is1=1,nsymb                                               4d15s24
            i3xddp(is1,is2,is3)=jsbuff                                  4d16s24
            jsbuff=jsbuff+j3xddp(is1,is2,is3)                           4d16s24
            nscheck=nscheck+j3xddp(is1,is2,is3)                         4d17s24
            j3xddp(is1,is2,is3)=i3xddp(is1,is2,is3)                     4d15s24
           end do                                                       4d15s24
          end do                                                        4d15s24
         end do                                                         4d15s24
         inewtop=jsbuff                                                 4d17s24
        else                                                            4d16s24
         ibcoff=inewtop                                                 4d17s24
        end if                                                          4d15s24
       end do                                                           4d15s24
       do isc=1,nsymb                                                   4d12s24
        do isb=1,nsymb                                                  4d12s24
         do isa=1,nsymb                                                 4d15s24
          j3xddp(isa,isb,isc)=i3xddp(isa,isb,isc)                       4d12s24
          nupdate(isa,isb,isc)=0                                        4d12s24
         end do                                                         4d12s24
        end do                                                          4d12s24
       end do                                                           4d12s24
       do isb=1,nsymb                                                   4d10s24
        irefo=irefo8(isb)
        nbpb=nbasisp(isb)*ncomp                                         4d10s24
        iuum(isb)=ibcoff                                                4d15s24
        ibcoff=iuum(isb)+irefo*nbpb                                     4d15s24
       end do                                                           4d10s24
       call enough('dorb4v.uun',bc,ibc)                                 4d10s24
       nuum=ibcoff-iuum(1)                                              4d16s24
       do iz=iuum(1),ibcoff-1                                           4d10s24
        bc(iz)=0d0                                                      4d10s24
       end do                                                           4d10s24
       igoal=2625922
       ifirst=1                                                         3d19s24
       return
      end if                                                            3d19s24
      if(iflag.eq.0)then                                                4d2s24
      if(npassa.eq.-1)then
       is1=ipassa                                                       4d10s24
       is2=ipassb                                                       4d10s24
       is12=multh(is1,is2)
       is3=npassb
       is4=multh(is12,is3)
c     isa and isb can be any thing.
c     for density, isa ge isb.
       nb=irefo8(is3)                                                   4d10s24
       nbpb=nbasisp(is3)*ncomp                                          4d10s24
       if(is1.eq.is2)then                                               4d10s24
        isw=0                                                           4d10s24
        nrow=(nbasisp(is1)*(nbasisp(is1)+1))/2                          4d10s24
       else                                                             4d10s24
        isw=1                                                           4d10s24
        nrow=nbasisp(is1)*nbasisp(is2)                                  4d10s24
       end if                                                           4d10s24
       in12=nocc/ncomp                                                    4d11s24
       in34=nocc-ncomp*in12                                               4d11s24
       icolx=0                                                          4d11s24
       ncolx=nb*nvirt(is4)*ncomp                                        4d11s24
       rms=0d0                                                          4d11s24
       nab0=nsza0*nszb0                                                 4d15s24
       if(is1.ge.is2)then                                               4d11s24
        do i2=0,nszb0-1                                                 4d11s24
         i22=ibc(iptrb(is2)+i2)                                         4d11s24
         do i1=0,nsza0-1                                                4d11s24
          i11=ibc(iptra(is1)+i1)                                        4d11s24
          if(i22.le.i11.or.isw.eq.1)then                                4d11s24
           do i4=0,nvirt(is4)-1                                         4d11s24
            i4p=i4+nvirt(is4)*in12                                      4d11s24
            i4pp=i4+noc(is4)                                            4d15s24
            irow=nb*i4p+j3xddp(is1,is2,is3)+ncolx*icolx                 4d11s24
            iadd4x0=ieraw(1,1)+i1+nsza0*(i2+nszb0*nbasisp(is3)*i4pp)    4d15s24
            do ia=0,nb-1                                                 4d11s24
             irowa=irow+ia                                              4d16s24
             juu=iuum(is3)+nbasisp(is3)*in34+nbpb*ia                    4d15s24
             iadd4x=iadd4x0                                             4d15s24
             do iv=0,nbasisp(is3)-1                                     4d15s24
              bc(juu+iv)=bc(juu+iv)+bc(irowa)*bc(iadd4x)                4d16s24
              iadd4x=iadd4x+nab0                                        4d15s24
             end do                                                     4d15s24
            end do                                                      4d11s24
           end do                                                       4d11s24
           icolx=icolx+1                                                4d11s24
          end if                                                        4d11s24
         end do                                                         4d11s24
        end do                                                          4d11s24
        nupdate(is1,is2,is3)=icolx                                      4d12s24
       else if(ja4.ne.jb4)then                                          4d11s24
        do i2=0,nszb0-1                                                 4d11s24
         do i1=0,nsza0-1                                                4d11s24
          do i4=0,nvirt(is4)-1                                          4d11s24
           i4p=i4+nvirt(is4)*in12                                       4d11s24
           i4pp=i4+noc(is4)                                             4d15s24
           irow=nb*i4p+j3xddp(is1,is2,is3)+ncolx*icolx                  4d16s24
           iadd4x0=ieraw(1,1)+i1+nsza0*(i2+nszb0*nbasisp(is3)*i4pp)     4d15s24
           do ia=0,nb-1                                                 4d11s24
            irowa=irow+ia                                               4d16s24
            juu=iuum(is3)+nbasisp(is3)*in34+nbpb*ia                     4d15s24
            iadd4x=iadd4x0                                              4d15s24
            do iv=0,nbasisp(is3)-1                                      4d15s24
             bc(juu+iv)=bc(juu+iv)+bc(irowa)*bc(iadd4x)                 4d16s24
             iadd4x=iadd4x+nab0                                         4d15s24
            end do                                                      4d15s24
           end do                                                       4d11s24
          end do                                                        4d11s24
          icolx=icolx+1                                                 4d11s24
         end do                                                         4d11s24
        end do                                                          4d11s24
        nupdate(is1,is2,is3)=icolx                                      4d16s24
       end if                                                           4d11s24
      else if(npassa.eq.-2)then                                         4d12s24
       is3=ipassb                                                       4d12s24
       is4=npassb                                                       4d12s24
       is34=multh(is3,is4)                                              4d12s24
       nb=irefo8(is3)                                                   4d10s24
       ncolx=nb*nvirt(is4)*ncomp                                        4d11s24
       do isb=1,nsymb                                                   4d12s24
        isa=multh(isb,is34)                                             4d12s24
         nadd=nupdate(isa,isb,is3)*ncolx                                4d12s24
         j3xddp(isa,isb,is3)=j3xddp(isa,isb,is3)+nadd                   4d12s24
         nupdate(isa,isb,is3)=0                                         4d12s24
       end do                                                           4d12s24
      end if                                                            3d4s24
      else
       jtthalf=itthalf                                                  4d17s24
       do isb=1,nsymb
         nbpb=nbasisp(isb)*ncomp                                        4d5s24
         irefo=irefo8(isb)
         if(min(nbpb,irefo).gt.0)then
         iuu=iuum(isb)
         itmp1=ibcoff                                                   4d5s24
         itmp2=itmp1+nbpb*nvirt(isb)                                    4d5s24
         ibcoff=itmp2+nvirt(isb)*irefo                                  4d5s24
         call enough('dorb4v.tmp12',bc,ibc)                             4d5s24
         do iv=0,nvirt(isb)-1
          ivp=iv+noc(isb)
          do ic=0,nbpb-1
           icv=iorb(isb)+ic+nbpb*ivp
           ivc=itmp1+iv+nvirt(isb)*ic
           bc(ivc)=bc(icv)
          end do
         end do
         call dgemm('n','n',nvirt(isb),irefo,nbpb,1d0,
     $        bc(itmp1),nvirt(isb),bc(iuu),nbpb,0d0,
     $        bc(itmp2),nvirt(isb),'dorb4v.tmp2')
         do icpy=0,nvirt(isb)*irefo-1                                   4d17s24
          bc(jtthalf+icpy)=bc(itmp2+icpy)                               4d17s24
         end do                                                         4d17s24
         jtthalf=jtthalf+nvirt(isb)*irefo                               4d17s24
         ibcoff=itmp1                                                   4d16s24
        end if
       end do
      end if
      return
      end
