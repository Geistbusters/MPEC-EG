c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine convert2tor(ihd,nff2,iff2,ncsf,ncsf2,isymmrci,nvirt,   7d22s21
     $     nroot,mdon,mdoop,nsymb,multh,nsing,vnew,nff20,ism,ismv,irel, 7d22s21
     $     irelv,iff20,norb,irefo,kff20,vxr,veci,ncsfti,nff0,iff0,bc,   11d14s22
     $     ibc)                                                         11d14s22
c
c     for testing, take singles vector and turn it into a reference
c     vector
c
      implicit real*8 (a-h,o-z)
      integer*8 ihd(mdoop,*),i18,i28,i38,i48
      dimension nff2(mdoop,nsymb,2),iff2(*),ncsf(*),nvirt(*),multh(8,8),7d16s21
     $     vnew(nsing,*),nff20(mdoop,3),ism(*),ismv(*),irel(*),irelv(*),7d15s21
     $     irefo(*),iptb(200,8),kff20(*),vxr(*),veci(ncsfti,*),ncsf2(*),7d22s21
     $     nff0(mdoop,3),iff0(*)                                        7d19s21
      include "common.store"
      write(6,*)('vxr in convert2tor: '),ibcoff,iff20
      call prntm2(vxr,1,nroot,1)
      do i=1,norb
       ismv(i)=ism(i)
       irelv(i)=irel(i)
      end do
      next=norb+1
      do isb=1,nsymb                                                    7d22s21
       if(nvirt(isb).gt.200)stop 'convert2tor iptb'
       do iv=1,nvirt(isb)                                               7d15s21
        ismv(next)=isb                                                  7d15s21
        irelv(next)=iv+irefo(isb)                                       7d15s21
        iptb(iv,isb)=next                                               7d15s21
        next=next+1
       end do
      end do
      write(6,*)('next at end '),next
      ibcoff=iff20                                                      7d15s21
      ivoff=1                                                           7d15s21
      ifoff=0                                                           7d15s21
      jff20=1                                                           7d15s21
      jff2=1                                                            7d15s21
      i18=1
      do ii=mdon+1,mdoop                                                7d23s21
       nff20(ii,1)=0                                                    7d23s21
       nff20(ii,2)=-1
      end do                                                            7d23s21
      icntr=0
      do ii=mdon+1,mdoop                                                7d16s21
       write(6,*)('for ii = '),ii
       write(6,*)(nff20(jj,1),jj=mdon+1,mdoop)
       iarg=ii-mdon
       nff20(ii,1)=nff20(ii,1)+nff0(ii,1)                               7d23s21
       write(6,*)('nff0: '),nff0(ii,1)
       if(nff20(ii,2).lt.0)then
        nff20(ii,2)=jff20                                                7d16s21
        nff20(ii,3)=ivoff                                                7d16s21
       end if
       jf=nff0(ii,2)                                                    7d19s21
       jv=nff0(ii,3)                                                    7d19s21
       write(6,*)('zero vd!!')
       do if=1,nff0(ii,1)                                               7d19s21
        kff20(jff20)=iff0(jf)                                           7d19s21
        jff20=jff20+1                                                   7d19s21
        jf=jf+1                                                         7d19s21
        kff20(jff20)=iff0(jf)                                           7d19s21
        jff20=jff20+1                                                   7d19s21
        jf=jf+1                                                         7d19s21
        do i=1,ncsf(iarg)                                               7d19s21
         do ir=1,nroot                                                   7d19s21
          vnew(ivoff,ir)=veci(jv,ir)                                      7d19s21
         end do                                                         7d19s21
         ivoff=ivoff+1                                                  7d19s21
         jv=jv+1                                                        7d19s21
        end do                                                          7d19s21
       end do
       do ipass=1,2                                                     7d23s21
        do isb=1,nsymb                                                    7d15s21
         isbv12=multh(isb,isymmrci)                                      7d22s21
         nvisv=0                                                         7d22s21
         nvnotv=0                                                        7d22s21
         do isbv1=1,nsymb                                                7d22s21
          isbv2=multh(isbv1,isbv12)                                      7d22s21
          if(isbv2.ge.isbv1)then                                         7d22s21
           if(isbv1.eq.isbv2)then                                        7d22s21
            nvisv=nvisv+nvirt(isbv1)                                     7d22s21
            nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        7d22s21
           else                                                          7d22s21
            nvv=nvirt(isbv1)*nvirt(isbv2)                                7d22s21
           end if                                                        7d22s21
           nvnotv=nvnotv+nvv                                             7d22s21
          end if                                                         7d22s21
         end do                                                          7d22s21
         nn=nvisv+nvnotv                                                 7d22s21
         if(min(nff2(ii,isb,1),nn).gt.0)then                             7d22s21
          ncol=nff2(ii,isb,1)                                            7d22s21
          nrow=nroot*(ncsf2(iarg)*nvisv+ncsf(iarg)*nvnotv)               7d22s21
          ibcoff=jff20+nff2(ii,isb,1)*nn+iff20-1                         7d22s21
          itmp=ibcoff
          ibcoff=itmp+nrow*ncol
          call enough('convert2tor.  1',bc,ibc)
          i28=nrow
          i38=ncol
          call ddi_get(bc,ibc,ihd(ii,isb),i18,i28,i18,i38,bc(itmp))     11d15s22
          if(ipass.eq.1)then                                            7d23s21
           nff20(ii,1)=nff20(ii,1)+nff2(ii,isb,1)*nvnotv                  7d22s21
          else                                                          7d23s21
           nff20(ii+1,1)=nff20(ii+1,1)+nff2(ii,isb,1)*nvisv             7d23s21
           if(nff20(ii+1,2).lt.0)then                                   7d23s21
            nff20(ii+1,2)=jff20                                         7d23s21
            nff20(ii+1,3)=ivoff                                         7d23s21
           end if                                                       7d23s21
          end if                                                        7d23s21
          write(6,*)('for '),ii,isb
          jff2=nff2(ii,isb,2)
          jtmp=itmp-1                                                    7d22s21
          do if1=1,nff2(ii,isb,1)
           jff2p=jff2+1
           do isbv1=1,nsymb                                              7d22s21
            isbv2=multh(isbv1,isbv12)                                    7d22s21
            if(isbv2.ge.isbv1)then                                       7d22s21
             if(isbv1.eq.isbv2)then                                      7d22s21
              if(ipass.eq.2)then                                        7d23s21
               do iv=1,nvirt(isbv1)                                       7d22s21
                kff20(jff20)=ibset(iff2(jff2),iptb(iv,isbv1))             7d22s21
                jff20=jff20+1
                kff20(jff20)=iff2(jff2p)
                jff20=jff20+1                                                7d15s21
                do i=0,ncsf2(iarg)-1                                      7d22s21
                 iad=jtmp+nroot*(i+ncsf2(iarg)*(iv-1))
                 do ir=1,nroot                                               7d15s21
                   bc(iad+ir)=0d0
                  vnew(ivoff,ir)=bc(iad+ir)
                 end do                                                      7d15s21
                 ivoff=ivoff+1                                               7d15s21
                end do                                                       7d15s21
               end do                                                        7d15s21
              end if                                                     7d23s21
              jtmp=jtmp+nroot*nvirt(isbv1)*ncsf2(iarg)                   7d22s21
              nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                      7d22s21
              isw=0                                                      7d22s21
             else                                                        7d22s21
              isw=1                                                      7d22s21
              nvv=nvirt(isbv1)*nvirt(isbv2)                              7d22s21
             end if                                                      7d22s21
             if(ipass.eq.1)then                                         7d23s21
              do iv2=0,nvirt(isbv2)-1                                     7d22s21
               iv2p=iv2+1                                                 7d22s21
               itop=(iv2+isw*(nvirt(isbv1)-iv2))-1                        7d22s21
               do iv1=0,itop                                              7d22s21
                iv1p=iv1+1                                                7d22s21
                itri=((iv2*(iv2-1))/2)+iv1                                7d22s21
                irec=iv1+nvirt(isbv1)*iv2                                 7d22s21
                itri=itri+isw*(irec-itri)                                 7d22s21
                kff20(jff20)=iff2(jff2)
                jff20=jff20+1
                kff20(jff20)=ibset(iff2(jff2p),iptb(iv1p,isbv1))            7d22s21
                kff20(jff20)=ibset(kff20(jff20),iptb(iv2p,isbv2))            7d22s21
                jff20=jff20+1                                                7d15s21
                do i=0,ncsf(iarg)-1                                         7d22s21
                 iad=jtmp+nroot*(i+ncsf(iarg)*itri)
                 do ir=1,nroot                                            7d22s21
                   bc(iad+ir)=0d0
                  vnew(ivoff,ir)=bc(iad+ir)                               7d22s21
                 end do                                                   7d22s21
                 ivoff=ivoff+1                                            7d22s21
                end do                                                    7d22s21
               end do                                                     7d22s21
              end do                                                      7d22s21
             end if                                                     7d23s21
             jtmp=jtmp+nroot*nvv*ncsf(iarg)                              7d22s21
            end if                                                       7d22s21
           end do                                                        7d22s21
           jff2=jff2p+1                                                  7d15s21
          end do
          call ddi_put(bc,ibc,ihd(ii,isb),i18,i28,i18,i38,bc(itmp))     11d15s22
          ibcoff=itmp
         end if
        end do                                                           7d15s21
       end do                                                           7d23s21
      end do                                                            7d15s21
      write(6,*)('ivoff at end  = '),ivoff,ibcoff,jff20,iff20
      return
      end
