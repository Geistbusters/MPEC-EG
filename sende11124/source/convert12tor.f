c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine convert12tor(ihs,ihd,nff1,iff1,nff2,iff2,ncsf,ncsf2,   7d27s21
     $     isymmrci,nvirt,nroot,mdon,mdoop,nsymb,multh,nsing,vnew,nff10,7d27s21
     $     ism,ismv,irel,irelv,iff10,norb,irefo,kff10,vxr,veci,ncsfti,  7d27s21
     $     nff0,iff0,bc,ibc)                                            11d14s22
c
c     for testing, take singles vector and doubles vectors and turn it
c     into a reference vector
c
      implicit real*8 (a-h,o-z)
      integer*8 ihd(mdoop,*),ihs(mdoop,*),i18,i28,i38,i48
      dimension nff1(mdoop,nsymb,2),iff1(*),ncsf(*),nvirt(*),multh(8,8),7d16s21
     $     vnew(nsing,*),nff10(mdoop,3),ism(*),ismv(*),irel(*),irelv(*),7d15s21
     $     irefo(*),iptb(200,8),kff10(*),vxr(*),veci(ncsfti,*),         7d19s21
     $     nff0(mdoop,3),iff0(*),nff2(mdoop,nsymb,2),iff2(*),ncsf2(*)   7d27s21
      include "common.store"
      data loop/0/
      data loopx/100000/
      write(6,*)('vxr in convert12tor: '),loop
      call prntm2(vxr,1,nroot,1)
      do i=1,norb
       ismv(i)=ism(i)
       irelv(i)=irel(i)
      end do
      next=norb+1
      do isb=1,nsymb
       if(nvirt(isb).gt.200)stop 'convert1tor iptb'
       do iv=1,nvirt(isb)                                               7d15s21
        ismv(next)=isb                                                  7d15s21
        irelv(next)=iv+irefo(isb)                                       7d15s21
        iptb(iv,isb)=next                                               7d15s21
        next=next+1
       end do
      end do
      write(6,*)('next at end '),next,loop
      vdx=4.3032307190660778d-002
      vdxx=vdx*0.9d0
      do ii=mdon+1,mdoop                                                7d23s21
       nff10(ii,1)=0                                                    7d23s21
       nff10(ii,2)=-1
      end do                                                            7d23s21
      ibcoff=iff10                                                      7d15s21
      ivoff=1                                                           7d15s21
      ifoff=0                                                           7d15s21
      jff10=1                                                           7d15s21
      jff1=1                                                            7d15s21
      i18=1
      thresd=1d-2
      do ii=mdon+1,mdoop                                                7d16s21
       write(6,*)('for nclop = '),ii,loop
       iarg=ii-mdon
       write(6,*)('nff10(ii,2: '),nff10(ii,2),loop
       nff10(ii,1)=nff10(ii,1)+nff0(ii,1)                               7d27s21
       if(nff10(ii,2).lt.0)then
        nff10(ii,2)=jff10                                                7d16s21
        nff10(ii,3)=ivoff                                                7d16s21
       end if
       jf=nff0(ii,2)                                                    7d19s21
       jv=nff0(ii,3)                                                    7d19s21
       write(6,*)('copy internal part '),nff0(ii,1),loop
       do if=1,nff0(ii,1)                                               7d19s21
        kff10(jff10)=iff0(jf)                                           7d19s21
        jff10=jff10+1                                                   7d19s21
        jf=jf+1                                                         7d19s21
        kff10(jff10)=iff0(jf)                                           7d19s21
        jff10=jff10+1                                                   7d19s21
        jf=jf+1                                                         7d19s21
        do i=1,ncsf(iarg)                                               7d19s21
         do ir=1,nroot                                                   7d19s21
          vnew(ivoff,ir)=veci(jv,ir)                                      7d19s21
         end do                                                         7d19s21
         ivoff=ivoff+1                                                  7d19s21
         jv=jv+1                                                        7d19s21
        end do                                                          7d19s21
       end do
       write(6,*)('moving on to singles part '),loop
       do isb=1,nsymb                                                    7d15s21
        isbv=multh(isb,isymmrci)                                         7d15s21
        write(6,*)('for isb,isbv: '),isb,isbv,loop
        if(min(nff1(ii,isb,1),nvirt(isbv)).gt.0)then
         write(6,*)('adding '),loop
         ncol=ncsf(iarg)*nff1(ii,isb,1)
         nrow=nroot*nvirt(isbv)
         ibcoff=jff10+nff1(ii,isb,1)*nvirt(isbv)+iff10-1
         itmp=ibcoff
         ibcoff=itmp+nrow*ncol
         call enough('convert12tor.  1',bc,ibc)
         i28=nrow
         i38=ncol
         write(6,*)('getting '),ihs(ii,isb),nrow,ncol,loop
         call ddi_get(bc,ibc,ihs(ii,isb),i18,i28,i18,i38,bc(itmp))      11d15s22
         nff10(ii,1)=nff10(ii,1)+nff1(ii,isb,1)*nvirt(isbv)             7d16s21
         jff1=nff1(ii,isb,2)
         nkeep=0
         do if1=1,nff1(ii,isb,1)
          jff1p=jff1+1
          do iv=1,nvirt(isbv)
           kff10(jff10)=iff1(jff1)
           jff10=jff10+1
           kff10(jff10)=ibset(iff1(jff1p),iptb(iv,isbv))                7d15s21
           jff10=jff10+1                                                7d15s21
           do i=0,ncsf(iarg)-1                                          7d15s21
            iad=itmp+nrow*(i+ncsf(iarg)*(if1-1))-1+nroot*(iv-1)         7d16s21
            do ir=1,nroot                                               7d15s21
             vnew(ivoff,ir)=bc(iad+ir)
            end do                                                      7d15s21
            ivoff=ivoff+1                                               7d15s21
           end do                                                       7d15s21
          end do                                                        7d15s21
          jff1=jff1p+1                                                  7d15s21
         end do
         write(6,*)('done! '),loop
         if(nkeep.gt.0)then
          nnxn=nroot*nvirt(isbv)*ncsf(iarg)
          call prntm2(bc(itmp),nnxn,nff1(ii,isb,1),nnxn)
         end if
         write(6,*)('re-saving singles vector '),loop
         call ddi_put(bc,ibc,ihs(ii,isb),i18,i28,i18,i38,bc(itmp))      11d15s22
         ibcoff=itmp
        end if
       end do                                                           7d15s21
       write(6,*)('moving on to doubles '),loop
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
          ibcoff=jff10+nff2(ii,isb,1)*nn+iff10-1                         7d22s21
          itmp=ibcoff
          ibcoff=itmp+nrow*ncol
          call enough('convert12tor.  2',bc,ibc)
          i28=nrow
          i38=ncol
          call ddi_get(bc,ibc,ihd(ii,isb),i18,i28,i18,i38,bc(itmp))     11d15s22
          if(ipass.eq.1)then                                            7d23s21
           nff10(ii,1)=nff10(ii,1)+nff2(ii,isb,1)*nvnotv                  7d22s21
          else                                                          7d23s21
           nff10(ii+1,1)=nff10(ii+1,1)+nff2(ii,isb,1)*nvisv             7d23s21
           if(nff10(ii+1,2).lt.0)then                                   7d23s21
            nff10(ii+1,2)=jff10                                         7d23s21
            nff10(ii+1,3)=ivoff                                         7d23s21
           end if                                                       7d23s21
          end if                                                        7d23s21
          jff2=nff2(ii,isb,2)
          jtmp=itmp-1                                                    7d22s21
          do if1=1,nff2(ii,isb,1)
           jff2p=jff2+1
           do isbv1=1,nsymb                                              7d22s21
            isbv2=multh(isbv1,isbv12)                                    7d22s21
            if(isbv2.ge.isbv1)then                                       7d22s21
             if(isbv1.eq.isbv2)then                                      7d22s21
               do iv=1,nvirt(isbv1)                                       7d22s21
                if(ipass.eq.2)then
                 kff10(jff10)=ibset(iff2(jff2),iptb(iv,isbv1))             7d22s21
                 jff10=jff10+1
                 kff10(jff10)=iff2(jff2p)
                 jff10=jff10+1                                                7d15s21
                end if
                do i=0,ncsf2(iarg)-1                                      7d22s21
                 iad=jtmp+nroot*(i+ncsf2(iarg)*(iv-1))
                 do ir=1,nroot                                               7d15s21
                  if(ipass.eq.2)vnew(ivoff,ir)=bc(iad+ir)
                 end do                                                      7d15s21
                 if(ipass.eq.2)ivoff=ivoff+1                                               7d15s21
                end do                                                       7d15s21
               end do                                                        7d15s21
              jtmp=jtmp+nroot*nvirt(isbv1)*ncsf2(iarg)                   7d22s21
              nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                      7d22s21
              isw=0                                                      7d22s21
             else                                                        7d22s21
              isw=1                                                      7d22s21
              nvv=nvirt(isbv1)*nvirt(isbv2)                              7d22s21
             end if                                                      7d22s21
              do iv2=0,nvirt(isbv2)-1                                     7d22s21
               iv2p=iv2+1                                                 7d22s21
               itop=(iv2+isw*(nvirt(isbv1)-iv2))-1                        7d22s21
               do iv1=0,itop                                              7d22s21
                iv1p=iv1+1                                                7d22s21
                itri=((iv2*(iv2-1))/2)+iv1                                7d22s21
                irec=iv1+nvirt(isbv1)*iv2                                 7d22s21
                itri=itri+isw*(irec-itri)                                 7d22s21
                if(ipass.eq.1)then
                 kff10(jff10)=iff2(jff2)
                 jff10=jff10+1
                 kff10(jff10)=ibset(iff2(jff2p),iptb(iv1p,isbv1))            7d22s21
                 kff10(jff10)=ibset(kff10(jff10),iptb(iv2p,isbv2))            7d22s21
                 jff10=jff10+1                                                7d15s21
                end if
                do i=0,ncsf(iarg)-1                                         7d22s21
                 iad=jtmp+nroot*(i+ncsf(iarg)*itri)
                 do ir=1,nroot                                            7d22s21
                  if(ipass.eq.1)vnew(ivoff,ir)=bc(iad+ir)                               7d22s21
                 end do                                                   7d22s21
                 if(ipass.eq.1)ivoff=ivoff+1                                            7d22s21
                end do                                                    7d22s21
               end do                                                     7d22s21
              end do                                                      7d22s21
             jtmp=jtmp+nroot*nvv*ncsf(iarg)                              7d22s21
            end if                                                       7d22s21
           end do                                                        7d22s21
           jff2=jff2p+1                                                  7d15s21
          end do
          if(ipass.eq.1)then
          call ddi_put(bc,ibc,ihd(ii,isb),i18,i28,i18,i38,bc(itmp))     11d15s22
          end if
          ibcoff=itmp
         end if
        end do                                                           7d15s21
       end do                                                           7d23s21
      end do                                                            7d15s21
      write(6,*)('vdx = '),vdx
      write(6,*)('returning from convert12tor! '),loop
      return
      end
c mpec2.1 version zeta copyright u.s. government
      subroutine looploopx(loop,loopx)
      if(loop.ge.loopx)then
       call dws_synca
       call dws_finalize
       stop 'looploopx'                                                 1d19s23
      end if
      loop=loop+1
      return
      end
