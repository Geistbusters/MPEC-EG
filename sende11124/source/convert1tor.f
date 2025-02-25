c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine convert1tor(ihs,nff1,iff1,ncsf,isymmrci,nvirt,nroot,   7d15s21
     $     mdon,mdoop,nsymb,multh,nsing,vnew,nff10,ism,ismv,irel,irelv, 7d15s21
     $     iff10,norb,irefo,kff10,vxr,veci,ncsfti,nff0,iff0,bc,ibc)     11d14s22
c
c     for testing, take singles vector and turn it into a reference
c     vector
c
      implicit real*8 (a-h,o-z)
      integer*8 ihs(mdoop,*),i18,i28,i38,i48
      dimension nff1(mdoop,nsymb,2),iff1(*),ncsf(*),nvirt(*),multh(8,8),7d16s21
     $     vnew(nsing,*),nff10(mdoop,3),ism(*),ismv(*),irel(*),irelv(*),7d15s21
     $     irefo(*),iptb(200,8),kff10(*),vxr(*),veci(ncsfti,*),         7d19s21
     $     nff0(mdoop,3),iff0(*)                                        7d19s21
      data icall/0/
      save icall
      include "common.store"
      icall=icall+1
      write(6,*)('icall in convert1tor: '),icall
      write(6,*)('vxr in convert1tor: ')
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
      write(6,*)('next at end '),next
      ibcoff=iff10                                                      7d15s21
      ivoff=1                                                           7d15s21
      ifoff=0                                                           7d15s21
      jff10=1                                                           7d15s21
      jff1=1                                                            7d15s21
      i18=1
      ishit=0
      ntkeep=1
      vcuti=-0.8d0
      vcutiu=1.2d0
      vcut=-3d-2
      vcutu=3d2
      do ii=mdon+1,mdoop                                                7d16s21
       iarg=ii-mdon
       nff10(ii,1)=nff0(ii,1)                                           7d19s21
       nff10(ii,2)=jff10                                                7d16s21
       nff10(ii,3)=ivoff                                                7d16s21
       jf=nff0(ii,2)                                                    7d19s21
       jv=nff0(ii,3)                                                    7d19s21
       do if=1,nff0(ii,1)                                               7d19s21
        kff10(jff10)=iff0(jf)                                           7d19s21
        jff10=jff10+1                                                   7d19s21
        jf=jf+1                                                         7d19s21
        kff10(jff10)=iff0(jf)                                           7d19s21
        jff10=jff10+1                                                   7d19s21
        jf=jf+1                                                         7d19s21
        do i=1,ncsf(iarg)                                               7d19s21
         do ir=1,nroot                                                   7d19s21
          if(abs(veci(jv,ir)).ge.vcuti.and.abs(veci(jv,ir)).le.vcutiu)
     $         then
           write(6,*)('keep vint '),veci(jv,ir)
          else
           veci(jv,ir)=0d0
          end if
          vnew(ivoff,ir)=veci(jv,ir)                                      7d19s21
         end do                                                         7d19s21
         ivoff=ivoff+1                                                  7d19s21
         jv=jv+1                                                        7d19s21
        end do                                                          7d19s21
       end do
       do isb=1,nsymb                                                    7d15s21
        isbv=multh(isb,isymmrci)                                         7d15s21
        if(min(nff1(ii,isb,1),nvirt(isbv)).gt.0)then
         ncol=ncsf(iarg)*nff1(ii,isb,1)
         nrow=nroot*nvirt(isbv)
         ibcoff=jff10+nff1(ii,isb,1)*nvirt(isbv)+iff10-1
         itmp=ibcoff
         ibcoff=itmp+nrow*ncol
         call enough('convert1tor.  1',bc,ibc)
         i28=nrow
         i38=ncol
         call ddi_get(bc,ibc,ihs(ii,isb),i18,i28,i18,i38,bc(itmp))      11d15s22
         nff10(ii,1)=nff10(ii,1)+nff1(ii,isb,1)*nvirt(isbv)             7d16s21
         jff1=nff1(ii,isb,2)
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
             if(abs(bc(iad+ir)).gt.vcut.and.abs(bc(iad+ir)).lt.vcutu)
     $            then
              vnew(ivoff,ir)=bc(iad+ir)
              write(6,*)('keep '),bc(iad+ir),iv,if1,ii,i,ir
             else
              vnew(ivoff,ir)=0d0
             end if
            end do                                                      7d15s21
            ivoff=ivoff+1                                               7d15s21
           end do                                                       7d15s21
          end do                                                        7d15s21
          jff1=jff1p+1                                                  7d15s21
         end do
         call ddi_put(bc,ibc,ihs(ii,isb),i18,i28,i18,i38,bc(itmp))      11d15s22
         ibcoff=itmp
        end if
       end do                                                           7d15s21
      end do                                                            7d15s21
      write(6,*)('ntkeep='),ntkeep
      if(ntkeep.eq.0)then
       call dws_synca
       call dws_finalize
       stop
      end if
      write(6,*)('ivoff at end  = '),ivoff
      return
      end
