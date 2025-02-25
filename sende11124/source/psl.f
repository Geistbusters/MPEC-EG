c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine psl(nps,ips,hdig,pham,numb,iaorb,nalpha,
     $                 iborb,nbeta,itmpa,itmpb,ixlzz,                   12d20s19
     $                 itmpc,itmpd,itmpe,itmpf,norb,                    12d22s19
     $                 nsymb,nsbeta,islz,multh,vecs,nspina,lambda,nlzz, 12d31s19
     $     lwrite,bc,ibc,ielecs)                                               11d9s22
      implicit real*8 (a-h,o-z)
      integer*2 ips(4,1)                                                8d22s06
      integer*1 iaorb,iborb,itmpc,itmpd,itmpe,itmpf                     11d14s05
      logical lprint,lwrite                                             12d31s19
      include "common.cas"
      include "common.store"                                            8d15s06
c
c     explicitly form p-space lz^2 matrix elements
c
      dimension hdig(1),pham(*),iaorb(nalpha,1),                        1d18s23
     $          iborb(nbeta,1),itmpa(1),itmpb(1),                       12d22s19
     $          itmp2(2,2),itmpc(1),itmpd(1),                           12d22s19
     $          itmpe(1),itmpf(1),ncsym(8),nsbeta(8),                   12d22s19
     $     ixlzz(8,*),multh(8,8),vecs(*),islz(*),igoalx(4)                        12d31s19
      data icall/0/
      data loopit/0/
      save icall,loopit                                                 12d22s19
      npass=nlzz/2                                                      12d31s19
      ione=npass+1                                                      12d31s19
      lprint=.false.
      icall=icall+1
      ism=ibcoff                                                        8d15s06
      irelo=ism+norb                                                    8d15s06
      ibcoff=irelo+norb                                                 8d15s06
      call enough('psl.  1',bc,ibc)
      jsm=ism-1                                                         8d15s06
      jrelo=irelo-1                                                     8d15s06
      do i=1,nsymb                                                      8d15s06
       do j=1,iacto(i)                                                  8d15s06
        ibc(jsm+j)=i                                                    8d15s06
        ibc(jrelo+j)=j                                                  8d15s06
  153   format(3i5)                                                     8d15s06
       end do                                                           8d15s06
       jsm=jsm+iacto(i)                                                 8d15s06
       jrelo=jrelo+iacto(i)                                             8d15s06
      end do                                                            8d15s06
      jsm=ism-1                                                         8d15s06
      jrelo=irelo-1                                                     8d15s06
      do i=1,nps
       ii=i+nps*(i-1)                                                   1d18s23
       pham(ii)=hdig(i)                                                 1d18s23
      end do
      do i=1,nps-1
       ipa=ips(1,i)                                                     8d22s06
       ipb=ips(2,i)                                                     8d22s06
       isbu=ips(3,i)                                                    8d22s06
       nasum2=0                                                         8d15s06
       do isb=1,isbu-1                                                  8d15s06
        nasum2=nasum2+nadet(isb)                                        8d15s06
       end do                                                           8d15s06
       nbsum2=0                                                         8d15s06
       do isb=1,nsbeta(isbu)-1                                          8d15s06
        nbsum2=nbsum2+nbdet(isb)                                        8d15s06
       end do                                                           8d15s06
       ia=ipa+nasum2                                                    8d15s06
       ib=ipb+nbsum2                                                    8d15s06
       do j=i+1,nps
        jpa=ips(1,j)                                                    8d22s06
        jpb=ips(2,j)                                                    8d22s06
        jsbu=ips(3,j)                                                   8d22s06
        nasum1=0                                                        8d15s06
        do isb=1,jsbu-1                                                 8d15s06
         nasum1=nasum1+nadet(isb)                                       8d15s06
        end do                                                          8d15s06
        nbsum1=0                                                        8d15s06
        do isb=1,nsbeta(jsbu)-1                                         8d15s06
         nbsum1=nbsum1+nbdet(isb)                                       8d15s06
        end do                                                          8d15s06
        ja=jpa+nasum1                                                   8d15s06
        jb=jpb+nbsum1                                                   8d15s06
        if(lprint)write(6,*)('ia,ib,ja,jb '),ia,ib,ja,jb
        isuma=0
        do i1=1,norb                                                    8d15s06
         itmpc(i1)=0
         itmpd(i1)=0
         itmpe(i1)=0
         itmpf(i1)=0
        end do
        do i1=1,nalpha
         itmpc(iaorb(i1,ja))=1
         itmpd(iaorb(i1,ia))=1
        end do
        do i1=1,norb                                                    8d15s06
         itmpa(i1)=itmpc(i1)-itmpd(i1)
         isuma=isuma+iabs(itmpa(i1))
        end do
        isumb=0
        do i1=1,nbeta
         itmpe(iborb(i1,jb))=1
         itmpf(iborb(i1,ib))=1
        end do
        do i1=1,norb                                                    8d15s06
         itmpb(i1)=itmpe(i1)-itmpf(i1)
         isumb=isumb+iabs(itmpb(i1))
        end do
        isum=isumb+isuma
         ans=0d0
         if(lprint)write(6,*)('no. difference = '),isum
         loopit=loopit+1
        if(isum.eq.2)then
         if(lprint)write(6,*)('single ')
         if(isuma.eq.0)then
          if(lprint)write(6,*)('beta single'),i,j
          do i1=1,norb                                                  8d15s06
           if(itmpb(i1).eq.1)nfrom=i1
           if(itmpb(i1).eq.-1)nto=i1
          end do
          nb=0
          do i1=min(nfrom,nto)+1,max(nfrom,nto)-1
           if(itmpe(i1).eq.1)nb=nb+1
          end do
          if(mod(nb,2).eq.0)then
           phs=1d0
          else
           phs=-1d0
          end if
          isbx=ibc(jsm+nto)                                             8d15s06
          isby=ibc(jsm+nfrom)                                           8d15s06
          if(isbx.ne.isby)then                                          8d15s06
           write(6,*)('single of wrong symmetry '),isbx,isby,nto,nfrom  8d15s06
           stop
          end if                                                        8d15s06
          into=ibc(jrelo+nto)                                           8d15s06
          infrom=ibc(jrelo+nfrom)                                       8d15s06
          if(lprint)write(6,*)('from to '),infrom,into,isbx             12d22s19
          iad1=ixlzz(isbx,ione)+into-1+iacto(isbx)*(infrom-1)           12d31s19
          part=bc(iad1)                                                 12d22s19
          if(lprint)write(6,*)('beta '),ib
          do i1=1,nbeta                                                 8d15s06
           isx=ibc(jsm+iborb(i1,ib))                                    8d15s06
           iisx=ibc(jrelo+iborb(i1,ib))                                 8d15s06
           do ipass=1,npass                                             12d31s19
            if(multh(isbx,isx).eq.islz(ipass))then                      12d31s19
             iad1=ixlzz(isbx,ipass)+into-1+iacto(isbx)*(iisx-1)         12d31s19
             iad2=ixlzz(isbx,ipass)+infrom-1+iacto(isbx)*(iisx-1)       12d31s19
             term=-2d0*bc(iad1)*bc(iad2)                                 12d22s19
             part=part-term                                              12d22s19
            end if
           end do                                                       12d31s19
          end do                                                        8d15s06
          ans=part*phs
         else
          do i1=1,norb                                                  8d15s06
           if(itmpa(i1).eq.1)nfrom=i1
           if(itmpa(i1).eq.-1)nto=i1
          end do
          nb=0
          do i1=min(nfrom,nto)+1,max(nfrom,nto)-1
           if(itmpc(i1).eq.1)nb=nb+1
          end do
          if(mod(nb,2).eq.0)then
           phs=1d0
          else
           phs=-1d0
          end if
          isbx=ibc(jsm+nto)                                             8d15s06
          isby=ibc(jsm+nfrom)                                           8d15s06
          if(isbx.ne.isby)then                                          8d15s06
           write(6,*)('symmetry of single does not match '),isbx,isby,
     $         nto,nfrom
           stop
          end if                                                        8d15s06
          into=ibc(jrelo+nto)                                           8d15s06
          infrom=ibc(jrelo+nfrom)                                       8d15s06
          iad1=ixlzz(isbx,ione)+into-1+iacto(isbx)*(infrom-1)           12d31s19
          part=bc(iad1)                                                 12d22s19
          do i1=1,nalpha                                                8d15s06
           isx=ibc(jsm+iaorb(i1,ia))                                    8d15s06
           iisx=ibc(jrelo+iaorb(i1,ia))                                 8d15s06
           do ipass=1,npass                                             12d31s19
            if(multh(isbx,isx).eq.islz(ipass))then                      12d31s19
             iad1=ixlzz(isbx,ipass)+into-1+iacto(isbx)*(iisx-1)         12d31s19
             iad2=ixlzz(isbx,ipass)+infrom-1+iacto(isbx)*(iisx-1)       12d31s19
             term=-2d0*bc(iad1)*bc(iad2)                                 12d22s19
             part=part-term                                              12d22s19
            end if                                                       12d22s19
           end do                                                       12d31s19
          end do                                                        8d15s06
          ans=part*phs
         end if
        else if(isum.eq.4)then
         if(lprint)write(6,*)('double ')
         if(isuma.eq.0)then
          if(lprint)write(6,*)('bb')
          ip=1
          im=1
          do i1=1,norb                                                  8d15s06
           if(itmpb(i1).eq.1)then
            itmp2(ip,1)=i1
            ip=ip+1
           end if
           if(itmpb(i1).eq.-1)then
            itmp2(im,2)=i1
            im=im+1
           end if
          end do
          call phsdet(itmpe,itmpb,itmp2,iphsb,norb)
          call phsdet(itmpf,itmpb,itmp2(1,2),iphsk,norb)
          iphs=iphsb+iphsk
          if(mod(iphs,2).eq.0)then
           phs=1d0
          else
           phs=-1d0
          end if
          isx1=ibc(jsm+itmp2(1,1))                                      8d15s06
          isx2=ibc(jsm+itmp2(1,2))                                      8d15s06
          isx3=ibc(jsm+itmp2(2,1))                                      8d15s06
          isx4=ibc(jsm+itmp2(2,2))                                      8d15s06
          iisx1=ibc(jrelo+itmp2(1,1))                                   8d15s06
          iisx2=ibc(jrelo+itmp2(1,2))                                   8d15s06
          iisx3=ibc(jrelo+itmp2(2,1))                                   8d15s06
          iisx4=ibc(jrelo+itmp2(2,2))                                   8d15s06
          val1=0d0                                                      12d22s19
          do ipass=1,npass                                              12d31s19
           if(multh(isx1,isx2).eq.islz(ipass))then                      12d31s19
            iad1=ixlzz(isx1,ipass)+iisx1-1+iacto(isx1)*(iisx2-1)        12d31s19
            iad2=ixlzz(isx3,ipass)+iisx3-1+iacto(isx3)*(iisx4-1)        12d31s19
            val1=val1-2d0*bc(iad1)*bc(iad2)                             12d31s19
           end if                                                        12d22s19
          end do                                                        12d31s19
          isx1=ibc(jsm+itmp2(2,1))                                      8d15s06
          isx2=ibc(jsm+itmp2(1,2))                                      8d15s06
          isx3=ibc(jsm+itmp2(1,1))                                      8d15s06
          isx4=ibc(jsm+itmp2(2,2))                                      8d15s06
          iisx1=ibc(jrelo+itmp2(2,1))                                   8d15s06
          iisx2=ibc(jrelo+itmp2(1,2))                                   8d15s06
          iisx3=ibc(jrelo+itmp2(1,1))                                   8d15s06
          iisx4=ibc(jrelo+itmp2(2,2))                                   8d15s06
          val2=0d0                                                      12d22s19
          do ipass=1,npass                                              12d31s19
           if(multh(isx1,isx2).eq.islz(ipass))then                      12d31s19
            iad1=ixlzz(isx1,ipass)+iisx1-1+iacto(isx1)*(iisx2-1)        12d31s19
            iad2=ixlzz(isx3,ipass)+iisx3-1+iacto(isx3)*(iisx4-1)        12d31s19
            val2=val2-2d0*bc(iad1)*bc(iad2)                             12d31s19
           end if                                                        12d22s19
          end do                                                        12d31s19
          ans=phs*(val1-val2)                                           12d22s19
         else if(isumb.eq.0)then
          if(lprint)write(6,*)('aa')
          ip=1
          im=1
          do i1=1,norb                                                  8d15s06
           if(itmpa(i1).eq.1)then
            itmp2(ip,1)=i1
            ip=ip+1
           end if
           if(itmpa(i1).eq.-1)then
            itmp2(im,2)=i1
            im=im+1
           end if
          end do
          call phsdet(itmpc,itmpa,itmp2,iphsb,norb)                     8d15s06
          call phsdet(itmpd,itmpa,itmp2(1,2),iphsk,norb)                8d15s06
          iphs=iphsb+iphsk
          if(mod(iphs,2).eq.0)then
           phs=1d0
          else
           phs=-1d0
          end if
          if(lprint)write(6,*)('phases '),iphsb,iphsk,iphs,phs
          isx1=ibc(jsm+itmp2(1,1))                                      8d15s06
          isx2=ibc(jsm+itmp2(1,2))                                      8d15s06
          isx3=ibc(jsm+itmp2(2,1))                                      8d15s06
          isx4=ibc(jsm+itmp2(2,2))                                      8d15s06
          iisx1=ibc(jrelo+itmp2(1,1))                                   8d15s06
          iisx2=ibc(jrelo+itmp2(1,2))                                   8d15s06
          iisx3=ibc(jrelo+itmp2(2,1))                                   8d15s06
          iisx4=ibc(jrelo+itmp2(2,2))                                   8d15s06
          val1=0d0                                                      12d22s19
          do ipass=1,npass                                              12d31s19
           if(multh(isx1,isx2).eq.islz(ipass))then                      12d31s19
            iad1=ixlzz(isx1,ipass)+iisx1-1+iacto(isx1)*(iisx2-1)        12d31s19
            iad3=ixlzz(isx3,ipass)+iisx3-1+iacto(isx3)*(iisx4-1)        12d31s19
            val1=val1-2d0*bc(iad1)*bc(iad3)                             12d31s19
           end if                                                        12d22s19
          end do                                                        12d31s19
          isx1=ibc(jsm+itmp2(2,1))                                      8d15s06
          isx2=ibc(jsm+itmp2(1,2))                                      8d15s06
          isx3=ibc(jsm+itmp2(1,1))                                      8d15s06
          isx4=ibc(jsm+itmp2(2,2))                                      8d15s06
          iisx1=ibc(jrelo+itmp2(2,1))                                   8d15s06
          iisx2=ibc(jrelo+itmp2(1,2))                                   8d15s06
          iisx3=ibc(jrelo+itmp2(1,1))                                   8d15s06
          iisx4=ibc(jrelo+itmp2(2,2))                                   8d15s06
          val2=0d0                                                      12d22s19
          do ipass=1,npass                                              12d31s19
           if(multh(isx1,isx2).eq.islz(ipass))then                      12d31s19
            iad1=ixlzz(isx1,ipass)+iisx1-1+iacto(isx1)*(iisx2-1)        12d31s19
            iad3=ixlzz(isx3,ipass)+iisx3-1+iacto(isx3)*(iisx4-1)        12d31s19
            val2=val2-2d0*bc(iad1)*bc(iad3)                             12d31s19
           end if                                                        12d22s19
          end do                                                        12d31s19
          ans=phs*(val1-val2)                                           12d22s19
          if(lprint)write(6,*)('integrals '),val1,val2
         else if(isuma.eq.2.and.isumb.eq.2)then
          if(lprint)write(6,*)('ab')
          do i1=1,norb                                                  8d15s06
           if(itmpa(i1).eq.1)nfroma=i1
           if(itmpa(i1).eq.-1)ntoa=i1
          end do
          nb=0
          do i1=min(nfroma,ntoa)+1,max(nfroma,ntoa)-1
           if(itmpc(i1).eq.1)nb=nb+1
          end do
          do i1=1,norb                                                  8d15s06
           if(itmpb(i1).eq.1)nfromb=i1
           if(itmpb(i1).eq.-1)ntob=i1
          end do
          do i1=min(nfromb,ntob)+1,max(nfromb,ntob)-1
           if(itmpe(i1).eq.1)nb=nb+1
          end do
          if(mod(nb,2).eq.0)then
           phs=1d0
          else
           phs=-1d0
          end if
          isfa=ibc(jsm+nfroma)                                          8d15s06
          ista=ibc(jsm+ntoa)                                            8d15s06
          isfb=ibc(jsm+nfromb)                                          8d15s06
          istb=ibc(jsm+ntob)                                            8d15s06
          infroma=ibc(jrelo+nfroma)                                     8d15s06
          intoa=ibc(jrelo+ntoa)                                         8d15s06
          infromb=ibc(jrelo+nfromb)                                     8d15s06
          intob=ibc(jrelo+ntob)                                         8d15s06
          term=0d0                                                      12d22s19
          do ipass=1,npass                                              12d31s19
           if(multh(isfa,ista).eq.islz(ipass))then                      12d31s19
            iada=ixlzz(isfa,ipass)+infroma-1+iacto(isfa)*(intoa-1)      12d31s19
            iadb=ixlzz(isfb,ipass)+infromb-1+iacto(isfb)*(intob-1)      12d31s19
            term=term-2d0*bc(iada)*bc(iadb)                             12d31s19
           end if                                                        12d22s19
          end do                                                        12d31s19
          ans=phs*term                                                  12d22s19
         end if
        end if
        ij=i+nps*(j-1)                                                  1d18s23
        pham(ij)=ans                                                    1d18s23
        ji=j+nps*(i-1)                                                  1d18s23
        pham(ji)=ans                                                    1d18s23
       end do
      end do
      ibcoff=ism                                                        2d16s07
      itmp=ibcoff                                                       12d31s19
      ibcoff=itmp+nps*nspina                                            12d31s19
      call enough('psl.  2',bc,ibc)
      if(lprint)then
       write(6,*)('L^2 matrix in ps basis ')
       call prntm2(pham,nps,nps,nps)
       ieig=ibcoff
       ivec=ieig+nps
       isym=ivec+nps*nps
       ibcoff=isym+nps
       call diagx(nps,pham,bc(ieig),bc(ivec),ibc(isym),bc,ibc)
       write(6,*)('eigenvalues in ps basis ')
       call prntm2(bc(ieig),1,nps,1)
       write(6,*)('eigenvectors in ps basis')
       call prntm2(bc(ivec),nps,nps,nps)
       ibcoff=ieig
      end if
      nrow=nps                                                          12d31s19
      do ipass=1,2                                                      12d31s19
       call dgemm('n','n',nrow,nspina,nps,1d0,pham,nrow,vecs,nps,0d0,   12d31s19
     $      bc(itmp),nrow,                                              12d31s19
     d' psl.  1')
       do i=0,nspina-1                                                  12d31s19
        do j=0,nrow-1                                                   12d31s19
         ji=itmp+j+nrow*i                                               12d31s19
         ij=i+1+nspina*j                                                12d31s19
         pham(ij)=bc(ji)                                                1d18s23
        end do                                                          12d31s19
       end do                                                           12d31s19
       nrow=nspina                                                      12d31s19
      end do                                                            12d31s19
      if(lprint)then
       write(6,*)('L^2 matrix in ps spin basis ')
       call prntm2(pham,nrow,nrow,nrow)
      end if
      ibcoff=itmp                                                       12d31s19
      ieig=ibcoff
      ivec=ieig+nspina                                                  12d31s19
      isym=ivec+nspina*nspina                                           12d31s19
      ibcoff=isym+nspina                                                12d31s19
      call enough('psl.  3',bc,ibc)
      call diagx(nspina,pham,bc(ieig),bc(ivec),ibc(isym),bc,ibc)        11d14s22
      if(lwrite)then                                                    1s1s20
       if(npass.eq.1)then                                               1s1s20
        write(6,*)('looking for lambda = '),lambda                      12d31s19
       else                                                             1s1s20
        write(6,*)('looking for angular momentum = '),lambda            1s1s20
       end if                                                           1s1s20
      end if                                                            1s1s20
      if(lprint)then
       write(6,*)('want: '),lambda,npass
       write(6,*)('got: ')
       call prntm2(bc(ieig),1,nspina,1)
      end if
      nkeep=0                                                           12d31s19
      nbad=0                                                            12d31s19
      do irt=1,nspina                                                   12d31s19
       ir=irt-1
       jvec=ivec+nspina*ir-1                                            12d31s19
       if(npass.eq.1)then                                               12d31s19
        try=sqrt(abs(bc(ieig+ir)))                                        12d24s19
        itry=nint(try)                                                   12d24s19
        delta=abs(bc(ieig+ir)-dfloat(itry*itry))                          12d24s19
       else                                                             12d31s19
        try=sqrt(abs(bc(ieig+ir))+0.25d0)-0.5d0                         12d31s19
        itry=nint(try)                                                   12d24s19
        delta=abs(bc(ieig+ir)-dfloat(itry*(itry+1)))                    12d31s19
       end if                                                           12d31s19
       if(delta.gt.1d-10)then                                           8d4s22
        write(6,*)('bad eigenvalue: '),bc(ieig+ir),delta                8d4s22
        write(6,*)('vector ')                                           8d4s22
        call prntm2(bc(jvec+1),nspina,1,nspina)                         8d4s22
        nbad=nbad+1                                                     8d4s22
       end if                                                           8d4s22
       if(delta.gt.1d-10.or.itry.eq.lambda)then                         4d17s21
        kvec=ivec-1+nspina*nkeep                                        12d31s19
        do j=1,nspina                                                   12d31s19
         bc(kvec+j)=bc(jvec+j)                                          12d31s19
        end do                                                          12d31s19
        nkeep=nkeep+1                                                   12d31s19
       end if                                                           12d31s19
      end do
      if(lwrite)then                                                    12d31s19
       write(6,*)('number of keepers = '),nkeep                         12d31s19
       if(nbad.gt.0)then                                                12d31s19
        write(6,*)('of which '),nbad,(' were not pure')                 12d31s19
       end if                                                           12d31s19
       end if                                                           12d31s19
       ratio=dfloat(nbad)/dfloat(nkeep)                                 4d25s21
       xnkeep=dfloat(nkeep)                                             3d4s21
       bc(ibcoff)=xnkeep                                                  4d25s21
       bc(ibcoff+1)=ratio                                               4d25s21
       call dws_bcast(bc(ibcoff),2)                                     4d25s21
       xnkeep=bc(ibcoff)                                                4d25s21
       ratio=bc(ibcoff+1)                                               4d25s21
       nkeep=nint(xnkeep)                                               3d4s21
       if(ratio.gt.0.d0)then                                            8d3s22
        if(lwrite)then                                                  4d25s21
         write(6,*)('there are too many unpure eigenvalues')            4d25s21
         write(6,*)('increase the number of p-space functions to fix '),4d25s21
     $        ('this')                                                  4d25s21
         stop
        end if
       end if                                                           4d25s21
       call dws_bcast(bc(ivec),nspina*nkeep)                            3d4s21
       itmp=ibcoff                                                      12d31s19
       ibcoff=itmp+nps*nkeep                                            12d31s19
       call enough('psl.  4',bc,ibc)
       call dgemm('n','n',nps,nkeep,nspina,1d0,vecs,nps,bc(ivec),nspina,12d31s19
     $      0d0,bc(itmp),nps,                                           12d31s19
     d' psl.  2')
       jtmp=itmp-1                                                      12d31s19
       do i=1,nps*nkeep                                                 12d31s19
        vecs(i)=bc(jtmp+i)                                              12d31s19
       end do                                                           12d31s19
       nspina=nkeep                                                     12d31s19
       ibcoff=ieig                                                      12d31s19
      ibcoff=ieig
      return
      end
