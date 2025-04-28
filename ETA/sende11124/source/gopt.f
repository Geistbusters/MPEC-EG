c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gopt(igopt,idata,nunique,nder,bc,ibc)                  4d11s23
      implicit real*8 (a-h,o-z)                                         7d28s22
      character*1 cphas                                                 7d28s22
      character*4 olab(3)                                               4d15s22
      logical lprint,ldebug
      integer*1 idogrado1(4)                                            4d11s23
      equivalence (idogrado4,idogrado1)                                 4d11s23
      include "common.basis"                                            7d28s22
      include "common.input"                                            7d28s22
      include "common.store"                                            7d28s22
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      common/singcm/iuse,nff
      data olab/'d/dx','d/dy','d/dz'/                                   4d15s22
      save                                                              7d28s22
      idogrado4=idogrado                                                4d11s23
      if(idata.eq.0)then                                                7d28s22
       lprint=mynowprog.eq.0                                            3d14s23
       ldebug=.false.                                                   3d14s23
       if(lprint)then                                                   3d14s23
        write(6,*)(' ')
        write(6,*)('setting up search directions in gopt')
       end if
       iderfact=ibcoff                                                  3d22s23
       ibcoff=iderfact+natom                                            3d22s23
       idata=ibcoff                                                     7d28s22
       ider=ibcoff                                                      7d28s22
       natom3=natom*3                                                   3d13s23
       ibcoff=ider+natom*natom3                                         3d13s23
       call enough('gopt.  1',bc,ibc)
       do iz=ider,ibcoff-1                                              7d28s22
        bc(iz)=0d0                                                      7d28s22
       end do                                                           7d28s22
       nder=0                                                           7d28s22
       do ixyz=1,3                                                       7d28s22
        do ia=1,natom                                                    7d28s22
         if(iapair(1,ia).ge.0)then
          ipuse=ipropsym(ixyz)                                           2d3s16
          if(iapair(1,ia).eq.0)then                                      2d3s16
           npass=1                                                       2d3s16
          else                                                           2d3s16
           npass=2
          end if                                                         2d3s16
          do ipass=1,npass
           if(ipuse.eq.1)then                                            7d28s22
            jder=ider+natom*(ixyz-1+3*nder)-1                           7d28s22
            nder=nder+1                                                 7d28s22
            bc(jder+ia)=1d0                                             7d28s22
            if(iapair(1,ia).eq.0)then                                      2d3s16
             bc(iderfact+nder-1)=1d0                                    3d22s23
             if(lprint)then
              write(6,*)('atom '),ia
              write(6,*)('symmetrically unique ')
             end if
            else                                                           2d3s16
             bc(iderfact+nder-1)=0.5d0                                  3d22s23
             if(lprint)then                                             3d14s23
              write(6,*)('atom '),ia                                    3d14s23
              write(6,*)('symmetrically related to atom '),iapair(1,ia) 3d14s23
             end if                                                     3d14s23
            end if
            if(npass.eq.1)then
             cphas=' '                                                    2d3s16
            else if(ipass.eq.1)then                                       2d3s16
             cphas='+'                                                    2d3s16
             bc(jder+iapair(1,ia))=1d0                                  7d28s22
            else
             cphas='-'                                                    2d3s16
             bc(jder+iapair(1,ia))=-1d0                                 7d28s22
            end if                                                        2d3s16
            if(lprint)then                                              3d14s23
             if(iapair(1,ia).eq.0)then                                    4d28s22
              write(6,1)olab(ixyz),ia,ipuse                               4d28s22
             else                                                         4d28s22
              write(6,11)olab(ixyz),ia,cphas,iapair(1,ia),ipuse           4d28s22
   11        format('for operator ',a4,i2,a1,i2,' with symmetry ',i1)   7d28s22
             end if                                                       4d28s22
            end if                                                      3d14s23
    1       format('for operator ',a4,i2,' with symmetry ',i1)           4d28s22
           end if                                                        7d28s22
           if(npass.eq.2)ipuse=multh(ipuse,iapair(2,ia))                4d19s23
          end do
         end if
        end do                                                           7d28s22
       end do                                                            7d28s22
       if(ldebug)then
        write(6,*)('ider matrix ')                                       7d28s22
        call prntm2(bc(ider),natom3,nder,natom3)                         7d28s22
       end if
       if(lprint)write(6,*)('eliminate cm motion ')                     3d22s23
       itrans=ibcoff                                                    7d28s22
       ibcoff=itrans+nder*nder                                          7d28s22
       nunique=0                                                        7d28s22
       ipos=ibcoff                                                      7d28s22
       in2=ipos+nder                                                    7d28s22
       ibcoff=in2+nder                                                  7d28s22
       call enough('gopt.  2',bc,ibc)
       do ixyz=1,3                                                      7d28s22
        if(ldebug)write(6,*)('for ixyz = '),ixyz
        jder=ider+natom*(ixyz-1)                                        7d28s22
        n1=0                                                            7d28s22
        n2=0                                                            7d28s22
        do i=1,nder                                                     7d28s22
         if(ldebug)write(6,*)('for der '),i
         kder=jder+natom3*(i-1)                                         7d28s22
         if(ldebug)call prntm2(bc(kder),natom,1,natom)
         dot=0d0                                                        7d28s22
         sz=0d0                                                         7d28s22
         do ia=0,natom-1                                                7d28s22
          dot=dot+bc(kder+ia)*xmaxs(ia+1)                               3d22s23
          sz=sz+bc(kder+ia)**2                                          7d28s22
         end do                                                         7d28s22
         if(sz.gt.1d-10)then                                            7d28s22
          if(ldebug)write(6,*)('dot with translation: '),dot
          if(abs(dot).lt.1d-10)then                                     7d28s22
           ibc(ipos+i-1)=1                                              7d28s22
           n1=n1+1                                                      7d28s22
          else                                                          7d28s22
           ibc(ipos+i-1)=2                                              7d28s22
           ibc(in2+n2)=i-1                                              7d28s22
           n2=n2+1                                                      7d28s22
          end if                                                        7d28s22
         else                                                           7d28s22
          ibc(ipos+i-1)=0
         end if                                                         7d28s22
        end do                                                          7d28s22
        if(ldebug)                                                      3d14s23
     $    write(6,*)('what we have for ipos: '),(ibc(ipos+i),i=0,nder-1)3d14s23
        if(n1.gt.0)then                                                 7d28s22
         do i=0,nder-1                                                  7d28s22
          if(ibc(ipos+i).eq.1)then                                      7d28s22
           if(lprint)write(6,*)('keep der '),i+1,('as is')              3d14s23
           jtrans=itrans+nder*nunique                                   7d28s22
           do j=0,nder-1                                                7d28s22
            bc(jtrans+j)=0d0                                            7d28s22
           end do                                                       7d28s22
           bc(jtrans+i)=1d0                                             7d28s22
           nunique=nunique+1                                            7d28s22
          end if                                                        7d28s22
         end do                                                         7d28s22
        end if                                                          7d28s22
        if(n2.gt.0)then                                                 7d28s22
         write(6,*)('itrans so far ')
         call prntm2(bc(itrans),nder,nunique,nder)
         ifmat=ibcoff                                                   7d28s22
         icoef=ifmat+natom*n2                                           7d28s22
         iwgt=icoef+n2                                                  7d28s22
         idat=iwgt+natom                                                7d28s22
         iscr=idat+natom                                                7d28s22
         ibcoff=iscr+2*n2+n2*n2+2*n2*natom+natom                        7d28s22
         call enough('gopt.  3',bc,ibc)
         jfmat=ifmat                                                    7d28s22
         do i=0,n2-1                                                    7d28s22
          kder=jder+natom3*ibc(in2+i)                                   7d28s22
          do j=0,natom-1                                                7d28s22
           bc(jfmat+j)=bc(kder+j)                                       7d28s22
          end do                                                        7d28s22
          jfmat=jfmat+natom                                             7d28s22
         end do                                                         7d28s22
         do i=0,natom-1                                                 7d28s22
          bc(idat+i)=xmaxs(i+1)                                            3d22s23
          bc(iwgt+i)=1d0                                                7d28s22
         end do                                                         7d28s22
         iuse=0                                                         7d28s22
         if(ldebug)then                                                 3d14s23
          iout=6                                                        3d14s23
         else                                                           3d14s23
          iout=0                                                        3d14s23
         end if                                                         3d14s23
         call lsqfit2(bc(ifmat),natom,n2,bc(idat),natom,1,natom,        7d28s22
     $        bc(icoef),n2,bc(iscr),bc(iwgt),iout,rms,bc,ibc)           3d27s23
         if(rms.gt.1d-10.or.idogrado1(4).ne.2)then                      4d11s23
          if(lprint)                                                    3d14s23
     $         write(6,*)('we failed to fit translation, so keep as is')3d14s23
          do i=0,n2-1                                                   7d28s22
           jtrans=itrans+nder*nunique                                   7d28s22
           do j=0,nder-1                                                7d28s22
            bc(jtrans+j)=0d0                                            7d28s22
           end do                                                       7d28s22
           bc(jtrans+ibc(in2+i))=1d0
           nunique=nunique+1                                             7d28s22
          end do                                                        7d28s22
          ibcoff=ifmat                                                  7d28s22
         else                                                           7d28s22
          if(lprint)then                                                3d14s23
           write(6,*)('linear combination that is translation: ')
           write(6,*)('derivatives '),(ibc(in2+i),i=0,n2-1)             3d14s23
           call prntm2(bc(icoef),1,n2,1)                                 7d28s22
          end if                                                        3d14s23
          do i=0,n2-1                                                   7d28s22
           bc(ifmat+i)=bc(icoef+i)                                      7d28s22
          end do                                                        7d28s22
          iseed=12353299
          call randset(iseed)                                           7d28s22
          do j=1,n2-1                                                   7d28s22
           do k=0,n2-1                                                  7d28s22
            iad=ifmat+k+n2*j                                            7d28s22
            bc(iad)=2d0*(randdws()-0.5d0)                               7d28s22
           end do                                                       7d28s22
          end do                                                        7d28s22
          if(ldebug)then
           write(6,*)('trial vector matrix ')
           call prntm2(bc(ifmat),n2,n2,n2)
          end if                                                        3d14s23
          do i=0,n2-1                                                   7d28s22
           ivec=ifmat+n2*i                                              7d28s22
           do j=0,i-1                                                   7d28s22
            jvec=ifmat+n2*j                                             7d28s22
            dot=0d0                                                     7d28s22
            do k=0,n2-1                                                 7d28s22
             dot=dot+bc(ivec+k)*bc(jvec+k)                              7d28s22
            end do                                                      7d28s22
            do k=0,n2-1                                                 7d28s22
             bc(ivec+k)=bc(ivec+k)-dot*bc(jvec+k)                       7d28s22
            end do                                                      7d28s22
           end do                                                       7d28s22
           dot=0d0                                                      7d28s22
           do k=0,n2-1                                                  7d28s22
            dot=dot+bc(ivec+k)**2                                       7d28s22
           end do                                                       7d28s22
           dot=1d0/sqrt(dot)                                            7d28s22
           do k=0,n2-1                                                  7d28s22
            bc(ivec+k)=bc(ivec+k)*dot                                   7d28s22
           end do                                                       7d28s22
          end do                                                        7d28s22
          if(ldebug)then
           write(6,*)('orthogonalized vectors: ')
           call prntm2(bc(ifmat),n2,n2,n2)
          end if                                                        3d14s23
          jfmat=ifmat+n2                                                7d28s22
          do i=1,n2-1                                                   7d28s22
           jtrans=itrans+nder*nunique                                   7d28s22
           do j=0,nder-1                                                7d28s22
            bc(jtrans+j)=0d0                                            7d28s22
           end do                                                       7d28s22
           do k=0,n2-1                                                  7d28s22
            bc(jtrans+ibc(in2+k))=bc(jfmat+k)*bc(iderfact+ibc(in2+k))   3d23s23
           end do                                                       7d28s22
           jfmat=jfmat+n2                                               7d28s22
           nunique=nunique+1                                             7d28s22
          end do                                                        7d28s22
          ibcoff=ifmat
         end if                                                         7d28s22
        end if                                                          7d28s22
       end do                                                           7d28s22
       if(ldebug)then                                                   3d14s23
        write(6,*)('transformation matrix: '),itrans,ibcoff
        call prntm2(bc(itrans),nder,nunique,nder)                         7d28s22
       end if                                                           3d14s23
       if(lprint)then                                                   3d14s23
        write(6,*)('number of unique directions: '),nunique             3d14s23
        write(6,*)(' ')                                                 3d14s23
       end if                                                           3d14s23
      else                                                              7d28s22
      end if                                                            7d28s22
      return                                                            7d28s22
      end                                                               7d28s22
