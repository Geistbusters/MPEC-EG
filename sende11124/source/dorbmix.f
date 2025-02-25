c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dorbmix(nsymb,idoubo,irefo,noc,nvirt,nsdlk,isblk,      6d24s24
     $     nsdlk1,isblk1,iden,nh0av,ioooo,ionex,bc,ibc,itt,nbasdws,     6d17s24
     $     jmats,kmats,nsdlkk,isblkk,i3x)                               10d25s23
      implicit real*8 (a-h,o-z)                                         7d20s23
      dimension idoubo(*),irefo(*),noc(*),nvirt(*),                     6d24s24
     $     isblk(4,*),ioooo(*),iden(*),nh0av(*),isblk1(4,*),            6d17s24
     $     ionex(*),itt(*),nbasdws(*),jmats(*),kmats(*),i3x(*),         10d25s23
     $     isblkk(4,*)                                                  10d25s23
      include "common.store"                                            7d20s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data loopx/11000/
      loop=0
      igoal=1
      igoul=itt(1)+7*3
      do is=1,nsdlk                                                     7d19s23
       if(isblk(1,is).eq.isblk(2,is))then                               10d26s23
        isb=isblk(1,is)
        nrow=(noc(isb)*(noc(isb)+1))/2                                  7d19s23
        jsb=isblk(3,is)                                                 7d19s23
        if(idoubo(isb).gt.0)then                                        10d25s23
         call ilimts(nvirt(jsb),nvirt(jsb),mynprocg,mynowprog,il,ih,     10d25s23
     $       i1s,i1e,i2s,i2e)                                           10d25s23
         i1n=nvirt(jsb)                                                 10d25s23
         i10=i1s                                                        10d25s23
         jj=jmats(is)                                                   10d25s23
         do i2=i2s,i2e                                                  10d25s23
          i2p=i2-1+irefo(jsb)                                           10d25s23
          if(i2.eq.i2e)i1n=i1e                                          10d25s23
          do i1=i10,i1n                                                 10d25s23
           i1p=i1-1+irefo(jsb)                                          10d25s23
           iadd=iden(jsb)+i1p+nh0av(jsb)*i2p                            10d25s23
           ff=4d0*bc(iadd)                                              10d25s23
           do ie=0,noc(isb)-1                                           10d25s23
            jtt=itt(isb)+nbasdws(isb)*ie                                10d25s23
            do id=0,idoubo(isb)-1                                       10d25s23
             ix=max(id,ie)                                              10d25s23
             in=min(id,ie)                                              10d25s23
             itri=((ix*(ix+1))/2)+in+jj                                 10d25s23
             bc(jtt+id)=bc(jtt+id)+ff*bc(itri)                          10d25s23
            end do                                                      10d25s23
           end do                                                       10d25s23
           jj=jj+nrow                                                   10d25s23
          end do                                                        10d25s23
          i10=1                                                         10d25s23
         end do                                                         10d25s23
         nhere2=i2e+1-i2s                                               10d25s23
         if(nhere2.gt.0)then                                            10d26s23
          itmp=ibcoff                                                    10d25s23
          ibcoff=itmp+nvirt(jsb)*nhere2                                  10d25s23
          call enough('dorbmix.goulJb',bc,ibc)                           10d25s23
          do iz=itmp,ibcoff-1                                            10d25s23
           bc(iz)=0d0                                                    10d25s23
          end do                                                         10d25s23
          i10=i1s                                                        10d25s23
          i1n=nvirt(jsb)                                                 10d25s23
          jj=jmats(is)-1                                                 10d25s23
          do i2=i2s,i2e                                                  10d25s23
           jtmp=itmp+nvirt(jsb)*(i2-i2s)-1                               10d25s23
           if(i2.eq.i2e)i1n=i1e                                          10d25s23
           do i1=i10,i1n                                                 10d25s23
            do id=1,idoubo(isb)                                          10d25s23
             iad=jj+((id*(id+1))/2)                                      10d25s23
             bc(jtmp+i1)=bc(jtmp+i1)+bc(iad)                             10d25s23
            end do                                                       10d25s23
            jj=jj+nrow                                                   10d25s23
           end do                                                        10d25s23
           i10=1                                                         10d25s23
          end do                                                         10d25s23
          jden=iden(jsb)+nh0av(jsb)*irefo(jsb)                           10d25s23
          jtt=itt(jsb)+idoubo(jsb)+nbasdws(jsb)*(noc(jsb)+i2s-1)         10d25s23
          call dgemm('n','n',nh0av(jsb),nhere2,nvirt(jsb),4d0,           10d25s23
     $         bc(jden),nh0av(jsb),bc(itmp),nvirt(jsb),1d0,              10d26s23
     $         bc(jtt),nbasdws(jsb),'dorbmix.goulJb2')                   10d26s23
          ibcoff=itmp                                                    10d26s23
         end if                                                         10d26s23
        end if                                                          10d25s23
c
        call ilimts(noc(jsb),noc(jsb),mynprocg,mynowprog,il,ih,i1s,     7d20s23
     $        i1e,i2s,i2e)                                              7d19s23
        if(min(irefo(jsb),idoubo(isb)).gt.0)then                        7d20s23
         nhere=ih+1-il                                                  7d19s23
         itmp=ibcoff                                                    7d20s23
         ibcoff=itmp+noc(jsb)*irefo(jsb)                                7d20s23
         call enough('dorbmix.1',bc,ibc)                                7d20s23
         do iz=itmp,ibcoff-1                                            7d20s23
          bc(iz)=0d0                                                    7d20s23
         end do                                                         7d20s23
         jtmp=itmp                                                      7d20s23
         do i2=0,irefo(jsb)-1                                           7d20s23
          i2p=i2+idoubo(jsb)                                            7d20s23
          do i1=0,noc(jsb)-1                                            7d20s23
           icol=i1+1+noc(jsb)*i2p                                       7d20s23
           if(icol.ge.il.and.icol.le.ih)then                            7d20s23
            icol0=icol
            icol=ioooo(is)+nrow*(icol-il)-1                             7d20s23
            sum=0d0                                                     7d20s23
            do i34=1,idoubo(isb)                                        7d20s23
             iad=icol+((i34*(i34+1))/2)                                 7d20s23
             sum=sum+bc(iad)                                            7d20s23
            end do                                                      7d20s23
            jtmp=itmp+i2+irefo(jsb)*i1
            bc(jtmp)=sum                                                10d26s23
           end if                                                       7d20s23
          end do                                                        7d20s23
          jtmp=jtmp+noc(jsb)                                            7d20s23
         end do                                                         7d20s23
         jtt=itt(jsb)+idoubo(jsb)                                       10d26s23
         call dgemm('n','n',nh0av(jsb),noc(jsb),irefo(jsb),4d0,         10d26s23
     $        bc(iden(jsb)),nh0av(jsb),bc(itmp),irefo(jsb),1d0,            7d20s23
     $        bc(jtt),nbasdws(jsb),'dorbmix.b')                           7d20s23
         ibcoff=itmp                                                    7d20s23
        end if                                                          7d20s23
        do i2=0,idoubo(jsb)-1                                             7d19s23
         do ia=0,noc(jsb)-1                                             7d20s23
          icol=ia+1+noc(jsb)*i2                                         7d20s23
          if(icol.ge.il.and.icol.le.ih)then                             7d19s23
           icol0=icol
           icol=ioooo(is)+nrow*(icol-il)                                7d19s23
           sum=0d0                                                      7d19s23
           iadd=iden(isb)                                               7d19s23
           do i4=0,irefo(isb)-1                                         7d19s23
            i4p=i4+idoubo(isb)                                          7d19s23
            do i3=0,irefo(isb)-1                                        7d19s23
             i3p=i3+idoubo(isb)                                         7d19s23
             ix=max(i3p,i4p)                                            7d19s23
             in=min(i3p,i4p)                                            7d19s23
             iad=icol+((ix*(ix+1))/2)+in                                7d19s23
             sum=sum+bc(iadd+i3)*bc(iad)                                7d19s23
            end do                                                      7d19s23
            iadd=iadd+nh0av(isb)                                        7d19s23
           end do                                                       7d19s23
           jtt=itt(jsb)+i2+nbasdws(jsb)*ia                              7d20s23
           bc(jtt)=bc(jtt)+4d0*sum                                      7d20s23
          end if                                                        7d19s23
         end do                                                         7d19s23
        end do                                                          7d19s23
c
       end if                                                           7d19s23
       if(isblk(1,is).eq.isblk(4,is))then                               10d26s23
        isb=isblk(1,is)                                                 10d26s23
        jsb=isblk(2,is)                                                 7d19s23
        if(isblk(1,is).eq.isblk(2,is))then                              7d20s23
         nrow12=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2               7d20s23
         iswitch=0                                                      7d20s23
        else                                                            7d20s23
         nrow12=noc(isblk(1,is))*noc(isblk(2,is))                       7d20s23
         iswitch=1                                                      7d20s23
        end if                                                          7d20s23
        if(idoubo(jsb).gt.0)then                                        10d25s23
         call ilimts(nvirt(jsb),nvirt(isb),mynprocg,mynowprog,il,ih,    10d25s23
     $       i1s,i1e,i2s,i2e)                                           10d25s23
         jj=jmats(is)                                                   10d25s23
         i10=i1s                                                        10d25s23
         i1n=nvirt(jsb)                                                 10d25s23
         do i2=i2s,i2e                                                  10d25s23
          if(i2.eq.i2e)i1n=i1e                                          10d25s23
          iadd=iden(isb)+nh0av(isb)*(irefo(isb)+i2-1)                   10d25s23
          do i1=i10,i1n                                                 10d25s23
           do id=0,idoubo(jsb)-1                                        10d25s23
            jtt=itt(jsb)+id+nbasdws(jsb)*(noc(jsb)+i1-1)                10d25s23
            do ia=0,irefo(isb)-1                                        10d25s23
             iap=ia+idoubo(isb)                                         10d25s23
             irec=iap+noc(isb)*id                                       10d25s23
             itri=((iap*(iap+1))/2)+id                                  10d25s23
             itri=itri+iswitch*(irec-itri)+jj                           10d25s23
             bc(jtt)=bc(jtt)-2d0*bc(iadd+ia)*bc(itri)                   10d25s23
            end do                                                      10d25s23
           end do                                                       10d25s23
           jj=jj+nrow12                                                 10d25s23
          end do                                                        10d25s23
          i10=1                                                         10d25s23
         end do                                                         10d25s23
        end if                                                          10d25s23
        call ilimts(noc(jsb),noc(isb),mynprocg,mynowprog,il,ih,i1s,     7d20s23
     $        i1e,i2s,i2e)                                              7d19s23
        nhere=ih+1-il                                                   7d20s23
        if(nhere.gt.0)then                                              7d20s23
         do i2=0,noc(isb)-1                                             7d20s23
          do ia=0,irefo(jsb)-1                                           7d20s23
           iadd=iden(jsb)+nh0av(jsb)*ia                                  7d20s23
           iap=ia+idoubo(jsb)                                            7d20s23
           icol=iap+1+noc(jsb)*i2                                        7d20s23
           if(icol.ge.il.and.icol.le.ih)then                             7d20s23
            icol0=icol                                                   7d20s23
            icol=ioooo(is)+nrow12*(icol-il)                              7d20s23
            do i4=0,irefo(jsb)-1                                         7d20s23
             ff=-2d0*bc(iadd+i4)                                         7d20s23
             i4p=i4+idoubo(jsb)                                          7d20s23
             do i3=0,idoubo(isb)-1                                      7d20s23
              jtt=itt(isb)+i3+nbasdws(isb)*i2                           7d20s23
              ix=max(i3,i4p)                                             7d20s23
              in=min(i3,i4p)                                             7d20s23
              itri=((ix*(ix+1))/2)+in                                    7d20s23
              irec=i3+noc(isb)*i4p                                       7d20s23
              itri=itri+iswitch*(irec-itri)                              7d20s23
              iad=icol+itri                                              7d20s23
              bc(jtt)=bc(jtt)+ff*bc(iad)                                 7d20s23
             end do                                                      7d20s23
            end do                                                       7d20s23
           end if                                                        7d20s23
          end do                                                         7d20s23
         end do                                                          7d20s23
         if(isb.ne.jsb)then                                             7d20s23
          do ia=0,irefo(isb)-1                                             7d19s23
           iadd=iden(isb)+nh0av(isb)*ia                                  7d19s23
           iap=ia+idoubo(isb)                                            7d19s23
           do i2=0,noc(jsb)-1                                           7d20s23
            icol=i2+1+noc(jsb)*iap                                       7d19s23
            if(icol.ge.il.and.icol.le.ih)then                            7d19s23
             icol0=icol
             icol=ioooo(is)+nrow12*(icol-il)                             7d19s23
             do i4=0,idoubo(jsb)-1                                      7d20s23
              jtt=itt(jsb)+i4+nbasdws(jsb)*i2                           7d20s23
              do i3=0,irefo(isb)-1                                       7d19s23
               i3p=i3+idoubo(isb)                                        7d19s23
               ix=max(i3p,i4)                                           7d19s23
               in=min(i3p,i4)                                           7d19s23
               itri=((ix*(ix+1))/2)+in                                   7d19s23
               irec=i3p+noc(isb)*i4                                     7d19s23
               itri=itri+iswitch*(irec-itri)                             7d19s23
               iad=icol+itri                                             7d19s23
               bc(jtt)=bc(jtt)-2d0*bc(iadd+i3)*bc(iad)                   7d19s23
              end do                                                     7d19s23
             end do                                                      7d19s23
            end if                                                       7d19s23
           end do                                                        7d19s23
          end do                                                         7d19s23
         end if                                                         7d20s23
         if(min(irefo(isb),idoubo(jsb)).gt.0)then                       7d20s23
          itmp=ibcoff                                                   7d20s23
          ibcoff=itmp+irefo(isb)*noc(isb)                               7d20s23
          call enough('dorbmix.k',bc,ibc)                               7d20s23
          do iz=itmp,ibcoff-1                                           7d20s23
           bc(iz)=0d0                                                   7d20s23
          end do                                                        7d20s23
          jtmp=itmp                                                     7d20s23
          do ia=0,noc(isb)-1                                            7d20s23
           do i1=0,idoubo(jsb)-1                                        7d20s23
            icol=1+i1+noc(jsb)*ia                                       7d20s23
            if(icol.ge.il.and.icol.le.ih)then                           7d20s23
             icol=ioooo(is)+nrow12*(icol-il)                            7d25s23
c     d|d eta)                                                        (a
             do ial=0,irefo(isb)-1                                      7d20s23
              ialp=ial+idoubo(isb)                                      7d20s23
              ix=max(ialp,i1)                                           7d20s23
              in=min(ialp,i1)                                           7d20s23
              itri=((ix*(ix+1))/2)+in                                   7d20s23
              irec=ialp+noc(isb)*i1                                     7d20s23
              itri=itri+iswitch*(irec-itri)                             7d20s23
              iad=icol+itri                                             7d20s23
              bc(jtmp+ial)=bc(jtmp+ial)+bc(iad)                         7d20s23
             end do                                                     7d20s23
            end if                                                      7d20s23
           end do                                                       7d20s23
           jtmp=jtmp+irefo(isb)                                         7d20s23
          end do                                                        7d20s23
          jtt=itt(isb)+idoubo(isb)                                      7d20s23
          call dgemm('n','n',nh0av(isb),noc(isb),irefo(isb),-2d0,       7d20s23
     $         bc(iden(isb)),nh0av(isb),bc(itmp),irefo(isb),1d0,        7d20s23
     $         bc(jtt),nbasdws(isb),'dorbmix.tmp2k')                    7d20s23
          ibcoff=itmp                                                   7d20s23
         end if                                                         7d20s23
        end if                                                          7d20s23
       else if(isblk(2,is).eq.isblk(4,is).and.                          10d26s23
     $       isblk(1,is).ne.isblk(2,is))then                            7d21s23
        isb=isblk(2,is)                                                 10d26s23
        jsb=isblk(1,is)                                                 7d19s23
        nrow12=noc(isblk(1,is))*noc(isblk(2,is))                        7d20s23
        if(idoubo(jsb).gt.0)then                                        10d25s23
         call ilimts(nvirt(jsb),nvirt(isb),mynprocg,mynowprog,il,ih,    10d25s23
     $       i1s,i1e,i2s,i2e)                                           10d25s23
         jj=jmats(is)                                                   10d25s23
         i10=i1s                                                        10d25s23
         i1n=nvirt(jsb)                                                 10d25s23
         do i2=i2s,i2e                                                  10d25s23
          if(i2.eq.i2e)i1n=i1e                                          10d25s23
          iadd=iden(isb)+nh0av(isb)*(irefo(isb)+i2-1)                   10d25s23
          do i1=i10,i1n                                                 10d25s23
           jtt=itt(jsb)+nbasdws(jsb)*(noc(jsb)+i1-1)                    10d25s23
           do ia=0,irefo(isb)-1                                         10d25s23
            do id=0,idoubo(jsb)-1                                        10d25s23
             iap=ia+idoubo(isb)                                         10d25s23
             irec=id+noc(jsb)*iap+jj                                    10d25s23
             bc(jtt+id)=bc(jtt+id)-2d0*bc(iadd+ia)*bc(irec)             10d25s23
            end do                                                      10d25s23
           end do                                                       10d25s23
           jj=jj+nrow12                                                 10d25s23
          end do                                                        10d25s23
          i10=1                                                         10d25s23
         end do                                                         10d25s23
        end if                                                          10d25s23
        call ilimts(noc(jsb),noc(isb),mynprocg,mynowprog,il,ih,i1s,     10d24s23
     $        i1e,i2s,i2e)                                              7d19s23
        nhere=ih+1-il                                                   7d20s23
        if(min(nhere,idoubo(jsb)).gt.0)then                             7d20s23
         itmp=ibcoff                                                    7d20s23
         ibcoff=itmp+noc(isb)*irefo(isb)                                7d20s23
         call enough('dorbmix.kb',bc,ibc)                               7d20s23
         do iz=itmp,ibcoff-1                                            7d20s23
          bc(iz)=0d0                                                    7d20s23
         end do                                                         7d20s23
         do ia=0,noc(isb)-1                                             10d24s23
          jtmp=itmp+irefo(isb)*ia                                       10d24s23
          do i1=0,idoubo(jsb)-1                                          7d20s23
           icol=1+i1+noc(jsb)*ia                                        10d24s23
           icol0=icol                                                   10d23s23
           if(icol.ge.il.and.icol.le.ih)then                            7d20s23
            icol=ioooo(is)+nrow12*(icol-il)                             7d20s23
c     a|eta d) -2.                                                    (d
            do ial=0,irefo(isb)-1                                       7d20s23
             ialp=ial+idoubo(isb)                                       7d20s23
             irec=i1+noc(jsb)*ialp                                      7d20s23
             iad=icol+irec                                              7d20s23
             bc(jtmp+ial)=bc(jtmp+ial)+bc(iad)                          7d20s23
            end do                                                      7d20s23
           end if                                                       7d20s23
          end do                                                        7d20s23
         end do                                                         7d20s23
         jtt=itt(isb)+idoubo(isb)                                       7d20s23
         call dgemm('n','n',nh0av(isb),noc(isb),irefo(isb),-2d0,        7d20s23
     $        bc(iden(isb)),nh0av(isb),bc(itmp),irefo(isb),1d0,         10d26s23
     $        bc(jtt),nbasdws(isb),'dorbmix.tmp2kb')                    10d26s23
         ibcoff=itmp                                                     7d20s23
        end if                                                          7d20s23
       end if                                                           7d19s23
      end do                                                            7d19s23
      do is=1,nsdlk1                                                    7d21s23
       if(isblk1(1,is).eq.isblk1(2,is))then                             10d26s23
        isb=isblk1(3,is)                                                10d26s23
        jsb=isblk1(1,is)                                                7d21s23
        call ilimts(noc(isblk1(3,is)),nvirt(isblk1(4,is)),mynprocg,     7d21s23
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           7d21s23
        nrow=(noc(jsb)*(noc(jsb)+1))/2                                  7d21s23
        nhere=ih+1-il                                                   7d21s23
        if(nhere.gt.0)then                                              7d21s23
         if(idoubo(isb).ne.0)then                                       10d26s23
          nrow3=(nvirt(jsb)*(nvirt(jsb)+1))/2                           10d26s23
          idtmp=ibcoff                                                  10d26s23
          ibcoff=idtmp+nrow3                                            10d26s23
          call enough('dorbmix.idtmp',bc,ibc)                           10d26s23
          jdtmp=idtmp                                                   10d26s23
          do i3=0,nvirt(jsb)-1                                          10d26s23
           iadd=iden(jsb)+irefo(jsb)+nh0av(jsb)*(i3+irefo(jsb))         10d26s23
           do i4=0,i3-1                                                 10d26s23
            bc(jdtmp+i4)=bc(iadd+i4)*8d0                                10d26s23
           end do                                                       10d26s23
           jdtmp=jdtmp+i3                                               10d26s23
           bc(jdtmp)=bc(iadd+i3)*4d0                                    10d26s23
           jdtmp=jdtmp+1                                                10d26s23
          end do                                                        10d26s23
          do i2=i2s,i2e                                                 10d26s23
           i2p=i2-1+noc(isb)                                            10d26s23
           do i1=1,idoubo(isb)                                          10d26s23
            icol=i1+noc(isb)*(i2-1)                                     10d26s23
            if(icol.ge.il.and.icol.le.ih)then                           10d26s23
             jtt=itt(isb)+i1-1+nbasdws(isb)*i2p                         10d26s23
             j3x=i3x(is)+nrow3*(icol-il)                                10d26s23
             do i=0,nrow3-1
              bc(jtt)=bc(jtt)+bc(j3x+i)*bc(idtmp+i)                     10d26s23
             end do                                                     10d26s23
            end if                                                      10d26s23
           end do                                                       10d26s23
          end do                                                        10d26s23
          ibcoff=idtmp                                                  10d26s23
         end if                                                         10d26s23
         if(min(irefo(isb),idoubo(jsb)).gt.0)then                       7d20s23
          itmp=ibcoff                                                   7d20s23
          ibcoff=itmp+irefo(isb)*nvirt(isb)                             7d20s23
          call enough('dorbmix.tmp1a',bc,ibc)                           7d20s23
          do iz=itmp,ibcoff-1                                           7d20s23
           bc(iz)=0d0                                                   7d20s23
          end do                                                        7d20s23
          do iv=0,nvirt(isb)-1                                          7d20s23
           do i1=0,irefo(isb)-1                                         7d20s23
            iad=itmp+i1+irefo(isb)*iv                                   7d20s23
            i1p=i1+idoubo(isb)                                          7d20s23
            icol=i1p+1+noc(isb)*iv                                      7d20s23
            if(icol.ge.il.and.icol.le.ih)then                           7d20s23
             icol=ionex(is)+nrow*(icol-il)-1                            7d20s23
             do id=1,idoubo(jsb)                                        7d20s23
              itri=((id*(id+1))/2)+icol                                 7d20s23
              bc(iad)=bc(iad)+bc(itri)                                  7d20s23
             end do                                                     7d20s23
            end if                                                      7d20s23
           end do                                                       7d20s23
          end do                                                        7d20s23
          jtt=itt(isb)+idoubo(isb)+nbasdws(isb)*noc(isb)                7d20s23
          call dgemm('n','n',nh0av(isb),nvirt(isb),irefo(isb),4d0,      10d24s23
     $         bc(iden(isb)),nh0av(isb),bc(itmp),irefo(isb),1d0,        10d26s23
     $         bc(jtt),nbasdws(isb),'dorbmix.tmp21a')                   10d26s23
          ibcoff=itmp                                                   7d20s23
         end if                                                         7d20s23
         if(min(noc(isb),idoubo(jsb)).gt.0)then                         10d24s23
          itmp=ibcoff                                                   7d20s23
          ibcoff=itmp+nvirt(isb)*noc(isb)                               10d24s23
          call enough('dorbmix.tmp1a',bc,ibc)                           7d20s23
          do iz=itmp,ibcoff-1                                           10d24s23
           bc(iz)=0d0                                                   7d20s23
          end do                                                        7d20s23
          do iv=0,nvirt(isb)-1                                          7d20s23
           do i1=0,noc(isb)-1                                           10d24s23
            iad=itmp+iv+nvirt(isb)*i1                                   10d24s23
            icol=i1+1+noc(isb)*iv                                       10d24s23
            if(icol.ge.il.and.icol.le.ih)then                           7d20s23
             icol=ionex(is)+nrow*(icol-il)-1                            7d20s23
             do id=1,idoubo(jsb)                                        7d20s23
              itri=((id*(id+1))/2)+icol                                 7d20s23
              bc(iad)=bc(iad)+bc(itri)                                  7d20s23
             end do                                                     7d20s23
            end if                                                      7d20s23
           end do                                                       7d20s23
          end do                                                        7d20s23
          jden=iden(isb)+nh0av(isb)*irefo(isb)                          10d24s23
          jtt=itt(isb)+idoubo(isb)                                      10d24s23
          call dgemm('n','n',nh0av(isb),noc(isb),nvirt(isb),4d0,        10d24s23
     $         bc(jden),nh0av(isb),bc(itmp),nvirt(isb),1d0,             10d26s23
     $         bc(jtt),nbasdws(isb),'dorbmix.tmp21b')                   10d26s23
          ibcoff=itmp                                                   7d20s23
         end if                                                         7d20s23
         if(idoubo(jsb).gt.0)then                                       10d24s23
          do iv=0,nvirt(isb)-1                                          10d24s23
           ivp=iv+irefo(isb)                                            10d24s23
           do ia=0,irefo(isb)-1                                         10d24s23
            iap=ia+idoubo(isb)                                          10d24s23
            icol=iap+1+noc(isb)*iv                                      10d24s23
            if(icol.ge.il.and.icol.le.ih)then                           10d24s23
             iad=iden(isb)+ia+nh0av(isb)*ivp                            10d24s23
             ff=8d0*bc(iad)                                             10d24s23
             icol=ionex(is)+nrow*(icol-il)                              10d24s23
             do id=0,idoubo(jsb)-1                                      10d24s23
              jtt=itt(jsb)+id                                           10d24s23
              do ie=0,noc(jsb)-1                                        10d24s23
               ix=max(ie,id)                                            10d24s23
               in=min(ie,id)                                            10d24s23
               itri=((ix*(ix+1))/2)+in+icol                             10d24s23
               bc(jtt)=bc(jtt)+ff*bc(itri)                              10d24s23
               jtt=jtt+nbasdws(jsb)                                     10d24s23
              end do                                                    10d24s23
             end do                                                     10d24s23
            end if                                                      10d24s23
           end do                                                       10d24s23
          end do                                                        10d24s23
         end if                                                         10d24s23
         if(min(irefo(jsb),idoubo(isb)).gt.0)then                        7d21s23
          do iv=0,nvirt(isb)-1                                          7d21s23
           ivp=iv+noc(isb)                                              7d20s23
           do idd=0,idoubo(isb)-1                                       7d21s23
            jtt=itt(isb)+idd+nbasdws(isb)*ivp                           7d20s23
            icol=idd+1+noc(isb)*iv                                      7d21s23
            if(icol.ge.il.and.icol.le.ih)then                           7d21s23
             icol=ionex(is)+nrow*(icol-il)                              7d21s23
             iadd=iden(jsb)                                             7d21s23
             sum=0d0                                                    7d21s23
             do i4=0,irefo(jsb)-1                                       7d21s23
              i4p=i4+idoubo(jsb)                                        7d21s23
              do i3=0,irefo(jsb)-1                                      7d21s23
               i3p=i3+idoubo(jsb)                                       7d21s23
               ix=max(i3p,i4p)                                          7d21s23
               in=min(i3p,i4p)                                          7d21s23
               itri=((ix*(ix+1))/2)+in                                  7d21s23
               sum=sum+bc(iadd+i3)*bc(icol+itri)                        7d21s23
              end do                                                    7d21s23
              iadd=iadd+nh0av(jsb)                                      7d21s23
             end do                                                     7d21s23
             bc(jtt)=bc(jtt)+4d0*sum                                    7d20s23
            end if                                                      7d21s23
           end do                                                       7d21s23
          end do                                                        7d21s23
         end if                                                         7d21s23
        end if                                                          7d21s23
       end if                                                           7d21s23
       if(isblk1(1,is).eq.isblk1(4,is))then                             10d26s23
        isb=isblk1(1,is)                                                10d26s23
        jsb=isblk1(2,is)                                                7d21s23
        call ilimts(noc(isblk1(3,is)),nvirt(isblk1(4,is)),mynprocg,     7d21s23
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           7d21s23
        if(isb.eq.jsb)then                                              7d21s23
         nrow=(noc(jsb)*(noc(jsb)+1))/2                                  7d21s23
         nrow3=(nvirt(jsb)*(nvirt(jsb)+1))/2                            10d26s23
         iswitch=0                                                      7d21s23
        else                                                            7d21s23
         nrow=noc(isb)*noc(jsb)                                         7d21s23
         nrow3=nvirt(isb)*nvirt(jsb)                                    10d26s23
         iswitch=1                                                      7d21s23
        end if                                                          7d21s23
        nhere=ih+1-il                                                   7d21s23
        if(nhere.gt.0)then                                              7d21s23
         if(idoubo(jsb).ne.0)then                                       10d26s23
          do i2=i2s,i2e                                                 10d26s23
           iadd=iden(isb)+irefo(isb)+nh0av(isb)*(i2-1+irefo(isb))       10d26s23
           do id=1,idoubo(jsb)                                          10d26s23
            icol=id+noc(jsb)*(i2-1)                                     10d26s23
            if(icol.ge.il.and.icol.le.ih)then                           10d26s23
             j3x=i3x(is)+nrow3*(icol-il)                                10d26s23
             do i4=0,nvirt(jsb)-1                                       10d26s23
              jtt=itt(jsb)+id-1+nbasdws(jsb)*(i4+noc(jsb))              10d26s23
              do i3=0,nvirt(isb)-1                                      10d26s23
               irec=i3+nvirt(isb)*i4                                    10d26s23
               ix=max(i3,i4)                                            10d26s23
               in=min(i3,i4)
               itri=((ix*(ix+1))/2)+in                                  10d26s23
               itri=itri+iswitch*(irec-itri)+j3x                        10d26s23
               bc(jtt)=bc(jtt)-2d0*bc(iadd+i3)*bc(itri)                 10d26s23
              end do                                                    10d26s23
             end do                                                     10d26s23
            end if                                                      10d26s23
           end do                                                       10d26s23
          end do                                                        10d26s23
         end if                                                         10d26s23
         if(idoubo(jsb).gt.0)then                                       10d25s23
          do iv=0,nvirt(isb)-1                                           10d25s23
           ivp=iv+irefo(isb)                                            10d25s23
           jden=iden(isb)+nh0av(isb)*ivp                                10d25s23
           do ie=0,noc(jsb)-1                                            10d25s23
            icol=ie+1+noc(jsb)*iv                                        10d25s23
            if(icol.ge.il.and.icol.le.ih)then                            10d25s23
             jtt=itt(jsb)+nbasdws(jsb)*ie                               10d25s23
             icol=ionex(is)+nrow*(icol-il)                              10d25s23
             do id=0,idoubo(jsb)-1                                      10d25s23
              do ia=0,irefo(isb)-1                                      10d25s23
               iap=ia+idoubo(isb)                                       10d25s23
               irec=iap+noc(isb)*id                                     10d25s23
               itri=((iap*(iap+1))/2)+id                                10d25s23
               itri=itri+iswitch*(irec-itri)+icol                       10d25s23
               bc(jtt+id)=bc(jtt+id)-2d0*bc(itri)*bc(jden+ia)           10d25s23
              end do                                                    10d25s23
             end do                                                     10d25s23
            end if                                                      10d25s23
           end do                                                       10d25s23
           do id=0,idoubo(jsb)-1                                        10d25s23
            icol=id+1+noc(jsb)*iv                                       10d25s23
            icol0=icol
            if(icol.ge.il.and.icol.le.ih)then                           10d25s23
             icol=ionex(is)+nrow*(icol-il)                              10d25s23
             do ie=0,noc(jsb)-1                                         10d25s23
              jtt=itt(jsb)+id+nbasdws(jsb)*ie                           10d25s23
              do ia=0,irefo(isb)-1                                      10d25s23
               iap=ia+idoubo(isb)                                       10d25s23
               irec=iap+noc(isb)*ie                                     10d25s23
               ix=max(ie,iap)                                           10d25s23
               in=min(ie,iap)                                           10d25s23
               itri=((ix*(ix+1))/2)+in                                  10d25s23
               itri=itri+iswitch*(irec-itri)+icol                       10d25s23
               bc(jtt)=bc(jtt)-2d0*bc(jden+ia)*bc(itri)                 10d25s23
              end do                                                    10d25s23
             end do                                                     10d25s23
            end if                                                      10d25s23
           end do                                                       10d25s23
          end do                                                         10d25s23
         end if                                                         10d25s23
         if(min(irefo(isb),idoubo(jsb)).gt.0)then                       7d20s23
          itmp=ibcoff                                                   7d20s23
          ibcoff=itmp+irefo(isb)*nvirt(isb)                             7d20s23
          call enough('dorbmix.tmp1b',bc,ibc)                           7d20s23
          do iz=itmp,ibcoff-1                                           10d25s23
           bc(iz)=0d0                                                   10d25s23
          end do                                                        10d25s23
          do iv=0,nvirt(isb)-1                                          7d20s23
           jtmp=itmp+irefo(isb)*iv                                      7d20s23
           do id=0,idoubo(jsb)-1                                        7d20s23
            icol=id+1+noc(jsb)*iv                                       7d20s23
            if(icol.ge.il.and.icol.le.ih)then                           7d20s23
             icol=ionex(is)+nrow*(icol-il)                              7d20s23
             do ia=0,irefo(isb)-1                                       7d20s23
              iap=ia+idoubo(isb)                                        7d20s23
              irec=iap+noc(isb)*id                                      7d20s23
              itri=((iap*(iap+1))/2)+id                                 7d20s23
              itri=itri+iswitch*(irec-itri)+icol                        7d20s23
              bc(jtmp+ia)=bc(jtmp+ia)+bc(itri)                          7d20s23
             end do                                                     7d20s23
            end if                                                      7d20s23
           end do                                                       7d20s23
          end do                                                        7d20s23
          jtt=itt(isb)+idoubo(isb)+nbasdws(isb)*noc(isb)                7d20s23
          call dgemm('n','n',nh0av(isb),nvirt(isb),irefo(isb),-2d0,     10d24s23
     $         bc(iden(isb)),nh0av(isb),bc(itmp),irefo(isb),1d0,        10d26s23
     $         bc(jtt),nbasdws(isb),'dorbmix.tmp21b')                   10d26s23
          ibcoff=itmp                                                   7d20s23
         end if                                                         7d20s23
         if(min(irefo(jsb),idoubo(isb)).gt.0)then                       7d21s23
          do iv=0,nvirt(isb)-1                                          7d21s23
           ivp=iv+noc(isb)                                              7d20s23
           jtt=itt(isb)+nbasdws(isb)*ivp                                7d20s23
           iadd=iden(jsb)                                               7d21s23
           do i1=0,irefo(jsb)-1                                         7d21s23
            i1p=i1+idoubo(jsb)                                          7d21s23
            icol=i1p+1+noc(jsb)*iv                                      7d21s23
            if(icol.ge.il.and.icol.le.ih)then                           7d21s23
             icol=ionex(is)+nrow*(icol-il)                              7d21s23
             do i4=0,irefo(jsb)-1                                       7d21s23
              i4p=i4+idoubo(jsb)                                        7d21s23
              do i3=0,idoubo(isb)-1                                     7d21s23
               ix=max(i4p,i3)                                           7d21s23
               in=min(i4p,i3)                                           7d21s23
               itri=((ix*(ix+1))/2)+in                                  7d21s23
               irec=i3+noc(isb)*i4p                                     7d21s23
               itri=itri+iswitch*(irec-itri)                            7d21s23
               term=-2d0*bc(iadd+i4)*bc(icol+itri)                      7d20s23
               bc(jtt+i3)=bc(jtt+i3)+term                               7d20s23
              end do                                                    7d21s23
             end do                                                     7d21s23
            end if                                                      7d21s23
            iadd=iadd+nh0av(jsb)                                        7d21s23
           end do                                                       7d21s23
          end do                                                        7d21s23
         end if                                                         7d21s23
         if(min(noc(isb),idoubo(jsb)).gt.0)then                         10d24s23
          itmp=ibcoff                                                   10d24s23
          ibcoff=itmp+nvirt(isb)*noc(isb)                               10d24s23
          call enough('dorbmix.gg',bc,ibc)                              10d24s23
          do iz=itmp,ibcoff-1                                           10d24s23
           bc(iz)=0d0                                                   10d24s23
          end do                                                        10d24s23
          do iv=0,nvirt(isb)-1                                          10d24s23
           do id=0,idoubo(jsb)-1                                        10d24s23
            icol=1+id+noc(jsb)*iv                                       10d24s23
            if(icol.ge.il.and.icol.le.ih)then                           10d24s23
             icol=ionex(is)+nrow*(icol-il)                              10d24s23
             jtmp=itmp+iv                                               10d24s23
             do ie=0,noc(isb)-1                                         10d24s23
              irec=ie+noc(isb)*id                                       10d24s23
              ix=max(id,ie)                                             10d24s23
              in=min(id,ie)                                             10d24s23
              itri=((ix*(ix+1))/2)+in                                   10d24s23
              itri=itri+iswitch*(irec-itri)+icol                        10d24s23
              bc(jtmp)=bc(jtmp)+bc(itri)                                10d24s23
              jtmp=jtmp+nvirt(isb)                                      10d24s23
             end do                                                     10d24s23
            end if                                                      10d24s23
           end do                                                       10d24s23
          end do                                                        10d24s23
          jden=iden(isb)+nh0av(isb)*irefo(isb)                          10d24s23
          jtt=itt(isb)+idoubo(isb)                                      10d24s23
          call dgemm('n','n',nh0av(isb),noc(isb),nvirt(isb),-2d0,       10d24s23
     $         bc(jden),nh0av(isb),bc(itmp),nvirt(isb),1d0,             10d26s23
     $         bc(jtt),nbasdws(isb),'dorbmix.gg2')                      10d26s23
          ibcoff=itmp                                                   10d24s23
         end if                                                         10d24s23
        end if                                                          7d21s23
       else if(isblk1(2,is).eq.isblk1(4,is).and.                        10d26s23
     $       isblk1(1,is).ne.isblk1(2,is))then                          10d26s23
        isb=isblk1(2,is)                                                10d26s23
        jsb=isblk1(1,is)                                                7d21s23
        call ilimts(noc(isblk1(3,is)),nvirt(isblk1(4,is)),mynprocg,     7d21s23
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           7d21s23
        nrow=noc(isb)*noc(jsb)                                          7d21s23
        nrow3=nvirt(isb)*nvirt(jsb)                                     10d26s23
        nhere=ih+1-il                                                   7d21s23
        if(nhere.gt.0)then                                              7d21s23
         if(idoubo(jsb).ne.0)then                                       10d26s23
          do i2=i2s,i2e                                                 10d26s23
           iadd=iden(isb)+irefo(isb)+nh0av(isb)*(i2-1+irefo(isb))       10d26s23
           do id=1,idoubo(jsb)                                          10d26s23
            icol=id+noc(jsb)*(i2-1)                                     10d26s23
            if(icol.ge.il.and.icol.le.ih)then                           10d26s23
             j3x=i3x(is)+nrow3*(icol-il)                                10d26s23
             do i3=0,nvirt(isb)-1                                       10d26s23
              ff=-2d0*bc(iadd+i3)                                       10d26s23
              jtt=itt(jsb)+id-1+nbasdws(jsb)*noc(jsb)                   10d26s23
              do i4=0,nvirt(jsb)-1                                       10d26s23
               irec=i4+nvirt(jsb)*i3+j3x                                10d26s23
               bc(jtt)=bc(jtt)+ff*bc(irec)                              10d26s23
               jtt=jtt+nbasdws(jsb)                                     10d26s23
              end do                                                    10d26s23
             end do                                                     10d26s23
            end if                                                      10d26s23
           end do                                                       10d26s23
          end do                                                        10d26s23
         end if                                                         10d26s23
         if(idoubo(jsb).gt.0)then                                       10d25s23
          do iv=0,nvirt(isb)-1                                           10d25s23
           ivp=iv+irefo(isb)                                            10d25s23
           jden=iden(isb)+nh0av(isb)*ivp                                10d25s23
           do ie=0,noc(jsb)-1                                            10d25s23
            icol=ie+1+noc(jsb)*iv                                        10d25s23
            if(icol.ge.il.and.icol.le.ih)then                            10d25s23
             jtt=itt(jsb)+nbasdws(jsb)*ie                               10d25s23
             icol=ionex(is)+nrow*(icol-il)                              10d25s23
             do ia=0,irefo(isb)-1                                       10d25s23
              iap=ia+idoubo(isb)                                        10d25s23
              do id=0,idoubo(jsb)-1                                      10d25s23
               irec=id+noc(jsb)*iap+icol                                10d25s23
               bc(jtt+id)=bc(jtt+id)-2d0*bc(irec)*bc(jden+ia)           10d25s23
              end do                                                    10d25s23
             end do                                                     10d25s23
            end if                                                      10d25s23
           end do                                                       10d25s23
           do id=0,idoubo(jsb)-1                                        10d25s23
            icol=id+1+noc(jsb)*iv                                       10d25s23
            if(icol.ge.il.and.icol.le.ih)then                           10d25s23
             icol=ionex(is)+nrow*(icol-il)                              10d25s23
             do ia=0,irefo(isb)-1                                       10d25s23
              iap=ia+idoubo(isb)                                        10d25s23
              ff=-2d0*bc(jden+ia)                                       10d25s23
              jtt=itt(jsb)+id                                           10d25s23
              do ie=0,noc(jsb)-1                                         10d25s23
               irec=ie+noc(jsb)*iap+icol                                10d25s23
               bc(jtt)=bc(jtt)+ff*bc(irec)                              10d25s23
               jtt=jtt+nbasdws(jsb)                                     10d25s23
              end do                                                    10d25s23
             end do                                                     10d25s23
            end if                                                      10d25s23
           end do                                                       10d25s23
          end do                                                         10d25s23
         end if                                                         10d25s23
         if(min(irefo(isb),idoubo(jsb)).gt.0)then                       7d20s23
          itmp=ibcoff                                                   7d20s23
          ibcoff=itmp+irefo(isb)*nvirt(isb)                             7d20s23
          call enough('dorbmix.tmp1b',bc,ibc)                           7d20s23
          do iz=itmp,ibcoff-1                                           10d25s23
           bc(iz)=0d0                                                   10d25s23
          end do                                                        10d25s23
          do iv=0,nvirt(isb)-1                                          7d20s23
           jtmp=itmp+irefo(isb)*iv                                      7d20s23
           do id=0,idoubo(jsb)-1                                        7d20s23
            icol=id+1+noc(jsb)*iv                                       7d20s23
            if(icol.ge.il.and.icol.le.ih)then                           7d20s23
             icol=ionex(is)+nrow*(icol-il)                              7d20s23
             do ia=0,irefo(isb)-1                                       7d20s23
              iap=ia+idoubo(isb)                                        7d20s23
              irec=id+noc(jsb)*iap+icol                                 10d24s23
              bc(jtmp+ia)=bc(jtmp+ia)+bc(irec)                          10d24s23
             end do                                                     7d20s23
            end if                                                      7d20s23
           end do                                                       7d20s23
          end do                                                        7d20s23
          jtt=itt(isb)+idoubo(isb)+nbasdws(isb)*noc(isb)                7d20s23
          call dgemm('n','n',nh0av(isb),nvirt(isb),irefo(isb),-2d0,     10d24s23
     $         bc(iden(isb)),nh0av(isb),bc(itmp),irefo(isb),1d0,        10d26s23
     $         bc(jtt),nbasdws(isb),'dorbmix.tmp21b')                   10d26s23
          ibcoff=itmp                                                   7d20s23
         end if                                                         7d20s23
         if(min(irefo(jsb),idoubo(isb)).gt.0)then                       7d21s23
          do iv=0,nvirt(isb)-1                                          7d21s23
           iadd=iden(jsb)                                               7d21s23
           ivp=iv+noc(isb)                                              7d20s23
           jtt=itt(isb)+nbasdws(isb)*ivp                                7d20s23
           do i1=0,irefo(jsb)-1                                         7d21s23
            i1p=i1+idoubo(jsb)                                          7d21s23
            icol=i1p+1+noc(jsb)*iv                                      7d21s23
            if(icol.ge.il.and.icol.le.ih)then                           7d21s23
             icol=ionex(is)+nrow*(icol-il)                              7d21s23
             do i3=0,idoubo(isb)-1                                      7d21s23
              do i4=0,irefo(jsb)-1                                       7d21s23
               i4p=i4+idoubo(jsb)                                       7d21s23
               irec=i4p+noc(jsb)*i3                                     7d21s23
               term=-2d0*bc(iadd+i4)*bc(icol+irec)                      7d20s23
               bc(jtt+i3)=bc(jtt+i3)+term                               7d20s23
              end do                                                    7d21s23
             end do                                                     7d21s23
            end if                                                      7d21s23
            iadd=iadd+nh0av(jsb)                                        7d21s23
           end do                                                       7d21s23
          end do                                                        7d21s23
         end if                                                         7d21s23
         if(min(noc(isb),idoubo(jsb)).gt.0)then                         10d24s23
          itmp=ibcoff                                                   10d24s23
          ibcoff=itmp+nvirt(isb)*noc(isb)                               10d24s23
          call enough('dorbmix.gg',bc,ibc)                              10d24s23
          do iz=itmp,ibcoff-1                                           10d24s23
           bc(iz)=0d0                                                   10d24s23
          end do                                                        10d24s23
          do iv=0,nvirt(isb)-1                                          10d24s23
           do id=0,idoubo(jsb)-1                                        10d24s23
            icol=1+id+noc(jsb)*iv                                       10d24s23
            if(icol.ge.il.and.icol.le.ih)then                           10d24s23
             icol=ionex(is)+nrow*(icol-il)                              10d24s23
             jtmp=itmp+iv                                               10d24s23
             do ie=0,noc(isb)-1                                         10d24s23
              irec=id+noc(jsb)*ie+icol                                  10d24s23
              bc(jtmp)=bc(jtmp)+bc(irec)                                10d24s23
              jtmp=jtmp+nvirt(isb)                                      10d24s23
             end do                                                     10d24s23
            end if                                                      10d24s23
           end do                                                       10d24s23
          end do                                                        10d24s23
          jden=iden(isb)+nh0av(isb)*irefo(isb)                          10d24s23
          jtt=itt(isb)+idoubo(isb)                                      10d24s23
          call dgemm('n','n',nh0av(isb),noc(isb),nvirt(isb),-2d0,       10d24s23
     $         bc(jden),nh0av(isb),bc(itmp),nvirt(isb),1d0,             10d26s23
     $         bc(jtt),nbasdws(isb),'dorbmix.gg2')                      10d26s23
          ibcoff=itmp                                                   10d24s23
         end if                                                         10d24s23
        end if                                                          7d21s23
       end if                                                           7d21s23
      end do                                                            7d21s23
      do is=1,nsdlkk                                                     10d25s23
       if(isblkk(3,is).eq.isblkk(4,is))then                             10d26s23
        isb=isblkk(3,is)
        jsb=isblkk(1,is)                                                10d25s23
        nrow=noc(jsb)*noc(jsb)                                          10d25s23
        call ilimts(nvirt(isb),nvirt(isb),mynprocg,mynowprog,il,ih,i1s, 10d25s23
     $       i1e,i2s,i2e)                                               10d25s23
        if(idoubo(jsb).ne.0)then                                        10d25s23
         kk=kmats(is)                                                   10d25s23
         i10=i1s                                                        10d25s23
         i1n=nvirt(isb)                                                 10d25s23
         nhere2=i2e+1-i2s                                               10d25s23
         if(nhere2.gt.0)then                                            10d26s23
          itmp=ibcoff                                                    10d25s23
          ibcoff=itmp+nvirt(isb)*nhere2                                  10d25s23
          call enough('dorbmix.goulKd',bc,ibc)                           10d25s23
          do iz=itmp,ibcoff-1                                            10d25s23
           bc(iz)=0d0                                                    10d25s23
          end do                                                         10d25s23
          do i2=i2s,i2e                                                  10d25s23
           i2p=i2-1+irefo(isb)                                           10d25s23
           if(i2.eq.i2e)i1n=i1e                                          10d25s23
           do i1=i10,i1n                                                 10d25s23
            i1p=i1-1+irefo(isb)
            iadd=iden(isb)+i1p+nh0av(isb)*i2p                            10d25s23
            ff=-2d0*bc(iadd)                                             10d25s23
            jtmp=itmp+i1-1+nvirt(isb)*(i2-i2s)                           10d25s23
            do ie=0,noc(jsb)-1                                           10d25s23
             jtt=itt(jsb)+nbasdws(jsb)*ie                                10d25s23
             do id=0,idoubo(jsb)-1                                        10d25s23
              irec=id+noc(jsb)*ie+kk                                     10d25s23
              bc(jtt+id)=bc(jtt+id)+ff*bc(irec)                          10d25s23
             end do                                                      10d25s23
            end do                                                       10d25s23
            do id=0,idoubo(jsb)-1                                        10d25s23
             irec=id+noc(jsb)*id+kk                                      10d25s23
             bc(jtmp)=bc(jtmp)+bc(irec)                                  10d25s23
            end do                                                       10d25s23
            kk=kk+nrow                                                   10d25s23
           end do                                                        10d25s23
           i10=1                                                         10d25s23
          end do                                                         10d25s23
          jden=iden(isb)+nh0av(isb)*irefo(isb)                           10d25s23
          jtt=itt(isb)+idoubo(isb)+nbasdws(isb)*(noc(isb)+i2s-1)           10d26s23
          call dgemm('n','n',nh0av(isb),nhere2,nvirt(isb),-2d0,          10d25s23
     $        bc(jden),nh0av(isb),bc(itmp),nvirt(isb),1d0,              10d26s23
     $        bc(jtt),nbasdws(isb),'dorbmix.goulKdb')                   10d26s23
          ibcoff=itmp                                                    10d25s23
         end if                                                         10d26s23
        end if                                                          10d25s23
       end if                                                           10d25s23
       if(isblkk(2,is).eq.isblkk(4,is))then                             10d26s23
        isb=isblkk(2,is)
        jsb=isblkk(1,is)                                                10d25s23
        nrow=noc(isb)*noc(jsb)                                          10d25s23
        if(idoubo(jsb).gt.0)then                                        10d25s23
         call ilimts(nvirt(jsb),nvirt(isb),mynprocg,mynowprog,il,ih,    10d25s23
     $       i1s,i1e,i2s,i2e)                                           10d25s23
         kk=kmats(is)                                                   10d25s23
         i10=i1s                                                        10d25s23
         i1n=nvirt(jsb)                                                 10d25s23
         do i2=i2s,i2e                                                  10d25s23
          i2p=i2-1+irefo(isb)                                           10d25s23
          iadd=iden(isb)+nh0av(isb)*i2p                                 10d25s23
          if(i2.eq.i2e)i1n=i1e                                          10d25s23
          do i1=i10,i1n                                                 10d25s23
           jtt=itt(jsb)+nbasdws(jsb)*(i1-1+noc(jsb))                    10d25s23
           do ia=0,irefo(isb)-1                                         10d25s23
            iap=ia+idoubo(isb)                                          10d25s23
            do id=0,idoubo(jsb)-1                                       10d26s23
             irec=id+noc(jsb)*iap+kk                                    10d26s23
             bc(jtt+id)=bc(jtt+id)-2d0*bc(iadd+ia)*bc(irec)             10d25s23
            end do                                                      10d25s23
           end do                                                       10d25s23
           kk=kk+nrow                                                   10d25s23
          end do                                                        10d25s23
          i10=1                                                         10d25s23
         end do                                                         10d25s23
        end if                                                          10d25s23
       end if                                                           10d25s23
       if(isblkk(1,is).eq.isblkk(4,is))then                             10d26s23
        isb=isblkk(1,is)
        jsb=isblkk(2,is)                                                10d25s23
        nrow=noc(isb)*noc(jsb)                                          10d25s23
        if(idoubo(jsb).gt.0)then                                        10d25s23
         call ilimts(nvirt(jsb),nvirt(isb),mynprocg,mynowprog,il,ih,    10d25s23
     $       i1s,i1e,i2s,i2e)                                           10d25s23
         nhere=ih+1-il                                                  10d26s23
         if(nhere.gt.0)then                                             10d26s23
          kk=kmats(is)                                                   10d25s23
          i10=i1s                                                        10d25s23
          i1n=nvirt(jsb)                                                 10d25s23
          do i2=i2s,i2e                                                  10d25s23
           if(i2.eq.i2e)i1n=i1e                                          10d25s23
           iadd=iden(isb)+nh0av(isb)*(i2+irefo(isb)-1)                   10d25s23
           do i1=i10,i1n                                                 10d25s23
            jtt=itt(jsb)+nbasdws(jsb)*(i1+noc(jsb)-1)                    10d25s23
            do id=0,idoubo(jsb)-1                                        10d25s23
             do ia=0,irefo(isb)-1                                         10d25s23
              iap=ia+idoubo(isb)                                          10d25s23
              irec=iap+noc(isb)*id+kk                                    10d26s23
              bc(jtt+id)=bc(jtt+id)+8d0*bc(iadd+ia)*bc(irec)             10d25s23
             end do
            end do
            kk=kk+nrow                                                   10d25s23
           end do                                                        10d25s23
           i10=1                                                         10d25s23
          end do                                                         10d25s23
         end if                                                         10d26s23
        end if                                                          10d25s23
       end if                                                           10d25s23
      end do                                                             10d25s23
      return
      end
