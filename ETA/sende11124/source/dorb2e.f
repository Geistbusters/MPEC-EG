c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dorb2e(nsymb,idoubo,irefo,noc,nvirt,isblk,nsdlk,       6d24s24
     $     isblk1,nsdlk1,ionex,ioooo,bc,ibc,itt,nbasdws,igoal)                6d20s24
      implicit real*8 (a-h,o-z)                                         7d20s23
      dimension idoubo(*),irefo(*),noc(*),nvirt(*),                     6d24s24
     $     isblk(4,*),isblk1(4,*),ionex(*),ioooo(*),itt(*),             6d20s24
     $     nbasdws(*)                                                   10d19s23
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      do is=1,nsdlk                                                     5d1s24
       if(isblk(1,is).eq.isblk(2,is))then                               5d1s24
        nrow=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                  5d1s24
        isw=0                                                           5d1s24
       else                                                             5d1s24
        nrow=noc(isblk(1,is))*noc(isblk(2,is))                          5d1s24
        isw=1                                                           5d1s24
       end if                                                           5d1s24
       ncol=noc(isblk(3,is))*noc(isblk(3,is))                           5d1s24
       if(min(nrow,ncol).gt.0)then                                      5d1s24
        if(isblk(1,is).eq.isblk(2,is))then                              5d1s24
         isb=isblk(1,is)                                                5d1s24
         jsb=isblk(3,is)                                                10d19s23
         if(min(idoubo(jsb),idoubo(isb)).gt.0)then                      5d1s24
          call ilimts(noc(jsb),noc(jsb),mynprocg,mynowprog,il,ih,i1s,   10d19s23
     $           i1e,i2s,i2e)                                           10d19s23
          nhere=ih+1-il                                                 10d19s23
          if(nhere.gt.0)then                                            5d1s24
           do idd=0,idoubo(jsb)-1                                       10d19s23
            icol=idd+1+noc(jsb)*idd                                     10d19s23
            if(icol.ge.il.and.icol.le.ih)then                           10d19s23
             icol0=icol                                                 10d19s23
             icol=ioooo(is)+nrow*(icol-il)                              10d19s23
             jff=itt(isb)                                               10d19s23
             do i4=0,noc(isb)-1                                         10d19s23
              iadd=icol+((i4*(i4+1))/2)                                 10d19s23
              do i3=0,idoubo(isb)-1                                     10d19s23
               ix=max(i3,i4)                                            8d16s24
               in=min(i3,i4)                                            8d16s24
               iadd=icol+((ix*(ix+1))/2)+in                             8d16s24
               bc(jff+i3)=bc(jff+i3)+8d0*bc(iadd)                       8d16s24
              end do                                                    10d19s23
              jff=jff+nbasdws(isb)                                      10d19s23
             end do                                                     10d19s23
            end if                                                      10d19s23
           end do                                                       10d19s23
          end if                                                        5d1s24
         end if                                                         5d1s24
        end if                                                          5d1s24
        if(isblk(1,is).eq.isblk(4,is))then                              5d1s24
         isb=isblk(1,is)                                                5d1s24
         jsb=isblk(2,is)                                                5d1s24
         if(min(idoubo(jsb),idoubo(isb)).gt.0)then                      5d1s24
          call ilimts(noc(jsb),noc(isb),mynprocg,mynowprog,il,ih,i1s,   10d19s23
     $           i1e,i2s,i2e)                                           10d19s23
          nhere=ih+1-il                                                 10d19s23
          if(nhere.gt.0)then                                            5d1s24
           do i2=0,noc(isb)-1                                           5d1s24
            jtt=itt(isb)+nbasdws(isb)*i2                                5d1s24
            do i1=1,idoubo(jsb)                                         5d1s24
             icol=i1+noc(jsb)*i2                                        5d1s24
             if(icol.ge.il.and.icol.le.ih)then                          5d1s24
              i1m=i1-1                                                  5d1s24
              icol=ioooo(is)+nrow*(icol-il)                             5d1s24
              do i4=0,idoubo(isb)-1                                     5d1s24
               irec=i4+noc(isb)*i1m                                     5d1s24
               ix=max(i4,i1m)                                           5d1s24
               in=min(i4,i1m)                                           5d1s24
               itri=((ix*(ix+1))/2)+in                                  5d1s24
               itri=itri+isw*(irec-itri)+icol                           5d1s24
               bc(jtt+i4)=bc(jtt+i4)-4d0*bc(itri)                       5d1s24
              end do                                                    5d1s24
             end if                                                     5d1s24
            end do                                                      5d1s24
           end do                                                       5d1s24
           if(isb.ne.jsb)then                                           5d1s24
            do i2=0,idoubo(isb)-1                                       5d1s24
             do i1=1,noc(jsb)                                           5d1s24
              icol=i1+noc(jsb)*i2                                        5d1s24
              if(icol.ge.il.and.icol.le.ih)then                          5d1s24
               i1m=i1-1                                                  5d1s24
               icol=ioooo(is)+nrow*(icol-il)                             5d1s24
               jtt=itt(jsb)+nbasdws(jsb)*i1m                            5d1s24
               do i3=0,idoubo(jsb)-1                                     5d1s24
                irec=i2+noc(isb)*i3+icol                                5d1s24
                bc(jtt+i3)=bc(jtt+i3)-4d0*bc(irec)                      5d1s24
               end do                                                    5d1s24
              end if                                                     5d1s24
             end do                                                      5d1s24
            end do                                                       5d1s24
           end if                                                       5d1s24
          end if                                                        5d1s24
         end if                                                         5d1s24
        end if                                                          5d1s24
       end if
      end do
      do isb=1,nsymb                                                    7d18s23
       if(idoubo(isb).gt.0)then                                         10d19s23
        nrow=((noc(isb)*(noc(isb)+1))/2)                                10d19s23
        do is=1,nsdlk                                                   10d19s23
         if(isblk(1,is).eq.isblk(2,is).and.isblk(1,is).eq.isb)then      10d19s23
          jsb=isblk(3,is)                                               10d19s23
          if(idoubo(jsb).gt.0)then                                      10d19s23
          end if                                                        10d19s23
         end if                                                         10d19s23
         if(isblk(1,is).eq.isblk(4,is).and.isblk(1,is).eq.isb)then      10d19s23
          jsb=isblk(2,is)                                               10d19s23
          if(idoubo(jsb).gt.0)then                                      10d19s23
           call ilimts(noc(jsb),noc(isb),mynprocg,mynowprog,il,ih,i1s,  10d19s23
     $           i1e,i2s,i2e)                                           10d19s23
           if(isblk(1,is).eq.isblk(2,is))then                           10d19s23
            nrow12=(noc(isb)*(noc(isb)+1))/2                            10d19s23
            iswitch=0                                                   10d19s23
           else                                                         10d19s23
            nrow12=noc(isb)*noc(jsb)                                    10d19s23
            iswitch=1                                                   10d19s23
           end if                                                       10d19s23
           nhere=ih+1-il                                                10d19s23
          end if                                                        10d19s23
         else if(isblk(1,is).ne.isblk(2,is).and.                        10d19s23
     $          isblk(1,is).eq.isblk(3,is).and.isblk(2,is).eq.isb)then  10d19s23
          jsb=isblk(1,is)                                               10d19s23
          if(idoubo(jsb).gt.0)then                                      10d19s23
           call ilimts(noc(isb),noc(jsb),mynprocg,mynowprog,il,ih,i1s,  10d19s23
     $           i1e,i2s,i2e)                                           10d19s23
           nrow12=noc(isb)*noc(jsb)                                     10d19s23
           nhere=ih+1-il                                                10d19s23
          end if                                                        10d19s23
         end if                                                         10d19s23
        end do                                                          10d19s23
        if(nvirt(isb).gt.0)then                                         7d19s23
         do is=1,nsdlk1                                                  7d18s23
          if(isblk1(1,is).eq.isblk1(2,is).and.isblk1(3,is).eq.isb)then     7d18s23
           jsb=isblk1(1,is)                                              7d18s23
           nrow=(noc(jsb)*(noc(jsb)+1))/2                               7d19s23
           if(idoubo(jsb).gt.0)then                                     7d18s23
            call ilimts(noc(isb),nvirt(isb),mynprocg,mynowprog,il,ih,   7d19s23
     $           i1s,i1e,i2s,i2e)                                       7d19s23
            do iv=0,nvirt(isb)-1                                        7d19s23
             ivp=iv+noc(isb)                                            10d19s23
             do idd=0,idoubo(isb)-1                                     7d19s23
              icol=idd+1+noc(isb)*iv                                    7d19s23
              if(icol.ge.il.and.icol.le.ih)then                          7d18s23
               icol=ionex(is)+nrow*(icol-il)-1                          7d20s23
               jtt=itt(isb)+idd+nbasdws(isb)*ivp                        10d19s23
               sum=0d0                                                  7d19s23
               do i34=1,idoubo(jsb)                                     7d19s23
                ii=(i34*(i34+1))/2                                      7d19s23
                iadd=icol+ii                                            7d19s23
                sum=sum+bc(iadd)                                        7d19s23
               end do                                                   7d18s23
               bc(jtt)=bc(jtt)+8d0*sum                                  10d19s23
              end if                                                     7d18s23
             end do                                                     7d19s23
            end do                                                      7d18s23
           end if                                                       7d18s23
          end if                                                        7d18s23
          if(isblk1(1,is).eq.isblk1(4,is).and.isblk1(1,is).eq.isb)then  7d19s23
           jsb=isblk1(2,is)                                              7d18s23
           if(idoubo(jsb).gt.0)then                                     7d18s23
            call ilimts(noc(jsb),nvirt(isb),mynprocg,mynowprog,il,ih,   7d19s23
     $           i1s,i1e,i2s,i2e)                                       7d19s23
            if(isblk1(1,is).eq.isblk1(2,is))then                        7d19s23
             nrow12=(noc(isb)*(noc(isb)+1))/2                           7d19s23
             iswitch=0                                                  7d19s23
            else                                                        7d19s23
             nrow12=noc(isb)*noc(jsb)                                   7d19s23
             iswitch=1                                                  7d19s23
            end if                                                      7d19s23
            nhere=ih+1-il
            do iv=0,nvirt(isb)-1                                        7d19s23
             ivp=iv+noc(isb)                                            10d19s23
             do i1=0,idoubo(jsb)-1                                      7d19s23
              icol=i1+1+noc(jsb)*iv                                     7d19s23
              if(icol.ge.il.and.icol.le.ih)then                          7d18s23
               icol0=icol
               icol=ionex(is)+nrow12*(icol-il)                          7d19s23
               jtt=itt(isb)+nbasdws(isb)*ivp                            10d19s23
               do idd=0,idoubo(isb)-1                                    7d19s23
                irec=idd+noc(isb)*i1                                     7d19s23
                ix=max(idd,i1)                                          7d19s23
                in=min(idd,i1)                                          7d19s23
                itri=((ix*(ix+1))/2)+in                                 7d21s23
                itri=itri+iswitch*(irec-itri)                           7d19s23
                iadd=icol+itri                                          7d19s23
                bc(jtt+idd)=bc(jtt+idd)-4d0*bc(iadd)                    10d19s23
               end do                                                   7d18s23
              end if                                                    7d19s23
             end do                                                     7d19s23
            end do                                                      7d18s23
           end if                                                       7d18s23
          else if(isblk1(2,is).eq.isblk1(4,is).and.                     7d21s23
     $         isblk1(2,is).eq.isb.and.isblk1(1,is).ne.isblk1(2,is))then7d21s23
           jsb=isblk1(1,is)                                              7d18s23
           if(idoubo(jsb).gt.0)then                                     7d18s23
            call ilimts(noc(jsb),nvirt(isb),mynprocg,mynowprog,il,ih,   7d19s23
     $           i1s,i1e,i2s,i2e)                                       7d19s23
            nrow12=noc(isb)*noc(jsb)                                    7d21s23
            nhere=ih+1-il
            do iv=0,nvirt(isb)-1                                        7d19s23
             ivp=iv+noc(isb)                                            10d19s23
             do i1=0,idoubo(jsb)-1                                      7d19s23
              icol=i1+1+noc(jsb)*iv                                     7d19s23
              if(icol.ge.il.and.icol.le.ih)then                          7d18s23
               icol0=icol
               icol=ionex(is)+nrow12*(icol-il)                          7d19s23
               jtt=itt(isb)+nbasdws(isb)*ivp                            10d19s23
               do idd=0,idoubo(isb)-1                                    7d19s23
                irec=i1+noc(jsb)*idd                                    7d21s23
                iadd=icol+irec                                          7d21s23
                bc(jtt+idd)=bc(jtt+idd)-4d0*bc(iadd)                    10d19s23
               end do                                                   7d18s23
              end if                                                    7d19s23
             end do                                                     7d19s23
            end do                                                      7d18s23
           end if                                                       7d18s23
          end if                                                        7d18s23
         end do                                                         7d18s23
        end if                                                          7d18s23
       end if                                                           7d18s23
      end do                                                            7d18s23
      return                                                            7d20s23
      end                                                               7d20s23
