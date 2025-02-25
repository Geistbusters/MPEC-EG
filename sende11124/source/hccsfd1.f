c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hccsfd1(vecx,ncsft,ibasis,ncsf,nfcn,iden,              3d3s21
     $     nrootz,mdon,mysym,iptrbit,ixw1,ixw2,mdoo,nh0av,norbu,ismu,   3d18s22
     $     irelu,bc,ibc)                                                11d10s22
      implicit real*8 (a-h,o-z)                                         7d11s19
c
c     1-e density
c
      external second                                                   8d1s19
      logical ltest,lnew                                                11d16s20
      integer*1 nab1(2),nab2(2)                                         3d3s21
      integer*8 ipack,itesta,itestb,gandcc,gandco,gandcb                2d6s23
      integer*2 ipack2(4)                                               12d1s19
      equivalence (ipack,ipack2)                                        12d1s19
      dimension vecx(ncsft,nrootz),iden(*),                             3d3s21
     $ ibasis(3,*),ncsf(*),iptrbit(2,mdoo+1,*),nab4(2,3),nh0av(*),      3d18s22
     $     ismu(*),irelu(*)                                             3d18s22
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      include "common.store"                                            7d11s19
      include "common.mrci"                                             6d19s19
      ltest=.false.                                                     11d16s20
      lnew=.true.                                                       11d16s20
      if(ltest)write(6,*)('hi, my name is hccsfd1')
      ngcode=0
      ips=1                                                             3d12s21
      loopit=0                                                          7d11s19
      isum=ibcoff                                                       3d3s21
      ibcoff=isum+nrootz                                                3d3s21
      call enough('hccsfd1.  1',bc,ibc)
      sumx=0d0                                                          3d3s21
      do if=1,nfcn                                                      7d11s19
       nclo=ibasis(1,if)                                                7d11s19
       nclop=nclo+1                                                     7d11s19
       iarg=nclop-mdon                                                  7d11s19
       iic=iptrbit(1,nclop,mysym)+ibasis(2,if)-1                        11d16s20
       iio=iptrbit(2,nclop,mysym)+ibasis(3,if)-1                        11d16s20
       jps=ips                                                          12d11s19
       nopen=popcnt(ibc(iio))                                           3d3s21
       do jf=if,nfcn                                                    12d9s19
        ncloj=ibasis(1,jf)                                              7d11s19
        if(iabs(ncloj-nclo).le.1)then                                   12d9s19
         nclopj=ncloj+1                                                  6d11s19
         jarg=ncloj+1-mdon                                               6d11s19
         jp=jps+1
         jjo=iptrbit(2,nclopj,mysym)+ibasis(3,jf)-1                     11d16s20
         jjc=iptrbit(1,nclopj,mysym)+ibasis(2,jf)-1                     11d16s20
         nopenj=popcnt(ibc(jjo))                                        3d3s21
         gandcc=ieor(ibc(jjc),ibc(iic))                                 2d6s23
         gandco=ieor(ibc(jjo),ibc(iio))                                 2d6s23
         gandcb=ior(gandcc,gandco)                                      2d6s23
         ndifb=popcnt(gandcb)                                           2d6s23
         if(ndifb.le.2)then                                             2d6s23
          ndifs=popcnt(gandco)                                          2d6s23
          ndifd=popcnt(gandcc)                                          2d6s23
          if(mod(loopit,mynprocg).eq.mynowprog)then                      12d9s19
           loopit=mynowprog                                             9d10s24
           ibcsav=ibcoff                                                  12d9s19
           if(ndifd.eq.0.and.ndifs.eq.0)then                            2d6s23
            do ir=1,nrootz                                               3d3s21
             irm=ir-1                                                    3d3s21
             sum=0d0                                                     3d3s21
             do ii=0,ncsf(iarg)-1                                       3d12s21
              sum=sum+vecx(ips+ii,ir)**2                                 3d3s21
             end do                                                      3d3s21
             sumx=sumx+sum                                               3d3s21
             bc(isum+irm)=sum                                            3d3s21
            end do                                                       3d3s21
            do i=1,norbu                                                 3d3s21
             fact=0d0                                                   3d3s21
             if(btest(ibc(jjc),i))then                                  3d3s21
              fact=2d0                                                  3d3s21
             else if(btest(ibc(jjo),i))then                             3d3s21
              fact=1d0                                                  3d3s21
             end if                                                     3d3s21
             if(fact.ne.0d0)then                                        3d3s21
              jbs=ismu(i)                                                3d3s21
              jbg=irelu(i)-1                                             3d3s21
              jden=iden(jbs)+jbg*(1+nh0av(jbs))                         3d3s21
              do ir=0,nrootz-1                                          3d3s21
               sum=bc(isum+ir)                                          3d3s21
               bc(jden)=bc(jden)+sum*fact                               3d3s21
               jden=jden+nh0av(jbs)*nh0av(jbs)                          3d25s21
              end do                                                    3d3s21
             end if                                                     3d3s21
            end do                                                      3d3s21
           else if(ndifs.eq.2.and.ndifb.eq.2)then                       2d6s23
            do i=1,norb                                                 2d6s23
             if(btest(gandco,i))then                                    2d6s23
              if((btest(ibc(jjo),i).and..not.btest(ibc(iic),i)).or.     2d6s23
     $            (btest(ibc(jjc),i).and.btest(ibc(iio),i)))then        2d6s23
               nab4(1,1)=i                                              2d6s23
              else                                                      2d6s23
               nab4(2,1)=i                                              2d6s23
              end if                                                    2d6s23
             end if                                                     2d6s23
            end do                                                      2d6s23
            call gandc(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,      3d3s21
     $            nopen,jarg,iarg,ncsf,norbu,ixw1,ixw2,nnot1,nab1,iwpb1, 11d16s20
     $            iwpk1,ncsfmid1,bc,ibc)                                11d14s22
            jprod=ibcoff                                                  11d13s20
            ibcoff=jprod+ncsf(jarg)*nrootz                               3d3s21
            call enough('hccsfd1.  2',bc,ibc)
            call xtimesn(ncsf(jarg),ncsf(iarg),ncsfmid1,nrootz,iwpb1,   3d3s21
     $            iwpk1,vecx(ips,1),ncsft,bc(jprod),ncsf(jarg),1d0,0d0, 11d15s22
     $           bc,ibc)                                                11d15s22
            jbs=ismu(nab1(1))                                            3d3s21
            jbg=irelu(nab1(1))-1                                         3d3s21
            jkg=irelu(nab1(2))-1                                         3d3s21
            jden=iden(jbs)+jbg+nh0av(jbs)*jkg                           3d3s21
            kden=iden(jbs)+jkg+nh0av(jbs)*jbg                           3d3s21
            do ir=1,nrootz                                              3d3s21
             jprodx=jprod+ncsf(jarg)*(ir-1)                             3d3s21
             sum=0d0                                                    3d3s21
             do j=0,ncsf(jarg)-1                                        3d3s21
              sum=sum+vecx(jps+j,ir)*bc(jprodx+j)                       3d4s21
             end do                                                     3d3s21
             bc(jden)=bc(jden)+sum                                      3d3s21
             bc(kden)=bc(kden)+sum                                      3d3s21
             jden=jden+nh0av(jbs)*nh0av(jbs)                            3d3s21
             kden=kden+nh0av(jbs)*nh0av(jbs)                            3d3s21
            end do                                                      3d3s21
            ibcoff=jprod                                                3d3s21
           end if
          end if                                                         12d9s19
          loopit=loopit+1                                                12d9s19
         end if                                                          12d9s19
        end if                                                          12d9s19
        jps=jps+ncsf(jarg)
       end do                                                           6d11s19
       ips=ips+ncsf(iarg)                                               3d3s21
      end do                                                            3d3s21
      ibcoff=isum                                                       3d3s21
      return
      end                                                               7d11s19
