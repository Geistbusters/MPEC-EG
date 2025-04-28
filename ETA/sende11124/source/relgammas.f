c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine relgammas(nec,mdon,mdoo,irw0,irw1,irw2,spinx,bc,ibc)   11d9s22
      implicit real*8 (a-h,o-z)
      include "common.store"
      integer*1 icode(64),ipackc1(8)                                    9d16s21
      integer*2 ipack28(4)                                              9d16s21
      integer*4 ipack48(2)
      integer*8 ipackc,ipack8                                           9d16s21
      equivalence (ipackc,ipackc1),(ipack8,ipack48),(ipack8,ipack28)    9d16s21
c
c     addressing ixw1. We will be given io and icase and we need to
c     point to index specifying full or compacted, then address if full
c     or no. words, then address if compacted.
c     this is essential triangle storage, so ((io*(io-1))/2)+icase      12d3s20
c     would do it.
c
      nopenx=nec-mdon*2                                                 5d8s20
      nopenn=nec-mdoo*2                                                 5d8s20
      nopenn=max(mod(nec,2),nopenn-2)                                   8d8s22
      nxx=nopenx+1                                                      5d8s20
      nxxa=(nxx*(nxx+1))/2                                              12d3s20
      loop=0
      igoal=5770914
c
c     addressing: piggy backing on non-rel code,                        9d16s21
c     iooa=ixw1+((io*(io-1))/2)+icase-1   ?                              9d16s21
c     let ixw1+io+delta point to origin, delta is 0 or 1 for normal or
c     transposed.
c     then use below to get spin address.
c     then use same for ket spins. at that address give bra&ket ncsf+
c     ntype, then compute icase offset.
c     will point to info for that io and icase.
c     then given is2b and ms2b
c     is2b ms2b (ms+in) is2b ms2b (ms+in) #          is2b orig   sum i0,in-1 2i+1 = (in-1+1-i0)(2in-
c     0    0     0       1    -1     0    1            0    1                           =(2in+2i0)*(
c     2   -2     0       1    +1     1    2            1    1
c     2    0     1       3    -3     0    3            2    2                      even:i0=0, get 2i
c     2   +2     2       3    -1     1    4            3    3                      2in       need to
c     4   -4     0       3    +1     2    5            4    5                       2 2*2/4=1
c     4   -2     1       3    +3     3    6            5    7                       4 4*4/4=4
c     4    0     2       5    -5     0    7            6   10                       6 6*6/4=9
c     4   +2     3       5    -3     1    8            7   13                      even:i1=1/2, get
c     4   +4     4       5    -1     2    9                                        2in         need
c     6   -6     0       5    +1     3   10                                        3   4*2/4=2
c     6   -4     1       5    +3     4   11                                        5   6*4/4=6
c     6   -2     2       5    +5     5   12                                        7   8*6/4=12
c     6    0     3       7    -7     0   13                                        (2ms+2in)/2
c     6   +2     4       7    -5     1   14
c
c     is2b ms2b  is2b ms2b  #         is2b  orig    sum i0,in-1 (i-i0)+1 = (in-1+1-i0)(in-1-i0+1+i0-
c     0    0      1    +1   1          0     1
c     2    0      3    +1   2          1     1                 even:i0=0, get in*(in+1)/2=2in*(2in+2
c     2   +2      3    +3   3          2     2                 2in need to add 1
c     4    0      5    +1   4          3     2                  2  2*4/8=1
c     4   +2      5    +3   5          4     4                  4  4*6/8=3
c     4   +4      5    +5   6          5     4                  6  6*8/8=6
c     6    0      7    +1   7          6     7                 odd:2i0=1, get (2in-1)*(2in+1)/8
c     6   +2      7    +3   8          7     7                  3  2*4/8=1  need to add 1.
c                                                               5  4*6/=3
c     now we need to know how to deal with ket spins. here we can reset the min
c     ntype = 1 or 2 < doesn't depend on case.
c     ms triangle rule iabs(m2sb-m2sk) .le. 2
c     is triangle rule iabs(i2sb-i2sk) .le. 2
c     next: compress, then ||
      islo=mod(nopenx,2)                                                9d16s21
      ishix=nint(spinx*2d0)                                             9d20s21
      irw0=ibcoff                                                       9d20s21
      irw1=irw0+1+(nopenx-islo)/2                                       9d20s21
      irw2=irw1+1+(nopenx-islo)/2                                         9d16s21
      ibcoff=irw2+1+(nopenx-islo)/2                                       9d16s21
      ibctop=ibcoff                                                     12d5s20
      do io=nopenn,nopenx,2                                             5d11s20
       ishi=min(io,ishix)                                               9d20s21
       iooa=irw0+(io-islo)/2                                            9d20s21
       ibc(iooa)=ibcoff                                                 9d16s21
       iooa=ibcoff                                                      9d16s21
       nspinx=ishi+2                                                    9d20s21
       nspinx=(nspinx-islo)*(nspinx-islo+2)                             9d21s21
       nspinx=nspinx/8                                                  9d21s21
       nspinx=nspinx+1                                                  9d17s21
       ibcoff=ibcoff+9*nspinx                                           9d16s21
       nx=io                                                            9d20s21
       do ms2b=islo,ishi,2                                                 9d21s21
        do ms2k=-ishi,ishi,2                                            9d17s21
         if(iabs(ms2b-ms2k).le.2)then                                   9d16s21
          idelm=(ms2b-ms2k)/2                                           9d16s21
          idelm=idelm+1                                                 9d16s21
          do is2b=iabs(ms2b),ishi,2                                        9d16s21
           ioffb=(is2b-islo)*(is2b-islo+2)                              9d21s21
           ioffb=ioffb/8                                                9d16s21
           ioffb=ioffb+((islo+ms2b)/2)                                  9d21s21
           do is2k=iabs(ms2k),ishi,2                                       9d16s21
            if(iabs(is2b-is2k).le.2)then                                9d16s21
             idels=(is2b-is2k)/2                                        9d16s21
             idels=idels+1                                              9d16s21
             ioff=iooa+idels+3*(idelm+3*ioffb)                          9d16s21
             ibc(ioff)=ibcoff                                           9d16s21
             call wfetch(io,mdon,idum,is2b,ms2b,ndetb,ncsfb,ivecb,      9d16s21
     $             iaorbb,idum,bc,ibc)                                  11d9s22
             nbb=ibcoff-ibc(ioff)                                       9d17s21
             call wfetch(io,mdon,idum,is2k,ms2k,ndetk,ncsfk,iveck,      9d16s21
     $             iaorbk,idum,bc,ibc)                                  11d9s22
             nbk=ibcoff-ibc(ioff)-nbb                                   9d17s21
             nn=ncsfb*ncsfk+1                                           9d20s21
             mxsto=nn*2*io                                              9d17s21
             iaorbbn=ibc(ioff)+mxsto                                    9d17s21
             iaorbkn=iaorbbn+nbb                                        9d17s21
             ibcoff=iaorbkn+nbk                                         9d17s21
             do i=nbk-1,0,-1                                            9d17s21
              ibc(iaorbkn+i)=ibc(iaorbk+i)                              9d17s21
             end do                                                     9d17s21
             do i=nbb-1,0,-1                                            9d17s21
              ibc(iaorbbn+i)=ibc(iaorbb+i)                              9d17s21
             end do                                                     9d17s21
             iaorbk=iaorbkn                                             9d17s21
             iaorbb=iaorbbn                                             9d17s21
             idata=ibc(ioff)                                            9d17s21
             ipack48(1)=ncsfb                                           9d16s21
             ipack48(2)=ncsfk                                           9d16s21
             ibc(idata)=ipack8                                          9d17s21
             idata=idata+1                                              9d17s21
             itype=idata                                                9d17s21
             imats=itype+1                                              9d16s21
             do i=1,io                                                        5d8s20
              icode(i)=3                                                      5d8s20
             end do                                                           5d8s20
             ibc0=ibcoff                                                9d20s21
             call cupo1rp(icode,nx,bc(ivecb),ndetb,ncsfb,ibc(iaorbb),   3d23s22
     $              io,bc(iveck),ndetk,ncsfk,ibc(iaorbk),io,iout,nnot,   9d16s21
     $              n1den,ntype,bc,ibc)                                 11d14s22
             ibc(itype)=ntype                                           9d20s21
             imat=iout+n1den                                             9d16s21
             jmats=imats                                                9d20s21
             nnbk=ncsfb*ncsfk                                           3d23s22
             loop=loop+1
             do i=0,ntype-1
              ipackc=ibc(iout+i)                                          9d16s21
              ipack48(1)=ipackc1(1)                                     9d20s21
              ipack48(2)=ipackc1(2)                                     9d20s21
              ibc(jmats)=ipack8                                         9d20s21
              jmats=jmats+1                                             9d20s21
              jmats0=jmats                                              3d23s22
              ibc(jmats)=0                                              3d23s22
              jmats=jmats+1                                             3d23s22
              jmat=imat+ncsfb*ncsfk*i                                     9d16s21
               icmp1=ibcoff
               icmp2=icmp1+ncsfb                                        3d23s22
               icmp3=icmp2+ncsfk*ncsfb
               ibcoff=icmp3+ncsfk*ncsfb
               call enough('relgammas.  1',bc,ibc)
               call cmpvcsf(ncsfb,ncsfk,ibc(icmp1),ibc(icmp2),
     $              bc(icmp3),bc(jmat),nused)                                  12d1s20
                nusedi=nused/2
                if(nusedi*2.ne.nused)nusedi=nusedi+1                      3d23s22
                nused4=ncsfb+nused+nusedi+1                             3d23s22
                if(nused4.lt.nnbk)then                                   3d23s22
                 ibc(jmats0)=1                                          3d23s22
                end if                                                  3d23s22
               ibcoff=icmp1                                              3d23s22
              if(ibc(jmats0).eq.0)then                                   3d23s22
               do ii=0,ncsfb*ncsfk-1                                     9d20s21
                bc(jmats+ii)=bc(jmat+ii)                                 9d20s21
               end do                                                    9d20s21
               jmats=jmats+ncsfb*ncsfk                                   9d20s21
              else                                                      3d23s22
               ibc(jmats)=nused                                         3d23s22
               jmats=jmats+1                                            3d23s22
               do ii=0,ncsfb+nusedi-1                                   3d23s22
                ibc(jmats+ii)=ibc(icmp1+ii)                             3d23s22
               end do                                                   3d23s22
               jmats=jmats+ncsfb+nusedi                                 3d23s22
               do ii=0,nused-1                                          3d23s22
                bc(jmats+ii)=bc(icmp3+ii)                               3d23s22
               end do                                                   3d23s22
               jmats=jmats+nused                                        3d23s22
              end if                                                    3d23s22
             end do
             imats=jmats                                                9d20s21
             ibcoff=imats                                               9d17s21
            end if                                                      9d16s21
           end do                                                        9d16s21
          end do                                                         9d16s21
         end if                                                         9d16s21
        end do                                                          9d16s21
       end do                                                           9d16s21
      end do                                                            9d16s21
      do io=nopenn,nopenx,2                                             5d11s20
       ishi=min(io,ishix)                                               9d20s21
       iooa=irw1+(io-islo)/2                                            9d20s21
       ibc(iooa)=ibcoff                                                 9d16s21
       iooa=ibcoff                                                      9d16s21
       nspinx=ishi+2                                                    9d20s21
       nspinx=(nspinx-islo)*(nspinx-islo+2)                             9d21s21
       nspinx=nspinx/8                                                  9d16s21
       nspinx=nspinx+1                                                  9d17s21
       ibcoff=ibcoff+9*nspinx                                           9d16s21
       nx=io+1
       do ms2b=islo,ishi,2                                              9d21s21
        do ms2k=-ishi,ishi,2                                            9d17s21
         if(iabs(ms2b-ms2k).le.2)then                                   9d16s21
          idelm=(ms2b-ms2k)/2                                           9d16s21
          idelm=idelm+1                                                 9d16s21
          do is2b=iabs(ms2b),ishi,2                                        9d16s21
           ioffb=(is2b-islo)*(is2b-islo+2)                              9d21s21
           ioffb=ioffb/8                                                9d16s21
           ioffb=ioffb+((mod(is2b,2)+ms2b)/2)                           9d21s21
           do is2k=iabs(ms2k),ishi,2                                       9d16s21
            if(iabs(is2b-is2k).le.2)then                                9d16s21
             idels=(is2b-is2k)/2                                        9d16s21
             idels=idels+1                                              9d16s21
             ioff=iooa+idels+3*(idelm+3*ioffb)                          9d16s21
             ibc(ioff)=ibcoff                                           9d16s21
             call wfetch(io,mdon,idum,is2b,ms2b,ndetb,ncsfb,ivecb,      9d16s21
     $             iaorbb,idum,bc,ibc)                                  11d9s22
             nbb=ibcoff-ibc(ioff)                                       9d17s21
             call wfetch(io,mdon,idum,is2k,ms2k,ndetk,ncsfk,iveck,      9d16s21
     $             iaorbk,idum,bc,ibc)                                  11d9s22
             nbk=ibcoff-ibc(ioff)-nbb                                   9d17s21
             nn=2*ncsfb*ncsfk+2                                         3d25s22
             mxsto=nn*2*io+io+1                                         3d25s22
             iaorbbn=ibc(ioff)+mxsto                                    9d17s21
             iaorbkn=iaorbbn+nbb                                        9d17s21
             ibcoff=iaorbkn+nbk                                         9d17s21
             do i=nbk-1,0,-1                                            9d17s21
              ibc(iaorbkn+i)=ibc(iaorbk+i)                              9d17s21
             end do                                                     9d17s21
             do i=nbb-1,0,-1                                            9d17s21
              ibc(iaorbbn+i)=ibc(iaorbb+i)                              9d17s21
             end do                                                     9d17s21
             iaorbk=iaorbkn                                             9d17s21
             iaorbb=iaorbbn                                             9d17s21
             idata=ibc(ioff)                                            9d17s21
             ipack48(1)=ncsfb                                           9d16s21
             ipack48(2)=ncsfk                                           9d16s21
             ibc(idata)=ipack8                                          9d17s21
             idata=idata+1                                              9d17s21
             ibc(idata)=io                                              3d25s22
             itoc=idata+1                                               3d25s22
             itocm=itoc-1                                               3d25s22
             idata=itoc+io                                              3d25s22
             itype=idata                                                9d17s21
             imats=itype+1                                              9d16s21
             do i=1,io                                                        5d8s20
              icode(i)=3                                                      5d8s20
             end do                                                           5d8s20
             icode(nx)=2                                                      5d8s20
             do icase=1,io                                                    5d8s20
              icode(icase)=1                                                  5d8s20
              ibc0=ibcoff                                               9d16s21
              call cupo1rp(icode,nx,bc(ivecb),ndetb,ncsfb,ibc(iaorbb),  3d23s22
     $              io,bc(iveck),ndetk,ncsfk,ibc(iaorbk),io,iout,nnot,   9d16s21
     $              n1den,ntype,bc,ibc)                                 11d14s22
              nnbk=ncsfb*ncsfk                                          3d25s22
              nnbk2=nnbk*2                                              3d25s22
              ibc(itype)=ntype                                          9d16s21
              imat=iout+n1den                                             9d16s21
              jmats=imats                                               9d16s21
              ibc(itocm+icase)=jmats                                    3d25s22
              do i=0,ntype-1
               ipackc=ibc(iout+i)                                          9d16s21
               ipack48(1)=ipackc1(1)                                    9d16s21
               ipack48(2)=ipackc1(2)                                    9d16s21
               ibc(jmats)=ipack8                                        9d16s21
               jmats=jmats+1                                            9d16s21
               jmats0=jmats                                             3d25s22
               ibc(jmats)=0                                             3d25s22
               jmats=jmats+1                                            9d16s21
               jmat=imat+ncsfb*ncsfk*i                                     9d16s21
               icmp1=ibcoff
               icmp2=icmp1+ncsfb                                        3d25s22
               icmp3=icmp2+ncsfk*ncsfb
               itrans=icmp3+ncsfk*ncsfb                                 3d25s22
               ibcoff=itrans+nnbk                                       3d25s22
               call enough('relgammas.  2',bc,ibc)
               call cmpvcsf(ncsfb,ncsfk,ibc(icmp1),ibc(icmp2),
     $              bc(icmp3),bc(jmat),nused)                                  12d1s20
               ibcoff=icmp1                                             3d25s22
               nusedi=nused/2
               if(nusedi*2.ne.nused)nusedi=nusedi+1                      3d23s22
               nused4=ncsfb+ncsfk+2*(nused+nusedi)                      3d25s22
               if(nused4.lt.+nnbk2)then                                 3d28s22
                ibc(jmats0)=1                                           3d25s22
                ipack48(1)=nused                                        3d25s22
                ipack48(2)=ncsfb+nused+nusedi                           3d25s22
                ibc(jmats)=ipack8                                       3d25s22
                jmats=jmats+1                                           3d25s22
               end if                                                   3d25s22
               if(ibc(jmats0).eq.0)then                                 3d25s22
                do ii=0,ncsfb*ncsfk-1                                    9d16s21
                 bc(jmats+ii)=bc(jmat+ii)                                9d16s21
                end do                                                   9d16s21
                jmats=jmats+ncsfb*ncsfk                                  9d16s21
               else                                                     3d25s22
                do ii=0,ncsfb+nusedi-1                                   3d23s22
                 ibc(jmats+ii)=ibc(icmp1+ii)                             3d23s22
                end do                                                   3d23s22
                jmats=jmats+ncsfb+nusedi                                 3d23s22
                do ii=0,nused-1                                          3d23s22
                 bc(jmats+ii)=bc(icmp3+ii)                               3d23s22
                end do                                                   3d23s22
                jmats=jmats+nused                                        3d23s22
               end if                                                   3d25s22
               if(ibc(jmats0).eq.0)then                                 3d25s22
                do k=0,ncsfk-1                                           9d16s21
                 do ib=0,ncsfb-1                                         9d16s21
                  ik=jmat+ib+ncsfb*k                                     9d16s21
                  ki=jmats+k+ncsfk*ib                                    9d16s21
                  bc(ki)=bc(ik)                                          9d16s21
                 end do                                                  9d16s21
                end do                                                   9d16s21
                jmats=jmats+ncsfb*ncsfk                                  9d16s21
               else                                                     3d25s22
                icmp1=ibcoff
                icmp2=icmp1+ncsfk                                        3d25s22
                icmp3=icmp2+ncsfk*ncsfb
                itrans=icmp3+ncsfk*ncsfb                                 3d25s22
                ibcoff=itrans+nnbk                                       3d25s22
                call enough('relgammas.  3',bc,ibc)
                do k=0,ncsfk-1                                           9d16s21
                 do ib=0,ncsfb-1                                         9d16s21
                  ik=jmat+ib+ncsfb*k                                     9d16s21
                  ki=itrans+k+ncsfk*ib                                  3d25s22
                  bc(ki)=bc(ik)                                          9d16s21
                 end do                                                  9d16s21
                end do                                                   9d16s21
                call cmpvcsf(ncsfk,ncsfb,ibc(icmp1),ibc(icmp2),         3d25s22
     $              bc(icmp3),bc(itrans),nused)                         3d25s22
                do ii=0,ncsfk+nusedi-1                                   3d23s22
                 ibc(jmats+ii)=ibc(icmp1+ii)                             3d23s22
                end do                                                   3d23s22
                jmats=jmats+ncsfk+nusedi                                 3d23s22
                do ii=0,nused-1                                          3d23s22
                 bc(jmats+ii)=bc(icmp3+ii)                               3d23s22
                end do                                                   3d23s22
                jmats=jmats+nused                                        3d23s22
               end if                                                   3d25s22
               ibcoff=icmp1                                             3d25s22
              end do
              imats=jmats                                               9d16s21
              icode(icase)=3                                                   5d8s20
              ibcoff=ibc0                                               9d17s21
             end do                                                            5d8s20
             ibcoff=imats                                               9d17s21
            end if                                                      9d16s21
           end do                                                        9d16s21
          end do                                                         9d16s21
         end if                                                         9d16s21
        end do                                                          9d16s21
       end do                                                           9d16s21
      end do                                                            9d16s21
      do io=nopenn,nopenx-1,2                                           5d11s20
       iop=io+2
       ishi=min(iop,ishix)                                              9d20s21
       iooa=irw2+(io-islo)/2                                            9d20s21
       ibc(iooa)=ibcoff                                                 9d16s21
       iooa=ibcoff                                                      9d16s21
       nspinb=ishi+2                                                    9d20s21
       nspinb=(nspinb-islo)*(nspinb-islo+2)                             9d21s21
       nspinb=nspinb/8                                                  9d21s21
       nspinb=nspinb+1                                                  9d17s21
       ibcoff=ibcoff+nspinb*9                                           9d16s21
       nx=iop                                                           5d11s20
       nxm=nx-1                                                         5d11s20
       do ms2b=islo,ishi,2                                              9d21s21
        do ms2k=-min(ishi,io),min(ishi,io),2                            9d20s21
         if(iabs(ms2b-ms2k).le.2)then                                   9d16s21
          idelm=(ms2b-ms2k)/2                                           9d16s21
          idelm=idelm+1                                                 9d16s21
          do is2b=iabs(ms2b),ishi,2                                     10d12s21
           do is2k=iabs(ms2k),min(io,ishi),2                            10d12s21
            if(iabs(is2k-is2b).le.2)then                                9d16s21
             idels=(is2b-is2k)/2                                        9d16s21
             idels=idels+1                                              9d16s21
             ioffb=(is2b-islo)*(is2b-islo+2)                            9d21s21
             ioffb=ioffb/8                                              9d21s21
             ioffb=ioffb+((mod(is2b,2)+ms2b)/2)                         9d21s21
             ioff=iooa+idels+3*(idelm+3*ioffb)                          9d16s21
             ibc(ioff)=ibcoff                                           9d16s21
             call wfetch(iop,mdon,idum,is2b,ms2b,ndetb,ncsfb,ivecb,       9d16s21
     $             iaorbb,idum,bc,ibc)                                  11d9s22
             nbb=ibcoff-ibc(ioff)                                       9d17s21
             call wfetch(io,mdon,idum,is2k,ms2k,ndetk,ncsfk,iveck,        9d16s21
     $             iaorbk,idum,bc,ibc)                                  11d9s22
             nbk=ibcoff-ibc(ioff)-nbb                                   9d17s21
             nn=2*ncsfb*ncsfk+2                                         3d25s22
             mxsto=nn*2*nxm+nxm+1                                       3d25s22
             iaorbbn=ibc(ioff)+mxsto                                    9d17s21
             iaorbkn=iaorbbn+nbb                                        9d17s21
             ibcoff=iaorbkn+nbk                                         9d17s21
             do i=nbk-1,0,-1                                            9d17s21
              ibc(iaorbkn+i)=ibc(iaorbk+i)                              9d17s21
             end do                                                     9d17s21
             do i=nbb-1,0,-1                                            9d17s21
              ibc(iaorbbn+i)=ibc(iaorbb+i)                              9d17s21
             end do                                                     9d17s21
             iaorbk=iaorbkn                                             9d17s21
             iaorbb=iaorbbn                                             9d17s21
             idata=ibc(ioff)                                            9d17s21
             ipack48(1)=ncsfb                                           9d16s21
             ipack48(2)=ncsfk                                           9d16s21
             ibc(idata)=ipack8                                          9d17s21
             idata=idata+1                                              9d17s21
             ibc(idata)=nxm                                             3d25s22
             itoc=idata+1                                               3d25s22
             itocm=itoc-1                                               3d25s22
             idata=itoc+nxm                                             3d25s22
             itype=idata                                                9d17s21
             imats=idata+1                                              9d17s21
             do i=1,nx                                                        5d11s20
              icode(i)=3                                                      5d11s20
             end do                                                           5d11s20
             icode(nx)=1                                                      5d11s20
             do icase=1,nxm                                                   5d11s20
              icode(icase)=7                                                  5d11s20
              ibctop=ibcoff                                             9d17s21
              call cupo1rp(icode,nx,bc(ivecb),ndetb,ncsfb,ibc(iaorbb),  3d23s22
     $             iop,bc(iveck),ndetk,ncsfk,ibc(iaorbk),io,iout,nnot,  9d16s21
     $             n1den,ntype,bc,ibc)                                  11d14s22
              nnbk=ncsfb*ncsfk                                          3d25s22
              nnbk2=nnbk*2                                              3d25s22
              ibc(itype)=ntype                                          9d16s21
              imat=iout+n1den                                             9d16s21
              jmats=imats                                               9d16s21
              ibc(itocm+icase)=jmats                                    3d25s22
              do i=0,ntype-1
               ipackc=ibc(iout+i)                                         9d16s21
               ipack48(1)=ipackc1(1)                                    9d16s21
               ipack48(2)=ipackc1(2)                                    9d16s21
               ibc(jmats)=ipack8                                        9d16s21
               jmats=jmats+1                                            9d16s21
               jmats0=jmats                                             3d23s22
               ibc(jmats0)=0                                            3d23s22
               jmats=jmats+1                                            3d25s22
               jmat=imat+ncsfb*ncsfk*i                                    9d16s21
               icmp1=ibcoff
               icmp2=icmp1+ncsfb                                        3d23s22
               icmp3=icmp2+ncsfk*ncsfb
               itrans=icmp3+ncsfk*ncsfb                                 3d25s22
               ibcoff=itrans+nnbk                                       3d25s22
               call enough('relgammas.  4',bc,ibc)
               call cmpvcsf(ncsfb,ncsfk,ibc(icmp1),ibc(icmp2),
     $              bc(icmp3),bc(jmat),nused)                                  12d1s20
               ibcoff=icmp1                                             3d25s22
               nusedi=nused/2
               if(nusedi*2.ne.nused)nusedi=nusedi+1                      3d23s22
               nused4=ncsfb+ncsfk+2*(nused+nusedi)                      3d25s22
               if(nused4.lt.+nnbk2)then                                   3d23s22
                ibc(jmats0)=1                                           3d25s22
                ipack48(1)=nused                                        3d25s22
                ipack48(2)=ncsfb+nused+nusedi                           3d25s22
                ibc(jmats)=ipack8                                       3d25s22
                jmats=jmats+1                                           3d25s22
               end if                                                   3d25s22
               if(ibc(jmats0).eq.0)then                                 3d25s22
                do ii=0,ncsfb*ncsfk-1                                    9d16s21
                 bc(jmats+ii)=bc(jmat+ii)                                9d16s21
                end do                                                   9d16s21
                jmats=jmats+ncsfb*ncsfk                                  9d16s21
               else                                                     3d25s22
                do ii=0,ncsfb+nusedi-1                                   3d23s22
                 ibc(jmats+ii)=ibc(icmp1+ii)                             3d23s22
                end do                                                   3d23s22
                jmats=jmats+ncsfb+nusedi                                 3d23s22
                do ii=0,nused-1                                          3d23s22
                 bc(jmats+ii)=bc(icmp3+ii)                               3d23s22
                end do                                                   3d23s22
                jmats=jmats+nused                                        3d23s22
               end if                                                   3d25s22
               if(ibc(jmats0).eq.0)then                                 3d25s22
                do k=0,ncsfk-1                                           9d16s21
                 do ib=0,ncsfb-1                                         9d16s21
                  ik=jmat+ib+ncsfb*k                                     9d16s21
                  ki=jmats+k+ncsfk*ib                                    9d16s21
                  bc(ki)=bc(ik)                                          9d16s21
                 end do                                                  9d16s21
                end do                                                   9d16s21
                jmats=jmats+ncsfb*ncsfk                                  9d16s21
               else                                                     3d25s22
                icmp1=ibcoff
                icmp2=icmp1+ncsfk                                       3d25s22
                icmp3=icmp2+ncsfk*ncsfb
                itrans=icmp3+ncsfk*ncsfb                                 3d25s22
                ibcoff=itrans+nnbk                                       3d25s22
                call enough('relgammas.  5',bc,ibc)
                do k=0,ncsfk-1                                           9d16s21
                 do ib=0,ncsfb-1                                         9d16s21
                  ik=jmat+ib+ncsfb*k                                     9d16s21
                  ki=itrans+k+ncsfk*ib                                  3d25s22
                  bc(ki)=bc(ik)                                          9d16s21
                 end do                                                  9d16s21
                end do                                                   9d16s21
                call cmpvcsf(ncsfk,ncsfb,ibc(icmp1),ibc(icmp2),         3d25s22
     $              bc(icmp3),bc(itrans),nused)                         3d25s22
                do ii=0,ncsfk+nusedi-1                                   3d23s22
                 ibc(jmats+ii)=ibc(icmp1+ii)                             3d23s22
                end do                                                   3d23s22
                jmats=jmats+ncsfk+nusedi                                 3d23s22
                do ii=0,nused-1                                          3d23s22
                 bc(jmats+ii)=bc(icmp3+ii)                               3d23s22
                end do                                                   3d23s22
                jmats=jmats+nused                                        3d23s22
               end if                                                   3d25s22
               ibcoff=icmp1                                             3d25s22
              end do
              ibcoff=ibctop                                             9d17s21
              imats=jmats                                               9d16s21
              icode(icase)=3                                              9d16s21
             end do                                                       9d16s21
             ibcoff=imats                                               9d17s21
            end if                                                      9d16s21
           end do                                                        9d16s21
          end do                                                         9d16s21
         end if                                                         9d16s21
        end do                                                          9d16s21
       end do                                                           9d16s21
      end do
      return                                                            9d16s21
      end                                                               9d16s21
