c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine parseterm(line,nroot,irxinfos,nrxinfos,maxxref)        4d24s21
      character*(*) line                                                4d24s21
      dimension irxinfos(5,*),isymb(5),lz2(5)                           8d23s21
      write(6,*)('hi, my name is parseterm "'),line,('"'),              3d21s22
     $     nroot,nrxinfos,maxxref                                       3d21s22
      nlzz=0                                                            1d23s23
      lam=0                                                             2d10s23
      if(line(1:1).eq.'S')then                                          3d21s22
c                                                                       4d24s21
c     Sj (no extra symmetry), Sig+-gu (Cinfv or Dinfv), S or So (atom)    4d24s21
c                                                                       4d24s21
       if(len(line).eq.1)then                                           4d24s21
        nlzz=6                                                          4d24s21
        lam=0                                                           4d24s21
        isymb(1)=1
        lz2(1)=0                                                        8d23s21
       else if(len(line).eq.2)then                                      4d24s21
        if(line(2:2).eq.'o')then                                        4d24s21
         nlzz=6                                                          4d24s21
         lam=0                                                           4d24s21
         isymb(1)=8
         lz2(1)=0                                                       8d23s21
        else                                                            4d24s21
         nlzz=0                                                         4d24s21
         lam=0                                                          4d24s21
         read(line(2:2),*)isymb(1)                                      4d24s21
        end if                                                          4d24s21
       else                                                             4d24s21
        do iq=len(line),1,-1                                            2d10s22
         if(.not.(line(iq:iq).eq.' '.or.line(iq:iq).eq.','))then        2d18s22
          ntail=iq                                                      2d10s22
          go to 1510                                                    2d10s22
         end if                                                         2d10s22
        end do                                                          2d10s22
 1510   continue                                                        2d10s22
        nlzz=2                                                          4d24s21
        lam=0                                                           4d24s21
        igu=0                                                           4d24s21
        if(line(len(line):len(line)).eq.'g')then                        4d24s21
         ntail=ntail-1                                                  4d24s21
         igu=1                                                          4d24s21
        else if(line(len(line):len(line)).eq.'u')then                   4d24s21
         ntail=ntail-1                                                  4d24s21
         igu=-1                                                         4d24s21
        end if                                                          4d24s21
        if(line(ntail:ntail).eq.'+')then                                4d24s21
         if(igu.ge.0)then                                               4d24s21
          isymb(1)=1                                                    4d24s21
         else                                                           4d24s21
          isymb(1)=5                                                    4d24s21
         end if                                                         4d24s21
        else if(line(ntail:ntail).eq.'-')then                           4d24s21
         if(igu.ge.0)then                                               4d24s21
          isymb(1)=4                                                    4d24s21
         else                                                           4d24s21
          isymb(1)=8                                                    4d24s21
         end if                                                         4d24s21
        else                                                            4d24s21
         write(6,*)('unable to parsea '),line                            4d24s21
         stop
        end if                                                          4d24s21
       end if                                                           4d24s21
      else if(line(1:1).eq.'P')then                                     4d24s21
c
c     Pi gu or Phi gu (Cinfv or Dinfv), P or Po (atom)
c
       if(len(line).eq.1)then                                           4d24s21
        nlzz=6                                                          4d24s21
        lam=1                                                           4d24s21
        isymb(1)=4
        isymb(2)=6
        isymb(3)=7
        lz2(1)=0                                                        8d23s21
        lz2(2)=1                                                        8d23s21
        lz2(3)=1                                                        8d23s21
       else if(len(line).eq.2)then                                      4d24s21
        if(line(2:2).eq.'i')then                                        4d24s21
         nlzz=2                                                         4d24s21
         lam=1                                                          4d24s21
         isymb(1)=2                                                     4d24s21
         isymb(2)=3                                                     4d24s21
        else if(line(2:2).eq.'o')then                                   4d24s21
         nlzz=6                                                          4d24s21
         lam=1                                                           4d24s21
         isymb(1)=2
         isymb(2)=3
         isymb(3)=5
         lz2(1)=1                                                       8d23s21
         lz2(2)=1                                                       8d23s21
         lz2(3)=0                                                       8d23s21
        else                                                            4d24s21
         write(6,*)('unable to parseb '),line                            4d24s21
         stop                                                           4d24s21
        end if
       else                                                             4d24s21
        nlzz=2                                                          4d24s21
        if(line(2:2).eq.'i')then                                        4d24s21
         lam=1                                                          4d24s21
        else if(line(2:2).eq.'h')then                                   4d24s21
         lam=3                                                          4d24s21
        else                                                            4d24s21
         write(6,*)('unable to parsec '),line                            4d24s21
         stop                                                           4d24s21
        end if                                                          4d24s21
        igu=0                                                           4d24s21
        if(line(len(line):len(line)).eq.'g')then                        4d24s21
         igu=1                                                          4d24s21
        else if(line(len(line):len(line)).eq.'u')then                   4d24s21
         igu=-1                                                         4d24s21
        end if                                                          4d24s21
        if(igu.le.0)then                                                4d24s21
         isymb(1)=2                                                     4d24s21
         isymb(2)=3                                                     4d24s21
        else                                                            4d24s21
         isymb(1)=6                                                     4d24s21
         isymb(2)=7                                                     4d24s21
        end if                                                          4d24s21
       end if                                                           4d24s21
      else if(line(1:1).eq.'D')then                                     4d24s21
c
c     Del gu (Cinfv or Dinfv), D or Do (atom)
c
       lam=2                                                            4d24s21
       if(len(line).eq.1)then                                           4d24s21
        nlzz=6                                                          4d24s21
        isymb(1)=1                                                       4d24s21
        isymb(2)=1
        isymb(3)=4
        isymb(4)=6                                                       4d24s21
        isymb(5)=7                                                       4d24s21
        lz2(1)=0                                                        8d23s21
        lz2(2)=2                                                        8d23s21
        lz2(3)=2                                                        8d23s21
        lz2(4)=1                                                        8d23s21
        lz2(5)=1                                                        10d25s21
       else if(len(line).eq.2)then                                      4d24s21
        if(line(2:2).eq.'o')then                                        4d24s21
         nlzz=6                                                          4d24s21
         isymb(5)=5                                                     8d24s21
         isymb(4)=3                                                     4d24s21
         isymb(3)=2                                                     4d24s21
         isymb(2)=8                                                     5d14s21
         isymb(1)=8                                                     4d24s21
         lz2(5)=2                                                       8d23s21
         lz2(4)=1                                                       8d23s21
         lz2(3)=1                                                       8d23s21
         lz2(2)=2                                                       8d23s21
         lz2(1)=0                                                       8d23s21
        else                                                            4d24s21
         write(6,*)('unable to parsed '),line                            4d24s21
         stop                                                           4d24s21
        end if
       else                                                             4d24s21
        nlzz=2                                                          4d24s21
        igu=0                                                           4d24s21
        if(line(len(line):len(line)).eq.'g')then                        4d24s21
         igu=1                                                          4d24s21
        else if(line(len(line):len(line)).eq.'u')then                   4d24s21
         igu=-1                                                         4d24s21
        end if                                                          4d24s21
        if(igu.ge.0)then                                                4d24s21
         isymb(1)=1                                                     4d24s21
         isymb(2)=4                                                     4d24s21
        else                                                            4d24s21
         isymb(1)=5                                                     4d24s21
         isymb(2)=8                                                     4d24s21
        end if                                                          4d24s21
       end if                                                           4d24s21
      else if(line(1:1).eq.'G')then                                     4d24s21
c
c     Gam gu
c
       nlzz=2                                                           4d24s21
       lam=4                                                            4d24s21
       igu=0                                                            4d24s21
       if(line(len(line):len(line)).eq.'g')then                         4d24s21
        igu=1                                                           4d24s21
       else if(line(len(line):len(line)).eq.'u')then                    4d24s21
        igu=-1                                                          4d24s21
       end if                                                           4d24s21
       if(igu.le.0)then                                                 4d24s21
        isymb(1)=1                                                      12d30s21
        isymb(2)=4                                                      12d30s21
       else                                                             4d24s21
        isymb(1)=5                                                      12d30s21
        isymb(2)=8                                                      12d31s21
       end if
      else if(line(1:1).eq.'A')then                                     9d8s22
       nlzz=0                                                           9d8s22
       if(len(line).gt.2)then                                           1d23s23
        if(line(2:3).eq.'''''')then                                     1d23s23
         isymb(1)=2                                                      9d8s22
        else                                                            1d23s23
         write(6,*)('don''t know how to parse "'),line(1:3),('"')       1d23s23
        end if                                                          1d23s23
       else                                                             1d23s23
        if(line(2:2).eq.'"')then                                        1d23s23
         isymb(1)=2                                                      9d8s22
        else if(line(2:2).eq.'''')then                                   9d27s22
         isymb(1)=1
        else if(line(2:2).eq.'g')then                                    9d8s22
         isymb(1)=1                                                      9d8s22
        else if(line(2:2).eq.'u')then                                    9d8s22
         isymb(1)=8                                                      9d8s22
        else if(line(2:2).eq.'1')then                                    9d8s22
         isymb(1)=1                                                      9d8s22
        else if(line(2:2).eq.'2')then                                    9d8s22
         isymb(1)=4                                                      9d8s22
        else                                                             9d8s22
         write(6,*)('don''t know how to parse "'),line(1:2),('"')       1d23s23
        end if                                                          1d23s23
       end if                                                           9d8s22
      else if(line(1:1).eq.'B')then                                     9d8s22
       if(len(line).gt.2)then                                           1d23s23
        if(line(2:3).eq.'3u')then                                        9d8s22
         isymb(1)=2                                                      9d8s22
        else if(line(2:3).eq.'2u')then                                   9d8s22
         isymb(1)=3                                                      9d8s22
        else if(line(2:3).eq.'1g')then                                   9d8s22
         isymb(1)=4                                                      9d8s22
        else if(line(2:3).eq.'1u')then                                   9d8s22
         isymb(1)=5                                                      9d8s22
        else if(line(2:3).eq.'2g')then                                   9d8s22
         isymb(1)=6                                                      9d8s22
        else if(line(2:3).eq.'3g')then                                   9d8s22
         isymb(1)=7                                                      9d8s22
        else                                                             9d8s22
         write(6,*)('don''t know how to parse "'),line(1:3),('"')        9d8s22
        end if                                                          1d23s23
       else                                                             1d23s23
        if(line(2:2).eq.'1')then                                        1d23s23
         isymb(1)=3                                                      9d8s22
        else if(line(2:2).eq.'2')then                                    9d8s22
         isymb(1)=2                                                      9d8s22
        end if                                                          1d23s23
       end if                                                           9d8s22
      else                                                              9d8s22
c
c     ?
c
       write(6,*)('unable to parsee '),line                              4d24s21
       return                                                           3d10s22
      end if                                                            4d24s21
      if(nlzz.eq.0.or.lam.eq.0)then                                     4d24s21
       nhere=1                                                          4d24s21
      else if(nlzz.eq.2)then                                            4d24s21
       nhere=2                                                          4d24s21
      else                                                              4d24s21
       nhere=2*lam+1                                                    4d24s21
      end if                                                            4d24s21
      do ihere=1,nhere                                                  4d24s21
       nrxinfos=nrxinfos+1                                              4d24s21
       if(nrxinfos.ge.maxxref)then                                      4d24s21
        write(6,*)('tooooo many xref states: '),nrxinfos,maxxref,nhere  1d18s23
        stop                                                            4d24s21
       end if                                                           4d24s21
       irxinfos(1,nrxinfos)=isymb(ihere)                                4d24s21
       irxinfos(2,nrxinfos)=nroot                                       4d24s21
       irxinfos(3,nrxinfos)=nlzz                                        4d24s21
       irxinfos(4,nrxinfos)=lam                                         4d24s21
       irxinfos(5,nrxinfos)=lz2(ihere)                                  8d23s21
      end do                                                            4d24s21
      return                                                            4d24s21
      end                                                               4d24s21
