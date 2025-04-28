      subroutine randset(iseed)
      implicit real*8 (a-h,o-z)
c
c     set up data for random number generator.
c     reference: d.e. knuth, "the art of computer programming",
c     vol. 2 of "seminumerical algorithms" (addison-wesley, 1981),
c     2nd edition.
c
c     first generate a sequence of random numbers using the
c     multipicative congruential method, then shuffle it (pg. 32).
c
c     to continue a sequence, set iseed to a negitive number.
c     in this case the program will try to read in data from unit 50.
c
      dimension x(100),iseed3(4)
      common/rancom/x,iseed3
      save /rancom/
      data ib/65536/
      if(iseed.lt.0)then
       write(6,1)
    1  format(/1x,33hcontinuing random number sequence)
       write(6,2)
    2  format(/1x,20hreading from unit 50)
       read(50)iseed3,x
       write(6,3)iseed3
    3  format(/1x,24hnext seed in base 2**16:,/(1x,4i20))
      else
       is=iseed
       if(mod(is,2).eq.0)is=is-1
       do 4 i=1,4
        iseed3(i)=mod(is,ib)
        is=is/ib
    4  continue
   10  format(/1x,36hinitiating a random number sequence.,
     $        /1x,17hseed in base 10 =,i20,
     $        /1x,20hseed in base 2**16 =,4i20)
       do 5 i=1,100
        x(i)=rand1(iseed3)
    5  continue
       iseed=0
      end if
      return
      end
      subroutine randget(iseed)
      implicit real*8 (a-h,o-z)
c
c     save the date necessary to restart a random number sequence.
c     this information is saved on unit 50.
c
c     reference: d.e. knuth, "the art of computer programming",
c     vol. 2 of "seminumerical algorithms" (addison-wesley, 1981),
c     2nd edition.
c
c     first generate a sequence of random numbers using the
c     multipicative congruential method, then shuffle it (pg. 32).
c
      dimension x(100),iseed3(4)
      common/rancom/x,iseed3
      save /rancom/
      rewind 50
      write(50)iseed3,x
      write(6,3)iseed3
    3 format(1x,46hsave random no. seq. on unit 50 - next seed is,
     $       1x,4i20)
      return
      end
      function randdws(idum)
      implicit real*8 (a-h,o-z)
c
c     get the next random number.
c
c     reference: d.e. knuth, "the art of computer programming",
c     vol. 2 of "seminumerical algorithms" (addison-wesley, 1981),
c     2nd edition.
c
c     first generate a sequence of random numbers using the
c     multipicative congruential method, then shuffle it (pg. 32).
c
      dimension x(100),iseed3(4)
      common/rancom/x,iseed3
      save /rancom/
      j=int(99d0*x(100))+1
      randdws=x(100)
      x(100)=x(j)
      x(j)=rand1(iseed3)
      return
      end
      function rand1(iseed)
      implicit real*8 (a-h,o-z)
c
c     random number generatation using the multipicative congruential
c     method.
c     if x(n) is the current seed, then the next one x(n+1) is given by
c     x(n+1)=(a*x(n)+c)mod(m)
c     and the random number is x(n+1)/m.
c     this routine hits interger overflow at 2**32.  the arithmatic is
c     carried out using 4 "digit" numbers of base 2**16, thus iseed is
c     4 element array with each element a "digit". the units "digit" is
c     the first element of iseed.
c     parameters:
c     m=2**64=(0)(0)(0)(0)(1) base 2**16
c     a=6,364,136,223,846,793,005 base 10 =
c          (322557)(19605)(62509)(22609) base 2**16
c     c=0.
c     these parameters are from d.e. knuth, "the art of computer
c     programming", vol. 2 of "seminumerical algorithms" (addison-wesley
c     reading, mass. 1981) 2nd eddition, and are atributed to
c     michel lavaux and frank janssens.
c
      dimension iseed(4),ia(4),ic(4),id(8)
      data ia/32557,19605,62509,22609/
      data ic/0,0,0,0/
      data ib/65536/
      data bi/1.52587890625d-5/
      b=1d0/bi
c
c     id will equal iseed*ia+ic.
c     set it equal to ic
c
      do 1 i=1,4
       id(i)=ic(i)
       id(i+4)=0
    1 continue
c
c     form ia*iseed+ic
c
      do 2 j=1,4
c       do 2 i=1,4
       do 2 i=1,5-j
        k=j+i-1
c
c      ip is unnormalized product of k digit of result so far.
c      it should never exceed ib*(ib-1).
c
        ip=ia(j)*iseed(i)
c
c      now normalize the result (carry forword to next digits if
c      necessary)
c
    3   continue
        ip=ip+id(k)
        id(k)=mod(ip,ib)
        ip=ip/ib
c        if(ip.eq.0)go to 20
        if(ip.eq.0.or.k.eq.4)go to 20
        k=k+1
        go to 3
   20  continue
c
    2 continue
c
c     now perform modular arithmatic
c
      do 4 i=1,4
       iseed(i)=id(i)
    4 continue
c
c     now determine floating point random number
c
      rand1=dfloat(iseed(1))
      do 5 i=2,4
       rand1=dfloat(iseed(i))+rand1*bi
    5 continue
      rand1=rand1*bi
      return
      end
