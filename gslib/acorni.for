      double precision function acorni(idum)
c-----------------------------------------------------------------------
c
c Fortran implementation of ACORN random number generator of order less
c than or equal to 12 (higher orders can be obtained by increasing the
c parameter value MAXORD).
c
c
c NOTES: 1. The variable idum is a dummy variable. The common block
c           IACO is used to transfer data into the function.
c
c        2. Before the first call to ACORN the common block IACO must
c           be initialised by the user, as follows. The values of
c           variables in the common block must not subsequently be
c           changed by the user.
c
c             KORDEI - order of generator required ( must be =< MAXORD)
c
c             MAXINT - modulus for generator, must be chosen small
c                      enough that 2*MAXINT does not overflow
c
c             ixv(1) - seed for random number generator
c                      require 0 < ixv(1) < MAXINT
c
c             (ixv(I+1),I=1,KORDEI)
c                    - KORDEI initial values for generator
c                      require 0 =< ixv(I+1) < MAXINT
c
c        3. After initialisation, each call to ACORN generates a single
c           random number between 0 and 1.
c
c        4. An example of suitable values for parameters is
c
c             KORDEI   = 10
c             MAXINT   = 2**30
c             ixv(1)   = an odd integer in the (approximate) range 
c                        (0.001 * MAXINT) to (0.999 * MAXINT)
c             ixv(I+1) = 0, I=1,KORDEI
c
c
c
c Author: R.S.Wikramaratna,                           Date: October 1990
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
      common/iaco/ ixv(MAXOP1)

      do i=1,KORDEI
            ixv(i+1)=(ixv(i+1)+ixv(i))
            if(ixv(i+1).ge.MAXINT) ixv(i+1)=ixv(i+1)-MAXINT
      end do
      acorni=dble(ixv(KORDEI+1))/MAXINT
      return
      end

      double precision function acorniopt(idum)

      implicit double precision (a-h,o-z)
      parameter (KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
      common/iaco/ ixv(MAXOP1)

      do i=1,KORDEI,6

            ixv(i+1)=(ixv(i+1)+ixv(i))
            if(ixv(i+1).ge.MAXINT) ixv(i+1)=ixv(i+1)-MAXINT

            ixv(i+2)=(ixv(i+2)+ixv(i+1))
            if(ixv(i+2).ge.MAXINT) ixv(i+2)=ixv(i+2)-MAXINT

            ixv(i+3)=(ixv(i+3)+ixv(i+2))
            if(ixv(i+3).ge.MAXINT) ixv(i+3)=ixv(i+3)-MAXINT

            ixv(i+4)=(ixv(i+4)+ixv(i+3))
            if(ixv(i+4).ge.MAXINT) ixv(i+4)=ixv(i+4)-MAXINT

            ixv(i+5)=(ixv(i+5)+ixv(i+4))
            if(ixv(i+5).ge.MAXINT) ixv(i+5)=ixv(i+5)-MAXINT

            ixv(i+6)=(ixv(i+6)+ixv(i+5))
            if(ixv(i+6).ge.MAXINT) ixv(i+6)=ixv(i+6)-MAXINT

c            ixv(i+7)=(ixv(i+7)+ixv(i+6))
c            if(ixv(i+7).ge.MAXINT) ixv(i+7)=ixv(i+7)-MAXINT
c
c            ixv(i+8)=(ixv(i+8)+ixv(i+7))
c            if(ixv(i+8).ge.MAXINT) ixv(i+8)=ixv(i+8)-MAXINT
c
c            ixv(i+9)=(ixv(i+9)+ixv(i+8))
c            if(ixv(i+9).ge.MAXINT) ixv(i+9)=ixv(i+9)-MAXINT
c
c            ixv(i+10)=(ixv(i+10)+ixv(i+9))
c            if(ixv(i+10).ge.MAXINT) ixv(i+10)=ixv(i+10)-MAXINT
c
c            ixv(i+11)=(ixv(i+11)+ixv(i+10))
c            if(ixv(i+11).ge.MAXINT) ixv(i+11)=ixv(i+11)-MAXINT
c
c            ixv(i+12)=(ixv(i+12)+ixv(i+11))
c            if(ixv(i+12).ge.MAXINT) ixv(i+12)=ixv(i+12)-MAXINT


      end do
      acorniopt=dble(ixv(KORDEI+1))/MAXINT
c      acorniopt=dble(shiftl(ixv(KORDEI+1),30))
      return

      end

c      double precision function acorniopt2(idum)
c
c      implicit double precision (a-h,o-z)
c      parameter (KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
cc      common/iaco/ ixv(MAXOP1)
c      common/iaco2/ ixv0,ixv1,ixv2,ixv3,ixv4,ixv5,ixv6,
c     +ixv7,ixv8,ixv9,ixv10,ixv11,ixv12
c
cc      do i=1,KORDEI,6
c
c            ixv1=(ixv1+ixv0)
c            if(ixv1.ge.MAXINT) ixv1=ixv1-MAXINT
c
c            ixv2=(ixv2+ixv1)
c            if(ixv2.ge.MAXINT) ixv2=ixv2-MAXINT
c
c            ixv3=(ixv3+ixv2)
c            if(ixv3.ge.MAXINT) ixv3=ixv3-MAXINT
c
c            ixv4=(ixv4+ixv3)
c            if(ixv4.ge.MAXINT) ixv4=ixv4-MAXINT
c
c            ixv5=(ixv5+ixv4)
c            if(ixv5.ge.MAXINT) ixv5=ixv5-MAXINT
c
c            ixv6=(ixv6+ixv5)
c            if(ixv6.ge.MAXINT) ixv6=ixv6-MAXINT
c
c            ixv7=(ixv7+ixv6)
c            if(ixv7.ge.MAXINT) ixv7=ixv7-MAXINT
c
c            ixv8=(ixv8+ixv7)
c            if(ixv8.ge.MAXINT) ixv8=ixv8-MAXINT
c
c            ixv9=(ixv9+ixv8)
c            if(ixv9.ge.MAXINT) ixv9=ixv9-MAXINT
c
c            ixv10=(ixv10+ixv9)
c            if(ixv10.ge.MAXINT) ixv10=ixv10-MAXINT
c
c            ixv11=(ixv11+ixv10)
c            if(ixv11.ge.MAXINT) ixv11=ixv11-MAXINT
c
c            ixv12=(ixv12+ixv11)
c            if(ixv12.ge.MAXINT) ixv12=ixv12-MAXINT
c
c
cc      end do
cc      acorniopt2=dble(ixv(KORDEI+1))/MAXINT
c      acorniopt2=dble(ixv12)/MAXINT
cc      acorniopt=dble(shiftl(ixv(KORDEI+1),30))
c      return
c
c      end

      double precision function acornilocal(MAXOP1,KORDEI,MAXINT,
     + ixv,idum)

      implicit double precision (a-h,o-z)
c      parameter (KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
c      common/iaco/ ixv(MAXOP1)
      integer MAXOP1,KORDEI,MAXINT
      integer ixv(MAXOP1)

      do i=1,KORDEI
            ixv(i+1)=(ixv(i+1)+ixv(i))
            if(ixv(i+1).ge.MAXINT) ixv(i+1)=ixv(i+1)-MAXINT
      end do
      acornilocal=dble(ixv(KORDEI+1))/MAXINT
      return
      end


