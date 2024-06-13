      subroutine cova3(x1,y1,z1,x2,y2,z2,ivarg,nst,MAXNST,c0,it,cc,aa,
     +                 irot,MAXROT,rotmat,cmax,cova)
c-----------------------------------------------------------------------
c
c                    Covariance Between Two Points
c                    *****************************
c
c This subroutine calculated the covariance associated with a variogram
c model specified by a nugget effect and nested varigoram structures.
c The anisotropy definition can be different for each nested structure.
c
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         coordinates of first point
c   x2,y2,z2         coordinates of second point
c   nst(ivarg)       number of nested structures (maximum of 4)
c   ivarg            variogram number (set to 1 unless doing cokriging
c                       or indicator kriging)
c   MAXNST           size of variogram parameter arrays
c   c0(ivarg)        isotropic nugget constant
c   it(i)            type of each nested structure:
c                      1. spherical model of range a;
c                      2. exponential model of parameter a;
c                           i.e. practical range is 3a
c                      3. gaussian model of parameter a;
c                           i.e. practical range is a*sqrt(3)
c                      4. power model of power a (a must be gt. 0  and
c                           lt. 2).  if linear model, a=1,c=slope.
c                      5. hole effect model
c   cc(i)            multiplicative factor of each nested structure.
c                      (sill-c0) for spherical, exponential,and gaussian
c                      slope for linear model.
c   aa(i)            parameter "a" of each nested structure.
c   irot             index of the rotation matrix for the first nested 
c                    structure (the second nested structure will use
c                    irot+1, the third irot+2, and so on)
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c
c
c OUTPUT VARIABLES:
c
c   cmax             maximum covariance
c   cova             covariance between (x1,y1,z1) and (x2,y2,z2)
c
c
c
c EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
c                      rotmat    computes rotation matrix for distance
c-----------------------------------------------------------------------
      parameter(PI=3.14159265,PMX=999.,EPSLON=1.e-5)
      integer   nst(*),it(*)
      real      c0(*),cc(*),aa(*)
      real*8    rotmat(MAXROT,3,3),hsqd,sqdist

c      real haainv1,haainv2,haainv3,h1,h2,h3,hsqd1,hsqd2,hsqd3



c
c Calculate the maximum covariance value (used for zero distances and
c for power model covariance):
c


      istart = 1 + (ivarg-1)*MAXNST
      cmax   = c0(ivarg)
      do is=1,nst(ivarg)
            ist = istart + is - 1
            if(it(ist).eq.4) then
                  cmax = cmax + PMX
            else
                  cmax = cmax + cc(ist)
            endif
      end do
c
c Check for "zero" distance, return with cmax if so:
c
      hsqd = sqdist(x1,y1,z1,x2,y2,z2,irot,MAXROT,rotmat)
      if(real(hsqd).lt.EPSLON) then
            cova = cmax
            return
      endif
c
c Loop over all the structures:
c
      cova = 0.0
      do is=1,nst(ivarg)
            ist = istart + is - 1
c
c Compute the appropriate distance:
c
            if(ist.ne.1) then
                  ir = min((irot+is-1),MAXROT)
                  hsqd=sqdist(x1,y1,z1,x2,y2,z2,ir,MAXROT,rotmat)
            end if
            h = real(dsqrt(hsqd))
c
c Spherical Variogram Model?
c
            if(it(ist).eq.1) then
                  hr = h/aa(ist)
                  if(hr.lt.1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
c
c Exponential Variogram Model?
c
            else if(it(ist).eq.2) then
                  cova = cova + cc(ist)*exp(-3.0*h/aa(ist))
c                  cova = cova + cc(ist)*exp(-3.0*h*aa(ist))

c
c Gaussian Variogram Model?
c
            else if(it(ist).eq.3) then
                  cova = cova + cc(ist)*exp(-3.*(h/aa(ist))*(h/aa(ist)))
c                  cova = cova + cc(ist)*exp(-3.*(h*aainv(ist))*(h*aa(ist)))

c
c Power Variogram Model?
c
            else if(it(ist).eq.4) then
                  cova = cova + cmax - cc(ist)*(h**aa(ist))
c
c Hole Effect Model?
c
            else if(it(ist).eq.5) then
c                 d = 10.0 * aa(ist)
c                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
                  cova = cova + cc(ist)*cos(h/aa(ist)*PI)
            endif
      end do

c
c Finished:
c
      return
      end

      subroutine cova3O1(x1,y1,z1,x2,y2,z2,ivarg,nst,MAXNST,c0,it,cc,aa,
     +                 irot,MAXROT,rotmat,cmax,cova)
c-----------------------------------------------------------------------
c
c                    Covariance Between Two Points
c                    *****************************
c
c This subroutine calculated the covariance associated with a variogram
c model specified by a nugget effect and nested varigoram structures.
c The anisotropy definition can be different for each nested structure.
c
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         coordinates of first point
c   x2,y2,z2         coordinates of second point
c   nst(ivarg)       number of nested structures (maximum of 4)
c   ivarg            variogram number (set to 1 unless doing cokriging
c                       or indicator kriging)
c   MAXNST           size of variogram parameter arrays
c   c0(ivarg)        isotropic nugget constant
c   it(i)            type of each nested structure:
c                      1. spherical model of range a;
c                      2. exponential model of parameter a;
c                           i.e. practical range is 3a
c                      3. gaussian model of parameter a;
c                           i.e. practical range is a*sqrt(3)
c                      4. power model of power a (a must be gt. 0  and
c                           lt. 2).  if linear model, a=1,c=slope.
c                      5. hole effect model
c   cc(i)            multiplicative factor of each nested structure.
c                      (sill-c0) for spherical, exponential,and gaussian
c                      slope for linear model.
c   aa(i)            parameter "a" of each nested structure.
c   irot             index of the rotation matrix for the first nested 
c                    structure (the second nested structure will use
c                    irot+1, the third irot+2, and so on)
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c
c
c OUTPUT VARIABLES:
c
c   cmax             maximum covariance
c   cova             covariance between (x1,y1,z1) and (x2,y2,z2)
c
c
c
c EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
c                      rotmat    computes rotation matrix for distance
c-----------------------------------------------------------------------
      parameter(PI=3.14159265,PMX=999.,EPSLON=1.e-5)
      integer   nst(*),it(*)
      real      c0(*),cc(*),aa(*)
      real*8    rotmat(MAXROT,3,3),hsqd,sqdist

c      real haainv1,haainv2,haainv3,h1,h2,h3,hsqd1,hsqd2,hsqd3



c
c Calculate the maximum covariance value (used for zero distances and
c for power model covariance):
c


      istart = 1 + (ivarg-1)*MAXNST
      cmax   = c0(ivarg)
      do is=1,nst(ivarg)
            ist = istart + is - 1
            if(it(ist).eq.4) then
                  cmax = cmax + PMX
            else
                  cmax = cmax + cc(ist)
            endif
      end do
c
c Check for "zero" distance, return with cmax if so:
c
      hsqd = sqdist(x1,y1,z1,x2,y2,z2,irot,MAXROT,rotmat)
      if(real(hsqd).lt.EPSLON) then
            cova = cmax
            return
      endif
c
c Loop over all the structures:
c
      cova = 0.0
c      do is=1,nst(ivarg)
            is=1
            ist = istart + is - 1
c
c Compute the appropriate distance:
c
c            if(ist.ne.1) then
c                  ir = min((irot+is-1),MAXROT)
c                  hsqd=sqdist(x1,y1,z1,x2,y2,z2,ir,MAXROT,rotmat)
c            end if
            h = real(dsqrt(hsqd))
c
c Spherical Variogram Model?
c
c            if(it(ist).eq.1) then
                  hr = h/aa(ist)
                  if(hr.lt.1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
c
c Exponential Variogram Model?
c


            is=2
            ist = istart + is - 1
c
c Compute the appropriate distance:
c
c            if(ist.ne.1) then
                  ir = min((irot+is-1),MAXROT)
                  hsqd=sqdist(x1,y1,z1,x2,y2,z2,ir,MAXROT,rotmat)
c            end if
            h = real(dsqrt(hsqd))

c            else if(it(ist).eq.2) then
                  cova = cova + cc(ist)*exp(-3.0*h/aa(ist))
c                  cova = cova + cc(ist)*exp(-3.0*h*aa(ist))

cc
cc Gaussian Variogram Model?
cc
c            else if(it(ist).eq.3) then
c                  cova = cova + cc(ist)*exp(-3.*(h/aa(ist))*(h/aa(ist)))
cc                  cova = cova + cc(ist)*exp(-3.*(h*aainv(ist))*(h*aa(ist)))
c
cc
cc Power Variogram Model?
cc
c            else if(it(ist).eq.4) then
c                  cova = cova + cmax - cc(ist)*(h**aa(ist))
cc
cc Hole Effect Model?
cc
c            else if(it(ist).eq.5) then
cc                 d = 10.0 * aa(ist)
cc                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
c                  cova = cova + cc(ist)*cos(h/aa(ist)*PI)
c            endif
c      end do

c
c Finished:
c
      return
      end








      subroutine covaopt3(x1,y1,z1,x2,y2,z2,ivarg,nst,MAXNST,c0,it,cc,
     +                 aa,aainv,
     +                 irot,MAXROT,rotmat,cmax,cova)
c-----------------------------------------------------------------------
c
c                    Covariance Between Two Points
c                    *****************************
c
c This subroutine calculated the covariance associated with a variogram
c model specified by a nugget effect and nested varigoram structures.
c The anisotropy definition can be different for each nested structure.
c
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         coordinates of first point
c   x2,y2,z2         coordinates of second point
c   nst(ivarg)       number of nested structures (maximum of 4)
c   ivarg            variogram number (set to 1 unless doing cokriging
c                       or indicator kriging)
c   MAXNST           size of variogram parameter arrays
c   c0(ivarg)        isotropic nugget constant
c   it(i)            type of each nested structure:
c                      1. spherical model of range a;
c                      2. exponential model of parameter a;
c                           i.e. practical range is 3a
c                      3. gaussian model of parameter a;
c                           i.e. practical range is a*sqrt(3)
c                      4. power model of power a (a must be gt. 0  and
c                           lt. 2).  if linear model, a=1,c=slope.
c                      5. hole effect model
c   cc(i)            multiplicative factor of each nested structure.
c                      (sill-c0) for spherical, exponential,and gaussian
c                      slope for linear model.
c   aa(i)            parameter "a" of each nested structure.
c   irot             index of the rotation matrix for the first nested 
c                    structure (the second nested structure will use
c                    irot+1, the third irot+2, and so on)
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c
c
c OUTPUT VARIABLES:
c
c   cmax             maximum covariance
c   cova             covariance between (x1,y1,z1) and (x2,y2,z2)
c
c
c
c EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
c                      rotmat    computes rotation matrix for distance
c-----------------------------------------------------------------------
      parameter(PI=3.14159265,PMX=999.,EPSLON=1.e-5)
      integer   nst(*),it(*)
      real      c0(*),cc(*),aa(*),aainv(*)
      real*8    rotmat(MAXROT,3,3),hsqd,sqdist

c      real haainv1,haainv2,haainv3,h1,h2,h3,hsqd1,hsqd2,hsqd3



c
c Calculate the maximum covariance value (used for zero distances and
c for power model covariance):
c


      istart = 1 + (ivarg-1)*MAXNST
      cmax   = c0(ivarg)
      do is=1,nst(ivarg)
            ist = istart + is - 1
            if(it(ist).eq.4) then
                  cmax = cmax + PMX
            else
                  cmax = cmax + cc(ist)
            endif
      end do
c
c Check for "zero" distance, return with cmax if so:
c
      hsqd = sqdist(x1,y1,z1,x2,y2,z2,irot,MAXROT,rotmat)
      if(real(hsqd).lt.EPSLON) then
            cova = cmax
            return
      endif
c
c Loop over all the structures:
c
      cova = 0.0
      do is=1,nst(ivarg)
            ist = istart + is - 1
c
c Compute the appropriate distance:
c
            if(ist.ne.1) then
                  ir = min((irot+is-1),MAXROT)
                  hsqd=sqdist(x1,y1,z1,x2,y2,z2,ir,MAXROT,rotmat)
            end if
c            h = real(dsqrt(hsqd))

c
c Spherical Variogram Model?
c
            if(it(ist).eq.1) then
c                  hr = h/aa(ist)
                  hr = (real(dsqrt(hsqd)))*aainv(ist)
                  if(hr.lt.1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
c
c Exponential Variogram Model?
c
            else if(it(ist).eq.2) then
                  cova = cova + cc(ist)*exp(-3.0*(real(dsqrt(hsqd)))*
     +                                      aainv(ist))
c                  cova = cova + cc(ist)*exp(-3.0*h/aa(ist))

c
c Gaussian Variogram Model?
c
            else if(it(ist).eq.3) then
                  cova = cova + cc(ist)*exp(-3.*(real(hsqd)*aainv(ist)*
     + aainv(ist)))
c                  cova = cova + cc(ist)*exp(-3.*(h/aa(ist))*(h/aa(ist)))

c
c Power Variogram Model?
c
            else if(it(ist).eq.4) then
                  cova = cova + cmax - cc(ist)*((real(dsqrt(hsqd)))**
     +                                            aa(ist))
c                  cova = cova + cmax - cc(ist)*(h**aa(ist))

c
c Hole Effect Model?
c
            else if(it(ist).eq.5) then
c                 d = 10.0 * aa(ist)
c                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
c                  cova = cova + cc(ist)*cos(h/aa(ist)*PI)
                  cova = cova + cc(ist)*cos((real(dsqrt(hsqd)))*
     +                                      aainv(ist)*PI)
            endif
      end do

c
c Finished:
c
      return
      end

      subroutine covaopt4(x1,y1,z1,x2,y2,z2,ivarg,nst,MAXNST,c0,it,cc,
     +                 aa,aainv,
     +                 irot,MAXROT,rotmat,cmax,cova)
c-----------------------------------------------------------------------
c
c                    Covariance Between Two Points
c                    *****************************
c
c This subroutine calculated the covariance associated with a variogram
c model specified by a nugget effect and nested varigoram structures.
c The anisotropy definition can be different for each nested structure.
c
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         coordinates of first point
c   x2,y2,z2         coordinates of second point
c   nst(ivarg)       number of nested structures (maximum of 4)
c   ivarg            variogram number (set to 1 unless doing cokriging
c                       or indicator kriging)
c   MAXNST           size of variogram parameter arrays
c   c0(ivarg)        isotropic nugget constant
c   it(i)            type of each nested structure:
c                      1. spherical model of range a;
c                      2. exponential model of parameter a;
c                           i.e. practical range is 3a
c                      3. gaussian model of parameter a;
c                           i.e. practical range is a*sqrt(3)
c                      4. power model of power a (a must be gt. 0  and
c                           lt. 2).  if linear model, a=1,c=slope.
c                      5. hole effect model
c   cc(i)            multiplicative factor of each nested structure.
c                      (sill-c0) for spherical, exponential,and gaussian
c                      slope for linear model.
c   aa(i)            parameter "a" of each nested structure.
c   irot             index of the rotation matrix for the first nested 
c                    structure (the second nested structure will use
c                    irot+1, the third irot+2, and so on)
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c
c
c OUTPUT VARIABLES:
c
c   cmax             maximum covariance
c   cova             covariance between (x1,y1,z1) and (x2,y2,z2)
c
c
c
c EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
c                      rotmat    computes rotation matrix for distance
c-----------------------------------------------------------------------

      implicit none
      real, intent(in) :: x1,y1,z1,x2,y2,z2
      integer,intent(in) :: ivarg
      integer,intent(in) :: nst(ivarg)
      integer,intent(in) :: MAXNST
      real,intent(in) :: c0(ivarg)
      integer,intent(in) :: it(MAXNST)
      real,intent(in) :: cc(MAXNST),aa(MAXNST),aainv(MAXNST)
      integer,intent(in)::irot,MAXROT
      real*8,intent(in) :: rotmat(MAXROT,3,3)
      real,intent(inout) :: cmax,cova


c      parameter(PI=3.14159265,PMX=999.,EPSLON=1.e-5)
c      integer   nst(*),it(*)
c      real      c0(*),cc(*),aa(*),aainv(*)
c      real*8    rotmat(MAXROT,3,3),hsqd,sqdist

      real PI, EPSLON
      integer PMX
      real*8 hsqd,sqdist
      integer istart,is,ist,ir
      real hr


c      real haainv1,haainv2,haainv3,h1,h2,h3,hsqd1,hsqd2,hsqd3

      PI=3.14159265
      PMX=999
      EPSLON=1.e-5


c
c Calculate the maximum covariance value (used for zero distances and
c for power model covariance):
c


      istart = 1 + (ivarg-1)*MAXNST
      cmax   = c0(ivarg)
      do is=1,nst(ivarg)
            ist = istart + is - 1
            if(it(ist).eq.4) then
                  cmax = cmax + PMX
            else
                  cmax = cmax + cc(ist)
            endif
      end do
c
c Check for "zero" distance, return with cmax if so:
c
      hsqd = sqdist(x1,y1,z1,x2,y2,z2,irot,MAXROT,rotmat)
      if(real(hsqd).lt.EPSLON) then
            cova = cmax
            return
      endif
c
c Loop over all the structures:
c
      cova = 0.0
      do is=1,nst(ivarg)
            ist = istart + is - 1
c
c Compute the appropriate distance:
c
            if(ist.ne.1) then
                  ir = min((irot+is-1),MAXROT)
                  hsqd=sqdist(x1,y1,z1,x2,y2,z2,ir,MAXROT,rotmat)
            end if
c            h = real(dsqrt(hsqd))

c
c Spherical Variogram Model?
c
            if(it(ist).eq.1) then
c                  hr = h/aa(ist)
                  hr = (real(dsqrt(hsqd)))*aainv(ist)
                  if(hr.lt.1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
c
c Exponential Variogram Model?
c
            else if(it(ist).eq.2) then
                  cova = cova + cc(ist)*exp(-3.0*(real(dsqrt(hsqd)))*
     +                                      aainv(ist))
c                  cova = cova + cc(ist)*exp(-3.0*h/aa(ist))

c
c Gaussian Variogram Model?
c
            else if(it(ist).eq.3) then
                  cova = cova + cc(ist)*exp(-3.*(real(hsqd)*aainv(ist)*
     + aainv(ist)))
c                  cova = cova + cc(ist)*exp(-3.*(h/aa(ist))*(h/aa(ist)))

c
c Power Variogram Model?
c
            else if(it(ist).eq.4) then
                  cova = cova + cmax - cc(ist)*((real(dsqrt(hsqd)))**
     +                                            aa(ist))
c                  cova = cova + cmax - cc(ist)*(h**aa(ist))

c
c Hole Effect Model?
c
            else if(it(ist).eq.5) then
c                 d = 10.0 * aa(ist)
c                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
c                  cova = cova + cc(ist)*cos(h/aa(ist)*PI)
                  cova = cova + cc(ist)*cos((real(dsqrt(hsqd)))*
     +                                      aainv(ist)*PI)
            endif
      end do

c
c Finished:
c
      return
      end

      subroutine covaopt5(sqdistance,ivarg,nst,MAXNST,c0,it,cc,
     +                 aa,aainv,
     +                 irot,MAXROT,rotmat,cmax,cova)
c-----------------------------------------------------------------------
c
c                    Covariance Between Two Points
c                    *****************************
c
c This subroutine calculated the covariance associated with a variogram
c model specified by a nugget effect and nested varigoram structures.
c The anisotropy definition can be different for each nested structure.
c
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         coordinates of first point
c   x2,y2,z2         coordinates of second point
c   nst(ivarg)       number of nested structures (maximum of 4)
c   ivarg            variogram number (set to 1 unless doing cokriging
c                       or indicator kriging)
c   MAXNST           size of variogram parameter arrays
c   c0(ivarg)        isotropic nugget constant
c   it(i)            type of each nested structure:
c                      1. spherical model of range a;
c                      2. exponential model of parameter a;
c                           i.e. practical range is 3a
c                      3. gaussian model of parameter a;
c                           i.e. practical range is a*sqrt(3)
c                      4. power model of power a (a must be gt. 0  and
c                           lt. 2).  if linear model, a=1,c=slope.
c                      5. hole effect model
c   cc(i)            multiplicative factor of each nested structure.
c                      (sill-c0) for spherical, exponential,and gaussian
c                      slope for linear model.
c   aa(i)            parameter "a" of each nested structure.
c   irot             index of the rotation matrix for the first nested 
c                    structure (the second nested structure will use
c                    irot+1, the third irot+2, and so on)
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c
c
c OUTPUT VARIABLES:
c
c   cmax             maximum covariance
c   cova             covariance between (x1,y1,z1) and (x2,y2,z2)
c
c
c
c EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
c                      rotmat    computes rotation matrix for distance
c-----------------------------------------------------------------------

      implicit none
      real, intent(in) :: sqdistance
      integer,intent(in) :: ivarg
      integer,intent(in) :: nst(ivarg)
      integer,intent(in) :: MAXNST
      real,intent(in) :: c0(ivarg)
      integer,intent(in) :: it(MAXNST)
      real,intent(in) :: cc(MAXNST),aa(MAXNST),aainv(MAXNST)
      integer,intent(in)::irot,MAXROT
      real*8,intent(in) :: rotmat(MAXROT,3,3)
      real,intent(inout) :: cmax,cova


c      parameter(PI=3.14159265,PMX=999.,EPSLON=1.e-5)
c      integer   nst(*),it(*)
c      real      c0(*),cc(*),aa(*),aainv(*)
c      real*8    rotmat(MAXROT,3,3),hsqd,sqdist

      real PI, EPSLON
      integer PMX
      real*8 hsqd,sqdist
      integer istart,is,ist,ir
      real hr


c      real haainv1,haainv2,haainv3,h1,h2,h3,hsqd1,hsqd2,hsqd3

      PI=3.14159265
      PMX=999
      EPSLON=1.e-5


c
c Calculate the maximum covariance value (used for zero distances and
c for power model covariance):
c


      istart = 1 + (ivarg-1)*MAXNST
      cmax   = c0(ivarg)
      do is=1,nst(ivarg)
            ist = istart + is - 1
            if(it(ist).eq.4) then
                  cmax = cmax + PMX
            else
                  cmax = cmax + cc(ist)
            endif
      end do
c
c Check for "zero" distance, return with cmax if so:
c
c      hsqd = sqdist(x1,y1,z1,x2,y2,z2,irot,MAXROT,rotmat)
      hsqd = sqdistance
      if(real(hsqd).lt.EPSLON) then
            cova = cmax
            return
      endif
c
c Loop over all the structures:
c
      cova = 0.0
      do is=1,nst(ivarg)
            ist = istart + is - 1
c
c Compute the appropriate distance:
c
            if(ist.ne.1) then
                  ir = min((irot+is-1),MAXROT)
                  hsqd = sqdistance
c                  hsqd=sqdist(x1,y1,z1,x2,y2,z2,ir,MAXROT,rotmat)
            end if
c            h = real(dsqrt(hsqd))

c
c Spherical Variogram Model?
c
            if(it(ist).eq.1) then
c                  hr = h/aa(ist)
                  hr = (real(dsqrt(hsqd)))*aainv(ist)
                  if(hr.lt.1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
c
c Exponential Variogram Model?
c
            else if(it(ist).eq.2) then
                  cova = cova + cc(ist)*exp(-3.0*(real(dsqrt(hsqd)))*
     +                                      aainv(ist))
c                  cova = cova + cc(ist)*exp(-3.0*h/aa(ist))

c
c Gaussian Variogram Model?
c
            else if(it(ist).eq.3) then
                  cova = cova + cc(ist)*exp(-3.*(real(hsqd)*aainv(ist)*
     + aainv(ist)))
c                  cova = cova + cc(ist)*exp(-3.*(h/aa(ist))*(h/aa(ist)))

c
c Power Variogram Model?
c
            else if(it(ist).eq.4) then
                  cova = cova + cmax - cc(ist)*((real(dsqrt(hsqd)))**
     +                                            aa(ist))
c                  cova = cova + cmax - cc(ist)*(h**aa(ist))

c
c Hole Effect Model?
c
            else if(it(ist).eq.5) then
c                 d = 10.0 * aa(ist)
c                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
c                  cova = cova + cc(ist)*cos(h/aa(ist)*PI)
                  cova = cova + cc(ist)*cos((real(dsqrt(hsqd)))*
     +                                      aainv(ist)*PI)
            endif
      end do

c
c Finished:
c
      return
      end


      subroutine cova3_opt02(x1,y1,z1,x2,y2,z2,ivarg,nst,MAXNST,c0,it,
     +                 cc,aa,aainv,
     +                 irot,MAXROT,rotmat,cmax,cova)
c-----------------------------------------------------------------------
c
c                    Covariance Between Two Points
c                    *****************************
c
c This subroutine calculated the covariance associated with a variogram
c model specified by a nugget effect and nested varigoram structures.
c The anisotropy definition can be different for each nested structure.
c
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         coordinates of first point
c   x2,y2,z2         coordinates of second point
c   nst(ivarg)       number of nested structures (maximum of 4)
c   ivarg            variogram number (set to 1 unless doing cokriging
c                       or indicator kriging)
c   MAXNST           size of variogram parameter arrays
c   c0(ivarg)        isotropic nugget constant
c   it(i)            type of each nested structure:
c                      1. spherical model of range a;
c                      2. exponential model of parameter a;
c                           i.e. practical range is 3a
c                      3. gaussian model of parameter a;
c                           i.e. practical range is a*sqrt(3)
c                      4. power model of power a (a must be gt. 0  and
c                           lt. 2).  if linear model, a=1,c=slope.
c                      5. hole effect model
c   cc(i)            multiplicative factor of each nested structure.
c                      (sill-c0) for spherical, exponential,and gaussian
c                      slope for linear model.
c   aa(i)            parameter "a" of each nested structure.
c   irot             index of the rotation matrix for the first nested 
c                    structure (the second nested structure will use
c                    irot+1, the third irot+2, and so on)
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c
c
c OUTPUT VARIABLES:
c
c   cmax             maximum covariance
c   cova             covariance between (x1,y1,z1) and (x2,y2,z2)
c
c
c
c EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
c                      rotmat    computes rotation matrix for distance
c-----------------------------------------------------------------------
      parameter(PI=3.14159265,PMX=999.,EPSLON=1.e-5)
      integer   nst(*),it(*)
      real      c0(*),cc(*),aa(*),aainv(*)
c      real*8    rotmat(MAXROT,3,3),hsqd,sqdist
      real*8    rotmat(3,3,MAXROT),hsqd,sqdist

c      real haainv1,haainv2,haainv3,h1,h2,h3,hsqd1,hsqd2,hsqd3



c
c Calculate the maximum covariance value (used for zero distances and
c for power model covariance):
c


      istart = 1 + (ivarg-1)*MAXNST
      cmax   = c0(ivarg)
      do is=1,nst(ivarg)
            ist = istart + is - 1
            if(it(ist).eq.4) then
                  cmax = cmax + PMX
            else
                  cmax = cmax + cc(ist)
            endif
      end do
c
c Check for "zero" distance, return with cmax if so:
c
c      hsqd = sqdist(x1,y1,z1,x2,y2,z2,irot,MAXROT,rotmat)
      hsqd = sqdist_opt01(x1,y1,z1,x2,y2,z2,irot,MAXROT,rotmat)
      if(real(hsqd).lt.EPSLON) then
            cova = cmax
            return
      endif
c
c Loop over all the structures:
c
      cova = 0.0
      do is=1,nst(ivarg)
            ist = istart + is - 1
c
c Compute the appropriate distance:
c
            if(ist.ne.1) then
                  ir = min((irot+is-1),MAXROT)
c                  hsqd=sqdist(x1,y1,z1,x2,y2,z2,ir,MAXROT,rotmat)
                  hsqd=sqdist_opt01(x1,y1,z1,x2,y2,z2,ir,MAXROT,rotmat)
            end if
c            h = real(dsqrt(hsqd))

c
c Spherical Variogram Model?
c
            if(it(ist).eq.1) then
c                  hr = h/aa(ist)
                  hr = (real(dsqrt(hsqd)))*aainv(ist)
                  if(hr.lt.1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
c
c Exponential Variogram Model?
c
            else if(it(ist).eq.2) then
                  cova = cova + cc(ist)*exp(-3.0*(real(dsqrt(hsqd)))*
     +                                      aainv(ist))
c                  cova = cova + cc(ist)*exp(-3.0*h/aa(ist))

c
c Gaussian Variogram Model?
c
            else if(it(ist).eq.3) then
                  cova = cova + cc(ist)*exp(-3.*(real(hsqd)*aainv(ist)*
     + aainv(ist)))
c                  cova = cova + cc(ist)*exp(-3.*(h/aa(ist))*(h/aa(ist)))

c
c Power Variogram Model?
c
            else if(it(ist).eq.4) then
                  cova = cova + cmax - cc(ist)*((real(dsqrt(hsqd)))**
     +                                            aa(ist))
c                  cova = cova + cmax - cc(ist)*(h**aa(ist))

c
c Hole Effect Model?
c
            else if(it(ist).eq.5) then
c                 d = 10.0 * aa(ist)
c                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
c                  cova = cova + cc(ist)*cos(h/aa(ist)*PI)
                  cova = cova + cc(ist)*cos((real(dsqrt(hsqd)))*
     +                                      aainv(ist)*PI)
            endif
      end do

c
c Finished:
c
      return
      end


