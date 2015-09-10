      real*8 function sqdist(x1,y1,z1,x2,y2,z2,ind,MAXROT,rotmat)
c-----------------------------------------------------------------------
c
c    Squared Anisotropic Distance Calculation Given Matrix Indicator
c    ***************************************************************
c
c This routine calculates the anisotropic distance between two points
c  given the coordinates of each point and a definition of the
c  anisotropy.
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         Coordinates of first point
c   x2,y2,z2         Coordinates of second point
c   ind              The rotation matrix to use
c   MAXROT           The maximum number of rotation matrices dimensioned
c   rotmat           The rotation matrices
c
c
c
c OUTPUT VARIABLES:
c
c   sqdist           The squared distance accounting for the anisotropy
c                      and the rotation of coordinates (if any).
c
c
c NO EXTERNAL REFERENCES
c
c
c-----------------------------------------------------------------------
      real*8 rotmat(MAXROT,3,3),cont,dx,dy,dz

c      real*8 cont1,cont2,cont3
c      real*8 rmi11,rmi12,rmi13,
c     +       rmi21,rmi22,rmi23,
c     +       rmi31,rmi32,rmi33

c
c Compute component distance vectors and the squared distance:
c
      dx = dble(x1 - x2)
      dy = dble(y1 - y2)
      dz = dble(z1 - z2)
      sqdist = 0.0
      do i=1,3
            cont   = rotmat(ind,i,1) * dx
     +             + rotmat(ind,i,2) * dy
     +             + rotmat(ind,i,3) * dz
            sqdist = sqdist + cont * cont
      end do

c      do i=1,3
c            cont   = rotmat(ind,1,1) * dx
c     +             + rotmat(ind,1,2) * dy
c     +             + rotmat(ind,1,3) * dz
c            sqdist = sqdist + cont * cont
c            cont   = rotmat(ind,2,1) * dx
c     +             + rotmat(ind,2,2) * dy
c     +             + rotmat(ind,2,3) * dz
c            sqdist = sqdist + cont * cont
c            cont   = rotmat(ind,3,1) * dx
c     +             + rotmat(ind,3,2) * dy
c     +             + rotmat(ind,3,3) * dz
c            sqdist = sqdist + cont * cont
c      end do

c      rmi11=rotmat(ind,1,1)
c      rmi21=rotmat(ind,2,1)
c      rmi31=rotmat(ind,3,1)
c      rmi12=rotmat(ind,1,2)
c      rmi22=rotmat(ind,2,2)
c      rmi32=rotmat(ind,3,2)
c      rmi13=rotmat(ind,1,3)
c      rmi23=rotmat(ind,2,3)
c      rmi33=rotmat(ind,3,3)
c
c      cont1   = rmi11 * dx
c     +       + rmi12 * dy
c     +       + rmi13 * dz
c
c      cont2   = rmi21 * dx
c     +       + rmi22 * dy
c     +       + rmi23 * dz
c
c      cont3   = rmi31 * dx
c     +       + rmi32 * dy
c     +       + rmi33 * dz
c
c      sqdist = cont1 * cont1 + cont2 * cont2 + cont3 * cont3 

      return
      end



      real*8 function sqdist_opt01(x1,y1,z1,x2,y2,z2,ind,MAXROT,rotmat)
c-----------------------------------------------------------------------
c
c    Squared Anisotropic Distance Calculation Given Matrix Indicator
c    ***************************************************************
c
c This routine calculates the anisotropic distance between two points
c  given the coordinates of each point and a definition of the
c  anisotropy.
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         Coordinates of first point
c   x2,y2,z2         Coordinates of second point
c   ind              The rotation matrix to use
c   MAXROT           The maximum number of rotation matrices dimensioned
c   rotmat           The rotation matrices
c
c
c
c OUTPUT VARIABLES:
c
c   sqdist           The squared distance accounting for the anisotropy
c                      and the rotation of coordinates (if any).
c
c
c NO EXTERNAL REFERENCES
c
c
c-----------------------------------------------------------------------
      real*8 rotmat(3,3,MAXROT),cont,dx,dy,dz

c      real*8 cont1,cont2,cont3
c      real*8 rmi11,rmi12,rmi13,
c     +       rmi21,rmi22,rmi23,
c     +       rmi31,rmi32,rmi33

c
c Compute component distance vectors and the squared distance:
c
      dx = dble(x1 - x2)
      dy = dble(y1 - y2)
      dz = dble(z1 - z2)
      sqdist = 0.0
      do i=1,3
            cont   = rotmat(1,i,ind) * dx
     +             + rotmat(2,i,ind) * dy
     +             + rotmat(3,i,ind) * dz
            sqdist = sqdist + cont * cont
      end do

c      do i=1,3
c            cont   = rotmat(ind,1,1) * dx
c     +             + rotmat(ind,1,2) * dy
c     +             + rotmat(ind,1,3) * dz
c            sqdist = sqdist + cont * cont
c            cont   = rotmat(ind,2,1) * dx
c     +             + rotmat(ind,2,2) * dy
c     +             + rotmat(ind,2,3) * dz
c            sqdist = sqdist + cont * cont
c            cont   = rotmat(ind,3,1) * dx
c     +             + rotmat(ind,3,2) * dy
c     +             + rotmat(ind,3,3) * dz
c            sqdist = sqdist + cont * cont
c      end do

c      rmi11=rotmat(ind,1,1)
c      rmi21=rotmat(ind,2,1)
c      rmi31=rotmat(ind,3,1)
c      rmi12=rotmat(ind,1,2)
c      rmi22=rotmat(ind,2,2)
c      rmi32=rotmat(ind,3,2)
c      rmi13=rotmat(ind,1,3)
c      rmi23=rotmat(ind,2,3)
c      rmi33=rotmat(ind,3,3)
c
c      cont1   = rmi11 * dx
c     +       + rmi12 * dy
c     +       + rmi13 * dz
c
c      cont2   = rmi21 * dx
c     +       + rmi22 * dy
c     +       + rmi23 * dz
c
c      cont3   = rmi31 * dx
c     +       + rmi32 * dy
c     +       + rmi33 * dz
c
c      sqdist = cont1 * cont1 + cont2 * cont2 + cont3 * cont3 

      return
      end
