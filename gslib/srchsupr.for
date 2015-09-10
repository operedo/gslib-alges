
      subroutine srchsupr(xloc,yloc,zloc,radsqd,irot,MAXROT,
     +                    rotmat,
     +                    nsbtosr,ixsbtosr,iysbtosr,izsbtosr,noct,nd,
     +                    x,y,z,tmp,nisb,nxsup,xmnsup,xsizsup,
     +                    nysup,ymnsup,ysizsup,nzsup,zmnsup,zsizsup,
     +                    nclose,close,infoct)
c-----------------------------------------------------------------------
c
c              Search Within Super Block Search Limits
c              ***************************************
c
c
c This subroutine searches through all the data that have been tagged in
c the super block subroutine.  The close data are passed back in the
c index array "close".  An octant search is allowed.
c
c
c
c INPUT VARIABLES:
c
c   xloc,yloc,zloc   location of point being estimated/simulated
c   radsqd           squared search radius
c   irot             index of the rotation matrix for searching
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c   nsbtosr          Number of super blocks to search
c   ixsbtosr         X offsets for super blocks to search
c   iysbtosr         Y offsets for super blocks to search
c   izsbtosr         Z offsets for super blocks to search
c   noct             If >0 then data will be partitioned into octants
c   nd               Number of data
c   x(nd)            X coordinates of the data
c   y(nd)            Y coordinates of the data
c   z(nd)            Z coordinates of the data
c   tmp(nd)          Temporary storage to keep track of the squared
c                      distance associated with each data
c   nisb()                Array with cumulative number of data in each
c                           super block.
c   nxsup,xmnsup,xsizsup  Definition of the X super block grid
c   nysup,ymnsup,ysizsup  Definition of the X super block grid
c   nzsup,zmnsup,zsizsup  Definition of the X super block grid
c
c
c
c OUTPUT VARIABLES:
c
c   nclose           Number of close data
c   close()          Index of close data
c   infoct           Number of informed octants (only computes if
c                      performing an octant search)
c
c
c
c EXTERNAL REFERENCES:
c
c   sqdist           Computes anisotropic squared distance
c   sortem           Sorts multiple arrays in ascending order
c
c
c
c-----------------------------------------------------------------------
      real    x(*),y(*),z(*),tmp(*),close(*)
      real*8  rotmat(MAXROT,3,3),hsqd,sqdist
c      real*8  rotmat(3,3,MAXROT),hsqd,sqdist
      integer nisb(*),inoct(8)
      integer ixsbtosr(*),iysbtosr(*),izsbtosr(*)
      logical inflag
c
c Determine the super block location of point being estimated:
c
      call getindx(nxsup,xmnsup,xsizsup,xloc,ix,inflag)
      call getindx(nysup,ymnsup,ysizsup,yloc,iy,inflag)
      call getindx(nzsup,zmnsup,zsizsup,zloc,iz,inflag)
c
c Loop over all the possible Super Blocks:
c
      nclose = 0
      do 1 isup=1,nsbtosr
c
c Is this super block within the grid system:
c
            ixsup = ix + ixsbtosr(isup)
            iysup = iy + iysbtosr(isup)
            izsup = iz + izsbtosr(isup)
            if(ixsup.le.0.or.ixsup.gt.nxsup.or.
     +         iysup.le.0.or.iysup.gt.nysup.or.
     +         izsup.le.0.or.izsup.gt.nzsup) go to 1
c
c Figure out how many samples in this super block:
c
            ii = ixsup + (iysup-1)*nxsup + (izsup-1)*nxsup*nysup
            if(ii.eq.1) then
                  nums = nisb(ii)
                  i    = 0
            else
c                  print *,ii
                  nums = nisb(ii) - nisb(ii-1)
                  i    = nisb(ii-1)
            endif
c
c Loop over all the data in this super block:
c
            do 2 ii=1,nums
                  i = i + 1
c
c Check squared distance:
c
                  hsqd = sqdist(xloc,yloc,zloc,x(i),y(i),z(i),irot,
c                  hsqd = sqdist_opt01(xloc,yloc,zloc,x(i),y(i),z(i),
     +                          MAXROT,rotmat)
                  if(real(hsqd).gt.radsqd) go to 2
c
c Accept this sample:
c
                  nclose = nclose + 1
                  close(nclose) = real(i)
                  tmp(nclose)  = real(hsqd)
 2          continue
 1    continue

c      print *,'srchsupr:',nclose

c
c Sort the nearby samples by distance to point being estimated:
c
      call sortem(1,nclose,tmp,1,close,c,d,e,f,g,h)
c      call sortem2(1,nclose,tmp,1,close)
c
c If we aren't doing an octant search then just return:
c
      if(noct.le.0) return
c
c PARTITION THE DATA INTO OCTANTS:
c
      do i=1,8
            inoct(i) = 0
      end do
c
c Now pick up the closest samples in each octant:
c
      nt = 8*noct
      na = 0
      do j=1,nclose
            i  = int(close(j))
            h  = tmp(j)
            dx = x(i) - xloc
            dy = y(i) - yloc
            dz = z(i) - zloc
            if(dz.lt.0.) go to 5
            iq=4
            if(dx.le.0.0 .and. dy.gt.0.0) iq=1
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=2
            if(dx.lt.0.0 .and. dy.le.0.0) iq=3
            go to 6
 5          iq=8
            if(dx.le.0.0 .and. dy.gt.0.0) iq=5
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=6
            if(dx.lt.0.0 .and. dy.le.0.0) iq=7
 6          continue
            inoct(iq) = inoct(iq) + 1
c
c Keep this sample if the maximum has not been exceeded:
c
            if(inoct(iq).le.noct) then
                  na = na + 1
                  close(na) = i
                  tmp(na)   = h
                  if(na.eq.nt) go to 7
            endif
      end do
c
c End of data selection. Compute number of informed octants and return:
c
 7    nclose = na
      infoct = 0
      do i=1,8
            if(inoct(i).gt.0) infoct = infoct + 1
      end do
c
c Finished:
c
      return
      end

      subroutine srchsuprO1(xloc,yloc,zloc,radsqd,irot,MAXROT,rotmat,
     +                    nsbtosr,ixsbtosr,iysbtosr,izsbtosr,noct,nd,
     +                    x,y,z,tmp,nisb,nxsup,xmnsup,xsizsup,
     +                    nysup,ymnsup,ysizsup,nzsup,zmnsup,zsizsup,
     +                    nclose,close,infoct)
c-----------------------------------------------------------------------
c
c              Search Within Super Block Search Limits
c              ***************************************
c
c
c This subroutine searches through all the data that have been tagged in
c the super block subroutine.  The close data are passed back in the
c index array "close".  An octant search is allowed.
c
c
c
c INPUT VARIABLES:
c
c   xloc,yloc,zloc   location of point being estimated/simulated
c   radsqd           squared search radius
c   irot             index of the rotation matrix for searching
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c   nsbtosr          Number of super blocks to search
c   ixsbtosr         X offsets for super blocks to search
c   iysbtosr         Y offsets for super blocks to search
c   izsbtosr         Z offsets for super blocks to search
c   noct             If >0 then data will be partitioned into octants
c   nd               Number of data
c   x(nd)            X coordinates of the data
c   y(nd)            Y coordinates of the data
c   z(nd)            Z coordinates of the data
c   tmp(nd)          Temporary storage to keep track of the squared
c                      distance associated with each data
c   nisb()                Array with cumulative number of data in each
c                           super block.
c   nxsup,xmnsup,xsizsup  Definition of the X super block grid
c   nysup,ymnsup,ysizsup  Definition of the X super block grid
c   nzsup,zmnsup,zsizsup  Definition of the X super block grid
c
c
c
c OUTPUT VARIABLES:
c
c   nclose           Number of close data
c   close()          Index of close data
c   infoct           Number of informed octants (only computes if
c                      performing an octant search)
c
c
c
c EXTERNAL REFERENCES:
c
c   sqdist           Computes anisotropic squared distance
c   sortem           Sorts multiple arrays in ascending order
c
c
c
c-----------------------------------------------------------------------

      implicit none

      real xloc,yloc,zloc,radsqd
      integer irot, MAXROT
      real*8  rotmat(MAXROT,3,3)
      integer nsbtosr
      integer ixsbtosr(*),iysbtosr(*),izsbtosr(*)
      integer noct,nd
      real    x(*),y(*),z(*),tmp(*)
      integer nisb(*)
      integer nxsup
      real xmnsup,xsizsup
      integer nysup
      real ymnsup,ysizsup
      integer nzsup
      real zmnsup,zsizsup
      integer nclose
      real    close(*)
      integer infoct
      

      integer inoct(8)
      logical inflag
      real*8  hsqd,sqdist
      real c,d,e,f,g,h
      real dx,dy,dz
      integer i,j,ii,isup,nums,ixsup,iysup,izsup,ix,iy,iz
      integer ixsup1,iysup1,izsup1
      integer ixsup2,iysup2,izsup2
      integer ixsup3,iysup3,izsup3
      integer na,iq,nt,liminf
c
c Determine the super block location of point being estimated:
c
      call getindx(nxsup,xmnsup,xsizsup,xloc,ix,inflag)
      call getindx(nysup,ymnsup,ysizsup,yloc,iy,inflag)
      call getindx(nzsup,zmnsup,zsizsup,zloc,iz,inflag)
c
c Loop over all the possible Super Blocks:
c
      nclose = 0


      do 1 isup=1,nsbtosr

c
c Is this super block within the grid system:
c
            ixsup = ix + ixsbtosr(isup)
c            ixsup1 = ix + ixsbtosr(isup+1)
c            ixsup2 = ix + ixsbtosr(isup+2)
c            ixsup3 = ix + ixsbtosr(isup+3)

            iysup = iy + iysbtosr(isup)
c            iysup1 = iy + iysbtosr(isup+1)
c            iysup2 = iy + iysbtosr(isup+2)
c            iysup3 = iy + iysbtosr(isup+3)

            izsup = iz + izsbtosr(isup)
c            izsup1 = iz + izsbtosr(isup+1)
c            izsup2 = iz + izsbtosr(isup+2)
c            izsup3 = iz + izsbtosr(isup+3)

            if(ixsup.le.0.or.ixsup.gt.nxsup.or.
     +         iysup.le.0.or.iysup.gt.nysup.or.
     +         izsup.le.0.or.izsup.gt.nzsup) go to 1

c
c Figure out how many samples in this super block:
c
            ii = ixsup + (iysup-1)*nxsup + (izsup-1)*nxsup*nysup

c            i = transfer(ii.ne.1,i) * nisb(ii-1)
c            nums = nisb(ii)- i

            if(ii.eq.1) then
                  nums = nisb(ii)
                  i    = 0
            else
                  nums = nisb(ii) - nisb(ii-1)
                  i    = nisb(ii-1)
            endif

c
c Loop over all the data in this super block:
c
            do 2 ii=1,nums
                  i = i + 1
c
c Check squared distance:
c
                  hsqd = sqdist(xloc,yloc,zloc,x(i),y(i),z(i),irot,
     +                          MAXROT,rotmat)
                  if(real(hsqd).gt.radsqd) go to 2
c
c Accept this sample:
c
                  nclose = nclose + 1
                  close(nclose) = real(i)
                  tmp(nclose)  = real(hsqd)
 2          continue

c            end if

cc
cc Is this super block within the grid system:
cc
c 11         if(ixsup1.le.0.or.ixsup1.gt.nxsup.or.
c     +         iysup1.le.0.or.iysup1.gt.nysup.or.
c     +         izsup1.le.0.or.izsup1.gt.nzsup) go to 1
c
cc
cc Figure out how many samples in this super block:
cc
c            ii = ixsup1 + (iysup1-1)*nxsup + (izsup1-1)*nxsup*nysup
c
c            i = transfer(ii.ne.1,i) * nisb(ii-1)
c            nums = nisb(ii)- i
c
cc            if(ii.eq.1) then
cc                  nums = nisb(ii)
cc                  i    = 0
cc            else
cc                  nums = nisb(ii) - nisb(ii-1)
cc                  i    = nisb(ii-1)
cc            endif
c
cc
cc Loop over all the data in this super block:
cc
c            do 22 ii=1,nums
c                  i = i + 1
cc
cc Check squared distance:
cc
c                  hsqd = sqdist(xloc,yloc,zloc,x(i),y(i),z(i),irot,
c     +                          MAXROT,rotmat)
c                  if(real(hsqd).gt.radsqd) go to 22
cc
cc Accept this sample:
cc
c                  nclose = nclose + 1
c                  close(nclose) = real(i)
c                  tmp(nclose)  = real(hsqd)
c 22         continue
c
cc
cc Is this super block within the grid system:
cc
c 12         if(ixsup2.le.0.or.ixsup2.gt.nxsup.or.
c     +         iysup2.le.0.or.iysup2.gt.nysup.or.
c     +         izsup2.le.0.or.izsup2.gt.nzsup) go to 13
c
cc
cc Figure out how many samples in this super block:
cc
c            ii = ixsup2 + (iysup2-1)*nxsup + (izsup2-1)*nxsup*nysup
c
c            i = transfer(ii.ne.1,i) * nisb(ii-1)
c            nums = nisb(ii)- i
c
cc            if(ii.eq.1) then
cc                  nums = nisb(ii)
cc                  i    = 0
cc            else
cc                  nums = nisb(ii) - nisb(ii-1)
cc                  i    = nisb(ii-1)
cc            endif
c
cc
cc Loop over all the data in this super block:
cc
c            do 23 ii=1,nums
c                  i = i + 1
cc
cc Check squared distance:
cc
c                  hsqd = sqdist(xloc,yloc,zloc,x(i),y(i),z(i),irot,
c     +                          MAXROT,rotmat)
c                  if(real(hsqd).gt.radsqd) go to 23
cc
cc Accept this sample:
cc
c                  nclose = nclose + 1
c                  close(nclose) = real(i)
c                  tmp(nclose)  = real(hsqd)
c 23         continue
c
cc
cc Is this super block within the grid system:
cc
c 13         if(ixsup3.le.0.or.ixsup3.gt.nxsup.or.
c     +         iysup3.le.0.or.iysup3.gt.nysup.or.
c     +         izsup3.le.0.or.izsup3.gt.nzsup) go to 1
c
cc
cc Figure out how many samples in this super block:
cc
c            ii = ixsup3 + (iysup3-1)*nxsup + (izsup3-1)*nxsup*nysup
c
c            i = transfer(ii.ne.1,i) * nisb(ii-1)
c            nums = nisb(ii)- i
c
cc            if(ii.eq.1) then
cc                  nums = nisb(ii)
cc                  i    = 0
cc            else
cc                  nums = nisb(ii) - nisb(ii-1)
cc                  i    = nisb(ii-1)
cc            endif
c
cc
cc Loop over all the data in this super block:
cc
c            do 24 ii=1,nums
c                  i = i + 1
cc
cc Check squared distance:
cc
c                  hsqd = sqdist(xloc,yloc,zloc,x(i),y(i),z(i),irot,
c     +                          MAXROT,rotmat)
c                  if(real(hsqd).gt.radsqd) go to 24
cc
cc Accept this sample:
cc
c                  nclose = nclose + 1
c                  close(nclose) = real(i)
c                  tmp(nclose)  = real(hsqd)
c 24         continue

 1    continue

c      liminf=isup
c
c      do 111 isup=liminf,nsbtosr
cc
cc Is this super block within the grid system:
cc
c            ixsup = ix + ixsbtosr(isup)
c            iysup = iy + iysbtosr(isup)
c            izsup = iz + izsbtosr(isup)
c
c            if(ixsup.le.0.or.ixsup.gt.nxsup.or.
c     +         iysup.le.0.or.iysup.gt.nysup.or.
c     +         izsup.le.0.or.izsup.gt.nzsup) go to 111
cc
cc Figure out how many samples in this super block:
cc
c            ii = ixsup + (iysup-1)*nxsup + (izsup-1)*nxsup*nysup
c
c            i = transfer(ii.ne.1,i) * nisb(ii-1)
c            nums = nisb(ii)- i
c
cc            if(ii.eq.1) then
cc                  nums = nisb(ii)
cc                  i    = 0
cc            else
cc                  nums = nisb(ii) - nisb(ii-1)
cc                  i    = nisb(ii-1)
cc            endif
c
cc
cc Loop over all the data in this super block:
cc
c            do 222 ii=1,nums
c                  i = i + 1
cc
cc Check squared distance:
cc
c                  hsqd = sqdist(xloc,yloc,zloc,x(i),y(i),z(i),irot,
c     +                          MAXROT,rotmat)
c                  if(real(hsqd).gt.radsqd) go to 222
cc
cc Accept this sample:
cc
c                  nclose = nclose + 1
c                  close(nclose) = real(i)
c                  tmp(nclose)  = real(hsqd)
c 222        continue
c
cc            end if
c
c 111  continue






c      print *,'srchsupr:',nclose

c
c Sort the nearby samples by distance to point being estimated:
c
c      call sortem(1,nclose,tmp,1,close,c,d,e,f,g,h)
      call sortem2(1,nclose,tmp,1,close)
c
c If we aren't doing an octant search then just return:
c
      if(noct.le.0) return
c
c PARTITION THE DATA INTO OCTANTS:
c
      do i=1,8
            inoct(i) = 0
      end do
c
c Now pick up the closest samples in each octant:
c
      nt = 8*noct
      na = 0
      do j=1,nclose
            i  = int(close(j))
            h  = tmp(j)
            dx = x(i) - xloc
            dy = y(i) - yloc
            dz = z(i) - zloc
            if(dz.lt.0.) go to 5
            iq=4
            if(dx.le.0.0 .and. dy.gt.0.0) iq=1
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=2
            if(dx.lt.0.0 .and. dy.le.0.0) iq=3
            go to 6
 5          iq=8
            if(dx.le.0.0 .and. dy.gt.0.0) iq=5
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=6
            if(dx.lt.0.0 .and. dy.le.0.0) iq=7
 6          continue
            inoct(iq) = inoct(iq) + 1
c
c Keep this sample if the maximum has not been exceeded:
c
            if(inoct(iq).le.noct) then
                  na = na + 1
                  close(na) = i
                  tmp(na)   = h
                  if(na.eq.nt) go to 7
            endif
      end do
c
c End of data selection. Compute number of informed octants and return:
c
 7    nclose = na
      infoct = 0
      do i=1,8
            if(inoct(i).gt.0) infoct = infoct + 1
      end do
c
c Finished:
c
      return
      end

      subroutine srchsuprO2(xloc,yloc,zloc,radsqd,irot,MAXROT,rotmat,
     +                    nsbtosr,ixsbtosr,iysbtosr,izsbtosr,noct,nd,
     +                    x,y,z,tmp,nisb,nxsup,xmnsup,xsizsup,
     +                    nysup,ymnsup,ysizsup,nzsup,zmnsup,zsizsup,
     +                    nclose,close,infoct)
c-----------------------------------------------------------------------
c
c              Search Within Super Block Search Limits
c              ***************************************
c
c
c This subroutine searches through all the data that have been tagged in
c the super block subroutine.  The close data are passed back in the
c index array "close".  An octant search is allowed.
c
c
c
c INPUT VARIABLES:
c
c   xloc,yloc,zloc   location of point being estimated/simulated
c   radsqd           squared search radius
c   irot             index of the rotation matrix for searching
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c   nsbtosr          Number of super blocks to search
c   ixsbtosr         X offsets for super blocks to search
c   iysbtosr         Y offsets for super blocks to search
c   izsbtosr         Z offsets for super blocks to search
c   noct             If >0 then data will be partitioned into octants
c   nd               Number of data
c   x(nd)            X coordinates of the data
c   y(nd)            Y coordinates of the data
c   z(nd)            Z coordinates of the data
c   tmp(nd)          Temporary storage to keep track of the squared
c                      distance associated with each data
c   nisb()                Array with cumulative number of data in each
c                           super block.
c   nxsup,xmnsup,xsizsup  Definition of the X super block grid
c   nysup,ymnsup,ysizsup  Definition of the X super block grid
c   nzsup,zmnsup,zsizsup  Definition of the X super block grid
c
c
c
c OUTPUT VARIABLES:
c
c   nclose           Number of close data
c   close()          Index of close data
c   infoct           Number of informed octants (only computes if
c                      performing an octant search)
c
c
c
c EXTERNAL REFERENCES:
c
c   sqdist           Computes anisotropic squared distance
c   sortem           Sorts multiple arrays in ascending order
c
c
c
c-----------------------------------------------------------------------

      implicit none

      real xloc,yloc,zloc,radsqd
      integer irot, MAXROT
      real*8  rotmat(MAXROT,3,3)
      integer nsbtosr
      integer ixsbtosr(*),iysbtosr(*),izsbtosr(*)
      integer noct,nd
      real    x(*),y(*),z(*),tmp(*)
      integer nisb(*)
      integer nxsup
      real xmnsup,xsizsup
      integer nysup
      real ymnsup,ysizsup
      integer nzsup
      real zmnsup,zsizsup
      integer nclose
      real    close(*)
      integer infoct
      

      integer inoct(8)
      logical inflag
      real*8  hsqd,sqdist
      real c,d,e,f,g,h
      real dx,dy,dz
      integer i,j,ii,isup,nums,ixsup,iysup,izsup,ix,iy,iz
      integer na,iq,nt
c
c Determine the super block location of point being estimated:
c
      call getindx(nxsup,xmnsup,xsizsup,xloc,ix,inflag)
      call getindx(nysup,ymnsup,ysizsup,yloc,iy,inflag)
      call getindx(nzsup,zmnsup,zsizsup,zloc,iz,inflag)
c
c Loop over all the possible Super Blocks:
c
      nclose = 0
      do 1 isup=1,nsbtosr
c
c Is this super block within the grid system:
c
            ixsup = ix + ixsbtosr(isup)
            iysup = iy + iysbtosr(isup)
            izsup = iz + izsbtosr(isup)
c            if(ixsup.le.0.or.ixsup.gt.nxsup.or.
c     +         iysup.le.0.or.iysup.gt.nysup.or.
c     +         izsup.le.0.or.izsup.gt.nzsup) go to 1

            if(ixsup.gt.0.and.ixsup.le.nxsup.and.
     +         iysup.gt.0.and.iysup.le.nysup.and.
     +         izsup.gt.0.and.izsup.le.nzsup) then

c
c Figure out how many samples in this super block:
c
            ii = ixsup + (iysup-1)*nxsup + (izsup-1)*nxsup*nysup

            i = transfer(ii.ne.1,i) * nisb(ii-1) 
            nums = nisb(ii) - i 
c            if(ii.eq.1) then
c                  nums = nisb(ii)
c                  i    = 0
c            else
c                  nums = nisb(ii) - nisb(ii-1)
c                  i    = nisb(ii-1)
c            endif

c
c Loop over all the data in this super block:
c
            do 2 ii=1,nums
                  i = i + 1
c
c Check squared distance:
c
                  hsqd = sqdist(xloc,yloc,zloc,x(i),y(i),z(i),irot,
     +                          MAXROT,rotmat)
                  if(real(hsqd).gt.radsqd) go to 2
c
c Accept this sample:
c
                  nclose = nclose + 1
                  close(nclose) = real(i)
                  tmp(nclose)  = real(hsqd)
 2          continue

            end if

 1    continue

c      print *,'srchsupr:',nclose

c
c Sort the nearby samples by distance to point being estimated:
c
c      call sortem(1,nclose,tmp,1,close,c,d,e,f,g,h)
      call sortem2(1,nclose,tmp,1,close)
c
c If we aren't doing an octant search then just return:
c
      if(noct.le.0) return
c
c PARTITION THE DATA INTO OCTANTS:
c
      do i=1,8
            inoct(i) = 0
      end do
c
c Now pick up the closest samples in each octant:
c
      nt = 8*noct
      na = 0
      do j=1,nclose
            i  = int(close(j))
            h  = tmp(j)
            dx = x(i) - xloc
            dy = y(i) - yloc
            dz = z(i) - zloc
            if(dz.lt.0.) go to 5
            iq=4
            if(dx.le.0.0 .and. dy.gt.0.0) iq=1
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=2
            if(dx.lt.0.0 .and. dy.le.0.0) iq=3
            go to 6
 5          iq=8
            if(dx.le.0.0 .and. dy.gt.0.0) iq=5
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=6
            if(dx.lt.0.0 .and. dy.le.0.0) iq=7
 6          continue
            inoct(iq) = inoct(iq) + 1
c
c Keep this sample if the maximum has not been exceeded:
c
            if(inoct(iq).le.noct) then
                  na = na + 1
                  close(na) = i
                  tmp(na)   = h
                  if(na.eq.nt) go to 7
            endif
      end do
c
c End of data selection. Compute number of informed octants and return:
c
 7    nclose = na
      infoct = 0
      do i=1,8
            if(inoct(i).gt.0) infoct = infoct + 1
      end do
c
c Finished:
c
      return
      end






c      subroutine srchsuprO2(xloc,yloc,zloc,radsqd,irot,MAXROT,rotmat,
c     +                    nsbtosr,ixyzsbtosr,noct,nd,
c     +                    x,y,z,tmp,nisb,nxsup,xmnsup,xsizsup,
c     +                    nysup,ymnsup,ysizsup,nzsup,zmnsup,zsizsup,
c     +                    nclose,close,infoct)
cc-----------------------------------------------------------------------
cc
cc              Search Within Super Block Search Limits
cc              ***************************************
cc
cc
cc This subroutine searches through all the data that have been tagged in
cc the super block subroutine.  The close data are passed back in the
cc index array "close".  An octant search is allowed.
cc
cc
cc
cc INPUT VARIABLES:
cc
cc   xloc,yloc,zloc   location of point being estimated/simulated
cc   radsqd           squared search radius
cc   irot             index of the rotation matrix for searching
cc   MAXROT           size of rotation matrix arrays
cc   rotmat           rotation matrices
cc   nsbtosr          Number of super blocks to search
cc   ixsbtosr         X offsets for super blocks to search
cc   iysbtosr         Y offsets for super blocks to search
cc   izsbtosr         Z offsets for super blocks to search
cc   noct             If >0 then data will be partitioned into octants
cc   nd               Number of data
cc   x(nd)            X coordinates of the data
cc   y(nd)            Y coordinates of the data
cc   z(nd)            Z coordinates of the data
cc   tmp(nd)          Temporary storage to keep track of the squared
cc                      distance associated with each data
cc   nisb()                Array with cumulative number of data in each
cc                           super block.
cc   nxsup,xmnsup,xsizsup  Definition of the X super block grid
cc   nysup,ymnsup,ysizsup  Definition of the X super block grid
cc   nzsup,zmnsup,zsizsup  Definition of the X super block grid
cc
cc
cc
cc OUTPUT VARIABLES:
cc
cc   nclose           Number of close data
cc   close()          Index of close data
cc   infoct           Number of informed octants (only computes if
cc                      performing an octant search)
cc
cc
cc
cc EXTERNAL REFERENCES:
cc
cc   sqdist           Computes anisotropic squared distance
cc   sortem           Sorts multiple arrays in ascending order
cc
cc
cc
cc-----------------------------------------------------------------------
c
c      implicit none
c
c      real xloc,yloc,zloc,radsqd
c      integer irot, MAXROT
c      real*8  rotmat(MAXROT,3,3)
c      integer nsbtosr
cc      integer ixsbtosr(*),iysbtosr(*),izsbtosr(*)
c      integer ixyzsbtosr(*)
c      integer noct,nd
c      real    x(*),y(*),z(*),tmp(*)
c      integer nisb(*)
c      integer nxsup
c      real xmnsup,xsizsup
c      integer nysup
c      real ymnsup,ysizsup
c      integer nzsup
c      real zmnsup,zsizsup
c      integer nclose
c      real    close(*)
c      integer infoct
c      
c
c      integer inoct(8)
c      logical inflag
c      real*8  hsqd,sqdist
c      real c,d,e,f,g,h
c      real dx,dy,dz
c      integer i,j,ii,isup,nums,ixsup,iysup,izsup,ix,iy,iz
c      integer na,iq,nt
cc
cc Determine the super block location of point being estimated:
cc
c      call getindx(nxsup,xmnsup,xsizsup,xloc,ix,inflag)
c      call getindx(nysup,ymnsup,ysizsup,yloc,iy,inflag)
c      call getindx(nzsup,zmnsup,zsizsup,zloc,iz,inflag)
cc
cc Loop over all the possible Super Blocks:
cc
c      nclose = 0
c      do 1 isup=1,nsbtosr
cc
cc Is this super block within the grid system:
cc
cc            ixsup = ix + ixsbtosr(isup)
cc            iysup = iy + iysbtosr(isup)
cc            izsup = iz + izsbtosr(isup)
c            ixsup = ix + ixyzsbtosr(3*(isup-1)+1)
c            iysup = iy + ixyzsbtosr(3*(isup-1)+2)
c            izsup = iz + ixyzsbtosr(3*(isup-1)+3)
c
c            if(ixsup.le.0.or.ixsup.gt.nxsup.or.
c     +         iysup.le.0.or.iysup.gt.nysup.or.
c     +         izsup.le.0.or.izsup.gt.nzsup) go to 1
cc
cc Figure out how many samples in this super block:
cc
c            ii = ixsup + (iysup-1)*nxsup + (izsup-1)*nxsup*nysup
c            if(ii.eq.1) then
c                  nums = nisb(ii)
c                  i    = 0
c            else
c                  nums = nisb(ii) - nisb(ii-1)
c                  i    = nisb(ii-1)
c            endif
c
cc
cc Loop over all the data in this super block:
cc
c            do 2 ii=1,nums
c                  i = i + 1
cc
cc Check squared distance:
cc
c                  hsqd = sqdist(xloc,yloc,zloc,x(i),y(i),z(i),irot,
c     +                          MAXROT,rotmat)
c                  if(real(hsqd).gt.radsqd) go to 2
cc
cc Accept this sample:
cc
c                  nclose = nclose + 1
c                  close(nclose) = real(i)
c                  tmp(nclose)  = real(hsqd)
c 2          continue
c 1    continue
c
cc      print *,'srchsupr:',nclose
c
cc
cc Sort the nearby samples by distance to point being estimated:
cc
c      call sortem(1,nclose,tmp,1,close,c,d,e,f,g,h)
cc
cc If we aren't doing an octant search then just return:
cc
c      if(noct.le.0) return
cc
cc PARTITION THE DATA INTO OCTANTS:
cc
c      do i=1,8
c            inoct(i) = 0
c      end do
cc
cc Now pick up the closest samples in each octant:
cc
c      nt = 8*noct
c      na = 0
c      do j=1,nclose
c            i  = int(close(j))
c            h  = tmp(j)
c            dx = x(i) - xloc
c            dy = y(i) - yloc
c            dz = z(i) - zloc
c            if(dz.lt.0.) go to 5
c            iq=4
c            if(dx.le.0.0 .and. dy.gt.0.0) iq=1
c            if(dx.gt.0.0 .and. dy.ge.0.0) iq=2
c            if(dx.lt.0.0 .and. dy.le.0.0) iq=3
c            go to 6
c 5          iq=8
c            if(dx.le.0.0 .and. dy.gt.0.0) iq=5
c            if(dx.gt.0.0 .and. dy.ge.0.0) iq=6
c            if(dx.lt.0.0 .and. dy.le.0.0) iq=7
c 6          continue
c            inoct(iq) = inoct(iq) + 1
cc
cc Keep this sample if the maximum has not been exceeded:
cc
c            if(inoct(iq).le.noct) then
c                  na = na + 1
c                  close(na) = i
c                  tmp(na)   = h
c                  if(na.eq.nt) go to 7
c            endif
c      end do
cc
cc End of data selection. Compute number of informed octants and return:
cc
c 7    nclose = na
c      infoct = 0
c      do i=1,8
c            if(inoct(i).gt.0) infoct = infoct + 1
c      end do
cc
cc Finished:
cc
c      return
c      end


