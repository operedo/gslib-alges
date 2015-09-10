C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                                                      %
C Copyright (C) 2003, Statios Software and Services Incorporated.  All %
C rights reserved.                                                     %
C                                                                      %
C This program has been modified from the one distributed in 1996 (see %
C below).  This version is also distributed in the hope that it will   %
C be useful, but WITHOUT ANY WARRANTY. Compiled programs based on this %
C code may be redistributed without restriction; however, this code is %
C for one developer only. Each developer or user of this source code   %
C must purchase a separate copy from Statios.                          %
C                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                                                      %
C Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %
C Junior University.  All rights reserved.                             %
C                                                                      %
C The programs in GSLIB are distributed in the hope that they will be  %
C useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
C responsibility to anyone for the consequences of using them or for   %
C whether they serve any particular purpose or work at all, unless he  %
C says so in writing.  Everyone is granted permission to copy, modify  %
C and redistribute the programs in GSLIB, but only under the condition %
C that this notice and the above copyright notice remain intact.       %
C                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      program main
c-----------------------------------------------------------------------
c
c               Variogram of Irregularly Spaced 3-D Data
c               ****************************************
c
c This is a template driver program for GSLIB's "gamv" subroutine. The
c input data must be entered with coordinates in a GEOEAS format file.
c The User's Guide can be referred to for more details.
c
c The program is executed with no command line arguments.  The user
c will be prompted for the name of a parameter file.  The parameter
c file is described in the documentation (see the example gamv.par)
c
c
c
c The output file will contain each directional variogram ordered by
c direction and then variogram (the directions cycle fastest then the
c variogram number).  For each variogram there will be a one line
c description and then "nlag" lines with the following:
c
c        a) lag number (increasing from 1 to nlag)
c        b) separation distance
c        c) the "variogram" value
c        d) the number of pairs for the lag
c        e) the mean of the data contributing to the tail
c        f) the mean of the data contributing to the head
c
c
c
c AUTHOR: Clayton V. Deutsch                             DATE: 1989-1999
c OpenMP/MPI Modifications: Oscar F. Peredo              DATE: 2014-2015
c-----------------------------------------------------------------------

#ifdef _OPENMP
      use omp_lib
#endif
c#ifdef TRACE
c      use extrae_module
c#endif
#ifdef USE_MPI
      include 'mpif.h'
#endif

c geostat module      
      real,allocatable         :: x(:),y(:),z(:),vr(:,:),azm(:),atol(:),
     +                            bandwh(:),dip(:),dtol(:),bandwd(:)
      real*8,allocatable       :: sills(:),dis(:),gam(:),hm(:),
     +                            tm(:),hv(:),tv(:),np(:)
      integer,allocatable      :: ivtail(:),ivhead(:),ivtype(:)
      character*12,allocatable :: names(:)

      real      EPSLON,VERSION
      real      xlag,xltol,tmin,tmax
      integer   nd,nlag,ndir,nvarg,isill,test
      character outfl*512

c readparam variables
      parameter(MV=500)

      real      var(MV),cut(MV)
      real*8    avg(MV),ssq(MV)
      integer   ivar(MV),num(MV),ivc(MV),indflag(MV)
      character datafl*512,str*512
      logical   testfl,testdat
      real,allocatable :: vrmin(:),vrmax(:)
      data      lin/1/,ncut/0/

c check variables
      character title*80

c gamv variables
      parameter(PI=3.14159265)
      real      uvxazm(100),uvyazm(100),uvzdec(100),
     +          uvhdec(100),csatol(100),csdtol(100)
      logical   omni
      integer   threadId,numThreads
      real*8,allocatable ::    reducedVariables(:,:,:)
      integer extractValue
      integer i,j,id,ii,il,it,iv,jj
      real dx,dy,dz,dxs,dys,dzs,hs,h
      integer lagbeg,lagend,ilag
      real band,dcazm,dcdec,dxy,gamma,vrh,vrhpr,vrt,vrtpr
      real xi,yi,zi 
      integer liminf,limsup,iinf,isup,k
      real xlaginv
#ifdef USE_MPI
      real*8,allocatable :: nptmp(:),distmp(:),gamtmp(:),hmtmp(:),
     + tmtmp(:),hvtmp(:),tvtmp(:)
      integer ierr, myrank,mpisize,blocksize,nblocks
#endif


c writeout variables
      character titlewr*132
      data      lout/1/

c#ifdef TRACE
c      call extrae_init()
c#endif

#ifdef USE_MPI
      call MPI_INIT(ierr) 
      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
      print *,'Process ',myrank
#else
      myrank=1
#endif

#ifdef _OPENMP
c$omp parallel
      numThreads = OMP_get_num_threads() 
c$omp end parallel
#else
      numThreads = 1
#endif

      EPSLON  = 1.0e-20
      VERSION = 3.000

ccccccccccccccc readparms cccccccccccccccccccccc

c
c Note VERSION number:
c
      write(*,9999) VERSION
 9999 format(/' GAMV Version: ',f5.3/)
c
c Get the name of the parameter file - try the default name if no input:
c
      do i=1,512
            str(i:i) = ' '
      end do
      call getarg(1,str)
      if(str(1:1).eq.' ')then
            write(*,*) 'Which parameter file do you want to use?'
            read (*,'(a)') str
      end if
      if(str(1:1).eq.' ') str(1:20) = 'gamv.par            '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'gamv.par            ') then
                  write(*,*) '        creating a blank parameter file'
                  call makepar
                  write(*,*)
            end if
            stop
      endif
      open(lin,file=str,status='OLD')
c
c Find Start of Parameters:
c
 1    read(lin,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
c
c Read Input Parameters:
c
      read(lin,'(a512)',err=98) datafl
      call chknam(datafl,512)
      write(*,*) ' data file = ',datafl(1:40)

      read(lin,*,err=98) ixl,iyl,izl
      write(*,*) ' columns for X,Y,Z = ',ixl,iyl,izl

      read(lin,*,err=98) nvar
      write(*,*) ' number of variables = ',nvar
      backspace lin

      read(lin,*,err=98) j,(ivar(i),i=1,nvar)
      write(*,*) ' columns = ',(ivar(i),i=1,nvar)

      read(lin,*,err=98) tmin,tmax
      write(*,*) ' trimming limits = ',tmin,tmax

      read(lin,'(a512)',err=98) outfl
      call chknam(outfl,512)
      write(*,*) ' output file = ',outfl(1:40)

      read(lin,*,err=98) nlag
      write(*,*) ' number of lags = ',nlag
      if(nlag.lt.1)       stop 'nlag is too small: check parameters'

      read(lin,*,err=98) xlag
      write(*,*) ' lag distance = ',xlag
      if(xlag.le.0.0) stop 'xlag is too small: check parameter file'

      read(lin,*,err=98) xltol
      write(*,*) ' lag tolerance = ',xltol

      read(lin,*,err=98) ndir
      write(*,*) ' number of directions = ',ndir
      
      allocate (azm(ndir),stat = test)   
      allocate (atol(ndir),stat = test)
      allocate (bandwh(ndir),stat = test)
      allocate (dip(ndir),stat = test)
      allocate (dtol(ndir),stat = test)
      allocate (bandwd(ndir),stat = test)      

      if(ndir.lt.1)       stop 'ndir is too small: check parameters'

      do i=1,ndir
            read(lin,*,err=98) azm(i),atol(i),bandwh(i),
     +                         dip(i),dtol(i),bandwd(i)
            write(*,*) ' azm, atol, bandwh = ',azm(i),atol(i),bandwh(i)
            write(*,*) ' dip, dtol, bandwd = ',dip(i),dtol(i),bandwd(i)
            if(bandwh(i).lt.0.0) then
                  write(*,*) ' Horizontal bandwidth is too small!'
                  stop
            endif
            if(bandwd(i).lt.0.0) then
                  write(*,*) ' Vertical bandwidth is too small!'
                  stop
            endif
      end do

      read(lin,*,err=98) isill
      write(*,*) ' flag to standardize sills = ',isill

      read(lin,*,err=98) nvarg
      write(*,*) ' number of variograms = ',nvarg
      if(nvarg.lt.1)      stop 'nvarg is too small: check parameters'   

      mxdlv=ndir*(nlag+2)*nvarg

      allocate(reducedVariables(7,mxdlv,numThreads),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      

      allocate (dis(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (gam(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (hm(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (tm(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (hv(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (tv(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (np(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (ivtail(nvarg),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (ivhead(nvarg),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (ivtype(nvarg),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (names(nvar+nvarg),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
      
      ncut = 0
      do i=1,nvarg
            read(lin,*,err=98) ivtail(i),ivhead(i),ivtype(i)
            write(*,*) ' tail,head,type = ',
     +                   ivtail(i),ivhead(i),ivtype(i)
            if(ivtype(i).eq.9.or.ivtype(i).eq.10) then
                   ncut = ncut + 1
                   if(tmin.gt.0.0)stop'tmin interferes with indicators!'
                   if(tmax.le.1.0)stop'tmax interferes with indicators!'
                   backspace lin
                   read(lin,*,err=98) ii,jj,kk,cut(ncut)
                   if(ivtype(i).eq.9)  indflag(ncut) = 1
                   if(ivtype(i).eq.10) indflag(ncut) = 0
                   ivc(ncut) = ivtail(i)
                   ivtail(i) = nvar + ncut
                   ivhead(i) = nvar + ncut
                   write(names(nvar+ncut),140) ncut
 140               format('Indicator ',i2)
                   write(*,*) ' indicator threshold = ',cut(ncut)
            endif
      end do
      write(*,*)
      close(lin)
      MAXVAR = nvar + ncut
c
c Perform some quick error checking:
c
      if(xltol.le.0.0) then
            write(*,*) 'xltol is too small: resetting to xlag/2'
            xltol = 0.5*xlag
      endif
c
c Check to make sure the data file exists, then either read in the
c data or write an error message and stop:
c
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR data file ',datafl,' does not exist!'
            stop
      endif
c
c The data file exists so open the file and read in the header
c information. Initialize the storage that will be used to summarize
c the data found in the file:
c
      open(lin,file=datafl,status='OLD')
c
      read(lin,*,err=99)
      read(lin,*,err=99)       nvari
      do i=1,nvari
            read(lin,*)
      end do
      maxdat = 0
 22   read(lin,*,end=44,err=99)(var(j),j=1,nvari)
      maxdat = maxdat + 1
      go to 22
 44   continue
      write(*,*)'maxdat = ',maxdat
c
      allocate (vr(maxdat,MAXVAR),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
c
      allocate (vrmin(MAXVAR),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
c
      allocate (vrmax(MAXVAR),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
c
      allocate (sills(MAXVAR),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
c
      allocate (x(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
c
      allocate (y(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
c
      allocate (z(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
c
      rewind(lin) 
      read(lin,'(a40)',err=99) str
      read(lin,*,err=99)       nvari
      do i=1,nvari
            read(lin,'(a40)',err=99) str
            do iv=1,nvar
                  j=ivar(iv)
                  if(i.eq.j) names(iv) = str(1:12)
            end do
            num(i) = 0
            avg(i) = 0.0
            ssq(i) = 0.0
      end do
c
c Read all the data until the end of the file:
c
      nd = 0
 2    continue
      read(lin,*,end=9,err=99) (var(j),j=1,nvari)
      testdat = .false.
      do iv=1,nvar
            j=ivar(iv)
            if(var(j).ge.tmin.and.var(j).lt.tmax) testdat = .true.
      end do
      if(.not.testdat) go to 2
      nd = nd + 1
c
c Acceptable data, make sure there are not too many data:
c
      do iv=1,nvar
            j=ivar(iv)
            vr(nd,iv) = var(j)
            if(var(j).ge.tmin.and.var(j).lt.tmax) then
                  num(iv) = num(iv) + 1
                  avg(iv) = avg(iv) + dble(var(j))
                  ssq(iv) = ssq(iv) + dble(var(j)*var(j))
            endif
      end do
      if(ixl.le.0) then
            x(nd) = 0.0
      else
            x(nd) = var(ixl)
      endif
      if(iyl.le.0) then
            y(nd) = 0.0
      else
            y(nd) = var(iyl)
      endif
      if(izl.le.0) then
            z(nd) = 0.0
      else
            z(nd) = var(izl)
      endif
      go to 2
 9    continue
      close(lin)
c
c Compute the averages and variances as an error check for the user:
c
      do iv=1,nvar
            sills(iv) = -999.
            if(num(iv).gt.0) then
                  avg(iv)   = avg(iv)/dble(num(iv))
                  ssq(iv)   =(ssq(iv)/dble(num(iv)))-avg(iv)*avg(iv)
                  sills(iv) = ssq(iv)
                  write(*,*) 'Variable number ',iv
                  write(*,*) '  Number   = ',num(iv)
                  write(*,*) '  Average  = ',real(avg(iv))
                  write(*,*) '  Variance = ',real(ssq(iv))
            endif
      end do
c
c Construct Indicator Variables if necessary:
c
      do ic=1,ncut
            iv   = ivc(ic)
            jv   = nvar + ic
            ptot = 0.0
            p1   = 0.0
            do id=1,nd
                  if(vr(id,iv).le.tmin.or.vr(id,iv).gt.tmax) then
                        vr(id,jv) = tmin - EPSLON
                  else
                        if(indflag(ic).eq.1) then
                              if(vr(id,iv).lt.cut(ic)) then
                                    vr(id,jv) = 0.0
                              else
                                    vr(id,jv) = 1.0
                              endif
                              p1   = p1   + vr(id,jv)
                              ptot = ptot + 1.0
                        else
                              vr(id,jv) = 0.0
                              if(int(vr(id,iv)+0.5).eq.int(cut(ic)+0.5))
     +                        vr(id,jv) = 1.0
                              p1   = p1   + vr(id,jv)
                              ptot = ptot + 1.0
                        end if
                  end if
            end do
            p1        = p1 / max(ptot,1.0)
            sills(jv) = dble (p1*(1.0-p1))
      end do
c
c Establish minimums and maximums:
c
      do i=1,MAXVAR
            vrmin(i) =  1.0e21
            vrmax(i) = -1.0e21
      end do
      do id=1,nd
            do iv=1,nvar+ncut
                  if(vr(id,iv).ge.tmin.and.vr(id,iv).lt.tmax) then
                        if(vr(id,iv).lt.vrmin(iv)) vrmin(iv) = vr(id,iv)
                        if(vr(id,iv).gt.vrmax(iv)) vrmax(iv) = vr(id,iv)
                  end if
            end do
      end do
c
c Check on the variogams that were requested:
c
c      call check(vrmin,vrmax)

c
c Loop over all the variograms to be computed:
c
      write(*,*)
      do iv=1,nvarg
c
c Note the variogram type and the variables being used:
c
      it = abs(ivtype(iv))
      if(it.eq. 1) title(1:24) = 'Semivariogram          :'
      if(it.eq. 2) title(1:24) = 'Cross Semivariogram    :'
      if(it.eq. 3) title(1:24) = 'Covariance             :'
      if(it.eq. 4) title(1:24) = 'Correlogram            :'
      if(it.eq. 5) title(1:24) = 'General Relative       :'
      if(it.eq. 6) title(1:24) = 'Pairwise Relative      :'
      if(it.eq. 7) title(1:24) = 'Variogram of Logarithms:'
      if(it.eq. 8) title(1:24) = 'Semimadogram           :'
      if(it.eq. 9) title(1:24) = 'Indicator 1/2 Variogram:'
      if(it.eq.10) title(1:24) = 'Indicator 1/2 Variogram:'
      write(title(25:64),100) names(ivtail(iv)),names(ivhead(iv))
 100  format('  tail=',a12,' head=',a12)
      write(*,101) iv,title(1:64)
 101  format(' Variogram ',i2,1x,a64)
c
c Check for possible errors or inconsistencies:
c
      if(it.eq.2) then
            if(ivtail(iv).eq.ivhead(iv)) write(*,201)
 201        format('  WARNING: cross variogram with the same variable!')
      else if(it.eq.5) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,501)
            if(vrmin(ivtail(iv)).lt.0.0.and.vrmax(ivtail(iv)).gt.0.0)
     +            write(*,502)
            if(vrmin(ivhead(iv)).lt.0.0.and.vrmax(ivhead(iv)).gt.0.0)
     +            write(*,502)
 501        format('  WARNING: cross general relative variogram are',
     +             ' difficult to interpret!')
 502        format('  WARNING: there are both positive and negative',
     +             ' values - lag mean could be zero!')
      else if(it.eq.6) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,601)
            if(vrmin(ivtail(iv)).lt.0.0.and.vrmax(ivtail(iv)).gt.0.0)
     +            write(*,602)
            if(vrmin(ivhead(iv)).lt.0.0.and.vrmax(ivhead(iv)).gt.0.0)
     +            write(*,602)
 601        format('  WARNING: cross pairwise relative variogram are',
     +             ' difficult to interpret!')
 602        format('  WARNING: there are both positive and negative',
     +             ' values - pair means could be zero!')
      else if(it.eq.7) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,701)
            if(vrmin(ivtail(iv)).lt.0.0.or.vrmin(ivhead(iv)).lt.0.0)
     +      write(*,702)
 701        format('  WARNING: cross logarithmic variograms may be',
     +             ' difficult to interpret!')
 702        format('  WARNING: there are zero or negative',
     +             ' values - logarithm undefined!')
      else if(it.eq.8) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,801)
 801        format('  WARNING: cross rodograms may be difficult to',
     +             ' interpret!')
      else if(it.eq.9) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,901)
 901        format('  WARNING: cross madograms may be difficult to',
     +             ' interpret!')
      endif
c
c END Loop over all variograms:
c
      end do


      goto 1001
c
c Error in an Input File Somewhere:
c
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'

1001  print *,'Parameters ok.'

ccccccccccccccc gamv cccccccccccccccccccccc

c
c Define the distance tolerance if it isn't already:
c
      if(xltol.le.0.0) xltol = 0.5 * xlag
c
c Define the angles and tolerance for each direction:
c
      do id=1,ndir
c
c The mathematical azimuth is measured counterclockwise from EW and
c not clockwise from NS as the conventional azimuth is:
c
            azmuth     = (90.0-azm(id))*PI/180.0
            uvxazm(id) = cos(azmuth)
            uvyazm(id) = sin(azmuth)
            if(atol(id).le.0.0) then
                  csatol(id) = cos(45.0*PI/180.0)
            else
                  csatol(id) = cos(atol(id)*PI/180.0)
            endif
c
c The declination is measured positive down from vertical (up) rather
c than negative down from horizontal:
c
            declin     = (90.0-dip(id))*PI/180.0
            uvzdec(id) = cos(declin)
            uvhdec(id) = sin(declin)
            if(dtol(id).le.0.0) then
                  csdtol(id) = cos(45.0*PI/180.0)
            else
                  csdtol(id) = cos(dtol(id)*PI/180.0)
            endif
      end do
c
c Initialize the arrays for each direction, variogram, and lag:
c
      nsiz = ndir*nvarg*(nlag+2)
      do i=1,nsiz
            np(i)  = 0.
            dis(i) = 0.0
            gam(i) = 0.0
            hm(i)  = 0.0
            tm(i)  = 0.0
            hv(i)  = 0.0
            tv(i)  = 0.0
      end do
      dismxs = ((real(nlag) + 0.5 - EPSLON) * xlag) ** 2
c
c MAIN LOOP OVER ALL PAIRS:
c
      irepo = max(1,min((nd/10),1000))

#ifdef USE_MPI

      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, ierr)

      allocate(nptmp(mxdlv))
      allocate(distmp(mxdlv))
      allocate(gamtmp(mxdlv))
      allocate(hmtmp(mxdlv))
      allocate(tmtmp(mxdlv))
      allocate(hvtmp(mxdlv))
      allocate(tvtmp(mxdlv))

      if(myrank.eq.0)then
         nptmp(:)=0.0
         distmp(:)=0.0
         gamtmp(:)=0.0
         hmtmp(:)=0.0
         tmtmp(:)=0.0
         hvtmp(:)=0.0
         tvtmp(:)=0.0
      end if

#endif

      xlaginv=1.0/xlag
c$omp parallel default(firstprivate)
c$omp& shared(x,y,z,reducedVariables,vr)
#ifdef _OPENMP
      threadId = int(OMP_get_thread_num())+1 
      print *,'Thread ',threadId
#else
      threadId = 1
#endif

#ifdef USE_MPI
      blocksize=USE_MPI*numThreads
      nblocks =ceiling(real(nd-blocksize+1)/real(blocksize))
      do k=1,nblocks
         if( mod(k,mpisize).eq.myrank )then
            iinf=(k-1)*blocksize+1
            isup=(k)*blocksize
            if(isup.gt.nd)then
               isup=nd
            endif
c$omp do schedule(runtime) 
            do i=iinf,isup
#else
c$omp do schedule(runtime) 
      do i=1,nd
#endif

        if((int(i/irepo)*irepo).eq.i) write(*,103) i,nd,threadId,myrank
 103  format('   currently on ',i9,' of ',i9,' with thread ',i3,
     + ' and process ',i3)

c            xi=x(i)
c            yi=y(i)
c            zi=z(i)

      do 4 j=i,nd
c
c Definition of the lag corresponding to the current pair:
c
            dx  = x(j) - x(i)
            dy  = y(j) - y(i)
            dz  = z(j) - z(i)
c            dx  = x(j) - xi
c            dy  = y(j) - yi
c            dz  = z(j) - zi
            dxs = dx*dx
            dys = dy*dy
            dzs = dz*dz
            hs  = dxs + dys + dzs
            if(hs.gt.dismxs) go to 4
            if(hs.lt.0.0) hs = 0.0
            h   = sqrt(hs)
c
c Determine which lag this is and skip if outside the defined distance
c tolerance:
c
            if(h.le.EPSLON) then
                  lagbeg = 1
                  lagend = 1
            else
                  lagbeg = -1
                  lagend = -1
c                  do ilag=2,nlag+2
c                        if(h.ge.(xlag*real(ilag-2)-xltol).and.
c     +                     h.le.(xlag*real(ilag-2)+xltol)) then
c                              if(lagbeg.lt.0) lagbeg = ilag 
c                              lagend = ilag 
c                        end if
c                  end do
                  liminf=(ceiling((h-xltol)*xlaginv)+2)
                  limsup=(floor((h+xltol)*xlaginv)+2)
                  do ilag=liminf,limsup
                     if(lagbeg.lt.0) lagbeg = ilag
                     lagend = ilag
                  end do
                  if(lagend.lt.0) go to 4
            endif

c
c Definition of the direction corresponding to the current pair. All
c directions are considered (overlapping of direction tolerance cones
c is allowed):
c
            do 5 id=1,ndir
c
c Check for an acceptable azimuth angle:
c
                  dxy = sqrt(max((dxs+dys),0.0))
                  if(dxy.lt.EPSLON) then
                        dcazm = 1.0
                  else
                        dcazm = (dx*uvxazm(id)+dy*uvyazm(id))/dxy
                  endif
                  if(abs(dcazm).lt.csatol(id)) go to 5

c
c Check the horizontal bandwidth criteria (maximum deviation 
c perpendicular to the specified direction azimuth):
c
                  band = uvxazm(id)*dy - uvyazm(id)*dx
                  if(abs(band).gt.bandwh(id)) go to 5
c
c Check for an acceptable dip angle:
c
                  if(dcazm.lt.0.0) dxy = -dxy
                  if(lagbeg.eq.1) then
                        dcdec = 0.0
                  else
                        dcdec = (dxy*uvhdec(id)+dz*uvzdec(id))/h
                        if(abs(dcdec).lt.csdtol(id)) go to 5
                  endif
c
c Check the vertical bandwidth criteria (maximum deviation perpendicular
c to the specified dip direction):
c
                  band = uvhdec(id)*dz - uvzdec(id)*dxy
                  if(abs(band).gt.bandwd(id)) go to 5

c
c Check whether or not an omni-directional variogram is being computed:
c
                  omni = .false.
                  if(atol(id).ge.90.0) omni = .true.
c
c This direction is acceptable - go ahead and compute all variograms:
c

                  do 6 iv=1,nvarg
c
c For this variogram, sort out which is the tail and the head value:
c
                      it = ivtype(iv)
                      if(dcazm.ge.0.0.and.dcdec.ge.0.0) then
                            ii = ivtail(iv)
                            vrh   = vr(i,ii)
                            ii = ivhead(iv)
                            vrt   = vr(j,ii)
                            if(omni.or.it.eq.2) then
                                  ii    = ivhead(iv)
                                  vrtpr = vr(i,ii)
                                  ii    = ivtail(iv)
                                  vrhpr = vr(j,ii)
                            endif
                      else
                            ii = ivtail(iv)
                            vrh   = vr(j,ii)
                            ii = ivhead(iv)
                            vrt   = vr(i,ii)
                            if(omni.or.it.eq.2) then
                                  ii    = ivhead(iv)
                                  vrtpr = vr(j,ii)
                                  ii    = ivtail(iv)
                                  vrhpr = vr(i,ii)
                            endif
                      endif
c
c Reject this pair on the basis of missing values:
c
                      if(vrt.lt.tmin.or.vrh.lt.tmin.or.
     +                   vrt.gt.tmax.or.vrh.gt.tmax) go to 6
                      if(it.eq.2.and.(vrtpr.lt.tmin.or.vrhpr.lt.tmin.or.
     +                                vrtpr.gt.tmax.or.vrhpr.gt.tmax))
     +                                               go to 6
c
c             COMPUTE THE APPROPRIATE "VARIOGRAM" MEASURE
c
c The Semivariogram:
c
      if(it.eq.1.or.it.eq.5.or.it.ge.9) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               np(ii)  = np(ii)  + 1.
               dis(ii) = dis(ii) + dble(h)
               tm(ii)  = tm(ii)  + dble(vrt)
               hm(ii)  = hm(ii)  + dble(vrh)
               gam(ii) = gam(ii) + dble((vrh-vrt)*(vrh-vrt))
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           gam(ii) = gam(ii) + dble((vrhpr-vrtpr)*
     +                                              (vrhpr-vrtpr))
                     endif
               endif
            end do
c
c The Traditional Cross Semivariogram:
c
      else if(it.eq.2) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               np(ii)  = np(ii)  + 1.
               dis(ii) = dis(ii) + dble(h)
               tm(ii)  = tm(ii)  + dble(0.5*(vrt+vrtpr))
               hm(ii)  = hm(ii)  + dble(0.5*(vrh+vrhpr))
               gam(ii) = gam(ii) + dble((vrhpr-vrh)*(vrt-vrtpr))
            end do
c
c The Covariance:
c
      else if(abs(it).eq.3) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               np(ii)  = np(ii)  + 1.
               dis(ii) = dis(ii) + dble(h)
               tm(ii)  = tm(ii)  + dble(vrt)
               hm(ii)  = hm(ii)  + dble(vrh)
               gam(ii) = gam(ii) + dble(vrh*vrt)
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           gam(ii) = gam(ii) + dble(vrhpr*vrtpr)
                     endif
               endif
            end do
c
c The Correlogram:
c
      else if(abs(it).eq.4) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               np(ii)  = np(ii)  + 1.
               dis(ii) = dis(ii) + dble(h)
               tm(ii)  = tm(ii)  + dble(vrt)
               hm(ii)  = hm(ii)  + dble(vrh)
               hv(ii)  = hv(ii)  + dble(vrh*vrh)
               tv(ii)  = tv(ii)  + dble(vrt*vrt)
               gam(ii) = gam(ii) + dble(vrh*vrt)
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           hv(ii)  = hv(ii)  + dble(vrhpr*vrhpr)
                           tv(ii)  = tv(ii)  + dble(vrtpr*vrtpr)
                           gam(ii) = gam(ii) + dble(vrhpr*vrtpr)
                     endif
               endif
            end do
c
c The Pairwise Relative:
c
      else if(it.eq.6) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               if(abs(vrt+vrh).gt.EPSLON) then
                     np(ii)  = np(ii)  + 1.
                     dis(ii) = dis(ii) + dble(h)
                     tm(ii)  = tm(ii)  + dble(vrt)
                     hm(ii)  = hm(ii)  + dble(vrh)
                     gamma   = 2.0*(vrt-vrh)/(vrt+vrh)
                     gam(ii) = gam(ii) + dble(gamma*gamma)
               endif
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                     if(abs(vrtpr+vrhpr).gt.EPSLON) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           gamma   = 2.0*(vrt-vrh)/(vrt+vrh)
                           gam(ii) = gam(ii) + dble(gamma*gamma)
                     endif
                     endif
               endif
            enddo
c
c Variogram of Logarithms:
c
      else if(it.eq.7) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               if(vrt.gt.EPSLON.and.vrh.gt.EPSLON) then
                     np(ii)  = np(ii)  + 1.
                     dis(ii) = dis(ii) + dble(h)
                     tm(ii)  = tm(ii)  + dble(vrt)
                     hm(ii)  = hm(ii)  + dble(vrh)
                     gamma   = alog(vrt)-alog(vrh)
                     gam(ii) = gam(ii) + dble(gamma*gamma)
               endif
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                     if(vrtpr.gt.EPSLON.and.vrhpr.gt.EPSLON) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           gamma   = alog(vrt)-alog(vrh)
                           gam(ii) = gam(ii) + dble(gamma*gamma)
                     endif
                     endif
               endif
            end do
c
c Madogram:
c
      else if(it.eq.8) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               np(ii)  = np(ii)  + 1.
               dis(ii) = dis(ii) + dble(h)
               tm(ii)  = tm(ii)  + dble(vrt)
               hm(ii)  = hm(ii)  + dble(vrh)
               gam(ii) = gam(ii) + dble(abs(vrh-vrt))
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           gam(ii) = gam(ii) + dble(abs(vrhpr-vrtpr))
                     endif
               endif
            end do
      endif
c
c Finish loops over variograms, directions, and the double data loops:
c
 6              continue
 5          continue
 4    continue

      end do
c$omp end do nowait


#ifdef USE_MPI
         end if
      end do
#endif


#ifdef _OPENMP

c$omp barrier
      reducedVariables(1,:,threadId)=dis(:)
      reducedVariables(2,:,threadId)=gam(:)
      reducedVariables(3,:,threadId)=np(:)
      reducedVariables(4,:,threadId)=hm(:)
      reducedVariables(5,:,threadId)=tm(:)
      reducedVariables(6,:,threadId)=hv(:)
      reducedVariables(7,:,threadId)=tv(:)
#endif


c$omp end parallel

#ifdef _OPENMP
      dis(:)=0.0 
      gam(:)=0.0 
      np(:)=0.0  
      hm(:)=0.0  
      tm(:)=0.0  
      hv(:)=0.0  
      tv(:)=0.0       
      do ii=1,numThreads 
         do jj=1,mxdlv
            dis(jj) = dis(jj) + reducedVariables(1,jj,ii)
            gam(jj) = gam(jj) + reducedVariables(2,jj,ii)
            np(jj)  = np(jj)  + reducedVariables(3,jj,ii)
            hm(jj)  = hm(jj)  + reducedVariables(4,jj,ii)
            tm(jj)  = tm(jj)  + reducedVariables(5,jj,ii)
            hv(jj)  = hv(jj)  + reducedVariables(6,jj,ii)
            tv(jj)  = tv(jj)  + reducedVariables(7,jj,ii)
         end do
      end do
#endif

#ifdef USE_MPI
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(dis, distmp, mxdlv, MPI_DOUBLE_PRECISION, 
     + MPI_SUM,0, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(gam, gamtmp, mxdlv, MPI_DOUBLE_PRECISION, 
     + MPI_SUM,0, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(np, nptmp, mxdlv, MPI_DOUBLE_PRECISION, 
     + MPI_SUM,0, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(hm, hmtmp, mxdlv, MPI_DOUBLE_PRECISION, 
     + MPI_SUM,0, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(tm, tmtmp, mxdlv, MPI_DOUBLE_PRECISION, 
     + MPI_SUM,0, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(hv, hvtmp, mxdlv, MPI_DOUBLE_PRECISION, 
     + MPI_SUM,0, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(tv, tvtmp, mxdlv, MPI_DOUBLE_PRECISION, 
     + MPI_SUM,0, MPI_COMM_WORLD, ierr)

      if(myrank.eq.0)then
         do jj=1,mxdlv
            dis(jj) = distmp(jj) 
            gam(jj) = gamtmp(jj) 
            np(jj)  = nptmp(jj)  
            hm(jj)  = hmtmp(jj)  
            tm(jj)  = tmtmp(jj)  
            hv(jj)  = hvtmp(jj)  
            tv(jj)  = tvtmp(jj)  
         end do
      end if
#endif

#ifdef USE_MPI
      if(myrank.eq.0)then
#endif
c
c Get average values for gam, hm, tm, hv, and tv, then compute
c the correct "variogram" measure:
c
      do 7 id=1,ndir
      do 7 iv=1,nvarg
      do 7 il=1,nlag+2
            i = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
            if(np(i).le.0.) go to 7
            rnum   = np(i)
            dis(i) = dis(i) / dble(rnum)
            gam(i) = gam(i) / dble(rnum)
            hm(i)  = hm(i)  / dble(rnum)
            tm(i)  = tm(i)  / dble(rnum)
            hv(i)  = hv(i)  / dble(rnum)
            tv(i)  = tv(i)  / dble(rnum)
            it = ivtype(iv)
c
c Attempt to standardize:
c
            if(isill.eq.1) then
                  if(ivtail(iv).eq.ivhead(iv)) then
                        iii = ivtail(iv)
                        if((it.eq.1.or.it.ge.9).and.sills(iii).gt.0.0)
     +                    gam(i) = gam(i) / sills(iii)
                  end if
            end if
c
c  1. report the semivariogram rather than variogram
c  2. report the cross-semivariogram rather than variogram
c  3. the covariance requires "centering"
c  4. the correlogram requires centering and normalizing
c  5. general relative requires division by lag mean
c  6. report the semi(pairwise relative variogram)
c  7. report the semi(log variogram)
c  8. report the semi(madogram)
c
            if(it.eq.1.or.it.eq.2) then
                  gam(i) = 0.5 * gam(i)
            else if(abs(it).eq.3) then
                  gam(i) = gam(i) - hm(i)*tm(i)
                  if(it.lt.0) then
                        if(sills(ivtail(iv)).lt.0.0.or.
     +                     sills(ivhead(iv)).lt.0.0) then
                              gam(i) = -999.0
                        else
                              variance = ( sqrt(sills(ivtail(iv)))
     +                                 *   sqrt(sills(ivhead(iv))) )
                              gam(i) = variance - gam(i)
                        end if
                  end if
            else if(abs(it).eq.4) then
                  hv(i)  = hv(i)-hm(i)*hm(i)
                  if(hv(i).lt.0.0) hv(i) = 0.0
                  hv(i)  = dsqrt(hv(i))
                  tv(i)  = tv(i)-tm(i)*tm(i)
                  if(tv(i).lt.0.0) tv(i) = 0.0
                  tv(i)  = dsqrt(tv(i))
                  if((hv(i)*tv(i)).lt.EPSLON) then
                        gam(i) = 0.0
                  else
                        gam(i) =(gam(i)-hm(i)*tm(i))/(hv(i)*tv(i))
                  endif
                          if(it.lt.0) gam(i) = 1.0 - gam(i)
c
c Square "hv" and "tv" so that we return the variance:
c
                  hv(i)  = hv(i)*hv(i)
                  tv(i)  = tv(i)*tv(i)
            else if(it.eq.5) then
                  htave  = 0.5*(hm(i)+tm(i))
                  htave  = htave   *   htave
                  if(htave.lt.EPSLON) then
                        gam(i) = 0.0
                  else
                        gam(i) = gam(i)/dble(htave)
                  endif
            else if(it.ge.6) then
                  gam(i) = 0.5 * gam(i)
            endif
 7    continue

ccccccccccccccc writeout cccccccccccccccccccccc

c
c Loop over all the variograms that have been computed:
c
      open(lout,file=outfl,status='UNKNOWN')
      do iv=1,nvarg
c
c Construct a title that reflects the variogram type and the variables
c that were used to calculate the variogram:
c
      it = abs(ivtype(iv))
      if(it.eq. 1) titlewr(1:24) = 'Semivariogram           '
      if(it.eq. 2) titlewr(1:24) = 'Cross Semivariogram     '
      if(it.eq. 3) titlewr(1:24) = 'Covariance              '
      if(it.eq. 4) titlewr(1:24) = 'Correlogram             '
      if(it.eq. 5) titlewr(1:24) = 'General Relative        '
      if(it.eq. 6) titlewr(1:24) = 'Pairwise Relative       '
      if(it.eq. 7) titlewr(1:24) = 'Variogram of Logarithms '
      if(it.eq. 8) titlewr(1:24) = 'Semimadogram            '
      if(it.eq. 9) titlewr(1:24) = 'Indicator 1/2 Variogram '
      if(it.eq.10) titlewr(1:24) = 'Indicator 1/2 Variogram '
      write(titlewr(25:62),1100) names(ivtail(iv)),names(ivhead(iv))
 1100 format('tail:',a12,' head:',a12)
c
c Loop over all the directions (note the direction in the title):
c
      do id=1,ndir
            write(titlewr(62:74),1101) id
 1101       format('direction ',i2)
            write(lout,'(a74)') titlewr(1:74)
c
c Write out all the lags:
c
            do il=1,nlag+2
                  i = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
                  nump = int(np(i))
                  if(it.eq.4) then
                        write(lout,1102) il,dis(i),gam(i),nump,
     +                                  hm(i),tm(i),hv(i),tv(i)
                  else
                        write(lout,1102) il,dis(i),gam(i),nump,
     +                                  hm(i),tm(i)
                  endif
 1102             format(1x,i3,1x,f12.3,1x,f12.5,1x,i8,4(1x,f14.5))
            end do
      end do
c
c End loop over variograms
c
      end do
c
c Finished:
c
      close(lout)

#ifdef USE_MPI
      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr);
      deallocate(nptmp)
      deallocate(distmp)
      deallocate(gamtmp)
      deallocate(hmtmp)
      deallocate(tmtmp)
      deallocate(hvtmp)
      deallocate(tvtmp)

#endif

c
c Finished:
c

      write(*,9998) VERSION
 9998 format(/' GAMV Version: ',f5.3, ' Finished'/)

c#ifdef TRACE
c      call extrae_fini()
c#endif

#ifdef USE_MPI
      call MPI_FINALIZE(ierr) 
#endif
      stop
      end
 
 


      subroutine makepar
c-----------------------------------------------------------------------
c
c                      Write a Parameter File
c                      **********************
c
c
c
c-----------------------------------------------------------------------
      lun = 99
      open(lun,file='gamv.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for GAMV',/,
     +       '                  *******************',/,/,
     +       'START OF PARAMETERS:')

      write(lun,11)
 11   format('../data/cluster.dat               ',
     +       '-file with data')
      write(lun,12)
 12   format('1   2   0                         ',
     +       '-   columns for X, Y, Z coordinates')
      write(lun,13)
 13   format('2   3   4                         ',
     +       '-   number of variables,col numbers')
      write(lun,14)
 14   format('-1.0e21     1.0e21                ',
     +       '-   trimming limits')
      write(lun,15)
 15   format('gamv.out                          ',
     +       '-file for variogram output')
      write(lun,16)
 16   format('10                                ',
     +       '-number of lags')
      write(lun,17)
 17   format('5.0                               ',
     +       '-lag separation distance')
      write(lun,18)
 18   format('3.0                               ',
     +       '-lag tolerance')
      write(lun,19)
 19   format('3                                 ',
     +       '-number of directions')
      write(lun,20)
 20   format('0.0  90.0 50.0   0.0  90.0  50.0  ',
     +       '-azm,atol,bandh,dip,dtol,bandv')
      write(lun,21)
 21   format('0.0  22.5 25.0   0.0  22.5  25.0  ',
     +       '-azm,atol,bandh,dip,dtol,bandv')
      write(lun,22)
 22   format('90.  22.5 25.0   0.0  22.5  25.0  ',
     +       '-azm,atol,bandh,dip,dtol,bandv')
      write(lun,23)
 23   format('0                                 ',
     +       '-standardize sills? (0=no, 1=yes)')
      write(lun,24)
 24   format('3                                 ',
     +       '-number of variograms')
      write(lun,25)
 25   format('1   1   1                         ',
     +       '-tail var., head var., variogram type')
      write(lun,26)
 26   format('1   2   2                         ',
     +       '-tail var., head var., variogram type')
      write(lun,27)
 27   format('2   2   1                         ',
     +       '-tail var., head var., variogram type')
      write(lun,40)
 40   format(//,'type 1 = traditional semivariogram',/,
     +          '     2 = traditional cross semivariogram',/,
     +          '     3 = covariance',/,
     +          '     4 = correlogram',/,
     +          '     5 = general relative semivariogram',/,
     +          '     6 = pairwise relative semivariogram',/,
     +          '     7 = semivariogram of logarithms',/,
     +          '     8 = semimadogram',/,
     +          '     9 = indicator semivariogram - continuous',/,
     +          '     10= indicator semivariogram - categorical')

      close(lun)
      return
      end
