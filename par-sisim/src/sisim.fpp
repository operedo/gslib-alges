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
c           Conditional Simulation of a 3-D Rectangular Grid
c           ************************************************
c
c The output file will be a GEOEAS file containing the simulated values
c The file is ordered by x,y,z, and then simulation (i.e., x cycles
c fastest, then y, then z, then simulation number).
c
c
c
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

      include  'sisim.inc'



      real,allocatable   :: x(:),y(:),z(:),vr(:,:),close(:),actloc(:),
     +                      tmpdat(:),order(:),gcut(:),sim(:),tmp(:),   
     +                      gcdf(:),covtab(:,:,:,:),cnodex(:),cnodey(:),
     +                      cnodez(:),cnodev(:),cnodet(:),vra(:),
     +                      thres(:),cdf(:),ccdf(:),ccdfo(:),beez(:),
     +                      c0(:),cmax(:),cc(:),aa(:),ang1(:),ang2(:),
     +                      ang3(:),anis1(:),anis2(:),aviol(:),xviol(:),
     +                      aainv(:)
      real*8,allocatable  :: r(:),rr(:),s(:),a(:)
      integer,allocatable :: ixnode(:),iynode(:),iznode(:),nisb(:),
     +                       icnode(:),ixsbtosr(:),iysbtosr(:),
     +                       izsbtosr(:),it(:),nst(:),nviol(:)

      integer, allocatable:: aclose(:)

      real*8,allocatable  :: rotmat(:,:,:)
      logical,allocatable :: atnode(:),softdat(:)


      integer nodmax
      integer nctx,ncty,nctz,nlooku,ncnode,maxsec


      parameter(MV=20)
      real      var(MV)
      real*8    p,acorni
      integer,allocatable :: ivrs(:)
      integer   test
      character datafl*512,tabfl*512,softfl*512,outfl*512,dbgfl*512,
     +          str*512,title*80
      logical   testfl
      real      ntviol,atviol


      integer finLoop

c#ifdef _OPENMP
      integer MAXTHREADS
      parameter (MAXTHREADS=60)
      integer BUFFSIZE
      parameter(BUFFSIZE=2)
      integer threadId,numThreads,ithread
      integer loutThreads(MAXTHREADS)
      integer ldbgThreads(MAXTHREADS)
      character outflThreads(43,MAXTHREADS) , outfltmp*43
      character dbgflThreads(43,MAXTHREADS) , dbgfltmp*43
c#endif

c MPI variables
      integer ierr, myrank,mpisize,nloopIni,nloopFin

#ifdef USE_MPI
      call MPI_INIT(ierr) 
      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, ierr)
      print *,'Hello from process ',myrank,'/',mpisize
#else
      myrank=0
      mpisize=1
#endif



c#ifdef TRACE
c      call extrae_init()
c#endif

      numThreads=1
#ifdef _OPENMP
c$omp parallel
      numThreads = OMP_get_num_threads()
c$omp end parallel
#endif


c
c Fortran unit numbers needed:
c
      lin  = 1
      lout = 2
#ifdef _OPENMP
      do i=1,MAXTHREADS
#ifdef USE_MPI
            loutThreads(i)=lout*100+i+myrank*MAXTHREADS
#else
            loutThreads(i)=lout*100+i
#endif
      end do
#else
      loutThreads(1)=lout
#endif
      ldbg = 3
#ifdef _OPENMP
      do i=1,MAXTHREADS
#ifdef USE_MPI
            ldbgThreads(i)=ldbg*100+i+myrank*MAXTHREADS
#else
            ldbgThreads(i)=ldbg*100+i
#endif
      end do
#else
      ldbgThreads(1)=ldbg
#endif

ccccccccccccccc readparms cccccccccccccccccccccc

c
c Note VERSION number:
c
      write(*,9999) VERSION
 9999 format(/' SISIM Version: ',f5.3/)

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
      if(str(1:1).eq.' ') str(1:20) = 'sisim.par           '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'sisim.par           ') then
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

      read(lin,*,err=98) ivtype
      write(*,*) ' variable type (1=continuous, 0=categorical)= ',ivtype

      read(lin,*,err=98) ncut
      write(*,*) ' number of thresholds / categories = ',ncut
c
c Find the needed parameters:
c
      MAXCUT = ncut
      MAXROT = MAXCUT * MAXNST + 1
      MXCUT = MAXCUT + 1
c
c Allocate the needed memory:
c34
      allocate(thres(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 34: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c35
      allocate(cdf(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 35: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c36
      allocate(ccdf(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 36: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c37
      allocate(ccdfo(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 37: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c38
      allocate(beez(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 38: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c39
      allocate(c0(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 39: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c40
      allocate(cmax(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 40: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c41
      allocate(cc(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 41: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c42
      allocate(aa(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 42: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if

      allocate(aainv(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 42: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if

c43
      allocate(ang1(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 43: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c44
      allocate(ang2(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 44: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c45
      allocate(ang3(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 45: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c46
      allocate(anis1(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 46: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c47
      allocate(anis2(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 47: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c48
      allocate(aviol(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 48: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c49
      allocate(xviol(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 49: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c50
      allocate(it(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 50: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c51
      allocate(nst(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 51: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c52
      allocate(nviol(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 52: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c53
      allocate(rotmat(MAXROT,3,3),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 53: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c54
      allocate(ivrs(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 54: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c
      read(lin,*,err=98) (thres(i),i=1,ncut)
      write(*,*) ' thresholds / categories = ',(thres(i),i=1,ncut)

      read(lin,*,err=98) (cdf(i),i=1,ncut)
      write(*,*) ' global cdf / pdf        = ',(cdf(i),i=1,ncut)

      read(lin,'(a512)',err=98) datafl
      call chknam(datafl,512)
      write(*,*) ' data file = ',datafl(1:40)

      read(lin,*,err=98) ixl,iyl,izl,ivrl
      write(*,*) ' input columns = ',ixl,iyl,izl,ivrl

      read(lin,'(a512)',err=98) softfl
      call chknam(softfl,512)
      write(*,*) ' soft data file = ',softfl(1:40)
      inquire(file=softfl,exist=testfl)

      if(testfl) then
            read(lin,*,err=98) ixs,iys,izs,(ivrs(i),i=1,ncut)
            write(*,*) ' columns = ',ixs,iys,izs,(ivrs(i),i=1,ncut)
            read(lin,*,err=98) imbsim
            write(*,*) ' Markov-Bayes simulation = ',imbsim
            if(imbsim.eq.1) then
                  read(lin,*,err=98) (beez(i),i=1,ncut)
            else
                  read(lin,*,err=98)
            end if
      else
            read(lin,*,err=98)
            read(lin,*,err=98)
            read(lin,*,err=98)
      end if

      read(lin,*,err=98) tmin,tmax
      write(*,*) ' trimming limits      ',tmin,tmax

      read(lin,*,err=98) zmin,zmax
      write(*,*) ' data limits (tails)  ',zmin,zmax

      read(lin,*,err=98) ltail,ltpar
      write(*,*) ' lower tail = ',ltail,ltpar

      read(lin,*,err=98) middle,mpar
      write(*,*) ' middle = ',middle,mpar

      read(lin,*,err=98) utail,utpar
      write(*,*) ' upper tail = ',utail,utpar

      read(lin,'(a512)',err=98) tabfl
      call chknam(tabfl,512)
      write(*,*) ' file for tab. quant. ',tabfl(1:40)

      read(lin,*,err=98) itabvr,itabwt
      write(*,*) ' columns for vr wt = ',itabvr,itabwt

      read(lin,*,err=98) idbg
      write(*,*) ' debugging level = ',idbg

      read(lin,'(a512)',err=98) dbgfl
      call chknam(dbgfl,512)
      write(*,*) ' debugging file = ',dbgfl(1:40)

#ifdef _OPENMP
      do i=1,numThreads
            write(dbgfltmp,"(A40,I3)") trim(dbgfl(1:40)),
     + ldbgThreads(i)
            dbgfltmp = adjustl(trim(dbgfltmp))
            write(*,*) ' debugging file (threads) = ',
     + dbgfltmp

      end do
#endif

      read(lin,'(a512)',err=98) outfl
      call chknam(outfl,512)
      write(*,*) ' output file = ',outfl(1:40)

#ifdef _OPENMP
      do i=1,numThreads
            write(outfltmp,"(A40,I3)") trim(outfl(1:40)),
     + loutThreads(i)
            outfltmp = adjustl(trim(outfltmp))
            write(*,*) ' output file (threads) = ',
     + outfltmp

      end do
#endif

      read(lin,*,err=98) nsim
      write(*,*) ' number of simulations = ',nsim

      read(lin,*,err=98) nx,xmn,xsiz
      write(*,*) ' X grid specification = ',nx,xmn,xsiz

      read(lin,*,err=98) ny,ymn,ysiz
      write(*,*) ' Y grid specification = ',ny,ymn,ysiz

      read(lin,*,err=98) nz,zmn,zsiz
      write(*,*) ' Z grid specification = ',nz,zmn,zsiz
      nxy  = nx*ny
      nxyz = nx*ny*nz

      read(lin,*,err=98) ixv(1)
      write(*,*) ' random number seed = ',ixv(1)
      do i=1,1000
             p = acorni(idum)
      end do

      read(lin,*,err=98) ndmax
      write(*,*) ' ndmax = ',ndmax

      read(lin,*,err=98) nodmax
      write(*,*) ' max prev sim nodes = ',nodmax

      read(lin,*,err=98) maxsec
      write(*,*) ' max soft indicator data = ',maxsec

      read(lin,*,err=98) sstrat
      write(*,*) ' search strategy = ',sstrat

      read(lin,*,err=98) mults,nmult
      write(*,*) ' multiple grid search flag = ',mults,nmult

      read(lin,*,err=98) noct
      write(*,*) ' max per octant = ',noct

      read(lin,*,err=98) radius,radius1,radius2
      write(*,*) ' search radii = ',radius,radius1,radius2
      if(radius.lt.EPSLON) stop 'radius must be greater than zero'
      radsqd = radius  * radius
      sanis1 = radius1 / radius
      sanis2 = radius2 / radius

      read(lin,*,err=98) sang1,sang2,sang3
      write(*,*) ' search anisotropy angles = ',sang1,sang2,sang3

      read(lin,*,err=98) mxctx,mxcty,mxctz
      write(*,*) ' size of covariance lookup = ',mxctx,mxcty,mxctz
      
      read(lin,*,err=98) mik,cutmik
      write(*,*) ' median IK switch = ',mik,cutmik

      read(lin,*,err=98) ktype
      write(*,*) ' kriging type switch = ',ktype

c
c Output now goes to debugging file:
c
      open(ldbg,file=dbgfl,status='UNKNOWN')

#ifdef _OPENMP
c each thread will write in the files file.dbg401, file.dbg402, so on.
c up to MAXTHREADS
      do i=1,numThreads
            write(dbgfltmp,"(A40,I3)") trim(dbgfl(1:40)),
     + ldbgThreads(i)
            dbgfltmp = adjustl(trim(dbgfltmp))
            open(ldbgThreads(i),file=dbgfltmp,status='UNKNOWN')
       end do
#endif


      do i=1,ncut
            read(lin,*,err=98) nst(i),c0(i)
#ifdef DEBUG
            if(ivtype.eq.0)
     +      write(ldbg,100)  i,thres(i),cdf(i),nst(i),c0(i)
            if(ivtype.eq.1)
     +      write(ldbg,101)  i,thres(i),cdf(i),nst(i),c0(i)
#endif
            if(nst(i).gt.MAXNST) stop 'nst is too big'
            istart = 1 + (i-1)*MAXNST
            do j=1,nst(i)
                  index = istart + j - 1
                  read(lin,*,err=98) it(index),cc(index),ang1(index),
     +                               ang2(index),ang3(index)
                  if(it(index).eq.3) STOP 'Gaussian Model Not Allowed!'
                  read(lin,*,err=98) aa(index),aa1,aa2

                  aainv(index)=1.0/aa(index)

#ifdef DEBUG
                  write(ldbg,102)  j,it(index),aa(index),cc(index)
#endif
                  anis1(index) = aa1 / max(EPSLON,aa(index))
                  anis2(index) = aa2 / max(EPSLON,aa(index))
#ifdef DEBUG
                  write(ldbg,103) ang1(index),ang2(index),ang3(index),
     +                            anis1(index),anis2(index)
#endif
            end do
      end do
      close(lin)
c
c Find the needed parameters:
c
      MAXX = nx
      MAXY = ny
      MAXZ = nz
      MXYZ = MAXX * MAXY * MAXZ
      MAXCTX = mxctx
      MAXCTY = mxcty
      MAXCTZ = mxctz
      MAXCXY = MAXCTX * MAXCTY
      MAXXYZ = MAXCTX * MAXCTY * MAXCTZ
      MAXSAM = ndmax 
      MAXNOD = nodmax
      MAXKR1 = 2 * MAXNOD + 2 * MAXSAM + 1
      MAXKR2 = MAXKR1 * MAXKR1
      MAXSBX = 1
      if(nx.gt.1)then
            MAXSBX = int(nx/2)
            if(MAXSBX.gt.50)MAXSBX=50
      end if
c
      MAXSBY = 1
      if(ny.gt.1)then
            MAXSBY = int(ny/2)
            if(MAXSBY.gt.50)MAXSBY=50
      end if
c
      MAXSBZ = 1
      if(nz.gt.1)then
            MAXSBZ = int(nz/2)
            if(MAXSBZ.gt.50)MAXSBZ=50
      end if
c
      MAXSB = MAXSBX * MAXSBY * MAXSBZ
      MXSXY = 4 * MAXSBX * MAXSBY
      MXSX = 2 * MAXSBX
      av = 0.0
      ss = 0.0
c
c Find the paramater MAXDAT:
c
      MAXDAT = 1
      inquire(file=datafl,exist=testfl)
      if(testfl) then
            open(lin,file=datafl,status='OLD')
            read(lin,*,err=99)
            read(lin,*,err=99)       nvari
            do i=1,nvari
            read(lin,'()',err=99)
            end do
            MAXDAT = 0
 55         read(lin,*,end=66,err=99) (var(j),j=1,nvari)
            if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax)go to 55
            MAXDAT = MAXDAT + 1
            go to 55
 66         continue
            rewind(lin)
            close(lin)
      end if
      inquire(file=softfl,exist=testfl)
      if(testfl) then
            open(lin,file=softfl,status='OLD')
            read(lin,*,err=97)
            read(lin,*,err=97) nvari
            do i=1,nvari
            read(lin,'()',err=99)
            end do
 56         read(lin,*,end=67,err=99) (var(j),j=1,nvari)
            MAXDAT = MAXDAT + 1
            go to 56
 67         continue
            close(lin)
      endif
c
c Find the paramater MAXTAB:
c
      MAXTAB = ncut
      inquire(file=tabfl,exist=testfl)
      if(testfl) then
            open(lin,file=tabfl,status='OLD')
            read(lin,*,err=99)
            read(lin,*,err=99)       nvari
            do i=1,nvari
                  read(lin,'()',err=99)
            end do
            MAXTAB = 0
 77         read(lin,*,end=88,err=99) (var(j),j=1,nvari)
            if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax)go to 77
            MAXTAB = MAXTAB + 1
            go to 77
 88         continue
            rewind(lin)
            close(lin)
      end if
c
c Allocate the needed memory:
c1
      allocate(x(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 1: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c2
      allocate(y(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 2: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c3
      allocate(z(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 3: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c4
      allocate(vr(MAXDAT,MXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 4: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c5

      allocate(aclose(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 5: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if

      allocate(close(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 5: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c6
      allocate(actloc(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 6: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c7
      allocate(tmpdat(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 7: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c8
      allocate(sim(MXYZ),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 8: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if

c9
      MAXORD = MXYZ
      if(MXYZ.lt.MAXCXY) MAXORD=MAXCXY
      allocate(order(MAXORD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 9: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c10
      allocate(tmp(MAXORD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 10: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c11
      allocate(gcut(MAXTAB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 11: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c12
      allocate(gcdf(MAXTAB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 12: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c13
      allocate(covtab(MAXCTX,MAXCTY,MAXCTZ,MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 13: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c14
      allocate(cnodex(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 14: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c15
      allocate(cnodey(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 15: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c16
      allocate(cnodez(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 16: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c17
      allocate(cnodev(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 17: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c18
      allocate(cnodet(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 18: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c19
      allocate(vra(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 19: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c20
      allocate(r(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 20: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c21
      allocate(rr(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 21: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c22
      allocate(s(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 22: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c23
      allocate(a(MAXKR2),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 23: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c24
      allocate(ixnode(MAXXYZ),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 24: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c25
      allocate(iynode(MAXXYZ),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 25: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c26
      allocate(iznode(MAXXYZ),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 26: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if

c27
      allocate(nisb(MAXSB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 27: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c28
      allocate(icnode(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 28: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c29
      allocate(ixsbtosr(8*MAXSB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 29: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c30
      allocate(iysbtosr(8*MAXSB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 30: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c31
      allocate(izsbtosr(8*MAXSB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 31: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c32
      allocate(atnode(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 32: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c33
      allocate(softdat(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 33: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c
 100  format(/,' Category  number ',i2,' = ',f12.3,/,
     +         '           global prob value = ',f8.4,/,
     +         '           number of structures = ',i3,/,
     +         '           nugget effect        = ',f8.4)
 101  format(/,' Threshold number ',i2,' = ',f12.3,/,
     +         '           global prob value = ',f8.4,/,
     +         '           number of structures = ',i3,/,
     +         '           nugget effect        = ',f8.4)
 102  format(  '           type of structure ',i3,' = ',i3,/,
     +         '           aa parameter         = ',f12.4,/,
     +         '           cc parameter         = ',f12.4)
 103  format(  '           ang1, ang2, ang3     = ',3f6.2,/,
     +         '           anis1, anis2         = ',2f12.4)
c
c Perform some quick error checking:
c
      if(nx.gt.MAXX) stop 'nx is too big - modify .inc file'
      if(ny.gt.MAXY) stop 'ny is too big - modify .inc file'
      if(nz.gt.MAXZ) stop 'nz is too big - modify .inc file'
c
c Check to make sure the data file exists, then either read in the
c data or write a warning:
c
      title = 'SISIM SIMULATIONS:                      '//
     +        '                                        '
      nd = 0
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,113) datafl
 113        format('WARNING data file ',a40,' does not exist!',/,
     +             ' Hope your intention was to create an',
     +             ' unconditional simulation.')
      else
c
c The data file exists so open the file and read in the header
c information.
c
            write(*,*) 'Reading input data'
            av = 0.0
            ss = 0.0
            open(lin,file=datafl,status='OLD')
            read(lin,'(a60)',err=99) title(21:80)
            read(lin,*,err=99)       nvari
            do i=1,nvari
                  read(lin,'()',err=99)
            end do
c
c Read all the data until the end of the file:
c
 5          read(lin,*,end=6,err=99) (var(j),j=1,nvari)
            vrt = var(ivrl)
            if(vrt.lt.tmin.or.vrt.ge.tmax) go to 5
            nd  = nd + 1
            x(nd) = xmn
            y(nd) = ymn
            z(nd) = zmn
            if(ixl.gt.0) x(nd) = var(ixl)
            if(iyl.gt.0) y(nd) = var(iyl)
            if(izl.gt.0) z(nd) = var(izl)
            av = av + vrt
            ss = ss + vrt*vrt
c
c The indicator data are constructed knowing the thresholds and the
c data value.
c
            if(ivtype.eq.0) then
                  do ic=1,ncut
                        vr(nd,ic) = 0.0
                        if(int(vrt+0.5).eq.int(thres(ic)+0.5)) 
     +                  vr(nd,ic) = 1.0
                  end do
            else
                  do ic=1,ncut
                        vr(nd,ic) = 1.0
                        if(vrt.gt.thres(ic)) vr(nd,ic) = 0.0
                  end do
            end if
            vr(nd,MXCUT) = vrt
            go to 5
 6          close(lin)
c     
c Compute the averages and variances as an error check for the user:
c
            xd = max(real(nd),1.0)
            av = av / xd
            ss =(ss / xd ) - av * av
            write(*,120)    ivrl,nd,av,ss
#ifdef DEBUG
            write(ldbg,120) ivrl,nd,av,ss
#endif
 120        format(/,'  Data for SISIM: Variable number ',i2,
     +             /,'  Number of acceptable data  = ',i8,
     +             /,'  Equal Weighted Average     = ',f12.4,
     +             /,'  Equal Weighted Variance    = ',f12.4,/)
c
c Check to make sure that the grid is compatible with the data:
c
            if(ixl.le.0.and.nx.gt.1) then
               write(*,*) 'ERROR there is no X coordinate in data file'
               write(*,*) '      nx must be set to 1'
               stop
            end if
            if(iyl.le.0.and.ny.gt.1) then
               write(*,*) 'ERROR there is no Y coordinate in data file'
               write(*,*) '      ny must be set to 1'
               stop
            end if
            if(izl.le.0.and.nz.gt.1) then
               write(*,*) 'ERROR there is no Z coordinate in data file'
               write(*,*) '      nz must be set to 1'
               stop
            end if
      endif
c
c Now, if required, read in the tabulated values for details of the dist
c
      gcut(:)=0.0
      gcdf(:)=0.0
      if(ltail.eq.3.or.middle.eq.3.or.utail.eq.3) then
            ng = 0
            inquire(file=tabfl,exist=testfl)
            if(.not.testfl) stop 'ERROR tabfl does not exist'
            open(lin,file=tabfl,status='OLD')
            read(lin,*,err=97)
            read(lin,*,err=97) nvari
            do i=1,nvari
                  read(lin,*,err=97)
            end do
            tcdf = 0.0
            ng   = 0
 21         read(lin,*,end=22,err=97) (var(j),j=1,nvari)
            if(var(itabvr).lt.tmin.or.var(itabvr).ge.tmax) go to 21
            ng = ng + 1
            gcut(ng) = var(itabvr)
            gcdf(ng) = 1.0
            if(itabwt.gt.0) gcdf(ng) = var(itabwt)
            tcdf = tcdf + gcdf(ng)
            go to 21
 22         close(lin)
c
c Sort in ascending order and keep track of where the tabulated values
c switch classes:
c
            if(tcdf.le.0.0) then
                  write(*,*) 'ERROR: either the weights are zero or'
                  write(*,*) '       there are no tabulated data.'
                  stop
            endif
            call sortem(1,ng,gcut,1,gcdf,c,d,e,f,g,h)
c
c Set up gcdf for tabulated quantiles:
c
            oldcp = 0.0
            cp    = 0.0
            tcdf  = 1.0 / tcdf
            do i=1,ng
                  cp      = cp + gcdf(i) * tcdf
                  gcdf(i) =(cp + oldcp) * 0.5
                  oldcp   = cp
            end do
      end if
c
c Direct input of indicator data:
c
      nhd = nd
      inquire(file=softfl,exist=testfl)
      if(testfl) then
            write(*,*)
            write(*,*) 'Reading soft indicator data'
            open(lin,file=softfl,status='OLD')
            read(lin,*,err=97)
            read(lin,*,err=97) nvari
            if(ivrs(ncut).gt.nvari) then
                  write(*,*) ' ERROR: too few variables in ',softfl
                  write(*,*) '        inconsistent with parameters'
                  stop
            end if
            do i=1,nvari
                  read(lin,*,err=97)
            end do
 12         read(lin,*,end=13,err=96) (var(j),j=1,nvari)
c
c Don't keep soft data co-located with hard data:
c
            xx = xmn
            yy = ymn
            zz = zmn
            if(ixs.gt.0) xx = var(ixs)
            if(iys.gt.0) yy = var(iys)
            if(izs.gt.0) zz = var(izs)
            do i=1,nhd
                  test = abs(xx-x(i)) + abs(yy-y(i)) + abs(zz-z(i))
                  if(test.le.EPSLON) go to 12
            end do
c
c Accept this data:
c
            nd = nd + 1
            x(nd) = xx
            y(nd) = yy
            z(nd) = zz
            do j=1,ncut
                  i = ivrs(j)
                  vr(nd,j) = var(i)
                  ccdf(j)  = var(i)
            end do
c
c Draw a value for this soft distribution (in case the distribution is
c co-located with a grid node and Markov-Bayes is not used):
c

            cdfval = real(acorni(idum))
c            cdfval = 0.5
            call ordrel(ivtype,ncut,ccdf,ccdfo,nviol,aviol,xviol)
            zval = UNEST
            call beyond(ivtype,ncut,thres,ccdfo,ng,gcut,gcdf,zmin,zmax,
     +                  ltail,ltpar,middle,mpar,utail,utpar,zval,
     +                  cdfval,ierr)
            vr(nd,MXCUT) = zval
c
c If performing median IK then check for missing values:
c
            if(mik.eq.1) then
                  do ic=1,ncut
                        if(vr(nd,ic).lt.0.0) then
                              write(*,150) softfl
                              stop
                        endif
                  end do
 150              format(' Since the median IK approach is being',
     +                   ' considered no missing values are',
     +                   ' allowed',/,' Check file ',a40)
            endif
            go to 12
 13         close(lin)
            write(*,*) ' finished reading soft data'
      endif
c
c Load the right variogram as the first one if performing median IK:
c
      if(mik.eq.1) then
            icut = 1
            clos = abs(cutmik-thres(1))
            do ic=2,ncut
                  test = abs(cutmik-thres(ic))
                  if(test.lt.clos) then
                        icut = ic
                        clos = test
                  end if
            end do
            c0(1)   = c0(icut)
            nst(1)  = nst(icut)
            istart1 = 1
            istarti = 1 + (icut-1)*MAXNST
            do ist=1,nst(1)
                  index1        = istart1 + ist - 1
                  indexi        = istarti + ist - 1
                  it(index1)    = it(indexi)
                  aa(index1)    = aa(indexi)
                  aainv(index1)    = aainv(indexi)
                  cc(index1)    = cc(indexi)
                  ang1(index1)  = ang1(indexi)
                  ang2(index1)  = ang2(indexi)
                  ang3(index1)  = ang3(indexi)
                  anis1(index1) = anis1(indexi)
                  anis2(index1) = anis2(indexi)
            end do
      end if
c
c Open the output file and return:
c

#ifdef UNFORMATTED
      open(lout,file=outfl,status='UNKNOWN',form='UNFORMATTED')
#else
      open(lout,file=outfl,status='UNKNOWN')
#endif

#ifdef _OPENMP
c each thread will write in the files file.out401, file.out402, so on.
c up to MAXTHREADS
      do i=1,numThreads
            write(outfltmp,"(A40,I3)") trim(outfl(1:40)),
     + loutThreads(i)
            outfltmp = adjustl(trim(outfltmp))

#ifdef UNFORMATTED
            open(loutThreads(i),file=outfltmp,status='UNKNOWN',
     + form='UNFORMATTED')
#else
            open(loutThreads(i),file=outfltmp,status='UNKNOWN')
#endif

       end do
#endif

#ifndef UNFORMATTED
      write(lout,104) title
 104  format(a80)
      write(lout,105) 1,nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz,nsim
 105  format(i2,3(1x,i4),3(1x,g14.8),3(1x,g12.6),i4) 
      write(lout,106)
 106  format('Simulated Value')
#endif

      go to 1111
c
c Error in an Input File Somewhere:
c
 96   stop 'ERROR in soft data file!'
 97   stop 'ERROR in table look up file!'
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'



 1111 print *,'END READING PARAMETERS'
      print *,'START SIMULATION'


ccccccccccccccc sisim cccccccccccccccccccccc

c
c Set up the rotation/anisotropy matrices that are needed for the
c variogram and search:
c
      write(*,*) 'Setting up rotation matrices for variogram and search'
      do ic=1,ncut
      do is=1,nst(ic)
            ind = is + (ic-1)*MAXNST
            call setrot(ang1(ind),ang2(ind),ang3(ind),anis1(ind),
     +                  anis2(ind),ind,MAXROT,rotmat)
      end do
      end do
      isrot = MAXNST*MAXCUT + 1
      call setrot(sang1,sang2,sang3,sanis1,sanis2,isrot,MAXROT,rotmat)
c
c Set up for super block searching:
c
      if(sstrat.eq.0) then
            write(*,*) 'Setting up super block search strategy'
            do i=1,nd
                  actloc(i) = real(i)
            end do
            nsec = 0
            call setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z,
     +             actloc,tmp,nsec,sec1,sec2,sec3,MAXSBX,MAXSBY,MAXSBZ,
     +             nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,nzsup,
     +             zmnsup,zsizsup)
            call picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,
     +             isrot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr,
     +             iysbtosr,izsbtosr)
      end if
c
c Set up the covariance table and the spiral search:
c
      call ctable(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXROT,
     + MAXCUT,MAXORD,MAXXYZ,covtab,tmp,order,iznode,iynode,ixnode,
     + MAXNST,aa,c0,cc,cmax,idbg,isrot,it,ldbg,ncut,nlooku,nodmax,
     + nst,nx,ny,nz,xsiz,ysiz,zsiz,radsqd,rotmat,nctx,ncty,nctz
     + )
c
c MAIN LOOP OVER ALL THE SIMULAUTIONS:
c

#ifdef USE_MPI
      nsim =ceiling(real(nsim-mpisize+1)/real(mpisize))
#endif

c$omp parallel default(firstprivate)
c$omp& shared(x,y,z,vr,covtab,ixsbtosr,iysbtosr,izsbtosr)
#ifdef _OPENMP
      threadId = int(OMP_get_thread_num())
      if(numThreads>1)then

#ifdef USE_MPI
         idum=(threadId+1+myrank*MAXTHREADS)*ixv(1)
#else
         idum=(threadId+1)*ixv(1)
#endif
         ixv(:)=0
         ixv(1)=idum 

         do i=1,1000
            p = acorni(idum)
         end do

      end if
#else
      threadId = 0
#endif

      finLoop=ceiling(real(nsim-numThreads+1)/real(numThreads))

      print *,'Simulations per thread=',finLoop

      do isim=1,finLoop


c
c Work out a random path for this realization:
c
 
            do ind=1,nxyz
                  sim(ind)   = real(acorni(idum))
c                  sim(ind)   = real(ind)/real(nxyz)
                  order(ind) = ind
            end do

c
c The multiple grid search works with multiples of 4 (yes, that is
c somewhat arbitrary):
c
            if(mults.eq.1) then
                  do imult=1,nmult
                        nnz = max(1,nz/(imult*4))
                        nny = max(1,ny/(imult*4))
                        nnx = max(1,nx/(imult*4))
                        jz  = 1
                        jy  = 1
                        jx  = 1
                        do iz=1,nnz
                           if(nnz.gt.1) jz = iz*imult*4
                           do iy=1,nny
                              if(nny.gt.1) jy = iy*imult*4
                              do ix=1,nnx
                                 if(nnx.gt.1) jx = ix*imult*4
                                 index = jx + (jy-1)*nx + (jz-1)*nxy
                                 sim(index) = sim(index) - imult
                              end do
                           end do
                        end do
                  end do
            end if
c            call sortem(1,nxyz,sim,1,order,c,d,e,f,g,h)
            call sortem2(1,nxyz,sim,1,order)


c
c Initialize the simulation:
c
            do i=1,nxyz
                  sim(i) = UNEST
                  tmp(i) = 0.0
            end do
            write(*,*)
            write(*,*) ' Working on realization number: ',isim
c
c Assign the data to the closest grid node:
c
            TINY = 0.0001
            do id=1,nd
                  call getindx(nx,xmn,xsiz,x(id),ix,testind)
                  call getindx(ny,ymn,ysiz,y(id),iy,testind)
                  call getindx(nz,zmn,zsiz,z(id),iz,testind)
                  ind  = ix + (iy-1)*nx + (iz-1)*nxy
                  xx   = xmn + real(ix-1)*xsiz
                  yy   = ymn + real(iy-1)*ysiz
                  zz   = zmn + real(iz-1)*zsiz
                  test = abs(xx-x(id)) + abs(yy-y(id)) + abs(zz-z(id))
c
c Assign this data to the node (unless there is a closer data):
c
                  atnode(id) = .false.
                  if(sstrat.eq.1)                  atnode(id) = .true.
                  if(sstrat.eq.0.and.test.le.TINY) atnode(id) = .true.
                  if(atnode(id)) then
                        if(sim(ind).ge.0.0) then
                              id2 = int(sim(ind)+0.5)
                              test2 = abs(xx-x(id2)) + abs(yy-y(id2))
     +                                               + abs(zz-z(id2))
                              if(test.le.test2) sim(ind) = real(id)
#ifdef DEBUG
                              if(idbg.ge.2) 
     + write(ldbgThreads(threadId+1),1002) id,id2
#endif
                        else
                              sim(ind) = real(id)
                        end if
                  end if
            end do
 1002       format(' WARNING data values ',2i5,' are both assigned to ',
     +           /,'         the same node - taking the closest')
c
c Now, enter the hard data values into the "sim" array and keep the
c data number in the "tmp" array (to be reset when a hard value
c is assigned to that node):
            do i=1,nxyz
                  id = int(sim(i)+0.5)
                  if(id.gt.0) then
                        if(id.le.nhd) then
                              sim(i) = vr(id,MXCUT)
                        else
                              tmp(i) = sim(i)
                              sim(i) = UNEST
                        end if
                  end if
            end do
c
c Accumulate the number and magnitude of order relations violations:
c
            nclose = 0
            irepo  = max(1,min((nxyz/10),10000))
            ntviol = 0.0
            atviol = 0.0
            do icut=1,ncut
                  nviol(icut) =  0
                  aviol(icut) =  0.0
                  xviol(icut) = -1.0
            end do
c
c MAIN LOOP OVER ALL THE NODES:
c

            do in=1,nxyz

c                  if((int(in/irepo)*irepo).eq.in) write(*,1004) in
c 1004             format('   currently on node ',i9)

                  index = int(order(in)+0.5)
c
c Do we really need to simulate this grid node location?
c

                  if(sim(index).ne.UNEST) go to 20
                  if(imbsim.eq.0.and.tmp(index).ne.0.0) then
                        id = int(tmp(index)+0.5)
                        sim(index) = vr(id,MXCUT)
                        go to 20
                  end if
c
c Location of the node we are currently working on:
c
                  iz = int((index-1)/nxy) + 1
                  iy = int((index-(iz-1)*nxy-1)/nx) + 1
                  ix = index - (iz-1)*nxy - (iy-1)*nx
                  xx = xmn + real(ix-1)*xsiz
                  yy = ymn + real(iy-1)*ysiz
                  zz = zmn + real(iz-1)*zsiz
#ifdef DEBUG
                  if(idbg.ge.3)
     +            write(ldbgThreads(threadId+1),*) 
     + 'Working on grid index:',ix,iy,iz
#endif
c
c Now, we'll simulate the point ix,iy,iz.  First, get the close data
c and make sure that there are enough to actually simulate a value,
c we'll only keep the closest "ndmax" data, and look for previously
c simulated grid nodes:
c

                  if(sstrat.eq.0) then
                        call srchsupr(xx,yy,zz,radsqd,isrot,MAXROT,
     +                       rotmat,nsbtosr,ixsbtosr,iysbtosr,
     +                       izsbtosr,noct,nd,x,y,z,tmpdat,nisb,nxsup,
     +                       xmnsup,xsizsup,nysup,ymnsup,ysizsup,
     +                       nzsup,zmnsup,zsizsup,nclose,close,
     +                       infoct)
                        if(nclose.gt.ndmax) nclose = ndmax
c                       do i=1,nclose
c                             iii = int(close(i)+0.5)
c                             close(i) = real(actloc(iii))
c                       end do
                  endif

c                  call srchnd(ix,iy,iz,ncnodeaux,in)

                  call srchnd(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
     +        maxsec,nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz,UNEST,
     +        xmn,ymn,zmn,xsiz,ysiz,zsiz,icnode,ixnode,iynode,iznode,
     +        cnodex,cnodey,cnodez,cnodev,cnodet,tmp,sim)


c                  ncnode=ncnodeaux
#ifdef DEBUG
                  if(idbg.ge.3)
     +            write(ldbgThreads(threadId+1),*) 
     + '  there are ',nclose,' close data and ',ncnode,' close nodes'
#endif
c
c What cdf value are we looking for?
c
                  zval   = UNEST
                  cdfval = real(acorni(idum))
c                  cdfval = real(in)/real(nxyz)
c
c Use the global distribution?
c
                  if((nclose+ncnode).le.0) then
                        call beyond(ivtype,ncut,thres,cdf,ng,gcut,gcdf,
     +                              zmin,zmax,ltail,ltpar,middle,mpar,
     +                              utail,utpar,zval,cdfval,ierr)
                  else
c
c Estimate the local distribution by indicator kriging:
c
                        do ic=1,ncut
                              aclose(:)=0

              call krige(ix,iy,iz,ic,MAXCTX,MAXCTY,MAXCTZ,MAXKR1,
     + MAXROT,MAXDAT,MXCUT,MAXKR2,MAXNOD,MAXCUT,MAXXYZ,MAXNST,idbg,
     + imbsim,ivtype,ktype,ldbg,mik,nclose,ncnode,nctx,ncty,nctz, 
     + xx,yy,zz,cdf(ic),ccdf(ic),atnode,softdat,vr,vra,x,y,z,icnode,
     + iznode,iynode,ixnode,it,nst,aclose,cnodex,cnodey,cnodez,
     + cnodev,cnodet,covtab,thres,aa,aainv,c0,cc,cmax,close,beez,r,rr,
     + s,a,rotmat,MAXTHREADS,ldbgThreads,threadId)

                        end do

c
c Correct order relations:
c

                        call ordrel(ivtype,ncut,ccdf,ccdfo,nviol,aviol,
     +                              xviol)

c
c Draw from the local distribution:
c
                        call beyond(ivtype,ncut,thres,ccdfo,ng,gcut,
     +                              gcdf,zmin,zmax,ltail,ltpar,middle,
     +                              mpar,utail,utpar,zval,cdfval,ierr)

c
c Write some debugging information:
c

#ifdef DEBUG
                        if(idbg.ge.3) then
                              do ic=1,ncut
                              write(ldbgThreads(threadId+1),202) 
     + ccdf(ic),ccdfo(ic)
 202                          format('  CDF (original and fixed)',2f7.4)
                              end do
                        endif
#endif
                  endif
                  sim(index) = zval

c
c END MAIN LOOP OVER NODES:
c
 20         continue

            tmp(index) = 0.0

           end do

c
c Write this simulation to the output file:
c
            nxysim = 0
            do ic=1,ncut
                  ccdf(ic) = 0.0
            end do


            do ind=1,nxyz
c
c Calculate the cdf of the simulated values (for error checking):
c

#ifdef DEBUG
                  if(sim(ind).gt.UNEST) then
                        nxysim = nxysim + 1
                        do ic=1,ncut
                              if(ivtype.eq.0) then
                                    if(sim(ind).eq.thres(ic))
     +                                ccdf(ic)=ccdf(ic)+1.0
                              else
                                    if(sim(ind).le.thres(ic))
     +                                ccdf(ic)=ccdf(ic)+1.0
                              end if
                        end do
                  endif
#endif

#ifndef BUFFERED

#ifdef UNFORMATTED
                  write(loutThreads(threadId+1)) sim(ind)
#else
                  write(loutThreads(threadId+1),'(g14.8)') sim(ind)
#endif

#endif

            end do

c
c Report on the reproduction of the cdf and the number and magnitude
c of order relations violations:
c 
#ifdef DEBUG
c            write(*,203)    isim
            write(ldbgThreads(threadId+1),203) isim
            do icut=1,ncut
                  ccdf(icut) = ccdf(icut) / max(real(nxysim),1.0)
c                  write(*,204)    icut,cdf(icut),ccdf(icut)
                  write(ldbgThreads(threadId+1),204) 
     + icut,cdf(icut),ccdf(icut)
            end do
 203        format(/,' Finished simulation ',i2)
 204        format('     threshold ',i3,' input cdf = ',f6.4,
     +                 ' realization cdf = ',f6.4)
c            write(*,   300)
            write(ldbgThreads(threadId+1),300)
 300        format(/,' Summary of order relations: ')
            ntot = 0
            atot = 0.0
            do icut=1,ncut
               ntot = ntot + nviol(icut)
               atot = atot + aviol(icut)
               aviol(icut) = aviol(icut) / real(max(1,nviol(icut)))
c               write(*,302) icut,nviol(icut),aviol(icut),xviol(icut)
               write(ldbgThreads(threadId+1),302) 
     + icut,nviol(icut),aviol(icut),xviol(icut)
 302           format('     threshold',i2,' number = ',i6,
     +                ' average = ',f8.4,' maximum = ',f8.4)
            end do
            atot = atot / real(max(1,ntot))
            btot =(ntot / real(ncut*nxysim)) * 100.0
            write(ldbgThreads(threadId+1),303) btot,atot
c            write(*,   303) btot,atot
 303        format(/,' total of ',f18.6,'% with average of ',f8.4)
#endif

c
c END MAIN LOOP OVER SIMULATIONS:
c

#ifdef BUFFERED

#ifdef UNFORMATTED
            write(loutThreads(threadId+1)) sim
#else
            write(loutThreads(threadId+1),*) sim
#endif

#endif

      end do
c$omp end parallel

c
c Return to the main program:
c


c
c Finished:
c
      close(lout)
#ifdef _OPENMP
      do i=1,numThreads
            close(loutThreads(i))
      end do
#endif
      close(ldbg)
#ifdef _OPENMP
      do i=1,numThreads
            close(ldbgThreads(i))
      end do
#endif
      write(*,9998) VERSION
 9998 format(/' SISIM Version: ',f5.3, ' Finished'/)

c#ifdef TRACE
c      call extrae_fini()
c#endif

#ifdef USE_MPI
      call MPI_FINALIZE(ierr) 
#endif

      stop
      end
 

      subroutine ctable(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXROT,
     + MAXCUT,MAXORD,MAXXYZ,covtab,tmp,order,iznode,iynode,ixnode,
     + MAXNST,aa,c0,cc,cmax,idbg,isrot,it,ldbg,ncut,nlooku,nodmax,
     + nst,nx,ny,nz,xsiz,ysiz,zsiz,radsqd,rotmat,nctx,ncty,nctz
     + )
c-----------------------------------------------------------------------
c
c               Establish the Covariance Look up Table
c               **************************************
c
c The idea is to establish a 3-D network that contains the covariance
c value for a range of grid node offsets that should be at as large
c as twice the search radius in each direction.  The reason it has to
c be twice as large as the search radius is because we want to use it
c to compute the data covariance matrix as well as the data-block
c covariance matrix.
c
c Secondly, we want to establish a search for nearby nodes that 
c in order of closeness as defined by the variogram.
c
c
c
c INPUT VARIABLES:
c
c   xsiz,ysiz,zsiz  Definition of the grid being considered
c   MAXCTX,Y,Z      Number of blocks in covariance table
c
c   covariance table parameters
c
c
c
c OUTPUT VARIABLES:  covtab()         Covariance table
c
c EXTERNAL REFERENCES:
c
c   sqdist          Computes 3-D anisotropic squared distance
c   sortem          Sorts multiple arrays in ascending order
c   cova3           Computes the covariance according to a 3-D model
c
c
c
c-----------------------------------------------------------------------
c      use geostat

      implicit none

      integer MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXROT,
     + MAXCUT,MAXORD,MAXXYZ
      real    covtab(MAXCTX,MAXCTY,MAXCTZ,MAXCUT),tmp(MAXORD),
     + order(MAXORD)
      integer iznode(MAXXYZ),iynode(MAXXYZ),ixnode(MAXXYZ)
      integer MAXNST
      real    aa(MAXCUT*MAXNST),c0(MAXCUT),cc(MAXCUT*MAXNST)  
      real    cmax(MAXCUT) 
      integer idbg,isrot
      integer it(MAXCUT*MAXNST) 
      integer ldbg,ncut,nlooku,nodmax
      integer nst(MAXCUT)
      integer nx,ny,nz
      real    xsiz,ysiz,zsiz
      real    radsqd
      real*8  rotmat(MAXROT,3,3) 
      integer nctx,ncty,nctz

      real TINY
      parameter(TINY=1.0e-10)
c      include  'sisim.inc'
      real*8    sqdist,hsqd
      integer i,j,k,ilooku,icut,irot,ic,jc,kc
      integer il,loc,ix,iy,iz
      real    xx,yy,zz 
      real    c(1),d(1),e(1),f(1),g(1),h(1)
      real    cbb 
      real    cov
c
c Size of the look-up table:
c
      nctx = min(((MAXCTX-1)/2),(nx-1))
      ncty = min(((MAXCTY-1)/2),(ny-1))
      nctz = min(((MAXCTZ-1)/2),(nz-1))
c
c Initialize the covariance subroutine and cbb at the same time:
c
      call cova3(0.,0.,0.,0.,0.,0.,1,nst,MAXNST,c0,it,cc,aa,1,
     +           MAXROT,rotmat,cmax,cbb)
c
c Now, set up the table and keep track of the node offsets that are
c within the search radius:
c
      ilooku = max((ncut/2),1)
      nlooku = 0
      do icut=1,ncut
      irot = 1 + (icut-1)*MAXNST
      do i=-nctx,nctx
      xx = i * xsiz
      ic = nctx + 1 + i
      do j=-ncty,ncty
      yy = j * ysiz
      jc = ncty + 1 + j
      do k=-nctz,nctz
      zz = k * zsiz
      kc = nctz + 1 + k
            call cova3(0.,0.,0.,xx,yy,zz,icut,nst,MAXNST,c0,it,cc,aa,
     +                 irot,MAXROT,rotmat,cmax,cov)
            covtab(ic,jc,kc,icut) = cov
            if(icut.eq.ilooku) then
               hsqd = sqdist(0.0,0.0,0.0,xx,yy,zz,isrot,MAXROT,rotmat)
               if(real(hsqd).le.radsqd) then
                  nlooku = nlooku + 1
c
c We subtract the covariance from a large value so that the ascending
c sort subroutine will accomplish the sort we want.  Furthermore, a
c fraction of the distance is also taken off so that we search by
c anisotropic distance once we are beyond the range:
c
                  tmp(nlooku)   =-(covtab(ic,jc,kc,icut)-TINY*hsqd)
                  order(nlooku) =real((kc-1)*MAXCXY+(jc-1)*MAXCTX+ic)
               endif 
            endif
      end do
      end do
      end do
      end do

c
c Finished setting up the look-up table, now order the nodes such
c that the closest ones, according to variogram distance, are searched
c first. Note: the "loc" array is used because I didn't want to make 
c special allowance for 2 byte integers in the sorting subroutine:
c
      call sortem(1,nlooku,tmp,1,order,c,d,e,f,g,h)
      do il=1,nlooku
            loc = int(order(il))
            iz  = int((loc-1)/MAXCXY) + 1
            iy  = int((loc-(iz-1)*MAXCXY-1)/MAXCTX) + 1
            ix  = loc-(iz-1)*MAXCXY - (iy-1)*MAXCTX
            iznode(il) = iz
            iynode(il) = iy
            ixnode(il) = ix
      end do
      if(nodmax.gt.MAXNOD) then
#ifdef DEBUG
            write(ldbg,*)
            write(ldbg,*) 'The maximum number of close nodes = ',nodmax
            write(ldbg,*) 'this was reset from your specification due '
            write(ldbg,*) 'to storage limitations.'
#endif
            nodmax = MAXNOD
      endif
c
c Debugging output if requested:
c
      if(idbg.le.2) return
#ifdef DEBUG
      write(ldbg,*)
      write(ldbg,*) 'There are ',nlooku,' nearby nodes that will be '
      write(ldbg,*) 'checked until enough close data are found.'
      write(ldbg,*)
#endif
      if(idbg.lt.4) return
      do i=1,nlooku
            xx = (ixnode(i) - nctx - 1) * xsiz
            yy = (iynode(i) - ncty - 1) * ysiz
            zz = (iznode(i) - nctz - 1) * zsiz
#ifdef DEBUG
            write(ldbg,100) i,xx,yy,zz
#endif
      end do
 100  format('Point ',i6,' at ',3f18.6)
c
c All finished:
c
      return
      end



      subroutine srchnd(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
     +        maxsec,nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz,UNEST,
     +        xmn,ymn,zmn,xsiz,ysiz,zsiz,icnode,ixnode,iynode,iznode,
     +        cnodex,cnodey,cnodez,cnodev,cnodet,tmp,sim)
c-----------------------------------------------------------------------
c
c               Search for nearby Simulated Grid nodes
c               **************************************
c
c The idea is to spiral away from the node being simulated and note all
c the nearby nodes that have been simulated.
c
c
c
c INPUT VARIABLES:
c
c   ix,iy,iz        index of the point currently being simulated
c   sim(nx,ny,nz)   the simulation so far
c   nodmax          the maximum number of nodes that we want
c   nlooku          the number of nodes in the look up table
c   i[x,y,z]node    the relative indices of those nodes.
c   [x,y,z]mn       the origin of the global grid netwrok
c   [x,y,z]siz      the spacing of the grid nodes.
c
c
c
c OUTPUT VARIABLES:
c
c   ncnode          the number of close nodes
c   icnode()        the number in the look up table
c   cnode[x,y,z]()  the location of the nodes
c   cnodev()        the values at the nodes
c
c
c
c-----------------------------------------------------------------------
c      use geostat
c      include  'sisim.inc'

      implicit none

      integer ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,maxsec,
     +        nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz 
      real    UNEST,xmn,ymn,zmn,xsiz,ysiz,zsiz
      integer icnode(MAXNOD),ixnode(MAXXYZ),iynode(MAXXYZ),
     +        iznode(MAXXYZ)
      real    cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD),
     + cnodev(MAXNOD),cnodet(MAXNOD),tmp(MAXORD),sim(MXYZ)

      integer i,j,k,il,indexloc,ncsec
c
c Consider all the nearby nodes until enough have been found:
c
      ncnode = 0
      ncsec  = 0


      do il=1,nlooku

            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il))-nctx-1)
            if(i.lt.1.or.i.gt.nx) go to 1
            j = iy + (int(iynode(il))-ncty-1)
            if(j.lt.1.or.j.gt.ny) go to 1
            k = iz + (int(iznode(il))-nctz-1)
            if(k.lt.1.or.k.gt.nz) go to 1

c
c Check this potentially informed grid node:
c
            indexloc = (k-1)*nx*ny + (j-1)*nx + i


            if(sim(indexloc).gt.UNEST.or.
     + tmp(indexloc).gt.0.5) then
                  if(sim(indexloc).le.UNEST.and.
     + tmp(indexloc).gt.0.5) then
                        ncsec  = ncsec + 1
                        if(ncsec.gt.maxsec) go to 1
                  end if
                  ncnode         = ncnode + 1
                  icnode(ncnode) = il
                  cnodex(ncnode) = xmn + real(i-1)*xsiz
                  cnodey(ncnode) = ymn + real(j-1)*ysiz
                  cnodez(ncnode) = zmn + real(k-1)*zsiz
                  cnodev(ncnode) = sim(indexloc)
                  cnodet(ncnode) = tmp(indexloc)

            endif
 1          continue
      end do
c
c Return to calling program:
c
      return
      end



      subroutine krige(ix,iy,iz,icut,MAXCTX,MAXCTY,MAXCTZ,MAXKR1,
     + MAXROT,MAXDAT,MXCUT,MAXKR2,MAXNOD,MAXCUT,MAXXYZ,MAXNST,idbg,
     + imbsim,ivtype,ktype,ldbg,mik,nclose,ncnode,nctx,ncty,nctz, 
     + xx,yy,zz,gmean,cmean,atnode,softdat,vr,vra,x,y,z,icnode,
     + iznode,iynode,ixnode,it,nst,aclose,cnodex,cnodey,cnodez,
     + cnodev,cnodet,covtab,thres,aa,aainv,c0,cc,cmax,close,beez,r,
     + rr,s,a,
     + rotmat,MAXTHREADS,ldbgThreads,threadId)
c-----------------------------------------------------------------------
c
c            Builds and Solves the SK or OK Kriging System
c            *********************************************
c
c INPUT VARIABLES:
c
c   ix,iy,iz        index of the point currently being simulated
c   xx,yy,zz        location of the point currently being simulated
c   icut            cutoff number to use for either the covariance look
c                     up table or the covariance calculation
c
c
c
c OUTPUT VARIABLES:
c
c   cmean           kriged estimate
c
c
c 
c EXTERNAL REFERENCES: ksol   Gaussian elimination system solution
c
c
c NOTE: 1. the array "aclose" is used to flag those samples which exist
c          at the cutoff currently being kriged.
c
c
c-----------------------------------------------------------------------
c      use      geostat
c      include 'sisim.inc'

      implicit none

      integer ix,iy,iz,icut,MAXCTX,MAXCTY,MAXCTZ,MAXKR1,MAXROT,MAXDAT,
     +        MXCUT,MAXKR2,MAXNOD,MAXCUT,MAXXYZ,MAXNST,idbg,imbsim,
     +        ivtype,ktype,ldbg,mik,nclose,ncnode,nctx,ncty,nctz 
      real    xx,yy,zz,gmean,cmean
      logical atnode(MAXDAT),softdat(MAXKR1)
      real    vr(MAXDAT,MXCUT),vra(MAXKR1),x(MAXDAT),y(MAXDAT),z(MAXDAT)
      integer icnode(MAXNOD),iznode(MAXXYZ),iynode(MAXXYZ),
     +        ixnode(MAXXYZ),it(MAXCUT*MAXNST),nst(MAXCUT),
     +        aclose(MAXKR1)
      real    cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD),
     +        cnodev(MAXNOD),cnodet(MAXNOD),
     +        covtab(MAXCTX,MAXCTY,MAXCTZ,MAXCUT),thres(MAXCUT),
     +        aa(MAXCUT*MAXNST),aainv(MAXCUT*MAXNST),c0(MAXCUT),
     +        cc(MAXCUT*MAXNST),
     +        cmax(MAXCUT),close(MAXDAT),beez(MAXCUT)  
      real*8  r(MAXKR1),rr(MAXKR1),s(MAXKR1),a(MAXKR2),
     +        rotmat(MAXROT,3,3) 
      integer MAXTHREADS
      integer ldbgThreads(MAXTHREADS)
      integer threadId


      logical krig,somesoft,bothsoft
      integer i,j,index,irot,in,j1,iii,ind,ix1,iy1,iz1,ie,ii,jj,is,ix2,
     +        iy2,iz2,kk,ising,mclose,na,neq,nhd,jnew
      real    x1,y1,z1,x2,y2,z2,sumwt,cov
c
c Size of the kriging system:  Some of the data values may be missing
c which would correspond to a constraint interval.  Note that there
c should not be any missing values if the median approximation is being
c considered.  The variable ``krig'' is used
c to flag whether kriging is to be done or if the previous weights are
c to be used.
c

      somesoft = .false.
      krig     = .true.
      if(mik.eq.1.and.icut.gt.1) krig = .false.
      if(krig) then
            mclose = 0
            do i=1,nclose
                  index     =  int(close(i))
                  if(.not.atnode(index).and.vr(index,icut).ge.0.0) then
                        mclose = mclose + 1
                        aclose(mclose) = index
                  endif
            end do
            na  = mclose + ncnode
            neq = na + ktype
      endif
c
c There are no data yet:
c
      irot   = 1 + (icut-1)*MAXNST
c
c Set up kriging matrices:
c
      in = 0
      j1 = 0
      do 1 j=1,na
            softdat(j) = .false.
c
c Sort out the actual location of point "j"
c
            if(j.le.mclose) then
                  index  = aclose(j)
                  vra(j) = vr(index,icut)
                  x1     = x(index)
                  y1     = y(index)
                  z1     = z(index)
                  if(index.gt.nhd) softdat(j) = .true.
            else
c
c It is a previously simulated node (keep index for table look-up):
c
                  index  = j-mclose
                  x1     = cnodex(index)
                  y1     = cnodey(index)
                  z1     = cnodez(index)
c
c Is this node informed by a hard datum or a soft datum?
c
                  if(cnodet(index).le.0.5) then
                        if(ivtype.eq.0) then
                           vra(j) = 0.0
                           if(int(cnodev(index)+0.5).eq.
     +                        int(thres(icut)+0.5)) vra(j) = 1.0
                        else
                           vra(j) = 1.0
                           if(cnodev(index).gt.thres(icut)) vra(j) = 0.0
                        end if
                  else
                        iii    = int(cnodet(index)+0.5)
                        vra(j) = vr(iii,icut)
                        softdat(j) = .true.
                  end if
                  ind    = icnode(index)
                  ix1    = ix + (int(ixnode(ind))-nctx-1)
                  iy1    = iy + (int(iynode(ind))-ncty-1)
                  iz1    = iz + (int(iznode(ind))-nctz-1)
            endif
c
c Only set up the matrix and the RHS if kriging:
c
            if(krig) then
               do 2 i=1,j
c
c Sort out the actual location of point "i"
c
                  if(i.le.mclose) then
                        index  = aclose(i)
                        x2     = x(index)
                        y2     = y(index)
                        z2     = z(index)
                  else
c
c It is a previously simulated node (keep index for table look-up):
c
                        index  = i-mclose
                        x2     = cnodex(index)
                        y2     = cnodey(index)
                        z2     = cnodez(index)
                        ind    = icnode(index)
                        ix2    = ix + (int(ixnode(ind))-nctx-1)
                        iy2    = iy + (int(iynode(ind))-ncty-1)
                        iz2    = iz + (int(iznode(ind))-nctz-1)
                  endif
c
c Now, get the covariance value:
c
                  in = in + 1
c
c Decide whether or not to use the covariance look-up table:
c
                  if(j.le.mclose.or.i.le.mclose) then
                        call cova3(x1,y1,z1,x2,y2,z2,icut,nst,MAXNST,
     +                       c0,it,cc,aa,irot,MAXROT,rotmat,cmax,cov)
c                        call covaopt3(x1,y1,z1,x2,y2,z2,icut,nst,MAXNST,
c     +                    c0,it,cc,aa,aainv,irot,MAXROT,rotmat,cmax,cov)
                        a(in) = dble(cov)
                  else
c
c Try to use the covariance look-up (if the distance is in range):
c
                        ii = nctx + 1 + (ix1 - ix2)
                        jj = ncty + 1 + (iy1 - iy2)
                        kk = nctz + 1 + (iz1 - iz2)
                        if(ii.lt.1.or.ii.gt.MAXCTX.or.
     +                     jj.lt.1.or.jj.gt.MAXCTY.or.
     +                     kk.lt.1.or.kk.gt.MAXCTZ) then
                              call cova3(x1,y1,z1,x2,y2,z2,icut,nst,
     +                                   MAXNST,c0,it,cc,aa,irot,MAXROT,
     +                                   rotmat,cmax,cov)
c                          call covaopt3(x1,y1,z1,x2,y2,z2,icut,nst,
c     +                            MAXNST,c0,it,cc,aa,aainv,irot,MAXROT,
c     +                                   rotmat,cmax,cov)
                              a(in) = dble(cov)
                        else
                              a(in) = dble(covtab(ii,jj,kk,icut))
                        endif
                  endif
 2          continue
c
c Get the RHS value (possibly with covariance look-up table):
c
            if(j.le.mclose) then
                  call cova3(xx,yy,zz,x1,y1,z1,icut,nst,MAXNST,
     +                       c0,it,cc,aa,irot,MAXROT,rotmat,cmax,cov)
c               call covaopt3(xx,yy,zz,x1,y1,z1,icut,nst,MAXNST,
c     +                  c0,it,cc,aa,aainv,irot,MAXROT,rotmat,cmax,cov)
                  r(j) = dble(cov)
            else
c
c Try to use the covariance look-up (if the distance is in range):
c
                  ii = nctx + 1 + (ix - ix1)
                  jj = ncty + 1 + (iy - iy1)
                  kk = nctz + 1 + (iz - iz1)
                  if(ii.lt.1.or.ii.gt.MAXCTX.or.
     +               jj.lt.1.or.jj.gt.MAXCTY.or.
     +               kk.lt.1.or.kk.gt.MAXCTZ) then
                        call cova3(xx,yy,zz,x1,y1,z1,icut,nst,MAXNST,
     +                       c0,it,cc,aa,irot,MAXROT,rotmat,cmax,cov)
c                     call covaopt3(xx,yy,zz,x1,y1,z1,icut,nst,MAXNST,
c     +                    c0,it,cc,aa,aainv,irot,MAXROT,rotmat,cmax,cov)
                        r(j) = dble(cov)
                  else
                        r(j) = dble(covtab(ii,jj,kk,icut))
                  endif
            endif
            rr(j) = r(j)
c
c End ``if'' block (true if kriging)
c
         endif
c
c End loop over all of the nearby data
c
      if(softdat(j)) somesoft = .true.
 1    continue

c
c If we are doing Markov-Bayes are there are soft data we need to
c correct some of the covariance values in the kriging matrix:
c
      if(imbsim.eq.1.and.somesoft) then
            in = 0
            do j=1,na
                  do i=1,j
                        in = in + 1
                        bothsoft = .false.
                        if(softdat(j).and.softdat(i)) bothsoft = .true.
c
c Correct for soft-soft covariance or soft-hard covariance:
c
                        if(bothsoft) then
                              a(in) = a(in)*dble(beez(icut))
                              if(i.ne.j) a(in) = a(in)*dble(beez(icut))
                        else
                              if(softdat(j).or.softdat(i))
     +                        a(in) = a(in)*dble(beez(icut))
                        end if
                  end do
c
c Correct the right hand side for soft-hard covariance:
c
                  if(softdat(j)) then
                        r(j)  = r(j)*dble(beez(icut))
                        rr(j) = r(j)
                  end if
            end do
      end if
c
c Addition of OK constraint:
c
      if(krig.and.ktype.eq.1) then
            do i=1,na
                  in    = in + 1
                  a(in) = 1.0
            end do
            in      = in + 1
            a(in)   = 0.0
            r(neq)  = 1.0
            rr(neq) = 1.0
      endif
c
c Write out the kriging Matrix if Seriously Debugging:
c
#ifdef DEBUG
      if(krig.and.idbg.ge.4) then
            write(ldbgThreads(threadId+1),101) ix,iy,iz
            is = 1
            do i=1,neq
                  ie = is + i - 1
                  write(ldbgThreads(threadId+1),102) 
     + i,r(i),(a(j),j=is,ie)
                  is = is + i
            end do
 101        format(/,'Kriging Matrices for Node: ',3i4,
     +               ' RHS first')
 102        format('    r(',i2,') =',f7.4,'  a= ',9(10f7.4))
      endif
#endif
c
c Solve the Kriging System:
c
      if(krig) then
            if(neq.eq.1.and.ktype.eq.0) then
                  s(1)  = r(1) / a(1)
                  ising = 0
            else
c                  call ksol(1,neq,1,a,r,s,ising)
                  call ksolO1(1,neq,1,a,r,s,ising)
            endif
      endif
c
c Write a warning if the matrix is singular:
c
      if(ising.ne.0) then
#ifdef DEBUG
            if(idbg.ge.1) then
                  write(ldbgThreads(threadId+1),*) 
     + 'WARNING SISIM: singular matrix'
                  write(ldbgThreads(threadId+1),*) 
     + '              for block',ix,iy,iz
            endif
#endif
            cmean  = 0.0
            return
      endif
c
c Write out the kriging Matrix if Seriously Debugging:
c
#ifdef DEBUG
      if(krig.and.idbg.ge.4) then
            do i=1,na
                  write(ldbgThreads(threadId+1),140) i,s(i)
            end do
 140        format(' Kriging weight for data: ',i4,' = ',f8.4)
      endif
#endif
c
c Compute the estimate, the sum of weights, correct for SK, and return:
c
      cmean = 0.0
      sumwt = 0.0
      do i=1,na
            cmean = cmean + real(s(i)) * vra(i)
            sumwt = sumwt + real(s(i))
      end do
      if(ktype.eq.0) cmean = cmean + (1.0-sumwt)*gmean
      return
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
      open(lun,file='sisim.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for SISIM',/,
     +       '                  ********************',/,/,
     +       'START OF PARAMETERS:')

      write(lun,11)
 11   format('1                             ',
     +       '-1=continuous(cdf), 0=categorical(pdf)')
      write(lun,12)
 12   format('5                             ',
     +       '-number thresholds/categories')
      write(lun,13)
 13   format('0.5   1.0   2.5   5.0   10.0  ',
     +       '-   thresholds / categories')
      write(lun,14)
 14   format('0.12  0.29  0.50  0.74  0.88  ',
     +       '-   global cdf / pdf')
      write(lun,15)
 15   format('../data/cluster.dat           ',
     +       '-file with data')
      write(lun,16)
 16   format('1   2   0   3                 ',
     +       '-   columns for X,Y,Z, and variable')
      write(lun,17)
 17   format('direct.ik                     ',
     +       '-file with soft indicator input')
      write(lun,18)
 18   format('1   2   0   3 4 5 6 7         ',
     +       '-   columns for X,Y,Z, and indicators')
      write(lun,19)
 19   format('0                             ',
     +       '-   Markov-Bayes simulation (0=no,1=yes)')
      write(lun,20)
 20   format('0.61  0.54  0.56  0.53  0.29  ',
     +       '-      calibration B(z) values')
      write(lun,21)
 21   format('-1.0e21    1.0e21             ',
     +       '-trimming limits')
      write(lun,22)
 22   format('0.0   30.0                    ',
     +       '-minimum and maximum data value')
      write(lun,23)
 23   format('1      0.0                    ',
     +       '-   lower tail option and parameter')
      write(lun,24)
 24   format('1      1.0                    ',
     +       '-   middle     option and parameter')
      write(lun,25)
 25   format('1     30.0                    ',
     +       '-   upper tail option and parameter')
      write(lun,26)
 26   format('cluster.dat                   ',
     +       '-   file with tabulated values')
      write(lun,27)
 27   format('3   0                         ',
     +       '-      columns for variable, weight')
      write(lun,28)
 28   format('0                             ',
     +       '-debugging level: 0,1,2,3')
      write(lun,29)
 29   format('sisim.dbg                     ',
     +       '-file for debugging output')
      write(lun,30)
 30   format('sisim.out                     ',
     +       '-file for simulation output')
      write(lun,31)
 31   format('1                             ',
     +       '-number of realizations')
      write(lun,32)
 32   format('50   0.5    1.0               ',
     +       '-nx,xmn,xsiz')
      write(lun,33)
 33   format('50   0.5    1.0               ',
     +       '-ny,ymn,ysiz')
      write(lun,34)
 34   format('1    1.0   10.0               ',
     +       '-nz,zmn,zsiz')
      write(lun,35)
 35   format('69069                         ',
     +       '-random number seed')
      write(lun,36)
 36   format('12                            ',
     +       '-maximum original data  for each kriging')
      write(lun,37)
 37   format('12                            ',
     +       '-maximum previous nodes for each kriging')
      write(lun,38)
 38   format('1                             ',
     +       '-maximum soft indicator nodes for kriging')
      write(lun,39)
 39   format('1                             ',
     +       '-assign data to nodes? (0=no,1=yes)')
      write(lun,40)
 40   format('0     3                       ',
     +       '-multiple grid search? (0=no,1=yes),num')
      write(lun,41)
 41   format('0                             ',
     +       '-maximum per octant    (0=not used)')
      write(lun,42)
 42   format('20.0  20.0  20.0              ',
     +       '-maximum search radii')
      write(lun,43)
 43   format(' 0.0   0.0   0.0              ',
     +       '-angles for search ellipsoid')
      write(lun,44)
 44   format('51    51    11                ',
     +       '-size of covariance lookup table')
      write(lun,47)
 47   format('0    2.5                      ',
     +       '-0=full IK, 1=median approx. (cutoff)')
      write(lun,48)
 48   format('0                             ',
     +       '-0=SK, 1=OK')
      write(lun,49)
 49   format('1    0.15                     ',
     +       '-One   nst, nugget effect')
      write(lun,50)
 50   format('1    0.85 0.0   0.0   0.0     ',
     +       '-      it,cc,ang1,ang2,ang3')
      write(lun,51)
 51   format('         10.0  10.0  10.0     ',
     +       '-      a_hmax, a_hmin, a_vert')
      write(lun,52)
 52   format('1    0.10                     ',
     +       '-Two   nst, nugget effect')
      write(lun,53)
 53   format('1    0.90 0.0   0.0   0.0     ',
     +       '-      it,cc,ang1,ang2,ang3')
      write(lun,54)
 54   format('         10.0  10.0  10.0     ',
     +       '-      a_hmax, a_hmin, a_vert')
      write(lun,55)
 55   format('1    0.10                     ',
     +       '-Three nst, nugget effect')
      write(lun,56)
 56   format('1    0.90 0.0   0.0   0.0     ',
     +       '-      it,cc,ang1,ang2,ang3')
      write(lun,57)
 57   format('         10.0  10.0  10.0     ',
     +       '-      a_hmax, a_hmin, a_vert')
      write(lun,58)
 58   format('1    0.10                     ',
     +       '-Four  nst, nugget effect')
      write(lun,59)
 59   format('1    0.90 0.0   0.0   0.0     ',
     +       '-      it,cc,ang1,ang2,ang3')
      write(lun,60)
 60   format('         10.0  10.0  10.0     ',
     +       '-      a_hmax, a_hmin, a_vert')
      write(lun,61)
 61   format('1    0.15                     ',
     +       '-Five  nst, nugget effect')
      write(lun,62)
 62   format('1    0.85 0.0   0.0   0.0     ',
     +       '-      it,cc,ang1,ang2,ang3')
      write(lun,63)
 63   format('         10.0  10.0  10.0     ',
     +       '-      a_hmax, a_hmin, a_vert')

      close(lun)
      return
      end
