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
c                Sequential Gaussian Simulation
c                ******************************
c
c The program is executed with no command line arguments.  The user
c will be prompted for the name of a parameter file.  The parameter
c file is described in the documentation (see the example sgsim.par)
c
c The output file will be a GEOEAS file containing the simulated values
c The file is ordered by x,y,z, and then simulation (i.e., x cycles
c fastest, then y, then z, then simulation number).  The values will be
c backtransformed to the original data values if a normal scores
c transform was performed.
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

      include  'sgsim.inc'

      real,allocatable      :: x(:),y(:),z(:),vr(:),wt(:),
     +          vrtr(:),vrgtr(:),close(:),sec(:),sim(:),lvm(:),
     +          tmp(:),covtab(:,:,:),cnodex(:),
     +          cnodey(:),cnodez(:),cnodev(:),vra(:),vrea(:),
     +          simbuffer(:)
      real*8,allocatable    :: r(:),rr(:),s(:),a(:)
      integer,allocatable   :: nisb(:),icnode(:),order(:)
      integer*2,allocatable :: ixnode(:),iynode(:),iznode(:),
     +          ixsbtosr(:),iysbtosr(:),izsbtosr(:)

c variables from sgsim.inc
      real cmax(1) 
      integer ndmax,mxctx,mxcty,mxctz
      integer nctx,ncty,nctz,nlooku,ncnode,nodmax
      integer KORDEI,MAXOP1,MAXINT
      parameter(KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
      integer ixv(MAXOP1)


c variables from readparams
      real      var(50)
      real*8    p,acorni,cp,oldcp,w
      real*8    acornilocal
      character transfl*512,smthfl*512,tmpfl*512,datafl*512,outfl*512,
     +          dbgfl*512,lvmfl*512,str*512
      logical   testfl,testind,trans

c#ifdef _OPENMP
      integer MAXTHREADS,BUFFSIZE
      parameter (MAXTHREADS=60)
      parameter (BUFFSIZE=4)
      integer threadId,numThreads,ithread
      integer loutThreads(MAXTHREADS)
      integer ldbgThreads(MAXTHREADS)
      character outflThreads(43,MAXTHREADS) , outfltmp*43
      character dbgflThreads(43,MAXTHREADS) , dbgfltmp*43
      integer buffcounter
c#endif

      integer finLoop
      real t1,t2,rate 
      integer c1,c2,cr,cm

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
c Input/Output units used:
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
      llvm = 4

ccccccccccccccc readparms cccccccccccccccccccccc

c
c Note VERSION number:
c
      write(*,9999) VERSION
 9999 format(/' SGSIM Version: ',f5.3/)
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
      if(str(1:1).eq.' ') str(1:20) = 'sgsim.par           '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'sgsim.par           ') then
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

      read(lin,*,err=98) ixl,iyl,izl,ivrl,iwt,isecvr
      write(*,*) ' input columns = ',ixl,iyl,izl,ivrl,iwt,isecvr

      read(lin,*,err=98) tmin,tmax
      write(*,*) ' trimming limits = ',tmin,tmax

      read(lin,*,err=98) itrans
      write(*,*) ' transformation flag = ',itrans

      read(lin,'(a512)',err=98) transfl
      call chknam(transfl,512)
      write(*,*) ' transformation file = ',transfl(1:40)

      read(lin,*,err=98) ismooth
      write(*,*) ' consider smoothed distribution (1=yes) = ',ismooth

      read(lin,'(a512)',err=98) smthfl
      call chknam(smthfl,512)
      write(*,*) ' file with smoothed distribution = ',smthfl(1:40)

      read(lin,*,err=98) isvr,iswt
      write(*,*) ' columns = ',isvr,iswt

      read(lin,*,err=98) zmin,zmax
      write(*,*) ' data limits (tails) = ',zmin,zmax

      read(lin,*,err=98) ltail,ltpar
      write(*,*) ' lower tail = ',ltail,ltpar

      read(lin,*,err=98) utail,utpar
      write(*,*) ' upper tail = ',utail,utpar

      read(lin,*,err=98) idbg
      write(*,*) ' debugging level = ',idbg

      read(lin,'(a512)',err=98) dbgfl
      call chknam(dbgfl,512)
      write(*,*) ' debugging file = ',dbgfl(1:40)
      open(ldbg,file=dbgfl,status='UNKNOWN')

#ifdef _OPENMP
      do i=1,numThreads
            write(dbgfltmp,"(A40,I3)") trim(dbgfl(1:40)),
     + ldbgThreads(i)
            dbgfltmp = adjustl(trim(dbgfltmp))
            write(*,*) ' debuging file (threads) = ',
     + dbgfltmp
      end do
#endif

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

      read(lin,'(a512)',err=98) outfl
      call chknam(outfl,512)
      write(*,*) ' output file ',outfl(1:40)

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
      write(*,*) ' number of realizations = ',nsim

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
c             p = acorni(idum)
             p = acornilocal(MAXOP1,KORDEI,MAXINT,ixv,idum)
      end do

      read(lin,*,err=98) ndmin,ndmax
      write(*,*) ' min and max data = ',ndmin,ndmax

      read(lin,*,err=98) nodmax
      write(*,*) ' maximum previous nodes = ',nodmax

      read(lin,*,err=98) sstrat
      write(*,*) ' two-part search flag = ',sstrat
      if(sstrat.eq.1) ndmax = 0

      read(lin,*,err=98) mults,nmult
      write(*,*) ' multiple grid search flag = ',mults,nmult

      read(lin,*,err=98) noct
      write(*,*) ' number of octants = ',noct

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
      
      read(lin,*,err=98) ktype
      write(*,*) ' kriging type = ',ktype
      
      trans = .true.
      if(ktype.lt.0) then
            trans = .false.
            ktype = abs(ktype)
      end if

      colocorr = 0.0
      if(ktype.eq.4) then
            backspace lin
            read(lin,*,err=98) i,colocorr
            varred = 1.0
            backspace lin
            read(lin,*,err=9990) i,xx,varred
 9990       continue
            write(*,*) ' correlation coefficient = ',colocorr
            write(*,*) ' secondary variable varred = ',varred
      end if

      read(lin,'(a512)',err=98) lvmfl
      call chknam(lvmfl,512)
      write(*,*) ' secondary model file = ',lvmfl(1:40)

      read(lin,*,err=98) icollvm
      write(*,*) ' column in secondary model file = ',icollvm

      read(lin,*,err=98) nst(1),c0(1)
      sill = c0(1)
      write(*,*) ' nst, c0 = ',nst(1),c0(1)

      if(nst(1).le.0) then
            write(*,9997) nst(1)
 9997       format(' nst must be at least 1, it has been set to ',i4,/,
     +             ' The c or a values can be set to zero')
            stop
      endif

      do i=1,nst(1)
            read(lin,*,err=98) it(i),cc(i),ang1(i),ang2(i),ang3(i)
            read(lin,*,err=98) aa(i),aa1,aa2

            aainv(i)=1.0/aa(i)

            anis1(i) = aa1 / max(aa(i),EPSLON)
            anis2(i) = aa2 / max(aa(i),EPSLON)
            sill     = sill + cc(i)
            if(it(i).eq.4) then
                  write(*,*) ' A power model is NOT allowed '
                  write(*,*) ' Choose a different model and re start '
                  stop
            endif
            write(*,*) ' it,cc,ang[1,2,3]; ',it(i),cc(i),
     +                   ang1(i),ang2(i),ang3(i)
            write(*,*) ' a1 a2 a3: ',aa(i),aa1,aa2
      end do
      write(*,*)
      close(lin)
c
c Find the needed parameters:
c
      MAXCTX = mxctx
      MAXCTY = mxcty
      MAXCTZ = mxctz
      MAXCXY = MAXCTX * MAXCTY
      MAXXYZ = MAXCTX * MAXCTY * MAXCTZ
      MAXX   = nx
      MAXY   = ny
      MAXZ   = nz
      MXYZ   = MAXX * MAXY * MAXZ
      if(MXYZ.lt.100) MXYZ = 100
      MAXNOD = nodmax
      MAXSAM = ndmax
      MAXKR1 = MAXNOD + MAXSAM + 1
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
      MAXSB = MAXSBX*MAXSBY*MAXSBZ
c
c Find MAXDAT:
c
      MAXDAT = 100
      inquire(file=datafl,exist=testfl)
      if(testfl)then
            open(lin,file=datafl,status='UNKNOWN')
            read(lin,*,err=98)
            read(lin,*,err=99) nvari
            do i=1,nvari
                  read(lin,*)
            end do
            MAXDAT = 0
 33         read(lin,*,end=66,err=98)(var(j),j=1,nvari)
            MAXDAT = MAXDAT + 1
            go to 33
 66         continue
            rewind(lin)
            close(lin)
      end if
c
      MAXTMP = 1
      inquire(file=smthfl,exist=testfl)
      if(testfl)then
            open(lin,file=smthfl,status='UNKNOWN')
            read(lin,*,err=98)
            read(lin,*,err=97) nvari
            do i=1,nvari
                  read(lin,*)
            end do
            MAXTMP = 0
 22         read(lin,*,end=55,err=97)(var(j),j=1,nvari)
            MAXTMP = MAXTMP + 1
            go to 22
 55         continue
            rewind(lin)
            close(lin)
      end if
      if(MAXTMP.gt.MAXDAT)MAXDAT = MAXTMP
c
c Allocate the needed memory:
c
      allocate(x(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 1: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(y(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 2: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(z(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 3: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(vr(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 4: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(wt(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 5: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(vrtr(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 6: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(vrgtr(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 7: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(close(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 8: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(sec(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 9: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(sim(MXYZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 10: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
      allocate(simbuffer(MXYZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 10: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if


c
      allocate(lvm(MXYZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 11: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(tmp(MAXXYZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 12: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      MAXORD = MXYZ
      if(MXYZ.lt.MAXCXY) MAXORD=MAXCXY
      allocate(order(MAXORD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 13: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(covtab(MAXCTX,MAXCTY,MAXCTZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 14: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if

c
      allocate(cnodex(MAXNOD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 15: Allocation failed',
     +                  ' due to insufficient memory.'
                        stop
            end if
c
      allocate(cnodey(MAXNOD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 16: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(cnodez(MAXNOD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 17: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(cnodev(MAXNOD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 18: Allocation failed',
     +                  ' due to insufficient memory.'
                        stop
            end if
c
      allocate(vra(MAXKR1),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 19: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(vrea(MAXKR1),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 20: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(r(MAXKR1),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 21: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(rr(MAXKR1),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 22: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(s(MAXKR1),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 23: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(a(MAXKR2),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 24: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(nisb(MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 25: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(icnode(MAXNOD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 26: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(ixnode(MAXXYZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 27: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(iynode(MAXXYZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 28: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(iznode(MAXXYZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 29: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(ixsbtosr(8*MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 30: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(iysbtosr(8*MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 31: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(izsbtosr(8*MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 32: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
c Warn the user if the sill is different than 1.0:
c
      if(sill.gt.(1.0+EPSLON).or.sill.lt.(1.0-EPSLON)) then
            write(*,*) 'WARNING the sill of your variogram is not 1.0!'
            write(*,*) '        the sill = ',sill
            write(*,*)
      end if
c
c Perform some quick error checking:
c
      testfl = .false.
      if(nx.gt.MAXX.or.ny.gt.MAXY.or.nz.gt.MAXZ) then
            write(*,*) 'ERROR: available grid size: ',MAXX,MAXY,MAXZ
            write(*,*) '       you have asked for : ',nx,ny,nz
            testfl = .true.
      end if
      if(ltail.ne.1.and.ltail.ne.2) then
            write(*,*) 'ERROR invalid lower tail option ',ltail
            write(*,*) '      only allow 1 or 2 - see manual '
            testfl = .true.
      endif
      if(utail.ne.1.and.utail.ne.2.and.utail.ne.4) then
            write(*,*) 'ERROR invalid upper tail option ',ltail
            write(*,*) '      only allow 1,2 or 4 - see manual '
            testfl = .true.
      endif
      if(utail.eq.4.and.utpar.lt.1.0) then
            write(*,*) 'ERROR invalid power for hyperbolic tail',utpar
            write(*,*) '      must be greater than 1.0!'
            testfl = .true.
      endif
      if(ltail.eq.2.and.ltpar.lt.0.0) then
            write(*,*) 'ERROR invalid power for power model',ltpar
            write(*,*) '      must be greater than 0.0!'
            testfl = .true.
      endif
      if(utail.eq.2.and.utpar.lt.0.0) then
            write(*,*) 'ERROR invalid power for power model',utpar
            write(*,*) '      must be greater than 0.0!'
            testfl = .true.
      endif
      if(testfl) stop
c
c Check to make sure the data file exists:
c
      nd = 0
      av = 0.0
      ss = 0.0
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'WARNING data file ',datafl,' does not exist!'
            write(*,*) '   - Hope your intention was to create an ',
     +                       'unconditional simulation'
            write(*,*) '   - Resetting ndmin, ndmax, and itrans  to 0 '
            write(*,*) '   - Resetting sstrat to 1 '
            ndmin  = 0
            ndmax  = 0
            sstrat = 1
      end if
c
c Establish the reference histogram for the simulation (provided that
c we have data, and we are transforming the data):
c
      if(itrans.eq.1) then
            write(*,*) 'Setting up transformation table'
c
c Decide which file to use for establishing the transformation table:
c
            if(ismooth.eq.1) then
                  tmpfl  = smthfl
                  icolvr = isvr
                  icolwt = iswt
            else
                  tmpfl  = datafl
                  icolvr = ivrl
                  icolwt = iwt
            end if
            inquire(file=tmpfl,exist=testfl)
            if(.not.testfl) then
                  write(*,*) 'ERROR: ',tmpfl,' does not exist'
                  write(*,*) '       this file is needed! '
                  stop
            endif
c
c Open up the file with reference distribution:
c
            open(lin,file=tmpfl,status='UNKNOWN')
            read(lin,'(a40)',err=98) str(1:40)
            read(lin,*,err=99) nvari
            do i=1,nvari
                  read(lin,*,err=98)
            end do
c
c Now, read in the actual data:
c
            nt     = 0
            ntr    = 0
            twt    = 0.0
 3          read(lin,*,end=4,err=99) (var(j),j=1,nvari)
c
c Trim this data?
c
            if(var(icolvr).lt.tmin.or.var(icolvr).ge.tmax) then
                  nt = nt + 1
                  go to 3
            endif
            ntr = ntr + 1
c
c Exceeded available storage?
c
            if(icolvr.gt.nvari.or.icolwt.gt.nvari) then
                  write(*,*) ' ERROR: too few columns in ref data '
                  stop
            endif
c
c Keep this data: Assign the data value and coordinate location:
c
            vrtr(ntr) = var(icolvr)
            if(icolwt.le.0) then
                  vrgtr(ntr) = 1.0
            else
                  vrgtr(ntr) = var(icolwt)
            endif
            if(vrgtr(ntr).le.0.0) then
                  ntr = ntr - 1
                  nt  = nt  + 1
                  go to 3
            end if
            twt = twt + vrgtr(ntr)
c
c Go back for another datum:
c
            go to 3
 4          close(lin)
            if(ntr.le.1) then
                  write(*,*) 'ERROR: too few data for transformation'
                  stop
            endif
c
c Write transformation table:
c
            open(lout,file=transfl,status='UNKNOWN')
c
c Sort data by value:
c
            istart = 1
            iend   = ntr
            call sortem(istart,iend,vrtr,1,vrgtr,c,d,e,f,g,h)
c
c Compute the cumulative probabilities and write transformation table
c
            twt   = max(twt,EPSLON)
            oldcp = 0.0
            cp    = 0.0
            do j=istart,iend
                  cp =  cp + dble(vrgtr(j)/twt)
                  w  = (cp + oldcp)*0.5
                  call gauinv(w,vrg,ierr)
                  if(ierr.eq.1) vrg = UNEST
                  write(lout,201) vrtr(j),vrg
 201              format(f12.5,1x,f12.5)
                  oldcp =  cp
c
c Now, reset the weight to the normal scores value:
c
                  vrgtr(j) = vrg
            end do
            close(lout)
      end if
c
c Now, read the data if the file exists:
c
      inquire(file=datafl,exist=testfl)
      if(testfl) then
            write(*,*) 'Reading input data'
            open(lin,file=datafl,status='OLD')
            read(lin,*,err=99)
            read(lin,*,err=99) nvari
            do i=1,nvari
                  read(lin,*,err=99)
            end do
            if(ixl.gt.nvari.or.iyl.gt.nvari.or.izl.gt.nvari.or.
     +         ivrl.gt.nvari.or.isecvr.gt.nvari.or.iwt.gt.nvari) then
                  write(*,*) 'ERROR: you have asked for a column number'
                  write(*,*) '       greater than available in file'
                  stop
            end if
c
c Read all the data until the end of the file:
c
            twt = 0.0
            nd  = 0
            nt  = 0
 5          read(lin,*,end=6,err=99) (var(j),j=1,nvari)
            if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax) then
                  nt = nt + 1
                  go to 5
            end if
            nd = nd + 1
c
c Acceptable data, assign the value, X, Y, Z coordinates, and weight:
c
            vr(nd) = var(ivrl)
            if(ixl.le.0) then
                  x(nd) = xmn
            else
                  x(nd) = var(ixl)
            endif
            if(iyl.le.0) then
                  y(nd) = ymn
            else
                  y(nd) = var(iyl)
            endif
            if(izl.le.0) then
                  z(nd) = zmn
            else
                  z(nd) = var(izl)
            endif
            if(iwt.le.0) then
                  wt(nd) = 1.0
            else
                  wt(nd) = var(iwt)
            endif
            if(isecvr.le.0) then
                  sec(nd) = UNEST
            else
                  sec(nd) = var(isecvr)
            endif
c
c Normal scores transform?
c
            if(itrans.eq.1) then
                  vrr = vr(nd)
                  call locate(vrtr,ntr,1,ntr,vrr,j)
                  j   = min(max(1,j),(ntr-1))
                  vrg = powint(vrtr(j),vrtr(j+1),vrgtr(j),vrgtr(j+1),
     +                         vrr,1.0)
                  if(vrg.lt.vrgtr(1)  ) vrg = vrgtr(1)
                  if(vrg.gt.vrgtr(ntr)) vrg = vrgtr(ntr)
                  vr(nd) = vrg
            end if
            twt = twt + wt(nd)
            av  = av  + var(ivrl)*wt(nd)
            ss  = ss  + var(ivrl)*var(ivrl)*wt(nd)
            go to 5
 6          close(lin)
c
c Compute the averages and variances as an error check for the user:
c
            av = av / max(twt,EPSLON)
            ss =(ss / max(twt,EPSLON)) - av * av
#ifdef DEBUG
            write(ldbg,111) nd,nt,av,ss
            write(*,   111) nd,nt,av,ss
#endif
 111  format(/,' Data for SGSIM: Number of acceptable data  = ',i8,/,
     +         '                 Number trimmed             = ',i8,/,
     +         '                 Weighted Average           = ',f12.4,/,
     +         '                 Weighted Variance          = ',f12.4,/)
      endif
c
c Read secondary attribute model if necessary:
c
      if(ktype.ge.2) then
            write(*,*) 'Reading secondary attribute file'
            inquire(file=lvmfl,exist=testfl)
            if(.not.testfl) then
                  write(*,104) lvmfl
 104              format('WARNING secondary attribute file ',a40,
     +             ' does not exist!')
                  stop
            end if
            open(llvm,file=lvmfl,status='OLD')
            read(llvm,*,err=97)
            read(llvm,*,err=97) nvaril
            do i=1,nvaril
                  read(llvm,*,err=97)
            end do
            index = 0
             
            av = 0.0
            ss = 0.0
            ns = 0
            do iz=1,nz
                  do iy=1,ny
                        do ix=1,nx
                           index = index + 1
                           read(llvm,*,err=97) (var(j),j=1,nvaril)
                           vrr = var(icollvm)
                           lvm(index) = vrr
                           sim(index) = real(index)
c
c Do we to transform the secondary variable for a local mean?
c
                           if(trans.and.ktype.eq.2.and.itrans.eq.1) then
                                 if(vrr.le.tmin.or.vrr.ge.tmax) then
                                       lvm(index) = -1.0e21
                                 else   
                                 call locate(vrtr,ntr,1,ntr,vrr,j)
                                 j   =min(max(1,j),(ntr-1))
                                 vrg =powint(vrtr(j),vrtr(j+1),vrgtr(j),
     +                                       vrgtr(j+1),vrr,1.0)
                                 if(vrg.lt.vrgtr(1)  ) vrg = vrgtr(1)
                                 if(vrg.gt.vrgtr(ntr)) vrg = vrgtr(nd)
                                 lvm(index) = vrg
                                 end if
                           end if
                           if(vrr.ge.tmin.or.vrr.le.tmax) then
                                 av = av + vrr
                                 ss = ss + vrr*vrr
                                 ns = ns + 1
                           end if   
                        end do
                  end do
            end do
            ns = max(ns,1)
            av = av / real(ns)
            ss =(ss / real(ns)) - av * av
#ifdef DEBUG
            write(ldbg,112) ns,av,ss
            write(*,   112) ns,av,ss
#endif
 112  format(/,' Secondary Data: Number of data             = ',i8,/,
     +         '                 Equal Weighted Average     = ',f12.4,/,
     +         '                 Equal Weighted Variance    = ',f12.4,/)
c
c Do we need to work with data residuals? (Locally Varying Mean)
c
            if(ktype.eq.2) then
                  do i=1,nd
                        call getindx(nx,xmn,xsiz,x(i),ix,testind)
                        call getindx(ny,ymn,ysiz,y(i),iy,testind)
                        call getindx(nz,zmn,zsiz,z(i),iz,testind)
                        index = ix + (iy-1)*nx + (iz-1)*nxy
                        sec(i) = lvm(index)
c
c Calculation of residual moved to krige subroutine: vr(i)=vr(i)-sec(i)
c
                  end do
            end if
c
c Do we need to get an external drift attribute for the data?
c
            if(ktype.eq.3) then
                  do i=1,nd
                        if(sec(i).eq.UNEST) then
                              call getindx(nx,xmn,xsiz,x(i),ix,testind)
                              call getindx(ny,ymn,ysiz,y(i),iy,testind)
                              call getindx(nz,zmn,zsiz,z(i),iz,testind)
                              index = ix + (iy-1)*nx + (iz-1)*nxy
                              sec(i) = lvm(index)
                        end if
                  end do
            end if
c
c Transform the secondary attribute to normal scores?
c
            if(trans.and.ktype.eq.4) then
#ifdef DEBUG
                  write(ldbg,113) varred
#endif
 113              format(/,' Transforming Secondary Data with',
     +                     ' variance reduction of ',f12.4,/)
                  write(*,*) 'Transforming secondary variable'
                  write(*,*)
                  call sortem(1,nxyz,lvm,1,sim,c,d,e,f,g,h)
                  oldcp = 0.0
                  cp    = 0.0
                  do i=1,nxyz
                        if(lvm(i).gt.tmin.and.lvm(i).le.tmax) then
                              cp =  cp + dble(1.0/real(ns))
                              w  = (cp + oldcp)/2.0
                              call gauinv(w,lvm(i),ierr)
                              lvm(i) = lvm(i) * varred
                              oldcp  =  cp
                        end if      
                  end do
                  call sortem(1,nxyz,sim,1,lvm,c,d,e,f,g,h)
            end if
      end if
c
c Open the output file:
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
      write(lout,210)
 210  format('SGSIM Realizations')
      write(lout,211) 1,nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz,nsim
 211  format(i2,3(1x,i4),3(1x,g14.8),3(1x,g12.6),i4) 
      write(lout,212)
 212  format('value')
#endif
      
      go to 1111 
c
c Error in an Input File Somewhere:
c
 97   stop 'ERROR in secondary data file!'
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'



 1111 print *,'END READING PARAMETERS'
      print *,'START SIMULATION'


ccccccccccccccc sgsim cccccccccccccccccccccc

c
c Set up the rotation/anisotropy matrices that are needed for the
c variogram and search.
c
      write(*,*) 'Setting up rotation matrices for variogram and search'
      do is=1,nst(1)
            call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is),
     +                  is,MAXROT,rotmat)
      end do
      isrot = MAXNST + 1
      call setrot(sang1,sang2,sang3,sanis1,sanis2,isrot,MAXROT,rotmat)
c
c Set up the super block search:
c
      if(sstrat.eq.0) then
            write(*,*) 'Setting up super block search strategy'
            nsec = 1
            call setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z,
     +                   vr,wt,nsec,sec,sec2,sec3,MAXSBX,MAXSBY,MAXSBZ,
     +                   nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,
     +                   nzsup,zmnsup,zsizsup)
            call picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,
     +                   isrot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr,
     +                   iysbtosr,izsbtosr)
      end if
c
c Set up the covariance table and the spiral search:
c
c      call ctable(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ)

      call ctable(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXROT,
     + MAXORD,MAXXYZ,covtab,tmp,order,iznode,iynode,ixnode,
     + MAXNST,aa,c0,cc,cmax,idbg,isrot,it,ldbg,nlooku,nodmax,
     + nst,nx,ny,nz,xsiz,ysiz,zsiz,radsqd,rotmat,nctx,ncty,nctz,
     + cbb)


c
c MAIN LOOP OVER ALL THE SIMULAUTIONS:
c

      buffcounter = 1

      call system_clock(count_rate=cr)
      call system_clock(count_max=cm)
      rate = REAL(cr)

      call system_clock(c1)

#ifdef USE_MPI
      nsim =ceiling(real(nsim-mpisize+1)/real(mpisize))
#endif

c$omp parallel default(firstprivate) shared(x,y,z,vr)
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
c            p = acorni(idum)
            p = acornilocal(MAXOP1,KORDEI,MAXINT,ixv,idum)
         end do

      end if

#else
      threadId = 0
#endif

      finLoop=ceiling(real(nsim-numThreads+1)/real(numThreads))

      write(*,*) 'Simulations per thread=',finLoop

      do isim=1,finLoop

c
c Read in the secondary data distribution for this realization:
c
            if(isim.gt.1.and.ktype.eq.4) then
                  write(*,*)
                  write(*,*) ' Reading next secondary model'
                  index = 0
                  do iz=1,nz
                        do iy=1,ny
                              do ix=1,nx
                                 index = index + 1
                                 read(llvm,*,end=977)(var(j),j=1,nvaril)
                                 lvm(index) = var(icollvm)
                                 sim(index) = real(index)
                              end do
                        end do
                  end do
                  write(*,*) ' Building CDF from  secondary model'
                  call sortem(1,nxyz,lvm,1,sim,c,d,e,f,g,h)
                  oldcp = 0.0
                  cp    = 0.0
                  do i=1,nxyz
                        cp =  cp + dble(1.0/real(nxyz))
                        w  = (cp + oldcp)/2.0
                        call gauinv(w,lvm(i),ierr)
                        lvm(i) = lvm(i) * varred
                        oldcp  =  cp
                  end do
                  write(*,*) ' Restoring order of secondary model'
                  call sortem(1,nxyz,sim,1,lvm,c,d,e,f,g,h)
 977              continue
                  write(*,*)
            end if
c
c Work out a random path for this realization:
c
            do ind=1,nxyz
c                  sim(ind)   = real(acorni(idum))
                  sim(ind)   = real(acornilocal(MAXOP1,KORDEI,MAXINT,
     +                             ixv,idum))
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
            do ind=1,nxyz
                  sim(ind) = UNEST
            end do
            write(*,*)
            write(*,*) 'Working on realization number ',isim
c
c Assign the data to the closest grid node:
c
            TINY = 0.0001
            do id=1,nd
                  call getindx(nx,xmn,xsiz,x(id),ix,testind)
                  call getindx(ny,ymn,ysiz,y(id),iy,testind)
                  call getindx(nz,zmn,zsiz,z(id),iz,testind)
                  ind = ix + (iy-1)*nx + (iz-1)*nxy
                  xx  = xmn + real(ix-1)*xsiz
                  yy  = ymn + real(iy-1)*ysiz
                  zz  = zmn + real(iz-1)*zsiz
                  test = abs(xx-x(id)) + abs(yy-y(id)) + abs(zz-z(id))
c
c Assign this data to the node (unless there is a closer data):
c
                  if(sstrat.eq.1) then
                        if(sim(ind).ge.0.0) then
                              id2 = int(sim(ind)+0.5)
                              test2 = abs(xx-x(id2)) + abs(yy-y(id2))
     +                                               + abs(zz-z(id2))
                              if(test.le.test2) sim(ind) = real(id)
#ifdef DEBUG
                              write(ldbgThreads(threadId+1),102) id,id2
#endif
                        else
                              sim(ind) = real(id)
                        end if
                  end if
c
c Assign a flag so that this node does not get simulated:
c
                  if(sstrat.eq.0.and.test.le.TINY) sim(ind)=10.0*UNEST
            end do
 102        format(' WARNING data values ',2i5,' are both assigned to ',
     +           /,'         the same node - taking the closest')
c
c Now, enter data values into the simulated grid:
c
            do ind=1,nxyz
                  id = int(sim(ind)+0.5)
                  if(id.gt.0) sim(ind) = vr(id)
            end do
            irepo = max(1,min((nxyz/10),10000))
c
c MAIN LOOP OVER ALL THE NODES:
c

            do in=1,nxyz

c                  if((int(in/irepo)*irepo).eq.in) write(*,103) in
c 103              format('   currently on node ',i9)

c
c Figure out the location of this point and make sure it has
c not been assigned a value already:
c
                  index = order(in)
                  if(sim(index).gt.(UNEST+EPSLON).or.
     +               sim(index).lt.(UNEST*2.0)) go to 50
                  iz = int((index-1)/nxy) + 1
                  iy = int((index-(iz-1)*nxy-1)/nx) + 1
                  ix = index - (iz-1)*nxy - (iy-1)*nx
                  xx = xmn + real(ix-1)*xsiz
                  yy = ymn + real(iy-1)*ysiz
                  zz = zmn + real(iz-1)*zsiz

c
c Now, we'll simulate the point ix,iy,iz.  First, get the close data
c and make sure that there are enough to actually simulate a value,
c we'll only keep the closest "ndmax" data, and look for previously
c simulated grid nodes:
c
                  if(sstrat.eq.0) then
                        close(:)=0.0
                        call srchsupr(xx,yy,zz,radsqd,isrot,MAXROT,
     +                          rotmat,nsbtosr,ixsbtosr,iysbtosr,
     +                          izsbtosr,noct,nd,x,y,z,wt,nisb,nxsup,
     +                          xmnsup,xsizsup,nysup,ymnsup,ysizsup,
     +                          nzsup,zmnsup,zsizsup,nclose,close,
     +                          infoct)
                        if(nclose.lt.ndmin) go to 50
                        if(nclose.gt.ndmax) nclose = ndmax
c                        call srchsuprO1(xx,yy,zz,radsqd,isrot,MAXROT,
c     +                          rotmat,nsbtosr,ixsbtosr,iysbtosr,
c     +                          izsbtosr,noct,nd,x,y,z,wt,nisb,nxsup,
c     +                          xmnsup,xsizsup,nysup,ymnsup,ysizsup,
c     +                          nzsup,zmnsup,zsizsup,nclose,close,
c     +                          infoct)
c                        if(nclose.lt.ndmin) go to 50
c                        if(nclose.gt.ndmax) nclose = ndmax
                  endif

c                  call srchnd(ix,iy,iz)
      call srchnd(ix,iy,iz,MAXNOD,MAXXYZ,MXYZ,ncnode,
     +        nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz,nxy,noct,
     +        UNEST,xmn,ymn,zmn,xsiz,ysiz,zsiz,icnode,
     +        ixnode,iynode,iznode,cnodex,cnodey,cnodez,cnodev,sim)

c
c Calculate the conditional mean and standard deviation.  This will be
c done with kriging if there are data, otherwise, the global mean and
c standard deviation will be used:
c
                  if(ktype.eq.2) then
                        gmean = lvm(index)
                  else
                        gmean = 0.0
                  end if
                  if((nclose+ncnode).lt.1) then
                        cmean  = gmean
                        cstdev = 1.0
                  else
c
c Perform the kriging.  Note that if there are fewer than four data
c then simple kriging is prefered so that the variance of the
c realization does not become artificially inflated:
c
                        lktype = ktype
                        if(ktype.eq.1.and.(nclose+ncnode).lt.4) lktype=0
c                        call krige(ix,iy,iz,xx,yy,zz,lktype,gmean,
c     +                             cmean,cstdev,MAXCTX,MAXCTY,MAXCTZ)


                        call krige(ix,iy,iz,lktype,MAXCTX,MAXCTY,MAXCTZ,
     +                        MAXKR1,MAXROT,MAXDAT,MAXKR2,MAXNOD,MAXXYZ,
     +                        MAXNST,MXYZ,idbg,ldbg,nclose,ncnode,nctx,
     +                        ncty,nctz,nx,nxy,xx,yy,zz,gmean,cmean,
     +                        cstdev,colocorr,EPSLON,vr,vra,vrea,sec,
     +                        x,y,z,lvm,icnode,iznode,iynode,ixnode,it,
     +                        nst,cnodex,cnodey,cnodez,cnodev,covtab,
     +                        aa,aainv,c0,cc,cmax,close,r,rr,s,a,rotmat,
     +                        cbb,MAXTHREADS,ldbgThreads,threadId)


                  endif
c
c Draw a random number and assign a value to this node:
c

c                  p = acorni(idum)
                  p = acornilocal(MAXOP1,KORDEI,MAXINT,
     +                             ixv,idum)
                  call gauinv(p,xp,ierr)
                  sim(index) = xp * cstdev + cmean
#ifdef DEBUG
                  if(idbg.ge.3) write(ldbgThreads(threadId+1),141) 
     + p,sim(index)
#endif
 141              format(' random number ',f6.4,' realization ',f7.4)
c
c Quick check for far out results:
c
                  if(abs(cmean).gt.5.0.or.abs(cstdev).gt.5.0.or.
     +               abs(sim(index)).gt.6.0) then
#ifdef DEBUG
                  write(ldbgThreads(threadId+1),1004) 
     + ix,iy,iz,cmean,cstdev,sim(index)
#endif
 1004             format('WARNING: grid node location: ',3i5,/,
     +                   '         conditional mean:   ',f12.5,/,
     +                   '         conditional stdev:  ',f12.5,/,
     +                   '         simulated value:    ',f12.5)
                  endif
c
c END MAIN LOOP OVER NODES:
c

 50               continue
            end do

c
c Do we need to reassign the data to the grid nodes?
c
            if(sstrat.eq.0) then
                  do id=1,nd
                        call getindx(nx,xmn,xsiz,x(id),ix,testind)
                        call getindx(ny,ymn,ysiz,y(id),iy,testind)
                        call getindx(nz,zmn,zsiz,z(id),iz,testind)
                        xx  = xmn + real(ix-1)*xsiz
                        yy  = ymn + real(iy-1)*ysiz
                        zz  = zmn + real(iz-1)*zsiz
                        ind = ix + (iy-1)*nx + (iz-1)*nxy
                        test=abs(xx-x(id))+abs(yy-y(id))+abs(zz-z(id))
                        if(test.le.TINY) sim(ind) = vr(id)
                  end do
            end if



c
c Back transform each value and write results:
c
            ne = 0
            av = 0.0
            ss = 0.0
            do ind=1,nxyz
                  simval = sim(ind)
                  if(simval.gt.-9.0.and.simval.lt.9.0) then
                        ne = ne + 1
                        av = av + simval
                        ss = ss + simval*simval
                  end if
                  if(itrans.eq.1.and.simval.gt.(UNEST+EPSLON)) then
                        simval = backtr(simval,ntr,vrtr,vrgtr,zmin,
     +                                  zmax,ltail,ltpar,utail,utpar)
                        if(simval.lt.zmin) simval = zmin
                        if(simval.gt.zmax) simval = zmax
                  end if

#ifdef BUFFERED
                  simbuffer(ind) = simval
#else

#ifdef UNFORMATTED
                  write(loutThreads(threadId+1)) simval
#else
                  write(loutThreads(threadId+1),'(g14.8)') simval
#endif

#endif

            end do
            av = av / max(real(ne),1.0)
            ss =(ss / max(real(ne),1.0)) - av * av
#ifdef DEBUG
            write(ldbgThreads(threadId+1),1011) isim,ne,av,ss
c            write(*,   1011) isim,ne,av,ss
#endif
 1011       format(/,' Realization ',i3,': number   = ',i8,/,
     +               '                  mean     = ',f12.4,
     +               ' (close to 0.0?)',/,
     +               '                  variance = ',f12.4,
     +               ' (close to gammabar(V,V)? approx. 1.0)',/)

c
c END MAIN LOOP OVER SIMULATIONS:
c

#ifdef BUFFERED

#ifdef UNFORMATTED
            write(loutThreads(threadId+1)) simbuffer
#else
            write(loutThreads(threadId+1),*) simbuffer
#endif

#endif

      end do
c$omp end parallel 

      call system_clock(c2)
      write(*,*) "simulation loop elapsed time : ",(c2 - c1)/rate

c
c Return to the main program:
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
c
c Finished:
c
      write(*,9998) VERSION
 9998 format(/' SGSIM Version: ',f5.3, ' Finished'/)

c#ifdef TRACE
c      call extrae_fini()
c#endif

#ifdef USE_MPI
      call MPI_FINALIZE(ierr) 
#endif

      stop
      end


      subroutine ctable(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXROT,
     + MAXORD,MAXXYZ,covtab,tmp,order,iznode,iynode,ixnode,
     + MAXNST,aa,c0,cc,cmax,idbg,isrot,it,ldbg,nlooku,nodmax,
     + nst,nx,ny,nz,xsiz,ysiz,zsiz,radsqd,rotmat,nctx,ncty,nctz,
     + cbb)

c-----------------------------------------------------------------------
c
c               Establish the Covariance Look up Table
c               **************************************
c
c The idea is to establish a 3-D network that contains the covariance
c value for a range of grid node offsets that should be at as large
c as twice the search radius in each direction.  The reason it has to
c be twice as large as the search radius is because we want to use it
c to compute the data covariance matrix as well as the data-point
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
c      use       geostat
c      parameter(TINY=1.0e-10)
c      include  'sgsim.inc'
c      real*8    hsqd,sqdist


      implicit none

      integer MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXROT,
     + MAXORD,MAXXYZ
      real    covtab(MAXCTX,MAXCTY,MAXCTZ),tmp(MAXORD)
      integer order(MAXORD)
      integer*2 iznode(MAXXYZ),iynode(MAXXYZ),ixnode(MAXXYZ)
      integer MAXNST
      real    aa(MAXNST),c0(1),cc(MAXNST)  
      real    cmax(1) 
      integer idbg,isrot
      integer it(MAXNST) 
      integer ldbg,nlooku,nodmax
      integer nst(1)
      integer nx,ny,nz
      real    xsiz,ysiz,zsiz
      real    radsqd
      real*8  rotmat(MAXROT,3,3) 
      integer nctx,ncty,nctz
      real    cbb 


      real TINY
      parameter(TINY=1.0e-10)
c      include  'sisim.inc'
      real*8    sqdist,hsqd
      integer i,j,k,ic,jc,kc
      integer il,loc,ix,iy,iz
      real    xx,yy,zz 
      real    c(1),d(1),e(1),f(1),g(1),h(1)



c
c Size of the look-up table:
c
      nctx = min(((MAXCTX-1)/2),(nx-1))
      ncty = min(((MAXCTY-1)/2),(ny-1))
      nctz = min(((MAXCTZ-1)/2),(nz-1))
c
c Debugging output:
c
#ifdef DEBUG
      write(ldbg,*)
      write(ldbg,*) 'Covariance Look up table and search for previously'
      write(ldbg,*) 'simulated grid nodes.  The maximum range in each '
      write(ldbg,*) 'coordinate direction for covariance look up is:'
      write(ldbg,*) '          X direction: ',nctx*xsiz
      write(ldbg,*) '          Y direction: ',ncty*ysiz
      write(ldbg,*) '          Z direction: ',nctz*zsiz
      write(ldbg,*) 'Node Values are not searched beyond this distance!'
      write(ldbg,*)
#endif
c
c NOTE: If dynamically allocating memory, and if there is no shortage
c       it would a good idea to go at least as far as the radius and
c       twice that far if you wanted to be sure that all covariances
c       in the left hand covariance matrix are within the table look-up.
c
c Initialize the covariance subroutine and cbb at the same time:
c
      call cova3(0.0,0.0,0.0,0.0,0.0,0.0,1,nst,MAXNST,c0,it,cc,aa,
     +           1,MAXROT,rotmat,cmax,cbb)
c
c Now, set up the table and keep track of the node offsets that are
c within the search radius:
c
      nlooku = 0
      do i=-nctx,nctx
      xx = i * xsiz
      ic = nctx + 1 + i
      do j=-ncty,ncty
      yy = j * ysiz
      jc = ncty + 1 + j
      do k=-nctz,nctz
      zz = k * zsiz
      kc = nctz + 1 + k
            call cova3(0.0,0.0,0.0,xx,yy,zz,1,nst,MAXNST,c0,it,cc,aa,
     +                 1,MAXROT,rotmat,cmax,covtab(ic,jc,kc))
            hsqd = sqdist(0.0,0.0,0.0,xx,yy,zz,isrot,
     +                          MAXROT,rotmat)
            if(real(hsqd).le.radsqd) then
                  nlooku         = nlooku + 1
c
c We want to search by closest variogram distance (and use the
c anisotropic Euclidean distance to break ties:
c
                  tmp(nlooku)   = - (covtab(ic,jc,kc)-TINY*real(hsqd))
c                  order(nlooku) = real((kc-1)*MAXCXY+(jc-1)*MAXCTX+ic)
                  order(nlooku) = (kc-1)*MAXCXY+(jc-1)*MAXCTX+ic
            endif
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
            iznode(il) = int(iz)
            iynode(il) = int(iy)
            ixnode(il) = int(ix)
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
      if(idbg.lt.2) return
#ifdef DEBUG
      write(ldbg,*)
      write(ldbg,*) 'There are ',nlooku,' nearby nodes that will be '
      write(ldbg,*) 'checked until enough close data are found.'
      write(ldbg,*)
#endif
      if(idbg.lt.14) return
      do i=1,nlooku
            xx = (ixnode(i) - nctx - 1) * xsiz
            yy = (iynode(i) - ncty - 1) * ysiz
            zz = (iznode(i) - nctz - 1) * zsiz
#ifdef DEBUG
            write(ldbg,100) i,xx,yy,zz
#endif
      end do
 100  format('Point ',i3,' at ',3f12.4)
c
c All finished:
c
      return
      end



c      subroutine srchnd(ix,iy,iz)

      subroutine srchnd(ix,iy,iz,MAXNOD,MAXXYZ,MXYZ,ncnode,
     +        nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz,nxy,noct,
     +        UNEST,
     +        xmn,ymn,zmn,xsiz,ysiz,zsiz,icnode,ixnode,iynode,iznode,
     +        cnodex,cnodey,cnodez,cnodev,sim)



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
c   sim             the realization so far
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
c      use       geostat
c      include  'sgsim.inc'


      implicit none

      integer ix,iy,iz,MAXNOD,MAXXYZ,MXYZ,ncnode,
     +        nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz,nxy,noct 
      real    UNEST,xmn,ymn,zmn,xsiz,ysiz,zsiz
      integer icnode(MAXNOD)
      integer*2 ixnode(MAXXYZ),iynode(MAXXYZ),
     +        iznode(MAXXYZ)
      real    cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD),
     + cnodev(MAXNOD),sim(MXYZ)

      integer i,j,k,il,ind,idx,idy,idz,iq

      integer   ninoct(8)
c
c Consider all the nearby nodes until enough have been found:
c
      ncnode = 0
      if(noct.gt.0) then
            do i=1,8
                  ninoct(i) = 0
            end do
      end if
      do 2 il=2,nlooku
            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il))-nctx-1)
            j = iy + (int(iynode(il))-ncty-1)
            k = iz + (int(iznode(il))-nctz-1)
            if(i.lt. 1.or.j.lt. 1.or.k.lt. 1) go to 2
            if(i.gt.nx.or.j.gt.ny.or.k.gt.nz) go to 2
            ind = i + (j-1)*nx + (k-1)*nxy
            if(sim(ind).gt.UNEST) then
c
c Check the number of data already taken from this octant:
c
                  if(noct.gt.0) then
                        idx = ix - i
                        idy = iy - j
                        idz = iz - k
                        if(idz.gt.0) then
                              iq = 4
                              if(idx.le.0 .and. idy.gt.0) iq = 1
                              if(idx.gt.0 .and. idy.ge.0) iq = 2
                              if(idx.lt.0 .and. idy.le.0) iq = 3
                        else
                              iq = 8
                              if(idx.le.0 .and. idy.gt.0) iq = 5
                              if(idx.gt.0 .and. idy.ge.0) iq = 6
                              if(idx.lt.0 .and. idy.le.0) iq = 7
                        end if
                        ninoct(iq) = ninoct(iq) + 1
                        if(ninoct(iq).gt.noct) go to 2
                  end if
                  ncnode = ncnode + 1
                  icnode(ncnode) = il
                  cnodex(ncnode) = xmn + real(i-1)*xsiz
                  cnodey(ncnode) = ymn + real(j-1)*ysiz
                  cnodez(ncnode) = zmn + real(k-1)*zsiz
                  cnodev(ncnode) = sim(ind)
            endif
 2    continue
c
c Return to calling program:
c
      return
      end



c      subroutine krige(ix,iy,iz,xx,yy,zz,lktype,gmean,cmean,cstdev,
c     +                  MAXCTX,MAXCTY,MAXCTZ)

      subroutine krige(ix,iy,iz,lktype,MAXCTX,MAXCTY,MAXCTZ,MAXKR1,
     + MAXROT,MAXDAT,MAXKR2,MAXNOD,MAXXYZ,MAXNST,MXYZ,idbg,
     + ldbg,nclose,ncnode,nctx,ncty,nctz,nx,nxy,
     + xx,yy,zz,gmean,cmean,cstdev,colocorr,EPSLON,vr,
     + vra,vrea,sec,x,y,z,lvm,icnode,
     + iznode,iynode,ixnode,it,nst,cnodex,cnodey,cnodez,
     + cnodev,covtab,aa,aainv,c0,cc,cmax,close,r,rr,s,a,
     + rotmat,cbb,MAXTHREADS,ldbgThreads,threadId)

      implicit none

      integer ix,iy,iz,lktype,MAXCTX,MAXCTY,MAXCTZ,MAXKR1,MAXROT,MAXDAT,
     +        MAXKR2,MAXNOD,MAXXYZ,MAXNST,MXYZ,idbg,
     +        ldbg,nclose,ncnode,nctx,ncty,nctz,nx,nxy 
      real    xx,yy,zz,gmean,cmean,cstdev,colocorr,EPSLON
      real    vr(MAXDAT),vra(MAXKR1),vrea(MAXKR1),sec(MAXDAT),
     + x(MAXDAT),y(MAXDAT),z(MAXDAT),lvm(MXYZ) 
      integer icnode(MAXNOD)
      integer*2 iznode(MAXXYZ),iynode(MAXXYZ),ixnode(MAXXYZ)
      integer it(MAXNST),nst(1)
      real    cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD),
     +        cnodev(MAXNOD),
     +        covtab(MAXCTX,MAXCTY,MAXCTZ),
     +        aa(MAXNST),aainv(MAXNST),
     +        c0(1),cc(MAXNST),
     +        cmax(1),close(MAXDAT)  
      real*8  r(MAXKR1),rr(MAXKR1),s(MAXKR1),a(MAXKR2),
     +        rotmat(MAXROT,3,3) 
      real    cbb
      integer MAXTHREADS
      integer ldbgThreads(MAXTHREADS)
      integer threadId


c-----------------------------------------------------------------------
c
c            Builds and Solves the SK or OK Kriging System
c            *********************************************
c
c INPUT VARIABLES:
c
c   ix,iy,iz        index of the point currently being simulated
c   xx,yy,zz        location of the point currently being simulated
c
c
c
c OUTPUT VARIABLES:
c
c   cmean           kriged estimate
c   cstdev          kriged standard deviation
c
c
c
c EXTERNAL REFERENCES: ksol   Gaussian elimination system solution
c
c
c
c-----------------------------------------------------------------------
c      use      geostat
c      include 'sgsim.inc'





      logical first
      integer i,j,index,in,ind,ix1,iy1,iz1,ie,ii,jj,is,ix2,
     +        iy2,iz2,kk,ising,na,neq
      real    x1,y1,z1,x2,y2,z2,sumwts,cov
      real*8  edmin,edmax,sfmin,sfmax


c
c Size of the kriging system:
c
      first = .false.
      na    = nclose + ncnode
 33   continue      
c      if(lktype.eq.0) neq = na
c      if(lktype.eq.1) neq = na + 1
c      if(lktype.eq.2) neq = na
c      if(lktype.eq.3) neq = na + 2
c      if(lktype.eq.4) neq = na + 1

      if(lktype.eq.0) then 
         neq = na
      else 
         if(lktype.eq.1) then 
            neq = na + 1
         else 
            if(lktype.eq.2) then 
               neq = na
            else 
               if(lktype.eq.3) then 
                  neq = na + 2
               else 
                  if(lktype.eq.4) neq = na + 1
               end if 
            end if 
         end if 
      end if 

      if(lktype.ge.3) then
            ind = ix + (iy-1)*nx + (iz-1)*nxy
            if(lvm(ind).le.-6.0.or.lvm(ind).ge.6.0) then
                  lktype = 0
                  go to 33
            end if      
      end if
c
c Set up kriging matrices:
c
      in=0
      do j=1,na
c
c Sort out the actual location of point "j"
c
            if(j.le.nclose) then
                  index  = int(close(j))
                  x1     = x(index)
                  y1     = y(index)
                  z1     = z(index)
                  vra(j) = vr(index)
                  vrea(j)= sec(index)
                  if(lktype.eq.2) vra(j) = vra(j) - vrea(j)
            else
c
c It is a previously simulated node (keep index for table look-up):
c
                  index  = j-nclose
                  x1     = cnodex(index)
                  y1     = cnodey(index)
                  z1     = cnodez(index)
                  vra(j) = cnodev(index)
                  ind    = icnode(index)
                  ix1    = ix + (int(ixnode(ind))-nctx-1)
                  iy1    = iy + (int(iynode(ind))-ncty-1)
                  iz1    = iz + (int(iznode(ind))-nctz-1)
                  index  = ix1 + (iy1-1)*nx + (iz1-1)*nxy
                  vrea(j)= lvm(index)
                  if(lktype.eq.2) vra(j) = vra(j) - vrea(j)
            endif
            do i=1,j
c
c Sort out the actual location of point "i"
c
                  if(i.le.nclose) then
                        index  = int(close(i))
                        x2     = x(index)
                        y2     = y(index)
                        z2     = z(index)
                  else
c
c It is a previously simulated node (keep index for table look-up):
c
                        index  = i-nclose
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
                  if(j.le.nclose.or.i.le.nclose) then
c                        call cova3(x1,y1,z1,x2,y2,z2,1,nst,MAXNST,c0,it,
c     +                             cc,aa,1,MAXROT,rotmat,cmax,cov)
                        call covaopt3(x1,y1,z1,x2,y2,z2,1,nst,MAXNST,c0,
     +                          it,cc,aa,aainv,1,MAXROT,rotmat,cmax,cov)
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
c                              call cova3(x1,y1,z1,x2,y2,z2,1,nst,MAXNST,
c     +                             c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
                              call covaopt3(x1,y1,z1,x2,y2,z2,1,nst,
     +                             MAXNST,c0,it,cc,aa,aainv,1,MAXROT,
     +                             rotmat,cmax,cov)
                        else
                              cov = covtab(ii,jj,kk)
                        endif
                        a(in) = dble(cov)
                  endif
            end do
c
c Get the RHS value (possibly with covariance look-up table):
c
            if(j.le.nclose) then
c                  call cova3(xx,yy,zz,x1,y1,z1,1,nst,MAXNST,c0,it,cc,aa,
c     +                       1,MAXROT,rotmat,cmax,cov)
                  call covaopt3(xx,yy,zz,x1,y1,z1,1,nst,MAXNST,c0,it,cc,
     +                       aa,aainv,1,MAXROT,rotmat,cmax,cov)
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
c                        call cova3(xx,yy,zz,x1,y1,z1,1,nst,MAXNST,c0,it,
c     +                             cc,aa,1,MAXROT,rotmat,cmax,cov)
                        call covaopt3(xx,yy,zz,x1,y1,z1,1,nst,MAXNST,c0,
     +                             it,cc,aa,aainv,1,MAXROT,rotmat,cmax,
     +                             cov)
                  else
                        cov = covtab(ii,jj,kk)
                  endif
                  r(j) = dble(cov)
            endif
            rr(j) = r(j)
      end do
c
c Addition of OK constraint:
c
      if(lktype.eq.1.or.lktype.eq.3) then
            do i=1,na
                  in    = in + 1
                  a(in) = 1.0
            end do
            in       = in + 1
            a(in)    = 0.0
            r(na+1)  = 1.0
            rr(na+1) = 1.0
      endif
c
c Addition of the External Drift Constraint:
c
      if(lktype.eq.3) then
            edmin =  999999.
            edmax = -999999.
            do i=1,na
                  in    = in + 1
                  a(in) = vrea(i)
                  if(a(in).lt.edmin) edmin = a(in)
                  if(a(in).gt.edmax) edmax = a(in)
            end do
            in       = in + 1
            a(in)    = 0.0
            in       = in + 1
            a(in)    = 0.0
            ind      = ix + (iy-1)*nx + (iz-1)*nxy
            r(na+2)  = dble(lvm(ind))
            rr(na+2) = r(na+2)
            if((edmax-edmin).lt.EPSLON) neq = neq - 1
      endif
c
c Addition of Collocated Cosimulation Constraint:
c
      if(lktype.eq.4) then
            sfmin =  1.0e21
            sfmax = -1.0e21
            do i=1,na
                  in    = in + 1
                  a(in) = dble(colocorr)*r(i)
                  if(a(in).lt.sfmin) sfmin = a(in)
                  if(a(in).gt.sfmax) sfmax = a(in)
            end do
            in    = in + 1
            a(in) = 1.0
            ii    = na + 1
            r(ii) = dble(colocorr)
            rr(ii)= r(ii)
c           if((sfmax-sfmin).lt.EPSLON) neq = neq - 1
      end if
c
c Write out the kriging Matrix if Seriously Debugging:
c
#ifdef DEBUG
      if(idbg.ge.3) then
            write(ldbgThreads(threadId+1),100) ix,iy,iz
            is = 1
            do i=1,neq
                  ie = is + i - 1
                  write(ldbgThreads(threadId+1),101) 
     + i,r(i),(a(j),j=is,ie)
                  is = is + i
            end do
 100        format(/,'Kriging Matrices for Node: ',3i4,' RHS first')
 101        format('    r(',i2,') =',f7.4,'  a= ',99f7.4)
      endif
#endif

c
c Solve the Kriging System:
c
      if(neq.eq.1.and.lktype.ne.3) then
            s(1)  = r(1) / a(1)
            ising = 0
      else
c            call ksol(1,neq,1,a,r,s,ising)
            call ksolO1(1,neq,1,a,r,s,ising)
      endif
c
c Write a warning if the matrix is singular:
c
      if(ising.ne.0) then
#ifdef DEBUG
            if(idbg.ge.1) then
                  write(ldbgThreads(threadId+1),*) 
     + 'WARNING SGSIM: singular matrix'
                  write(ldbgThreads(threadId+1),*) 
     + '               for node',ix,iy,iz
            endif
#endif
            cmean  = gmean
            cstdev = 1.0
            return
      endif
c
c Compute the estimate and kriging variance.  Recall that kriging type
c     0 = Simple Kriging:
c     1 = Ordinary Kriging:
c     2 = Locally Varying Mean:
c     3 = External Drift:
c     4 = Collocated Cosimulation:
c

      cmean  = 0.0
      cstdev = cbb
c      cstdev = 0.0
      sumwts = 0.0
      do i=1,na
            cmean  = cmean  + real(s(i))*vra(i)
            cstdev = cstdev - real(s(i)*rr(i))
            sumwts = sumwts + real(s(i))
      end do

      if(lktype.eq.1) then 
         cstdev = cstdev - real(s(na+1))
      else 
         if(lktype.eq.2) then 
            cmean  = cmean + gmean
         else 
            if(lktype.eq.4) then
               ind    = ix + (iy-1)*nx + (iz-1)*nxy
               cmean  = cmean  + real(s(na+1))*lvm(ind)
               cstdev = cstdev - real(s(na+1) *rr(na+1))
            end if
         end if
      end if

c
c Error message if negative variance:
c
      if(cstdev.lt.0.0) then
#ifdef DEBUG
            write(ldbgThreads(threadId+1),*) 
     + 'ERROR: Negative Variance: ',cstdev
#endif
            cstdev = 0.0
      endif
      cstdev = sqrt(max(cstdev,0.0))
c
c Write out the kriging Weights if Seriously Debugging:
c

#ifdef DEBUG
      if(idbg.ge.3) then
            do i=1,na
                  write(ldbgThreads(threadId+1),140) i,vra(i),s(i)
            end do
 140        format(' Data ',i4,' value ',f8.4,' weight ',f8.4)
            if(lktype.eq.4) write(ldbgThreads(threadId+1),141) 
     + lvm(ind),s(na+1)
 141        format(' Sec Data  value ',f8.4,' weight ',f8.4)
            write(ldbgThreads(threadId+1),142) gmean,cmean,cstdev
 142        format(' Global mean ',f8.4,' conditional ',f8.4,
     +             ' std dev ',f8.4)
      end if
#endif

c
c Finished Here:
c
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
      open(lun,file='sgsim.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for SGSIM',/,
     +       '                  ********************',/,/,
     +       'START OF PARAMETERS:')

      write(lun,11)
 11   format('../data/cluster.dat           ',
     +       '-file with data')
      write(lun,12)
 12   format('1  2  0  3  5  0              ',
     +       '-  columns for X,Y,Z,vr,wt,sec.var.')
      write(lun,13)
 13   format('-1.0       1.0e21             ',
     +       '-  trimming limits')
      write(lun,14)
 14   format('1                             ',
     +       '-transform the data (0=no, 1=yes)')
      write(lun,15)
 15   format('sgsim.trn                     ',
     +       '-  file for output trans table')
      write(lun,16)
 16   format('0                             ',
     +       '-  consider ref. dist (0=no, 1=yes)')
      write(lun,17)
 17   format('histsmth.out                  ',
     +       '-  file with ref. dist distribution')
      write(lun,18)
 18   format('1  2                          ',
     +       '-  columns for vr and wt')
      write(lun,19)
 19   format('0.0    15.0                   ',
     +       '-  zmin,zmax(tail extrapolation)')
      write(lun,20)
 20   format('1       0.0                   ',
     +       '-  lower tail option, parameter')
      write(lun,21)
 21   format('1      15.0                   ',
     +       '-  upper tail option, parameter')
      write(lun,22)
 22   format('1                             ',
     +       '-debugging level: 0,1,2,3')
      write(lun,23)
 23   format('sgsim.dbg                     ',
     +       '-file for debugging output')
      write(lun,24)
 24   format('sgsim.out                     ',
     +       '-file for simulation output')
      write(lun,25)
 25   format('1                             ',
     +       '-number of realizations to generate')
      write(lun,26)
 26   format('50    0.5    1.0              ',
     +       '-nx,xmn,xsiz')
      write(lun,27)
 27   format('50    0.5    1.0              ',
     +       '-ny,ymn,ysiz')
      write(lun,28)
 28   format('1     0.5    1.0              ',
     +       '-nz,zmn,zsiz')
      write(lun,29)
 29   format('69069                         ',
     +       '-random number seed')
      write(lun,30)
 30   format('0     8                       ',
     +       '-min and max original data for sim')
      write(lun,31)
 31   format('12                            ',
     +       '-number of simulated nodes to use')
      write(lun,32)
 32   format('1                             ',
     +       '-assign data to nodes (0=no, 1=yes)')
      write(lun,33)
 33   format('1     3                       ',
     +       '-multiple grid search (0=no, 1=yes),num')
      write(lun,34)
 34   format('0                             ',
     +       '-maximum data per octant (0=not used)')
      write(lun,35)
 35   format('10.0  10.0  10.0              ',
     +       '-maximum search radii (hmax,hmin,vert)')
      write(lun,36)
 36   format(' 0.0   0.0   0.0              ',
     +       '-angles for search ellipsoid')
      write(lun,37)
 37   format('51    51    11                ',
     +       '-size of covariance lookup table')
      write(lun,38)
 38   format('0     0.60   1.0              ',
     +       '-ktype: 0=SK,1=OK,2=LVM,3=EXDR,4=COLC')
      write(lun,39)
 39   format('../data/ydata.dat             ',
     +       '-  file with LVM, EXDR, or COLC variable')
      write(lun,40)
 40   format('4                             ',
     +       '-  column for secondary variable')
      write(lun,41)
 41   format('1    0.1                      ',
     +       '-nst, nugget effect')
      write(lun,42)
 42   format('1    0.9  0.0   0.0   0.0     ',
     +       '-it,cc,ang1,ang2,ang3')
      write(lun,43)
 43   format('         10.0  10.0  10.0     ',
     +       '-a_hmax, a_hmin, a_vert')

      close(lun)
      return
      end
