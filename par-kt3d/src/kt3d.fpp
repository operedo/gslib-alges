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
c             Kriging (SK,OK,KT) of a 3-D Rectangular Grid
c             ********************************************
c
c The program is executed with no command line arguments.  The user
c will be prompted for the name of a parameter file.  The parameter
c file is described in the documentation (see the example kt3d.par)
c and should contain the following information:
c
c
c
c AUTHOR: Clayton V. Deutsch                             DATE: 1989-1999
c OpenMP/MPI Modifications: Oscar F. Peredo              DATE: 2014-2015
c-----------------------------------------------------------------------

#ifdef _OPENMP
      use omp_lib
#endif
#ifdef USE_MPI
      include 'mpif.h'
#endif

      include  'kt3d.inc'


      integer,allocatable :: nisb(:),ixsbtosr(:),iysbtosr(:),izsbtosr(:)
      real,allocatable    :: x(:),y(:),z(:),vr(:),ve(:),dh(:),tmp(:),
     +          close(:),xa(:),ya(:),za(:),vra(:),vea(:),xdb(:),ydb(:),
     +          zdb(:),cut(:),cdf(:)
      real*8,allocatable  :: r(:),rr(:),s(:),a(:)
      real,allocatable :: buffer(:,:),buffertmp(:,:)
      
c readparm variables
      parameter(MV=100)
      real      var(MV)
      character datafl*512,jackfl*512,extfl*512,outfl*512,dbgfl*512,
     +          str*512,title*80
      logical   testfl


c kt3d variables
      real*8     cbb
      logical    first,fircon,accept
      data       fircon/.true./

      integer threadId,numThreads
      parameter (MAXTHREADS=60)
      integer ithread
      integer loutThreads(MAXTHREADS)
      character outflThreads(43,MAXTHREADS) , outfltmp*43

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

      numThreads=1
#ifdef _OPENMP
c$omp parallel
      numThreads = OMP_get_num_threads()
c$omp end parallel
#endif


ccccccccccccccc readparms cccccccccccccccccccccc

c
c FORTRAN Units:
c
      lin   = 1
      ldbg  = 3
      lout  = 4
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
      lext  = 7
      ljack = 8
c
c Note VERSION number:
c
      write(*,9999) VERSION
 9999 format(/' KT3D Version: ',f5.3/)
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
      if(str(1:1).eq.' ') str(1:20) = 'kt3d.par            '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'kt3d.par            ') then
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

      read(lin,*,err=98) idhl,ixl,iyl,izl,ivrl,iextv
      write(*,*) ' columns = ',idhl,ixl,iyl,izl,ivrl,iextv

      read(lin,*,err=98) tmin,tmax
      write(*,*) ' trimming limits = ',tmin,tmax

      read(lin,*,err=98) koption
      write(*,*) ' kriging option = ',koption

c
c This is an undocumented feature to have kt3d construct an IK-type
c distribution:
c
      iktype = 0
      if(koption.lt.0) then
            iktype  = 1
            koption = -koption
      end if
      if(iktype.eq.1) then

            read(lin,*,err=98) ncut
            write(*,*) ' number of cutoffs = ',ncut
c
c Find the needed parameter:
c
            MAXCUT = ncut
c
c Allocate the needed memory:
c21
            allocate(cut(MAXCUT),stat = test)
                  if(test.ne.0)then
                        write(*,*)'ERROR: Allocation failed due to',
     +                        ' insufficient memory.'
                        stop
                  end if
c22
            allocate(cdf(MAXCUT),stat = test)
                  if(test.ne.0)then
                        write(*,*)'ERROR: Allocation failed due to',
     +                        ' insufficient memory.'
                        stop
                  end if
c
            read(lin,*,err=98) (cut(i),i=1,ncut)
            write(*,*) ' cutoffs = ',(cut(i),i=1,ncut)

      end if

      read(lin,'(a512)',err=98) jackfl
      call chknam(jackfl,512)
      write(*,*) ' jackknife data file = ',jackfl(1:40)

      read(lin,*,err=98) ixlj,iylj,izlj,ivrlj,iextvj
      write(*,*) ' columns = ',ixlj,iylj,izlj,ivrlj,iextvj

      read(lin,*,err=98) idbg
      write(*,*) ' debugging level = ',idbg

      read(lin,'(a512)',err=98) dbgfl
      call chknam(dbgfl,512)
      write(*,*) ' debugging file = ',dbgfl(1:40)

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

      read(lin,*,err=98) nx,xmn,xsiz
      write(*,*) ' nx, xmn, xsiz = ',nx,xmn,xsiz

      read(lin,*,err=98) ny,ymn,ysiz
      write(*,*) ' ny, ymn, ysiz = ',ny,ymn,ysiz

      read(lin,*,err=98) nz,zmn,zsiz
      write(*,*) ' nz, zmn, zsiz = ',nz,zmn,zsiz

      read(lin,*,err=98) nxdis,nydis,nzdis
      write(*,*) ' block discretization:',nxdis,nydis,nzdis

      read(lin,*,err=98) ndmin,ndmax
      write(*,*) ' ndmin,ndmax = ',ndmin,ndmax

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

      read(lin,*,err=98) ktype,skmean
      write(*,*) ' ktype, skmean =',ktype,skmean

      read(lin,*,err=98) (idrif(i),i=1,9)
      write(*,*) ' drift terms = ',(idrif(i),i=1,9)

      read(lin,*,err=98) itrend
      write(*,*) ' itrend = ',itrend

      read(lin,'(a512)',err=98) extfl
      call chknam(extfl,40)
      write(*,*) ' external drift file = ',extfl(1:40)

      read(lin,*,err=98) iextve
      write(*,*) ' variable in external drift file = ',iextve

      read(lin,*,err=98) nst(1),c0(1)
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
            write(*,*) ' it,cc,ang[1,2,3]; ',it(i),cc(i),
     +                   ang1(i),ang2(i),ang3(i)
            write(*,*) ' a1 a2 a3: ',aa(i),aa1,aa2
            if(it(i).eq.4) then
                  if(aa(i).lt.0.0) stop ' INVALID power variogram'
                  if(aa(i).gt.2.0) stop ' INVALID power variogram'
            end if
      end do

      close(lin)
c
c Find the needed parameters:
c
      MAXDIS = nxdis*nydis*nzdis
      MAXSAM = ndmax + 1
      MAXEQ = MAXSAM + MAXDT + 2
      MAXSBX = 1
      if(nx.gt.1)then
            MAXSBX = int(nx/2.00)
            if(MAXSBX.gt.50)MAXSBX=50
      end if
c
      MAXSBY = 1
      if(ny.gt.1)then
            MAXSBY = int(ny/2.00)
            if(MAXSBY.gt.50)MAXSBY=50
      end if
c
      MAXSBZ = 1
      if(nz.gt.1)then
            MAXSBZ = int(nz/2.00)
            if(MAXSBZ.gt.50)MAXSBZ=50
      end if
c
      MAXSB = MAXSBX*MAXSBY*MAXSBZ
      MXSXY = 4 * MAXSBX * MAXSBY
      MXSX  = 2 * MAXSBX
c
c Allocate the needed memory:
c1
      allocate(nisb(MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c2
      allocate(ixsbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c3
      allocate(iysbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c4
      allocate(izsbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if

c13
      allocate(xa(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c14
      allocate(ya(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c15
      allocate(za(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c16
      allocate(vra(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c17
      allocate(vea(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c18
      allocate(xdb(MAXDIS),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c19
      allocate(ydb(MAXDIS),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c20
      allocate(zdb(MAXDIS),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c23
      allocate(r(MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c24
      allocate(rr(MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c25
      allocate(s(MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c26
      allocate(a(MAXEQ * MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c
c Perform some quick error checking:
c
      if(ndmax.gt.MAXSAM) stop 'ndmax is too big - modify .inc file'
      if(ktype.eq.3.and.iextv.le.0) stop 'must have external variable'
      if(ixl.le.0.and.nx.gt.1) write(*,*) ' WARNING: ixl=0 and nx>1 ! '
      if(iyl.le.0.and.ny.gt.1) write(*,*) ' WARNING: iyl=0 and ny>1 ! '
      if(izl.le.0.and.nz.gt.1) write(*,*) ' WARNING: izl=0 and nz>1 ! '
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
      title(1:22) = 'KT3D ESTIMATES WITH: '
      open(lin,file=datafl,status='OLD')
      read(lin,*)
      read(lin,*,err=99)       nvari
      do i=1,nvari
            read(lin,*)
      end do
      MAXDAT = 0
 22   read(lin,*,end=33,err=99) (var(j),j=1,nvari)
      if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax) go to 22
      MAXDAT = MAXDAT + 1
      go to 22
 33   continue
c
c Allocate the needed memory:
c5
      allocate(x(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c6
      allocate(y(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c7
      allocate(z(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c8
      allocate(vr(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c9
      allocate(ve(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c10
      allocate(dh(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c11
      allocate(tmp(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c12
      allocate(close(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if

      rewind(lin)
      read(lin,'(a58)') title(23:80)
      read(lin,*,err=99)       nvari
      nd = 0
      av = 0.0
      ss = 0.0
      do i=1,nvari
            read(lin,'(a40)',err=99) str
      end do
c
c Some tests on column numbers:
c
      if(ixl.gt.nvari.or.iyl.gt.nvari.or.izl.gt.nvari.or.ivrl.gt.nvari)
     +      then
            write(*,*) 'There are only ',nvari,' columns in input data'
            write(*,*) '  your specification is out of range'
            stop
      end if
c
c Read all the data until the end of the file:
c
 2    read(lin,*,end=3,err=99) (var(j),j=1,nvari)
      vrt = var(ivrl)
      if(vrt.lt.tmin.or.vrt.ge.tmax) go to 2
      nd = nd + 1
      if(nd.gt.MAXDAT) then
            write(*,*) ' ERROR: Exceeded available memory for data'
            stop
      end if
c
c Establish the location of this datum:
c
      if(idhl.le.0) then
            dh(nd) = -99
      else
            dh(nd) = var(idhl)
      endif
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
c
c Establish the external drift variable (if needed):
c
      ve(nd) = 1.0
      if(ktype.eq.3.or.ktype.eq.2) then
            ve(nd) = var(iextv)
            if(ve(nd).lt.tmin.or.ve(nd).ge.tmax) then
                  write(*,*) ' External drift variable must be present',
     +                       ' at all data locations!'
                  write(*,*) ' Encountered at data number ',nd
                  stop
            end if
      end if
      vr(nd) = vrt
      av     = av + vrt
      ss     = ss + vrt*vrt
      go to 2
 3    close(lin)
c
c Compute the averages and variances as an error check for the user:
c
      av = av / max(real(nd),1.0)
      ss =(ss / max(real(nd),1.0)) - av * av
      write(*,*) 'Data for KT3D: Variable number ',ivrl
      write(*,*) '  Number   = ',nd
      write(*,*) '  Average  = ',av
      write(*,*) '  Variance = ',ss
      if(nd.lt.1) then
            write(*,*) ' ERROR: there are no data'
            stop
      end if
c
c Open the debugging and output files:
c
      open(ldbg,file=dbgfl,status='UNKNOWN')

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
      write(lout,'(a80)') title

      if(iktype.eq.0.and.koption.eq.0) then
           write(lout,201) 2,nx,ny,nz
           write(lout,102)
 102       format('Estimate',/,'EstimationVariance')
      end if
      if(iktype.eq.0.and.koption.ge.1) then
           write(lout,201) 7
           write(lout,103)
 103       format('X',/,'Y',/,'Z',/,'True',/,'Estimate',/,
     +            'EstimationVariance',/,'Error: est-true')
      end if
 201  format(4(1x,i4))

      if(iktype.eq.1) then
            if(koption.eq.0) then
                  write(lout,201) ncut,nx,ny,nz
            else
                  write(lout,201) ncut+1
            end if
            do i=1,ncut
                  write(lout,104) i,cut(i)
 104              format('Threshold: ',i2,' = ',f12.5)
            end do
            if(koption.eq.1) write(lout,105)
 105        format('true value')
      end if
#endif

c
c Open the external drift file if needed and position it at the
c first grid node in the file:
c
      if((ktype.eq.2.or.ktype.eq.3).and.koption.eq.0) then
#ifdef _OPENMP
      write(*,*) 'ERROR: multi-core kt3d only supports kriging type'
      write(*,*) '       0=SK or 1=OK, combined with option 0=grid.'
      write(*,*) '       Re-run using single-core kt3d.'
            stop
#endif
            inquire(file=extfl,exist=testfl)
            if(.not.testfl) then
                  write(*,*) 'ERROR file ',extfl,' does not exist!'
                  stop
            endif
            open(lext,file=extfl,status='UNKNOWN')
            read(lext,'(a40)',err=97) str
            read(lext,*,err=97)       nvari
            do i=1,nvari
                  read(lext,'(a40)',err=97) str
            end do
#ifdef DEBUG
            if(idbg.ge.3) write(ldbg,100) iextve
 100        format('A secondary variable is being used.  The gridded '
     +             'file',/,'must have the same grid specifications '
     +             'as the grid you are kriging.',/,'The external '
     +             'drift variable was taken from column ',i2)
#endif
      endif
c
c Set up for cross validation:
c
      if(koption.eq.1) then
            jackfl = datafl
            idhlj  = idhl
            ixlj   = ixl
            iylj   = iyl
            izlj   = izl
            ivrlj  = ivrl
            iextvj = iextv
      end if
c
c Open the file with the jackknife data?
c
      if(koption.gt.0) then
#ifdef _OPENMP
      write(*,*) 'ERROR: multi-core kt3d only supports option 0=grid.'
      write(*,*) '       Re-run using single-core kt3d.'
            stop
#endif
            inquire(file=jackfl,exist=testfl)
            if(.not.testfl) then
                  write(*,*) 'ERROR file ',jackfl,' does not exist!'
                  stop
            endif
            open(ljack,file=jackfl,status='OLD')
            read(ljack,*,err=96)
            read(ljack,*,err=96) nvarij
            do i=1,nvarij
                  read(ljack,*,err=96)
            end do

      end if
c
c Finished here:
c
c      return
      go to 1111
c
c Error in an Input File Somewhere:
c
 96   stop 'ERROR in jackknife file!'
 97   stop 'ERROR in external drift file!'
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'


 1111 print *,'END READING PARAMETERS'
      print *,'START KRIGING'


ccccccccccccccc kt3d cccccccccccccccccccccc


c
c Set up the rotation/anisotropy matrices that are needed for the
c variogram and search.  Also compute the maximum covariance for
c the rescaling factor:
c
      write(*,*) 'Setting up rotation matrices for variogram and search'
      radsqd = radius * radius
      PMX    = 999.0
      covmax = c0(1)
      do is=1,nst(1)
            call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is),
     +                  is,MAXROT,rotmat)
            if(it(is).eq.4) then
                  covmax = covmax + PMX 
            else
                  covmax = covmax + cc(is)
            endif
      end do
      isrot = MAXNST + 1
      call setrot(sang1,sang2,sang3,sanis1,sanis2,isrot,MAXROT,rotmat)

c
c Finish computing the rescaling factor and stop if unacceptable:
c
      if(radsqd.lt.1.0) then
            resc = 2.0 * radius / max(covmax,0.0001)
      else
            resc =(4.0 * radsqd)/ max(covmax,0.0001)
      endif
      if(resc.le.0.0) then
            write(*,*) 'ERROR KT3D: The rescaling value is wrong ',resc
            write(*,*) '            Maximum covariance: ',covmax
            write(*,*) '            search radius:      ',radius
            stop
      endif
      resc = 1.0 / resc
c
c Set up for super block searching:
c
      write(*,*) 'Setting up super block search strategy'
      nsec = 2
      call setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z,
     +             vr,tmp,nsec,ve,dh,sec3,MAXSBX,MAXSBY,MAXSBZ,nisb,
     +             nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,nzsup,
     +             zmnsup,zsizsup)
      call picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,
     +             isrot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr,
     +             iysbtosr,izsbtosr)

c
c Compute the number of drift terms, if an external drift is being
c considered then it is one more drift term, if SK is being considered
c then we will set all the drift terms off and mdt to 0):
c
      mdt = 1
      do i=1,9
            if(ktype.eq.0.or.ktype.eq.2) idrif(i) = 0
            if(idrif(i).lt.0.or.idrif(i).gt.1) then
                  write(*,*) 'ERROR KT3D: invalid drift term',idrif(i)
                  stop
            endif
            mdt = mdt + idrif(i)
      end do
      if(ktype.eq.3) mdt = mdt + 1
      if(ktype.eq.0) mdt = 0
      if(ktype.eq.2) mdt = 0
c
c Set up the discretization points per block.  Figure out how many
c are needed, the spacing, and fill the xdb,ydb, and zdb arrays with
c the offsets relative to the block center (this only gets done once):
c
c In all cases the offsets are relative to the lower left corner.
c This is done for rescaling the drift terms in the kriging matrix.
c
      if(nxdis.lt.1) nxdis = 1
      if(nydis.lt.1) nydis = 1
      if(nzdis.lt.1) nzdis = 1
      ndb = nxdis * nydis * nzdis
      if(ndb.gt.MAXDIS) then
            write(*,*) 'ERROR KT3D: Too many discretization points',ndb
            write(*,*) '            Increase MAXDIS or lower n[xyz]dis'
            stop
      endif
      xdis = xsiz  / max(real(nxdis),1.0)
      ydis = ysiz  / max(real(nydis),1.0)
      zdis = zsiz  / max(real(nzdis),1.0)
      i    = 0
      xloc = -0.5*(xsiz+xdis)
      do ix =1,nxdis
            xloc = xloc + xdis
            yloc = -0.5*(ysiz+ydis)
            do iy=1,nydis
                  yloc = yloc + ydis
                  zloc = -0.5*(zsiz+zdis)
                  do iz=1,nzdis
                        zloc = zloc + zdis
                        i = i+1
                        xdb(i) = xloc + 0.5*xsiz
                        ydb(i) = yloc + 0.5*ysiz
                        zdb(i) = zloc + 0.5*zsiz
                  end do
            end do
      end do
c
c Initialize accumulators:
c
      nk    = 0
      xk    = 0.0
      vk    = 0.0
      xkmae = 0.0
      xkmse = 0.0
c
c Calculate Block Covariance. Check for point kriging.
c
      call cova3(xdb(1),ydb(1),zdb(1),xdb(1),ydb(1),zdb(1),1,nst,MAXNST,
     +           c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)

c
c Set the ``unbias'' variable so that the matrix solution is more stable
c
      unbias = cov
      cbb    = dble(cov)
      if(ndb.gt.1) then
            cbb = 0.0
            do i=1,ndb
               do j=1,ndb
                  call cova3(xdb(i),ydb(i),zdb(i),xdb(j),ydb(j),zdb(j),
     +               1,nst,MAXNST,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
                  if(i.eq.j) cov = cov - c0(1)
                  cbb = cbb + dble(cov)
               end do
            end do
            cbb = cbb/dble(real(ndb*ndb))
      end if
#ifdef DEBUG
      if(idbg.gt.1) then
            write(ldbg,*) ' '
            write(ldbg,*) 'Block Covariance: ',cbb
            write(ldbg,*) ' '
      end if
#endif
c
c Mean values of the drift functions:
c
      do i=1,9
            bv(i) = 0.0
      end do
      do i=1,ndb
            bv(1) = bv(1) + xdb(i)
            bv(2) = bv(2) + ydb(i)
            bv(3) = bv(3) + zdb(i)
            bv(4) = bv(4) + xdb(i)*xdb(i)
            bv(5) = bv(5) + ydb(i)*ydb(i)
            bv(6) = bv(6) + zdb(i)*zdb(i)
            bv(7) = bv(7) + xdb(i)*ydb(i)
            bv(8) = bv(8) + xdb(i)*zdb(i)
            bv(9) = bv(9) + ydb(i)*zdb(i)
      end do  
      do i=1,9
            bv(i) = (bv(i) / real(ndb)) * resc
      end do  
c
c Report on progress from time to time:
c
      if(koption.eq.0) then
            nxy   = nx*ny
            nxyz  = nx*ny*nz
            nloop = nxyz
            irepo = max(1,min((nxyz/10),10000))
      else
            nloop = 10000000
            irepo = max(1,min((nd/10),10000))
      end if
      ddh = 0.0
      write(*,*)
      write(*,*) 'Working on the kriging '

#ifdef USE_MPI
      nloopIni=(nloop/mpisize)*myrank +1 
      nloopFin=(nloop/mpisize)*(myrank+1)
      if(nloopFin.gt.nloop) then
         nloopFin=nloop
      end if
#endif

c
c MAIN LOOP OVER ALL THE BLOCKS IN THE GRID:
c

c$omp parallel default(firstprivate) 
c$omp& shared(x,y,z,vr,ixsbtosr,iysbtosr,
c$omp&        izsbtosr,rotmat,covtab)
#ifdef _OPENMP
      threadId = int(OMP_get_thread_num())
#else
      threadId = 0
#endif

c$omp do schedule(static)
#ifdef USE_MPI
      do index=nloopIni,nloopFin
#else
      do index=1,nloop
#endif

#ifndef _OPENMP

      if((int(index/irepo)*irepo).eq.index) write(*,1003) index
 1003 format('   currently on estimate ',i9)

#else

      if(threadId.eq.0.and.myrank.eq.0)then
      if((int(index/irepo)*irepo).eq.index) write(*,1003) 
     + (index*numThreads*mpisize)
 1003 format('   currently on estimate ',i9)
      endif

#endif


c
c Where are we making an estimate?
c
      if(koption.eq.0) then
            iz   = int((index-1)/nxy) + 1
            iy   = int((index-(iz-1)*nxy-1)/nx) + 1
            ix   = index - (iz-1)*nxy - (iy-1)*nx
            xloc = xmn + real(ix-1)*xsiz
            yloc = ymn + real(iy-1)*ysiz
            zloc = zmn + real(iz-1)*zsiz
      else
#ifndef _OPENMP
            read(ljack,*,err=96,end=2222) (var(i),i=1,nvarij)
            ddh  = 0.0
            xloc = xmn
            yloc = ymn
            zloc = zmn
            true = UNEST
            secj = UNEST
            if(idhlj.gt.0)  ddh    = var(idhlj)
            if(ixlj.gt.0)   xloc   = var(ixlj)
            if(iylj.gt.0)   yloc   = var(iylj)
            if(izlj.gt.0)   zloc   = var(izlj)
            if(ivrlj.gt.0)  true   = var(ivrlj)
            if(iextvj.gt.0) extest = var(iextvj)
            if(true.lt.tmin.or.true.ge.tmax) true = UNEST
#endif
      end if

c
c Read in the external drift variable for this grid node if needed:
c
      if(ktype.eq.2.or.ktype.eq.3) then
            if(koption.eq.0) then
                  read(lext,*) (var(i),i=1,iextve)
                  extest = var(iextve)
            end if
            if(extest.lt.tmin.or.extest.ge.tmax) then
                  est  = UNEST
                  estv = UNEST
                  go to 111
            end if
            resce  = covmax / max(extest,0.0001)
      endif
c
c Find the nearest samples:
c

      close(:)=0.0
      call srchsupr(xloc,yloc,zloc,radsqd,isrot,MAXROT,rotmat,nsbtosr,
     +              ixsbtosr,iysbtosr,izsbtosr,noct,nd,x,y,z,tmp,
     +              nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,
     +              nzsup,zmnsup,zsizsup,nclose,close,infoct)

c      call srchsuprO1(xloc,yloc,zloc,radsqd,isrot,MAXROT,rotmat,nsbtosr,
c     +              ixsbtosr,iysbtosr,izsbtosr,noct,nd,x,y,z,tmp,
c     +              nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,
c     +              nzsup,zmnsup,zsizsup,nclose,close,infoct)

c
c Load the nearest data in xa,ya,za,vra,vea:
c

c      vra(:)=0.0
c      vea(:)=0.0
      na = 0

      do i=1,nclose
            ind    = int(close(i)+0.5)
            accept = .true.
            if(koption.ne.0.and.
     +         (abs(x(ind)-xloc)+abs(y(ind)-yloc)+ abs(z(ind)-zloc))
     +                           .lt.EPSLON) accept = .false.
            if(koption.ne.0.and.
     +         (abs(dh(ind)-ddh)).lt.EPSLON) accept = .false.
            if(accept) then
                  if(na.lt.ndmax) then
                        na = na + 1
                        xa(na)  = x(ind) - xloc + 0.5*xsiz
                        ya(na)  = y(ind) - yloc + 0.5*ysiz
                        za(na)  = z(ind) - zloc + 0.5*zsiz
                        vra(na) = vr(ind)
                        vea(na) = ve(ind)
                  end if
            end if
      end do
c
c Test number of samples found:
c
      if(na.lt.ndmin) then
            est  = UNEST
            estv = UNEST
            go to 111
      end if
c
c Test if there are enough samples to estimate all drift terms:
c
      if(na.ge.1.and.na.le.mdt) then
            if(fircon) then
                  write(ldbg,999)
                  fircon = .false.
            end if
            est  = UNEST
            estv = UNEST
            go to 111
      end if
 999  format(' Encountered a location where there were too few data ',/,
     +       ' to estimate all of the drift terms but there would be',/,
     +       ' enough data for OK or SK.   KT3D currently leaves ',/,
     +       ' these locations unestimated.',/,
     +       ' This message is only written once - the first time.',/)
c
c There are enough samples - proceed with estimation.
c
      if(na.le.1) then
c
c Handle the situation of only one sample:
c
c            call cova3(xa(1),ya(1),za(1),xa(1),ya(1),za(1),1,nst,MAXNST,
c     +                 c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb1)
          call cova3O1(xa(1),ya(1),za(1),xa(1),ya(1),za(1),1,nst,MAXNST,
     +                 c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb1)
c
c Establish Right Hand Side Covariance:
c
            if(ndb.le.1) then
c                  call cova3(xa(1),ya(1),za(1),xdb(1),ydb(1),zdb(1),1,
c     +                 nst,MAXNST,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb)
                call cova3O1(xa(1),ya(1),za(1),xdb(1),ydb(1),zdb(1),1,
     +                 nst,MAXNST,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb)
            else
                  cb  = 0.0
                  do i=1,ndb
c                        call cova3(xa(1),ya(1),za(1),xdb(i),ydb(i),
c     +                             zdb(i),1,nst,MAXNST,c0,it,cc,aa,1,
c     +                             MAXROT,rotmat,cmax,cov)
                      call cova3O1(xa(1),ya(1),za(1),xdb(i),ydb(i),
     +                             zdb(i),1,nst,MAXNST,c0,it,cc,aa,1,
     +                             MAXROT,rotmat,cmax,cov)
                        cb = cb + cov
                        dx = xa(1) - xdb(i)
                        dy = ya(1) - ydb(i)
                        dz = za(1) - zdb(i)
                        if((dx*dx+dy*dy+dz*dz).lt.EPSLON) cb=cb-c0(1)
                  end do
                  cb = cb / real(ndb)
            end if
cvd
c
c Early bug - always did OK in presence of one data.
c
cvd
            if(ktype.eq.2) skmean = extest
            if(ktype.eq.0.or.ktype.eq.2) then
                  wt   = cb / cb1
                  est  = wt * vra(1) + (1.0-wt) * skmean
                  estv = real(cbb) - wt*cb
            else
                  est  = vra(1)
                  estv = real(cbb) - 2.0*cb + cb1
            end if
            nk   = nk + 1
            xk   = xk + est
            vk   = vk + est*est
            go to 111
      end if
c
c Go ahead and set up the OK portion of the kriging matrix:
c
      neq = mdt+na
c
c Initialize the main kriging matrix:
c
      first = .false.
      do i=1,neq*neq
            a(i) = 0.0
      end do
c
c Fill in the kriging matrix:
c
      do i=1,na
      do j=i,na
c            call cova3(xa(i),ya(i),za(i),xa(j),ya(j),za(j),1,nst,MAXNST,
c     +                 c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
          call cova3O1(xa(i),ya(i),za(i),xa(j),ya(j),za(j),1,nst,MAXNST,
     +                 c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
            a(neq*(i-1)+j) = dble(cov)
            a(neq*(j-1)+i) = dble(cov)
      end do
      end do
c
c Fill in the OK unbiasedness portion of the matrix (if not doing SK):
c
      if(neq.gt.na) then
            do i=1,na
                  a(neq*(i-1)+na+1) = dble(unbias)
                  a(neq*na+i)       = dble(unbias)
            end do
      endif
c
c Set up the right hand side:
c
      do i=1,na
            if(ndb.le.1) then
c                  call cova3(xa(i),ya(i),za(i),xdb(1),ydb(1),zdb(1),1,
c     +                 nst,MAXNST,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb)
                call cova3O1(xa(i),ya(i),za(i),xdb(1),ydb(1),zdb(1),1,
     +                 nst,MAXNST,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb)
            else
                  cb  = 0.0
                  do j=1,ndb
c                        call cova3(xa(i),ya(i),za(i),xdb(j),ydb(j),
c     +                             zdb(j),1,nst,MAXNST,c0,it,cc,aa,1,
c     +                             MAXROT,rotmat,cmax,cov)
                      call cova3O1(xa(i),ya(i),za(i),xdb(j),ydb(j),
     +                             zdb(j),1,nst,MAXNST,c0,it,cc,aa,1,
     +                             MAXROT,rotmat,cmax,cov)
                        cb = cb + cov
                        dx = xa(i) - xdb(j)
                        dy = ya(i) - ydb(j)
                        dz = za(i) - zdb(j)
                        if((dx*dx+dy*dy+dz*dz).lt.EPSLON) cb=cb-c0(1)
                  end do
                  cb = cb / real(ndb)
            end if
            r(i) = dble(cb)
      end do
      if(neq.gt.na) r(na+1) = dble(unbias)
c
c Add the additional unbiasedness constraints:
c
      im = na + 1
#ifndef _OPENMP
c
c First drift term (linear in "x"):
c
      if(idrif(1).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(xa(k)*resc)
                  a(neq*(k-1)+im) = dble(xa(k)*resc)
            end do
            r(im) = dble(bv(1))
      endif
c
c Second drift term (linear in "y"):
c
      if(idrif(2).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(ya(k)*resc)
                  a(neq*(k-1)+im) = dble(ya(k)*resc)
            end do
            r(im) = dble(bv(2))
      endif
c
c Third drift term (linear in "z"):
c
      if(idrif(3).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(za(k)*resc)
                  a(neq*(k-1)+im) = dble(za(k)*resc)
            end do
            r(im) = dble(bv(3))
      endif
c
c Fourth drift term (quadratic in "x"):
c
      if(idrif(4).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(xa(k)*xa(k)*resc)
                  a(neq*(k-1)+im) = dble(xa(k)*xa(k)*resc)
            end do
            r(im) = dble(bv(4))
      endif
c
c Fifth drift term (quadratic in "y"):
c
      if(idrif(5).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(ya(k)*ya(k)*resc)
                  a(neq*(k-1)+im) = dble(ya(k)*ya(k)*resc)
            end do
            r(im) = dble(bv(5))
      endif
c
c Sixth drift term (quadratic in "z"):
c
      if(idrif(6).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(za(k)*za(k)*resc)
                  a(neq*(k-1)+im) = dble(za(k)*za(k)*resc)
            end do
            r(im) = dble(bv(6))
      endif
c
c Seventh drift term (quadratic in "xy"):
c
      if(idrif(7).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(xa(k)*ya(k)*resc)
                  a(neq*(k-1)+im) = dble(xa(k)*ya(k)*resc)
            end do
            r(im) = dble(bv(7))
      endif
c
c Eighth drift term (quadratic in "xz"):
c
      if(idrif(8).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(xa(k)*za(k)*resc)
                  a(neq*(k-1)+im) = dble(xa(k)*za(k)*resc)
            end do
            r(im) = dble(bv(8))
      endif
c
c Ninth drift term (quadratic in "yz"):
c
      if(idrif(9).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(ya(k)*za(k)*resc)
                  a(neq*(k-1)+im) = dble(ya(k)*za(k)*resc)
            end do
            r(im) = dble(bv(9))
      endif
c
c External drift term (specified by external variable):
c
      if(ktype.eq.3) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(vea(k)*resce)
                  a(neq*(k-1)+im) = dble(vea(k)*resce)
            end do
            r(im) = dble(extest*resce)
      endif
#endif
c
c Copy the right hand side to compute the kriging variance later:
c
      do k=1,neq
            rr(k) = r(k)
      end do
      kadim = neq * neq
      ksdim = neq
      nrhs  = 1
      nv    = 1
c
c If estimating the trend then reset all the right hand side terms=0.0:
c
      if(itrend.ge.1) then
            do i=1,na
                  r(i)  = 0.0
                  rr(i) = 0.0
            end do
      endif
c
c Write out the kriging Matrix if Seriously Debugging:
c
#ifdef DEBUG
      if(idbg.eq.3) then
            write(ldbg,*) 'Estimating node index : ',ix,iy,iz
            is = 1 - neq
            do i=1,neq
                  is = 1 + (i-1)*neq
                  ie = is + neq - 1
                  write(ldbg,1000) i,r(i),(a(j),j=is,ie)
 1000             format('    r(',i2,') =',f7.4,'  a= ',9(10f7.4))
            end do
      endif
#endif
c
c Solve the kriging system:
c
      call ktsol(neq,nrhs,nv,a,r,s,ising,maxeq)
c
c Compute the solution:
c
      if(ising.ne.0) then
#ifdef DEBUG
            if(idbg.ge.3) write(ldbg,*) ' Singular Matrix ',ix,iy,iz
#endif
            est  = UNEST
            estv = UNEST
      else
            est  = 0.0
            estv = real(cbb)
            if(ktype.eq.2) skmean = extest
            do j=1,neq
                  estv = estv - real(s(j))*rr(j)
                  if(j.le.na) then
                        if(ktype.eq.0) then
                              est = est + real(s(j))*(vra(j)-skmean)
                        else if(ktype.eq.2) then
                              est = est + real(s(j))*(vra(j)-vea(j))
                        else
                              est = est + real(s(j))*vra(j)
                        endif
                  endif
            end do
            if(ktype.eq.0.or.ktype.eq.2) est = est + skmean
            nk   = nk + 1
            xk   = xk + est
            vk   = vk + est*est
c
c Write the kriging weights and data if debugging level is above 2:
c
#ifdef DEBUG
            if(idbg.ge.2) then
                  write(ldbg,*) '       '
                  write(ldbg,*) 'BLOCK: ',ix,iy,iz,' at ',xloc,yloc,zloc
                  write(ldbg,*) '       '
                  if(ktype.ne.0) 
     +            write(ldbg,*) '  Lagrange : ',s(na+1)*unbias
                  write(ldbg,*) '  BLOCK EST: x,y,z,vr,wt '
                  do i=1,na
                        xa(i) = xa(i) + xloc - 0.5*xsiz
                        ya(i) = ya(i) + yloc - 0.5*ysiz
                        za(i) = za(i) + zloc - 0.5*zsiz
                        write(ldbg,'(5f12.3)') xa(i),ya(i),za(i),
     +                                         vra(i),s(i)
                  end do
                  write(ldbg,*) '  estimate, variance  ',est,estv
            endif
#endif
      endif
c
c END OF MAIN KRIGING LOOP:
c
 111        continue
            if(iktype.eq.0) then
                  if(koption.eq.0) then
#ifdef UNFORMATTED
                        write(loutThreads(threadId+1)) est,estv
#else
                        write(loutThreads(threadId+1),
     + '(g14.8,1x,g14.8)') est,estv
#endif
                  else
                        err = UNEST
                        if(true.ne.UNEST.and.est.ne.UNEST) then
                              err=est-true
                              xkmae = xkmae + abs(err)
                              xkmse = xkmse + err*err
                        end if

#ifdef UNFORMATTED
                        write(loutThreads(threadId+1))
     + xloc,yloc,zloc,true,est,estv,err
#else
                        write(loutThreads(threadId+1),'(7(g14.8,1x))') 
     + xloc,yloc,zloc,true,est,estv,err
#endif

                  end if
            else
c
c Work out the IK-type distribution implicit to this data configuration
c and kriging weights:
c
                  do icut=1,ncut
                        cdf(icut) = -1.0
                  end do
                  wtmin = 1.0
                  do i=1,na
                        if(s(i).lt.wtmin) wtmin = s(i)
                  end do
                  sumwt = 0.0
                  do i=1,na
                        s(i)  = s(i) - wtmin
                        sumwt = sumwt + s(i)
                  end do
                  do i=1,na
                        s(i) = s(i) / max(0.00001,sumwt)
                  end do
                  if(na.gt.1.and.sumwt.gt.0.00001) then
                        do icut=1,ncut
                              cdf(icut) = 0.0
                              do i=1,na
                                    if(vra(i).le.cut(icut))
     +                              cdf(icut)=cdf(icut)+s(i)
                              end do
                        end do
                  end if
                  if(koption.eq.0) then
#ifdef UNFORMATTED
                        write(loutThreads(threadId+1)) 
     +                           (cdf(i),i=1,ncut)
#else
                        write(loutThreads(threadId+1),'(30(f8.4))') 
     +                           (cdf(i),i=1,ncut)
#endif
                  else
#ifdef UNFORMATTED
                        write(loutThreads(threadId+1)) 
     +                           (cdf(i),i=1,ncut),true
#else
                        write(loutThreads(threadId+1),'(30(f8.4))') 
     +                           (cdf(i),i=1,ncut),true
#endif
                  end if
            end if
      end do
c$omp end do nowait

c$omp end parallel


#ifndef _OPENMP
 2222 continue
      if(koption.gt.0) close(ljack)
#endif

c
c Write statistics of kriged values:
c
 
#ifdef DEBUG
      if(nk.gt.0.and.idbg.gt.0) then
            xk    = xk/real(nk)
            vk    = vk/real(nk) - xk*xk
            xkmae = xkmae/real(nk)
            xkmse = xkmse/real(nk)
            write(ldbg,1005) nk,xk,vk
            write(*,   1005) nk,xk,vk
 1005       format(/,'Estimated   ',i8,' blocks ',/,
     +               '  average   ',g14.8,/,'  variance  ',g14.8,/)
            if(koption.ne.0) then
                  write(*,1006) xkmae,xkmse
 1006             format(/,'  mean error',g14.8,/,'  mean sqd e',g14.8)
            end if
      endif
#endif
c
c All finished the kriging:
c
c      return
c 960  stop 'ERROR in jackknife file!'

c
c Finished:
c
      close(ldbg)
      close(lout)
#ifdef _OPENMP
      do i=1,numThreads
            close(loutThreads(i))
      end do
#endif
      write(*,9998) VERSION
 9998 format(/' KT3D Version: ',f5.3, ' Finished'/)

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
      open(lun,file='kt3d.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for KT3D',/,
     +       '                  *******************',/,/,
     +       'START OF PARAMETERS:')

      write(lun,11)
 11   format('../data/cluster.dat              ',
     +       '-file with data')
      write(lun,12)
 12   format('0  1  2  0  3  0                 ',
     +       '-   columns for DH,X,Y,Z,var,sec var')
      write(lun,13)
 13   format('-1.0e21   1.0e21                 ',
     +       '-   trimming limits')
      write(lun,14)
 14   format('0                                ',
     +       '-option: 0=grid, 1=cross, 2=jackknife')
      write(lun,15)
 15   format('xvk.dat                          ',
     +       '-file with jackknife data')
      write(lun,16)
 16   format('1   2   0    3    0              ',
     +       '-   columns for X,Y,Z,vr and sec var')
      write(lun,17)
 17   format('3                                ',
     +       '-debugging level: 0,1,2,3')
      write(lun,18)
 18   format('kt3d.dbg                         ',
     +       '-file for debugging output')
      write(lun,19)
 19   format('kt3d.out                         ',
     +       '-file for kriged output')
      write(lun,20)
 20   format('50   0.5    1.0                  ',
     +       '-nx,xmn,xsiz')
      write(lun,21)
 21   format('50   0.5    1.0                  ',
     +       '-ny,ymn,ysiz')
      write(lun,22)
 22   format('1    0.5    1.0                  ',
     +       '-nz,zmn,zsiz')
      write(lun,23)
 23   format('1    1      1                    ',
     +       '-x,y and z block discretization')
      write(lun,24)
 24   format('4    8                           ',
     +       '-min, max data for kriging')
      write(lun,25)
 25   format('0                                ',
     +       '-max per octant (0-> not used)')
      write(lun,26)
 26   format('20.0  20.0  20.0                 ',
     +       '-maximum search radii')
      write(lun,27)
 27   format(' 0.0   0.0   0.0                 ',
     +       '-angles for search ellipsoid')
      write(lun,28)
 28   format('0     2.302                      ',
     +       '-0=SK,1=OK,2=non-st SK,3=exdrift')
      write(lun,29)
 29   format('0 0 0 0 0 0 0 0 0                ',
     +       '-drift: x,y,z,xx,yy,zz,xy,xz,zy')
      write(lun,30)
 30   format('0                                ',
     +       '-0, variable; 1, estimate trend')
      write(lun,31)
 31   format('extdrift.dat                     ',
     +       '-gridded file with drift/mean')
      write(lun,32)
 32   format('4                                ',
     +       '-  column number in gridded file')
      write(lun,33)
 33   format('1    0.2                         ',
     +       '-nst, nugget effect')
      write(lun,34)
 34   format('1    0.8  0.0   0.0   0.0        ',
     +       '-it,cc,ang1,ang2,ang3')
      write(lun,35)
 35   format('         10.0  10.0  10.0        ',
     +       '-a_hmax, a_hmin, a_vert')

      close(lun)
      return
      end
