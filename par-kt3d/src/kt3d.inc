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
c-----------------------------------------------------------------------
c
c                  Kriging of a 3-D Rectangular Grid
c                  *********************************
c
c This file controls both the static dimensioning and allocates most of
c the storage in common blocks.
c
c   MAXSBX    maximum super block nodes in X direction
c   MAXSBY    maximum super block nodes in Y direction
c   MAXSBZ    maximum super block nodes in Z direction
c   MAXDAT    maximum number of data points
c   MAXDT     maximum number of drift terms
c   MAXSAM    maximum number of data points to use in one kriging system
c   MAXDIS    maximum number of discretization points per block
c   MAXNST    maximum number of nested structures
c
c-----------------------------------------------------------------------
c
c User Adjustable Parameters:
c
      parameter(MAXDT  =   9, MAXNST =   4)
c
c Fixed Parameters:
c
      parameter(UNEST = -999.0, EPSLON = 0.000001, MAXROT = MAXNST + 1,
     +          VERSION = 3.000)
c
c Static Array Dimensioning:
c
      integer   idrif(MAXDT),nst(1),it(MAXNST),iktype,ncut,test
      real      bv(9),c0(1),cc(MAXNST),aa(MAXNST),aainv(MAXNST),
     +          ang1(MAXNST),
     +          ang2(MAXNST),ang3(MAXNST),anis1(MAXNST),anis2(MAXNST)
      real*8    rotmat(MAXROT,3,3)
c      real*8    rotmat(3,3,MAXROT)
c
c The data and other input variables:
c
      common /datcom/ nd,tmin,tmax,nx,ny,nz,xmn,ymn,zmn,test,
     +                xsiz,ysiz,zsiz,ndmax,ndmin,radius,noct,nxdis,
     +                nydis,nzdis,idrif,itrend,ktype,skmean,koption,
     +                idbg,ldbg,iout,lout,iext,lext,iextve,ljack,
     +                idhlj,ixlj,iylj,izlj,ivrlj,iextvj,nvarij
c
c Variogram Parameters:
c
      common /vargdt/ nst,it,c0,cc,aa,ang1,ang2,ang3,anis1,anis2,rotmat
c
c Search variables and data for kriging:
c
      common /srccom/ sang1,sang2,sang3,sanis1,sanis2,isrot,na,ndb,
     +                bv,unbias,iktype,ncut
