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
c The following Parameters control static dimensioning for SGSIM:
c
c
c-----------------------------------------------------------------------
c
c Parameters:
c
      parameter(MAXNST=4,MAXROT=MAXNST+1,UNEST=-99.0,EPSLON=1.0e-20,
     +          VERSION=3.001)
c
c ACORN parameters:
c
c      parameter(KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)

c
c Variable Declaration:
c
      real      c0(1),cc(MAXNST),aa(MAXNST),aainv(MAXNST),
     +          ang1(MAXNST),ang2(MAXNST),
     +          ang3(MAXNST),anis1(MAXNST),anis2(MAXNST),ltpar
      real*8    rotmat(MAXROT,3,3)
      integer   ltail,utail,sstrat,mults,nst(1),it(MAXNST),test
c
c Common Blocks:
c
c      common /iaco/   ixv(MAXOP1)

      common /generl/ nsim,nx,ny,nz,xsiz,ysiz,zsiz,xmn,ymn,zmn,nxy,nxyz,
     +                nd,itrans,ntr,idbg,lin,lout,ldbg,llvm,icollvm,
     +                nvaril,ktype,colocorr,test,
     +                ltail,ltpar,utail,utpar,zmin,zmax,varred

      common /search/ radius,radsqd,sang1,sang2,sang3,sanis1,sanis2,
c     +                noct,nclose,ndmin,ndmax,sstrat,mults,nmult,
     +                noct,nclose,ndmin      ,sstrat,mults,nmult
c     +                mxctx,mxcty,mxctz

c      common /cova3d/ nst,it,cmax,c0,cc,aa,ang1,ang2,ang3,
      common /cova3d/ nst,it,ang1,ang2,ang3,
     +                anis1,anis2,isrot

c      common /clooku/ nctx,ncty,nctz,nlooku,ncnode,nodmax

      common /krigev/ rotmat,cbb
