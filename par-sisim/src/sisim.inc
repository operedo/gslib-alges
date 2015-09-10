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
c Parameters:
c
      parameter(MAXNST=4,EPSLON=1.0e-20,UNEST=-99.0,VERSION=3.000)
c
c ACORN parameters:
c
      parameter(KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
c
c Static Array Dimensioning:
c
      integer inoct(8),ltail,middle,utail,sstrat,mults,nmult,ivtype,
     +        imbsim
      real    ltpar,mpar,utpar
c
c Common Blocks:
c
      common /iaco/   ixv(MAXOP1)

      common /grid3d/ nx,ny,nz,nxy,nxyz,xsiz,ysiz,zsiz,xmn,ymn,zmn

      common /data3d/ ivtype,nhd,nd,zmin,zmax,ltail,middle,
     +                utail,ltpar,mpar,utpar,tmin,tmax,mxctx,
     +                mxcty,mxctz

      common /simula/ nsim,ktype,lin,lout,idbg,ldbg

      common /srch3d/ radius,radsqd,sang1,sang2,sang3,sanis1,sanis2,
c     +                noct,nclose,ndmax,mik,sstrat,nodmax,
     +                noct,nclose,ndmax,mik,sstrat,
     +                infoct,inoct,mults,nmult
c     +                infoct,inoct,mults,nmult,maxsec

      common /cova3d/ ncut,imbsim,ng,isrot,ivrot,cbb

c      common /clooku/ nctx,ncty,nctz,nlooku,ncnode
c      common /clooku/ nctx,ncty,nctz,nlooku
