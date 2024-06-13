      subroutine ksol(nright,neq,nsb,a,r,s,ising)
c-----------------------------------------------------------------------
c
c                Solution of a System of Linear Equations
c                ****************************************
c
c
c
c INPUT VARIABLES:
c
c   nright,nsb       number of columns in right hand side matrix.
c                      for KB2D: nright=1, nsb=1
c   neq              number of equations
c   a()              upper triangular left hand side matrix (stored 
c                      columnwise)
c   r()              right hand side matrix (stored columnwise)
c                      for kb2d, one column per variable
c
c
c
c OUTPUT VARIABLES:
c
c   s()              solution array, same dimension as  r  above.
c   ising            singularity indicator
c                      0,  no singularity problem
c                     -1,  neq .le. 1
c                      k,  a null pivot appeared at the kth iteration
c
c
c
c PROGRAM NOTES:
c
c   1. Requires the upper triangular left hand side matrix.
c   2. Pivots are on the diagonal.
c   3. Does not search for max. element for pivot.
c   4. Several right hand side matrices possible.
c   5. USE for ok and sk only, NOT for UK.
c
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8   a(*),r(*),s(*)
c
c If there is only one equation then set ising and return:
c
      if(neq.le.1) then
            ising = -1
            return
      endif
c
c Initialize:
c
      tol   = 0.1e-06
      ising = 0
      nn    = neq*(neq+1)/2
      nm    = nsb*neq
      m1    = neq-1
      kk    = 0
c
c Start triangulation:
c
      do k=1,m1
            kk=kk+k
            ak=a(kk)
            if(abs(ak).lt.tol) then
                  ising=k
                  return
            endif
            km1=k-1
            do iv=1,nright
                  nm1=nm*(iv-1)
                  ii=kk+nn*(iv-1)
                  piv=1./a(ii)
                  lp=0
                  do i=k,m1
                        ll=ii
                        ii=ii+i
                        ap=a(ii)*piv
                        lp=lp+1
                        ij=ii-km1
                        do j=i,m1
                              ij=ij+j
                              ll=ll+j
                              a(ij)=a(ij)-ap*a(ll)
                        end do
                        do llb=k,nm,neq
                              in=llb+lp+nm1
                              ll1=llb+nm1
                              r(in)=r(in)-ap*r(ll1)
                        end do
                  end do
            end do
      end do
c
c Error checking - singular matrix:
c
      ijm=ij-nn*(nright-1)
      if(abs(a(ijm)).lt.tol) then
            ising=neq
            return
      endif
c
c Finished triangulation, start solving back:
c
      do iv=1,nright
            nm1=nm*(iv-1)
            ij=ijm+nn*(iv-1)
            piv=1./a(ij)
            do llb=neq,nm,neq
                  ll1=llb+nm1
                  s(ll1)=r(ll1)*piv
            end do
            i=neq
            kk=ij
            do ii=1,m1
                  kk=kk-i
                  piv=1./a(kk)
                  i=i-1
                  do llb=i,nm,neq
                        ll1=llb+nm1
                        in=ll1
                        ap=r(in)
                        ij=kk
                        do j=i,m1
                              ij=ij+j
                              in=in+1
                              ap=ap-a(ij)*s(in)
                        end do
                        s(ll1)=ap*piv
                  end do
            end do
      end do
c
c Finished solving back, return:
c
      return
      end

      subroutine ksolopt(nright,neq,nsb,a,r,s,ising)
c-----------------------------------------------------------------------
c
c                Solution of a System of Linear Equations
c                ****************************************
c
c
c
c INPUT VARIABLES:
c
c   nright,nsb       number of columns in right hand side matrix.
c                      for KB2D: nright=1, nsb=1
c   neq              number of equations
c   a()              upper triangular left hand side matrix (stored 
c                      columnwise)
c   r()              right hand side matrix (stored columnwise)
c                      for kb2d, one column per variable
c
c
c
c OUTPUT VARIABLES:
c
c   s()              solution array, same dimension as  r  above.
c   ising            singularity indicator
c                      0,  no singularity problem
c                     -1,  neq .le. 1
c                      k,  a null pivot appeared at the kth iteration
c
c
c
c PROGRAM NOTES:
c
c   1. Requires the upper triangular left hand side matrix.
c   2. Pivots are on the diagonal.
c   3. Does not search for max. element for pivot.
c   4. Several right hand side matrices possible.
c   5. USE for ok and sk only, NOT for UK.
c
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8   a(*),r(*),s(*)
c
c If there is only one equation then set ising and return:
c
      if(neq.le.1) then
            ising = -1
            return
      endif
c
c Initialize:
c
      tol   = 0.1e-06
      ising = 0
      nn    = neq*(neq+1)/2
      nm    = nsb*neq
      m1    = neq-1
      kk    = 0
c
c Start triangulation:
c
      do k=1,m1
            kk=kk+k
            ak=a(kk)
            if(abs(ak).lt.tol) then
                  ising=k
                  return
            endif
            km1=k-1
c            do iv=1,nright
c                  nm1=nm*(iv-1)
                  nm1=0
c                  ii=kk+nn*(iv-1)
                  ii=kk
                  piv=1./a(ii)
                  lp=0
                  do i=k,m1
                        ll=ii
                        ii=ii+i
                        ap=a(ii)*piv
                        lp=lp+1
                        ij=ii-km1
                        do j=i,m1
                              ij=ij+j
                              ll=ll+j
                              a(ij)=a(ij)-ap*a(ll)
                        end do
                        do llb=k,nm,neq
                              in=llb+lp+nm1
                              ll1=llb+nm1
                              r(in)=r(in)-ap*r(ll1)
                        end do
                  end do
c            end do
      end do
c
c Error checking - singular matrix:
c
c      ijm=ij-nn*(nright-1)
      ijm=ij
      if(abs(a(ijm)).lt.tol) then
            ising=neq
            return
      endif
c
c Finished triangulation, start solving back:
c
c      do iv=1,nright
c            nm1=nm*(iv-1)
            nm1=0
c            ij=ijm+nn*(iv-1)
            ij=ijm
            piv=1./a(ij)
            do llb=neq,nm,neq
                  ll1=llb+nm1
                  s(ll1)=r(ll1)*piv
            end do
            i=neq
            kk=ij
            do ii=1,m1
                  kk=kk-i
                  piv=1./a(kk)
                  i=i-1
                  do llb=i,nm,neq
                        ll1=llb+nm1
                        in=ll1
                        ap=r(in)
                        ij=kk
                        do j=i,m1
                              ij=ij+j
                              in=in+1
                              ap=ap-a(ij)*s(in)
                        end do
                        s(ll1)=ap*piv
                  end do
            end do
c      end do
c
c Finished solving back, return:
c
      return
      end



      subroutine ksolunroll24(nright,neq,nsb,a,r,s,ising)
c-----------------------------------------------------------------------
c
c                Solution of a System of Linear Equations
c                ****************************************
c
c
c
c INPUT VARIABLES:
c
c   nright,nsb       number of columns in right hand side matrix.
c                      for KB2D: nright=1, nsb=1
c   neq              number of equations
c   a()              upper triangular left hand side matrix (stored 
c                      columnwise)
c   r()              right hand side matrix (stored columnwise)
c                      for kb2d, one column per variable
c
c
c
c OUTPUT VARIABLES:
c
c   s()              solution array, same dimension as  r  above.
c   ising            singularity indicator
c                      0,  no singularity problem
c                     -1,  neq .le. 1
c                      k,  a null pivot appeared at the kth iteration
c
c
c
c PROGRAM NOTES:
c
c   1. Requires the upper triangular left hand side matrix.
c   2. Pivots are on the diagonal.
c   3. Does not search for max. element for pivot.
c   4. Several right hand side matrices possible.
c   5. USE for ok and sk only, NOT for UK.
c
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8   a(*),r(*),s(*)
c
c If there is only one equation then set ising and return:
c
c      if(neq.le.1) then
c            ising = -1
c            return
c      endif
c
c Initialize:
c
      tol   = 0.1e-06
      ising = 0
      nn    = 12*(24+1)
      nm    = 24
      m1    = 24-1
      kk    = 0
c
c Start triangulation:
c
      do k=1,22,2
            kk=kk+k
            ak=a(kk)
            if(abs(ak).lt.tol) then
                  ising=k
                  return
            endif
            km1=k-1
            ii=kk
            piv=1./a(ii)
            lp=0
            do i=k,23
                  ll=ii
                  ii=ii+i
                  ap=a(ii)*piv
                  lp=lp+1
                  ij=ii-km1
                  do j=i,23
                        ij=ij+j
                        ll=ll+j
                        a(ij)=a(ij)-ap*a(ll)
                  end do
                  r(k+lp)=r(k+lp)-ap*r(k)
            end do

            kk=kk+k+1
            ak=a(kk)
            if(abs(ak).lt.tol) then
                  ising=k+1
                  return
            endif
            km1=k-1+1
            ii=kk
            piv=1./a(ii)
            lp=0
            do i=(k+1),23
                  ll=ii
                  ii=ii+i
                  ap=a(ii)*piv
                  lp=lp+1
                  ij=ii-km1
                  do j=i,23
                        ij=ij+j
                        ll=ll+j
                        a(ij)=a(ij)-ap*a(ll)
                  end do
                  r(k+1+lp)=r(k+1+lp)-ap*r(k+1)
            end do

c            kk=kk+k+2
c            ak=a(kk)
c            if(abs(ak).lt.tol) then
c                  ising=k+2
c                  return
c            endif
c            km1=k-1+2
c            ii=kk
c            piv=1./a(ii)
c            lp=0
c            do i=(k+2),23
c                  ll=ii
c                  ii=ii+i
c                  ap=a(ii)*piv
c                  lp=lp+1
c                  ij=ii-km1
c                  do j=i,23
c                        ij=ij+j
c                        ll=ll+j
c                        a(ij)=a(ij)-ap*a(ll)
c                  end do
c                  r(k+2+lp)=r(k+2+lp)-ap*r(k+2)
c            end do
c
c            kk=kk+k+3
c            ak=a(kk)
c            if(abs(ak).lt.tol) then
c                  ising=k+3
c                  return
c            endif
c            km1=k-1+3
c            ii=kk
c            piv=1./a(ii)
c            lp=0
c            do i=(k+3),23
c                  ll=ii
c                  ii=ii+i
c                  ap=a(ii)*piv
c                  lp=lp+1
c                  ij=ii-km1
c                  do j=i,23
c                        ij=ij+j
c                        ll=ll+j
c                        a(ij)=a(ij)-ap*a(ll)
c                  end do
c                  r(k+3+lp)=r(k+3+lp)-ap*r(k+3)
c            end do


      end do

c      kk=kk+21
c      ak=a(kk)
c      if(abs(ak).lt.tol) then
c            ising=k+21
c            return
c      endif
c      km1=21-1
c      ii=kk
c      piv=1./a(ii)
c      lp=0
c      do i=21,23
c            ll=ii
c            ii=ii+i
c            ap=a(ii)*piv
c            lp=lp+1
c            ij=ii-km1
c            do j=i,23
c                  ij=ij+j
c                  ll=ll+j
c                  a(ij)=a(ij)-ap*a(ll)
c            end do
c            r(21+lp)=r(21+lp)-ap*r(21)
c      end do
c
c
c      kk=kk+22
c      ak=a(kk)
c      if(abs(ak).lt.tol) then
c            ising=k+22
c            return
c      endif
c      km1=22-1
c      ii=kk
c      piv=1./a(ii)
c      lp=0
c      do i=22,23
c            ll=ii
c            ii=ii+i
c            ap=a(ii)*piv
c            lp=lp+1
c            ij=ii-km1
c            do j=i,23
c                  ij=ij+j
c                  ll=ll+j
c                  a(ij)=a(ij)-ap*a(ll)
c            end do
c            r(22+lp)=r(22+lp)-ap*r(22)
c      end do



      kk=kk+23
      ak=a(kk)
      if(abs(ak).lt.tol) then
            ising=23
            return
      endif
      km1=23-1
      ii=kk
      piv=1./a(ii)
      lp=0
      i=23
            ll=ii
            ii=ii+i
            ap=a(ii)*piv
            lp=lp+1
            ij=ii-km1
            j=23
                  ij=ij+j
                  ll=ll+j
                  a(ij)=a(ij)-ap*a(ll)
            
            r(23+lp)=r(23+lp)-ap*r(23)
      




c
c Error checking - singular matrix:
c
c      ijm=ij-nn*(nright-1)
      ijm=ij
      if(abs(a(ijm)).lt.tol) then
            ising=neq
            return
      endif
c
c Finished triangulation, start solving back:
c
c      do iv=1,nright
c            nm1=nm*(iv-1)
            nm1=0
c            ij=ijm+nn*(iv-1)
            ij=ijm
            piv=1./a(ij)
            do llb=neq,nm,neq
                  ll1=llb+nm1
                  s(ll1)=r(ll1)*piv
            end do
            i=neq
            kk=ij
            do ii=1,m1
                  kk=kk-i
                  piv=1./a(kk)
                  i=i-1
                  do llb=i,nm,neq
                        ll1=llb+nm1
                        in=ll1
                        ap=r(in)
                        ij=kk
                        do j=i,m1
                              ij=ij+j
                              in=in+1
                              ap=ap-a(ij)*s(in)
                        end do
                        s(ll1)=ap*piv
                  end do
            end do
c      end do
c
c Finished solving back, return:
c
      return
      end


      subroutine ksolO1(nright,neq,nsb,a,r,s,ising)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8   a(*),r(*),s(*)
c
c If there is only one equation then set ising and return:
c
      if(neq.le.1) then
            ising = -1
            return
      endif
c
c Initialize:
c
      tol   = 0.1e-06
      ising = 0
      nn    = neq*(neq+1)/2
c      nm    = nsb*neq
      nm    = neq
      m1    = neq-1
      kk    = 0
c
c Start triangulation:
c
      do k=1,m1
            kk=kk+k
            ak=a(kk)
            if(abs(ak).lt.tol) then
                  ising=k
                  return
            endif
            km1=k-1
            nm1=0
            ii=kk
            piv=1./a(ii)
            lp=0
c            do i=k,m1,2
cc first iteration
c                  ll=ii
c                  ii=ii+i
c                  ap=a(ii)*piv
c                  lp=lp+1
c                  ij=ii-km1
c                  do j=i,m1
c                        ij=ij+j
c                        ll=ll+j
c                        a(ij)=a(ij)-ap*a(ll)
c                  end do
c                  do llb=k,nm,neq
c                        in=llb+lp+nm1
c                        ll1=llb+nm1
c                        r(in)=r(in)-ap*r(ll1)
c                  end do
cc second iteration
c                  ll=ii
c                  ii=ii+i+1
c                  ap=a(ii)*piv
c                  lp=lp+1
c                  ij=ii-km1
c                  do j=(i+1),m1
c                        ij=ij+j
c                        ll=ll+j
c                        a(ij)=a(ij)-ap*a(ll)
c                  end do
c                  do llb=k,nm,neq
c                        in=llb+lp+nm1
c                        ll1=llb+nm1
c                        r(in)=r(in)-ap*r(ll1)
c                  end do
c            end do
c last iteration
c            ilast=i
            do i=k,m1
                  ll=ii
                  ii=ii+i
                  ap=a(ii)*piv
                  lp=lp+1
                  ij=ii-km1
                  do j=i,m1
                        ij=ij+j
                        ll=ll+j
                        a(ij)=a(ij)-ap*a(ll)
                  end do
                  do llb=k,nm,neq
                        in=llb+lp+nm1
                        ll1=llb+nm1
                        r(in)=r(in)-ap*r(ll1)
                  end do
            end do


      end do
c
c Error checking - singular matrix:
c
c      ijm=ij-nn*(nright-1)
      ijm=ij
      if(abs(a(ijm)).lt.tol) then
            ising=neq
            return
      endif
c
c Finished triangulation, start solving back:
c
c      do iv=1,nright
       iv=1

            nm1=0
            ij=ijm
            piv=1./a(ij)
            do llb=neq,nm,neq
                  ll1=llb+nm1
                  s(ll1)=r(ll1)*piv
            end do
            i=neq
            kk=ij
            do ii=1,m1
                  kk=kk-i
                  piv=1./a(kk)
                  i=i-1
                  do llb=i,nm,neq
                        ll1=llb+nm1
                        in=ll1
                        ap=r(in)
                        ij=kk
                        do j=i,m1
                              ij=ij+j
                              in=in+1
                              ap=ap-a(ij)*s(in)
                        end do
                        s(ll1)=ap*piv
                  end do
            end do
c      end do
c
c Finished solving back, return:
c
      return
      end


