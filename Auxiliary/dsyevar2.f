      SUBROUTINE dsyevar2(m1, AD, IndexD, m2, AOD, IndexODi,
     &                    IndexODj, tol, maxitr, n, maxn, nev,
     &                    maxnev, ncv, maxncv, itrdone,
     &                    matvecmultdone, infoaupd, infoeupd, Evals,
     &                    Evecs)
c********************************************************************
c     The subroutine dsyevar2 computes the nev lowest lying
c     eigenstates of a symmetric nxn matrix of type double precision.
c     This code is based on dssimp.f from arpack96.tar.gz.
c     The m1 non-zero elements of the diagonal of the symmetric
c     matrix are stored in AD and the m2 non-zero elements of the
c     upper half are stored in AOD.
c     The kth element in AD, AD(k), is originally located at row and
c     column number IndexD(k), and the kth element in AOD, AOD(k), is
c     originally located at row IndexODi(k) and column IndexODj(k).
c     IndexD, IndexODi, and IndexODj are assumed to be 0-based, since
c     this subroutine is frequently called from a C++ program.
c     Matrix-vector multiplication is done by the subroutine dmyav
c     defined below. This subroutine makes use of OpenMP.
c
c     Parameters:
c     m1               int, input
c                      length of AD
c     AD(m1)           double precision, input, unchanged
c                      linear array of the non-zero diagonal elements
c     IndexD(m1)       int, input
c                      the element AD(k) is located in the original
c                      matrix at row and column number IndexD(k),
c                      0-based index
c     m2               int, input
c                      length of AOD
c     AOD(m2)          double precision, input, unchanged
c                      linear array of the non-zero offdiagonal
c                      elements of the upper half
c     IndexODi(m2)     int, input
c                      the element AOD(k) is located in the original
c                      matrix at row number IndexODi(k),
c                      0-based index
c     IndexODj(m2)     int, input
c                      the element AOD(k) is located in the original
c                      matrix at column number IndexODj(k),
c                      0-based index
c     tol              double precision, input
c                      convergence criterion, if it is 0.0 macheps is
c                      used as tol
c     maxitr           int, input
c                      maximal number of Arnoldi iterations allowed,
c                      typical value is 300
c     n                int, input
c                      dimension of the matrix under consideration
c                      whose non-zero elements are in A
c     maxn             int, input
c                      maximal n, can be set equal to n
c     nev              int, input
c                      number of eigenvalues searched
c     maxnev           int, input
c                      maximal number of eigenvalues searched, can be
c                      set equal to nev
c     ncv              int, input
c                      number of Lanczos vectors computed, must be
c                      greater or equal to nev+1, typical value for
c                      computing only the ground state is 25
c     maxncv           int, input
c                      maximal ncv, can be set equal to ncv
c     itrdone          int, output
c                      resulting number of Arnoldi iterations done
c     matvecmultdone   int, output
c                      resulting number of matrix-vector
c                      multiplications done
c     infoaupd         int, output
c                      info of dsaupd, must be 0 on successful exit
c     infoeupd         int, output
c                      info of dseupd, must be 0 on successful exit
c     Evals(nev)       double precision, output
c                      the nev lowest lying eigenvalues
c     Evecs(n,nev)     double precision, output
c                      the nev lowest lying eigenstates
c
c     Caution:
c     The ARPACK subroutine dsaupd has problems with tol!
c
c     This software is part of the PEPS Library.
c     This software was programmed by Michael Lubasch when he was a
c     PhD student in the Theory Division of the Max Planck Institute
c     of Quantum Optics where his work was supervised by
c     Prof. Dr. J. Ignacio Cirac and Dr. Mari-Carmen Ba\~{n}uls.
c     This software is distributed under the
c     PEPS Library License v1.0 (see accompanying file LICENSE.txt).
c********************************************************************
      INTEGER m1, IndexD(m1), m2, IndexODi(m2), IndexODj(m2), maxitr
      INTEGER n, maxn, nev, maxnev, ncv, maxncv, itrdone
      INTEGER matvecmultdone, infoaupd, infoeupd
      DOUBLE PRECISION AD(m1), AOD(m2), tol, Evals(nev), Evecs(n,nev)
      LOGICAL select(maxncv), rvec
      CHARACTER bmat*1, which*2
      INTEGER ldv, iparam(11), ipntr(11), ido, lworkl, info, ierr, i
      INTEGER j, nx, ishfts, mode1, nconv
      DOUBLE PRECISION v(maxn,maxncv), workl(maxncv*(maxncv+8))
      DOUBLE PRECISION workd(3*maxn), d(maxncv,2), resid(maxn)
      DOUBLE PRECISION ax(maxn), sigma, zero, dnrm2
      external dnrm2, daxpy
      intrinsic abs
      zero=0.0d0
      ldv=maxn
      bmat='I'
      which='SA'
      lworkl=ncv*(ncv+8)
      info=0
      ido=0
      ishfts=1
      mode1=1
      iparam(1)=ishfts
      iparam(3)=maxitr
      iparam(7)=mode1
 10   continue
       call dsaupd(ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,
     &             ipntr,workd,workl,lworkl,info)
       if ((ido .eq. -1) .or. (ido .eq. 1)) then
        call dmyav(m1,AD,IndexD,m2,AOD,IndexODi,IndexODj,n,
     &            workd(ipntr(1)),workd(ipntr(2)))
        go to 10
       end if
      infoaupd=info
      if (info .lt. 0) then
       print *, ' '
       print *, ' Error in ARPACK routine dsyevar2: '
       print *, ' Error with _saupd, info = ', info
       print *, ' Check documentation in _saupd '
       print *, ' '
      else
       rvec=.true.
       call dseupd(rvec,'All',select,d,v,ldv,sigma,bmat,n,which,nev,
     &             tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,
     &             lworkl,ierr)
       infoeupd=ierr
       if (ierr .ne. 0) then
        print *, ' '
        print *, ' Error in ARPACK routine dsyevar2: '
        print *, ' Error with _seupd, info = ', ierr
        print *, ' Check the documentation of _seupd. '
        print *, ' '
       else
        nconv=iparam(5)
        do 20 j=1, nconv
         call dmyav(m1,AD,IndexD,m2,AOD,IndexODi,IndexODj,n,v(1,j),ax)
         call daxpy(n,-d(j,1),v(1,j),1,ax,1)
         d(j,2)=dnrm2(n,ax,1)
         d(j,2)=d(j,2)/abs(d(j,1))
 20     continue
        do 30 j=1, nconv
         Evals(j)=d(j,1)
         do 40 i=1, n
          Evecs(i,j)=v(i,j)
 40      continue
 30     continue
       end if
       itrdone=iparam(3)
       matvecmultdone=iparam(9)
      end if
      END

      SUBROUTINE dmyav(m1,AD,IndexD,m2,AOD,IndexODi,IndexODj,n,v,w)
      INTEGER m1, IndexD(m1), m2, IndexODi(m2), IndexODj(m2), n
      DOUBLE PRECISION AD(m1), AOD(m2), v(n), w(n)
      INTEGER i, j, k
c     Keep in mind that all indices are 0-based!
      !$OMP PARALLEL PRIVATE(i, j)
      !$OMP DO
      do 11 k=1, n
       w(k)=0.0d0
 11   continue
      !$OMP END DO
      !$OMP DO
      do 12 k=1, m1
       i=IndexD(k)+1
       w(i)=AD(k)*v(i)
 12   continue
      !$OMP END DO
      !$OMP DO
      do 13 k=1, m2
       i=IndexODi(k)+1
       j=IndexODj(k)+1
       w(i)=w(i)+AOD(k)*v(j)
       w(j)=w(j)+AOD(k)*v(i)
 13   continue
      !$OMP END DO
      !$OMP END PARALLEL
      END
