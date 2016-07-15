      subroutine dsyevar(n, A, tol, maxitr, maxn, maxnev, ncv,
     &                   maxncv, itrdone, matvecmultdone, infoaupd,
     &                   infoeupd, eval, evec)
c********************************************************************
c     The subroutine dsyevar computes the ground state of a symmetric
c     matrix of type double precision. This code is based on dssimp.f
c     from arpack96.tar.gz. The matrix-vector multiplication is done
c     by LAPACK's routine dsymv and hence makes no use of sparseness.
c
c     Parameters:
c     n                int, input
c                      dimension of A
c     A(n,n)           double precision, input, unchanged
c                      symmetric matrix of type double precision and
c                      dimension n
c     tol              double precision, input
c                      convergence criterion, if it is 0.0 macheps is
c                      used as tol
c     maxitr           int, input
c                      maximal number of Arnoldi iterations allowed,
c                      standard value is 300
c     maxn             int, input
c                      maximal n, can be set equal to n
c     maxnev           int, input
c                      maximal number of eigenvalues searched, can be
c                      set equal to 1
c     ncv              int, input
c                      number of Lanczos vectors computed, should be
c                      greater or equal to 2, standard value is 25
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
c     eval             double precision, output
c                      ground state energy
c     evec(n)          double precision, output
c                      ground state
c
c     This software is part of the PEPS Library.
c     This software was programmed by Michael Lubasch when he was a
c     PhD student in the Theory Division of the Max Planck Institute
c     of Quantum Optics where his work was supervised by
c     Prof. Dr. J. Ignacio Cirac and Dr. Mari-Carmen Ba\~{n}uls.
c     This software is distributed under the
c     PEPS Library License v1.0 (see accompanying file LICENSE.txt).
c********************************************************************
      integer          n, itrdone, matvecmultdone, infoaupd
      integer          infoeupd
      Double precision
     &                 A(n,n), eval, evec(n)
      integer          maxn, maxnev, maxncv, ldv
      Double precision
     &                 v(maxn,maxncv), workl(maxncv*(maxncv+8)),
     &                 workd(3*maxn), d(maxncv,2), resid(maxn),
     &                 ax(maxn)
      logical          select(maxncv)
      integer          iparam(11), ipntr(11)
      character        bmat*1, which*2
      integer          ido, nev, ncv, lworkl, info, ierr,
     &                 j, nx, ishfts, maxitr, mode1, nconv
      logical          rvec
      Double precision      
     &                 tol, sigma
      Double precision
     &                 zero
      parameter        (zero = 0.0D+0)
      Double precision           
     &                 dnrm2
      external         dnrm2, daxpy, dsymv
      intrinsic        abs
      ldv = maxn
      nev   = 1
      bmat  = 'I'
      which = 'SA'
      lworkl = ncv*(ncv+8)
      info = 0
      ido = 0
      ishfts = 1
      mode1 = 1
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode1
 10   continue
         call dsaupd ( ido, bmat, n, which, nev, tol, resid, 
     &                 ncv, v, ldv, iparam, ipntr, workd, workl,
     &                 lworkl, info )
         if (ido .eq. -1 .or. ido .eq. 1) then
            call dsymv('U', n, 1.0d0, A, n, workd(ipntr(1)),
     &                 1, 0.0d0, workd(ipntr(2)), 1)
            go to 10
         end if
      infoaupd = info
      if ( info .lt. 0 ) then
      else
          rvec = .true.
          call dseupd ( rvec, 'All', select, d, v, ldv, sigma, 
     &         bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &         iparam, ipntr, workd, workl, lworkl, ierr )
          infoeupd = ierr
          if ( ierr .ne. 0) then
          else
             nconv =  iparam(5)
             do 20 j=1, nconv
                call dsymv('U', n, 1.0d0, A, n, v(1,j), 1, 0.0d0,
     &                     ax, 1)
                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
 20          continue
             eval = d(1,1)
             do 30 i=1, n
              evec(i) = v(i,1)
 30          continue
          end if
      itrdone = iparam(3)
      matvecmultdone = iparam(9)
      end if
 9000 continue
      end subroutine dsyevar
