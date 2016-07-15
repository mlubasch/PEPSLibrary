      SUBROUTINE zrk4(y,h0,h1,h2,n,h,yout)
c********************************************************************
c     The subroutine zrk4 propagates an initial state y at time x to
c     a final state yout at time x+h with a possibly time-dependent
c     Hamiltonian H(t). The subroutine needs derivs(n,h,y,dydx)
c     defined below.
c     Given values for the variables y(1:n) use the fourth-order
c     Runge-Kutta method to advance the solution over an interval h
c     and return the incremented variables as yout(1:n), which need
c     not be a distinct array from y.
c
c     Parameters:
c     y(n)             complex*16, input
c                      initial state at x
c     h0(n,n)          complex*16, input
c                      Hamiltonian at time x
c     h1(n,n)          complex*16, input
c                      Hamiltonian at time x+0.5*h
c     h2(n,n)          complex*16, input
c                      Hamiltonian at time x+h
c     n                int, input
c                      dimension of y
c     h                double precision, input
c                      time step over which y is propagated
c     yout(n)          complex*16, output
c                      the resulting y propagated over time step h
c
c     This software is part of the PEPS Library.
c     This software was programmed by Michael Lubasch when he was a
c     PhD student in the Theory Division of the Max Planck Institute
c     of Quantum Optics where his work was supervised by
c     Prof. Dr. J. Ignacio Cirac and Dr. Mari-Carmen Ba\~{n}uls.
c     This software is distributed under the
c     PEPS Library License v1.0 (see accompanying file LICENSE.txt).
c********************************************************************
      INTEGER n
      DOUBLE PRECISION h
      COMPLEX*16 h0(n,n),h1(n,n),h2(n,n),y(n),yout(n)
      INTEGER i
      DOUBLE PRECISION h6,hh
      COMPLEX*16 dydx(n),dym(n),dyt(n),yt(n)
      call derivs(n,h0,y,dydx)
      hh=h*0.5d0
      h6=h/6.0d0
      do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
11    continue
      call derivs(n,h1,yt,dyt)
      do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
12    continue
      call derivs(n,h1,yt,dym)
      do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue
      call derivs(n,h2,yt,dyt)
      do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.0d0*dym(i))
14    continue
      return
      END

      SUBROUTINE derivs(n,h,y,dydx)
c********************************************************************
c     The subroutine derivs computes the time derivative of the state
c     vector y at time x. It is given by the Schroedinger equation
c     dydx = -i*H*y where H is the Hamiltonian at time x. The
c     subroutine uses LAPACK's ZSYMV for the matrix-vector
c     multiplication.
c
c     Parameters:
c     n                int, input
c                      dimension of y
c     h(n,n)           complex*16, input
c                      Hamiltonian at time x
c     y(n)             complex*16, input
c                      state vector at time x
c     dydx(n)          complex*16, output
c                      the derivative of y at time x
c********************************************************************
      INTEGER n
      COMPLEX*16 dydx(n),h(n,n),y(n)
      COMPLEX*16 imagUnit
      imagUnit=dcmplx(0.0d0,1.0d0)
      call zsymv('U', n, -imagUnit, h, n, y, 1, 0.0d0, dydx, 1)
      return
      END
