/// The required BLAS, LAPACK, and ARPACK functions.
/** This header file contains the required BLAS, LAPACK, and ARPACK functions.
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

extern "C"
{
// Various elementary routines:
 float snrm2_(int* n, float* X, int* incx);
 double dnrm2_(int* n, double* X, int* incx);
 void saxpy_(int* n, float* sa, float* Sx, int* incx, float* Sy, int* incy);
 void daxpy_(int* n, double* da, double* Dx, int* incx, double* Dy, int* incy);
// Random number generators:
 void slarnv_(int* idist, int* Iseed, int* n, float* X);
 void dlarnv_(int* idist, int* Iseed, int* n, double* X);
 void clarnv_(int* idist, int* Iseed, int* n, complex<float>* X);
 void zlarnv_(int* idist, int* Iseed, int* n, complex<double>* X);
// Matrix-vector multiplication:
 void sgemv_(char* trans, int* m, int* n, float* alpha, float* A, int* lda, float* X, int* incx,
             float* beta, float* Y, int* incy);
 void dgemv_(char* trans, int* m, int* n, double* alpha, double* A, int* lda, double* X, int* incx,
             double* beta, double* Y, int* incy);
 void cgemv_(char* trans, int* m, int* n, complex<float>* alpha, complex<float>* A, int* lda,
             complex<float>* X, int* incx, complex<float>* beta, complex<float>* Y, int* incy);
 void zgemv_(char* trans, int* m, int* n, complex<double>* alpha, complex<double>* A, int* lda,
             complex<double>* X, int* incx, complex<double>* beta, complex<double>* Y, int* incy);
// Matrix-matrix multiplication:
 void sgemm_(char* transa, char* transb, int* m, int* n, int* k, float* alpha, float* A,
             int* lda, float* B, int* ldb, float* beta, float* C, int* ldc);
 void dgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* A,
             int* lda, double* B, int* ldb, double* beta, double* C, int* ldc);
 void cgemm_(char* transa, char* transb, int* m, int* n, int* k, complex<float>* alpha,
             complex<float>* A, int* lda, complex<float>* B, int* ldb, complex<float>* beta,
             complex<float>* C, int* ldc);
 void zgemm_(char* transa, char* transb, int* m, int* n, int* k, complex<double>* alpha,
             complex<double>* A, int* lda, complex<double>* B, int* ldb, complex<double>* beta,
             complex<double>* C, int* ldc);
// QR decomposition:
 void sgeqrf_(int* m, int* n, float* A, int* lda, float* Tau, float* Work, int* lwork, int* info);
 void dgeqrf_(int* m, int* n, double* A, int* lda, double* Tau, double* Work, int* lwork, int* info);
 void cgeqrf_(int* m, int* n, complex<float>* A, int* lda, complex<float>* Tau, complex<float>* Work, int* lwork, int* info);
 void zgeqrf_(int* m, int* n, complex<double>* A, int* lda, complex<double>* Tau, complex<double>* Work, int* lwork, int* info);
// compute Q after QR decomposition with xgeqrf:
 void sormqr_(char* side, char* trans, int* m, int* n, int* k, float* A, int* lda, float* Tau, float* C, int* ldc, float* Work, int* lwork, int* info);
 void dormqr_(char* side, char* trans, int* m, int* n, int* k, double* A, int* lda, double* Tau, double* C, int* ldc, double* Work, int* lwork, int* info);
 void cunmqr_(char* side, char* trans, int* m, int* n, int* k, complex<float>* A, int* lda, complex<float>* Tau, complex<float>* C, int* ldc, complex<float>* Work, int* lwork, int* info);
 void zunmqr_(char* side, char* trans, int* m, int* n, int* k, complex<double>* A, int* lda, complex<double>* Tau, complex<double>* C, int* ldc, complex<double>* Work, int* lwork, int* info);
// Singular value decomposition:
 void sgesvd_(char* jobu, char* jobvt, int* m, int* n, float* A, int* lda, float* S, float* U,
              int* ldu, float* Vt, int* ldvt, float* Work, int* lwork, int* info);
 void dgesvd_(char* jobu, char* jobvt, int* m, int* n, double* A, int* lda, double* S,
              double* U, int* ldu, double* Vt, int* ldvt, double* Work, int* lwork, int* info);
 void cgesvd_(char* jobu, char* jobvt, int* m, int* n, complex<float>* A, int* lda,
              float* S, complex<float>* U, int* ldu, complex<float>* Vt, int* ldvt,
              complex<float>* Work, int* lwork, float* Rwork, int* info);
 void zgesvd_(char* jobu, char* jobvt, int* m, int* n, complex<double>* A, int* lda,
              double* S, complex<double>* U, int* ldu, complex<double>* Vt, int* ldvt,
              complex<double>* Work, int* lwork, double* Rwork, int* info);
// General eigenvalue problem:
 void sgeev_(char* jobvl, char* jobvr, int* n, float* A, int* lda, float* Wr, float* Wi, float* Vl,
             int* ldvl, float* Vr, int* ldvr, float* Work, int* lwork, int* info);
 void dgeev_(char* jobvl, char* jobvr, int* n, double* A, int* lda, double* Wr, double* Wi, double* Vl,
             int* ldvl, double* Vr, int* ldvr, double* Work, int* lwork, int* info);
 void cgeev_(char* jobvl, char* jobvr, int* n, complex<float>* A, int* lda, complex<float>* W,
             complex<float>* Vl, int* ldvl, complex<float>* Vr, int* ldvr, complex<float>* Work,
             int* lwork, float* Rwork, int* info);
 void zgeev_(char* jobvl, char* jobvr, int* n, complex<double>* A, int* lda, complex<double>* W,
             complex<double>* Vl, int* ldvl, complex<double>* Vr, int* ldvr, complex<double>* Work,
             int* lwork, double* Rwork, int* info);
// Hermitian eigenvalue problem:
 void ssyev_(char* jobz, char* uplo, int* n, float* A, int* lda, float* W, float* Work, int* lwork,
             int* info);
 void dsyev_(char* jobz, char* uplo, int* n, double* A, int* lda, double* W, double* Work, int* lwork,
             int* info);
 void cheev_(char* jobz, char* uplo, int* n, complex<float>* A, int* lda, float* W,
             complex<float>* Work, int* lwork, float* Rwork, int* info);
 void zheev_(char* jobz, char* uplo, int* n, complex<double>* A, int* lda, double* W,
             complex<double>* Work, int* lwork, double* Rwork, int* info);
// Linear least squares using divide-and-conquer SVD:
 void sgelsd_(int* m, int* n, int* nrhs, float* A, int* lda, float* B, int* ldb, float* S, float* rcond,
              int* rank, float* Work, int* lwork, int* Iwork, int* info);
 void dgelsd_(int* m, int* n, int* nrhs, double* A, int* lda, double* B, int* ldb, double* S, double* rcond,
              int* rank, double* Work, int* lwork, int* Iwork, int* info);
 void cgelsd_(int* m, int* n, int* nrhs, complex<float>* A, int* lda, complex<float>* B, int* ldb,
              float* S, float* rcond, int* rank, complex<float>* Work, int* lwork, float* Rwork,
              int* Iwork, int* info);
 void zgelsd_(int* m, int* n, int* nrhs, complex<double>* A, int* lda, complex<double>* B, int* ldb,
              double* S, double* rcond, int* rank, complex<double>* Work, int* lwork, double* Rwork,
              int* Iwork, int* info);
// Generalized symmetric definite eigenproblem using divide-and-conquer:
 void ssygvd_(int* itype, char* jobz, char* uplo, int* n, float* A, int* lda, float* B, int* ldb, float* W,
              float* Work, int* lwork, int* Iwork, int* liwork, int* info);
 void dsygvd_(int* itype, char* jobz, char* uplo, int* n, double* A, int* lda, double* B, int* ldb, double* W,
              double* Work, int* lwork, int* Iwork, int* liwork, int* info);
 void chegvd_(int* itype, char* jobz, char* uplo, int* n, complex<float>* A, int* lda, complex<float>* B, int* ldb,
              float* W, complex<float>* Work, int* lwork, float* Rwork, int* lrwork, int* Iwork, int* liwork, int* info);
 void zhegvd_(int* itype, char* jobz, char* uplo, int* n, complex<double>* A, int* lda, complex<double>* B, int* ldb,
              double* W, complex<double>* Work, int* lwork, double* Rwork, int* lrwork, int* Iwork, int* liwork, int* info);
// Generalized nonsymmetric eigenproblem:
// logicals are passed as integers which must be either 0 for false or 1 for true:
 void sggevx_(char* balanc, char* jobvl, char* jobvr, char* sense, int* n, float* A, int* lda, float* B, int* ldb,
              float* Alphar, float* Alphai, float* Beta, float* Vl, int* ldvl, float* Vr, int* ldvr, int* ilo,
              int* ihi, float* Lscale, float* Rscale, float* abnrm, float* bbnrm, float* Rconde, float* Rcondv,
              float* Work, int* lwork, int* Iwork, int* Bwork, int* info);
 void dggevx_(char* balanc, char* jobvl, char* jobvr, char* sense, int* n, double* A, int* lda, double* B, int* ldb,
              double* Alphar, double* Alphai, double* Beta, double* Vl, int* ldvl, double* Vr, int* ldvr, int* ilo,
              int* ihi, double* Lscale, double* Rscale, double* abnrm, double* bbnrm, double* Rconde, double* Rcondv,
              double* Work, int* lwork, int* Iwork, int* Bwork, int* info);
 void cggevx_(char* balanc, char* jobvl, char* jobvr, char* sense, int* n, complex<float>* A, int* lda,
              complex<float>* B, int* ldb, complex<float>* Alpha, complex<float>* Beta, complex<float>* Vl,
              int* ldvl, complex<float>* Vr, int* ldvr, int* ilo, int* ihi, float* Lscale, float* Rscale,
              float* abnrm, float* bbnrm, float* Rconde, float* Rcondv, complex<float>* Work, int* lwork,
              float* Rwork, int* Iwork, int* Bwork, int* info);
 void zggevx_(char* balanc, char* jobvl, char* jobvr, char* sense, int* n, complex<double>* A, int* lda,
              complex<double>* B, int* ldb, complex<double>* Alpha, complex<double>* Beta, complex<double>* Vl,
              int* ldvl, complex<double>* Vr, int* ldvr, int* ilo, int* ihi, double* Lscale, double* Rscale,
              double* abnrm, double* bbnrm, double* Rconde, double* Rcondv, complex<double>* Work, int* lwork,
              double* Rwork, int* Iwork, int* Bwork, int* info);
// ARPACK routines:
 void ssaupd_(int* ido, char* bmat, int* n, char* Which, int* nev, float* tol, float* Resid, int* ncv,
              float* V, int* ldv, int* Iparam, int* Ipntr, float* Workd, float* Workl, int* lworkl,
              int* info);
 void dsaupd_(int* ido, char* bmat, int* n, char* Which, int* nev, double* tol, double* Resid, int* ncv,
              double* V, int* ldv, int* Iparam, int* Ipntr, double* Workd, double* Workl, int* lworkl,
              int* info);
// logicals are passed as integers which must be either 0 for false or 1 for true:
 void sseupd_(int* rvec, char* howmny, int* Select, float* D, float* Z, int* ldz, float* sigma, char* bmat,
              int* n, char* Which, int* nev, float* tol, float* Resid, int* ncv, float* V, int* ldv,
              int* Iparam, int* Ipntr, float* Workd, float* Workl, int* lworkl, int* info);
// logicals are passed as integers which must be either 0 for false or 1 for true:
 void dseupd_(int* rvec, char* howmny, int* Select, double* D, double* Z, int* ldz, double* sigma,
              char* bmat, int* n, char* Which, int* nev, double* tol, double* Resid, int* ncv, double* V,
              int* ldv, int* Iparam, int* Ipntr, double* Workd, double* Workl, int* lworkl, int* info);
// Ground state computation:
 void ssyevar_(int* n, float* A, float* tol, int* maxitr, int* maxn, int* maxnev, int* ncv,
               int* maxncv, int* itrdone, int* matvecmultdone, int* infoaupd, int* infoeupd,
               float* eval, float* evec);
 void dsyevar_(int* n, double* A, double* tol, int* maxitr, int* maxn, int* maxnev, int* ncv,
               int* maxncv, int* itrdone, int* matvecmultdone, int* infoaupd, int* infoeupd,
               double* eval, double* evec);
// Lowest lying eigenstates:
 void ssyevar2_(int* m1, float* AD, int* IndexD, int* m2, float* AOD, int* IndexODi, int* IndexODj,
                float* tol, int* maxitr, int* n, int* maxn, int* nev, int* maxnev, int* ncv, int* maxncv,
                int* itrdone, int* matvecmultdone, int* infoaupd, int* infoeupd, float* Evals, float* Evecs);
 void dsyevar2_(int* m1, double* AD, int* IndexD, int* m2, double* AOD, int* IndexODi, int* IndexODj,
                double* tol, int* maxitr, int* n, int* maxn, int* nev, int* maxnev, int* ncv, int* maxncv,
                int* itrdone, int* matvecmultdone, int* infoaupd, int* infoeupd, double* Evals, double* Evecs);
// Runge-Kutta fourth order:
 void zrk4_(complex<double>* Y, complex<double>* H0, complex<double>* H1, complex<double>* H2,
            int* n, double* h, complex<double>* Yout);
};
