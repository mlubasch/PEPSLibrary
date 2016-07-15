/// Template class Matrix implements matrices.
/** The template class Matrix implements matrices and related operations by making use of BLAS and
    LAPACK.
    \param dim0 unsigned int, the dimension of the first index
    \param dim1 unsigned int, the dimension of the second index
    \param Type string, the type, must be "general" for a general matrix or "hermitian" for a hermitian
                matrix
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

// Declaration of template classes:
template<class T> class Tensor;
template<class T> class Matrix;
template<class T> class MPO;

// Declaration of friend functions:
template<class T> T expectationValueMatrix(const vector<T>& Vector0, const Matrix<T>& Matrix0);
template<class T> T scalarProductMatrixVector(const vector<T>& Vector0, const Matrix<T>& Matrix0,
                                              const vector<T>& Vector1);
template<class T> void getMPOFromLocalOperator(string BC, unsigned int N, unsigned int x,
                                               const Matrix<T>& O, MPO<T>& MPO0);
template<class T> void getMPOFromSumLocalOperator(string BC, unsigned int N, const Matrix<T>& O,
                                                  MPO<T>& MPO0);

template<class T> class Matrix : public Tensor<T>
{
 public:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  Matrix();

/// Constructor for Matrix with specific dimensions dim0 and dim1.
/** This constructor initializes a Matrix with a specific Shape.
    \param dim0 input: unsigned int, the dimension of the first index
    \param dim1 input: unsigned int, the dimension of the second index */
  Matrix(unsigned int dim0, unsigned int dim1);

/// Standard copy constructor.
/** The standard copy constructor copies the input Matrix into this.
    \param Matrix0 input: const Matrix<T>&, to be copied into this
    \sa Matrix<T>& operator=(const Matrix<T>& Matrix0) */
  Matrix(const Matrix<T>& Matrix0);

/// Standard destructor.
/** The standard destructor deletes the elements of the Matrix. */
  ~Matrix();

/// Assigns Matrix to this.
/** The operator= allows to assign a Matrix to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side Matrix.
    \param Matrix0 input: const Matrix<T>&, to be copied into this
    \return Matrix<T>&, a reference to the new this
    \sa Matrix(const Matrix<T>& Matrix0) */
  Matrix<T>& operator=(const Matrix<T>& Matrix0);

/// Returns dimension of first index of Matrix.
/** The returned dimension is the extent of the first index of the Matrix.
    \return unsigned int, the dimension of the first index of this Matrix */
  unsigned int getDim0() const { return this->Shape[0]; }

/// Returns dimension of second index of Matrix.
/** The returned dimension is the extent of the second index of this Matrix.
    \return unsigned int, the dimension of the second index of this Matrix */
  unsigned int getDim1() const { return this->Shape[1]; }

/// Allows to set and get elements from this Matrix by subscripts.
/** The operator() allows to set and get elements from this Matrix by two indices.
    \param index0 input: unsigned int, the row number
    \param index1 input: unsigned int, the column number
    \return T&, a reference to the element at the given position
    \sa inline void Tensor<T>::set(const vector<unsigned int>& Index, T element)
    \sa inline T Tensor<T>::get(const vector<unsigned int>& Index) const */
  inline T& operator()(unsigned int index0, unsigned int index1);

/// Allows to get elements from this Matrix by subscripts.
/** The operator() allows to get elements from this Matrix by two indices.
    \param index0 input: unsigned int, the row number
    \param index1 input: unsigned int, the column number
    \return const T&, a const reference to the element at the given position
    \sa inline void Tensor<T>::set(const vector<unsigned int>& Index, T element)
    \sa inline T Tensor<T>::get(const vector<unsigned int>& Index) const */
  inline const T& operator()(unsigned int index0, unsigned int index1) const;

/// Sets Type of Matrix.
/** This function sets the type of this Matrix.
    \param Type0 input: const string&, the type, must be "general" or "hermitian" */
  inline void setType(const string& Type0);

/// Returns Type of Matrix.
/** The returned type is either "general" or "hermitian".
    \param Type0 output: string&, the type, must be "general" or "hermitian" */
  inline void getType(string& Type0) const;

/// This is transposed.
/** This Matrix is transposed.
    \sa void transpose(Matrix<T>& Matrix0) const */
  void transpose();

/// Saves transpose of this Matrix in second Matrix.
/** The transpose of this Matrix is stored in the second Matrix0.
    \param Matrix0 output: Matrix<T>&, the resulting transpose
    \sa void transpose() */
  void transpose(Matrix<T>& Matrix0) const;

/// Adjoints this Matrix.
/** This Matrix is adjointed.
    \sa void Tensor<T>::complexConjugate()
    \sa void transpose()
    \sa void adjoint(Matrix<T>& Matrix0) const */
  void adjoint();

/// Returns adjoint of this Matrix.
/** The adjoint of this Matrix is returned as Matrix0.
    \param Matrix0 output: Matrix<T>&, the resulting adjoint of this Matrix
    \sa void Tensor<T>::complexConjugate()
    \sa void transpose()
    \sa void adjoint() */
  void adjoint(Matrix<T>& Matrix0) const;

/// Adds another Matrix to this and stores result in this.
/** Matrix0 is added to this and the result is stored in this.
    \param Matrix0 input: const Matrix<T>&, which is added to this, must have the same Shape as this
                          Matrix */
  void add(const Matrix<T>& Matrix0);

/// Multiplies this with scalar.
/** This is multiplied with a scalar x.
    \param x input: T, with which this Matrix is multiplied */
  void multiply(T x);

/// Multiplies this with vector to right.
/** This is multiplied with Vector0 to the right and the result is stored in Vector1. LAPACK's XGEMV
    is used.
    \param Vector0 input: const vector<T>&, with which this Matrix is right-multiplied
    \param Vector1 output: vector<T>&, the resulting vector, must have the correct form */
  void multiply(const vector<T>& Vector0, vector<T>& Vector1) const;

/// Multiplies this with second Matrix to right and stores result in this.
/** This is multiplied with Matrix0 to the right and the result is stored in this. LAPACK's XGEMM
    is used.
    \param Matrix0 input: const Matrix<T>&, with which this is right-multiplied
    \sa void multiply(const Matrix<T>& Matrix0, Matrix<T>& Matrix1) const */
  void multiply(const Matrix<T>& Matrix0);

/// Multiplies this with second Matrix to right and stores result in third Matrix.
/** This is multiplied with Matrix0 to the right and the result is stored in Matrix1. LAPACK's XGEMM
    is used.
    \param Matrix0 input: const Matrix<T>&, with which this is right-multiplied
    \param Matrix1 output: Matrix<T>&, the resulting Matrix, must be of dimension
                   this->getDim0() x Matrix0.getDim1()
    \sa void multiply(const Matrix<T>& Matrix0) */
  void multiply(const Matrix<T>& Matrix0, Matrix<T>& Matrix1) const;

/// Multiplies this with second Matrix to right as direct product and stores result in third Matrix.
/** This is multiplied with Matrix0 to the right as a direct product and the result is stored in
    Matrix1.
    \param Matrix0 input: const Matrix<T>&, with which this is right-multiplied as direct product
    \param Matrix1 output: Matrix<T>&, the resulting Matrix, on input it does not have to have a
                           specific shape, on output it will have a new shape and new elements */
  void multiplyDirectProduct(const Matrix<T>& Matrix0, Matrix<T>& Matrix1) const;

/// QR decomposes this Matrix A=Q*R.
/** The QR decomposition of this Matrix A=Q*R is computed, and Q is returned in this matrix and R in
    the input parameter.
    If dim0 >= dim1, then Q will be dim0 x dim1 and R dim1 x dim1, else Q will be dim0 x dim0 and
    R dim0 x dim1.
    LAPACK's XGEQRF and XORMQR/XUNMQR are used.
    \param R output: Matrix<T>&, the right matrix */
  void QRDecompose(Matrix<T>& R);

/// Singularvaluedecomposes this Matrix A=U*Sigma*Vt.
/** The singular value decomposition of this Matrix A=U*Sigma*Vt is computed and the results are
    stored in the input parameters. LAPACK's XGESVD is used.
   \param U output: Matrix<T>&, the left singular vectors, must be of dimension dim0 x dim0
   \param Sigma output: Matrix<T>&, the singular values, must be of dimension dim0 x dim1
   \param Vt output: Matrix<T>&, the right singular vectors, must be of dimension dim1 x dim1 */
  void singularValueDecompose(Matrix<T>& U, Matrix<T>& Sigma, Matrix<T>& Vt) const;

/// Pseudoinverts this Matrix.
/** Singularvaluedecomposing this Matrix and setting all singular values s_{i} to zero that fulfill
       s_{i}/s_{max} <= cutoff   ,
    the pseudoinverse of this Matrix is computed and returned in Matrix0.
    \param cutoff input: double, the cutoff, singular values s_{i} fulfilling s_{i}/s_{max} <= cutoff
                         are set to zero
    \param Matrix0 output: Matrix<T>&, the pseudoinverse of this Matrix */
  void pseudoinvert(double cutoff, Matrix<T>& Matrix0) const;

/// Diagonalizes this Matrix of Type "general".
/** This Matrix of Type "general" is diagonalized, eigenvalues and left and right eigenvectors are
    computed. LAPACK's XGEEV is used. This Matrix must be of type "general" and it must be square.
    \param W output: vector< complex<double> >&, the eigenvalues, must have size == dim0 == dim1
    \param Vr output: Matrix<T>&, the right eigenvectors, must be of dimension dim0 x dim0
    \param Vl output: Matrix<T>&, the left eigenvectors, must be of dimension dim0 x dim0 */
  void eigenDecompose(vector< complex<double> >& W, Matrix<T>& Vr, Matrix<T>& Vl) const;

/// Diagonalizes this Matrix of Type "hermitian".
/** This Matrix of Type "hermitian" is diagonalized, eigenvalues and eigenvectors are computed.
    LAPACK's XSYEV and XHEEV are used. This Matrix must be of type "hermitian" and it must be square.
    \param W output: vector<double>&, the eigenvalues, must have size == dim0 == dim1
    \param Vr output: Matrix<T>&, the eigenvectors, must be of dimension dim0 x dim0 */
  void eigenDecompose(vector<double>& W, Matrix<T>& Vr) const;

/// Solves generalized symmetric definite eigenproblem with divide-and-conquer.
/** The generalized eigenproblem A*x=w*B*x is solved, where A denotes this Matrix and B is Matrix0.
    LAPACK's XSYGVD and XHEGVD are used.
    This Matrix and Matrix0 must be of Type "hermitian" and Matrix0 must also be positive definite.
    On output this Matrix will contain the eigenvectors.
    \param Matrix0 input: Matrix<T>&, on input the right hand side matrix B,
                                      must be of Type "hermitian" and must be positive definite,
                                      on output it will be overwritten
    \param W output: vector<double>&, the eigenvalues, must fulfill (W.size() == dim0) */
  void eigenDecompose(Matrix<T>& Matrix0, vector<double>& W);

/// Solves generalized nonsymmetric eigenproblem.
/** The generalized eigenproblem A*x=w*B*x is solved, where A denotes this Matrix and B is Matrix0.
    LAPACK's XGGEVX are used.
    On input this matrix must be square, and on output it will be overwritten.
    \param Matrix0 input: Matrix<T>&, on input the right hand side matrix B,
                                      must have the same shape as this matrix,
                                      on output it will be overwritten
    \param balanc input: char, specifies the balance option,
                               (balanc == 'N'): no balancing,
                               (balanc == 'P'): permute only,
                               (balanc == 'S'): scale only,
                               (balanc == 'B'): both permute and scale
    \param jobvl input: char, (jobvl == 'N'): do not compute left eigenvectors,
                              (jobvl == 'V'): compute left eigenvectors
    \param jobvr input: char, (jobvr == 'N'): do not compute right eigenvectors,
                              (jobvr == 'V'): compute right eigenvectors
    \param sense input: char, determines which reciprocal condition numbers are computed,
                              (sense == 'N'): none,
                              (sense == 'E'): for eigenvalues only,
                              (sense == 'V'): for eigenvectors only,
                              (sense == 'B'): both for eigenvalues and eigenvectors
    \param Alpha output: vector< complex<double> >&, Alpha[k]/Beta[k] is eigenvalue k,
                         must fulfill (Alpha.size() == dim0)
    \param Beta output: vector< complex<double> >&, Alpha[k]/Beta[k] is eigenvalue k,
                        must fulfill (Beta.size() == dim0)
    \param Vl output: Matrix<T>&, (jobvl=='V'): the left eigenvectors,
                      must have the same shape as this matrix
    \param Vr output: Matrix<T>&, (jobvr=='V'): the right eigenvectors
                      must have the same shape as this matrix
    \param ilo output: int&, see LAPACK's documentation for XGGEVX
    \param ihi output: int&, see LAPACK's documentation for XGGEVX
    \param Lscale output: vector<double>&, see LAPACK's documentation for XGGEVX
                          must fulfill (Lscale.size() == dim0)
    \param Rscale output: vector<double>&, see LAPACK's documentation for XGGEVX
                          must fulfill (Lscale.size() == dim0)
    \param abnrm output: double&, one-norm of the balanced matrix A
    \param bbnrm output: double&, one-norm of the balanced matrix B
    \param Rconde output: vector<double>&, reciprocal condition numbers of the eigenvalues
    \param Rcondv output: vector<double>&, reciprocal condition numbers of the eigenvectors */
  void eigenDecompose(Matrix<T>& Matrix0, char balanc, char jobvl, char jobvr, char sense,
                      vector< complex<double> >& Alpha, vector< complex<double> >& Beta,
                      Matrix<T>& Vl, Matrix<T>& Vr,
                      int& ilo, int& ihi, vector<double>& Lscale, vector<double>& Rscale,
                      double& abnrm, double& bbnrm, vector<double>& Rconde, vector<double>& Rcondv);

/// Solves linear least squares problem.
/** The linear least squares problem A*x=b, where A denotes this Matrix, is solved.
    LAPACK's XGELSD is used.
    This Matrix must fulfill (this->Shape[0] >= this->Shape[1]) and it will be overwritten.
    \param rcond input: double, singular values <= rcond*S(1) are treated as zero
    \param b input/output: Matrix<T>&, on input the right hand side b,
                                       on output the solution x
    \param rank output: unsigned int&, the rank of this Matrix
    \param S output: vector<double>&, the singular values of this Matrix */
  void linearLeastSquares(double rcond, Matrix<T>& b, unsigned int& rank, vector<double>& S);

/// Computes ground state of this Matrix of Type "hermitian".
/** The ground state of this Matrix of Type "hermitian" is computed.
    ARPACK's XSAUPD, XSEUPD and LAPACK's XSYMV are used. This Matrix must be of type "hermitian", it
    must be square and real.
    \param e output: T&, the ground state energy
    \param V output: vector<T>&, the ground state, must have size dim0 == dim1 */
  void computeGroundState(T& e, vector<T>& V) const;

/// Computes matrix exponent exp(x * M).
/** This function computes the matrix exponent exp(x * M) from this Matrix M. It diagonalizes this
    Matrix with void Matrix<T>::eigenDecompose. This Matrix must be square, i.e.
    this->Shape[0]==this->Shape[1].
    \param x input: T, the prefactor in exp(x * M)
    \param MatrixExponent output: Matrix<T>&, the matrix exponent exp(x * M) of this Matrix M, on input
                                  it must have the same Shape as M, on output it is overwritten */
  void exponentialize(T x, Matrix<T>& MatrixExponent) const;

/// Returns singular value decomposed two body operator as MPO from the two body exponent exp(-delta*H).
/** This function computes a two body operator exp(-delta*H) from this two body Hamiltonian H. It
    performs a SVD on this operator Tensor=U*Sigma*Vt and returns two tensors TensorLeft=U*sqrt(Sigma)
    and TensorRight=sqrt(Sigma)*Vt representing the SVD. TensorLeft and TensorRight form the MPO
    representing exp(-delta*H). This function uses Tensor<T>::twoBodyOperatorSVD.
    This Matrix must be a two body Hamiltonian, i.e. it must be "hermitian" and this->Shape[0]==d0*d1
    and this->Shape[1]==d0*d1 where d0 and d1 are the physical dimensions of the two systems.
    \param delta input: T, the prefactor in exp(-delta*H)
    \param d0 input: unsigned int, the physical dimension of the left system
    \param d1 input: unsigned int, the physical dimension of the right system
    \param TensorLeft output: Tensor<T>&, the left tensor of the SVD, must have the correct form
    \param TensorRight output: Tensor<T>&, the right tensor of the SVD, must have the correct form
    \sa void Tensor<T>::twoBodyOperatorSVD(Tensor<T>& TensorLeft, Tensor<T>& TensorRight) */
  void twoBodyExponentSVD(T delta, unsigned int d0, unsigned int d1,
                          Tensor<T>& TensorLeft, Tensor<T>& TensorRight) const;

/// Computes lowest lying eigenstates of a real symmetric matrix.
/** The lowest lying eigenstates of a real symmetric matrix are computed. The non-zero diagonal elements shall
    be in AD and shall be indexed by IndexD, whereas the non-zero offdiagonal elements of the upper half shall
    be in AOD and shall be indexed by IndexODi and IndexODj. The kth element in AD, i.e. AD(k), is located at row and column
    IndexD(k) in the original matrix, and the kth element in AOD, i.e. AOD(k), is located at row IndexODi(k) and column IndexODj(k).
    ARPACK's XSAUPD, XSEUPD and a self-written matrix-vector multiplication subroutine are used. The latter
    subroutine makes use of symmetry, sparseness and OpenMP.
    T must be real.
    \param AD input: const vector<T>&, the non-zero diagonal elements
    \param IndexD input: const vector<int>&, the index to AD, must be 0-based
    \param AOD input: const vector<T>&, the non-zero offdiagonal elements of the upper half
    \param IndexODi input: const vector<int>&, the row index to AOD, must be 0-based
    \param IndexODj input: const vector<int>&, the column index to AOD, must be 0-based
    \param n input: int, the dimension of the original matrix
    \param nev input: int, the number of eigenvalues and eigenstates requested
    \param ncv input: int, the number of Arnoldi vectors used, ncv > nev must be satisfied
    \param Evals output: vector<T>&, the nev lowest lying eigenvalues
    \param Evecs output: Matrix<T>&, the nev lowest lying eigenstates in the columns */
  static void computeLowestEigenstates(const vector<T>& AD, const vector<int>& IndexD, const vector<T>& AOD,
                                       const vector<int>& IndexODi, const vector<int>& IndexODj,
                                       int n, int nev, int ncv, vector<T>& Evals, Matrix<T>& Evecs);

/// Computes low lying eigenstate for void GroundState<T>::computeLowEigenstate.
/** This function computes a low lying eigenstate for void GroundState<T>::computeLowEigenstate.
    The Hamiltonian is given by the product of tensors TensorLeft, TensorMiddle and TensorRight, and the
    resulting eigenstate specified by which is returned in Tensor0: which==0 implies the ground state,
    which==1 the first excited, and so on.
    ARPACK's XSAUPD, XSEUPD and a self-written tensor contraction are used. The latter needs
    2*d*D'*D^3+d^2*D'^2*D^2 operations for a matrix vector multiplication instead of d^2*D^4, where d is
    the physical dimension, D' is the MPO bond dimension and D is the MPS bond dimension.
    T must be real.
    \param TensorLeft input: const Tensor<T>&, the left tensor for the left part of the MPS MPO sandwich
    \param TensorMiddle input: const Tensor<T>&, the middle tensor of the MPO
    \param TensorRight input: const Tensor<T>&, the right tensor for the right part of the MPS MPO sandwich
    \param which input: unsigned int, specifies the returned eigenstate, which==0 implies the ground state,
                        which==1 the first excited, and so on
    \param Tensor0 output: Tensor<T>&, the desired eigenstate
    \sa void GroundState<T>::computeLowEigenstate(const Hamiltonian<T>& H0, unsigned int which, double eps,
                                                  unsigned int maxNumSweeps, double& energy,
                                                  unsigned int& numSweepsDone, MPS<T>& MPS0) const */
  static void computeArnoldi(const Tensor<T>& TensorLeft, const Tensor<T>& TensorMiddle,
                             const Tensor<T>& TensorRight, unsigned int which, Tensor<T>& Tensor0);

/// Computes expectation value of Matrix with Vector.
/** This function computes the expectation value of the input Matrix with the input Vector.
    \param Vector0 input: const vector<T>&, the vector
    \param Matrix0 input: const Matrix<T>&, the Matrix
    \return T, the resulting expectation value
    \sa friend T scalarProductMatrixVector<>(const vector<T>& Vector0, const Matrix<T>& Matrix0,
                                             const vector<T>& Vector1) */
  friend T expectationValueMatrix<>(const vector<T>& Vector0, const Matrix<T>& Matrix0);

/// Computes scalar product between Vector and Matrix applied to another Vector.
/** This friend function computes the scalar product of Vector0 as bra with Matrix0 applied to Vector1
    as ket.
    \param Vector0 input: const vector<T>&, the bra vector, must fulfill
                          Vector0.size()==Matrix0.getDim0()
    \param Matrix0 input: const Matrix<T>&, the Matrix
    \param Vector1 input: const vector<T>&, the ket vector, must fulfill
                          Vector1.size()==Matrix0.getDim1()
    \return T, the resulting scalar product
    \sa friend T expectationValueMatrix<>(const vector<T>& Vector0, const Matrix<T>& Matrix0)*/
  friend T scalarProductMatrixVector<>(const vector<T>& Vector0, const Matrix<T>& Matrix0,
                                       const vector<T>& Vector1);

/// Returns MPO representing local operator.
/** This friend function returns an MPO MPO0 that represents the local operator O at position x in a
    system of N spins.
    \param BC input: unsigned int, the boundary conditions
    \param N input: unsigned int, the number of spins
    \param x input: unsigned int, the position of the local operator O
    \param O input: const Matrix<T>&, the local operator
    \param MPO0 output: MPO<T>&, the resulting MPO representing the local operator O */
  friend void getMPOFromLocalOperator<>(string BC, unsigned int N, unsigned int x,
                                        const Matrix<T>& O, MPO<T>& MPO0);

/// Returns MPO representing sum of local operator.
/** This friend function returns an MPO MPO0 that represents the sum \sum_{l} O_{l} of a local operator O in a system of N spins.
    \param BC input: unsigned int, the boundary conditions
    \param N input: unsigned int, the number of spins
    \param O input: const Matrix<T>&, the local operator
    \param MPO0 output: MPO<T>&, the resulting MPO representing the sum of local operator O */
  friend void getMPOFromSumLocalOperator<>(string BC, unsigned int N, const Matrix<T>& O, MPO<T>& MPO0);

 protected:

/// Type of Matrix.
/** The type of this Matrix. It must be either "general" for a general matrix or "hermitian" for a
    hermitian matrix. By default it is "general". */
  string Type;
};

template<class T> Matrix<T>::Matrix() {}

template<class T> Matrix<T>::Matrix(unsigned int dim0, unsigned int dim1)
{
 this->rank = 2;
 vector<unsigned int> newShape(this->rank);
 newShape[0] = dim0; newShape[1] = dim1;
 this->Shape = newShape;
 this->size = dim0 * dim1;
 this->Elements = new T[this->size];
 this->Type = "general";
}

template<class T> Matrix<T>::Matrix(const Matrix<T>& Matrix0)
{
 this->rank = Matrix0.rank;
 this->Shape = Matrix0.Shape;
 this->size = Matrix0.size;
 this->Elements = new T[this->size];
 for (int i = 0; i < this->size; i++)
 {
  this->Elements[i] = Matrix0.Elements[i];
 }
 this->Type = Matrix0.Type;
}

template<class T> Matrix<T>::~Matrix() {}

template<class T> Matrix<T>& Matrix<T>::operator=(const Matrix<T>& Matrix0)
{
 if (this != &Matrix0)
 {
  this->rank = Matrix0.rank;
  this->Shape = Matrix0.Shape;
  this->size = Matrix0.size;
  delete[] this->Elements;
  this->Elements = new T[this->size];
  for (int i = 0; i < this->size; i++)
  {
   this->Elements[i] = Matrix0.Elements[i];
  }
  this->Type = Matrix0.Type;
 }
 return *this;
}

template<class T> inline T& Matrix<T>::operator()(unsigned int index0, unsigned int index1)
{
#ifdef DEBUG
 if ((index0 >= this->Shape[0]) || (index1 >= this->Shape[1]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> inline T& Matrix<T>::" <<
          "operator()(unsigned int index0, unsigned int index1): " <<
          "((index0 >= dim0) || (index1 >= dim1))." << endl;
  exit(1);
 }
#endif
 vector<unsigned int> Index(2);
 Index[0] = index0; Index[1] = index1;
 unsigned int position = this->getPosition(Index);
 return this->Elements[position];
}

template<class T> inline const T& Matrix<T>::operator()(unsigned int index0, unsigned int index1) const
{
#ifdef DEBUG
 if ((index0 >= this->Shape[0]) || (index1 >= this->Shape[1]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> inline const T& Matrix<T>::" <<
          "operator()(unsigned int index0, unsigned int index1) const: " <<
          "((index0 >= dim0) || (index1 >= dim1))." << endl;
  exit(1);
 }
#endif
 vector<unsigned int> Index(2);
 Index[0] = index0; Index[1] = index1;
 unsigned int position = this->getPosition(Index);
 return this->Elements[position];
}

template<class T> inline void Matrix<T>::setType(const string& Type0)
{
#ifdef DEBUG
 if ((Type0 != "general") && (Type0 != "hermitian"))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline void Matrix<T>::" <<
          "setType(const string& Type0)" << endl;
  exit(1);
 }
#endif
 this->Type = Type0;
}

template<class T> inline void Matrix<T>::getType(string& Type0) const
{
 Type0 = this->Type;
}

template<class T> void Matrix<T>::transpose()
{
#ifdef DEBUG
 if (this->size == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "transpose(): " <<
          "(this->size == 0)." << endl;
  exit(1);
 }
#endif
 unsigned int dim0 = this->Shape[0], dim1 = this->Shape[1];
 T* NewElements = new T[this->size];
 for (int col = 0; col < dim1; col++)
 {
  for (int row = 0; row < dim0; row++)
  {
   NewElements[row*dim1+col] = this->Elements[col*dim0+row];
  }
 }
 this->Shape[0] = dim1; this->Shape[1] = dim0;
 delete[] this->Elements;
 this->Elements = NewElements;
}

template<class T> void Matrix<T>::transpose(Matrix<T>& Matrix0) const
{
#ifdef DEBUG
 if (this->size == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "transpose(Matrix<T>& Matrix0) const: " <<
          "(this->size == 0)." << endl;
  exit(1);
 }
#endif
 Matrix0 = *this;
 Matrix0.transpose();
}

template<class T> void Matrix<T>::adjoint()
{
#ifdef DEBUG
 if (this->size == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "adjoint(): " <<
          "(this->size == 0)." << endl;
  exit(1);
 }
#endif
 this->transpose();
 this->complexConjugate();
}

template<class T> void Matrix<T>::adjoint(Matrix<T>& Matrix0) const
{
#ifdef DEBUG
 if (this->size == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "adjoint(Matrix<T>& Matrix0) const: " <<
          "(this->size == 0)." << endl;
  exit(1);
 }
#endif
 Matrix0 = *this;
 Matrix0.transpose();
 Matrix0.complexConjugate();
}

template<class T> void Matrix<T>::add(const Matrix<T>& Matrix0)
{
#ifdef DEBUG
 if ((this->size == 0) || (this->Shape[0] != Matrix0.Shape[0]) || (this->Shape[1] != Matrix0.Shape[1]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "add(const Matrix<T>& Matrix0): " <<
          "((this->size == 0) || (this->Shape[0] != Matrix0.Shape[0]) || " <<
          "(this->Shape[1] != Matrix0.Shape[1]))." << endl;
  exit(1);
 }
#endif
 for (int i = 0; i < this->size; i++)
 {
  this->Elements[i] += Matrix0.Elements[i];
 }
}

template<class T> void Matrix<T>::multiply(T x)
{
#ifdef DEBUG
 if (this->size == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "multiply(T x): " <<
          "(this->size == 0)." << endl;
  exit(1);
 }
#endif
 for (int i = 0; i < this->size; i++)
 {
  this->Elements[i] = x * this->Elements[i];
 }
}

template<class T> void Matrix<T>::multiply(const vector<T>& Vector0, vector<T>& Vector1) const
{
 int m = this->Shape[0]; int n = this->Shape[1];
#ifdef DEBUG
 if ((n != Vector0.size()) || (m != Vector1.size()))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "multiply(const vector<T>& Vector0, vector<T>& Vector1) const: " <<
          "((n != Vector0.size()) || (m != Vector1.size()))." << endl;
  exit(1);
 }
#endif
 char trans = 'N';
 T* X = new T[n]; T* Y = new T[m];
 for (int i = 0; i < n; i++)
  X[i] = Vector0[i];
 int incx = 1, incy = 1;
 if (typeid(T) == typeid(float))
 {
  float alpha = 1.0, beta = 0.0;
  sgemv_(&trans, &m, &n, &alpha, (float*)this->Elements, &m, (float*)X, &incx, &beta, (float*)Y, &incy);
 }
 else if (typeid(T) == typeid(double))
 {
  double alpha = 1.0, beta = 0.0;
  dgemv_(&trans, &m, &n, &alpha, (double*)this->Elements, &m, (double*)X, &incx, &beta, (double*)Y,
         &incy);
 }
 else if (typeid(T) == typeid(complex<float>))
 {
  complex<float> alpha = 1.0, beta = 0.0;
  cgemv_(&trans, &m, &n, &alpha, (complex<float>*)this->Elements, &m, (complex<float>*)X, &incx, &beta,
         (complex<float>*)Y, &incy);
 }
 else if (typeid(T) == typeid(complex<double>))
 {
  complex<double> alpha = 1.0, beta = 0.0;
  zgemv_(&trans, &m, &n, &alpha, (complex<double>*)this->Elements, &m, (complex<double>*)X, &incx,
         &beta, (complex<double>*)Y, &incy);
 }
 for (int i = 0; i < m; i++)
 {
  Vector1[i] = Y[i];
 }
 delete[] X; delete[] Y;
}

template<class T> void Matrix<T>::multiply(const Matrix<T>& Matrix0)
{
 int m = this->Shape[0]; int n = this->Shape[1];
 int m0 = Matrix0.Shape[0]; int n0 = Matrix0.Shape[1];
#ifdef DEBUG
 if (n != m0)
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void Matrix<T>::" <<
          "multiply(const Matrix<T>& Matrix0)" << endl;
  exit(1);
 }
#endif
 char transa = 'N', transb = 'N';
 unsigned int newSize = m * n0;
 T* newElements = new T[newSize];
 if (typeid(T) == typeid(float))
 {
  float alpha = 1.0, beta = 0.0;
  sgemm_(&transa, &transb, &m, &n0, &n, &alpha, (float*)this->Elements, &m,
         (float*)Matrix0.Elements, &n, &beta, (float*)newElements, &m);
 }
 else if (typeid(T) == typeid(double))
 {
  double alpha = 1.0, beta = 0.0;
  dgemm_(&transa, &transb, &m, &n0, &n, &alpha, (double*)this->Elements, &m,
         (double*)Matrix0.Elements, &n, &beta, (double*)newElements, &m);
 }
 else if (typeid(T) == typeid(complex<float>))
 {
  complex<float> alpha = 1.0, beta = 0.0;
  cgemm_(&transa, &transb, &m, &n0, &n, &alpha, (complex<float>*)this->Elements, &m,
         (complex<float>*)Matrix0.Elements, &n, &beta, (complex<float>*)newElements, &m);
 }
 else if (typeid(T) == typeid(complex<double>))
 {
  complex<double> alpha = 1.0, beta = 0.0;
  zgemm_(&transa, &transb, &m, &n0, &n, &alpha, (complex<double>*)this->Elements, &m,
         (complex<double>*)Matrix0.Elements, &n, &beta, (complex<double>*)newElements, &m);
 }
 this->Shape[0] = m; this->Shape[1] = n0;
 this->size = newSize;
 delete[] this->Elements;
 this->Elements = newElements;
}

template<class T> void Matrix<T>::multiply(const Matrix<T>& Matrix0, Matrix<T>& Matrix1) const
{
 int m = this->Shape[0]; int n = this->Shape[1];
 int m0 = Matrix0.Shape[0]; int n0 = Matrix0.Shape[1];
#ifdef DEBUG
 if ((n != m0) || (Matrix1.Shape[0] != m) || (Matrix1.Shape[1] != n0))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void Matrix<T>::" <<
          "multiply(const Matrix<T>& Matrix0, Matrix<T>& Matrix1) const" << endl;
  exit(1);
 }
#endif
 char transa = 'N', transb = 'N';
 if (typeid(T) == typeid(float))
 {
  float alpha = 1.0, beta = 0.0;
  sgemm_(&transa, &transb, &m, &n0, &n, &alpha, (float*)this->Elements, &m,
         (float*)Matrix0.Elements, &n, &beta, (float*)Matrix1.Elements, &m);
 }
 else if (typeid(T) == typeid(double))
 {
  double alpha = 1.0, beta = 0.0;
  dgemm_(&transa, &transb, &m, &n0, &n, &alpha, (double*)this->Elements, &m,
         (double*)Matrix0.Elements, &n, &beta, (double*)Matrix1.Elements, &m);
 }
 else if (typeid(T) == typeid(complex<float>))
 {
  complex<float> alpha = 1.0, beta = 0.0;
  cgemm_(&transa, &transb, &m, &n0, &n, &alpha, (complex<float>*)this->Elements, &m,
         (complex<float>*)Matrix0.Elements, &n, &beta, (complex<float>*)Matrix1.Elements, &m);
 }
 else if (typeid(T) == typeid(complex<double>))
 {
  complex<double> alpha = 1.0, beta = 0.0;
  zgemm_(&transa, &transb, &m, &n0, &n, &alpha, (complex<double>*)this->Elements, &m,
         (complex<double>*)Matrix0.Elements, &n, &beta, (complex<double>*)Matrix1.Elements, &m);
 }
}

template<class T> void Matrix<T>::multiplyDirectProduct(const Matrix<T>& Matrix0, Matrix<T>& Matrix1)
                                  const
{
 unsigned int dim0 = this->Shape[0]; unsigned int dim1 = this->Shape[1];
 unsigned int dim00 = Matrix0.Shape[0]; unsigned int dim10 = Matrix0.Shape[1];
 unsigned int dim01 = dim0 * dim00;
 unsigned int dim11 = dim1 * dim10;
 Matrix1 = Matrix<T>(dim01, dim11);
 int i1, j1;
 for (int j = 0; j < dim1; j++)
 {
  for (int i = 0; i < dim0; i++)
  {
   for (int j0 = 0; j0 < dim10; j0++)
   {
    for (int i0 = 0; i0 < dim00; i0++)
    {
     i1 = i * dim00 + i0;
     j1 = j * dim10 + j0;
     Matrix1(i1, j1) = this->Elements[j*dim0+i] * Matrix0.Elements[j0*dim00+i0];
    }
   }
  }
 }
}

template<class T> void Matrix<T>::QRDecompose(Matrix<T>& R)
{
 char side = 'L', trans = 'N';
 int info, s;
 int m = this->Shape[0], n = this->Shape[1];
 int k = min(m, n), lwork = n;
 int o = m*k;
 T* Tau = new T[k];
 T* Work = new T[lwork];
 Matrix<T> Q(m, k);
 for (int i = 0; i < o; i++)
  Q.Elements[i] = 0.0;
 for (int i = 0; i < k; i++)
  Q.Elements[(m+1)*i] = 1.0;
 R = Matrix<T>(k, n);
 if (typeid(T) == typeid(float))
 {
  sgeqrf_(&m, &n, (float*)this->Elements, &m, (float*)Tau, (float*)Work, &lwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "QRDecompose(Matrix<T>& R): " <<
           "INFO = " << info << " in LAPACK's SGEQRF." << endl;
   exit(1);
  }
  for (int j = 0; j < n; j++)
  {
   s = min(j, k-1);
   for (int i = 0; i <= s; i++)
    R.Elements[k*j+i] = this->Elements[m*j+i];
   for (int i = s+1; i < k; i++)
    R.Elements[k*j+i] = 0.0;
  }
  sormqr_(&side, &trans, &m, &k, &k, (float*)this->Elements, &m, (float*)Tau, (float*)Q.Elements, &m,
          (float*)Work, &lwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "QRDecompose(Matrix<T>& R): " <<
           "INFO = " << info << " in LAPACK's SORMQR." << endl;
   exit(1);
  }
 }
 else if (typeid(T) == typeid(double))
 {
  dgeqrf_(&m, &n, (double*)this->Elements, &m, (double*)Tau, (double*)Work, &lwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "QRDecompose(Matrix<T>& R): " <<
           "INFO = " << info << " in LAPACK's DGEQRF." << endl;
   exit(1);
  }
  for (int j = 0; j < n; j++)
  {
   s = min(j, k-1);
   for (int i = 0; i <= s; i++)
    R.Elements[k*j+i] = this->Elements[m*j+i];
   for (int i = s+1; i < k; i++)
    R.Elements[k*j+i] = 0.0;
  }
  dormqr_(&side, &trans, &m, &k, &k, (double*)this->Elements, &m, (double*)Tau, (double*)Q.Elements, &m,
          (double*)Work, &lwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "QRDecompose(Matrix<T>& R): " <<
           "INFO = " << info << " in LAPACK's DORMQR." << endl;
   exit(1);
  }
 }
 else if (typeid(T) == typeid(complex<float>))
 {
  cgeqrf_(&m, &n, (complex<float>*)this->Elements, &m, (complex<float>*)Tau, (complex<float>*)Work,
          &lwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "QRDecompose(Matrix<T>& R): " <<
           "INFO = " << info << " in LAPACK's CGEQRF." << endl;
   exit(1);
  }
  for (int j = 0; j < n; j++)
  {
   s = min(j, k-1);
   for (int i = 0; i <= s; i++)
    R.Elements[k*j+i] = this->Elements[m*j+i];
   for (int i = s+1; i < k; i++)
    R.Elements[k*j+i] = 0.0;
  }
  cunmqr_(&side, &trans, &m, &k, &k, (complex<float>*)this->Elements, &m, (complex<float>*)Tau,
          (complex<float>*)Q.Elements, &m, (complex<float>*)Work, &lwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "QRDecompose(Matrix<T>& R): " <<
           "INFO = " << info << " in LAPACK's CUNMQR." << endl;
   exit(1);
  }
 }
 else if (typeid(T) == typeid(complex<double>))
 {
  zgeqrf_(&m, &n, (complex<double>*)this->Elements, &m, (complex<double>*)Tau, (complex<double>*)Work,
          &lwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "QRDecompose(Matrix<T>& R): " <<
           "INFO = " << info << " in LAPACK's ZGEQRF." << endl;
   exit(1);
  }
  for (int j = 0; j < n; j++)
  {
   s = min(j, k-1);
   for (int i = 0; i <= s; i++)
    R.Elements[k*j+i] = this->Elements[m*j+i];
   for (int i = s+1; i < k; i++)
    R.Elements[k*j+i] = 0.0;
  }
  zunmqr_(&side, &trans, &m, &k, &k, (complex<double>*)this->Elements, &m, (complex<double>*)Tau,
          (complex<double>*)Q.Elements, &m, (complex<double>*)Work, &lwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "QRDecompose(Matrix<T>& R): " <<
           "INFO = " << info << " in LAPACK's ZUNMQR." << endl;
   exit(1);
  }
 }
 delete[] Tau;
 delete[] Work;
 *this = Q;
}

template<class T> void Matrix<T>::singularValueDecompose(Matrix<T>& U, Matrix<T>& Sigma, Matrix<T>& Vt) const
{
 int m = this->Shape[0]; int n = this->Shape[1];
#ifdef DEBUG
 if ((U.Shape[0] != m) || (U.Shape[1] != m) || (Sigma.Shape[0] != m) || (Sigma.Shape[1] != n) ||
     (Vt.Shape[0] != n) || (Vt.Shape[1] != n))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void Matrix<T>::" <<
          "singularValueDecompose(Matrix<T>& U, Matrix<T>& Sigma, Matrix<T>& Vt) const" << endl;
  exit(1);
 }
#endif
 char jobu = 'A', jobvt = 'A';
 T* A = new T[this->size];
 for (int i = 0; i < this->size; i++)
 {
  A[i] = this->Elements[i];
 }
 int info;
 Sigma.fillZeroes();
 if (typeid(T) == typeid(float))
 {
  float* S = new float[min(m,n)];
  int lwork = max(3*min(m,n)+max(m,n), 5*min(m,n));
  float* Work = new float[lwork];
  sgesvd_(&jobu, &jobvt, &m, &n, (float*)A, &m, S, (float*)U.Elements,
          &m, (float*)Vt.Elements, &n, Work, &lwork, &info);
  for (int i = 0; i < min(m,n); i++)
  {
   Sigma(i,i) = S[i];
  }
  delete[] S;
  delete[] Work;
 }
 else if (typeid(T) == typeid(double))
 {
  double* S = new double[min(m,n)];
  int lwork = max(3*min(m,n)+max(m,n), 5*min(m,n));
  double* Work = new double[lwork];
  dgesvd_(&jobu, &jobvt, &m, &n, (double*)A, &m, S, (double*)U.Elements,
          &m, (double*)Vt.Elements, &n, Work, &lwork, &info);
  for (int i = 0; i < min(m,n); i++)
  {
   Sigma(i,i) = S[i];
  }
  delete[] S;
  delete[] Work;
 }
 else if (typeid(T) == typeid(complex<float>))
 {
  float* S = new float[min(m,n)];
  int lwork = 2 * min(m,n) + max(m,n);
  complex<float>* Work = new complex<float>[lwork];
  float* Rwork = new float[5 * min(m,n)];
  cgesvd_(&jobu, &jobvt, &m, &n, (complex<float>*)A, &m, S, (complex<float>*)U.Elements,
          &m, (complex<float>*)Vt.Elements, &n, Work, &lwork, Rwork, &info);
  for (int i = 0; i < min(m,n); i++)
  {
   Sigma(i,i) = S[i];
  }
  delete[] S;
  delete[] Work;
  delete[] Rwork;
 }
 else if (typeid(T) == typeid(complex<double>))
 {
  double* S = new double[min(m,n)];
  int lwork = 2 * min(m,n) + max(m,n);
  complex<double>* Work = new complex<double>[lwork];
  double* Rwork = new double[5 * min(m,n)];
  zgesvd_(&jobu, &jobvt, &m, &n, (complex<double>*)A, &m, S, (complex<double>*)U.Elements,
          &m, (complex<double>*)Vt.Elements, &n, Work, &lwork, Rwork, &info);
  for (int i = 0; i < min(m,n); i++)
  {
   Sigma(i,i) = S[i];
  }
  delete[] S;
  delete[] Work;
  delete[] Rwork;
 }
 delete[] A;
}

template<class T> void Matrix<T>::pseudoinvert(double cutoff, Matrix<T>& Matrix0) const
{
#ifdef DEBUG
 if (this->size == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "pseudoinvert(double cutoff, Matrix<T>& Matrix0) const: " <<
          "(this->size == 0)." << endl;
  exit(1);
 }
#endif
 unsigned int dim0 = this->getDim0(), dim1 = this->getDim1();
 Matrix<T> U(dim0, dim0), Sigma(dim0, dim1), Vt(dim1, dim1);
 this->singularValueDecompose(U, Sigma, Vt);
 Matrix0 = Matrix<T>(dim1, dim0);
 Matrix0.fillZeroes();
 for (int k = 0; k < min(dim0, dim1); k++)
 {
  if (MathAuxiliary::convertToDouble(Sigma(k, k))/MathAuxiliary::convertToDouble(Sigma(0, 0)) > cutoff)
  {
   for (int j = 0; j < dim0; j++)
   {
    for (int i = 0; i < dim1; i++)
    {
     Matrix0(i, j) += MathAuxiliary::complexConjugate(Vt(k, i))*MathAuxiliary::complexConjugate(U(j, k))/Sigma(k, k);
    }
   }
  }
 }
}

template<class T> void Matrix<T>::eigenDecompose(vector< complex<double> >& W, Matrix<T>& Vr,
                                                 Matrix<T>& Vl) const
{
 int m = this->Shape[0]; int n = this->Shape[1];
#ifdef DEBUG
 if ((m != n) || (this->Type != "general") || (W.size() != m) || (Vr.Shape[0] != m) ||
     (Vr.Shape[1] != m) || (Vl.Shape[0] != m) || (Vl.Shape[1] != m))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void Matrix<T>::" <<
          "eigenDecompose(vector< complex<double> >& W, Matrix<T>& Vr, " <<
          "Matrix<T>& Vl) const" << endl;
  exit(1);
 }
#endif
 char jobvl = 'V', jobvr = 'V';
 T* A = new T[this->size];
 for (int i = 0; i < this->size; i++)
 {
  A[i] = this->Elements[i];
 }
 int info;
 if (typeid(T) == typeid(float))
 {
  float* Wr = new float[m];
  float* Wi = new float[m];
  int lwork = 4*m;
  float* Work = new float[lwork];
  sgeev_(&jobvl, &jobvr, &m, (float*)A, &m, Wr, Wi, (float*)Vl.Elements, &m, (float*)Vr.Elements,
         &m, Work, &lwork, &info);
  for (int i = 0; i < m; i++)
  {
   W[i] = complex<double>(Wr[i], Wi[i]);
  }
  delete[] Wr;
  delete[] Wi;
  delete[] Work;
 }
 else if (typeid(T) == typeid(double))
 {
  double* Wr = new double[m];
  double* Wi = new double[m];
  int lwork = 4*m;
  double* Work = new double[lwork];
  dgeev_(&jobvl, &jobvr, &m, (double*)A, &m, Wr, Wi, (double*)Vl.Elements, &m,
         (double*)Vr.Elements, &m, Work, &lwork, &info);
  for (int i = 0; i < m; i++)
  {
   W[i] = complex<double>(Wr[i], Wi[i]);
  }
  delete[] Wr;
  delete[] Wi;
  delete[] Work;
 }
 else if (typeid(T) == typeid(complex<float>))
 {
  complex<float>* W0 = new complex<float>[m];
  int lwork = 2*m;
  complex<float>* Work = new complex<float>[lwork];
  float* Rwork = new float[2*m];
  cgeev_(&jobvl, &jobvr, &m, (complex<float>*)A, &m, W0, (complex<float>*)Vl.Elements, &m,
         (complex<float>*)Vr.Elements, &m, Work, &lwork, Rwork, &info);
  for (int i = 0; i < m; i++)
  {
   W[i] = W0[i];
  }
  delete[] W0;
  delete[] Work;
  delete[] Rwork;
 }
 else if (typeid(T) == typeid(complex<double>))
 {
  complex<double>* W0 = new complex<double>[m];
  int lwork = 2*m;
  complex<double>* Work = new complex<double>[lwork];
  double* Rwork = new double[2*m];
  zgeev_(&jobvl, &jobvr, &m, (complex<double>*)A, &m, W0, (complex<double>*)Vl.Elements, &m,
         (complex<double>*)Vr.Elements, &m, Work, &lwork, Rwork, &info);
  for (int i = 0; i < m; i++)
  {
   W[i] = W0[i];
  }
  delete[] W0;
  delete[] Work;
  delete[] Rwork;
 }
 delete[] A;
}

template<class T> void Matrix<T>::eigenDecompose(vector<double>& W, Matrix<T>& Vr) const
{
 int m = this->Shape[0]; int n = this->Shape[1];
#ifdef DEBUG
 if ((m != n) || (this->Type != "hermitian") || (W.size() != m) || (Vr.Shape[0] != m) ||
     (Vr.Shape[1] != m))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void Matrix<T>::" <<
          "eigenDecompose(vector<double>& W, Matrix<T>& Vr) const" << endl;
  exit(1);
 }
#endif
 char jobz = 'V', uplo = 'U';
 T* A = new T[this->size];
 for (int i = 0; i < this->size; i++)
 {
  A[i] = this->Elements[i];
 }
 int info;
 if (typeid(T) == typeid(float))
 {
  float* W0 = new float[m];
  int lwork = 3*m-1;
  float* Work = new float[lwork];
  ssyev_(&jobz, &uplo, &m, (float*)A, &m, W0, Work, &lwork, &info);
  for (int i = 0; i < m; i++)
  {
   W[i] = (double)W0[i];
  }
  for (int i = 0; i < Vr.size; i++)
  {
   Vr.Elements[i] = A[i];
  }
  delete[] W0;
  delete[] Work;
 }
 else if (typeid(T) == typeid(double))
 {
  double* W0 = new double[m];
  int lwork = 3*m-1;
  double* Work = new double[lwork];
  dsyev_(&jobz, &uplo, &m, (double*)A, &m, W0, Work, &lwork, &info);
  for (int i = 0; i < m; i++)
  {
   W[i] = W0[i];
  }
  for (int i = 0; i < Vr.size; i++)
  {
   Vr.Elements[i] = A[i];
  }
  delete[] W0;
  delete[] Work;
 }
 else if (typeid(T) == typeid(complex<float>))
 {
  float* W0 = new float[m];
  int lwork = 2*m-1;
  complex<float>* Work = new complex<float>[lwork];
  float* Rwork = new float[3*n-2];
  cheev_(&jobz, &uplo, &m, (complex<float>*)A, &m, W0, Work, &lwork, Rwork, &info);
  for (int i = 0; i < m; i++)
  {
   W[i] = (double)W0[i];
  }
  for (int i = 0; i < Vr.size; i++)
  {
   Vr.Elements[i] = A[i];
  }
  delete[] W0;
  delete[] Work;
  delete[] Rwork;
 }
 else if (typeid(T) == typeid(complex<double>))
 {
  double* W0 = new double[m];
  int lwork = 2*m-1;
  complex<double>* Work = new complex<double>[lwork];
  double* Rwork = new double[3*n-2];
  zheev_(&jobz, &uplo, &m, (complex<double>*)A, &m, W0, Work, &lwork, Rwork, &info);
  for (int i = 0; i < m; i++)
  {
   W[i] = W0[i];
  }
  for (int i = 0; i < Vr.size; i++)
  {
   Vr.Elements[i] = A[i];
  }
  delete[] W0;
  delete[] Work;
  delete[] Rwork;
 }
 delete[] A;
}

template<class T> void Matrix<T>::eigenDecompose(Matrix<T>& Matrix0, vector<double>& W)
{
#ifdef DEBUG
 if ((this->size == 0) || (this->Shape[0] != this->Shape[1]) || (Matrix0.Shape[0] != this->Shape[0]) ||
     (Matrix0.Shape[1] != this->Shape[0]) || (W.size() != this->Shape[0]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "eigenDecompose(Matrix<T>& Matrix0, vector<double>& W): " <<
          "((this->size == 0) || (this->Shape[0] != this->Shape[1]) || (Matrix0.Shape[0] != this->Shape[0]) || " <<
           "(Matrix0.Shape[1] != this->Shape[0]) || (W.size() != this->Shape[0]))." << endl;
  exit(1);
 }
 else if ((this->Type != "hermitian") || (Matrix0.Type != "hermitian"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "eigenDecompose(Matrix<T>& Matrix0, vector<double>& W): " <<
          "((this->Type != hermitian) || (Matrix0.Type != hermitian))." << endl;
  exit(1);
 }
#endif
 int itype = 1;
 char jobz = 'V', uplo = 'U';
 int n = this->Shape[0];
 int lwork, info;
 int liwork = 3+5*n;
 int* Iwork = new int[liwork];
 if (typeid(T) == typeid(float))
 {
  float* W0 = new float[n];
  lwork = 1+6*n+2*n*n;
  float* Work = new float[lwork];
  ssygvd_(&itype, &jobz, &uplo, &n, (float*)this->Elements, &n, (float*)Matrix0.Elements, &n, W0,
          Work, &lwork, Iwork, &liwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "eigenDecompose(Matrix<T>& Matrix0, vector<double>& W): " <<
           "INFO = " << info << " in LAPACK's SSYGVD." << endl;
   exit(1);
  }
  for (int i = 0; i < n; i++)
  {
   W[i] = (double)W0[i];
  }
  delete[] W0;
  delete[] Work;
 }
 else if (typeid(T) == typeid(double))
 {
  double* W0 = new double[n];
  lwork = 1+6*n+2*n*n;
  double* Work = new double[lwork];
  dsygvd_(&itype, &jobz, &uplo, &n, (double*)this->Elements, &n, (double*)Matrix0.Elements, &n, W0,
          Work, &lwork, Iwork, &liwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "eigenDecompose(Matrix<T>& Matrix0, vector<double>& W): " <<
           "INFO = " << info << " in LAPACK's DSYGVD." << endl;
   exit(1);
  }
  for (int i = 0; i < n; i++)
  {
   W[i] = (double)W0[i];
  }
  delete[] W0;
  delete[] Work;
 }
 else if (typeid(T) == typeid(complex<float>))
 {
  float* W0 = new float[n];
  lwork = 2*n+n*n;
  complex<float>* Work = new complex<float>[lwork];
  int lrwork = 1+5*n+2*n*n;
  float* Rwork = new float[lrwork];
  chegvd_(&itype, &jobz, &uplo, &n, (complex<float>*)this->Elements, &n, (complex<float>*)Matrix0.Elements, &n, W0,
          Work, &lwork, Rwork, &lrwork, Iwork, &liwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "eigenDecompose(Matrix<T>& Matrix0, vector<double>& W): " <<
           "INFO = " << info << " in LAPACK's CHEGVD." << endl;
   exit(1);
  }
  for (int i = 0; i < n; i++)
  {
   W[i] = (double)W0[i];
  }
  delete[] W0;
  delete[] Work;
  delete[] Rwork;
 }
 else if (typeid(T) == typeid(complex<double>))
 {
  double* W0 = new double[n];
  lwork = 2*n+n*n;
  complex<double>* Work = new complex<double>[lwork];
  int lrwork = 1+5*n+2*n*n;
  double* Rwork = new double[lrwork];
  zhegvd_(&itype, &jobz, &uplo, &n, (complex<double>*)this->Elements, &n, (complex<double>*)Matrix0.Elements, &n, W0,
          Work, &lwork, Rwork, &lrwork, Iwork, &liwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "eigenDecompose(Matrix<T>& Matrix0, vector<double>& W): " <<
           "INFO = " << info << " in LAPACK's ZHEGVD." << endl;
   exit(1);
  }
  for (int i = 0; i < n; i++)
  {
   W[i] = (double)W0[i];
  }
  delete[] W0;
  delete[] Work;
  delete[] Rwork;
 }
 delete[] Iwork;
}

template<class T> void Matrix<T>::eigenDecompose(Matrix<T>& Matrix0, char balanc, char jobvl, char jobvr, char sense,
                                                 vector< complex<double> >& Alpha, vector< complex<double> >& Beta,
                                                 Matrix<T>& Vl, Matrix<T>& Vr,
                                                 int& ilo, int& ihi, vector<double>& Lscale, vector<double>& Rscale,
                                                 double& abnrm, double& bbnrm,
                                                 vector<double>& Rconde, vector<double>& Rcondv)
{
#ifdef DEBUG
 if ((this->size == 0) || (this->Shape[0] != this->Shape[1]) || (Matrix0.Shape[0] != this->Shape[0]) ||
     (Matrix0.Shape[1] != this->Shape[0]) || (Alpha.size() != this->Shape[0]) || (Beta.size() != this->Shape[0]) ||
     (Lscale.size() != this->Shape[0]) || (Rscale.size() != this->Shape[0]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "eigenDecompose(Matrix<T>& Matrix0, char balanc, char jobvl, char jobvr, char sense, " <<
                         "vector< complex<double> >& Alpha, vector< complex<double> >& Beta, " <<
                         "Matrix<T>& Vl, Matrix<T>& Vr, int& ilo, int& ihi, " <<
                         "vector<double>& Lscale, vector<double>& Rscale, double& abnrm, double& bbnrm, " <<
                         "vector<double>& Rconde, vector<double>& Rcondv): " <<
          "((this->size == 0) || (this->Shape[0] != this->Shape[1]) || (Matrix0.Shape[0] != this->Shape[0]) || " <<
           "(Matrix0.Shape[1] != this->Shape[0]) || (Alpha.size() != this->Shape[0]) || " <<
           "(Beta.size() != this->Shape[0]) || (Lscale.size() != this->Shape[0]) || " <<
           "(Rscale.size() != this->Shape[0]))." << endl;
  exit(1);
 }
 else if (((balanc != 'N') && (balanc != 'P') && (balanc != 'S') && (balanc != 'B')) ||
          ((jobvl != 'N') && (jobvl != 'V')) || ((jobvr != 'N') && (jobvr != 'V')) ||
          ((sense != 'N') && (sense != 'E') && (sense != 'V') && (sense != 'B')))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "eigenDecompose(Matrix<T>& Matrix0, char balanc, char jobvl, char jobvr, char sense, " <<
                         "vector< complex<double> >& Alpha, vector< complex<double> >& Beta, " <<
                         "Matrix<T>& Vl, Matrix<T>& Vr, int& ilo, int& ihi, " <<
                         "vector<double>& Lscale, vector<double>& Rscale, double& abnrm, double& bbnrm, " <<
                         "vector<double>& Rconde, vector<double>& Rcondv): " <<
          "(((balanc != 'N') && (balanc != 'P') && (balanc != 'S') && (balanc != 'B')) || " <<
           "((jobvl != 'N') && (jobvl != 'V')) || ((jobvr != 'N') && (jobvr != 'V')) || " <<
           "((sense != 'N') && (sense != 'E') && (sense != 'V') && (sense != 'B')))." << endl;
  exit(1);
 }
 else if ((jobvl == 'V') && ((Vl.Shape[0] != this->Shape[0]) || (Vl.Shape[1] != this->Shape[0])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "eigenDecompose(Matrix<T>& Matrix0, char balanc, char jobvl, char jobvr, char sense, " <<
                         "vector< complex<double> >& Alpha, vector< complex<double> >& Beta, " <<
                         "Matrix<T>& Vl, Matrix<T>& Vr, int& ilo, int& ihi, " <<
                         "vector<double>& Lscale, vector<double>& Rscale, double& abnrm, double& bbnrm, " <<
                         "vector<double>& Rconde, vector<double>& Rcondv): " <<
          "((jobvl == 'V') && ((Vl.Shape[0] != this->Shape[0]) || (Vl.Shape[1] != this->Shape[0])))." << endl;
  exit(1);
 }
 else if ((jobvr == 'V') && ((Vr.Shape[0] != this->Shape[0]) || (Vr.Shape[1] != this->Shape[0])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "eigenDecompose(Matrix<T>& Matrix0, char balanc, char jobvl, char jobvr, char sense, " <<
                         "vector< complex<double> >& Alpha, vector< complex<double> >& Beta, " <<
                         "Matrix<T>& Vl, Matrix<T>& Vr, int& ilo, int& ihi, " <<
                         "vector<double>& Lscale, vector<double>& Rscale, double& abnrm, double& bbnrm, " <<
                         "vector<double>& Rconde, vector<double>& Rcondv): " <<
          "((jobvr == 'V') && ((Vr.Shape[0] != this->Shape[0]) || (Vr.Shape[1] != this->Shape[0])))." << endl;
  exit(1);
 }
 else if (((sense == 'E') || (sense == 'B')) && (Rconde.size() != this->Shape[0]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "eigenDecompose(Matrix<T>& Matrix0, char balanc, char jobvl, char jobvr, char sense, " <<
                         "vector< complex<double> >& Alpha, vector< complex<double> >& Beta, " <<
                         "Matrix<T>& Vl, Matrix<T>& Vr, int& ilo, int& ihi, " <<
                         "vector<double>& Lscale, vector<double>& Rscale, double& abnrm, double& bbnrm, " <<
                         "vector<double>& Rconde, vector<double>& Rcondv): " <<
          "(((sense == 'E') || (sense == 'B')) && (Rconde.size() != this->Shape[0]))." << endl;
  exit(1);
 }
 else if (((sense == 'V') || (sense == 'B')) && (Rcondv.size() != this->Shape[0]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "eigenDecompose(Matrix<T>& Matrix0, char balanc, char jobvl, char jobvr, char sense, " <<
                         "vector< complex<double> >& Alpha, vector< complex<double> >& Beta, " <<
                         "Matrix<T>& Vl, Matrix<T>& Vr, int& ilo, int& ihi, " <<
                         "vector<double>& Lscale, vector<double>& Rscale, double& abnrm, double& bbnrm, " <<
                         "vector<double>& Rconde, vector<double>& Rcondv): " <<
          "(((sense == 'V') || (sense == 'B')) && (Rcondv.size() != this->Shape[0]))." << endl;
  exit(1);
 }
#endif
 int n = this->Shape[0];
 int lwork, info;
 int* Bwork = new int[n];
 if (typeid(T) == typeid(float))
 {
  float* Alphar0 = new float[n]; float* Alphai0 = new float[n]; float* Beta0 = new float[n];
  float* Lscale0 = new float[n]; float* Rscale0 = new float[n];
  float abnrm0, bbnrm0;
  float* Rconde0 = new float[n]; float* Rcondv0 = new float[n];
  lwork = 2*n*n+12*n+16;
  float* Work = new float[lwork]; int* Iwork = new int[n+6];
  sggevx_(&balanc, &jobvl, &jobvr, &sense, &n, (float*)this->Elements, &n, (float*)Matrix0.Elements, &n,
          Alphar0, Alphai0, Beta0, (float*)Vl.Elements, &n, (float*)Vr.Elements, &n, &ilo, &ihi,
          Lscale0, Rscale0, &abnrm0, &bbnrm0, Rconde0, Rcondv0, Work, &lwork, Iwork, Bwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "eigenDecompose(Matrix<T>& Matrix0, char balanc, char jobvl, char jobvr, char sense, " <<
                          "vector< complex<double> >& Alpha, vector< complex<double> >& Beta, " <<
                          "Matrix<T>& Vl, Matrix<T>& Vr, int& ilo, int& ihi, " <<
                          "vector<double>& Lscale, vector<double>& Rscale, double& abnrm, double& bbnrm, " <<
                          "vector<double>& Rconde, vector<double>& Rcondv): " <<
           "INFO = " << info << " in LAPACK's SGGEVX." << endl;
   exit(1);
  }
  for (int i = 0; i < n; i++)
  {
   Alpha[i] = complex<double>((double)Alphar0[i], (double)Alphai0[i]);
   Beta[i] = complex<double>((double)Beta0[i], 0.0);
   Lscale[i] = (double)Lscale0[i]; Rscale[i] = (double)Rscale0[i];
   Rconde[i] = (double)Rconde0[i]; Rcondv[i] = (double)Rcondv0[i];
  }
  abnrm = (double)abnrm0; bbnrm = (double)bbnrm0;
  delete[] Alphar0; delete[] Alphai0; delete[] Beta0;
  delete[] Lscale0; delete[] Rscale0;
  delete[] Rconde0; delete[] Rcondv0;
  delete[] Work; delete[] Iwork;
 }
 else if (typeid(T) == typeid(double))
 {
  double* Alphar0 = new double[n]; double* Alphai0 = new double[n]; double* Beta0 = new double[n];
  double* Lscale0 = new double[n]; double* Rscale0 = new double[n];
  double abnrm0, bbnrm0;
  double* Rconde0 = new double[n]; double* Rcondv0 = new double[n];
  lwork = 2*n*n+12*n+16;
  double* Work = new double[lwork]; int* Iwork = new int[n+6];
  dggevx_(&balanc, &jobvl, &jobvr, &sense, &n, (double*)this->Elements, &n, (double*)Matrix0.Elements, &n,
          Alphar0, Alphai0, Beta0, (double*)Vl.Elements, &n, (double*)Vr.Elements, &n, &ilo, &ihi,
          Lscale0, Rscale0, &abnrm0, &bbnrm0, Rconde0, Rcondv0, Work, &lwork, Iwork, Bwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "eigenDecompose(Matrix<T>& Matrix0, char balanc, char jobvl, char jobvr, char sense, " <<
                          "vector< complex<double> >& Alpha, vector< complex<double> >& Beta, " <<
                          "Matrix<T>& Vl, Matrix<T>& Vr, int& ilo, int& ihi, " <<
                          "vector<double>& Lscale, vector<double>& Rscale, double& abnrm, double& bbnrm, " <<
                          "vector<double>& Rconde, vector<double>& Rcondv): " <<
           "INFO = " << info << " in LAPACK's DGGEVX." << endl;
   exit(1);
  }
  for (int i = 0; i < n; i++)
  {
   Alpha[i] = complex<double>((double)Alphar0[i], (double)Alphai0[i]);
   Beta[i] = complex<double>((double)Beta0[i], 0.0);
   Lscale[i] = (double)Lscale0[i]; Rscale[i] = (double)Rscale0[i];
   Rconde[i] = (double)Rconde0[i]; Rcondv[i] = (double)Rcondv0[i];
  }
  abnrm = (double)abnrm0; bbnrm = (double)bbnrm0;
  delete[] Alphar0; delete[] Alphai0; delete[] Beta0;
  delete[] Lscale0; delete[] Rscale0;
  delete[] Rconde0; delete[] Rcondv0;
  delete[] Work; delete[] Iwork;
 }
 else if (typeid(T) == typeid(complex<float>))
 {
  complex<float>* Alpha0 = new complex<float>[n]; complex<float>* Beta0 = new complex<float>[n];
  float* Lscale0 = new float[n]; float* Rscale0 = new float[n];
  float abnrm0, bbnrm0;
  float* Rconde0 = new float[n]; float* Rcondv0 = new float[n];
  lwork = 2*n*n+2*n;
  complex<float>* Work = new complex<float>[lwork]; float* Rwork = new float[6*n]; int* Iwork = new int[n+2];
  cggevx_(&balanc, &jobvl, &jobvr, &sense, &n, (complex<float>*)this->Elements, &n,
          (complex<float>*)Matrix0.Elements, &n, Alpha0, Beta0, (complex<float>*)Vl.Elements, &n,
          (complex<float>*)Vr.Elements, &n, &ilo, &ihi, Lscale0, Rscale0, &abnrm0, &bbnrm0, Rconde0, Rcondv0,
          Work, &lwork, Rwork, Iwork, Bwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "eigenDecompose(Matrix<T>& Matrix0, char balanc, char jobvl, char jobvr, char sense, " <<
                          "vector< complex<double> >& Alpha, vector< complex<double> >& Beta, " <<
                          "Matrix<T>& Vl, Matrix<T>& Vr, int& ilo, int& ihi, " <<
                          "vector<double>& Lscale, vector<double>& Rscale, double& abnrm, double& bbnrm, " <<
                          "vector<double>& Rconde, vector<double>& Rcondv): " <<
           "INFO = " << info << " in LAPACK's CGGEVX." << endl;
   exit(1);
  }
  for (int i = 0; i < n; i++)
  {
   Alpha[i] = (complex<double>)Alpha0[i];
   Beta[i] = (complex<double>)Beta0[i];
   Lscale[i] = (double)Lscale0[i]; Rscale[i] = (double)Rscale0[i];
   Rconde[i] = (double)Rconde0[i]; Rcondv[i] = (double)Rcondv0[i];
  }
  abnrm = (double)abnrm0; bbnrm = (double)bbnrm0;
  delete[] Alpha0; delete[] Beta0;
  delete[] Lscale0; delete[] Rscale0;
  delete[] Rconde0; delete[] Rcondv0;
  delete[] Work; delete[] Rwork; delete[] Iwork;
 }
 else if (typeid(T) == typeid(complex<double>))
 {
  complex<double>* Alpha0 = new complex<double>[n]; complex<double>* Beta0 = new complex<double>[n];
  double* Lscale0 = new double[n]; double* Rscale0 = new double[n];
  double abnrm0, bbnrm0;
  double* Rconde0 = new double[n]; double* Rcondv0 = new double[n];
  lwork = 2*n*n+2*n;
  complex<double>* Work = new complex<double>[lwork]; double* Rwork = new double[6*n]; int* Iwork = new int[n+2];
  zggevx_(&balanc, &jobvl, &jobvr, &sense, &n, (complex<double>*)this->Elements, &n,
          (complex<double>*)Matrix0.Elements, &n, Alpha0, Beta0, (complex<double>*)Vl.Elements, &n,
          (complex<double>*)Vr.Elements, &n, &ilo, &ihi, Lscale0, Rscale0, &abnrm0, &bbnrm0, Rconde0, Rcondv0,
          Work, &lwork, Rwork, Iwork, Bwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "eigenDecompose(Matrix<T>& Matrix0, char balanc, char jobvl, char jobvr, char sense, " <<
                          "vector< complex<double> >& Alpha, vector< complex<double> >& Beta, " <<
                          "Matrix<T>& Vl, Matrix<T>& Vr, int& ilo, int& ihi, " <<
                          "vector<double>& Lscale, vector<double>& Rscale, double& abnrm, double& bbnrm, " <<
                          "vector<double>& Rconde, vector<double>& Rcondv): " <<
           "INFO = " << info << " in LAPACK's ZGGEVX." << endl;
   exit(1);
  }
  for (int i = 0; i < n; i++)
  {
   Alpha[i] = (complex<double>)Alpha0[i];
   Beta[i] = (complex<double>)Beta0[i];
   Lscale[i] = (double)Lscale0[i]; Rscale[i] = (double)Rscale0[i];
   Rconde[i] = (double)Rconde0[i]; Rcondv[i] = (double)Rcondv0[i];
  }
  abnrm = (double)abnrm0; bbnrm = (double)bbnrm0;
  delete[] Alpha0; delete[] Beta0;
  delete[] Lscale0; delete[] Rscale0;
  delete[] Rconde0; delete[] Rcondv0;
  delete[] Work; delete[] Rwork; delete[] Iwork;
 }
 delete[] Bwork;
}

template<class T> void Matrix<T>::linearLeastSquares(double rcond, Matrix<T>& b, unsigned int& rank,
                                                     vector<double>& S)
{
#ifdef DEBUG
 if ((this->size == 0) || (this->Shape[0] < this->Shape[1]) || (b.size == 0) ||
     (b.Shape[0] != this->Shape[0]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "linearLeastSquares(double rcond, Matrix<T>& b, unsigned int& rank, vector<double>& S): " <<
          "((this->size == 0) || (this->Shape[0] < this->Shape[1]) || (b.size == 0) || " <<
           "(b.Shape[0] != this->Shape[0]))." << endl;
  exit(1);
 }
#endif
 int m = this->Shape[0], n = this->Shape[1], nrhs = b.Shape[1];
 S = vector<double>(min(m,n));
 int rank0, lwork, info;
 int* Iwork = new int[50*min(m,n)];
 if (typeid(T) == typeid(float))
 {
  float* S0 = new float[min(m,n)];
  float rcond0 = rcond;
  lwork = (200+nrhs)*n+1000;
  float* Work = new float[lwork];
  sgelsd_(&m, &n, &nrhs, (float*)this->Elements, &m, (float*)b.Elements, &m, S0, &rcond0, &rank0, Work,
          &lwork, Iwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "linearLeastSquares(double rcond, Matrix<T>& b, unsigned int& rank, vector<double>& S): " <<
           "INFO = " << info << " in LAPACK's SGELSD." << endl;
   exit(1);
  }
  rank = rank0;
  for (int i = 0; i < min(m,n); i++)
  {
   S[i] = (double)S0[i];
  }
  delete[] S0;
  delete[] Work;
 }
 else if (typeid(T) == typeid(double))
 {
  double* S0 = new double[min(m,n)];
  double rcond0 = rcond;
  lwork = (200+nrhs)*n+1000;
  double* Work = new double[lwork];
  dgelsd_(&m, &n, &nrhs, (double*)this->Elements, &m, (double*)b.Elements, &m, S0, &rcond0, &rank0, Work,
          &lwork, Iwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "linearLeastSquares(double rcond, Matrix<T>& b, unsigned int& rank, vector<double>& S): " <<
           "INFO = " << info << " in LAPACK's DGELSD." << endl;
   exit(1);
  }
  rank = rank0;
  for (int i = 0; i < min(m,n); i++)
  {
   S[i] = (double)S0[i];
  }
  delete[] S0;
  delete[] Work;
 }
 else if (typeid(T) == typeid(complex<float>))
 {
  float* S0 = new float[min(m,n)];
  float rcond0 = rcond;
  lwork = (10+nrhs)*n;
  complex<float>* Work = new complex<float>[lwork];
  float* Rwork = new float[200*n+100*nrhs+1000];
  cgelsd_(&m, &n, &nrhs, (complex<float>*)this->Elements, &m, (complex<float>*)b.Elements, &m, S0, &rcond0,
          &rank0, Work, &lwork, Rwork, Iwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "linearLeastSquares(double rcond, Matrix<T>& b, unsigned int& rank, vector<double>& S): " <<
           "INFO = " << info << " in LAPACK's CGELSD." << endl;
   exit(1);
  }
  rank = rank0;
  for (int i = 0; i < min(m,n); i++)
  {
   S[i] = (double)S0[i];
  }
  delete[] S0;
  delete[] Work;
  delete[] Rwork;
 }
 else if (typeid(T) == typeid(complex<double>))
 {
  double* S0 = new double[min(m,n)];
  double rcond0 = rcond;
  lwork = (10+nrhs)*n;
  complex<double>* Work = new complex<double>[lwork];
  double* Rwork = new double[200*n+100*nrhs+1000];
  zgelsd_(&m, &n, &nrhs, (complex<double>*)this->Elements, &m, (complex<double>*)b.Elements, &m, S0,
          &rcond0, &rank0, Work, &lwork, Rwork, Iwork, &info);
  if (info != 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Matrix<T>::" <<
           "linearLeastSquares(double rcond, Matrix<T>& b, unsigned int& rank, vector<double>& S): " <<
           "INFO = " << info << " in LAPACK's ZGELSD." << endl;
   exit(1);
  }
  rank = rank0;
  for (int i = 0; i < min(m,n); i++)
  {
   S[i] = (double)S0[i];
  }
  delete[] S0;
  delete[] Work;
  delete[] Rwork;
 }
 delete[] Iwork;
}

template<class T> void Matrix<T>::computeGroundState(T& e, vector<T>& V) const
{
 int m = this->Shape[0]; int n = this->Shape[1];
#ifdef DEBUG
 if ((m != n) || (this->Type != "hermitian") || (V.size() != m))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "computeGroundState(T& e, vector<T>& V) const: " <<
          "((m != n) || (this->Type != hermitian) || (V.size() != m))." << endl;
  exit(1);
 }
 else if ((typeid(T) != typeid(float)) && (typeid(T) != typeid(double)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "computeGroundState(T& e, vector<T>& V) const: " <<
          "This function is not implemented yet for complex matrices." << endl;
  exit(1);
 }
#endif
 int maxitr = 300; int maxn = n; int maxnev = 1;
 int ncv, maxncv;
 if (n < 25)
 {
  ncv = n;
 }
 else
 {
  ncv = 25;
 }
 maxncv = ncv;
 int itrdone, matvecmultdone, infoaupd, infoeupd;
 T* e0 = new T;
 T* V0 = new T[n];
 if (typeid(T) == typeid(float))
 {
  float tol = 0.0;
  ssyevar_(&n, (float*)this->Elements, &tol, &maxitr, &maxn, &maxnev, &ncv, &maxncv, &itrdone,
           &matvecmultdone, &infoaupd, &infoeupd, (float*)e0, (float*)V0);
 }
 else if (typeid(T) == typeid(double))
 {
  double tol = 0.0;
  dsyevar_(&n, (double*)this->Elements, &tol, &maxitr, &maxn, &maxnev, &ncv, &maxncv, &itrdone,
           &matvecmultdone, &infoaupd, &infoeupd, (double*)e0, (double*)V0);
 }
 if (((infoaupd == 0) && (infoeupd == 0)) || ((infoaupd == 1) && (infoeupd == 0)))
 {
  e = *e0;
  for (int i = 0; i < n; i++)
  {
   V[i] = V0[i];
  }
 }
 else
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "computeGroundState(T& e, vector<T>& V) const: " <<
          "info = " << infoaupd << " in ARPACK's XSAUPD and " <<
          "info = " << infoeupd << " in ARPACK's XSEUPD." << endl;
  exit(1);
 }
 delete e0;
 delete[] V0;
}

template<class T> void Matrix<T>::exponentialize(T x, Matrix<T>& MatrixExponent) const
{
#ifdef DEBUG
 if ((this->size == 0) || (this->Shape[0] != this->Shape[1]) ||
     (this->Shape[0] != MatrixExponent.Shape[0]) ||
     (MatrixExponent.Shape[0] != MatrixExponent.Shape[1]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "exponentialize(T x, Matrix<T>& MatrixExponent) const: " <<
          "((this->size == 0) || (this->Shape[0] != this->Shape[1]) || " <<
            "(this->Shape[0] != MatrixExponent.Shape[0]) || " <<
            "(MatrixExponent.Shape[0] != MatrixExponent.Shape[1]))." << endl;
  exit(1);
 }
#endif
 unsigned int dim = this->Shape[0];
 Matrix<T> Vr(dim, dim);
 T element;
 if (this->Type == "hermitian")
 {
  vector<double> W(dim);
  this->eigenDecompose(W, Vr);
  for (int i = 0; i < dim; i++)
  {
   for (int j = 0; j < dim; j++)
   {
    element = 0.0;
    for (int k = 0; k < dim; k++)
    {
     element += exp(x * T(W[k])) * Vr(i, k) * MathAuxiliary::complexConjugate(Vr(j, k));
    }
    MatrixExponent(i, j) = element;
   }
  }
 }
 else if (this->Type == "general")
 {
  cerr << "The following function is not implemented yet for Matrices of Type general: " <<
          "template<class T> void Matrix<T>::" <<
          "exponentialize(T x, Matrix<T>& MatrixExponent) const." << endl;
  exit(1);
 }
}

template<class T> void Matrix<T>::twoBodyExponentSVD(T delta, unsigned int d0,
                                                     unsigned int d1, Tensor<T>& TensorLeft,
                                                     Tensor<T>& TensorRight) const
{
 vector<unsigned int> TensorLeftShape, TensorRightShape;
 TensorLeft.getShape(TensorLeftShape);
 TensorRight.getShape(TensorRightShape);
 unsigned int numSV = min(d0*d0, d1*d1);
#ifdef DEBUG
 if ((this->Type != "hermitian") || (this->Shape[0] != d0*d1) || (this->Shape[1] != d0*d1) ||
     (TensorLeftShape[0] != 1) || (TensorLeftShape[1] != numSV) || (TensorLeftShape[2] != d0) ||
     (TensorLeftShape[3] != d0) || (TensorRightShape[0] != numSV) || (TensorRightShape[1] != 1) ||
     (TensorRightShape[2] != d1) || (TensorRightShape[3] != d1))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Matrix<T>::" <<
          "twoBodyExponentSVD(T delta, unsigned int d0, unsigned int d1, Tensor<T>& TensorLeft, " <<
                             "Tensor<T>& TensorRight) const: " <<
          "This Matrix or the input arguments have an incorrect form." << endl;
  exit(1);
 }
#endif
 vector<unsigned int> Shape0(4);
 Shape0[0] = d0; Shape0[1] = d1; Shape0[2] = d0; Shape0[3] = d1;
 Tensor<T> TwoBodyOperator(Shape0);
 vector<unsigned int> Index(4);
 T element;
 vector<double> W(d0*d1);
 Matrix<T> Vr(d0*d1, d0*d1);
 this->eigenDecompose(W, Vr);
 for (int i = 0; i < d0*d1; i++)
 {
  Index[2] = i % d0; Index[3] = i / d0;
  for (int j = 0; j < d0*d1; j++)
  {
   Index[0] = j % d0; Index[1] = j / d0;
   element = 0.0;
   for (int k = 0; k < d0*d1; k++)
    element += exp(-delta * T(W[k])) * Vr(i, k) * MathAuxiliary::complexConjugate(Vr(j, k));
   TwoBodyOperator.set(Index, element);
  }
 }
 TwoBodyOperator.twoBodyOperatorSVD(TensorLeft, TensorRight);
}

template<class T> void Matrix<T>::computeLowestEigenstates(const vector<T>& AD, const vector<int>& IndexD,
                                                           const vector<T>& AOD, const vector<int>& IndexODi,
                                                           const vector<int>& IndexODj, int n, int nev, int ncv,
                                                           vector<T>& Evals, Matrix<T>& Evecs)
{
 int m1 = AD.size(); int m2 = AOD.size();
#ifdef DEBUG
 if ((IndexD.size() != m1) || (IndexODi.size() != m2) || (IndexODj.size() != m2) || (nev >= ncv))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void Matrix<T>::" <<
          "computeLowestEigenstates(const vector<T>& AD, const vector<int>& IndexD, " <<
                                   "const vector<T>& AOD, const vector<int>& IndexODi, const vector<int>& IndexODj, " <<
                                   "int n, int nev, int ncv, vector<T>& Evals, Matrix<T>& Evecs): " <<
          "((IndexD.size() != AD.size()) || (IndexODi.size() != AOD.size()) || (IndexODj.size() != AOD.size()) || " <<
           "(nev >= ncv))." << endl;
  exit(1);
 }
 else if ((typeid(T) != typeid(float)) && (typeid(T) != typeid(double)))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void Matrix<T>::" <<
          "computeLowestEigenstates(const vector<T>& AD, const vector<int>& IndexD, " <<
                                   "const vector<T>& AOD, const vector<int>& IndexODi, const vector<int>& IndexODj, " <<
                                   "int n, int nev, int ncv, vector<T>& Evals, Matrix<T>& Evecs): " <<
          "This function is not implemented yet for complex matrices." << endl;
  exit(1);
 }
#endif
 int maxitr = 300; int maxn = n; int maxnev = nev; int maxncv = ncv;
 int itrdone, matvecmultdone, infoaupd, infoeupd;
 Evals = vector<T>(nev);
 Evecs = Matrix<T>(n, nev);
 if (typeid(T) == typeid(float))
 {
  float tol = 0.0;
  ssyevar2_(&m1, (float*)&AD[0], (int*)&IndexD[0], &m2, (float*)&AOD[0], (int*)&IndexODi[0], (int*)&IndexODj[0], &tol, &maxitr, &n,
            &maxn, &nev, &maxnev, &ncv, &maxncv, &itrdone, &matvecmultdone, &infoaupd, &infoeupd, (float*)&Evals[0],
            (float*)Evecs.Elements);
 }
 else if (typeid(T) == typeid(double))
 {
  double tol = 0.0;
  dsyevar2_(&m1, (double*)&AD[0], (int*)&IndexD[0], &m2, (double*)&AOD[0], (int*)&IndexODi[0], (int*)&IndexODj[0], &tol, &maxitr, &n,
            &maxn, &nev, &maxnev, &ncv, &maxncv, &itrdone, &matvecmultdone, &infoaupd, &infoeupd, (double*)&Evals[0],
            (double*)Evecs.Elements);
 }
 if ((infoaupd != 0) || (infoeupd != 0))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void Matrix<T>::" <<
          "computeLowestEigenstates(const vector<T>& AD, const vector<int>& IndexD, " <<
                                   "const vector<T>& AOD, const vector<int>& IndexODi, const vector<int>& IndexODj, " <<
                                   "int n, int nev, int ncv, vector<T>& Evals, Matrix<T>& Evecs): " <<
          "info = " << infoaupd << " in ARPACK's XSAUPD and " <<
          "info = " << infoeupd << " in ARPACK's XSEUPD." << endl;
  exit(1);
 }
}

template<class T> void Matrix<T>::computeArnoldi(const Tensor<T>& TensorLeft, const Tensor<T>& TensorMiddle,
                                                 const Tensor<T>& TensorRight, unsigned int which,
                                                 Tensor<T>& Tensor0)
{
 vector<unsigned int> ShapeLeft, ShapeMiddle, ShapeRight;
 TensorLeft.getShape(ShapeLeft); TensorMiddle.getShape(ShapeMiddle); TensorRight.getShape(ShapeRight);
 unsigned int D1Left = ShapeLeft[3], D2Left = ShapeLeft[4], d = ShapeMiddle[2], D1Right = ShapeRight[3],
              D2Right = ShapeRight[4];
#ifdef DEBUG
 if ((ShapeLeft[0] != 1) || (ShapeLeft[1] != 1) || (ShapeLeft[2] != 1) || (ShapeLeft[5] != D1Left) ||
     (ShapeMiddle[0] != D2Left) || (ShapeMiddle[1] != D2Right) || (ShapeMiddle[3] != d) ||
     (ShapeRight[0] != 1) || (ShapeRight[1] != 1) || (ShapeRight[2] != 1) || (ShapeRight[5] != D1Right))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void Matrix<T>::" <<
          "computeArnoldi(const Tensor<T>& TensorLeft, const Tensor<T>& TensorMiddle, " <<
                         "const Tensor<T>& TensorRight, unsigned int which, Tensor<T>& Tensor0): " <<
          "TensorLeft, TensorMiddle or TensorRight do not have the correct Shapes." << endl;
  exit(1);
 }
 else if ((typeid(T) != typeid(float)) && (typeid(T) != typeid(double)))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void Matrix<T>::" <<
          "computeArnoldi(const Tensor<T>& TensorLeft, const Tensor<T>& TensorMiddle, " <<
                         "const Tensor<T>& TensorRight, unsigned int which, Tensor<T>& Tensor0): " <<
          "This function is not implemented yet for complex tensors." << endl;
  exit(1);
 }
#endif
 vector<unsigned int> Shape0(3);
 Shape0[0] = D1Left; Shape0[1] = D1Right; Shape0[2] = d;
 Tensor0 = Tensor<T>(Shape0);
 Tensor<T> Tensor1, Tensor2;
 vector<unsigned int> Indices1(1), Indices2(1), Indices12(2), Indices22(2);
 vector<unsigned int> Order(9);
 Order[0] = 3; Order[1] = 8; Order[2] = 4; Order[3] = 0; Order[4] = 1; Order[5] = 2; Order[6] = 5;
 Order[7] = 6; Order[8] = 7;
// parameters for ARPACK's xsaupd:
 int ido = 0; char bmat = 'I'; int n = D1Left*d*D1Right; char Which[] = {'S', 'A'}; int nev = which + 1;
 T tol = 0.0; T* Resid = new T[n]; int ncv = min(25, n); T* V = new T[n*ncv]; int ldv = n;
 int Iparam[11];
 int ishift = 1, mxiter = 300, mode = 1;
 Iparam[0] = ishift; Iparam[2] = mxiter; Iparam[6] = mode;
 int Ipntr[11]; T* Workd = new T[3*n]; T* Workl = new T[ncv*(ncv+8)]; int lworkl = ncv*(ncv+8);
 int info = 0;
// parameters for ARPACK's xseupd:
 int rvec = true; char howmny = 'A'; int* Select = new int[ncv]; T* D = new T[nev]; T sigma; int ierr;
 if (typeid(T) == typeid(float))
 {
  ssaupd_(&ido, &bmat, &n, Which, &nev, (float*)&tol, (float*)Resid, &ncv, (float*)V, &ldv, Iparam, Ipntr,
          (float*)Workd, (float*)Workl, &lworkl, &info);
  while ((ido == -1) || (ido == 1))
  {
/* Begin of the user defined matrix-vector multiplication where the action of the Hamiltonian on
   Workd[Ipntr[0]-1] is returned in Workd[Ipntr[1]-1]: */
   for (int i = 0; i < n; i++)
    Tensor0.set(i, Workd[Ipntr[0]-1+i]);
   Tensor1 = TensorLeft; Tensor2 = Tensor0;
   Indices1[0] = 3; Indices2[0] = 0;
   Tensor1.contract(Indices1, Tensor2, Indices2);
   Tensor2 = TensorMiddle;
   Indices12[0] = 3; Indices12[1] = 6; Indices22[0] = 0; Indices22[1] = 2;
   Tensor1.contract(Indices12, Tensor2, Indices22);
   Tensor2 = TensorRight;
   Indices12[0] = 4; Indices12[1] = 5; Indices22[0] = 3; Indices22[1] = 4;
   Tensor1.contract(Indices12, Tensor2, Indices22);
   Tensor1.permute(Order);
   for (int i = 0; i < n; i++)
    Workd[Ipntr[1]-1+i] = Tensor1.get(i);
/* End of the user defined matrix-vector multiplication */
   ssaupd_(&ido, &bmat, &n, Which, &nev, (float*)&tol, (float*)Resid, &ncv, (float*)V, &ldv, Iparam, Ipntr,
           (float*)Workd, (float*)Workl, &lworkl, &info);
  }
  if (info < 0)
  {
   cerr << "Program terminated because of error in static function " <<
           "template<class T> void Matrix<T>::" <<
           "computeArnoldi(const Tensor<T>& TensorLeft, const Tensor<T>& TensorMiddle, " <<
                          "const Tensor<T>& TensorRight, unsigned int which, Tensor<T>& Tensor0): " << endl;
   cerr << "Error with ssaupd, info = " << info << endl;
   cerr << "Check documentation in ssaupd." << endl;
   exit(1);
  }
  else
  {
   sseupd_(&rvec, &howmny, Select, (float*)D, (float*)V, &ldv, (float*)&sigma, &bmat, &n, Which, &nev,
           (float*)&tol, (float*)Resid, &ncv, (float*)V, &ldv, Iparam, Ipntr, (float*)Workd, (float*)Workl,
           &lworkl, &ierr);
   if (ierr != 0)
   {
    cerr << "Program terminated because of error in static function " <<
            "template<class T> void Matrix<T>::" <<
            "computeArnoldi(const Tensor<T>& TensorLeft, const Tensor<T>& TensorMiddle, " <<
                           "const Tensor<T>& TensorRight, unsigned int which, " <<
                           "Tensor<T>& Tensor0): " << endl;
    cerr << "Error with sseupd, info = " << ierr << endl;
    cerr << "Check the documentation of sseupd." << endl;
    exit(1);
   }
   else
   {
    for (int i = 0; i < n; i++)
     Tensor0.set(i, V[which*n+i]);
   }
  }
 }
 else if (typeid(T) == typeid(double))
 {
  dsaupd_(&ido, &bmat, &n, Which, &nev, (double*)&tol, (double*)Resid, &ncv, (double*)V, &ldv, Iparam,
          Ipntr, (double*)Workd, (double*)Workl, &lworkl, &info);
  while ((ido == -1) || (ido == 1))
  {
/* Begin of the user defined matrix-vector multiplication where the action of the Hamiltonian on
   Workd[Ipntr[0]-1] is returned in Workd[Ipntr[1]-1]: */
   for (int i = 0; i < n; i++)
    Tensor0.set(i, Workd[Ipntr[0]-1+i]);
   Tensor1 = TensorLeft; Tensor2 = Tensor0;
   Indices1[0] = 3; Indices2[0] = 0;
   Tensor1.contract(Indices1, Tensor2, Indices2);
   Tensor2 = TensorMiddle;
   Indices12[0] = 3; Indices12[1] = 6; Indices22[0] = 0; Indices22[1] = 2;
   Tensor1.contract(Indices12, Tensor2, Indices22);
   Tensor2 = TensorRight;
   Indices12[0] = 4; Indices12[1] = 5; Indices22[0] = 3; Indices22[1] = 4;
   Tensor1.contract(Indices12, Tensor2, Indices22);
   Tensor1.permute(Order);
   for (int i = 0; i < n; i++)
    Workd[Ipntr[1]-1+i] = Tensor1.get(i);
/* End of the user defined matrix-vector multiplication */
   dsaupd_(&ido, &bmat, &n, Which, &nev, (double*)&tol, (double*)Resid, &ncv, (double*)V, &ldv, Iparam,
           Ipntr, (double*)Workd, (double*)Workl, &lworkl, &info);
  }
  if (info < 0)
  {
   cerr << "Program terminated because of error in static function " <<
           "template<class T> void Matrix<T>::" <<
           "computeArnoldi(const Tensor<T>& TensorLeft, const Tensor<T>& TensorMiddle, " <<
                          "const Tensor<T>& TensorRight, unsigned int which, Tensor<T>& Tensor0): " << endl;
   cerr << "Error with dsaupd, info = " << info << endl;
   cerr << "Check documentation in dsaupd." << endl;
   exit(1);
  }
  else
  {
   dseupd_(&rvec, &howmny, Select, (double*)D, (double*)V, &ldv, (double*)&sigma, &bmat, &n, Which, &nev,
           (double*)&tol, (double*)Resid, &ncv, (double*)V, &ldv, Iparam, Ipntr, (double*)Workd,
           (double*)Workl, &lworkl, &ierr);
   if (ierr != 0)
   {
    cerr << "Program terminated because of error in static function " <<
            "template<class T> void Matrix<T>::" <<
            "computeArnoldi(const Tensor<T>& TensorLeft, const Tensor<T>& TensorMiddle, " <<
                           "const Tensor<T>& TensorRight, unsigned int which, " <<
                           "Tensor<T>& Tensor0): " << endl;
    cerr << "Error with dseupd, info = " << ierr << endl;
    cerr << "Check the documentation of dseupd." << endl;
    exit(1);
   }
   else
   {
    for (int i = 0; i < n; i++)
     Tensor0.set(i, V[which*n+i]);
   }
  }
 }
 delete[] Resid;
 delete[] V;
 delete[] Workd;
 delete[] Workl;
 delete[] Select;
 delete[] D;
}
