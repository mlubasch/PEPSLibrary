/// Template class MPO implements Matrix Product Operators.
/** The template class MPO implements Matrix Product Operators as a generalization of MPS. MPO consist
    of Tensors on a linear chain.
    \param BC string, the boundary conditions of the MPO
    \param N unsigned int, the number of tensors of the MPO
    \param d vector<unsigned int>, the physical dimensions of the MPO
    \param D unsigned int, the maximal virtual bond dimension of the MPO
    \param Tensors Tensor<T>*, the tensors of the MPO
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

// Declaration of template classes:
template<class T> class Matrix;
template<class T> class MPS;
template<class T> class MPO;

// Declaration of friend functions:
template<class T> void getMPOFromLocalOperator(string BC, unsigned int N, unsigned int x,
                                               const Matrix<T>& O, MPO<T>& MPO0);
template<class T> void getMPOFromSumLocalOperator(string BC, unsigned int N, const Matrix<T>& O,
                                                  MPO<T>& MPO0);
template<class T> T expectationValueMPO(const MPS<T>& MPS0, const MPO<T>& MPO0);
template<class T> T scalarProductMPOMPS(const MPS<T>& MPS0, const MPO<T>& MPO0, const MPS<T>& MPS1);
template<class T> T contractReverseMPOMPS(const MPS<T>& MPS0, const MPO<T>& MPO0, const MPS<T>& MPS1);
template<class T> void exactMultiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, MPS<T>& MPS1);
template<class T> double distanceMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, const MPS<T>& MPS1);
template<class T> void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                      unsigned int maxNumSweeps, double& errorAchieved,
                                      unsigned int& numSweepsDone, MPS<T>& MPS1);
template<class T> void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                      unsigned int maxNumSweeps, double& errorAchieved,
                                      unsigned int& numSweepsDone, MPS<T>& MPS1,
                                      double cutoff, unsigned int mode, unsigned int d2,
                                      double alphaStart, double x, double precision, unsigned int maxNumIter);
template<class T> void multiplyMPDOMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                          unsigned int maxNumSweeps, double& errorAchieved,
                                          unsigned int& numSweepsDone, MPS<T>& MPS1);
template<class T> void multiplyEvenOddMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, MPS<T>& MPS1);

template<class T> class MPO
{
 public:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  MPO();

/// Constructor for MPO with specific BC, N, d and D.
/** This constructor initializes a MPO of a specific form. The virtual bond dimension at each site is
    the maximal virtual bond dimension D.
    \param BC input: string, the boundary conditions
    \param N input: unsigned int, the number of tensors
    \param d input: const vector<unsigned int>&, the physical dimensions
    \param D input: unsigned int, the maximal virtual bond dimension */
  MPO(string BC, unsigned int N, const vector<unsigned int>& d, unsigned int D);

/// Constructor for MPO with specific BC, N, d and D.
/** This constructor initializes a MPO of a specific form with all physical dimensions equal to d. The
    virtual bond dimension at each site is the maximal virtual bond dimension D.
    \param BC input: string, the boundary conditions
    \param N input: unsigned int, the number of tensors
    \param d input: unsigned int, the physical dimension
    \param D input: unsigned int, the maximal virtual bond dimension */
  MPO(string BC, unsigned int N, unsigned int d, unsigned int D);

/// Standard copy constructor.
/** The standard copy constructor copies the input MPO into this.
    \param MPO0 input: const MPO<T>&, to be copied into this
    \sa MPO<T>& operator=(const MPO<T>& MPO0) */
  MPO(const MPO<T>& MPO0);

/// Standard destructor.
/** The standard destructor deletes the elements of the MPO. */
  ~MPO();

/// Assigns MPO to this.
/** The operator= allows to assign a MPO to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side MPO.
    \param MPO0 input: const MPO<T>&, to be copied into this
    \return MPO<T>&, a reference to the new this
    \sa MPO(const MPO<T>& MPO0) */
  MPO<T>& operator=(const MPO<T>& MPO0);

/// Returns boundary conditions of MPO.
/** The returned boundary conditions are either "open" or "periodic".
    \param BC0 output: string&, the boundary conditions of the MPO */
  void getBC(string& BC0) const { BC0 = this->BC; }

/// Returns number of tensors of MPO.
/** This function returns the number of tensors building up the MPO.
    \return unsigned int, the number of tensors of the MPO */
  unsigned int getN() const { return this->N; }

/// Sets physical dimension.
/** This function sets the physical dimension of this MPO at position.
    \param position input: unsigned int, the desired position
    \param d0 input: unsigned int, the physical dimension to be written at position */
  void setd(unsigned int position, unsigned int d0);

/// Returns physical dimensions of MPO.
/** This function returns the physical dimensions of the MPO.
    \param d0 output: vector<unsigned int>&, the physical dimensions of the MPO */
  void getd(vector<unsigned int>& d0) const { d0 = this->d; }

/// Sets maximal virtual bond dimension.
/** This function sets the maximal virtual bond dimension of this MPO.
    If the maximal virtual bond dimension is decreased, i.e. if D0 < this->D, then the dispensable
    elements will be simply discarded. If it is increased, i.e. if D0 > this->D, then the new elements
    will be set according to Filling either as "zeroes" or "random". Filling has the default value "zeroes".
    \param D0 input: unsigned int, the new maximal virtual bond dimension
    \param Filling input: const string&, the filling for D0 > this->D, has default value "zeroes", must
                          be either "zeroes" or "random" */
  void setD(unsigned int D0, const string& Filling = "zeroes");

/// Returns maximal virtual bond dimension of MPO.
/** This function returns the maximal virtual bond dimension of the MPO.
    \return unsigned int, the maximal virtual bond dimension of the MPO */
  unsigned int getD() const { return this->D; }

/// Sets Tensor at position.
/** Given a position the corresponding tensor in Tensors is set.
    \param position input: unsigned int, the desired position
    \param Tensor0 input: const Tensor<T>&, to be written at position */
  void set(unsigned int position, const Tensor<T>& Tensor0);

/// Returns Tensor at position.
/** Given a position the corresponding tensor in Tensors is returned.
    \param position input: unsigned int, the desired position
    \param Tensor0 output: Tensor<T>&, a copy of the tensor at position */
  void get(unsigned int position, Tensor<T>& Tensor0) const;

/// Returns this MPO as matrix.
/** This MPO is returned as a matrix in the standard basis.
    \param Matrix0 output: Matrix<T>&, this MPO as a matrix,
                           must fulfill Matrix0.dim0==d^{N} and Matrix0.dim1==d^{N} */
  void getMatrix(Matrix<T>& Matrix0) const;

/// This is filled with random entries.
/** This MPO is filled with uniformly distributed random entries. LAPACK's XLARNV is used.
    \param Seed input: const vector<unsigned int>&, must fulfill (Seed.size()==this->N) and contains the seed
                       values for each tensor
    \sa void fillZeroes() */
  void fillRandomly(const vector<unsigned int>& Seed);

/// Fills this MPO with zeroes.
/** This MPO is filled with zeroes.
    \sa void fillRandomly(const vector<unsigned int>& Seed) */
  void fillZeroes();

/// Normalizes this MPO.
/** This function normalizes this MPO by normalizing each tensor with
       T Tensor<T>::normalize()   .
    \return T, the product of the absolute values of the largest elements of all tensors
    \sa T Tensor<T>::normalize() */
  T simplifiedNormalize();

/// Transposes this MPO.
/** This MPO is transposed by swapping the physical indices of each tensor. */
  void transpose();

/// Returns adjoint of this MPO.
/** The adjoint of this MPO is returned in the input parameter.
    \param MPO0 output: MPO<T>&, the adjoint of this MPO */
  void adjoint(MPO<T>& MPO0) const;

/// Returns MPS representation of this MPO.
/** This function returns the MPS representation of this MPO.
    Hereby the physical index of the MPS comprises the two physical indices of this MPO in their
    standard order.
    \param MPS0 output: MPS<T>&, the MPS representation of this MPO */
  void getMPS(MPS<T>& MPS0) const;

/// Computes canonical form of this MPO.
/** This function computes the canonical form of this MPO and stores it in MPO0 and Lambdas.
    After obtaining the MPS representation of this MPO with
       void getMPS(MPS<T>& MPS0) const   ,
    this function uses
       void MPS<T>::canonicalize(double eps, MPS<T>& MPS0, vector< Matrix<T> >& Lambdas) const   .
    First the MPS is brought into normal form from right to left using
       void bringIntoNormalForm()   .
    Then each tensor is normalized with
       void normalizeTensor(unsigned int position, const string& Direction, Matrix<T>& X)
    from left to right, and X is singular value decomposed.
    \param eps input: double, the cutoff for the pseudoinverse, elements <= eps are assumed to be zero
    \param MPO0 output: MPO<T>&, the desired MPO in canonical form
    \param Lambdas output: vector< Matrix<T> >&, the lambdas
    \sa void getMPS(MPS<T>& MPS0) const
    \sa void MPS<T>::canonicalize(double eps, MPS<T>& MPS0, vector< Matrix<T> >& Lambdas) const
    \sa void decanonicalize(const vector< Matrix<T> >& Lambdas, MPO<T>& MPO0) const */
  void canonicalize(double eps, MPO<T>& MPO0, vector< Matrix<T> >& Lambdas) const;

/// Computes standard MPO from canonical form.
/** Given this MPO in canonical form with lambda-matrices Lambdas, its standard form, MPO0, is obtained by
    contracting each lambda with the tensor to its right.
    \param Lambdas input: const vector< Matrix<T> >&, the lambdas
    \param MPO0 output: MPO<T>&, the resulting MPO in standard form
    \sa void canonicalize(double eps, MPO<T>& MPO0, vector< Matrix<T> >& Lambdas) const */
  void decanonicalize(const vector< Matrix<T> >& Lambdas, MPO<T>& MPO0) const;

/// This is multiplied with another MPO exactly.
/** This MPO is multiplied with the first input MPO MPO0 exactly and the result is stored in the
    second input MPO MPO1.
    \param MPO0 input: const MPO<T>&, the MPO to be right-multiplied with this MPO
    \param MPO1 output: MPO<T>&, the resulting MPO */
  void multiply(const MPO<T>& MPO0, MPO<T>& MPO1) const;

/// Computes distance between MPO-MPO product and another MPO.
/** This function computes the norm of the difference between this MPO multiplied with the
    first input MPO MPO0 and the second input MPO MPO1:
       || this->MPO*MPO0 - MPO1 ||   .
    \param MPO0 input: const MPO<T>&, the MPO to be right-multiplied with this MPO
    \param MPO1 input: const MPO<T>&, the approximating MPO
    \return double, the resulting distance */
  double distance(const MPO<T>& MPO0, const MPO<T>& MPO1) const;

/// Multiplies this MPO with MPO0 and approximates result by MPO1.
/** This function approximates the product of this MPO with the input MPO MPO0 by a MPO MPO1 with
    a given bond dimension.
    It uses
       friend void multiplyMPOMPS<>(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                    unsigned int maxNumSweeps, double& errorAchieved,
                                    unsigned int& numSweepsDone, MPS<T>& MPS1)   .
    If (eps >= 0), this function evaluates the error of the approximation after each sweep and stops if the change
    in the error is below eps, i.e. in sweep n: eps >= abs(error(n) - error(n-1)). The function always stops if the
    number of sweeps exceeds maxNumSweeps. Hereby the error is defined as
       error := ||this->MPO*MPO0-MPO1||/||this->MPO*MPO0||   .
    If (eps < 0.0), the error is never computed and maxNumSweeps sweeps are performed.
    On input, if (MPO1.D >= this->D*MPO0.D) then the exact multiplication is done.
    On output, errorAchieved is the achieved approximation error, numSweepsDone is the number of sweeps done and
    MPO1 is the approximating MPO.
    \param MPO0 input: const MPO<T>&, the MPO on which this MPO acts
    \param eps input: double, if (eps >= 0.0) the convergence precision, i.e. convergence in sweep n if
                      eps >= abs(error(n) - error(n-1)),
                      if (eps < 0.0) no error computation is done
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param errorAchieved output: double&, the achieved approximation error, if (eps < 0.0) not addressed
    \param numSweepsDone output: unsigned int&, the number of sweeps done
    \param MPO1 input/output: MPO<T>&, on input the initial approximating MPO, on output the resulting
                              MPO, must have the correct form */
  void multiply(const MPO<T>& MPO0, double eps, unsigned int maxNumSweeps, double& errorAchieved,
                unsigned int& numSweepsDone, MPO<T>& MPO1) const;

/// Sets d2 of this MPDO.
/** Assuming this MPO to be a Matrix Product Density Operator (MPDO), this function sets the purification bond d2.
    It uses
       void canonicalize(double eps, MPO<T>& MPO0, vector< Matrix<T> >& Lambdas) const
    to obtain the canonical form of this MPDO. Then each purification bond is reduced independently to d2 by locally
    projecting the one-site reduced density matrix onto the subspace of largest eigenvalues as explained in
       [L. Wang and F. Verstraete, arXiv:1110.4362]   .
    On output, Eigenvalues captures the eigenvalues of the one-site reduced density matrices.
    \param d2 input: unsigned int, the new purification bond
    \param eps input: double, the cutoff for the pseudoinverse, elements <= eps are assumed to be zero
    \param Eigenvalues output: vector< vector<double> >&, the eigenvalues of the one-site reduced density matrices */
  void setd2(unsigned int d2, double eps, vector< vector<double> >& Eigenvalues);

/// Sets d2 of this MPDO.
/** Assuming this MPO to be a Matrix Product Density Operator (MPDO), this function sets the purification bond d2.
    If the purification bond is decreased, it uses
       friend void multiplyMPOMPS<>(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                    unsigned int maxNumSweeps, double& errorAchieved,
                                    unsigned int& numSweepsDone, MPS<T>& MPS1,
                                    double cutoff, unsigned int mode, unsigned int d2)   .
    \param d2 input: unsigned int, the new purification bond d2,
                     if d2 is larger than the current purification bond then the tensors are enlarged and the new
                     elements are zeroes,
                     if d2 is smaller than the current purification bond then the new tensors are found via
                     multiplyMPOMPS
    \param cutoff input: double, the cutoff for the pseudoinverse
    \param eps input: double, if (eps >= 0.0) the convergence precision, i.e. convergence in sweep n if
                      eps >= abs(error(n)-error(n-1)),
                      if (eps < 0.0) no error computation is done
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param errorAchieved output: double&, the achieved approximation error, if (eps < 0.0) not addressed
    \param numSweepsDone output: unsigned int&, the number of sweeps done */
  void setd2(unsigned int d2, double cutoff, double eps, unsigned int maxNumSweeps, double& errorAchieved,
             unsigned int& numSweepsDone);

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
/** This friend function returns an MPO MPO0 that represents the sum \sum_{l} O_{l} of a local operator     O in a system of N spins.
    \param BC input: unsigned int, the boundary conditions
    \param N input: unsigned int, the number of spins
    \param O input: const Matrix<T>&, the local operator
    \param MPO0 output: MPO<T>&, the resulting MPO representing the sum of local operator O */
  friend void getMPOFromSumLocalOperator<>(string BC, unsigned int N, const Matrix<T>& O, MPO<T>& MPO0);

/// Computes expectation value of MPO with MPS.
/** This friend function computes the expectation value of the input MPO with the input MPS.
    \param MPS0 input: const MPS<T>&, the desired MPS
    \param MPO0 input: const MPO<T>&, the desired MPO
    \return T, the resulting expectation value */
  friend T expectationValueMPO<>(const MPS<T>& MPS0, const MPO<T>& MPO0);

/// Computes scalar product between MPS and MPO applied to another MPS.
/** This friend function computes the scalar product of the first input MPS with the input MPO applied
    to the second input MPS.
    \param MPS0 input: const MPS<T>&, the desired bra MPS
    \param MPO0 input: const MPO<T>&, the desired MPO
    \param MPS1 input: const MPS<T>&, the desired ket MPS
    \return T, the resulting scalar product */
  friend T scalarProductMPOMPS<>(const MPS<T>& MPS0, const MPO<T>& MPO0, const MPS<T>& MPS1);

/// Contracts MPS in reverse order with MPO applied to another MPS.
/** This friend function computes the scalar product of MPS0 with MPO0 applied to MPS1, taking the tensors
    of MPS0 in reverse order and not complex conjugating them.
    It is used in the PEPS contraction.
    \param MPS0 input: const MPS<T>&, the bra MPS that is taken in reverse order and
                       that is not complex conjugated
    \param MPO0 input: const MPO<T>&, the MPO
    \param MPS1 input: const MPS<T>&, the ket MPS
    \return T, the resulting contraction value
    \sa T MPS<T>::contractReverse(const MPS<T>& MPS0) const */
  friend T contractReverseMPOMPS<>(const MPS<T>& MPS0, const MPO<T>& MPO0, const MPS<T>& MPS1);

/// Multiplies MPO with MPS exactly.
/** This function computes the product of the input MPO0 with the input MPS0 exactly and stores
    the result in MPS1.
    \param MPO0 input: const MPO<T>&, the MPO
    \param MPS0 input: const MPS<T>&, the MPS on which the MPO acts
    \param MPS1 output: MPS<T>&, the resulting MPS */
  friend void exactMultiplyMPOMPS<>(const MPO<T>& MPO0, const MPS<T>& MPS0, MPS<T>& MPS1);

/// Computes distance between MPO-MPS product and another MPS.
/** This friend function computes the norm of the difference between the input MPO multiplied with the
    first input MPS and the second input MPS, i.e. || MPO0|MPS0> - |MPS1> ||.
    \param MPO0 input: const MPO<T>&, the MPO multiplied with MPS0
    \param MPS0 input: const MPS<T>&, the MPS multiplied by MPO0
    \param MPS1 input: const MPS<T>&, the approximating MPS
    \return double, the resulting distance */
  friend double distanceMPOMPS<>(const MPO<T>& MPO0, const MPS<T>& MPS0, const MPS<T>& MPS1);

/// Multiplies MPS with MPO and approximates result by MPS.
/** This function approximates the product of the input MPS with the input MPO by a MPS with a
    given bond dimension. It evaluates the error of the approximation after each sweep and stops either
    if the change in the error is below eps (in sweep n: eps(n) = abs(error(n) - error(n-1))) or if the
    number of sweeps exceeds maxNumSweeps. On input, if ((eps == 0.0) && (maxNumSweeps == 0) &&
    (MPS1.D == MPO0.D * MPS0.D)) then the exact multiplication is done. On output, errorAchieved is the
    achieved approximation error, numSweepsDone is the number of sweeps done and MPS1 is the
    approximating MPS.
    \param MPO0 input: const MPO<T>&, the MPO
    \param MPS0 input: const MPS<T>&, the MPS on which the MPO acts
    \param eps input: double, the desired precision (in sweep n: eps(n) = abs(error(n+1) - error(n-1)))
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param errorAchieved output: unsigned int&, the achieved approximation error
    \param numSweepsDone output: unsigned int&, the number of sweeps done
    \param MPS1 input/output: MPS<T>&, on input the initial approximating MPS, on output the resulting
                              MPS, must have the correct form */
  friend void multiplyMPOMPS<>(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                               unsigned int maxNumSweeps, double& errorAchieved,
                               unsigned int& numSweepsDone, MPS<T>& MPS1);

  friend void multiplyMPDOMPOMPS<>(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                   unsigned int maxNumSweeps, double& errorAchieved,
                                   unsigned int& numSweepsDone, MPS<T>& MPS1);

/// Multiplies MPS with MPO and approximates result by structured MPS.
/** This function approximates the product of the input MPS MPS0 with the input MPO MPO0 by a MPS MPS1
    with a given bond dimension and with structured tensors.
    If (eps >= 0), it evaluates the error of the approximation after each sweep and
    stops if the change in the error is below eps, i.e. in sweep n: eps >= abs(error(n)-error(n-1)). The
    function always stops if the number of sweeps exceeds maxNumSweeps. Hereby the error is defined as
       error := ||MPO0*|MPS0>-|MPS1>||/||MPO0*|MPS0>||   .
    If (eps < 0.0), the error is never computed and maxNumSweeps sweeps are performed.
    On input, if (MPS1.D >= MPO0.D*MPS0.D) then the exact multiplication is done.
    On output, errorAchieved is the achieved approximation error, numSweepsDone is the number of sweeps
    done and MPS1 is the approximating MPS.
    This function uses
       void updateTensor(const string& Direction, Tensor<T>& NTensor, Tensor<T>& bTensor,
                         double cutoff, unsigned int mode, unsigned int d2)
    and
       void updateTensor(unsigned int position, const string& Direction,
                         Tensor<T>& NTensorLeft, Tensor<T>& NTensorRight, Tensor<T>& bTensor,
                         double cutoff, unsigned int mode, unsigned int d2)
    to update the tensors.
    \param MPO0 input: const MPO<T>&, the MPO
    \param MPS0 input: const MPS<T>&, the MPS on which the MPO acts
    \param eps input: double, if (eps >= 0.0) the convergence precision, i.e. convergence in sweep n if
                      eps >= abs(error(n) - error(n-1)),
                      if (eps < 0.0) no error computation is done
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param errorAchieved output: double&, the achieved approximation error, if (eps < 0.0) not addressed
    \param numSweepsDone output: unsigned int&, the number of sweeps done
    \param MPS1 input/output: MPS<T>&, on input the initial approximating MPS, on output the resulting
                              MPS, must have the correct form
    \param cutoff input: double, the cutoff for the pseudoinverse
    \param mode input: unsigned int, the mode,
                       if (mode == 0) then the norm matrix is not processed before a LAPACK routine
                       solves the system of linear equations,
                       if (mode == 1) the system of linear equations is solved via computation
                       of the pseudoinverse of the norm matrix
    \param d2 input: unsigned int, the maximal d2,
                     if (d2 == 0) then the tensors of MPS1 are not processed,
                     if (d2 != 0) then the tensors of MPS1 are made hermitian with rank d2
    \param alphaStart input: double, used only if (mode == 2)
    \param x input: double, used only if (mode == 2)
    \param precision input: double, used only if (mode == 2)
    \param maxNumIter input: unsigned int, used only if (mode == 2)
    \sa void updateTensor(const string& Direction, Tensor<T>& NTensor, Tensor<T>& bTensor,
                          double cutoff, unsigned int mode, unsigned int d2,
                          double alphaStart, double x, double precision, unsigned int maxNumIter)
    \sa void updateTensor(unsigned int position, const string& Direction,
                          Tensor<T>& NTensorLeft, Tensor<T>& NTensorRight, Tensor<T>& bTensor,
                          double cutoff, unsigned int mode, unsigned int d2,
                          double alphaStart, double x, double precision, unsigned int maxNumIter) */
  friend void multiplyMPOMPS<>(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                               unsigned int maxNumSweeps, double& errorAchieved,
                               unsigned int& numSweepsDone, MPS<T>& MPS1,
                               double cutoff, unsigned int mode, unsigned int d2,
                               double alphaStart, double x, double precision, unsigned int maxNumIter);

/// Multiplies even-odd MPO with MPS by means of TEBD.
/** This function approximates the product of an even-odd MPO with a MPS by a MPS with the same bond
    dimension. An even-odd MPO arises in the even-odd Trotter decomposition of the evolution operator.
    \param MPO0 input: const MPO<T>&, the even-odd MPO
    \param MPS0 input: const MPS<T>&, the MPS on which the even-odd MPO acts
    \param MPS1 output: MPS<T>&, the resulting MPS, must be of the same form as MPS0 */
  friend void multiplyEvenOddMPOMPS<>(const MPO<T>& MPO0, const MPS<T>& MPS0, MPS<T>& MPS1);

 protected:

/// Boundary conditions of MPO.
 string BC;

/// Number of tensors of MPO.
 unsigned int N;

/// Physical dimensions of MPO.
 vector<unsigned int> d;

/// Maximal virtual bond dimension of MPO.
 unsigned int D;

/// Tensors of MPO.
 Tensor<T>* Tensors;
};

template<class T> MPO<T>::MPO()
{
 this->N = 0;
 this->D = 0;
 this->Tensors = 0;
}

template<class T> MPO<T>::MPO(string BC, unsigned int N, const vector<unsigned int>& d, unsigned int D)
{
#ifdef DEBUG
 if (((BC != "open") && (BC != "periodic")) || (N != d.size()))
 {
  cerr << "Program terminated because of error in constructor: " <<
          "template<class T> MPO<T>::" <<
          "MPO(string BC, unsigned int N, const vector<unsigned int>& d, unsigned int D): " <<
          "(((BC != open) && (BC != periodic)) || (N != d.size()))" << endl;
  exit(1);
 }
#endif
 this->BC = BC;
 this->N = N;
 this->d = d;
 this->D = D;
 this->Tensors = new Tensor<T>[this->N];
 vector<unsigned int> Shape(4);
 if (this->BC == "open")
 {
  Shape[0] = 1; Shape[1] = this->D; Shape[2] = this->d[0]; Shape[3] = this->d[0];
  this->Tensors[0] = Tensor<T>(Shape);
  Shape[0] = this->D; Shape[1] = this->D;
  for (int i = 1; i < this->N-1; i++)
  {
   Shape[2] = this->d[i]; Shape[3] = this->d[i];
   this->Tensors[i] = Tensor<T>(Shape);
  }
  Shape[0] = this->D; Shape[1] = 1; Shape[2] = this->d[this->N-1]; Shape[3] = this->d[this->N-1];
  this->Tensors[N-1] = Tensor<T>(Shape);
 }
 else if (this->BC == "periodic")
 {
  Shape[0] = this->D; Shape[1] = this->D;
  for (int i = 0; i < this->N; i++)
  {
   Shape[2] = this->d[i]; Shape[3] = this->d[i];
   this->Tensors[i] = Tensor<T>(Shape);
  }
 }
}

template<class T> MPO<T>::MPO(string BC, unsigned int N, unsigned int d, unsigned int D)
{
#ifdef DEBUG
 if ((BC != "open") && (BC != "periodic"))
 {
  cerr << "Program terminated because of error in constructor: " <<
          "template<class T> MPO<T>::" <<
          "MPO(string BC, unsigned int N, unsigned int d, unsigned int D): " <<
          "The BC have to be either open or periodic." << endl;
  exit(1);
 }
#endif
 this->BC = BC;
 this->N = N;
 this->d = vector<unsigned int>(this->N);
 for (int i = 0; i < this->N; i++)
 {
  this->d[i] = d;
 }
 this->D = D;
 this->Tensors = new Tensor<T>[this->N];
 vector<unsigned int> Shape(4);
 Shape[2] = d; Shape[3] = d;
 if (this->BC == "open")
 {
  Shape[0] = 1; Shape[1] = this->D;
  this->Tensors[0] = Tensor<T>(Shape);
  Shape[0] = this->D; Shape[1] = this->D;
  for (int i = 1; i < this->N-1; i++)
  {
   this->Tensors[i] = Tensor<T>(Shape);
  }
  Shape[0] = this->D; Shape[1] = 1;
  this->Tensors[N-1] = Tensor<T>(Shape);
 }
 else if (this->BC == "periodic")
 {
  Shape[0] = this->D; Shape[1] = this->D;
  for (int i = 0; i < this->N; i++)
  {
   this->Tensors[i] = Tensor<T>(Shape);
  }
 }
}

template<class T> MPO<T>::MPO(const MPO<T>& MPO0)
{
 this->BC = MPO0.BC;
 this->N = MPO0.N;
 this->d = MPO0.d;
 this->D = MPO0.D;
 this->Tensors = new Tensor<T>[this->N];
 for (int i = 0; i < this->N; i++)
 {
  this->Tensors[i] = MPO0.Tensors[i];
 }
}

template<class T> MPO<T>::~MPO()
{
 delete[] this->Tensors;
}

template<class T> MPO<T>& MPO<T>::operator=(const MPO<T>& MPO0)
{
 if (this != &MPO0)
 {
  this->BC = MPO0.BC;
  this->N = MPO0.N;
  this->d = MPO0.d;
  this->D = MPO0.D;
  delete[] this->Tensors;
  this->Tensors = new Tensor<T>[this->N];
  for (int i = 0; i < this->N; i++)
  {
   this->Tensors[i] = MPO0.Tensors[i];
  }
 }
 return *this;
}

template<class T> void MPO<T>::setd(unsigned int position, unsigned int d0)
{
#ifdef DEBUG
 if ((position > this->N-1) || (d0 == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPO<T>::" <<
          "setd(unsigned int position, unsigned int d0): " <<
          "((position > this->N-1) || (d0 == 0))." << endl;
  exit(1);
 }
#endif
 this->d[position] = d0;
}

template<class T> void MPO<T>::setD(unsigned int D0, const string& Filling)
{
#ifdef DEBUG
 if ((this->N == 0) || (D0 == 0) || ((Filling != "zeroes") && (Filling != "random")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPO<T>::" <<
          "setD(unsigned int D0, const string& Filling): " <<
          "((this->N == 0) || (D0 == 0) || ((Filling != zeroes) && (Filling != random)))." << endl;
  exit(1);
 }
#endif
 if (D0 == this->D)
  return;
 vector<unsigned int> Shape0(4), Shape1(4), Index(4);
 Tensor<T> Tensor0, Tensor1;
 for (unsigned int position = 0; position < this->N; position++)
 {
  this->get(position, Tensor0);
  Tensor0.getShape(Shape0);
  Shape1[0] = D0; Shape1[1] = D0; Shape1[2] = Shape0[2]; Shape1[3] = Shape0[3];
  if ((this->BC == "open") && (position == 0))
   Shape1[0] = 1;
  if ((this->BC == "open") && (position == this->N-1))
   Shape1[1] = 1;
  Tensor1 = Tensor<T>(Shape1);
  if (Filling == "zeroes")
   Tensor1.fillZeroes();
  else if (Filling == "random")
   Tensor1.fillRandomly();
  for (int i3 = 0; i3 < min(Shape0[3], Shape1[3]); i3++)
  {
   Index[3] = i3;
   for (int i2 = 0; i2 < min(Shape0[2], Shape1[2]); i2++)
   {
    Index[2] = i2;
    for (int i1 = 0; i1 < min(Shape0[1], Shape1[1]); i1++)
    {
     Index[1] = i1;
     for (int i0 = 0; i0 < min(Shape0[0], Shape1[0]); i0++)
     {
      Index[0] = i0;
      Tensor1.set(Index, Tensor0.get(Index));
     }
    }
   }
  }
  this->set(position, Tensor1);
 }
 this->D = D0;
}

template<class T> void MPO<T>::set(unsigned int position, const Tensor<T>& Tensor0)
{
#ifdef DEBUG
 if ((position > this->N-1) || (Tensor0.getRank() != 4))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPO<T>::" <<
          "set(unsigned int position, const Tensor<T>& Tensor0): " <<
          "((position > this->N-1) || (Tensor0.getRank() != 4))." << endl;
  exit(1);
 }
#endif
 this->Tensors[position] = Tensor0;
}

template<class T> void MPO<T>::get(unsigned int position, Tensor<T>& Tensor0) const
{
#ifdef DEBUG
 if (position > this->N-1)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPO<T>::" <<
          "get(unsigned int position, Tensor<T>& Tensor0) const: " <<
          "(position > this->N-1)." << endl;
  exit(1);
 }
#endif
 Tensor0 = this->Tensors[position];
}

template<class T> void MPO<T>::getMatrix(Matrix<T>& Matrix0) const
{
 unsigned int dim = 1;
 for (int i = 0; i < this->N; i++)
 {
  dim *= d[i];
 }
 vector<unsigned int> Shape;
 Matrix0.getShape(Shape);
#ifdef DEBUG
 if ((Shape[0] != dim) || (Shape[1] != dim))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void MPO<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const: " <<
          "Matrix0 has an incorrect shape." << endl;
  exit(1);
 }
#endif
 unsigned int dimRest, resti, restj;
 vector<unsigned int> Indices(2), Index(2), Indices0(1), Indices1(1), Index2(2);
 Indices[0] = 0; Indices[1] = 1;
 Tensor<T> Tensor0, Tensor1;
 Indices0[0] = 1; Indices1[0] = 0;
 Index2[0] = 0; Index2[1] = 0;
 for (int i = 0; i < dim; i++)
 {
  for (int j = 0; j < dim; j++)
  {
   dimRest = dim / d[0];
   resti = i; restj = j;
   Index[1] = resti / dimRest;
   Index[0] = restj / dimRest;
   this->Tensors[0].getSubtensor(Indices, Index, Tensor0);
   resti %= dimRest; restj %= dimRest;
   for (int k = 1; k < this->N; k++)
   {
    dimRest /= d[k];
    Index[1] = resti / dimRest;
    Index[0] = restj / dimRest;
    this->Tensors[k].getSubtensor(Indices, Index, Tensor1);
    Tensor0.contract(Indices0, Tensor1, Indices1);
    resti %= dimRest; restj %= dimRest;
   }
   if (this->BC == "open")
   {
    Matrix0(i, j) = Tensor0.get(Index2);
   }
   else if (this->BC == "periodic")
   {}
  }
 }
}

template<class T> void MPO<T>::fillRandomly(const vector<unsigned int>& Seed)
{
#ifdef DEBUG
 if ((this->N == 0) || (Seed.size() != this->N))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPO<T>::" <<
          "fillRandomly(const vector<unsigned int>& Seed): " <<
          "((this->N == 0) || (Seed.size() != this->N))." << endl;
  exit(1);
 }
#endif
 for (int i = 0; i < this->N; i++)
 {
  this->Tensors[i].fillRandomly(Seed[i]);
 }
}

template<class T> void MPO<T>::fillZeroes()
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPO<T>::" <<
          "fillZeroes(): " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
 for (int i = 0; i < this->N; i++)
 {
  this->Tensors[i].fillZeroes();
 }
}

template<class T> T MPO<T>::simplifiedNormalize()
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T MPO<T>::" <<
          "simplifiedNormalize(): " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
 T result = 1.0;
 for (int position = 0; position < this->N; position++)
 {
  result *= this->Tensors[position].normalize();
 }
 return result;
}

template<class T> void MPO<T>::transpose()
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPO<T>::" <<
          "transpose(): " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
 vector<unsigned int> Order(4);
 Order[0] = 0; Order[1] = 1; Order[2] = 3; Order[3] = 2;
 for (unsigned int position = 0; position < this->N; position++)
  this->Tensors[position].permute(Order);
}

template<class T> void MPO<T>::adjoint(MPO<T>& MPO0) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPO<T>::" <<
          "adjoint(MPO<T>& MPO0): " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
 MPO0 = *this;
 vector<unsigned int> Order(4);
 Order[0] = 0; Order[1] = 1; Order[2] = 3; Order[3] = 2;
 for (unsigned int position = 0; position < this->N; position++)
 {
  MPO0.Tensors[position].permute(Order);
  MPO0.Tensors[position].complexConjugate();
 }
}

template<class T> void MPO<T>::getMPS(MPS<T>& MPS0) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPO<T>::" <<
          "getMPS(MPS<T>& MPS0) const: " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
 vector<unsigned int> d2(this->N), Shape3(3), Shape4(4);
 Tensor<T> Tensor0;
 for (int pos = 0; pos < this->N; pos++)
 {
  this->Tensors[pos].getShape(Shape4);
  d2[pos] = Shape4[2]*Shape4[3];
 }
 MPS0 = MPS<T>(this->BC, this->N, d2, this->D);
 for (int pos = 0; pos < this->N; pos++)
 {
  Tensor0 = this->Tensors[pos];
  Tensor0.getShape(Shape4);
  Shape3[0] = Shape4[0]; Shape3[1] = Shape4[1]; Shape3[2] = Shape4[2]*Shape4[3];
  Tensor0.reshape(Shape3);
  MPS0.set(pos, Tensor0);
 }
 if (this->BC == "open")
  MPS0.bringIntoNormalShape();
}

template<class T> void MPO<T>::canonicalize(double eps, MPO<T>& MPO0, vector< Matrix<T> >& Lambdas) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPO<T>::" <<
          "canonicalize(double eps, MPO<T>& MPO0, vector< Matrix<T> >& Lambdas) const: " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
 vector<unsigned> Shape, Shape0;
 Tensor<T> Tensor0;
 MPS<T> MPS0, MPS1;
 MPO0 = *this;
 MPO0.getMPS(MPS0);
 MPS0.canonicalize(eps, MPS1, Lambdas);
 for (unsigned int position = 0; position < MPO0.N; position++)
 {
  MPS1.get(position, Tensor0);
  MPO0.Tensors[position].getShape(Shape);
  Tensor0.getShape(Shape0);
  Shape[0] = Shape0[0]; Shape[1] = Shape0[1];
  Tensor0.reshape(Shape);
  MPO0.Tensors[position] = Tensor0;
 }
}

template<class T> void MPO<T>::decanonicalize(const vector< Matrix<T> >& Lambdas, MPO<T>& MPO0) const
{
#ifdef DEBUG
 if ((this->N == 0) || (this->BC == "periodic") || (Lambdas.size() != this->N-1))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPO<T>::" <<
          "decanonicalize(const vector< Matrix<T> >& Lambdas, MPO<T>& MPO0) const: " <<
          "((this->N == 0) || (this->BC == periodic) || (Lambdas.size() != this->N-1))." << endl;
  exit(1);
 }
#endif
 unsigned int position;
 vector<unsigned int> Indices0(1), Indices1(1);
 Tensor<T> Tensor0, Tensor1;
 MPO0 = *this;
 Indices0[0] = 1; Indices1[0] = 0;
 for (unsigned int position = 1; position < MPO0.N; position++)
 {
  Tensor0 = Lambdas[position-1]; Tensor1 = MPO0.Tensors[position];
  Tensor0.contract(Indices0, Tensor1, Indices1);
  MPO0.Tensors[position] = Tensor0;
 }
}

template<class T> void MPO<T>::multiply(const MPO<T>& MPO0, MPO<T>& MPO1) const
{
#ifdef DEBUG
 if ((this->BC != MPO0.BC) || (this->N != MPO0.N) || (this->N == 0) || (this->d != MPO0.d))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPO<T>::" <<
          "multiply(const MPO<T>& MPO0, MPO<T>& MPO1) const: " <<
          "This MPO and MPO0 do not have the same form or (this->N == 0)." << endl;
  exit(1);
 }
#endif
 MPO1 = MPO<T>(this->BC, this->N, this->d, this->D*MPO0.D);
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Indices0(1), Indices1(1), Order(6), Shape0, Shape1;
 Indices0[0] = 2; Indices1[0] = 3;
 Order[0] = 0; Order[1] = 3; Order[2] = 1; Order[3] = 4; Order[4] = 5; Order[5] = 2;
 for (unsigned int position = 0; position < this->N; position++)
 {
  Tensor0 = this->Tensors[position]; Tensor1 = MPO0.Tensors[position];
  Tensor0.getShape(Shape0); Tensor1.getShape(Shape1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order);
  Shape0[0] = Shape0[0]*Shape1[0]; Shape0[1] = Shape0[1]*Shape1[1]; Shape0[2] = Shape1[2];
  Tensor0.reshape(Shape0);
  MPO1.Tensors[position] = Tensor0;
 }
}

template<class T> double MPO<T>::distance(const MPO<T>& MPO0, const MPO<T>& MPO1) const
{
#ifdef DEBUG
 if ((this->BC != MPO0.BC) || (this->BC != MPO1.BC) || (this->N != MPO0.N) || (this->N != MPO1.N) ||
     (this->N == 0) || (this->d != MPO0.d) || (this->d != MPO1.d))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> double MPO<T>::" <<
          "distance(const MPO<T>& MPO0, const MPO<T>& MPO1) const: " <<
          "This MPO, MPO0, and MPO1 are not of the same form or (this->N == 0)." << endl;
  exit(1);
 }
#endif
// normSquaredMPOMPO0 := tr((MPO*MPO0)^{+}*(MPO*MPO0)):
 MPO<T> MPO2; MPS<T> MPS0;
 this->multiply(MPO0, MPO2);
 MPO2.getMPS(MPS0);
 double normSquaredMPOMPO0 = MPS0.scalarProduct();
// realScalarProductMPO1MPOMPO0 = real(tr(MPO1^{+}*(MPO*MPO0))):
 MPS<T> MPS1;
 MPO1.getMPS(MPS1);
 double realScalarProductMPO1MPOMPO0 = MathAuxiliary::convertToDouble(MPS1.scalarProduct(MPS0));
// normSquaredMPO1 = tr(MPO1^{+}*MPO1):
 double normSquaredMPO1 = MPS1.scalarProduct();
// returning the resulting distance:
 return sqrt(abs(normSquaredMPOMPO0 - 2.0*realScalarProductMPO1MPOMPO0 + normSquaredMPO1));
}

template<class T> void MPO<T>::multiply(const MPO<T>& MPO0, double eps, unsigned int maxNumSweeps,
                                        double& errorAchieved, unsigned int& numSweepsDone, MPO<T>& MPO1)
                                        const
{
#ifdef DEBUG
 if ((this->BC != MPO0.BC) || (this->BC != MPO1.BC) || (this->N != MPO0.N) || (this->N != MPO1.N) ||
     (this->N == 0) || (this->d != MPO0.d) || (this->d != MPO1.d))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPO<T>::" <<
          "multiply(const MPO<T>& MPO0, double eps, unsigned int maxNumSweeps, double& errorAchieved, " <<
                   "unsigned int& numSweepsDone, MPO<T>& MPO1): " <<
          "This MPO, MPO0, and MPO1 are not of the same form or (this->N == 0)." << endl;
  exit(1);
 }
#endif
 unsigned int position; T element;
 vector<unsigned int> d2(this->N), Shape3(3), Shape04(4), Shape4(4), Shape6(6), Index4(4), Index6(6);
 Tensor<T> Tensor0, Tensor1;
// d2 := this->d*this->d:
 for (int i = 0; i < this->N; i++)
  d2[i] = this->d[i]*this->d[i];
// define MPO2, MPS0, and MPS1 for multiplyMPOMPS:
// - MPO2 is the same as this MPO where each tensor gets two additional physical indices acting as a delta function
// - MPS0 results from reshaping the tensors in MPO0
// - MPS1 results from reshaping the tensors in MPO1
 MPO<T> MPO2(this->BC, this->N, d2, this->D);
 MPS<T> MPS0(MPO0.BC, MPO0.N, d2, MPO0.D), MPS1(MPO1.BC, MPO1.N, d2, MPO1.D);
 for (position = 0; position < this->N; position++)
 {
// - MPS0:
  Tensor0 = MPO0.Tensors[position];
  Tensor0.getShape(Shape04);
  Shape3[0] = Shape04[0]; Shape3[1] = Shape04[1]; Shape3[2] = Shape04[2]*Shape04[3];
  Tensor0.reshape(Shape3);
  MPS0.set(position, Tensor0);
// - MPO2:
  Tensor0 = this->Tensors[position];
  Tensor0.getShape(Shape4);
  Shape6[0] = Shape4[0]; Shape6[1] = Shape4[1]; Shape6[2] = Shape04[2]; Shape6[3] = Shape4[2];
  Shape6[4] = Shape04[2]; Shape6[5] = Shape4[3];
  Tensor1 = Tensor<T>(Shape6);
  Tensor1.fillZeroes();
  for (int i5 = 0; i5 < Shape6[5]; i5++)
  {
   Index4[3] = i5;
   Index6[5] = i5;
   for (int i4 = 0; i4 < Shape6[4]; i4++)
   {
    Index6[4] = i4; Index6[2] = i4;
    for (int i3 = 0; i3 < Shape6[3]; i3++)
    {
     Index4[2] = i3;
     Index6[3] = i3;
     for (int i1 = 0; i1 < Shape6[1]; i1++)
     {
      Index4[1] = i1;
      Index6[1] = i1;
      for (int i0 = 0; i0 < Shape6[0]; i0++)
      {
       Index4[0] = i0;
       Index6[0] = i0;
       element = Tensor0.get(Index4);
       Tensor1.set(Index6, element);
      }
     }
    }
   }
  }
  Shape4[2] = Shape04[2]*Shape4[2]; Shape4[3] = Shape04[2]*Shape4[3];
  Tensor1.reshape(Shape4);
  MPO2.set(position, Tensor1);
// - MPS1:
  Tensor0 = MPO1.Tensors[position];
  Tensor0.getShape(Shape4);
  Shape3[0] = Shape4[0]; Shape3[1] = Shape4[1]; Shape3[2] = Shape4[2]*Shape4[3];
  Tensor0.reshape(Shape3);
  MPS1.set(position, Tensor0);
 }
 MPS0.bringIntoNormalShape(); MPS1.bringIntoNormalShape();
 multiplyMPOMPS(MPO2, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
// reshape tensors in MPS1 for MPO1:
 for (position = 0; position < this->N; position++)
 {
  MPS1.get(position, Tensor0);
  MPO1.Tensors[position].getShape(Shape4);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1];
  Tensor0.reshape(Shape4);
  MPO1.set(position, Tensor0);
 }
}

template<class T> void MPO<T>::setd2(unsigned int d2, double eps, vector< vector<double> >& Eigenvalues)
{
#ifdef DEBUG
 if ((this->BC != "open") || (this->N == 0) || (d2 == 0) || (eps < 0.0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPO<T>::" <<
          "setd2(unsigned int d2, double eps, vector< vector<double> >& Eigenvalues): " <<
          "((this->BC != open) || (this->N == 0) || (d2 == 0) || (eps < 0.0))." << endl;
  exit(1);
 }
#endif
 Eigenvalues = vector< vector<double> >(this->N);
 unsigned int position, count = 0;
 for (position = 0; position < this->N; position++)
 {
  if (this->d[position] <= d2)
   count++;
 }
 if (count == this->N)
  return;
 unsigned int d0;
 vector<unsigned int> Indices01(1), Indices11(1), Indices03(3), Indices13(3);
 vector<double> W, W0;
 vector< Matrix<T> > Lambdas;
 Tensor<T> Tensor0, Tensor1;
 Matrix<T> Matrix0, Vr, P;
 MPO<T> MPO0;
 this->canonicalize(eps, MPO0, Lambdas);
// position 0:
 position = 0;
 d0 = this->d[position];
 if (d0 > d2)
 {
  MPO0.get(position, Tensor0); Tensor1 = Lambdas[position];
  Indices01[0] = 1; Indices11[0] = 0;
  Tensor0.contract(Indices01, Tensor1, Indices11);
  Tensor0.complexConjugate(Tensor1);
  Indices03[0] = 0; Indices03[1] = 1; Indices03[2] = 3;
  Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 3;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  Matrix0 = Matrix<T>(d0, d0);
  for (int i = 0; i < Tensor0.getSize(); i++)
   Matrix0.set(i, Tensor0.get(i));
  Matrix0.setType("hermitian");
  W = vector<double>(d0); Vr = Matrix<T>(d0, d0);
  Matrix0.eigenDecompose(W, Vr);
  W0 = W;
  for (int i = 0; i < d0; i++)
   W0[i] = W[d0-1-i];
  Eigenvalues[position] = W0;
  P = Matrix<T>(d2, d0);
  for (int i = 0; i < d2; i++)
  {
   for (int j = 0; j < d0; j++)
   {
    P(i, j) = MathAuxiliary::complexConjugate(Vr(j, d0-1-i));
   }
  }
  MPO0.get(position, Tensor0);
  Indices01[0] = 3; Indices11[0] = 1;
  Tensor0.contract(Indices01, P, Indices11);
  MPO0.d[position] = d2;
  MPO0.set(position, Tensor0);
 }
// 0 < position < this->N-1:
 Indices03[0] = 0; Indices03[1] = 2; Indices03[2] = 3;
 Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
 for (position = 1; position < this->N-1; position++)
 {
  d0 = this->d[position];
  if (d0 > d2)
  {
   MPO0.get(position, Tensor0); Tensor1 = Lambdas[position-1];
   Indices01[0] = 0; Indices11[0] = 1;
   Tensor0.contract(Indices01, Tensor1, Indices11);
   Tensor1 = Lambdas[position];
   Indices11[0] = 0;
   Tensor0.contract(Indices01, Tensor1, Indices11);
   Tensor0.complexConjugate(Tensor1);
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Matrix0 = Matrix<T>(d0, d0);
   for (int i = 0; i < Tensor0.getSize(); i++)
    Matrix0.set(i, Tensor0.get(i));
   Matrix0.setType("hermitian");
   W = vector<double>(d0); Vr = Matrix<T>(d0, d0);
   Matrix0.eigenDecompose(W, Vr);
   W0 = W;
   for (int i = 0; i < d0; i++)
    W0[i] = W[d0-1-i];
   Eigenvalues[position] = W0;
   P = Matrix<T>(d2, d0);
   for (int i = 0; i < d2; i++)
   {
    for (int j = 0; j < d0; j++)
    {
     P(i, j) = MathAuxiliary::complexConjugate(Vr(j, d0-1-i));
    }
   }
   MPO0.get(position, Tensor0);
   Indices01[0] = 3; Indices11[0] = 1;
   Tensor0.contract(Indices01, P, Indices11);
   MPO0.d[position] = d2;
   MPO0.set(position, Tensor0);
  }
 }
// position this->N-1:
 position = this->N-1;
 d0 = this->d[position];
 if (d0 > d2)
 {
  MPO0.get(position, Tensor0); Tensor1 = Lambdas[position-1];
  Indices01[0] = 0; Indices11[0] = 1;
  Tensor0.contract(Indices01, Tensor1, Indices11);
  Tensor0.complexConjugate(Tensor1);
  Indices03[0] = 0; Indices03[1] = 1; Indices03[2] = 3;
  Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 3;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  Matrix0 = Matrix<T>(d0, d0);
  for (int i = 0; i < Tensor0.getSize(); i++)
   Matrix0.set(i, Tensor0.get(i));
  Matrix0.setType("hermitian");
  W = vector<double>(d0); Vr = Matrix<T>(d0, d0);
  Matrix0.eigenDecompose(W, Vr);
  W0 = W;
  for (int i = 0; i < d0; i++)
   W0[i] = W[d0-1-i];
  Eigenvalues[position] = W0;
  P = Matrix<T>(d2, d0);
  for (int i = 0; i < d2; i++)
  {
   for (int j = 0; j < d0; j++)
   {
    P(i, j) = MathAuxiliary::complexConjugate(Vr(j, d0-1-i));
   }
  }
  MPO0.get(position, Tensor0);
  Indices01[0] = 3; Indices11[0] = 1;
  Tensor0.contract(Indices01, P, Indices11);
  MPO0.d[position] = d2;
  MPO0.set(position, Tensor0);
 }
 MPO0.decanonicalize(Lambdas, *this);
}

template<class T> void MPO<T>::setd2(unsigned int d2, double cutoff, double eps, unsigned int maxNumSweeps,
                                     double& errorAchieved, unsigned int& numSweepsDone)
{
#ifdef DEBUG
 if ((this->N == 0) || (d2 == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPO<T>::" <<
          "setd2(unsigned int d2, double cutoff, double eps, unsigned int maxNumSweeps, " <<
                "double& errorAchieved, unsigned int& numSweepsDone): " <<
          "((this->N == 0) || (d2 == 0))." << endl;
  exit(1);
 }
#endif
 vector<unsigned int> Shape0(4), Shape1(4), Index(4);
 Tensor<T> Tensor0, Tensor1;
// if d2 is larger or equal to the current purification bond:
 unsigned int count = 0;
 for (int pos = 0; pos < this->N; pos++)
 {
  this->get(pos, Tensor0);
  Tensor0.getShape(Shape0);
  if (Shape0[3] < d2)
  {
   count++;
   Shape1 = Shape0; Shape1[3] = d2;
   Tensor1 = Tensor<T>(Shape1);
   Tensor1.fillZeroes();
   for (int i3 = 0; i3 < Shape0[3]; i3++){
    Index[3] = i3;
   for (int i2 = 0; i2 < Shape0[2]; i2++){
    Index[2] = i2;
   for (int i1 = 0; i1 < Shape0[1]; i1++){
    Index[1] = i1;
   for (int i0 = 0; i0 < Shape0[0]; i0++){
    Index[0] = i0;
    Tensor1.set(Index, Tensor0.get(Index));
   }
   }
   }
   }
   this->set(pos, Tensor1);
  }
  else if (Shape0[3] == d2)
  {
   count++;
  }
 }
// if d2 is smaller than the current purification bond:
 if (count < this->N)
 {
// MPS representation MPS0 of this MPDO and MPS1 as initial approximation:
  MPS<T> MPS0, MPS1;
  MPO<T> MPO0(*this), MPO1, MPO2;
  MPO0.adjoint(MPO1);
  MPO1.multiply(MPO0, MPO2);
  MPO2.getMPS(MPS0);
  for (int pos = 0; pos < this->N; pos++)
  {
   MPO0.get(pos, Tensor0);
   Tensor0.getShape(Shape0);
   Shape1 = Shape0; Shape1[3] = d2;
   Tensor1 = Tensor<T>(Shape1);
   for (int i3 = 0; i3 < Shape1[3]; i3++){
    Index[3] = i3;
   for (int i2 = 0; i2 < Shape1[2]; i2++){
    Index[2] = i2;
   for (int i1 = 0; i1 < Shape1[1]; i1++){
    Index[1] = i1;
   for (int i0 = 0; i0 < Shape1[0]; i0++){
    Index[0] = i0;
    Tensor1.set(Index, Tensor0.get(Index));
   }
   }
   }
   }
   MPO0.set(pos, Tensor1);
  }
  MPO0.adjoint(MPO1);
  MPO1.multiply(MPO0, MPO2);
  MPO2.getMPS(MPS1);
// build IMPO:=1:
  vector<unsigned int> d0;
  MPS0.getd(d0);
  MPO<T> IMPO(this->BC, this->N, d0, 1);
  Shape0[0] = 1; Shape0[1] = 1;
  Index[0] = 0; Index[1] = 0;
  for (int pos = 0; pos < this->N; pos++)
  {
   Shape0[2] = d0[pos]; Shape0[3] = d0[pos];
   Tensor0 = Tensor<T>(Shape0);
   Tensor0.fillZeroes();
   for (int i = 0; i < d0[pos]; i++){
    Index[2] = i; Index[3] = i;
    Tensor0.set(Index, 1.0);
   }
   IMPO.set(pos, Tensor0);
  }
// multiplyMPOMPS:
  IMPO.setD(2);
  unsigned int mode = 2;
  multiplyMPOMPS(IMPO, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1, cutoff, mode, d2);
// set new tensors in this MPDO:
  unsigned int dim, d3;
  vector<unsigned int> Shape2(2), Shape3(3), Shape6(6), Order6(6), Index2(2);
  vector<double> W;
  Matrix<T> Matrix1, Vr;
  Order6[0] = 0; Order6[1] = 2; Order6[2] = 4; Order6[3] = 1; Order6[4] = 3; Order6[5] = 5;
  for (int pos = 0; pos < this->N; pos++)
  {
   MPS1.get(pos, Tensor1);
   Tensor1.getShape(Shape3);
   dim = sqrt(double(Shape3[0]*Shape3[1]*Shape3[2]));
   Shape6[0] = sqrt(double(Shape3[0])); Shape6[1] = Shape6[0]; Shape6[2] = sqrt(double(Shape3[1]));
   Shape6[3] = Shape6[2]; Shape6[4] = sqrt(double(Shape3[2])); Shape6[5] = Shape6[4];
   Tensor1.reshape(Shape6);
   Tensor1.permute(Order6);
   Shape2[0] = dim; Shape2[1] = dim;
   Tensor1.reshape(Shape2);
   Matrix1 = Matrix<T>(dim, dim);
   for (int j = 0; j < dim; j++){
    Index2[1] = j;
    for (int i = 0; i < dim; i++){
     Index2[0] = i;
     Matrix1(i, j) = Tensor1.get(Index2);
    }
   }
   Matrix1.setType("hermitian");
   W = vector<double>(dim); Vr = Matrix<T>(dim, dim);
   Matrix1.eigenDecompose(W, Vr);
   d3 = min(d2, dim);
   Shape0[0] = sqrt(double(Shape3[0])); Shape0[1] = sqrt(double(Shape3[1])); Shape0[2] = sqrt(double(Shape3[2]));
   Shape0[3] = d3;
   Tensor0 = Tensor<T>(Shape0);
   for (int j = dim-1; j >= int(dim-d3); j--){
    for (int i = 0; i < dim; i++){
     Tensor0.set(i+(dim-1-j)*dim, Vr(i, j)*sqrt(W[j]));
    }
   }
   this->set(pos, Tensor0);
  }
 }
}

template<class T> void getMPOFromLocalOperator(string BC, unsigned int N, unsigned int x,
                                               const Matrix<T>& O, MPO<T>& MPO0)
{
#ifdef DEBUG
 if (((BC != "open") && (BC != "periodic")) || (N == 0) || (x > N-1) || (O.size == 0) ||
     (O.getDim0() != O.getDim1()))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> void " <<
          "getMPOFromLocalOperator(string BC, unsigned int N, unsigned int x, " <<
                                  "const Matrix<T>& O, MPO<T>& MPO0): " <<
          "(((BC != open) && (BC != periodic)) || (N == 0) || (x > N-1) || (O.size == 0) || " <<
           "(O.getDim0() != O.getDim1()))." << endl;
  exit(1);
 }
#endif
 unsigned int d = O.getDim0(), D = 1;
 MPO0 = MPO<T>(BC, N, d, D);
 vector<unsigned int> Shape(4);
 Shape[0] = 1; Shape[1] = 1; Shape[2] = d; Shape[3] = d;
 Tensor<T> A0(Shape), A1(Shape);
 T element;
 vector<unsigned int> Index(4);
// the Identity tensor A0=Identity has to be set on every position but on position x:
 A0.fillZeroes();
 element = 1.0;
 Index[0] = 0; Index[1] = 0;
 for (int i = 0; i < d; i++)
 {
  Index[2] = i; Index[3] = i;
  A0.set(Index, element);
 }
 unsigned int position;
 for (position = 0; position < N; position++)
 {
  if (position != x)
  {
   MPO0.set(position, A0);
  }
 }
// the Operator tensor A1=O has to be set on position x:
 A1.fillZeroes();
 Index[0] = 0; Index[1] = 0;
 vector<unsigned int> Index2(2);
 for (int j = 0; j < d; j++)
 {
  for (int i = 0; i < d; i++)
  {
   Index2[0] = i; Index2[1] = j;
   element = O.get(Index2);
   Index[2] = j; Index[3] = i;
   A1.set(Index, element);
  }
 }
 position = x;
 MPO0.set(position, A1);
}

template<class T> void getMPOFromSumLocalOperator(string BC, unsigned int N, const Matrix<T>& O,
                                                  MPO<T>& MPO0)
{
#ifdef DEBUG
 if (((BC != "open") && (BC != "periodic")) || (N == 0) || (O.size == 0) ||
     (O.getDim0() != O.getDim1()))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> void " <<
          "getMPOFromSumLocalOperator(string BC, unsigned int N, const Matrix<T>& O, MPO<T>& MPO0): " <<
          "(((BC != open) && (BC != periodic)) || (N == 0) || (O.size == 0) || " <<
           "(O.getDim0() != O.getDim1()))." << endl;
  exit(1);
 }
#endif
 unsigned int d = O.getDim0(), D = 2;
 MPO0 = MPO<T>(BC, N, d, D);
 if (BC == "open")
 {
// the MPO operators Sigmas^{i}:
  vector<unsigned int> Shape(3);
  Shape[0] = 2; Shape[1] = d; Shape[2] = d;
  Tensor<T> Sigmas(Shape);
  Sigmas.fillZeroes();
  vector<unsigned int> Index(3);
// - Sigmas^{0} = Identity:
  T element = 1.0;
  Index[0] = 0;
  for (int i = 0; i < d; i++)
  {
   Index[1] = i; Index[2] = i;
   Sigmas.set(Index, element);
  }
// - Sigmas^{1} = O:
  Index[0] = 1;
  vector<unsigned int> Index2(2);
  for (int j = 0; j < d; j++)
  {
   for (int i = 0; i < d; i++)
   {
    Index2[0] = i; Index2[1] = j;
    element = O.get(Index2);
    Index[1] = j; Index[2] = i;
    Sigmas.set(Index, element);
   }
  }
// A for position 0:
  Shape[0] = 1; Shape[1] = D; Shape[2] = 2;
  Tensor<T> A(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 0; Index[1] = 1; Index[2] = 1;
  A.set(Index, element);
  Tensor<T> SigmasC(Sigmas);
  vector<unsigned int> IndexA(1), IndexSigmasC(1);
  IndexA[0] = 2; IndexSigmasC[0] = 0;
  A.contract(IndexA, SigmasC, IndexSigmasC);
  unsigned int position = 0;
  MPO0.set(position, A);
// A for position 1 <= l <= N-2:
  Shape[0] = D; Shape[1] = D; Shape[2] = 2;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 1; Index[1] = 1; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 0; Index[1] = 1; Index[2] = 1;
  A.set(Index, element);
  SigmasC = Sigmas;
  A.contract(IndexA, SigmasC, IndexSigmasC);
  for (position = 1; position < N-1; position++)
  {
   MPO0.set(position, A);
  }
// A for position N-1:
  Shape[0] = D; Shape[1] = 1; Shape[2] = 2;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 1; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 0; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  SigmasC = Sigmas;
  A.contract(IndexA, SigmasC, IndexSigmasC);
  position = N-1;
  MPO0.set(position, A);
 }
 else if (BC == "periodic")
 {
  cerr << "The following friend function is not implemented yet for periodic BC: " <<
          "template<class T> void " <<
          "getMPOFromSumLocalOperator(string BC, unsigned int N, const Matrix<T>& O, " <<
                                     "MPO<T>& MPO0)." << endl;
  exit(1);
 }
}
