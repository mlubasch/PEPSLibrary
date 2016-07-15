/// Template class MPS implements Matrix Product States.
/** The template class MPS implements Matrix Product States as Projected Entangled-Pair States in 1D.
    MPS consist of Tensors on a linear chain.
    If the MPS has "open" boundary conditions, it is always in normal shape meaning that the Shape of
    each Tensor is given by the normal form. If the MPS has "periodic" boundary conditions, all Tensors
    have the same Shape with virtual dimension D.
    \param BC string, the boundary conditions of the MPS, is "open" or "periodic"
    \param N unsigned int, the number of tensors of the MPS
    \param d vector<unsigned int>, the physical dimensions of the MPS
    \param D unsigned int, the maximal virtual bond dimension of the MPS
    \param Tensors Tensor<T>*, the tensors of the MPS
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

// Declaration of template classes:
template<class T> class Tensor;
template<class T> class MPS;
template<class T> class MPO;

// Declaration of friend functions:
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
template<class T> void getMPS(Tensor<T>& Tensor0, MPS<T>& MPS0);
template<class T> void getTensor(const MPS<T>& MPS0, Tensor<T>& Tensor0);

template<class T> class MPS
{
 public:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0.
    \sa MPS(string BC, unsigned int N, const vector<unsigned int>& d, unsigned int D)
    \sa MPS(string BC, unsigned int N, unsigned int d, unsigned int D) */
  MPS();

/// Constructor for MPS with specific BC, N, d and D.
/** This constructor initializes a MPS of a specific form.
    \param BC input: string, the boundary conditions, must be "open" or "periodic"
    \param N input: unsigned int, the number of tensors, must be > 1
    \param d input: const vector<unsigned int>&, the physical dimensions, must fulfill d.size()==N
    \param D input: unsigned int, the maximal virtual bond dimension, must be > 0
    \sa MPS()
    \sa MPS(string BC, unsigned int N, unsigned int d, unsigned int D) */
  MPS(string BC, unsigned int N, const vector<unsigned int>& d, unsigned int D);

/// Constructor for MPS with specific BC, N, d and D.
/** This constructor initializes a MPS of a specific form with all physical dimensions equal to d.
    \param BC input: string, the boundary conditions, must be "open" or "periodic"
    \param N input: unsigned int, the number of tensors, must be > 1
    \param d input: unsigned int, the physical dimension, must be > 1
    \param D input: unsigned int, the maximal virtual bond dimension, must be > 0
    \sa MPS()
    \sa MPS(string BC, unsigned int N, const vector<unsigned int>& d, unsigned int D) */
  MPS(string BC, unsigned int N, unsigned int d, unsigned int D);

/// Standard copy constructor.
/** The standard copy constructor copies the input MPS into this.
    \param MPS0 input: const MPS<T>&, to be copied into this
    \sa MPS<T>& operator=(const MPS<T>& MPS0) */
  MPS(const MPS<T>& MPS0);

/// Standard destructor.
/** The standard destructor deletes the elements of the MPS. */
  ~MPS();

/// Assigns MPS to this.
/** The operator= allows to assign a MPS to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side MPS.
    \param MPS0 input: const MPS<T>&, to be copied into this
    \return MPS<T>&, a reference to the new this
    \sa MPS(const MPS<T>& MPS0) */
  MPS<T>& operator=(const MPS<T>& MPS0);

/// Returns boundary conditions of MPS.
/** The returned boundary conditions are either "open" or "periodic".
    \param BC0 output: string&, the boundary conditions of this MPS */
  void getBC(string& BC0) const { BC0 = this->BC; }

/// Returns number of tensors of MPS.
/** This function returns the number of tensors building up this MPS.
    \return unsigned int, the number of tensors of this MPS */
  unsigned int getN() const { return this->N; }

/// Returns physical dimensions of MPS.
/** This function returns the physical dimensions of this MPS.
    \param d0 output: vector<unsigned int>&, the physical dimensions of this MPS */
  void getd(vector<unsigned int>& d0) const { d0 = this->d; }

/// Returns physical dimension d.
/** This function returns the physical dimension of this MPS at site 0. It is useful, if this
    MPS has all physical dimensions equal to each other.
    \return unsigned int, the physical dimension of this MPS at site 0 */
  unsigned int getd() const { return this->d[0]; }

/// Sets maximal virtual bond dimension.
/** This function changes the maximal virtual bond dimension of this MPS.
    If the maximal virtual bond dimension is decreased, D0 < this->D, then the dispensable elements are
    simply discarded. If it is increased, D0 > this->D, then the new elements are set as random numbers
    multiplied by element.
    Element has the default value 0.0, in which case the new elements are zeroes.
    The random numbers filling a tensor at position pos are seeded with time(0)+pos*13.
    \param D0 input: unsigned int, the new maximal virtual bond dimension, must be > 0
    \param element optional input: T, multiplies all random numbers */
  void setD(unsigned int D0, T element = 0.0);

/// Returns maximal virtual bond dimension of MPS.
/** This function returns the maximal virtual bond dimension of the MPS.
    \return unsigned int, the maximal virtual bond dimension of the MPS */
  unsigned int getD() const { return this->D; }

/// Sets Tensor at position.
/** Given a position the corresponding tensor in Tensors is set.
    Since every MPS with open boundary conditions is assumed to be in normal shape, we should always use
       void bringIntoNormalShape()
    after setting new Tensors in this MPS.
    \param position input: unsigned int, the position
    \param Tensor0 input: const Tensor<T>&, to be written at position */
  void set(unsigned int position, const Tensor<T>& Tensor0);

/// Returns Tensor at position.
/** Given a position the corresponding tensor in Tensors is returned.
    \param position input: unsigned int, the position
    \param Tensor0 output: Tensor<T>&, a copy of the tensor at position */
  void get(unsigned int position, Tensor<T>& Tensor0) const;

/// Writes this MPS to binary file.
/** Given a file name FileName, a new binary file is constructed into which this MPS is written.
    A MPS is represented in a binary file by:
    {BCSize, BC[0], ..., BC[BCSize-1], N, d[0], ..., d[N-1], D, Tensors[0], ..., Tensors[N-1]}
    where each Tensor is represented by:
    {rank, Shape[0], ..., Shape[rank-1], size, Elements[0], ..., Elements[size-1]}.
    \param FileName input: const string&, the name for the new binary file to which this MPS is
                           written
    \sa void read(const string& FileName) */
  void write(const string& FileName) const;

/// Reads MPS from binary file.
/** Given a binary file called FileName, this MPS is replaced by the MPS in FileName.
    A MPS is represented in a binary file by:
    {BCSize, BC[0], ..., BC[BCSize-1], N, d[0], ..., d[N-1], D, Tensors[0], ..., Tensors[N-1]}
    where each Tensor is represented by:
    {rank, Shape[0], ..., Shape[rank-1], size, Elements[0], ..., Elements[size-1]}.
    \param FileName input: const string&, the binary file from which this MPS is read
    \sa void write(const string& FileName) const */
  void read(const string& FileName);

/// Converts float MPS to complex<float>.
/** This function converts this MPS of type float to complex<float>. This MPS must have T==float.
    \param MPS0 output: MPS< complex<float> >&, the complex equivalent to this MPS */
  void convertToComplex(MPS< complex<float> >& MPS0) const;

/// Converts double MPS to complex<double>.
/** This function converts this MPS of type double to complex<double>. This MPS must have T==double.
    \param MPS0 output: MPS< complex<double> >&, the complex equivalent to this MPS */
  void convertToComplex(MPS< complex<double> >& MPS0) const;

/// Returns this MPS as vector.
/** This MPS is returned as a ket vector in the standard basis.
    \param Vector0 output: vector<T>&, this MPS as a vector, must fulfill Vector0.size()==d^N */
  void getVector(vector<T>& Vector0) const;

/// This is filled with random entries.
/** This MPS is filled with uniformly distributed random entries. LAPACK's XLARNV is used.
    \param Seed input: const vector<unsigned int>&, must have Seed.size()==this->N and contains the seed
                       values for each tensor
    \sa void MPS<T>::fillZeroes() */
  void fillRandomly(const vector<unsigned int>& Seed);

/// Fills this MPS with zeroes.
/** This MPS is filled with zeroes.
    \sa void MPS<T>::fillRandomly() */
  void fillZeroes();

/// This is initialized as separable many-qubit state.
/** This MPS is initialized as a separable many-qubit state |psi>=|phi>*N where
    |phi>=Coefficients[0]|0>+Coefficients[1]|1>+...+Coefficients[d-1]|d-1>.
    \param Coefficients input: const vector<T>&, the coefficients in the standard basis */
  void setSeparable(const vector<T>& Coefficients);

/// Multiplies this MPS with element.
/** This MPS is multiplied with element, by multiplying the first tensor with element.
    \param element input: T, the scalar with which this MPS is multiplied */
  void multiply(T element);

/// Normalizes this MPS.
/** This function normalizes this MPS by normalizing each tensor with
       T Tensor<T>::normalize()   .
    \return T, the product of the absolute values of the largest elements of all tensors
    \sa T Tensor<T>::normalize() */
  T simplifiedNormalize();

/// Normalizes Tensor in this.
/** This function normalizes a Tensor in the course of bringing this MPS into normal shape or form. If
    this MPS is not in normal shape, this function can change the Shape of a Tensor according to its
    SVD. LAPACK's XGESVD is used.
    \param position input: unsigned int, the desired Tensor
    \param Direction input: const string&, the sweep direction, must be "left" or "right"
    \param X output: Matrix<T>&, the resulting gauge matrix
    \sa void bringIntoNormalShape()
    \sa void bringIntoNormalForm()
    \sa T normalize() */
  void normalizeTensor(unsigned int position, const string& Direction, Matrix<T>& X);

/// Brings this MPS into normal shape.
/** This MPS is brought into normal shape by normalizing each Tensor first from left to right until the
    right bond is smaller than the product of physical and left bond and then from right to left until
    the left bond is smaller than the product of physical and right bond.
    Since every MPS with open boundary conditions is assumed to be in normal shape, this function is only used by
       friend void multiplyMPOMPS<>(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                    unsigned int maxNumSweeps, double& errorAchieved,
                                    unsigned int& numSweepsDone, MPS<T>& MPS1)
    in the case of exact multiplication. It should always be used after setting new Tensors in this MPS with
       void set(unsigned int position, const Tensor<T>& Tensor0)   .
    \sa void normalizeTensor(unsigned int position, const string& Direction, Matrix<T>& X)
    \sa void bringIntoNormalForm()
    \sa T normalize() */
  void bringIntoNormalShape();

/// Brings this MPS into normal form.
/** Assuming that this MPS is in normal shape, it is brought into normal form by normalizing each
    tensor from right to left, such that this MPS is in normal form from right to left up to the first
    tensor.
    \sa void normalizeTensor(unsigned int position, const string& Direction, Matrix<T>& X)
    \sa void bringIntoNormalShape()
    \sa T normalize() */
  void bringIntoNormalForm();

/// Normalizes this MPS.
/** This MPS is normalized by bringing it into normal form and then normalizing the first tensor. The
    resulting norm of this MPS is returned.
    \return T, the norm
    \sa void normalizeTensor(unsigned int position, const string& Direction, Matrix<T>& X)
    \sa void bringIntoNormalShape()
    \sa void bringIntoNormalForm() */
  T normalize();

/// Scales MPS in canonical form.
/** This static function scales a MPS in canonical form, MPS0 and D0, by dividing each tensor by the
    absolute value of its largest element and each lambda matrix by the sum of its elements.
    The factor by which the MPS got scaled is returned.
    \param MPS0 input/output: MPS<T>&, the MPS in canonical form
    \param D0 input/output: vector< Matrix<T> >&, the lambda-matrices to MPS0
    \return T, the scaling factor */
  static T scale(MPS<T>& MPS0, vector< Matrix<T> >& D0);

/// Concatenates two MPSs.
/** This MPS and MPS0 are concatenated in MPS1:
    |MPS1> = |this->MPS> \otimes |MPS0>.
    \param MPS0 input: const MPS<T>&, the MPS to be put right to this MPS
    \param MPS1 output: MPS<T>&, the resulting MPS */
  void concatenate(const MPS<T>& MPS0, MPS<T>& MPS1) const;

/// Takes scalar product of this MPS with itself.
/** This function computes the scalar product of this MPS with itself.
    \return double, the resulting value for the scalar product */
  double scalarProduct() const;

/// Takes scalar product of this MPS with another MPS.
/** This function computes the scalar product of this MPS as bra with MPS0 as ket.
    \param MPS0 input: const MPS<T>&, the ket MPS
    \return T, the resulting value for the scalar product
    \sa double scalarProduct() const
    \sa T contractReverse(const MPS<T>& MPS0) const */
  T scalarProduct(const MPS<T>& MPS0) const;

/// Contracts this MPS in reverse order with another MPS.
/** This function computes the scalar product of this MPS with MPS0, taking the tensors of this MPS in
    reverse order and not complex conjugating them.
    It is used in the PEPS contraction.
    \param MPS0 input: const MPS<T>&, the ket MPS
    \return T, the resulting contraction value
    \sa friend T contractReverseMPOMPS<>(const MPS<T>& MPS0, const MPO<T>& MPO0, const MPS<T>& MPS1) */
  T contractReverse(const MPS<T>& MPS0) const;

/// This MPS is changed to approximate a superposition of MPSs.
/** This function changes this MPS to approximate a superposition of MPSs given in the vector MPSs.
    It evaluates the error of the approximation after each sweep and stops either
    if the change in the error is below eps (in sweep n: eps > abs(error(n) - error(n-1))) or
    if the number of sweeps exceeds maxNumSweeps. On input, this MPS and all MPSs in the vector
    MPSs must have the same BC, N, and d. On output, errorAchieved is the achieved approximation
    error and numSweepsDone is the number of sweeps done.
    \param MPSs input: const vector< MPS<T> >&, the superposition of MPSs, must all have the same
                       BC, N, and d
    \param eps input: double, the convergence precision, i.e. convergence in sweep n if
                      eps > abs(error(n) - error(n-1))
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param errorAchieved output: double&, the achieved approximation error
    \param numSweepsDone output: unsigned int&, the number of sweeps done */
  void approximate(const vector< MPS<T> >& MPSs, double eps, unsigned int maxNumSweeps,
                   double& errorAchieved, unsigned int& numSweepsDone);

/// Computes canonical form of this MPS.
/** This function computes the canonical form of this MPS and stores it in MPS0 and Lambdas.
    First the MPS is brought into normal form from right to left using
       void bringIntoNormalForm()   .
    Then each tensor is normalized with
       void normalizeTensor(unsigned int position, const string& Direction, Matrix<T>& X)
    from left to right, and X is singular value decomposed.
    \param eps input: double, the cutoff for the pseudoinverse, elements <= eps are assumed to be zero
    \param MPS0 output: MPS<T>&, the desired MPS in canonical form
    \param Lambdas output: vector< Matrix<T> >&, the lambdas
    \sa void normalizeTensor(unsigned int position, const string& Direction, Matrix<T>& X)
    \sa void bringIntoNormalForm()
    \sa void canonicalize2(double eps, MPS<T>& MPS0, vector< Matrix<T> >& Lambdas) const
    \sa void canonicalizeTEBD(double eps, unsigned int& numSweeps, MPS<T>& MPS0, vector< Matrix<T> >& D0) const
    \sa void decanonicalize(const vector< Matrix<T> >& Lambdas, MPS<T>& MPS0) const */
  void canonicalize(double eps, MPS<T>& MPS0, vector< Matrix<T> >& Lambdas) const;

/// Computes canonical form of this MPS.
/** This function computes the canonical form of this MPS and stores it in MPS0 and Lambdas.
    We use the algorithm described in
       [Y.-Y. Shi, L.-M. Duan, and G. Vidal, PRA 74, 022320 (2006)]   .
    \param eps input: double, the cutoff for the pseudoinverse, elements <= eps are assumed to be zero
    \param MPS0 output: MPS<T>&, the desired MPS in canonical form
    \param Lambdas output: vector< Matrix<T> >&, the lambdas
    \sa void canonicalize(double eps, MPS<T>& MPS0, vector< Matrix<T> >& Lambdas) const
    \sa void canonicalizeTEBD(double eps, unsigned int& numSweeps, MPS<T>& MPS0, vector< Matrix<T> >& D0) const
    \sa void decanonicalize(const vector< Matrix<T> >& Lambdas, MPS<T>& MPS0) const */
  void canonicalize2(double eps, MPS<T>& MPS0, vector< Matrix<T> >& Lambdas) const;

/// Computes canonical form of this MPS with TEBD.
/** This function computes the canonical form of this MPS and stores it in MPS0 and D0.
    The canonical form is obtained from TEBD with unit time evolution gates, and the time evolution stops
    if the change in the distance between all lambdas is below eps.
    \param eps input: double, desired relative precision of all lambdas
    \param numSweeps output: unsigned int&, the number of sweeps done
    \param MPS0 output: MPS<T>&, the desired MPS in canonical form
    \param D0 output: vector< Matrix<T> >&, the lambdas
    \sa void canonicalize(double eps, MPS<T>& MPS0, vector< Matrix<T> >& Lambdas) const
    \sa void canonicalize2(double eps, MPS<T>& MPS0, vector< Matrix<T> >& Lambdas) const
    \sa void decanonicalize(const vector< Matrix<T> >& Lambdas, MPS<T>& MPS0) const */
  void canonicalizeTEBD(double eps, unsigned int& numSweeps, MPS<T>& MPS0, vector< Matrix<T> >& D0) const;

/// Computes standard MPS from canonical form.
/** Given this MPS in canonical form with lambda-matrices Lambdas, its standard form, MPS0, is obtained by
    contracting each lambda with the tensor to its right.
    \param Lambdas input: const vector< Matrix<T> >&, the lambdas
    \param MPS0 output: MPS<T>&, the resulting MPS in standard form
    \sa void canonicalize(double eps, MPS<T>& MPS0, vector< Matrix<T> >& Lambdas) const
    \sa void canonicalize2(double eps, MPS<T>& MPS0, vector< Matrix<T> >& Lambdas) const
    \sa void canonicalizeTEBD(double eps, unsigned int& numSweeps, MPS<T>& MPS0, vector< Matrix<T> >& D0) const */
  void decanonicalize(const vector< Matrix<T> >& Lambdas, MPS<T>& MPS0) const;

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
/** This function approximates the product of the input MPS MPS0 with the input MPO MPO0 by a MPS MPS1
    with a given bond dimension. If (eps >= 0), it evaluates the error of the approximation after each sweep and
    stops if the change in the error is below eps, i.e. in sweep n: eps >= abs(error(n) - error(n-1)). The
    function always stops if the number of sweeps exceeds maxNumSweeps. Hereby the error is defined as
       error := ||MPO0*|MPS0>-|MPS1>||/||MPO0*|MPS0>||   .
    If (eps < 0.0), the error is never computed and maxNumSweeps sweeps are performed.
    On input, if (MPS1.D >= MPO0.D*MPS0.D) then the exact multiplication is done.
    On output, errorAchieved is the achieved approximation error, numSweepsDone is the number of sweeps
    done and MPS1 is the approximating MPS.
    \param MPO0 input: const MPO<T>&, the MPO
    \param MPS0 input: const MPS<T>&, the MPS on which the MPO acts
    \param eps input: double, if (eps >= 0.0) the convergence precision, i.e. convergence in sweep n if
                      eps >= abs(error(n) - error(n-1)),
                      if (eps < 0.0) no error computation is done
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param errorAchieved output: double&, the achieved approximation error, if (eps < 0.0) not addressed
    \param numSweepsDone output: unsigned int&, the number of sweeps done
    \param MPS1 input/output: MPS<T>&, on input the initial approximating MPS, on output the resulting
                              MPS, must have the correct form */
  friend void multiplyMPOMPS<>(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                               unsigned int maxNumSweeps, double& errorAchieved,
                               unsigned int& numSweepsDone, MPS<T>& MPS1);

  void getPurificationFromMPS(MPS<T>& MPS0, vector<unsigned int>& d2) const;

  void getPurificationEigenvalues(vector< vector<double> >& PuriEvals) const;

  void getMPSFromPurification(const vector<unsigned int>& d2, MPS<T>& MPS0) const;

  void getMPDOTensor(const vector<unsigned int>& d2, unsigned int position, Tensor<T>& Tensor0) const;

  void updateTensor(Tensor<T>& Tensor0, Tensor<T>& PuriTensor) const;

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
                         double cutoff, unsigned int mode, unsigned int d2,
                         double alphaStart, double x, double precision, unsigned int maxNumIter)
    and
       void updateTensor(unsigned int position, const string& Direction,
                         Tensor<T>& NTensorLeft, Tensor<T>& NTensorRight, Tensor<T>& bTensor,
                         double cutoff, unsigned int mode, unsigned int d2,
                         double alphaStart, double x, double precision, unsigned int maxNumIter)
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
                       of the pseudoinverse of the norm matrix,
                       if (mode == 2) then MPS1 is assumed to be a MPDO and the cubic equations are
                       solved iteratively
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
 
/// Updates boundary tensor.
/** This function updates a boundary tensor and is used by
       friend void multiplyMPOMPS<>(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                    unsigned int maxNumSweeps, double& errorAchieved,
                                    unsigned int& numSweepsDone, MPS<T>& MPS1,
                                    double cutoff, unsigned int mode, unsigned int d2,
                                    double alphaStart, double x, double precision, unsigned int maxNumIter)   .
    \param Direction input: const string&, the sweep direction
    \param NormTensor input: const Tensor<T>&, the norm tensor
    \param MTensor input: const Tensor<T>&, the M tensor from the right-hand side of the system of linear equations
    \param cutoff input: double, the cutoff for the pseudoinverse
    \param mode input: unsigned int, the mode,
                       if (mode == 0) then the norm matrix is not processed before a LAPACK routine
                       solves the system of linear equations,
                       if (mode == 1) the system of linear equations is solved via computation
                       of the pseudoinverse of the norm matrix,
                       if (mode == 2) then MPS1 is assumed to be a MPDO and the cubic equations are
                       solved iteratively
    \param d2 input: unsigned int, the maximal d2,
                     if (d2 == 0) then the tensors of MPS1 are not processed,
                     if (d2 != 0) then the tensors of MPS1 are made hermitian with rank d2
    \param alphaStart input: double, used only if (mode == 2)
    \param x input: double, used only if (mode == 2)
    \param precision input: double, used only if (mode == 2)
    \param maxNumIter input: unsigned int, used only if (mode == 2)
    \sa friend void multiplyMPOMPS<>(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                     unsigned int maxNumSweeps, double& errorAchieved,
                                     unsigned int& numSweepsDone, MPS<T>& MPS1,
                                     double cutoff, unsigned int mode, unsigned int d2,
                                     double alphaStart, double x, double precision, unsigned int maxNumIter)
    \sa void updateTensor(unsigned int position, const string& Direction,
                          const Tensor<T>& NormTensorLeft, const Tensor<T>& NormTensorRight,
                          const Tensor<T>& MTensor,
                          double cutoff, unsigned int mode, unsigned int d2,
                          double alphaStart, double x, double precision, unsigned int maxNumIter) */
  void updateTensor(const string& Direction, const Tensor<T>& NormTensor, const Tensor<T>& MTensor,
                    double cutoff, unsigned int mode, unsigned int d2,
                    double alphaStart, double x, double precision, unsigned int maxNumIter);
 
/// Updates bulk tensor.
/** This function updates a bulk tensor and is used by
       friend void multiplyMPOMPS<>(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                    unsigned int maxNumSweeps, double& errorAchieved,
                                    unsigned int& numSweepsDone, MPS<T>& MPS1,
                                    double cutoff, unsigned int mode, unsigned int d2,
                                    double alphaStart, double x, double precision, unsigned int maxNumIter)   .
    \param position input: unsigned int, the position
    \param Direction input: const string&, the sweep direction
    \param NormTensorLeft input: const Tensor<T>&, the left part of the norm tensor
    \param NormTensorRight input: const Tensor<T>&, the right part of the norm tensor
    \param MTensor input: const Tensor<T>&, the M tensor from the right-hand side of the system of linear equations
    \param cutoff input: double, the cutoff for the pseudoinverse
    \param mode input: unsigned int, the mode,
                       if (mode == 0) then the norm matrix is not processed before a LAPACK routine
                       solves the system of linear equations,
                       if (mode == 1) the system of linear equations is solved via computation
                       of the pseudoinverse of the norm matrix,
                       if (mode == 2) then MPS1 is assumed to be a MPDO and the cubic equations are
                       solved iteratively
    \param d2 input: unsigned int, the maximal d2,
                     if (d2 == 0) then the tensors of MPS1 are not processed,
                     if (d2 != 0) then the tensors of MPS1 are made hermitian with rank d2
    \param alphaStart input: double, used only if (mode == 2)
    \param x input: double, used only if (mode == 2)
    \param precision input: double, used only if (mode == 2)
    \param maxNumIter input: unsigned int, used only if (mode == 2)
    \sa friend void multiplyMPOMPS<>(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                     unsigned int maxNumSweeps, double& errorAchieved,
                                     unsigned int& numSweepsDone, MPS<T>& MPS1,
                                     double cutoff, unsigned int mode, unsigned int d2,
                                     double alphaStart, double x, double precision, unsigned int maxNumIter)
    \sa void updateTensor(const string& Direction, Tensor<T>& NTensor, Tensor<T>& bTensor,
                          double cutoff, unsigned int mode, unsigned int d2,
                          double alphaStart, double x, double precision, unsigned int maxNumIter) */
  void updateTensor(unsigned int position, const string& Direction,
                    const Tensor<T>& NormTensorLeft, const Tensor<T>& NormTensorRight,
                    const Tensor<T>& MTensor,
                    double cutoff, unsigned int mode, unsigned int d2,
                    double alphaStart, double x, double precision, unsigned int maxNumIter);

/// Multiplies even-odd MPO with MPS by means of TEBD.
/** This function approximates the product of an even-odd MPO with a MPS by a MPS with the same bond
    dimension. An even-odd MPO arises in the even-odd Trotter decomposition of the evolution operator.
    \param MPO0 input: const MPO<T>&, the even-odd MPO
    \param MPS0 input: const MPS<T>&, the MPS on which the even-odd MPO acts
    \param MPS1 output: MPS<T>&, the resulting MPS, must be of the same form as MPS0 */
  friend void multiplyEvenOddMPOMPS<>(const MPO<T>& MPO0, const MPS<T>& MPS0, MPS<T>& MPS1);

/// Returns MPS representation.
/** This friend function computes the MPS MPS0 with open BC representing the Tensor Tensor0.
    It performs successive singular value decompositions on Tensor0 with
       void singularValueDecompose(const vector<unsigned int>& Indices0,
                                   const vector<unsigned int>& Indices1,
                                   unsigned int Dcut, Tensor<T>& Tensor0, Matrix<T>& Sigma,
                                   Tensor<T>& Tensor1)   .
    The exact MPS representation is returned.
    \param Tensor0 input: Tensor<T>&, the Tensor
    \param MPS0 output: MPS<T>&, the MPS with open BC representing Tensor0
    \sa void singularValueDecompose(const vector<unsigned int>& Indices0,
                                    const vector<unsigned int>& Indices1,
                                    unsigned int Dcut, Tensor<T>& Tensor0, Matrix<T>& Sigma,
                                    Tensor<T>& Tensor1)
    \sa void getTensor(const MPS<T>& MPS0, Tensor<T>& Tensor0) */
  friend void getMPS<>(Tensor<T>& Tensor0, MPS<T>& MPS0);

/// Returns Tensor representation.
/** This friend function computes the Tensor Tensor0 representing the MPS MPS0.
    The range 1 boundary indices are removed in the case of open boundary conditions.
    \param MPS0 input: const MPS<T>&, the MPS
    \param Tensor0 output: Tensor<T>&, the Tensor representing MPS0
    \sa void getMPS(Tensor<T>& Tensor0, MPS<T>& MPS0) */
  friend void getTensor<>(const MPS<T>& MPS0, Tensor<T>& Tensor0);

/// Returns Shape to position for open boundary conditions.
/** Given a position in Tensors of this MPS with open boundary conditions, the correct Shape is
    returned. The returned Shape is given by the normal shape of this MPS.
    \param position input: unsigned int, the desired position
    \param Shape0 output: vector<unsigned int>&, the resulting Shape, must fulfill Shape0.size()==3 */
  inline void getOpenBCShape(unsigned int position, vector<unsigned int>& Shape0) const;

 protected:

/// Boundary conditions BC.
 string BC;

/// Number of tensors N.
 unsigned int N;

/// Physical dimensions d.
 vector<unsigned int> d;

/// Maximal virtual bond dimension D.
 unsigned int D;

/// Tensors.
 Tensor<T>* Tensors;
};

template<class T> MPS<T>::MPS()
{
 this->N = 0;
 this->D = 0;
 this->Tensors = 0;
}

template<class T> MPS<T>::MPS(string BC, unsigned int N, const vector<unsigned int>& d, unsigned int D)
{
#ifdef DEBUG
 if (((BC != "open") && (BC != "periodic")) || (N != d.size()))
 {
  cerr << "Program terminated because of error in constructor: " <<
          "template<class T> MPS<T>::" <<
          "MPS(string BC, unsigned int N, const vector<unsigned int>& d, unsigned int D): " <<
          "(((BC != open) && (BC != periodic)) || (N != d.size()))." << endl;
  exit(1);
 }
#endif
 this->BC = BC;
 this->N = N;
 this->d = d;
 this->D = D;
 this->Tensors = new Tensor<T>[this->N];
 vector<unsigned int> Shape(3);
 if (this->BC == "open")
 {
  for (int i = 0; i < this->N; i++)
  {
   this->getOpenBCShape(i, Shape);
   this->Tensors[i] = Tensor<T>(Shape);
  }
 }
 else if (this->BC == "periodic")
 {
  Shape[0] = this->D; Shape[1] = this->D;
  for (int i = 0; i < this->N; i++)
  {
   Shape[2] = this->d[i];
   this->Tensors[i] = Tensor<T>(Shape);
  }
 }
}

template<class T> MPS<T>::MPS(string BC, unsigned int N, unsigned int d, unsigned int D)
{
#ifdef DEBUG
 if ((BC != "open") && (BC != "periodic"))
 {
  cerr << "Program terminated because of error in constructor: " <<
          "template<class T> MPS<T>::" <<
          "MPS(string BC, unsigned int N, unsigned int d, unsigned int D): " <<
          "((BC != open) && (BC != periodic))" << endl;
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
 vector<unsigned int> Shape(3);
 if (this->BC == "open")
 {
  for (int i = 0; i < this->N; i++)
  {
   this->getOpenBCShape(i, Shape);
   this->Tensors[i] = Tensor<T>(Shape);
  }
 }
 else if (this->BC == "periodic")
 {
  Shape[0] = this->D; Shape[1] = this->D; Shape[2] = d;
  for (int i = 0; i < this->N; i++)
  {
   this->Tensors[i] = Tensor<T>(Shape);
  }
 }
}

template<class T> MPS<T>::MPS(const MPS<T>& MPS0)
{
 this->BC = MPS0.BC;
 this->N = MPS0.N;
 this->d = MPS0.d;
 this->D = MPS0.D;
 this->Tensors = new Tensor<T>[this->N];
 for (int i = 0; i < this->N; i++)
 {
  this->Tensors[i] = MPS0.Tensors[i];
 }
}

template<class T> MPS<T>::~MPS()
{
 delete[] this->Tensors;
}

template<class T> MPS<T>& MPS<T>::operator=(const MPS<T>& MPS0)
{
 if (this != &MPS0)
 {
  this->BC = MPS0.BC;
  this->N = MPS0.N;
  this->d = MPS0.d;
  this->D = MPS0.D;
  delete[] this->Tensors;
  this->Tensors = new Tensor<T>[this->N];
  for (int i = 0; i < this->N; i++)
  {
   this->Tensors[i] = MPS0.Tensors[i];
  }
 }
 return *this;
}

template<class T> void MPS<T>::setD(unsigned int D0, T element)
{
#ifdef DEBUG
 if ((this->N == 0) || (D0 == 0))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void MPS<T>::" <<
          "setD(unsigned int D0, T element = 0.0): " <<
          "((this->N == 0) || (D0 == 0))." << endl;
  exit(1);
 }
#endif
 if (D0 == this->D)
 {
  return;
 }
 vector<unsigned int> Shape0(3), Shape1(3), Index(3);
 Tensor<T> Tensor0, Tensor1;
 T element0;
 unsigned int D1 = this->D;
 this->D = D0;
 if (this->BC == "open")
 {
  for (int position = 0; position < this->N; position++)
  {
   this->getOpenBCShape(position, Shape0);
   Tensor0 = Tensor<T>(Shape0);
   if (element == 0.0)
   {
    Tensor0.fillZeroes();
   }
   else if (element != 0.0)
   {
    Tensor0.fillRandomly(time(0)+position*13, element);
   }
   this->get(position, Tensor1);
   Tensor1.getShape(Shape1);
   for (int k = 0; k < this->d[position]; k++)
   {
    Index[2] = k;
    for (int j = 0; j < min(Shape0[1], Shape1[1]); j++)
    {
     Index[1] = j;
     for (int i = 0; i < min(Shape0[0], Shape1[0]); i++)
     {
      Index[0] = i;
      element0 = Tensor1.get(Index);
      Tensor0.set(Index, element0);
     }
    }
   }
   this->set(position, Tensor0);
  }
 }
 else if (this->BC == "periodic")
 {
  Shape0[0] = D0; Shape0[1] = D0;
  for (int position = 0; position < this->N; position++)
  {
   Shape0[2] = this->d[position];
   Tensor0 = Tensor<T>(Shape0);
   if (element == 0.0)
   {
    Tensor0.fillZeroes();
   }
   else if (element != 0.0)
   {
    Tensor0.fillRandomly(time(0)+position*13, element);
   }
   this->get(position, Tensor1);
   for (int k = 0; k < this->d[position]; k++)
   {
    Index[2] = k;
    for (int j = 0; j < min(D0, D1); j++)
    {
     Index[1] = j;
     for (int i = 0; i < min(D0, D1); i++)
     {
      Index[0] = i;
      element0 = Tensor1.get(Index);
      Tensor0.set(Index, element0);
     }
    }
   }
   this->set(position, Tensor0);
  }
 }
}

template<class T> void MPS<T>::set(unsigned int position, const Tensor<T>& Tensor0)
{
#ifdef DEBUG
 if ((position > this->N-1) || (Tensor0.getRank() != 3))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "set(unsigned int position, const Tensor<T>& Tensor0): " <<
          "((position > this->N-1) || (Tensor0.getRank() != 3))." << endl;
  exit(1);
 }
#endif
 this->Tensors[position] = Tensor0;
}

template<class T> void MPS<T>::get(unsigned int position, Tensor<T>& Tensor0) const
{
#ifdef DEBUG
 if (position > this->N-1)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "get(unsigned int position, Tensor<T>& Tensor0) const: " <<
          "(position > this->N-1)." << endl;
  exit(1);
 }
#endif
 Tensor0 = this->Tensors[position];
}

template<class T> void MPS<T>::write(const string& FileName) const
{
 ofstream File(FileName.c_str(), ios::out | ios::binary);
 if (File.is_open())
 {
  unsigned int BCSize = this->BC.size();
  File.write((char*)&BCSize, sizeof(BCSize));
  for (int i = 0; i < BCSize; i++)
  {
   File.write((char*)&(this->BC[i]), sizeof(this->BC[i]));
  }
  File.write((char*)&(this->N), sizeof(this->N));
  for (int i = 0; i < this->N; i++)
  {
   File.write((char*)&(this->d[i]), sizeof(this->d[i]));
  }
  File.write((char*)&(this->D), sizeof(this->D));
  Tensor<T> Tensor0;
  unsigned int rank0;
  vector<unsigned int> Shape0;
  unsigned int size0;
  T element0;
  for (int i = 0; i < this->N; i++)
  {
   Tensor0 = this->Tensors[i];
   rank0 = Tensor0.getRank();
   Tensor0.getShape(Shape0);
   size0 = Tensor0.getSize();
   File.write((char*)&rank0, sizeof(rank0));
   for (int j = 0; j < rank0; j++)
   {
    File.write((char*)&(Shape0[j]), sizeof(Shape0[j]));
   }
   File.write((char*)&size0, sizeof(size0));
   for (int j = 0; j < size0; j++)
   {
    element0 = Tensor0.get(j);
    File.write((char*)&element0, sizeof(element0));
   }
  }
  File.close();
 }
 else
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "write(const string& FileName) const: " <<
          "Binary output file " << FileName << " could not be opened." << endl;
  exit(1);
 }
}

template<class T> void MPS<T>::read(const string& FileName)
{
 ifstream File(FileName.c_str(), ios::in | ios::binary);
 if (File.is_open())
 {
  unsigned int BCSize;
  File.read((char*)&BCSize, sizeof(BCSize));
  this->BC = string(BCSize, ' ');
  for (int i = 0; i < BCSize; i++)
  {
   File.read((char*)&(this->BC[i]), sizeof(this->BC[i]));
  }
  File.read((char*)&(this->N), sizeof(this->N));
  this->d = vector<unsigned int>(this->N);
  for (int i = 0; i < this->N; i++)
  {
   File.read((char*)&(this->d[i]), sizeof(this->d[i]));
  }
  File.read((char*)&(this->D), sizeof(this->D));
  delete[] this->Tensors;
  this->Tensors = new Tensor<T>[this->N];
  unsigned int rank0;
  vector<unsigned int> Shape0;
  Tensor<T> Tensor0;
  unsigned int size0;
  T element0;
  for (int i = 0; i < this->N; i++)
  {
   File.read((char*)&rank0, sizeof(rank0));
   Shape0 = vector<unsigned int>(rank0);
   for (int j = 0; j < rank0; j++)
   {
    File.read((char*)&(Shape0[j]), sizeof(Shape0[j]));
   }
   Tensor0 = Tensor<T>(Shape0);
   File.read((char*)&size0, sizeof(size0));
   for (int j = 0; j < size0; j++)
   {
    File.read((char*)&element0, sizeof(element0));
    Tensor0.set(j, element0);
   }
   this->Tensors[i] = Tensor0;
  }
  File.close();
 }
 else
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "read(const string& FileName): " <<
          "Binary input file " << FileName << " could not be opened." << endl;
  exit(1);
 }
}

template<class T> void MPS<T>::convertToComplex(MPS< complex<float> >& MPS0) const
{
#ifdef DEBUG
 if (typeid(T) != typeid(float))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "convertToComplex(MPS< complex<float> >& MPS0) const: " <<
          "(typeid(T) != typeid(float))." << endl;
  exit(1);
 }
#endif
 MPS0 = MPS< complex<float> >(this->BC, this->N, this->d, this->D);
 Tensor< complex<float> > Tensor0;
 for (int position = 0; position < this->N; position++)
 {
  this->Tensors[position].convertToComplex(Tensor0);
  MPS0.set(position, Tensor0);
 }
}

template<class T> void MPS<T>::convertToComplex(MPS< complex<double> >& MPS0) const
{
#ifdef DEBUG
 if (typeid(T) != typeid(double))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "convertToComplex(MPS< complex<double> >& MPS0) const: " <<
          "(typeid(T) != typeid(double))." << endl;
  exit(1);
 }
#endif
 MPS0 = MPS< complex<double> >(this->BC, this->N, this->d, this->D);
 Tensor< complex<double> > Tensor0;
 for (int position = 0; position < this->N; position++)
 {
  this->Tensors[position].convertToComplex(Tensor0);
  MPS0.set(position, Tensor0);
 }
}

template<class T> void MPS<T>::getVector(vector<T>& Vector0) const
{
 unsigned int dim = 1;
 for (int i = 0; i < this->N; i++)
 {
  dim *= this->d[i];
 }
#ifdef DEBUG
 if (Vector0.size() != dim)
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void MPS<T>::" <<
          "getVector(vector<T>& Vector0) const: " <<
          "Vector0 has an incorrect size." << endl;
  exit(1);
 }
#endif
 unsigned int dimRest, rest;
 vector<unsigned int> Indices(2), Index(1), Indices0(1), Indices1(1), Index2(2);
 Indices[0] = 0; Indices[1] = 1;
 Tensor<T> Tensor0, Tensor1;
 Indices0[0] = 1; Indices1[0] = 0;
 Index2[0] = 0; Index2[1] = 0;
 for (int i = 0; i < dim; i++)
 {
  dimRest = dim / d[0];
  rest = i;
  Index[0] = rest / dimRest;
  this->Tensors[0].getSubtensor(Indices, Index, Tensor0);
  rest %= dimRest;
  for (int j = 1; j < this->N; j++)
  {
   dimRest /= d[j];
   Index[0] = rest / dimRest;
   this->Tensors[j].getSubtensor(Indices, Index, Tensor1);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   rest %= dimRest;
  }
  if (this->BC == "open")
  {
   Vector0[i] = Tensor0.get(Index2);
  }
  else if (this->BC == "periodic")
  {}
 }
}

template<class T> void MPS<T>::fillRandomly(const vector<unsigned int>& Seed)
{
#ifdef DEBUG
 if ((this->N == 0) || (Seed.size() != this->N))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
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

template<class T> void MPS<T>::fillZeroes()
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void MPS<T>::" <<
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

template<class T> void MPS<T>::setSeparable(const vector<T>& Coefficients)
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "setSeparable(const vector<T>& Coefficients): " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
 this->fillZeroes();
 vector<unsigned int> Index(3);
 Index[0] = 0; Index[1] = 0;
 for (int i = 0; i < this->N; i++)
 {
  for (int j = 0; j < this->d[i]; j++)
  {
   Index[2] = j;
   this->Tensors[i].set(Index, Coefficients[j]);
  }
 }
}

template<class T> void MPS<T>::multiply(T element)
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "multiply(T element): " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
 this->Tensors[0].multiply(element);
}

template<class T> T MPS<T>::simplifiedNormalize()
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T MPS<T>::" <<
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

template<class T> void MPS<T>::normalizeTensor(unsigned int position, const string& Direction,
                                               Matrix<T>& X)
{
#ifdef DEBUG
 if ((this->BC != "open") || (this->N == 0) || (position > this->N-1) ||
     ((Direction != "left") && (Direction != "right")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "normalizeTensor(unsigned int position, const string& Direction, Matrix<T>& X): " <<
          "((this->BC != open) || (this->N == 0) || (position > this->N-1) || " <<
           "((Direction != left) && (Direction != right)))." << endl;
  exit(1);
 }
#endif
 vector<unsigned int> Shape(3), Index(3);
 this->Tensors[position].getShape(Shape);
 unsigned int dim0, dim1;
 T element;
 if (Direction == "left")
 {
  dim0 = Shape[0]; dim1 = Shape[1]*Shape[2];
  Matrix<T> A(dim0, dim1);
  for (int k = 0; k < Shape[2]; k++)
  {
   for (int j = 0; j < Shape[1]; j++)
   {
    for (int i = 0; i < Shape[0]; i++)
    {
     Index[0] = i; Index[1] = j; Index[2] = k;
     A(i, j+k*Shape[1]) = this->Tensors[position].get(Index);
    }
   }
  }
  Matrix<T> U(dim0, dim0), Sigma(dim0, dim1), Vt(dim1, dim1);
  A.singularValueDecompose(U, Sigma, Vt);
// if (dim0 <= dim1) the Tensor's Shape is unchanged:
  if (dim0 <= dim1)
  {
   for (int j = 0; j < dim1; j++)
   {
    for (int i = 0; i < dim0; i++)
    {
     element = Vt(i, j);
     Index[0] = i; Index[1] = j%Shape[1]; Index[2] = j/Shape[1];
     this->Tensors[position].set(Index, element);
    }
   }
   X = Matrix<T>(dim0, dim0);
   for (int j = 0; j < dim0; j++)
   {
    for (int i = 0; i < dim0; i++)
    {
     X(i, j) = U(i, j) * Sigma(j, j);
    }
   }
  }
// else if (dim0 > dim1) the Tensor's Shape is changed:
  else if (dim0 > dim1)
  {
   vector<unsigned int> Shape0(3);
   Shape0[0] = dim1; Shape0[1] = Shape[1]; Shape0[2] = Shape[2];
   Tensor<T> Tensor0(Shape0);
   for (int j = 0; j < dim1; j++)
   {
    for (int i = 0; i < dim1; i++)
    {
     element = Vt(i, j);
     Index[0] = i; Index[1] = j%Shape[1]; Index[2] = j/Shape[1];
     Tensor0.set(Index, element);
    }
   }
   this->set(position, Tensor0);
   X = Matrix<T>(dim0, dim1);
   for (int j = 0; j < dim1; j++)
   {
    for (int i = 0; i < dim0; i++)
    {
     X(i, j) = U(i, j) * Sigma(j, j);
    }
   }
  }
 }
 else if (Direction == "right")
 {
  dim0 = Shape[0]*Shape[2]; dim1 = Shape[1];
  Matrix<T> A(dim0, dim1);
  for (int k = 0; k < Shape[2]; k++)
  {
   for (int j = 0; j < Shape[1]; j++)
   {
    for (int i = 0; i < Shape[0]; i++)
    {
     Index[0] = i; Index[1] = j; Index[2] = k;
     A(i+k*Shape[0], j) = this->Tensors[position].get(Index);
    }
   }
  }
  Matrix<T> U(dim0, dim0), Sigma(dim0, dim1), Vt(dim1, dim1);
  A.singularValueDecompose(U, Sigma, Vt);
// if (dim0 >= dim1) the Tensor's Shape is unchanged:
  if (dim0 >= dim1)
  {
   for (int j = 0; j < dim1; j++)
   {
    for (int i = 0; i < dim0; i++)
    {
     element = U(i, j);
     Index[0] = i%Shape[0]; Index[1] = j; Index[2] = i/Shape[0];
     this->Tensors[position].set(Index, element);
    }
   }
   X = Matrix<T>(dim1, dim1);
   for (int j = 0; j < dim1; j++)
   {
    for (int i = 0; i < dim1; i++)
    {
     X(i, j) = Sigma(i, i) * Vt(i, j);
    }
   }
  }
// else if (dim0 < dim1) the Tensor's Shape is changed:
  else if (dim0 < dim1)
  {
   vector<unsigned int> Shape0(3);
   Shape0[0] = Shape[0]; Shape0[1] = dim0; Shape0[2] = Shape[2];
   Tensor<T> Tensor0(Shape0);
   for (int j = 0; j < dim0; j++)
   {
    for (int i = 0; i < dim0; i++)
    {
     element = U(i, j);
     Index[0] = i%Shape[0]; Index[1] = j; Index[2] = i/Shape[0];
     Tensor0.set(Index, element);
    }
   }
   this->set(position, Tensor0);
   X = Matrix<T>(dim0, dim1);
   for (int j = 0; j < dim1; j++)
   {
    for (int i = 0; i < dim0; i++)
    {
     X(i, j) = Sigma(i, i) * Vt(i, j);
    }
   }
  }
 }
}

template<class T> void MPS<T>::bringIntoNormalShape()
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "bringIntoNormalShape(): " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
 Matrix<T> X;
 vector<unsigned int> Shape(3), Indices(1), IndicesX(1), Order(3);
 unsigned int position;
// starting at position 0 we normalize tensors from left to right until the product of the left bond
// and the physical index becomes larger or equal to the right bond:
 position = 0;
 this->Tensors[position].getShape(Shape);
 Indices[0] = 0; IndicesX[0] = 1;
 Order[0] = 2; Order[1] = 0; Order[2] = 1;
 while ((position < this->N/2) && (Shape[0]*Shape[2] < Shape[1]))
 {
  this->normalizeTensor(position, "right", X);
  this->Tensors[position+1].contract(Indices, X, IndicesX);
  this->Tensors[position+1].permute(Order);
  position++;
  this->Tensors[position].getShape(Shape);
 }
// starting at position N-1 we normalize tensors from right to left until the product of the right bond
// and the physical index becomes larger or equal to the left bond:
 position = this->N-1;
 this->Tensors[position].getShape(Shape);
 Indices[0] = 1; IndicesX[0] = 0;
 Order[0] = 0; Order[1] = 2; Order[2] = 1;
 while ((position > this->N/2) && (Shape[0] > Shape[1]*Shape[2]))
 {
  this->normalizeTensor(position, "left", X);
  this->Tensors[position-1].contract(Indices, X, IndicesX);
  this->Tensors[position-1].permute(Order);
  position--;
  this->Tensors[position].getShape(Shape);
 }
}

template<class T> void MPS<T>::bringIntoNormalForm()
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "bringIntoNormalForm(): " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
 Matrix<T> X;
 vector<unsigned int> Indices(1), IndicesX(1), Order(3);
 Indices[0] = 1; IndicesX[0] = 0;
 Order[0] = 0; Order[1] = 2; Order[2] = 1;
 for (int i = this->N-1; i > 0; i--)
 {
  this->normalizeTensor(i, "left", X);
  this->Tensors[i-1].contract(Indices, X, IndicesX);
  this->Tensors[i-1].permute(Order);
 }
}

template<class T> T MPS<T>::normalize()
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T MPS<T>::" <<
          "normalize(): " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
 Matrix<T> X;
 this->bringIntoNormalForm();
 this->normalizeTensor(0, "left", X);
 return X(0, 0);
}

template<class T> T MPS<T>::scale(MPS<T>& MPS0, vector< Matrix<T> >& D0)
{
#ifdef DEBUG
 if ((MPS0.N == 0) || (D0.size() != MPS0.N-1))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> T MPS<T>::" <<
          "scale(MPS<T>& MPS0, vector< Matrix<T> >& D0): " <<
          "((MPS0.N == 0) || (D0.size() != MPS0.N-1))." << endl;
  exit(1);
 }
#endif
 T maxAbsVal, sum;
 T scalingFactor = 1.0;
// loop over all tensors of MPS0:
 for (int position = 0; position < MPS0.N; position++)
 {
  maxAbsVal = 0.0;
  for (int i = 0; i < MPS0.Tensors[position].getSize(); i++)
  {
   if (abs(MPS0.Tensors[position].get(i)) > maxAbsVal)
    maxAbsVal = abs(MPS0.Tensors[position].get(i));
  }
  MPS0.Tensors[position].multiply(1.0/maxAbsVal);
  scalingFactor /= maxAbsVal;
 }
// loop over all lambda-matrices in D0:
 for (int position = 0; position < MPS0.N-1; position++)
 {
  sum = 0.0;
  for (int i = 0; i < D0[position].getDim0(); i++)
   sum += D0[position](i, i);
  D0[position].multiply(1.0/sum);
  scalingFactor /= sum;
 }
// we return the scaling factor:
 return scalingFactor;
}

template<class T> void MPS<T>::concatenate(const MPS<T>& MPS0, MPS<T>& MPS1) const
{
#ifdef DEBUG
 if ((this->BC != MPS0.BC) || (this->N == 0) || (MPS0.N == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "concatenate(const MPS<T>& MPS0, MPS<T>& MPS1) const: " <<
          "((this->BC != MPS0.BC) || (this->N == 0) || (MPS0.N == 0))." << endl;
  exit(1);
 }
#endif
 string BC1 = this->BC;
 unsigned int N1 = this->N + MPS0.N;
 vector<unsigned int> d1(N1);
 for (int i = 0; i < this->N; i++)
 {
  d1[i] = this->d[i];
 }
 for (int i = this->N; i < N1; i++)
 {
  d1[i] = MPS0.d[i-this->N];
 }
 unsigned int D1 = max(this->D, MPS0.D);
 MPS1 = MPS<T>(BC1, N1, d1, D1);
 vector<unsigned int> Shape(3), Index(3), Shape0(3);
 Tensor<T> A, Tensor0;
 T element;
 unsigned int position;
 if (this->BC == "open")
 {
// A for 0 <= position < this->N:
  for (position = 0; position < this->N; position++)
  {
   MPS1.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   Tensor0 = this->Tensors[position];
   Tensor0.getShape(Shape0);
   for (int k = 0; k < Shape0[2]; k++)
   {
    Index[2] = k;
    for (int j = 0; j < Shape0[1]; j++)
    {
     Index[1] = j;
     for (int i = 0; i < Shape0[0]; i++)
     {
      Index[0] = i;
      element = Tensor0.get(Index);
      A.set(Index, element);
     }
    }
   }
   MPS1.set(position, A);
  }
// A for this->N <= position < N1:
  for (position = this->N; position < N1; position++)
  {
   MPS1.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   Tensor0 = MPS0.Tensors[position-this->N];
   Tensor0.getShape(Shape0);
   for (int k = 0; k < Shape0[2]; k++)
   {
    Index[2] = k;
    for (int j = 0; j < Shape0[1]; j++)
    {
     Index[1] = j;
     for (int i = 0; i < Shape0[0]; i++)
     {
      Index[0] = i;
      element = Tensor0.get(Index);
      A.set(Index, element);
     }
    }
   }
   MPS1.set(position, A);
  }
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void MPS<T>::" <<
          "concatenate(const MPS<T>& MPS0, MPS<T>& MPS1) const." << endl;
  exit(1);
 }
}

template<class T> double MPS<T>::scalarProduct() const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> double MPS<T>::" <<
          "scalarProduct() const: " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Indices0(1), Indices1(1), Indices02(2), Indices12(2);
 Tensor0 = this->Tensors[0]; this->Tensors[0].complexConjugate(Tensor1);
 Indices0[0] = 2; Indices1[0] = 2;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = this->Tensors[1];
 Indices0[0] = 1; Indices1[0] = 0;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 this->Tensors[1].complexConjugate(Tensor1);
 Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 2;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 for (int i = 1; i < this->N-1; i++)
 {
  Tensor1 = this->Tensors[i+1];
  Indices0[0] = 2; Indices1[0] = 0;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  this->Tensors[i+1].complexConjugate(Tensor1);
  Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 2;
  Tensor0.contract(Indices02, Tensor1, Indices12);
 }
 if (this->BC == "open")
 {
  vector<unsigned int> Index(4);
  Index[0] = 0; Index[1] = 0; Index[2] = 0; Index[3] = 0;
  return MathAuxiliary::convertToDouble(Tensor0.get(Index));
 }
 else if (this->BC == "periodic")
 {}
}

template<class T> T MPS<T>::scalarProduct(const MPS<T>& MPS0) const
{
#ifdef DEBUG
 if ((this->BC != MPS0.BC) || (this->N != MPS0.N) || (this->N == 0) || (this->d != MPS0.d))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T MPS<T>::" <<
          "scalarProduct(const MPS<T>& MPS0): " <<
          "This MPS and MPS0 are not of the same form or (this->N == 0)." << endl;
  exit(1);
 }
#endif
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Indices0(1), Indices1(1), Indices02(2), Indices12(2);
 Tensor0 = MPS0.Tensors[0]; this->Tensors[0].complexConjugate(Tensor1);
 Indices0[0] = 2; Indices1[0] = 2;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = MPS0.Tensors[1];
 Indices0[0] = 1; Indices1[0] = 0;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 this->Tensors[1].complexConjugate(Tensor1);
 Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 2;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 for (int i = 1; i < this->N-1; i++)
 {
  Tensor1 = MPS0.Tensors[i+1];
  Indices0[0] = 2; Indices1[0] = 0;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  this->Tensors[i+1].complexConjugate(Tensor1);
  Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 2;
  Tensor0.contract(Indices02, Tensor1, Indices12);
 }
 if (this->BC == "open")
 {
  vector<unsigned int> Index(4);
  Index[0] = 0; Index[1] = 0; Index[2] = 0; Index[3] = 0;
  return Tensor0.get(Index);
 }
 else if (this->BC == "periodic")
 {}
}

template<class T> T MPS<T>::contractReverse(const MPS<T>& MPS0) const
{
#ifdef DEBUG
 if ((this->N == 0) || (this->BC != MPS0.BC) || (this->N != MPS0.N))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T MPS<T>::" <<
          "contractReverse(const MPS<T>& MPS0) const: " <<
          "((this->N == 0) || (this->BC != MPS0.BC) || (this->N != MPS0.N))." << endl;
  exit(1);
 }
 vector<unsigned int> d0(this->N);
 for (int i = 0; i < this->N; i++)
  d0[i] = this->d[N-1-i];
 if (d0 != MPS0.d)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T MPS<T>::" <<
          "contractReverse(const MPS<T>& MPS0) const: " <<
          "The physical dimensions of this MPS and MPS0 do not match." << endl;
  exit(1);
 }
#endif
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Indices0(1), Indices1(1), Indices02(2), Indices12(2);
 Tensor0 = MPS0.Tensors[0]; Tensor1 = this->Tensors[this->N-1];
 Indices0[0] = 2; Indices1[0] = 2;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = MPS0.Tensors[1];
 Indices0[0] = 1; Indices1[0] = 0;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = this->Tensors[this->N-2];
 Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 for (int i = 2; i < this->N; i++)
 {
  Tensor1 = MPS0.Tensors[i];
  Indices0[0] = 2; Indices1[0] = 0;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor1 = this->Tensors[this->N-1-i];
  Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
  Tensor0.contract(Indices02, Tensor1, Indices12);
 }
 if (this->BC == "open")
 {
  vector<unsigned int> Index(4);
  Index[0] = 0; Index[1] = 0; Index[2] = 0; Index[3] = 0;
  return Tensor0.get(Index);
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented for periodic BC: " <<
          "template<class T> T MPS<T>::" <<
          "contractReverse(const MPS<T>& MPS0) const." << endl;
  exit(1);
 }
}

template<class T> void MPS<T>::approximate(const vector< MPS<T> >& MPSs, double eps,
                                           unsigned int maxNumSweeps,
                                           double& errorAchieved,
                                           unsigned int& numSweepsDone)
{
 unsigned int M = MPSs.size();
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "approximate(const vector< MPS<T> >& MPSs, double eps, " <<
                      "unsigned int maxNumSweeps, double& errorAchieved, " <<
                      "unsigned int& numSweepsDone): " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
 for (int i = 0; i < M; i++)
 {
  if ((MPSs[i].BC != this->BC) || (MPSs[i].N != this->N) || (MPSs[i].d != this->d))
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void MPS<T>::" <<
           "approximate(const vector< MPS<T> >& MPSs, double eps, " <<
                       "unsigned int maxNumSweeps, double& errorAchieved, " <<
                       "unsigned int& numSweepsDone): " <<
            "This MPS and MPSs have different BC, N, or d." << endl;
   exit(1);
  }
 }
#endif
 double error0; double error1 = 0.0;
 numSweepsDone = 0;
// During the sweep we need normSquaredMPSs=\sum_{j,k}<MPSs[j]|MPS[k]> for computing the error:
 double normSquaredMPSs = 0.0;
 for (int i = 0; i < M; i++)
 {
  for (int j = i+1; j < M; j++)
  {
   normSquaredMPSs += 2.0*MathAuxiliary::convertToDouble(MPSs[i].scalarProduct(MPSs[j]));
  }
 }
 for (int i = 0; i < M; i++)
 {
  normSquaredMPSs += MPSs[i].scalarProduct();
 }
 double normSquaredMPS0, realScalarProductMPS0MPSs;
 if (this->BC == "open")
 {
// First we bring this MPS into normal form from right to left:
  this->bringIntoNormalForm();
// Then we construct the b for each site and each MPS:
  vector< vector< Tensor<T> > > Tensorsb(M, vector< Tensor<T> >(this->N));
  Tensor<T> Tensor0, Tensor1, Tensor2, Tensorb;
  vector<unsigned int> Indices0(1), Indices1(1), Order(4), Order0(3), Indices02(2), Indices12(2);
  for (int i = 0; i < M; i++)
  {
   Tensor0 = MPSs[i].Tensors[this->N-1];
   this->Tensors[this->N-1].complexConjugate(Tensor1);
   Indices0[0] = 2; Indices1[0] = 2;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Order[0] = 1; Order[1] = 3; Order[2] = 0; Order[3] = 2;
   Tensor0.permute(Order);
   Tensorsb[i][this->N-1] = Tensor0;
   for (int j = this->N-2; j > 0; j--)
   {
    Tensor1 = MPSs[i].Tensors[j];
    Indices0[0] = 2; Indices1[0] = 1;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    this->Tensors[j].complexConjugate(Tensor1);
    Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    Tensorsb[i][j] = Tensor0;
   }
  }
  string Direction = "right";
  vector<unsigned int> Indices03(3), Index02(2), Index04(4);
  Matrix<T> X;
  while (numSweepsDone < maxNumSweeps)
  {
   if (Direction == "right")
   {
// we compute the error:
    Tensor0 = this->Tensors[0];
    this->Tensors[0].complexConjugate(Tensor1);
    Indices02[0] = 1; Indices02[1] = 2; Indices12[0] = 1; Indices12[1] = 2;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    Indices02[0] = 0; Indices02[1] = 0;
    normSquaredMPS0 = MathAuxiliary::convertToDouble(Tensor0.get(Indices02));
    realScalarProductMPS0MPSs = 0.0;
    for (int i = 0; i < M; i++)
    {
     Tensor0 = Tensorsb[i][1];
     Tensor1 = MPSs[i].Tensors[0];
     Indices0[0] = 2; Indices1[0] = 1;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     this->Tensors[0].complexConjugate(Tensor1);
     Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     Index04[0] = 0; Index04[1] = 0; Index04[2] = 0; Index04[3] = 0;
     realScalarProductMPS0MPSs += MathAuxiliary::convertToDouble(Tensor0.get(Index04));
    }
    error0 = sqrt(abs(normSquaredMPS0 - 2.0*realScalarProductMPS0MPSs + normSquaredMPSs));
    errorAchieved = error0;
    if (abs(error1 - error0) <= eps)
    {
     return;
    }
// we prepare the sweep to the right:
    Tensor0 = Tensorsb[0][1];
    Tensor1 = MPSs[0].Tensors[0];
    Indices0[0] = 2; Indices1[0] = 1;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Indices03[0] = 1; Indices03[1] = 2; Indices03[2] = 4;
    Index02[0] = 0; Index02[1] = 0;
    Tensor0.getSubtensor(Indices03, Index02, Tensorb);
    for (int i = 1; i < M; i++)
    {
     Tensor0 = Tensorsb[i][1];
     Tensor1 = MPSs[i].Tensors[0];
     Indices0[0] = 2; Indices1[0] = 1;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Indices03[0] = 1; Indices03[1] = 2; Indices03[2] = 4;
     Index02[0] = 0; Index02[1] = 0;
     Tensor0.getSubtensor(Indices03, Index02, Tensor1);
     Tensorb.add(Tensor1);
    }
    this->Tensors[0] = Tensorb;
    this->normalizeTensor(0, Direction, X);
    for (int i = 0; i < M; i++)
    {
     Tensor0 = MPSs[i].Tensors[0];
     this->Tensors[0].complexConjugate(Tensor1);
     Indices0[0] = 2; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Order[0] = 0; Order[1] = 2; Order[2] = 1; Order[3] = 3;
     Tensor0.permute(Order);
     Tensorsb[i][0] = Tensor0;
    }
    for (int i = 1; i < this->N-1; i++)
    {
     Tensor0 = Tensorsb[0][i-1];
     Tensor1 = MPSs[0].Tensors[i];
     Indices0[0] = 2; Indices1[0] = 0;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Tensor1 = Tensorsb[0][i+1];
     Indices0[0] = 3; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Indices03[0] = 2; Indices03[1] = 6; Indices03[2] = 3;
     Index04[0] = 0; Index04[1] = 0; Index04[2] = 0; Index04[3] = 0;
     Tensor0.getSubtensor(Indices03, Index04, Tensorb);
     for (int j = 1; j < M; j++)
     {
      Tensor0 = Tensorsb[j][i-1];
      Tensor1 = MPSs[j].Tensors[i];
      Indices0[0] = 2; Indices1[0] = 0;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      Tensor1 = Tensorsb[j][i+1];
      Indices0[0] = 3; Indices1[0] = 2;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      Indices03[0] = 2; Indices03[1] = 6; Indices03[2] = 3;
      Index04[0] = 0; Index04[1] = 0; Index04[2] = 0; Index04[3] = 0;
      Tensor0.getSubtensor(Indices03, Index04, Tensor1);
      Tensorb.add(Tensor1);
     }
     this->Tensors[i] = Tensorb;
     this->normalizeTensor(i, Direction, X);
     for (int j = 0; j < M; j++)
     {
      Tensor0 = Tensorsb[j][i-1];
      Tensor1 = MPSs[j].Tensors[i];
      Indices0[0] = 2; Indices1[0] = 0;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      this->Tensors[i].complexConjugate(Tensor1);
      Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      Tensorsb[j][i] = Tensor0;
     }
    }
    Indices0[0] = 0; Indices1[0] = 1;
    this->Tensors[this->N-1].contract(Indices0, X, Indices1);
    Order0[0] = 2; Order0[1] = 0; Order0[2] = 1;
    this->Tensors[this->N-1].permute(Order0);
    Direction = "left";
   }
   else if (Direction == "left")
   {
// we compute the error:
    Tensor0 = this->Tensors[this->N-1];
    this->Tensors[this->N-1].complexConjugate(Tensor1);
    Indices02[0] = 0; Indices02[1] = 2; Indices12[0] = 0; Indices12[1] = 2;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    Indices02[0] = 0; Indices02[1] = 0;
    normSquaredMPS0 = MathAuxiliary::convertToDouble(Tensor0.get(Indices02));
    realScalarProductMPS0MPSs = 0.0;
    for (int i = 0; i < M; i++)
    {
     Tensor0 = Tensorsb[i][this->N-2];
     Tensor1 = MPSs[i].Tensors[this->N-1];
     Indices0[0] = 2; Indices1[0] = 0;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     this->Tensors[this->N-1].complexConjugate(Tensor1);
     Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     Index04[0] = 0; Index04[1] = 0; Index04[2] = 0; Index04[3] = 0;
     realScalarProductMPS0MPSs += MathAuxiliary::convertToDouble(Tensor0.get(Index04));
    }
    error1 = sqrt(abs(normSquaredMPS0 - 2.0*realScalarProductMPS0MPSs + normSquaredMPSs));
    errorAchieved = error1;
    if (abs(error0 - error1) <= eps)
    {
     return;
    }
// we prepare the sweep to the left:
    Tensor0 = Tensorsb[0][this->N-2];
    Tensor1 = MPSs[0].Tensors[this->N-1];
    Indices0[0] = 2; Indices1[0] = 0;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Indices03[0] = 2; Indices03[1] = 1; Indices03[2] = 4;
    Index02[0] = 0; Index02[1] = 0;
    Tensor0.getSubtensor(Indices03, Index02, Tensorb);
    for (int i = 1; i < M; i++)
    {
     Tensor0 = Tensorsb[i][this->N-2];
     Tensor1 = MPSs[i].Tensors[this->N-1];
     Indices0[0] = 2; Indices1[0] = 0;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Indices03[0] = 2; Indices03[1] = 1; Indices03[2] = 4;
     Index02[0] = 0; Index02[1] = 0;
     Tensor0.getSubtensor(Indices03, Index02, Tensor1);
     Tensorb.add(Tensor1);
    }
    this->Tensors[this->N-1] = Tensorb;
    this->normalizeTensor(this->N-1, Direction, X);
    for (int i = 0; i < M; i++)
    {
     Tensor0 = MPSs[i].Tensors[this->N-1];
     this->Tensors[this->N-1].complexConjugate(Tensor1);
     Indices0[0] = 2; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Order[0] = 1; Order[1] = 3; Order[2] = 0; Order[3] = 2;
     Tensor0.permute(Order);
     Tensorsb[i][this->N-1] = Tensor0;
    }
    for (int i = this->N-2; i > 0; i--)
    {
     Tensor0 = Tensorsb[0][i+1];
     Tensor1 = MPSs[0].Tensors[i];
     Indices0[0] = 2; Indices1[0] = 1;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Tensor1 = Tensorsb[0][i-1];
     Indices0[0] = 3; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Indices03[0] = 6; Indices03[1] = 2; Indices03[2] = 3;
     Index04[0] = 0; Index04[1] = 0; Index04[2] = 0; Index04[3] = 0;
     Tensor0.getSubtensor(Indices03, Index04, Tensorb);
     for (int j = 1; j < M; j++)
     {
      Tensor0 = Tensorsb[j][i+1];
      Tensor1 = MPSs[j].Tensors[i];
      Indices0[0] = 2; Indices1[0] = 1;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      Tensor1 = Tensorsb[j][i-1];
      Indices0[0] = 3; Indices1[0] = 2;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      Indices03[0] = 6; Indices03[1] = 2; Indices03[2] = 3;
      Index04[0] = 0; Index04[1] = 0; Index04[2] = 0; Index04[3] = 0;
      Tensor0.getSubtensor(Indices03, Index04, Tensor1);
      Tensorb.add(Tensor1);
     }
     this->Tensors[i] = Tensorb;
     this->normalizeTensor(i, Direction, X);
     for (int j = 0; j < M; j++)
     {
      Tensor0 = Tensorsb[j][i+1];
      Tensor1 = MPSs[j].Tensors[i];
      Indices0[0] = 2; Indices1[0] = 1;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      this->Tensors[i].complexConjugate(Tensor1);
      Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      Tensorsb[j][i] = Tensor0;
     }
    }
    Indices0[0] = 1; Indices1[0] = 0;
    this->Tensors[0].contract(Indices0, X, Indices1);
    Order0[0] = 0; Order0[1] = 2; Order0[2] = 1;
    this->Tensors[0].permute(Order0);
    Direction = "right";
   }
   numSweepsDone++;
  }
 }
 else if (this->BC == "periodic")
 {}
}

template<class T> void MPS<T>::canonicalize(double eps, MPS<T>& MPS0, vector< Matrix<T> >& Lambdas) const
{
#ifdef DEBUG
 if ((this->N == 0) || (this->BC == "periodic"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "canonicalize(double eps, MPS<T>& MPS0, vector< Matrix<T> >& Lambdas) const: " <<
          "((this->N == 0) || (this->BC == periodic))." << endl;
  exit(1);
 }
#endif
 MPS0 = *this;
 Lambdas = vector< Matrix<T> >(MPS0.N-1);
// Bring MPS0 into normal form from right to left:
 MPS0.bringIntoNormalForm();
// Canonicalize from left to right:
 string Direction = "right";
 unsigned int position, dim;
 Matrix<T> X, U, Sigma, Vt, SigmaInv;
 vector<unsigned int> IndicesR(1), IndicesL(1), OrderR(3), OrderL(3);
 IndicesR[0] = 1; IndicesL[0] = 0;
 OrderR[0] = 0; OrderR[1] = 2; OrderR[2] = 1;
 OrderL[0] = 2; OrderL[1] = 0; OrderL[2] = 1;
// - position 0:
 position = 0;
 MPS0.normalizeTensor(position, Direction, X);
 dim = X.getDim0();
 U = Matrix<T>(dim, dim); Sigma = Matrix<T>(dim, dim); Vt = Matrix<T>(dim, dim);
 X.singularValueDecompose(U, Sigma, Vt);
 MPS0.Tensors[position].contract(IndicesR, U, IndicesL);
 MPS0.Tensors[position].permute(OrderR);
 Lambdas[position] = Sigma;
 MPS0.Tensors[position+1].contract(IndicesL, Vt, IndicesR);
 MPS0.Tensors[position+1].permute(OrderL);
// - position 1 to N-2:
 for (position = 1; position < MPS0.N-1; position++)
 {
  SigmaInv = Sigma;
  for (int i = 0; i < SigmaInv.getDim0(); i++)
  {
   if (Sigma(i, i) > eps)
    SigmaInv(i, i) = 1.0/Sigma(i, i);
   else
    SigmaInv(i, i) = 0.0;
  }
  MPS0.Tensors[position].contract(IndicesL, Sigma, IndicesR);
  MPS0.Tensors[position].permute(OrderL);
  MPS0.normalizeTensor(position, Direction, X);
  MPS0.Tensors[position].contract(IndicesL, SigmaInv, IndicesR);
  MPS0.Tensors[position].permute(OrderL);
  dim = X.getDim0();
  U = Matrix<T>(dim, dim); Sigma = Matrix<T>(dim, dim); Vt = Matrix<T>(dim, dim);
  X.singularValueDecompose(U, Sigma, Vt);
  MPS0.Tensors[position].contract(IndicesR, U, IndicesL);
  MPS0.Tensors[position].permute(OrderR);
  Lambdas[position] = Sigma;
  MPS0.Tensors[position+1].contract(IndicesL, Vt, IndicesR);
  MPS0.Tensors[position+1].permute(OrderL);
 }
}

template<class T> void MPS<T>::canonicalize2(double eps, MPS<T>& MPS0, vector< Matrix<T> >& Lambdas) const
{
#ifdef DEBUG
 if ((this->N == 0) || (this->BC == "periodic"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "canonicalize2(double eps, MPS<T>& MPS0, vector< Matrix<T> >& Lambdas) const: " <<
          "((this->N == 0) || (this->BC == periodic))." << endl;
  exit(1);
 }
#endif
 MPS0 = *this;
 Lambdas = vector< Matrix<T> >(MPS0.N-1);
// compute all matrices from the right and store them in Matrices:
 vector< Tensor<T> > Matrices(MPS0.N);
 Tensor<T> Tensor0, Tensor1;
 Matrix<T> Matrix0, Vr, Xt, Xtinv, Y, Yinv, U, Sigma, Vt;
 vector<double> W;
 vector<unsigned int> Indices0(1), Indices1(1), Indices02(2), Indices12(2), Order(4), Shape(4), Index(4);
 Index[0] = 0; Index[1] = 0;
 vector<unsigned int> Order3(3);
 unsigned int dim0;
// for position N-1:
 unsigned int position = MPS0.N-1;
 Tensor0 = MPS0.Tensors[position];
 Tensor0.complexConjugate(Tensor1);
 Indices0[0] = 2; Indices1[0] = 2;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Order[0] = 1; Order[1] = 3; Order[2] = 0; Order[3] = 2;
 Tensor0.permute(Order);
 Matrices[position] = Tensor0;
// for positions N-2 to 1:
 for (position = MPS0.N-2; position > 0; position--)
 {
  Tensor1 = MPS0.Tensors[position];
  Indices0[0] = 2; Indices1[0] = 1;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  MPS0.Tensors[position].complexConjugate(Tensor1);
  Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Matrices[position] = Tensor0;
 }
// for position 0 come from the left:
 position = 0;
 Tensor0 = MPS0.Tensors[position];
 Tensor0.complexConjugate(Tensor1);
 Indices0[0] = 2; Indices1[0] = 2;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Order[0] = 0; Order[1] = 2; Order[2] = 1; Order[3] = 3;
 Tensor0.permute(Order);
// canonicalize from left to right:
 for (position = 0; position < MPS0.N-1; position++)
 {
// compute Xt and its pseudoinverse Xtinv:
  Tensor0.getShape(Shape);
  dim0 = Shape[2];
  Matrix0 = Matrix<T>(dim0, dim0);
  for (int j = 0; j < dim0; j++)
  {
   Index[3] = j;
   for (int i = 0; i < dim0; i++)
   {
    Index[2] = i;
    Matrix0(i, j) = Tensor0.get(Index);
   }
  }
  Matrix0.setType("hermitian");
  W = vector<double>(dim0); Vr = Matrix<T>(dim0, dim0);
  Matrix0.eigenDecompose(W, Vr);
  Xt = Matrix<T>(dim0, dim0);
  for (int j = 0; j < dim0; j++)
  {
   for (int i = 0; i < dim0; i++)
   {
    Xt(i, j) = sqrt(abs(W[i]))*Vr(j, i);
   }
  }
  U = Matrix<T>(dim0, dim0); Sigma = Matrix<T>(dim0, dim0); Vt = Matrix<T>(dim0, dim0);
  Xt.singularValueDecompose(U, Sigma, Vt);
  for (int i = 0; i < dim0; i++)
  {
   if (Sigma(i, i) <= eps)
   {
    Sigma(i, i) = 0.0;
   }
  }
  Xt = U; Xt.multiply(Sigma); Xt.multiply(Vt);
  Vt.transpose(); Vt.complexConjugate();
  for (int i = 0; i < dim0; i++)
  {
   if (Sigma(i, i) > 0.0)
    Sigma(i, i) = 1.0/Sigma(i, i);
  }
  U.transpose(); U.complexConjugate();
  Xtinv = Vt; Xtinv.multiply(Sigma); Xtinv.multiply(U);
// and Y and its pseudoinverse Yinv:
  Tensor1 = Matrices[position+1];
  for (int j = 0; j < dim0; j++)
  {
   Index[3] = j;
   for (int i = 0; i < dim0; i++)
   {
    Index[2] = i;
    Matrix0(i, j) = Tensor1.get(Index);
   }
  }
  Matrix0.eigenDecompose(W, Vr);
  Y = Matrix<T>(dim0, dim0);
  for (int j = 0; j < dim0; j++)
  {
   for (int i = 0; i < dim0; i++)
   {
    Y(i, j) = Vr(i, j)*sqrt(abs(W[j]));
   }
  }
  Y.singularValueDecompose(U, Sigma, Vt);
  for (int i = 0; i < dim0; i++)
  {
   if (Sigma(i, i) <= eps)
    Sigma(i, i) = 0.0;
  }
  Y = U; Y.multiply(Sigma); Y.multiply(Vt);
  Vt.transpose(); Vt.complexConjugate();
  for (int i = 0; i < dim0; i++)
  {
   if (Sigma(i, i) > 0.0)
    Sigma(i, i) = 1.0/Sigma(i, i);
  }
  U.transpose(); U.complexConjugate();
  Yinv = Vt; Yinv.multiply(Sigma); Yinv.multiply(U);
// compute XtY and store it in Xt:
  Xt.multiply(Y);
// and singular value decompose:
  Xt.singularValueDecompose(U, Sigma, Vt);
  Xtinv.multiply(U); Vt.multiply(Yinv);
  Indices0[0] = 1; Indices1[0] = 0;
  MPS0.Tensors[position].contract(Indices0, Xtinv, Indices1);
  Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
  MPS0.Tensors[position].permute(Order3);
  Lambdas[position] = Sigma;
  Indices0[0] = 0; Indices1[0] = 1;
  MPS0.Tensors[position+1].contract(Indices0, Vt, Indices1);
  Order3[0] = 2; Order3[1] = 0; Order3[2] = 1;
  MPS0.Tensors[position+1].permute(Order3);
// finally compute the new left tensor:
  Tensor1 = this->Tensors[position+1];
  Indices0[0] = 2; Indices1[0] = 0;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  this->Tensors[position+1].complexConjugate(Tensor1);
  Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 2;
  Tensor0.contract(Indices02, Tensor1, Indices12);
 }
}

template<class T> void MPS<T>::canonicalizeTEBD(double eps, unsigned int& numSweeps,
                                                MPS<T>& MPS0, vector< Matrix<T> >& D0) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "canonicalizeTEBD(double eps, unsigned int& numSweeps, MPS<T>& MPS0, " <<
                           "vector< Matrix<T> >& D0) const: " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
 MPS0 = *this;
 D0 = vector< Matrix<T> >(this->N-1);
// initial lambdas are unity such that the state remains unchanged:
 vector<unsigned int> Shape(3);
 unsigned int position;
 for (position = 0; position < MPS0.N-1; position++)
 {
  MPS0.Tensors[position].getShape(Shape);
  D0[position] = Matrix<T>(Shape[1], Shape[1]);
  D0[position].fillZeroes();
  for (int i = 0; i < D0[position].getDim0(); i++)
   D0[position](i, i) = 1.0;
 }
 vector< Matrix<T> > D1(D0);
// TEBD with unit time evolution operators until change in the sum of all lambdas is below eps:
 string Decomposition = "even";
 double distance = 1.0;
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Indices0(1), Indices1(1), Order(4), ShapeL(3), ShapeR(3), Index(4), Index03(3);
 vector<unsigned int> Index13(3);
 Order[0] = 0; Order[1] = 1; Order[2] = 3; Order[3] = 2;
 unsigned int dim0, dim1;
 Matrix<T> TwoBodyMatrix, U, Sigma, Vt;
 numSweeps = 0;
 while (distance > eps)
 {
  if (Decomposition == "even")
  {
   for (position = 0; position < MPS0.N-1; position += 2)
   {
    if (position == 0)
    {
     Tensor0 = MPS0.Tensors[position];
    }
    else
    {
     Tensor0 = D0[position-1]; Tensor1 = MPS0.Tensors[position];
     Indices0[0] = 1; Indices1[0] = 0;
     Tensor0.contract(Indices0, Tensor1, Indices1);
    }
    Tensor1 = D0[position];
    Indices0[0] = 1; Indices1[0] = 0;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = MPS0.Tensors[position+1];
    Indices0[0] = 2; Indices1[0] = 0;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    if (position != MPS0.N-2)
    {
     Tensor1 = D0[position+1];
     Indices0[0] = 2; Indices1[0] = 0;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Tensor0.permute(Order);
    }
    MPS0.Tensors[position].getShape(ShapeL);
    dim0 = ShapeL[0]*ShapeL[2];
    MPS0.Tensors[position+1].getShape(ShapeR);
    dim1 = ShapeR[1]*ShapeR[2];
    TwoBodyMatrix = Matrix<T>(dim0, dim1); U = Matrix<T>(dim0, dim0); Sigma = Matrix<T>(dim0, dim1);
    Vt = Matrix<T>(dim1, dim1);
    for (int i = 0; i < dim0; i++)
    {
     Index[0] = i % ShapeL[0]; Index[1] = i / ShapeL[0];
     for (int j = 0; j < dim1; j++)
     {
      Index[2] = j % ShapeR[1]; Index[3] = j / ShapeR[1];
      TwoBodyMatrix(i, j) = Tensor0.get(Index);
     }
    }
    TwoBodyMatrix.singularValueDecompose(U, Sigma, Vt);
    if (position == 0)
    {
     for (int i = 0; i < dim0; i++)
     {
      Index03[0] = i % ShapeL[0]; Index03[2] = i / ShapeL[0];
      for (int j = 0; j < dim1; j++)
      {
       Index13[1] = j % ShapeR[1]; Index13[2] = j / ShapeR[1];
       for (int k = 0; k < ShapeL[1]; k++)
       {
        Index03[1] = k; Index13[0] = k;
        MPS0.Tensors[position].set(Index03, U(i, k));
        MPS0.Tensors[position+1].set(Index13, Vt(k, j)/D0[position+1](j, j));
        D0[position](k, k) = Sigma(k, k);
       }
      }
     }
    }
    else if (position == MPS0.N-2)
    {
     for (int i = 0; i < dim0; i++)
     {
      Index03[0] = i % ShapeL[0]; Index03[2] = i / ShapeL[0];
      for (int j = 0; j < dim1; j++)
      {
       Index13[1] = j % ShapeR[1]; Index13[2] = j / ShapeR[1];
       for (int k = 0; k < ShapeL[1]; k++)
       {
        Index03[1] = k; Index13[0] = k;
        MPS0.Tensors[position].set(Index03, U(i, k)/D0[position-1](i, i));
        MPS0.Tensors[position+1].set(Index13, Vt(k, j));
        D0[position](k, k) = Sigma(k, k);
       }
      }
     }
    }
    else
    {
     for (int i = 0; i < dim0; i++)
     {
      Index03[0] = i % ShapeL[0]; Index03[2] = i / ShapeL[0];
      for (int j = 0; j < dim1; j++)
      {
       Index13[1] = j % ShapeR[1]; Index13[2] = j / ShapeR[1];
       for (int k = 0; k < ShapeL[1]; k++)
       {
        Index03[1] = k; Index13[0] = k;
        MPS0.Tensors[position].set(Index03, U(i, k)/D0[position-1](i, i));
        MPS0.Tensors[position+1].set(Index13, Vt(k, j)/D0[position+1](j, j));
        D0[position](k, k) = Sigma(k, k);
       }
      }
     }
    }
   }
   Decomposition = "odd";
  }
  else if (Decomposition == "odd")
  {
   for (position = 1; position < MPS0.N-1; position += 2)
   {
    Tensor0 = D0[position-1]; Tensor1 = MPS0.Tensors[position];
    Indices0[0] = 1; Indices1[0] = 0;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = D0[position];
    Indices0[0] = 1; Indices1[0] = 0;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = MPS0.Tensors[position+1];
    Indices0[0] = 2; Indices1[0] = 0;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    if (position != MPS0.N-2)
    {
     Tensor1 = D0[position+1];
     Indices0[0] = 2; Indices1[0] = 0;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Tensor0.permute(Order);
    }
    MPS0.Tensors[position].getShape(ShapeL);
    dim0 = ShapeL[0]*ShapeL[2];
    MPS0.Tensors[position+1].getShape(ShapeR);
    dim1 = ShapeR[1]*ShapeR[2];
    TwoBodyMatrix = Matrix<T>(dim0, dim1); U = Matrix<T>(dim0, dim0); Sigma = Matrix<T>(dim0, dim1);
    Vt = Matrix<T>(dim1, dim1);
    for (int i = 0; i < dim0; i++)
    {
     Index[0] = i % ShapeL[0]; Index[1] = i / ShapeL[0];
     for (int j = 0; j < dim1; j++)
     {
      Index[2] = j % ShapeR[1]; Index[3] = j / ShapeR[1];
      TwoBodyMatrix(i, j) = Tensor0.get(Index);
     }
    }
    TwoBodyMatrix.singularValueDecompose(U, Sigma, Vt);
    if (position == MPS0.N-2)
    {
     for (int i = 0; i < dim0; i++)
     {
      Index03[0] = i % ShapeL[0]; Index03[2] = i / ShapeL[0];
      for (int j = 0; j < dim1; j++)
      {
       Index13[1] = j % ShapeR[1]; Index13[2] = j / ShapeR[1];
       for (int k = 0; k < ShapeL[1]; k++)
       {
        Index03[1] = k; Index13[0] = k;
        MPS0.Tensors[position].set(Index03, U(i, k)/D0[position-1](i, i));
        MPS0.Tensors[position+1].set(Index13, Vt(k, j));
        D0[position](k, k) = Sigma(k, k);
       }
      }
     }
    }
    else
    {
     for (int i = 0; i < dim0; i++)
     {
      Index03[0] = i % ShapeL[0]; Index03[2] = i / ShapeL[0];
      for (int j = 0; j < dim1; j++)
      {
       Index13[1] = j % ShapeR[1]; Index13[2] = j / ShapeR[1];
       for (int k = 0; k < ShapeL[1]; k++)
       {
        Index03[1] = k; Index13[0] = k;
        MPS0.Tensors[position].set(Index03, U(i, k)/D0[position-1](i, i));
        MPS0.Tensors[position+1].set(Index13, Vt(k, j)/D0[position+1](j, j));
        D0[position](k, k) = Sigma(k, k);
       }
      }
     }
    }
   }
   Decomposition = "even";
  }
// distance:
  distance = 0.0;
  for (position = 0; position < MPS0.N-1; position++)
  {
   for (int i = 0; i < D0[position].getDim0(); i++)
    distance += (D1[position](i, i)-D0[position](i, i))*(D1[position](i, i)-D0[position](i, i));
  }
  distance = sqrt(distance);
  D1 = D0;
  numSweeps++;
 }
}

template<class T> void MPS<T>::decanonicalize(const vector< Matrix<T> >& Lambdas, MPS<T>& MPS0) const
{
#ifdef DEBUG
 if ((this->N == 0) || (this->BC == "periodic") || (Lambdas.size() != this->N-1))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "decanonicalize(const vector< Matrix<T> >& Lambdas, MPS<T>& MPS0) const: " <<
          "((this->N == 0) || (this->BC == periodic) || (Lambdas.size() != this->N-1))." << endl;
  exit(1);
 }
#endif
 MPS0 = *this;
 unsigned int position;
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Indices0(1), Indices1(1);
 Indices0[0] = 1; Indices1[0] = 0;
 for (position = 1; position < MPS0.N; position++)
 {
  Tensor0 = Lambdas[position-1]; Tensor1 = MPS0.Tensors[position];
  Tensor0.contract(Indices0, Tensor1, Indices1);
  MPS0.Tensors[position] = Tensor0;
 }
}

template<class T> T expectationValueMPO(const MPS<T>& MPS0, const MPO<T>& MPO0)
{
#ifdef DEBUG
 if ((MPS0.BC != MPO0.BC) || (MPS0.N != MPO0.N) || (MPS0.N == 0) || (MPS0.d != MPO0.d))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> T " <<
          "expectationValueMPO(const MPS<T>& MPS0, const MPO<T>& MPO0): " <<
          "MPS0 and MPO0 are not of the same form or (MPS0.N == 0)." << endl;
  exit(1);
 }
#endif
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Indices0(1), Indices1(1), Indices02(2), Indices12(2);
 Tensor0 = MPS0.Tensors[0]; Tensor1 = MPO0.Tensors[0];
 Indices0[0] = 2; Indices1[0] = 2;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 MPS0.Tensors[0].complexConjugate(Tensor1);
 Indices0[0] = 4; Indices1[0] = 2;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = MPS0.Tensors[1];
 Indices0[0] = 1; Indices1[0] = 0;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = MPO0.Tensors[1];
 Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 MPS0.Tensors[1].complexConjugate(Tensor1);
 Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 for (int i = 1; i < MPS0.N-1; i++)
 {
  Tensor1 = MPS0.Tensors[i+1];
  Indices0[0] = 3; Indices1[0] = 0;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor1 = MPO0.Tensors[i+1];
  Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  MPS0.Tensors[i+1].complexConjugate(Tensor1);
  Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
  Tensor0.contract(Indices02, Tensor1, Indices12);
 }
 if (MPS0.BC == "open")
 {
  vector<unsigned int> Index(6);
  Index[0] = 0; Index[1] = 0; Index[2] = 0; Index[3] = 0; Index[4] = 0; Index[5] = 0;
  return Tensor0.get(Index);
 }
 else if (MPS0.BC == "periodic")
 {}
}

template<class T> T scalarProductMPOMPS(const MPS<T>& MPS0, const MPO<T>& MPO0, const MPS<T>& MPS1)
{
#ifdef DEBUG
 if ((MPS0.BC != MPO0.BC) || (MPS1.BC != MPO0.BC) || (MPS0.N != MPO0.N) || (MPS1.N != MPO0.N) ||
     (MPS0.N == 0) || (MPS0.d != MPO0.d) || (MPS1.d != MPO0.d))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> T " <<
          "scalarProductMPOMPS(const MPS<T>& MPS0, const MPO<T>& MPO0, const MPS<T>& MPS1): " <<
          "MPS0, MPO0 and MPS1 are not of the same form or (MPS0.N == 0)." << endl;
  exit(1);
 }
#endif
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Indices0(1), Indices1(1), Indices02(2), Indices12(2);
 Tensor0 = MPS1.Tensors[0]; Tensor1 = MPO0.Tensors[0];
 Indices0[0] = 2; Indices1[0] = 2;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 MPS0.Tensors[0].complexConjugate(Tensor1);
 Indices0[0] = 4; Indices1[0] = 2;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = MPS1.Tensors[1];
 Indices0[0] = 1; Indices1[0] = 0;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = MPO0.Tensors[1];
 Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 MPS0.Tensors[1].complexConjugate(Tensor1);
 Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 for (int i = 1; i < MPS0.N-1; i++)
 {
  Tensor1 = MPS1.Tensors[i+1];
  Indices0[0] = 3; Indices1[0] = 0;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor1 = MPO0.Tensors[i+1];
  Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  MPS0.Tensors[i+1].complexConjugate(Tensor1);
  Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
  Tensor0.contract(Indices02, Tensor1, Indices12);
 }
 if (MPS0.BC == "open")
 {
  vector<unsigned int> Index(6);
  Index[0] = 0; Index[1] = 0; Index[2] = 0; Index[3] = 0; Index[4] = 0; Index[5] = 0;
  return Tensor0.get(Index);
 }
 else if (MPS0.BC == "periodic")
 {}
}

template<class T> T contractReverseMPOMPS(const MPS<T>& MPS0, const MPO<T>& MPO0, const MPS<T>& MPS1)
{
#ifdef DEBUG
 if ((MPS0.N == 0) || (MPS0.BC != MPO0.BC) || (MPS0.BC != MPS1.BC) ||
     (MPS0.N != MPO0.N) || (MPS0.N != MPS1.N))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> T contractReverseMPOMPS(const MPS<T>& MPS0, const MPO<T>& MPO0, " <<
                                                    "const MPS<T>& MPS1) const: " <<
          "((MPS0.N == 0) || (MPS0.BC != MPO0.BC) || (MPS0.BC != MPS1.BC) || " <<
           "(MPS0.N != MPO0.N) || (MPS0.N != MPS1.N))." << endl;
  exit(1);
 }
 vector<unsigned int> d0(MPS0.N);
 for (int i = 0; i < MPS0.N; i++)
  d0[i] = MPS0.d[MPS0.N-1-i];
 if ((d0 != MPO0.d) || (d0 != MPS1.d))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> T contractReverseMPOMPS(const MPS<T>& MPS0, const MPO<T>& MPO0, " <<
                                                    "const MPS<T>& MPS1) const: " <<
          "The physical dimensions of MPS0, MPO0, and MPS1 do not match." << endl;
  exit(1);
 }
#endif
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Indices0(1), Indices1(1), Indices02(2), Indices12(2);
 Tensor0 = MPS1.Tensors[0]; Tensor1 = MPO0.Tensors[0];
 Indices0[0] = 2; Indices1[0] = 2;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = MPS0.Tensors[MPS0.N-1];
 Indices0[0] = 4; Indices1[0] = 2;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = MPS1.Tensors[1];
 Indices0[0] = 1; Indices1[0] = 0;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = MPO0.Tensors[1];
 Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 Tensor1 = MPS0.Tensors[MPS0.N-2];
 Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 for (int i = 2; i < MPS0.N; i++)
 {
  Tensor1 = MPS1.Tensors[i];
  Indices0[0] = 3; Indices1[0] = 0;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor1 = MPO0.Tensors[i];
  Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor1 = MPS0.Tensors[MPS0.N-1-i];
  Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
  Tensor0.contract(Indices02, Tensor1, Indices12);
 }
 if (MPS0.BC == "open")
 {
  vector<unsigned int> Index(6);
  Index[0] = 0; Index[1] = 0; Index[2] = 0; Index[3] = 0; Index[4] = 0; Index[5] = 0;
  return Tensor0.get(Index);
 }
 else if (MPS0.BC == "periodic")
 {
  cerr << "The following friend function is not implemented for periodic BC: " <<
          "template<class T> T contractReverseMPOMPS(const MPS<T>& MPS0, const MPO<T>& MPO0, " <<
                                                    "const MPS<T>& MPS1) const." << endl;
  exit(1);
 }
}

template<class T> void exactMultiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, MPS<T>& MPS1)
{
#ifdef DEBUG
 if ((MPO0.BC != MPS0.BC) || (MPO0.N != MPS0.N) || (MPO0.N == 0) || (MPO0.d != MPS0.d))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> void " <<
          "exactMultiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, MPS<T>& MPS1): " <<
          "MPO0 and MPS0 are not of the same form or (MPO0.N == 0)." << endl;
  exit(1);
 }
#endif
// define MPS1:
 MPS1 = MPS<T>(MPO0.BC, MPO0.N, MPO0.d, MPO0.D*MPS0.D);
// compute tensors for MPS1:
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Indices0(1), Indices1(1), Order(5), ShapeMPO0(4), ShapeMPS0(3), Shape(3);
 Indices0[0] = 2; Indices1[0] = 2;
 Order[0] = 0; Order[1] = 2; Order[2] = 1; Order[3] = 3; Order[4] = 4;
 for (int pos = 0; pos < MPO0.N; pos++)
 {
  Tensor0 = MPS0.Tensors[pos]; Tensor1 = MPO0.Tensors[pos];
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order);
  MPS0.Tensors[pos].getShape(ShapeMPS0); MPO0.Tensors[pos].getShape(ShapeMPO0);
  Shape[0] = ShapeMPS0[0]*ShapeMPO0[0];
  Shape[1] = ShapeMPS0[1]*ShapeMPO0[1];
  Shape[2] = ShapeMPO0[3];
  Tensor0.reshape(Shape);
  MPS1.Tensors[pos] = Tensor0;
 }
// in the case of open boundary conditions recover the normal shape:
 if (MPO0.BC == "open")
  MPS1.bringIntoNormalShape();
}

template<class T> double distanceMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, const MPS<T>& MPS1)
{
#ifdef DEBUG
 if ((MPO0.BC != MPS0.BC) || (MPO0.BC != MPS1.BC) || (MPO0.N != MPS0.N) || (MPO0.N != MPS1.N) ||
     (MPO0.N == 0) || (MPO0.d != MPS0.d) || (MPO0.d != MPS1.d))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> double " <<
          "distanceMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, const MPS<T>& MPS1): " <<
          "MPO0, MPS0 and MPS1 are not of the same form or (MPO0.N == 0)." << endl;
  exit(1);
 }
#endif
// normSquaredMPO0MPS0 = <MPS0|MPO**{+}MPO|MPS0>:
 MPS<T> MPS2;
 exactMultiplyMPOMPS(MPO0, MPS0, MPS2);
 double normSquaredMPO0MPS0 = MPS2.scalarProduct();
// realScalarProductMPS1MPO0MPS0 = real(<MPS1|MPO0|MPS0>):
 double realScalarProductMPS1MPO0MPS0 = MathAuxiliary::convertToDouble(scalarProductMPOMPS(MPS1, MPO0,
                                                                                           MPS0));
// normSquaredMPS1 = <MPS1|MPS1>:
 double normSquaredMPS1 = MPS1.scalarProduct();
// returning the resulting distance:
 return sqrt(abs(normSquaredMPO0MPS0 - 2.0*realScalarProductMPS1MPO0MPS0 + normSquaredMPS1));
}

template<class T> void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                      unsigned int maxNumSweeps, double& errorAchieved,
                                      unsigned int& numSweepsDone, MPS<T>& MPS1)
{
#ifdef DEBUG
 if ((MPO0.BC != MPS0.BC) || (MPO0.BC != MPS1.BC) || (MPO0.N != MPS0.N) || (MPO0.N != MPS1.N) ||
     (MPO0.N == 0) || (MPO0.d != MPS0.d) || (MPO0.d != MPS1.d))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> void " <<
          "multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps, " <<
                         "unsigned int maxNumSweeps, double& errorAchieved, " <<
                         "unsigned int& numSweepsDone, MPS<T>& MPS1): " <<
          "MPO0, MPS0 and MPS1 are not of the same form or (MPO0.N == 0)." << endl;
  exit(1);
 }
#endif
// Variables for error computation:
 double normSquaredMPO0MPS0, realScalarProductMPS1MPO0MPS0, normSquaredMPS1;
 if (eps >= 0.0)
 {
// During the sweep we need normSquaredMPO0MPS0=<MPS0|MPO0^{+}MPO0|MPS0> for computing the error:
  MPO<T> MPO0Adjoint, MPO0Product;
  MPO0.adjoint(MPO0Adjoint);
  MPO0Adjoint.multiply(MPO0, MPO0Product);
  normSquaredMPO0MPS0 = MathAuxiliary::convertToDouble(expectationValueMPO(MPS0, MPO0Product));
 }
 if (MPS1.D >= MPO0.D*MPS0.D)
 {
  Tensor<T> Tensor0, Tensor1;
  vector<unsigned int> Indices0(1), Indices1(1), Order(5), ShapeMPO0(4), ShapeMPS0(3), Shape(3);
  Indices0[0] = 2; Indices1[0] = 2;
  Order[0] = 0; Order[1] = 2; Order[2] = 1; Order[3] = 3; Order[4] = 4;
  for (int i = 0; i < MPO0.N; i++)
  {
   Tensor0 = MPS0.Tensors[i]; Tensor1 = MPO0.Tensors[i];
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order);
   MPS0.Tensors[i].getShape(ShapeMPS0); MPO0.Tensors[i].getShape(ShapeMPO0);
   Shape[0] = ShapeMPS0[0]*ShapeMPO0[0];
   Shape[1] = ShapeMPS0[1]*ShapeMPO0[1];
   Shape[2] = ShapeMPO0[3];
   Tensor0.reshape(Shape);
   MPS1.Tensors[i] = Tensor0;
  }
// in the case of open boundary conditions we have to recover the normal shape:
  if (MPO0.BC == "open")
   MPS1.bringIntoNormalShape();
// we have to reset D:
  MPS1.setD(MPS1.D);
  if (eps >= 0.0)
   errorAchieved = distanceMPOMPS(MPO0, MPS0, MPS1)/sqrt(abs(normSquaredMPO0MPS0));
  numSweepsDone = 0;
 }
 else
 {
  double error0, error1;
  if (eps >= 0.0)
   error1 = 0.0;
  numSweepsDone = 0;
  if (MPO0.BC == "open")
  {
// First we bring MPS1 into normal form from right to left:
   MPS1.bringIntoNormalForm();
// Then we construct the b for each site:
   Tensor<T> Tensorsb[MPO0.N];
   Tensor<T> Tensor0, Tensor1, Tensor2;
   vector<unsigned int> Indices0(1), Indices1(1), Order(6), Indices02(2), Indices12(2);
   Tensor0 = MPS0.Tensors[MPO0.N-1]; Tensor1 = MPO0.Tensors[MPO0.N-1];
   Indices0[0] = 2; Indices1[0] = 2;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   MPS1.Tensors[MPO0.N-1].complexConjugate(Tensor1);
   Indices0[0] = 4; Indices1[0] = 2;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Order[0] = 1; Order[1] = 3; Order[2] = 5; Order[3] = 0; Order[4] = 2; Order[5] = 4;
   Tensor0.permute(Order);
   Tensorsb[MPO0.N-1] = Tensor0;
   for (int i = MPO0.N-2; i > 0; i--)
   {
    Tensor1 = MPS0.Tensors[i];
    Indices0[0] = 3; Indices1[0] = 1;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = MPO0.Tensors[i];
    Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    MPS1.Tensors[i].complexConjugate(Tensor1);
    Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    Tensorsb[i] = Tensor0;
   }
   string Direction = "right";
   vector<unsigned int> Indices03(3), Index04(4), Index06(6), Order0(3);
   Matrix<T> X;
   while (numSweepsDone < maxNumSweeps)
   {
    if (Direction == "right")
    {
     Tensor0 = Tensorsb[1]; Tensor1 = MPS0.Tensors[0];
     Indices0[0] = 3; Indices1[0] = 1;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Tensor1 = MPO0.Tensors[0];
     Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     if (eps >= 0.0)
     {
// we compute the error:
      Tensor1 = Tensor0;
      MPS1.Tensors[0].complexConjugate(Tensor2);
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
      Tensor1.contract(Indices02, Tensor2, Indices12);
      Index06[0] = 0; Index06[1] = 0; Index06[2] = 0; Index06[3] = 0; Index06[4] = 0; Index06[5] = 0;
      realScalarProductMPS1MPO0MPS0 = MathAuxiliary::convertToDouble(Tensor1.get(Index06));
      Tensor1 = MPS1.Tensors[0]; MPS1.Tensors[0].complexConjugate(Tensor2);
      Indices02[0] = 1; Indices02[1] = 2; Indices12[0] = 1; Indices12[1] = 2;
      Tensor1.contract(Indices02, Tensor2, Indices12);
      Indices02[0] = 0; Indices02[1] = 0;
      normSquaredMPS1 = MathAuxiliary::convertToDouble(Tensor1.get(Indices02));
      error0 = sqrt(abs(normSquaredMPO0MPS0 - 2.0*realScalarProductMPS1MPO0MPS0 + normSquaredMPS1));
      if (abs(error1 - error0)/sqrt(abs(normSquaredMPO0MPS0)) <= eps)
      {
       errorAchieved = error0/sqrt(abs(normSquaredMPO0MPS0));
       return;
      }
     }
// we prepare the sweep to the right:
     Indices03[0] = 2; Indices03[1] = 3; Indices03[2] = 6;
     Index04[0] = 0; Index04[1] = 0; Index04[2] = 0; Index04[3] = 0;
     Tensor0.getSubtensor(Indices03, Index04, Tensor1);
     MPS1.Tensors[0] = Tensor1;
     MPS1.normalizeTensor(0, Direction, X);
     Tensor0 = MPS0.Tensors[0]; Tensor1 = MPO0.Tensors[0];
     Indices0[0] = 2; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     MPS1.Tensors[0].complexConjugate(Tensor1);
     Indices0[0] = 4; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Order[0] = 0; Order[1] = 2; Order[2] = 4; Order[3] = 1; Order[4] = 3; Order[5] = 5;
     Tensor0.permute(Order);
     Tensorsb[0] = Tensor0;
     for (int i = 1; i < MPO0.N-1; i++)
     {
      Tensor1 = MPS0.Tensors[i];
      Indices0[0] = 3; Indices1[0] = 0;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      Tensor1 = MPO0.Tensors[i];
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      Tensor2 = Tensor0;
      Tensor1 = Tensorsb[i+1];
      Indices02[0] = 4; Indices02[1] = 5; Indices12[0] = 3; Indices12[1] = 4;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      Indices03[0] = 3; Indices03[1] = 8; Indices03[2] = 4;
      Index06[0] = 0; Index06[1] = 0; Index06[2] = 0; Index06[3] = 0; Index06[4] = 0; Index06[5] = 0;
      Tensor0.getSubtensor(Indices03, Index06, Tensor1);
      MPS1.Tensors[i] = Tensor1;
      MPS1.normalizeTensor(i, Direction, X);
      Tensor0 = Tensor2;
      MPS1.Tensors[i].complexConjugate(Tensor1);
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      Tensorsb[i] = Tensor0;
     }
     Indices0[0] = 0; Indices1[0] = 1;
     MPS1.Tensors[MPO0.N-1].contract(Indices0, X, Indices1);
     Order0[0] = 2; Order0[1] = 0; Order0[2] = 1;
     MPS1.Tensors[MPO0.N-1].permute(Order0);
     Direction = "left";
    }
    else if (Direction == "left")
    {
     Tensor0 = Tensorsb[MPO0.N-2]; Tensor1 = MPS0.Tensors[MPO0.N-1];
     Indices0[0] = 3; Indices1[0] = 0;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Tensor1 = MPO0.Tensors[MPO0.N-1];
     Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     if (eps >= 0.0)
     {
// we compute the error:
      Tensor1 = Tensor0;
      MPS1.Tensors[MPO0.N-1].complexConjugate(Tensor2);
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
      Tensor1.contract(Indices02, Tensor2, Indices12);
      Index06[0] = 0; Index06[1] = 0; Index06[2] = 0; Index06[3] = 0; Index06[4] = 0; Index06[5] = 0;
      realScalarProductMPS1MPO0MPS0 = MathAuxiliary::convertToDouble(Tensor1.get(Index06));
      Tensor1 = MPS1.Tensors[MPO0.N-1]; MPS1.Tensors[MPO0.N-1].complexConjugate(Tensor2);
      Indices02[0] = 0; Indices02[1] = 2; Indices12[0] = 0; Indices12[1] = 2;
      Tensor1.contract(Indices02, Tensor2, Indices12);
      Indices02[0] = 0; Indices02[1] = 0;
      normSquaredMPS1 = MathAuxiliary::convertToDouble(Tensor1.get(Indices02));
      error1 = sqrt(abs(normSquaredMPO0MPS0 - 2.0*realScalarProductMPS1MPO0MPS0 + normSquaredMPS1));
      if (abs(error0 - error1)/sqrt(abs(normSquaredMPO0MPS0)) <= eps)
      {
       errorAchieved = error1/sqrt(abs(normSquaredMPO0MPS0));
       return;
      }
     }
// we prepare the sweep to the left:
     Indices03[0] = 3; Indices03[1] = 2; Indices03[2] = 6;
     Index04[0] = 0; Index04[1] = 0; Index04[2] = 0; Index04[3] = 0;
     Tensor0.getSubtensor(Indices03, Index04, Tensor1);
     MPS1.Tensors[MPO0.N-1] = Tensor1;
     MPS1.normalizeTensor(MPO0.N-1, Direction, X);
     Tensor0 = MPS0.Tensors[MPO0.N-1]; Tensor1 = MPO0.Tensors[MPO0.N-1];
     Indices0[0] = 2; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     MPS1.Tensors[MPO0.N-1].complexConjugate(Tensor1);
     Indices0[0] = 4; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Order[0] = 1; Order[1] = 3; Order[2] = 5; Order[3] = 0; Order[4] = 2; Order[5] = 4;
     Tensor0.permute(Order);
     Tensorsb[MPO0.N-1] = Tensor0;
     for (int i = MPO0.N-2; i > 0; i--)
     {
      Tensor1 = MPS0.Tensors[i];
      Indices0[0] = 3; Indices1[0] = 1;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      Tensor1 = MPO0.Tensors[i];
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      Tensor2 = Tensor0;
      Tensor1 = Tensorsb[i-1];
      Indices02[0] = 4; Indices02[1] = 5; Indices12[0] = 3; Indices12[1] = 4;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      Indices03[0] = 8; Indices03[1] = 3; Indices03[2] = 4;
      Index06[0] = 0; Index06[1] = 0; Index06[2] = 0; Index06[3] = 0; Index06[4] = 0; Index06[5] = 0;
      Tensor0.getSubtensor(Indices03, Index06, Tensor1);
      MPS1.Tensors[i] = Tensor1;
      MPS1.normalizeTensor(i, Direction, X);
      Tensor0 = Tensor2;
      MPS1.Tensors[i].complexConjugate(Tensor1);
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      Tensorsb[i] = Tensor0;
     }
     Indices0[0] = 1; Indices1[0] = 0;
     MPS1.Tensors[0].contract(Indices0, X, Indices1);
     Order0[0] = 0; Order0[1] = 2; Order0[2] = 1;
     MPS1.Tensors[0].permute(Order0);
     Direction = "right";
    }
    numSweepsDone++;
   }
   if (eps >= 0.0)
    errorAchieved = distanceMPOMPS(MPO0, MPS0, MPS1)/sqrt(abs(normSquaredMPO0MPS0));
  }
  else if (MPO0.BC == "periodic")
  {
   cout << "The following friend function is not implemented for periodic boundary conditions: " <<
           "template<class T> void " <<
           "multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps, " <<
                          "unsigned int maxNumSweeps, double& errorAchieved, " <<
                          "unsigned int& numSweepsDone, MPS<T>& MPS1)." << endl;
  }
 }
}

template<class T> void MPS<T>::getPurificationFromMPS(MPS<T>& MPS0, vector<unsigned int>& d2) const
{
// define Purification-MPS MPS0:
 unsigned int sqrtD = sqrt(double(this->D));
 MPS0 = MPS<T>(this->BC, this->N, this->d, sqrtD);
 d2 = vector<unsigned int>(this->N);
// fill MPS0:
 unsigned int dim, sqrtDim;
 vector<unsigned int> Shape0(3), SqrtShape0(3), Shape6(6), Order6(6), Shape2(2), Index2(2);
 vector<unsigned int> AShape(4), Shape3(3);
 vector<double> W; Matrix<T> Vr;
 Tensor<T> Tensor0, ATensor;
 Matrix<T> Matrix0;
 Order6[0] = 0; Order6[1] = 2; Order6[2] = 4; Order6[3] = 1; Order6[4] = 3; Order6[5] = 5;
 for (unsigned int position = 0; position < this->N; position++)
 {
  this->get(position, Tensor0);
  Tensor0.getShape(Shape0);
  dim = Shape0[0]*Shape0[1]*Shape0[2];
  sqrtDim = sqrt(double(dim));
  SqrtShape0[0] = sqrt(double(Shape0[0])); SqrtShape0[1] = sqrt(double(Shape0[1]));
  SqrtShape0[2] = sqrt(double(Shape0[2]));
  Shape6[0] = SqrtShape0[0]; Shape6[1] = Shape6[0]; Shape6[2] = SqrtShape0[1]; Shape6[3] = Shape6[2];
  Shape6[4] = SqrtShape0[2]; Shape6[5] = Shape6[4];
  Tensor0.reshape(Shape6);
  Tensor0.permute(Order6);
  Shape2[0] = sqrtDim; Shape2[1] = sqrtDim;
  Tensor0.reshape(Shape2);
  Matrix0 = Matrix<T>(sqrtDim, sqrtDim);
  for (int j = 0; j < sqrtDim; j++)
  {
   Index2[1] = j;
   for (int i = 0; i < sqrtDim; i++)
   {
    Index2[0] = i;
    Matrix0(i, j) = Tensor0.get(Index2);
   }
  }
  Matrix0.setType("hermitian");
  W = vector<double>(sqrtDim); Vr = Matrix<T>(sqrtDim, sqrtDim);
  Matrix0.eigenDecompose(W, Vr);
  AShape[0] = SqrtShape0[0]; AShape[1] = SqrtShape0[1]; AShape[2] = SqrtShape0[2];
  AShape[3] = sqrtDim;
  ATensor = Tensor<T>(AShape);
  for (int j = 0; j < sqrtDim; j++)
  {
   for (int i = 0; i < sqrtDim; i++)
   {
    ATensor.set(i+j*sqrtDim, Vr(i, sqrtDim-1-j)*sqrt(abs(W[sqrtDim-1-j])));
   }
  }
  Shape3[0] = SqrtShape0[0]; Shape3[1] = SqrtShape0[1]; Shape3[2] = SqrtShape0[2]*sqrtDim;
  ATensor.reshape(Shape3);
  MPS0.set(position, ATensor);
  MPS0.d[position] = SqrtShape0[2]*sqrtDim;
  d2[position] = sqrtDim;
 }
}

template<class T> void MPS<T>::getPurificationEigenvalues(vector< vector<double> >& PuriEvals) const
{
// define PuriEvals:
 PuriEvals = vector< vector<double> >(this->N);
// fill PuriEvals:
 unsigned int dim, sqrtDim;
 vector<unsigned int> Shape0(3), SqrtShape0(3), Shape6(6), Order6(6), Shape2(2), Index2(2);
 vector<double> W; Matrix<T> Vr;
 Tensor<T> Tensor0;
 Matrix<T> Matrix0;
 Order6[0] = 0; Order6[1] = 2; Order6[2] = 4; Order6[3] = 1; Order6[4] = 3; Order6[5] = 5;
 for (unsigned int position = 0; position < this->N; position++)
 {
  this->get(position, Tensor0);
  Tensor0.getShape(Shape0);
  dim = Shape0[0]*Shape0[1]*Shape0[2];
  sqrtDim = sqrt(double(dim));
  SqrtShape0[0] = sqrt(double(Shape0[0])); SqrtShape0[1] = sqrt(double(Shape0[1]));
  SqrtShape0[2] = sqrt(double(Shape0[2]));
  Shape6[0] = SqrtShape0[0]; Shape6[1] = Shape6[0]; Shape6[2] = SqrtShape0[1]; Shape6[3] = Shape6[2];
  Shape6[4] = SqrtShape0[2]; Shape6[5] = Shape6[4];
  Tensor0.reshape(Shape6);
  Tensor0.permute(Order6);
  Shape2[0] = sqrtDim; Shape2[1] = sqrtDim;
  Tensor0.reshape(Shape2);
  Matrix0 = Matrix<T>(sqrtDim, sqrtDim);
  for (int j = 0; j < sqrtDim; j++)
  {
   Index2[1] = j;
   for (int i = 0; i < sqrtDim; i++)
   {
    Index2[0] = i;
    Matrix0(i, j) = Tensor0.get(Index2);
   }
  }
  Matrix0.setType("hermitian");
  W = vector<double>(sqrtDim); Vr = Matrix<T>(sqrtDim, sqrtDim);
  Matrix0.eigenDecompose(W, Vr);
  PuriEvals[position] = W;
 }
}

template<class T> void MPS<T>::getMPSFromPurification(const vector<unsigned int>& d2, MPS<T>& MPS0) const
{
// define MPDO-MPS MPS0:
 unsigned int D2 = this->D*this->D;
 MPS0 = MPS<T>(this->BC, this->N, this->d, D2);
// fill MPS0:
 vector<unsigned int> Shape0(3), Shape4(4), Indices0(1), Indices1(1), Order6(6);
 Tensor<T> Tensor0, Tensor1;
 Indices0[0] = 3; Indices1[0] = 3;
 Order6[0] = 0; Order6[1] = 3; Order6[2] = 1; Order6[3] = 4; Order6[4] = 2; Order6[5] = 5;
 for (unsigned int position = 0; position < this->N; position++)
 {
  this->get(position, Tensor0);
  Tensor0.getShape(Shape0);
  Shape4[0] = Shape0[0]; Shape4[1] = Shape0[1]; Shape4[2] = Shape0[2]/d2[position]; Shape4[3] = d2[position];
  Tensor0.reshape(Shape4);
  Tensor1 = Tensor0; Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order6);
  Shape0[0] = Shape4[0]*Shape4[0]; Shape0[1] = Shape4[1]*Shape4[1]; Shape0[2] = Shape4[2]*Shape4[2];
  Tensor0.reshape(Shape0);
  MPS0.set(position, Tensor0);
  MPS0.d[position] = Shape0[2];
 }
}

template<class T> void MPS<T>::getMPDOTensor(const vector<unsigned int>& d2, unsigned int position, Tensor<T>& Tensor0) const
{
 vector<unsigned int> Shape0(3), Shape4(4), Indices0(1), Indices1(1), Order6(6);
 Tensor<T> Tensor1;
 this->get(position, Tensor0);
 Tensor0.getShape(Shape0);
 Shape4[0] = Shape0[0]; Shape4[1] = Shape0[1]; Shape4[2] = Shape0[2]/d2[position]; Shape4[3] = d2[position];
 Tensor0.reshape(Shape4);
 Tensor1 = Tensor0; Tensor1.complexConjugate();
 Indices0[0] = 3; Indices1[0] = 3;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Order6[0] = 0; Order6[1] = 3; Order6[2] = 1; Order6[3] = 4; Order6[4] = 2; Order6[5] = 5;
 Tensor0.permute(Order6);
 Shape0[0] = Shape4[0]*Shape4[0]; Shape0[1] = Shape4[1]*Shape4[1]; Shape0[2] = Shape4[2]*Shape4[2];
 Tensor0.reshape(Shape0);
}

template<class T> void MPS<T>::updateTensor(Tensor<T>& Tensor0, Tensor<T>& PuriTensor) const
{
 unsigned int dim, sqrtDim;
 vector<unsigned int> Shape0(3), SqrtShape0(3), Shape6(6), Order6(6), Shape2(2), Index2(2);
 vector<unsigned int> PuriShape(3);
 Tensor0.getShape(Shape0);
 dim = Shape0[0]*Shape0[1]*Shape0[2];
 sqrtDim = sqrt(double(dim));
 SqrtShape0[0] = sqrt(double(Shape0[0])); SqrtShape0[1] = sqrt(double(Shape0[1]));
 SqrtShape0[2] = sqrt(double(Shape0[2]));
 Shape6[0] = SqrtShape0[0]; Shape6[1] = Shape6[0]; Shape6[2] = SqrtShape0[1]; Shape6[3] = Shape6[2];
 Shape6[4] = SqrtShape0[2]; Shape6[5] = Shape6[4];
 Tensor0.reshape(Shape6);
 Order6[0] = 0; Order6[1] = 2; Order6[2] = 4; Order6[3] = 1; Order6[4] = 3; Order6[5] = 5;
 Tensor0.permute(Order6);
 Shape2[0] = sqrtDim; Shape2[1] = sqrtDim;
 Tensor0.reshape(Shape2);
 Matrix<T> IMatrix(sqrtDim, sqrtDim);
 IMatrix.fillZeroes();
 for (int i = 0; i < sqrtDim; i++)
  IMatrix(i, i) = 1.0;
 Matrix<T> Matrix0(dim, dim);
 int i, j;
 for (int j1 = 0; j1 < sqrtDim; j1++)
 {
  for (int j0 = 0; j0 < sqrtDim; j0++)
  {
   Index2[1] = j0;
   j = j0+j1*sqrtDim;
   for (int i1 = 0; i1 < sqrtDim; i1++)
   {
    for (int i0 = 0; i0 < sqrtDim; i0++)
    {
     Index2[0] = i0;
     i = i0+i1*sqrtDim;
     Matrix0(j, i) = Tensor0.get(Index2)*IMatrix(i1, j1);
    }
   }
  }
 }
 Matrix0.setType("hermitian");
 vector<double> W(dim); Matrix<T> Vr(dim, dim);
 Matrix0.eigenDecompose(W, Vr);
 PuriShape[0] = SqrtShape0[0]; PuriShape[1] = SqrtShape0[1]; PuriShape[2] = SqrtShape0[2]*sqrtDim;
 PuriTensor = Tensor<T>(PuriShape);
 for (i = 0; i < dim; i++)
  PuriTensor.set(i, Vr(i, dim-1));
}

template<class T> void multiplyMPDOMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                          unsigned int maxNumSweeps, double& errorAchieved,
                                          unsigned int& numSweepsDone, MPS<T>& MPS1)
{
#ifdef DEBUG
 if ((MPO0.BC != MPS0.BC) || (MPO0.BC != MPS1.BC) || (MPO0.N != MPS0.N) || (MPO0.N != MPS1.N) ||
     (MPO0.N == 0) || (MPO0.d != MPS0.d) || (MPO0.d != MPS1.d))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> void " <<
          "multiplyMPDOMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps, " <<
                             "unsigned int maxNumSweeps, double& errorAchieved, " <<
                             "unsigned int& numSweepsDone, MPS<T>& MPS1): " <<
          "MPO0, MPS0 and MPS1 are not of the same form or (MPO0.N == 0)." << endl;
  exit(1);
 }
#endif
// Variables for error computation:
 double normSquaredMPO0MPS0, realScalarProductMPS1MPO0MPS0, normSquaredMPS1;
 if (eps >= 0.0)
 {
// During the sweep we need normSquaredMPO0MPS0=<MPS0|MPO0^{+}MPO0|MPS0> for computing the error:
  MPO<T> MPO0Adjoint, MPO0Product;
  MPO0.adjoint(MPO0Adjoint);
  MPO0Adjoint.multiply(MPO0, MPO0Product);
  normSquaredMPO0MPS0 = MathAuxiliary::convertToDouble(expectationValueMPO(MPS0, MPO0Product));
 }
 if (MPS1.D >= MPO0.D*MPS0.D)
 {
  Tensor<T> Tensor0, Tensor1;
  vector<unsigned int> Indices0(1), Indices1(1), Order(5), ShapeMPO0(4), ShapeMPS0(3), Shape(3);
  Indices0[0] = 2; Indices1[0] = 2;
  Order[0] = 0; Order[1] = 2; Order[2] = 1; Order[3] = 3; Order[4] = 4;
  for (int i = 0; i < MPO0.N; i++)
  {
   Tensor0 = MPS0.Tensors[i]; Tensor1 = MPO0.Tensors[i];
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order);
   MPS0.Tensors[i].getShape(ShapeMPS0); MPO0.Tensors[i].getShape(ShapeMPO0);
   Shape[0] = ShapeMPS0[0]*ShapeMPO0[0];
   Shape[1] = ShapeMPS0[1]*ShapeMPO0[1];
   Shape[2] = ShapeMPO0[3];
   Tensor0.reshape(Shape);
   MPS1.Tensors[i] = Tensor0;
  }
// in the case of open boundary conditions we have to recover the normal shape:
  if (MPO0.BC == "open")
   MPS1.bringIntoNormalShape();
// we have to reset D:
  MPS1.setD(MPS1.D);
  if (eps >= 0.0)
   errorAchieved = distanceMPOMPS(MPO0, MPS0, MPS1)/sqrt(abs(normSquaredMPO0MPS0));
  numSweepsDone = 0;
 }
 else
 {
// purification Puri1 of MPS1:
  MPS<T> Puri1; vector<unsigned int> d2;
  MPS1.getPurificationFromMPS(Puri1, d2);
  double error0, error1;
  if (eps >= 0.0)
   error1 = 0.0;
  numSweepsDone = 0;
  if (MPO0.BC == "open")
  {
// First we bring Puri1 into normal form from right to left:
   Puri1.bringIntoNormalForm();
// and obtain the corresponding MPS1:
   Puri1.getMPSFromPurification(d2, MPS1);
// Then we construct the b for each site:
   Tensor<T> Tensorsb[MPO0.N];
   Tensor<T> Tensor0, Tensor1, Tensor2, PuriTensor;
   vector<unsigned int> Indices0(1), Indices1(1), Order(6), Indices02(2), Indices12(2);
   Tensor0 = MPS0.Tensors[MPO0.N-1]; Tensor1 = MPO0.Tensors[MPO0.N-1];
   Indices0[0] = 2; Indices1[0] = 2;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   MPS1.Tensors[MPO0.N-1].complexConjugate(Tensor1);
   Indices0[0] = 4; Indices1[0] = 2;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Order[0] = 1; Order[1] = 3; Order[2] = 5; Order[3] = 0; Order[4] = 2; Order[5] = 4;
   Tensor0.permute(Order);
   Tensorsb[MPO0.N-1] = Tensor0;
   for (int i = MPO0.N-2; i > 0; i--)
   {
    Tensor1 = MPS0.Tensors[i];
    Indices0[0] = 3; Indices1[0] = 1;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = MPO0.Tensors[i];
    Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    MPS1.Tensors[i].complexConjugate(Tensor1);
    Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    Tensorsb[i] = Tensor0;
   }
   string Direction = "right";
   vector<unsigned int> Indices03(3), Index04(4), Index06(6), Order0(3);
   Matrix<T> X;
   while (numSweepsDone < maxNumSweeps)
   {
    if (Direction == "right")
    {
     Tensor0 = Tensorsb[1]; Tensor1 = MPS0.Tensors[0];
     Indices0[0] = 3; Indices1[0] = 1;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Tensor1 = MPO0.Tensors[0];
     Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     if (eps >= 0.0)
     {
// we compute the error:
      Tensor1 = Tensor0;
      MPS1.Tensors[0].complexConjugate(Tensor2);
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
      Tensor1.contract(Indices02, Tensor2, Indices12);
      Index06[0] = 0; Index06[1] = 0; Index06[2] = 0; Index06[3] = 0; Index06[4] = 0; Index06[5] = 0;
      realScalarProductMPS1MPO0MPS0 = MathAuxiliary::convertToDouble(Tensor1.get(Index06));
     }
// we prepare the sweep to the right:
     Indices03[0] = 2; Indices03[1] = 3; Indices03[2] = 6;
     Index04[0] = 0; Index04[1] = 0; Index04[2] = 0; Index04[3] = 0;
     Tensor0.getSubtensor(Indices03, Index04, Tensor1);
     MPS1.updateTensor(Tensor1, PuriTensor);
     Puri1.Tensors[0] = PuriTensor;
     Puri1.normalizeTensor(0, Direction, X);
     Puri1.getMPDOTensor(d2, 0, Tensor1);
     MPS1.Tensors[0] = Tensor1;
     Tensor0 = MPS0.Tensors[0]; Tensor1 = MPO0.Tensors[0];
     Indices0[0] = 2; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     MPS1.Tensors[0].complexConjugate(Tensor1);
     Indices0[0] = 4; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Order[0] = 0; Order[1] = 2; Order[2] = 4; Order[3] = 1; Order[4] = 3; Order[5] = 5;
     Tensor0.permute(Order);
     Tensorsb[0] = Tensor0;
     for (int i = 1; i < MPO0.N-1; i++)
     {
      Tensor1 = MPS0.Tensors[i];
      Indices0[0] = 3; Indices1[0] = 0;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      Tensor1 = MPO0.Tensors[i];
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      Tensor2 = Tensor0;
      Tensor1 = Tensorsb[i+1];
      Indices02[0] = 4; Indices02[1] = 5; Indices12[0] = 3; Indices12[1] = 4;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      Indices03[0] = 3; Indices03[1] = 8; Indices03[2] = 4;
      Index06[0] = 0; Index06[1] = 0; Index06[2] = 0; Index06[3] = 0; Index06[4] = 0; Index06[5] = 0;
      Tensor0.getSubtensor(Indices03, Index06, Tensor1);
      MPS1.updateTensor(Tensor1, PuriTensor);
      Puri1.Tensors[i] = PuriTensor;
      Puri1.normalizeTensor(i, Direction, X);
      Puri1.getMPDOTensor(d2, i, Tensor1);
      MPS1.Tensors[i] = Tensor1;
      Tensor0 = Tensor2;
      MPS1.Tensors[i].complexConjugate(Tensor1);
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      Tensorsb[i] = Tensor0;
     }
     Indices0[0] = 0; Indices1[0] = 1;
     Puri1.Tensors[MPO0.N-1].contract(Indices0, X, Indices1);
     Order0[0] = 2; Order0[1] = 0; Order0[2] = 1;
     Puri1.Tensors[MPO0.N-1].permute(Order0);
     Puri1.getMPDOTensor(d2, MPO0.N-1, Tensor1);
     MPS1.Tensors[MPO0.N-1] = Tensor1;
     Direction = "left";
    }
    else if (Direction == "left")
    {
     Tensor0 = Tensorsb[MPO0.N-2]; Tensor1 = MPS0.Tensors[MPO0.N-1];
     Indices0[0] = 3; Indices1[0] = 0;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Tensor1 = MPO0.Tensors[MPO0.N-1];
     Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     if (eps >= 0.0)
     {
// we compute the error:
      Tensor1 = Tensor0;
      MPS1.Tensors[MPO0.N-1].complexConjugate(Tensor2);
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
      Tensor1.contract(Indices02, Tensor2, Indices12);
      Index06[0] = 0; Index06[1] = 0; Index06[2] = 0; Index06[3] = 0; Index06[4] = 0; Index06[5] = 0;
      realScalarProductMPS1MPO0MPS0 = MathAuxiliary::convertToDouble(Tensor1.get(Index06));
     }
// we prepare the sweep to the left:
     Indices03[0] = 3; Indices03[1] = 2; Indices03[2] = 6;
     Index04[0] = 0; Index04[1] = 0; Index04[2] = 0; Index04[3] = 0;
     Tensor0.getSubtensor(Indices03, Index04, Tensor1);
     MPS1.updateTensor(Tensor1, PuriTensor);
     Puri1.Tensors[MPO0.N-1] = PuriTensor;
     Puri1.normalizeTensor(MPO0.N-1, Direction, X);
     Puri1.getMPDOTensor(d2, MPO0.N-1, Tensor1);
     MPS1.Tensors[MPO0.N-1] = Tensor1;
     Tensor0 = MPS0.Tensors[MPO0.N-1]; Tensor1 = MPO0.Tensors[MPO0.N-1];
     Indices0[0] = 2; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     MPS1.Tensors[MPO0.N-1].complexConjugate(Tensor1);
     Indices0[0] = 4; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Order[0] = 1; Order[1] = 3; Order[2] = 5; Order[3] = 0; Order[4] = 2; Order[5] = 4;
     Tensor0.permute(Order);
     Tensorsb[MPO0.N-1] = Tensor0;
     for (int i = MPO0.N-2; i > 0; i--)
     {
      Tensor1 = MPS0.Tensors[i];
      Indices0[0] = 3; Indices1[0] = 1;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      Tensor1 = MPO0.Tensors[i];
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      Tensor2 = Tensor0;
      Tensor1 = Tensorsb[i-1];
      Indices02[0] = 4; Indices02[1] = 5; Indices12[0] = 3; Indices12[1] = 4;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      Indices03[0] = 8; Indices03[1] = 3; Indices03[2] = 4;
      Index06[0] = 0; Index06[1] = 0; Index06[2] = 0; Index06[3] = 0; Index06[4] = 0; Index06[5] = 0;
      Tensor0.getSubtensor(Indices03, Index06, Tensor1);
      MPS1.updateTensor(Tensor1, PuriTensor);
      Puri1.Tensors[i] = PuriTensor;
      Puri1.normalizeTensor(i, Direction, X);
      Puri1.getMPDOTensor(d2, i, Tensor1);
      MPS1.Tensors[i] = Tensor1;
      Tensor0 = Tensor2;
      MPS1.Tensors[i].complexConjugate(Tensor1);
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      Tensorsb[i] = Tensor0;
     }
     Indices0[0] = 1; Indices1[0] = 0;
     Puri1.Tensors[0].contract(Indices0, X, Indices1);
     Order0[0] = 0; Order0[1] = 2; Order0[2] = 1;
     Puri1.Tensors[0].permute(Order0);
     Puri1.getMPDOTensor(d2, 0, Tensor1);
     MPS1.Tensors[0] = Tensor1;
     Direction = "right";
    }
    numSweepsDone++;
   }
   if (eps >= 0.0)
    errorAchieved = distanceMPOMPS(MPO0, MPS0, MPS1)/sqrt(abs(normSquaredMPO0MPS0));
  }
  else if (MPO0.BC == "periodic")
  {
   cout << "The following friend function is not implemented for periodic boundary conditions: " <<
           "template<class T> void " <<
           "multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps, " <<
                          "unsigned int maxNumSweeps, double& errorAchieved, " <<
                          "unsigned int& numSweepsDone, MPS<T>& MPS1)." << endl;
  }
 }
}

template<class T> void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                      unsigned int maxNumSweeps, double& errorAchieved,
                                      unsigned int& numSweepsDone, MPS<T>& MPS1,
                                      double cutoff, unsigned int mode, unsigned int d2,
                                      double alphaStart, double x, double precision, unsigned int maxNumIter)
{
#ifdef DEBUG
 if ((MPO0.BC != MPS0.BC) || (MPO0.BC != MPS1.BC) || (MPO0.N != MPS0.N) || (MPO0.N != MPS1.N) ||
     (MPO0.N == 0) || (MPO0.d != MPS0.d) || (MPO0.d != MPS1.d))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> void " <<
          "multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps, " <<
                         "unsigned int maxNumSweeps, double& errorAchieved, " <<
                         "unsigned int& numSweepsDone, MPS<T>& MPS1, " <<
                         "double cutoff, unsigned int mode, unsigned int d2, " <<
                         "double alphaStart, double x, double precision, unsigned int maxNumIter): " <<
          "MPO0, MPS0 and MPS1 are not of the same form or (MPO0.N == 0)." << endl;
  exit(1);
 }
#endif
 if (MPS1.D >= MPO0.D*MPS0.D)
 {
// exact MPO0-MPS0 product:
  unsigned int D1 = MPS1.D;
  exactMultiplyMPOMPS(MPO0, MPS0, MPS1);
// reset D:
  MPS1.setD(D1);
  errorAchieved = 0.0;
  numSweepsDone = 0;
 }
 else
 {
// approximate MPO0-MPS0 product:
// variables for error computation:
  double normSquaredMPO0MPS0, realScalarProductMPS1MPO0MPS0, normSquaredMPS1;
  if (eps >= 0.0)
  {
// normSquaredMPO0MPS0=<MPS0|MPO0^{+}MPO0|MPS0>:
   MPS<T> MPS2;
   exactMultiplyMPOMPS(MPO0, MPS0, MPS2);
   normSquaredMPO0MPS0 = MPS2.scalarProduct();
  }
  double error0, error1;
  if (eps >= 0.0)
   error1 = 0.0;
  numSweepsDone = 0;
  if (MPO0.BC == "open")
  {
   vector<unsigned int> Indices0(1), Indices1(1), Order4(4), Order6(6), Indices02(2), Indices12(2);
   Tensor<T> TensorsN[MPO0.N], Tensorsb[MPO0.N];
   Tensor<T> Tensor0, Tensor1, Tensor2;
// TensorsN from right to left:
   Tensor0 = MPS1.Tensors[MPO0.N-1]; MPS1.Tensors[MPO0.N-1].complexConjugate(Tensor1);
   Indices0[0] = 2; Indices1[0] = 2;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Order4[0] = 1; Order4[1] = 3; Order4[2] = 0; Order4[3] = 2;
   Tensor0.permute(Order4);
   TensorsN[MPO0.N-1] = Tensor0;
   Indices0[0] = 2; Indices1[0] = 1;
   Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
   for (int pos = MPO0.N-2; pos > 0; pos--)
   {
    Tensor1 = MPS1.Tensors[pos];
    Tensor0.contract(Indices0, Tensor1, Indices1);
    MPS1.Tensors[pos].complexConjugate(Tensor1);
    Tensor0.contract(Indices02, Tensor1, Indices12);
    TensorsN[pos] = Tensor0;
   }
// Tensorsb from right to left:
   Tensor0 = MPS0.Tensors[MPO0.N-1]; Tensor1 = MPO0.Tensors[MPO0.N-1];
   Indices0[0] = 2; Indices1[0] = 2;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   MPS1.Tensors[MPO0.N-1].complexConjugate(Tensor1);
   Indices0[0] = 4; Indices1[0] = 2;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Order6[0] = 1; Order6[1] = 3; Order6[2] = 5; Order6[3] = 0; Order6[4] = 2; Order6[5] = 4;
   Tensor0.permute(Order6);
   Tensorsb[MPO0.N-1] = Tensor0;
   Indices0[0] = 3; Indices1[0] = 1;
   Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
   for (int pos = MPO0.N-2; pos > 0; pos--)
   {
    Tensor1 = MPS0.Tensors[pos];
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = MPO0.Tensors[pos];
    Tensor0.contract(Indices02, Tensor1, Indices12);
    MPS1.Tensors[pos].complexConjugate(Tensor1);
    Tensor0.contract(Indices02, Tensor1, Indices12);
    Tensorsb[pos] = Tensor0;
   }
   string Direction = "right";
   vector<unsigned int> Index4(4), Index6(6);
   Index4[0] = 0; Index4[1] = 0; Index4[2] = 0; Index4[3] = 0;
   Index6[0] = 0; Index6[1] = 0; Index6[2] = 0; Index6[3] = 0; Index6[4] = 0; Index6[5] = 0;
   while (numSweepsDone < maxNumSweeps)
   {
    if (Direction == "right")
    {
// b for positon 0:
     Tensor0 = Tensorsb[1]; Tensor1 = MPS0.Tensors[0];
     Indices0[0] = 3; Indices1[0] = 1;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Tensor1 = MPO0.Tensors[0];
     Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     if (eps >= 0.0)
     {
// error computation:
      Tensor1 = Tensor0; MPS1.Tensors[0].complexConjugate(Tensor2);
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
      Tensor1.contract(Indices02, Tensor2, Indices12);
      realScalarProductMPS1MPO0MPS0 = MathAuxiliary::convertToDouble(Tensor1.get(Index6));
      Tensor1 = TensorsN[1]; Tensor2 = MPS1.Tensors[0];
      Indices0[0] = 2; Indices1[0] = 1;
      Tensor1.contract(Indices0, Tensor2, Indices1);
      MPS1.Tensors[0].complexConjugate(Tensor2);
      Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
      Tensor1.contract(Indices02, Tensor2, Indices12);
      normSquaredMPS1 = MathAuxiliary::convertToDouble(Tensor1.get(Index4));
      error0 = sqrt(abs(normSquaredMPO0MPS0 - 2.0*realScalarProductMPS1MPO0MPS0 + normSquaredMPS1));
      if (abs(error1-error0)/sqrt(normSquaredMPO0MPS0) < eps)
      {
       errorAchieved = error0/sqrt(normSquaredMPO0MPS0);
       return;
      }
     }
// update tensor on position 0:
     MPS1.updateTensor(Direction, TensorsN[1], Tensor0, cutoff, mode, d2, alphaStart, x, precision, maxNumIter);
// TensorsN[0]:
     Tensor0 = MPS1.Tensors[0]; MPS1.Tensors[0].complexConjugate(Tensor1);
     Indices0[0] = 2; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Order4[0] = 0; Order4[1] = 2; Order4[2] = 1; Order4[3] = 3;
     Tensor0.permute(Order4);
     TensorsN[0] = Tensor0;
// Tensorsb[0]:
     Tensor0 = MPS0.Tensors[0]; Tensor1 = MPO0.Tensors[0];
     Indices0[0] = 2; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     MPS1.Tensors[0].complexConjugate(Tensor1);
     Indices0[0] = 4; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Order6[0] = 0; Order6[1] = 2; Order6[2] = 4; Order6[3] = 1; Order6[4] = 3; Order6[5] = 5;
     Tensor0.permute(Order6);
     Tensorsb[0] = Tensor0;
     for (int pos = 1; pos < MPO0.N-1; pos++)
     {
// b for position pos:
      Tensor0 = Tensorsb[pos-1]; Tensor1 = MPS0.Tensors[pos];
      Indices0[0] = 3; Indices1[0] = 0;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      Tensor1 = MPO0.Tensors[pos];
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      Tensor2 = Tensor0;
      Tensor1 = Tensorsb[pos+1];
      Indices02[0] = 4; Indices02[1] = 5; Indices12[0] = 3; Indices12[1] = 4;
      Tensor0.contract(Indices02, Tensor1, Indices12);
// update tensor on position pos:
      MPS1.updateTensor(pos, Direction, TensorsN[pos-1], TensorsN[pos+1], Tensor0, cutoff, mode, d2, alphaStart, x, precision, maxNumIter);
// Tensorsb[pos]:
      MPS1.Tensors[pos].complexConjugate(Tensor1);
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
      Tensor2.contract(Indices02, Tensor1, Indices12);
      Tensorsb[pos] = Tensor2;
// TensorsN[pos]:
      Tensor0 = TensorsN[pos-1]; Tensor1 = MPS1.Tensors[pos];
      Indices0[0] = 2; Indices1[0] = 0;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      MPS1.Tensors[pos].complexConjugate(Tensor1);
      Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      TensorsN[pos] = Tensor0;
     }
     Direction = "left";
    }
    else if (Direction == "left")
    {
// b for positon N-1:
     Tensor0 = Tensorsb[MPO0.N-2]; Tensor1 = MPS0.Tensors[MPO0.N-1];
     Indices0[0] = 3; Indices1[0] = 0;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Tensor1 = MPO0.Tensors[MPO0.N-1];
     Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     if (eps >= 0.0)
     {
// error computation:
      Tensor1 = Tensor0; MPS1.Tensors[MPO0.N-1].complexConjugate(Tensor2);
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
      Tensor1.contract(Indices02, Tensor2, Indices12);
      realScalarProductMPS1MPO0MPS0 = MathAuxiliary::convertToDouble(Tensor1.get(Index6));
      Tensor1 = TensorsN[MPO0.N-2]; Tensor2 = MPS1.Tensors[MPO0.N-1];
      Indices0[0] = 2; Indices1[0] = 0;
      Tensor1.contract(Indices0, Tensor2, Indices1);
      MPS1.Tensors[MPO0.N-1].complexConjugate(Tensor2);
      Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 2;
      Tensor1.contract(Indices02, Tensor2, Indices12);
      normSquaredMPS1 = MathAuxiliary::convertToDouble(Tensor1.get(Index4));
      error1 = sqrt(abs(normSquaredMPO0MPS0 - 2.0*realScalarProductMPS1MPO0MPS0 + normSquaredMPS1));
      if (abs(error0-error1)/sqrt(normSquaredMPO0MPS0) < eps)
      {
       errorAchieved = error1/sqrt(normSquaredMPO0MPS0);
       return;
      }
     }
// update tensor on position N-1:
     MPS1.updateTensor(Direction, TensorsN[MPO0.N-2], Tensor0, cutoff, mode, d2, alphaStart, x, precision, maxNumIter);
// TensorsN[N-1]:
     Tensor0 = MPS1.Tensors[MPO0.N-1]; MPS1.Tensors[MPO0.N-1].complexConjugate(Tensor1);
     Indices0[0] = 2; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Order4[0] = 1; Order4[1] = 3; Order4[2] = 0; Order4[3] = 2;
     Tensor0.permute(Order4);
     TensorsN[MPO0.N-1] = Tensor0;
// Tensorsb[N-1]:
     Tensor0 = MPS0.Tensors[MPO0.N-1]; Tensor1 = MPO0.Tensors[MPO0.N-1];
     Indices0[0] = 2; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     MPS1.Tensors[MPO0.N-1].complexConjugate(Tensor1);
     Indices0[0] = 4; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Order6[0] = 1; Order6[1] = 3; Order6[2] = 5; Order6[3] = 0; Order6[4] = 2; Order6[5] = 4;
     Tensor0.permute(Order6);
     Tensorsb[MPO0.N-1] = Tensor0;
     for (int pos = MPO0.N-2; pos > 0; pos--)
     {
// b for position pos:
      Tensor0 = Tensorsb[pos+1]; Tensor1 = MPS0.Tensors[pos];
      Indices0[0] = 3; Indices1[0] = 1;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      Tensor1 = MPO0.Tensors[pos];
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      Tensor2 = Tensor0;
      Tensor1 = Tensorsb[pos-1];
      Indices02[0] = 4; Indices02[1] = 5; Indices12[0] = 3; Indices12[1] = 4;
      Tensor0.contract(Indices02, Tensor1, Indices12);
// update tensor on position pos:
      MPS1.updateTensor(pos, Direction, TensorsN[pos-1], TensorsN[pos+1], Tensor0, cutoff, mode, d2, alphaStart, x, precision, maxNumIter);
// Tensorsb[pos]:
      MPS1.Tensors[pos].complexConjugate(Tensor1);
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
      Tensor2.contract(Indices02, Tensor1, Indices12);
      Tensorsb[pos] = Tensor2;
// TensorsN[pos]:
      Tensor0 = TensorsN[pos+1]; Tensor1 = MPS1.Tensors[pos];
      Indices0[0] = 2; Indices1[0] = 1;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      MPS1.Tensors[pos].complexConjugate(Tensor1);
      Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      TensorsN[pos] = Tensor0;
     }
     Direction = "right";
    }
    numSweepsDone++;
   }
   if (eps >= 0.0)
   {
// final error computation:
    if (Direction == "right")
    {
     Tensor0 = Tensorsb[1]; Tensor1 = MPS0.Tensors[0];
     Indices0[0] = 3; Indices1[0] = 1;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Tensor1 = MPO0.Tensors[0];
     Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     MPS1.Tensors[0].complexConjugate(Tensor1);
     Tensor0.contract(Indices02, Tensor1, Indices12);
     realScalarProductMPS1MPO0MPS0 = MathAuxiliary::convertToDouble(Tensor0.get(Index6));
     Tensor0 = TensorsN[1]; Tensor1 = MPS1.Tensors[0];
     Indices0[0] = 2; Indices1[0] = 1;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     MPS1.Tensors[0].complexConjugate(Tensor1);
     Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     normSquaredMPS1 = MathAuxiliary::convertToDouble(Tensor0.get(Index4));
     errorAchieved = sqrt(abs(normSquaredMPO0MPS0-2.0*realScalarProductMPS1MPO0MPS0+normSquaredMPS1))/
                     sqrt(normSquaredMPO0MPS0);
    }
    else if (Direction == "left")
    {
     Tensor0 = Tensorsb[MPO0.N-2]; Tensor1 = MPS0.Tensors[MPO0.N-1];
     Indices0[0] = 3; Indices1[0] = 0;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Tensor1 = MPO0.Tensors[MPO0.N-1];
     Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     MPS1.Tensors[MPO0.N-1].complexConjugate(Tensor1);
     Tensor0.contract(Indices02, Tensor1, Indices12);
     realScalarProductMPS1MPO0MPS0 = MathAuxiliary::convertToDouble(Tensor0.get(Index6));
     Tensor0 = TensorsN[MPO0.N-2]; Tensor1 = MPS1.Tensors[MPO0.N-1];
     Indices0[0] = 2; Indices1[0] = 0;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     MPS1.Tensors[MPO0.N-1].complexConjugate(Tensor1);
     Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     normSquaredMPS1 = MathAuxiliary::convertToDouble(Tensor0.get(Index4));
     errorAchieved = sqrt(abs(normSquaredMPO0MPS0-2.0*realScalarProductMPS1MPO0MPS0+normSquaredMPS1))/
                     sqrt(normSquaredMPO0MPS0);
    }
   }
  }
  else if (MPO0.BC == "periodic")
  {
   cout << "The following friend function is not implemented for periodic boundary conditions: " <<
           "template<class T> void " <<
           "multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps, " <<
                          "unsigned int maxNumSweeps, double& errorAchieved, " <<
                          "unsigned int& numSweepsDone, MPS<T>& MPS1, " <<
                          "double cutoff, unsigned int mode, unsigned int d2, " <<
                          "double alphaStart, double x, double precision, unsigned int maxNumIter)." << endl;
  }
 }
}

template<class T> void MPS<T>::updateTensor(const string& Direction,
                                            const Tensor<T>& NormTensor, const Tensor<T>& MTensor,
                                            double cutoff, unsigned int mode, unsigned int d2,
                                            double alphaStart, double x, double precision, unsigned int maxNumIter)
{
 Tensor<T> NTensor(NormTensor), bTensor(MTensor);
 vector<unsigned int> Shape0;
 if (Direction == "right")
  this->Tensors[0].getShape(Shape0);
 else if (Direction == "left")
  this->Tensors[this->N-1].getShape(Shape0);
#ifdef DEBUG
 if ((this->BC != "open") || (this->N == 0) || ((Direction != "right") && (Direction != "left")) ||
     ((mode != 0) && (mode != 1) && (mode != 2)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "updateTensor(const string& Direction, const Tensor<T>& NormTensor, const Tensor<T>& MTensor, " <<
                       "double cutoff, unsigned int mode, unsigned int d2, " <<
                       "double alphaStart, double x, double precision, unsigned int maxNumIter): " <<
          "((this->BC != open) || (this->N == 0) || ((Direction != right) && (Direction != left)) || " <<
           "((mode != 0) && (mode != 1) && (mode != 2)))." << endl;
  exit(1);
 }
 vector<unsigned int> NShape, bShape;
 NTensor.getShape(NShape); bTensor.getShape(bShape);
 if ((NTensor.getRank() != 4) || (NShape[0] != 1) || (NShape[1] != 1) ||
     ((NShape[2] != Shape0[1]) && (Direction == "right")) || ((NShape[2] != Shape0[0]) && (Direction == "left")) ||
     ((NShape[3] != Shape0[1]) && (Direction == "right")) || ((NShape[3] != Shape0[0]) && (Direction == "left")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "updateTensor(const string& Direction, const Tensor<T>& NormTensor, const Tensor<T>& MTensor, " <<
                       "double cutoff, unsigned int mode, unsigned int d2, " <<
                       "double alphaStart, double x, double precision, unsigned int maxNumIter): " <<
          "NormTensor has wrong Shape." << endl;
  exit(1);
 }
 if ((bTensor.getRank() != 7) || (bShape[0] != 1) || (bShape[1] != 1) || (bShape[2] != 1) ||
     ((bShape[3] != Shape0[1]) && (Direction == "right")) || ((bShape[3] != Shape0[0]) && (Direction == "left")) ||
     (bShape[4] != 1) || (bShape[5] != 1) || (bShape[6] != Shape0[2]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "updateTensor(const string& Direction, const Tensor<T>& NormTensor, const Tensor<T>& MTensor, " <<
                       "double cutoff, unsigned int mode, unsigned int d2, " <<
                       "double alphaStart, double x, double precision, unsigned int maxNumIter): " <<
          "MTensor has wrong Shape." << endl;
  exit(1);
 }
#endif
 if (mode == 2)
 {
  double alpha, distance0, distance1;
  unsigned int dim = Shape0[0]*Shape0[1]*Shape0[2], sqrtDim = sqrt(double(dim)), i, j, k;
  vector<unsigned int> Indices1(1), Indices12(2), Indices13(3), Indices2(1), Indices22(2), Indices23(3);
  vector<unsigned int> Index4(4), Index6(6), Shape2(2), Index2(2), Shape4(4), Shape6(6), Order6(6);
  vector<double> W;
  Tensor<T> Tensor0, Tensor1, NTensor1, NTensor2, bTensor1, bTensor2;
  Matrix<T> NMatrix, NMatrixH, NMatrixInv, bMatrix;
  Matrix<T> Matrix1(sqrtDim, sqrtDim), Matrix1H, Vr;
// proper d2:
  d2 = min(d2, sqrtDim);
// SqrtShape0 = sqrt(Shape0):
  vector<unsigned int> SqrtShape0(3);
  SqrtShape0[0] = sqrt(double(Shape0[0])); SqrtShape0[1] = sqrt(double(Shape0[1]));
  SqrtShape0[2] = sqrt(double(Shape0[2]));
// define IMatrix:=1 on physical bond:
  Matrix<T> IMatrix(SqrtShape0[2], SqrtShape0[2]);
  IMatrix.fillZeroes();
  for (i = 0; i < SqrtShape0[2]; i++)
   IMatrix(i, i) = 1.0;
// declare purification-Tensor=:ATensor:
  vector<unsigned int> AShape(4);
  AShape[0] = SqrtShape0[0]; AShape[1] = SqrtShape0[1]; AShape[2] = SqrtShape0[2]; AShape[3] = d2;
  Tensor<T> ATensor(AShape), BTensor;
// reshape NTensor:
  if (Direction == "right"){
   Shape4[0] = SqrtShape0[1]; Shape4[1] = Shape4[0]; Shape4[2] = SqrtShape0[1]; Shape4[3] = Shape4[2];
  }
  else if (Direction == "left"){
   Shape4[0] = SqrtShape0[0]; Shape4[1] = Shape4[0]; Shape4[2] = SqrtShape0[0]; Shape4[3] = Shape4[2];
  }
  NTensor.reshape(Shape4);
// reshape bTensor:
  if (Direction == "right"){
   Shape4[0] = SqrtShape0[1]; Shape4[1] = Shape4[0]; Shape4[2] = SqrtShape0[2]; Shape4[3] = Shape4[2];
  }
  else if (Direction == "left"){
   Shape4[0] = SqrtShape0[0]; Shape4[1] = Shape4[0]; Shape4[2] = SqrtShape0[2]; Shape4[3] = Shape4[2];
  }
  bTensor.reshape(Shape4);
// initialize iteration and obtain ATensor:
  if (Direction == "right")
   Tensor0 = this->Tensors[0];
  else if (Direction == "left")
   Tensor0 = this->Tensors[this->N-1];
  Tensor1 = Tensor0;
  Shape6[0] = SqrtShape0[0]; Shape6[1] = Shape6[0]; Shape6[2] = SqrtShape0[1]; Shape6[3] = Shape6[2];
  Shape6[4] = SqrtShape0[2]; Shape6[5] = Shape6[4];
  Tensor1.reshape(Shape6);
  Order6[0] = 0; Order6[1] = 2; Order6[2] = 4; Order6[3] = 1; Order6[4] = 3; Order6[5] = 5;
  Tensor1.permute(Order6);
  Shape2[0] = sqrtDim; Shape2[1] = sqrtDim;
  Tensor1.reshape(Shape2);
  for (j = 0; j < sqrtDim; j++){
   Index2[1] = j;
   for (i = 0; i < sqrtDim; i++){
    Index2[0] = i;
    Matrix1(i, j) = Tensor1.get(Index2);
   }
  }
  Matrix1H = Matrix1;
  Matrix1H.adjoint();
  Matrix1.add(Matrix1H);
  Matrix1.multiply(0.5);
  Matrix1.setType("hermitian");
  W = vector<double>(sqrtDim); Vr = Matrix<T>(sqrtDim, sqrtDim);
  Matrix1.eigenDecompose(W, Vr);
  for (i = 0; i < sqrtDim; i++)
  {
   if (W[i]/abs(W[sqrtDim-1]) < cutoff)
    W[i] = 0.0;
  }
  for (k = sqrtDim-d2; k < sqrtDim; k++){
   for (i = 0; i < sqrtDim; i++){
    ATensor.set(i+(k-sqrtDim+d2)*sqrtDim, Vr(i, k)*sqrt(W[k]));
   }
  }
  Tensor<T> ATensor0(ATensor);
  T NVal, bVal;
  distance0 = 0.0; distance1 = 1.0e10;
  alpha = alphaStart; unsigned int iter = 0;
// iterate:
  while ((abs(distance0-distance1)/abs(distance1) > precision) && (iter < maxNumIter))
  {
// construct NMatrix:
   NMatrix = Matrix<T>(sqrtDim*d2, sqrtDim*d2);
   NTensor1 = NTensor;
   if (Direction == "right"){
// - contract NTensor:
    BTensor = ATensor;
    BTensor.complexConjugate();
    Indices1[0] = 1; Indices2[0] = 1;
    NTensor1.contract(Indices1, BTensor, Indices2);
    BTensor = ATensor;
    Indices12[0] = 2; Indices12[1] = 4; Indices22[0] = 1; Indices22[1] = 2;
    NTensor1.contract(Indices12, BTensor, Indices22);
    NTensor2 = NTensor1;
    BTensor = ATensor;
    Indices12[0] = 0; Indices12[1] = 3; Indices22[0] = 1; Indices22[1] = 3;
    NTensor2.contract(Indices12, BTensor, Indices22);
    BTensor = ATensor;
    BTensor.complexConjugate();
    Indices13[0] = 0; Indices13[1] = 3; Indices13[2] = 5; Indices23[0] = 1; Indices23[1] = 3; Indices23[2] = 2;
    NTensor2.contract(Indices13, BTensor, Indices23);
    NVal = NTensor2.get(0);
// - build NMatrix:
    for (int j3 = 0; j3 < d2; j3++){
     Index6[3] = j3;
    for (int j2 = 0; j2 < SqrtShape0[2]; j2++){
    for (int j1 = 0; j1 < SqrtShape0[1]; j1++){
     Index6[0] = j1;
    for (int j0 = 0; j0 < SqrtShape0[0]; j0++){
     Index6[2] = j0;
     j = j0 + j1*SqrtShape0[0] + j2*SqrtShape0[0]*SqrtShape0[1] + j3*SqrtShape0[0]*SqrtShape0[1]*SqrtShape0[2];
     for (int i3 = 0; i3 < d2; i3++){
      Index6[5] = i3;
     for (int i2 = 0; i2 < SqrtShape0[2]; i2++){
     for (int i1 = 0; i1 < SqrtShape0[1]; i1++){
      Index6[1] = i1;
     for (int i0 = 0; i0 < SqrtShape0[0]; i0++){
      Index6[4] = i0;
      i = i0 + i1*SqrtShape0[0] + i2*SqrtShape0[0]*SqrtShape0[1] + i3*SqrtShape0[0]*SqrtShape0[1]*SqrtShape0[2];
      NMatrix(i, j) = NTensor1.get(Index6)*IMatrix(i2, j2);
     }}}}
    }}}}
   }
   else if (Direction == "left"){
// - contract NTensor:
    BTensor = ATensor;
    BTensor.complexConjugate();
    Indices1[0] = 1; Indices2[0] = 0;
    NTensor1.contract(Indices1, BTensor, Indices2);
    BTensor = ATensor;
    Indices12[0] = 2; Indices12[1] = 4; Indices22[0] = 0; Indices22[1] = 2;
    NTensor1.contract(Indices12, BTensor, Indices22);
    NTensor2 = NTensor1;
    BTensor = ATensor;
    Indices12[0] = 0; Indices12[1] = 3; Indices22[0] = 0; Indices22[1] = 3;
    NTensor2.contract(Indices12, BTensor, Indices22);
    BTensor = ATensor;
    BTensor.complexConjugate();
    Indices13[0] = 0; Indices13[1] = 3; Indices13[2] = 5; Indices23[0] = 0; Indices23[1] = 3; Indices23[2] = 2;
    NTensor2.contract(Indices13, BTensor, Indices23);
    NVal = NTensor2.get(0);
// - build NMatrix:
    for (int j3 = 0; j3 < d2; j3++){
     Index6[3] = j3;
    for (int j2 = 0; j2 < SqrtShape0[2]; j2++){
    for (int j1 = 0; j1 < SqrtShape0[1]; j1++){
     Index6[2] = j1;
    for (int j0 = 0; j0 < SqrtShape0[0]; j0++){
     Index6[0] = j0;
     j = j0 + j1*SqrtShape0[0] + j2*SqrtShape0[0]*SqrtShape0[1] + j3*SqrtShape0[0]*SqrtShape0[1]*SqrtShape0[2];
     for (int i3 = 0; i3 < d2; i3++){
      Index6[5] = i3;
     for (int i2 = 0; i2 < SqrtShape0[2]; i2++){
     for (int i1 = 0; i1 < SqrtShape0[1]; i1++){
      Index6[4] = i1;
     for (int i0 = 0; i0 < SqrtShape0[0]; i0++){
      Index6[1] = i0;
      i = i0 + i1*SqrtShape0[0] + i2*SqrtShape0[0]*SqrtShape0[1] + i3*SqrtShape0[0]*SqrtShape0[1]*SqrtShape0[2];
      NMatrix(i, j) = NTensor1.get(Index6)*IMatrix(i2, j2);
     }}}}
    }}}}
   }
// construct bMatrix:
   bMatrix = Matrix<T>(sqrtDim*d2, 1);
   bTensor1 = bTensor;
   if (Direction == "right"){
// - contract bTensor:
    BTensor = ATensor;
    Indices12[0] = 1; Indices12[1] = 3; Indices22[0] = 1; Indices22[1] = 2;
    bTensor1.contract(Indices12, BTensor, Indices22);
    bTensor2 = bTensor1;
    BTensor = ATensor;
    BTensor.complexConjugate();
    Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 3; Indices23[0] = 1; Indices23[1] = 2; Indices23[2] = 3;
    bTensor2.contract(Indices13, BTensor, Indices23);
    bVal = bTensor2.get(0);
// - build bMatrix:
    for (int i3 = 0; i3 < d2; i3++){
     Index4[3] = i3;
    for (int i2 = 0; i2 < SqrtShape0[2]; i2++){
     Index4[1] = i2;
    for (int i1 = 0; i1 < SqrtShape0[1]; i1++){
     Index4[0] = i1;
    for (int i0 = 0; i0 < SqrtShape0[0]; i0++){
     Index4[2] = i0;
     i = i0 + i1*SqrtShape0[0] + i2*SqrtShape0[0]*SqrtShape0[1] + i3*SqrtShape0[0]*SqrtShape0[1]*SqrtShape0[2];
     bMatrix(i, 0) = bTensor1.get(Index4);
    }}}}
   }
   else if (Direction == "left"){
// - contract bTensor:
    BTensor = ATensor;
    Indices12[0] = 1; Indices12[1] = 3; Indices22[0] = 0; Indices22[1] = 2;
    bTensor1.contract(Indices12, BTensor, Indices22);
    bTensor2 = bTensor1;
    BTensor = ATensor;
    BTensor.complexConjugate();
    Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 3; Indices23[0] = 0; Indices23[1] = 2; Indices23[2] = 3;
    bTensor2.contract(Indices13, BTensor, Indices23);
    bVal = bTensor2.get(0);
// - build bMatrix:
    for (int i3 = 0; i3 < d2; i3++){
     Index4[3] = i3;
    for (int i2 = 0; i2 < SqrtShape0[2]; i2++){
     Index4[1] = i2;
    for (int i1 = 0; i1 < SqrtShape0[1]; i1++){
     Index4[2] = i1;
    for (int i0 = 0; i0 < SqrtShape0[0]; i0++){
     Index4[0] = i0;
     i = i0 + i1*SqrtShape0[0] + i2*SqrtShape0[0]*SqrtShape0[1] + i3*SqrtShape0[0]*SqrtShape0[1]*SqrtShape0[2];
     bMatrix(i, 0) = bTensor1.get(Index4);
    }}}}
   }
// solve system of linear equations via pseudoinverse:
   NMatrixH = NMatrix;
   NMatrixH.adjoint();
   NMatrix.add(NMatrixH);
   NMatrix.multiply(0.5);
   NMatrix.setType("hermitian");
   W = vector<double>(sqrtDim*d2); Vr = Matrix<T>(sqrtDim*d2, sqrtDim*d2);
   NMatrix.eigenDecompose(W, Vr);
   NMatrixInv = Matrix<T>(sqrtDim*d2, sqrtDim*d2);
   NMatrixInv.fillZeroes();
   for (k = 0; k < sqrtDim*d2; k++){
    if (W[k]/abs(W[sqrtDim*d2-1]) > cutoff){
     for (j = 0; j < sqrtDim*d2; j++){
      for (i = 0; i < sqrtDim*d2; i++){
       NMatrixInv(i, j) += Vr(i, k)*MathAuxiliary::complexConjugate(Vr(j, k))/W[k];
      }
     }
    }
   }
   NMatrixInv.multiply(bMatrix);
   for (i = 0; i < sqrtDim*d2; i++)
    ATensor.set(i, NMatrixInv.get(i));
   distance0 = distance1;
   distance1 = NVal-2.0*MathAuxiliary::convertToDouble(bVal);
   if (distance1 > distance0)
    alpha *= x;
   ATensor.multiply(alpha);
   ATensor0.multiply(1.0-alpha);
   ATensor.add(ATensor0);
   ATensor0 = ATensor;
   iter++;
  }
  BTensor = ATensor;
  BTensor.complexConjugate();
  Indices1[0] = 3; Indices2[0] = 3;
  ATensor.contract(Indices1, BTensor, Indices2);
  Order6[0] = 0; Order6[1] = 3; Order6[2] = 1; Order6[3] = 4; Order6[4] = 2; Order6[5] = 5;
  ATensor.permute(Order6);
  ATensor.reshape(Shape0);
  if (Direction == "right")
   this->Tensors[0] = ATensor;
  else if (Direction == "left")
   this->Tensors[this->N-1] = ATensor;
 }
 else
 {
// construct NMatrix and bMatrix:
  unsigned int dim = Shape0[0]*Shape0[1]*Shape0[2], i, j, k;
  Matrix<T> IMatrix(Shape0[2], Shape0[2]), NMatrix(dim, dim), bMatrix(dim, 1);
// - define IMatrix:=1 on physical indices:
  IMatrix.fillZeroes();
  for (i = 0; i < Shape0[2]; i++)
   IMatrix(i, i) = 1.0;
// - fill NMatrix and bMatrix:
  NMatrix.fillZeroes(); bMatrix.fillZeroes();
  vector<unsigned int> Index4(4), Index7(7);
  Index4[0] = 0; Index4[1] = 0;
  Index7[0] = 0; Index7[1] = 0; Index7[2] = 0; Index7[4] = 0; Index7[5] = 0;
  if (Direction == "right")
  {
   for (int i2 = 0; i2 < Shape0[2]; i2++)
   {
    for (int i1 = 0; i1 < Shape0[1]; i1++)
    {
     Index4[3] = i1;
     for (int j2 = 0; j2 < Shape0[2]; j2++)
     {
      for (int j1 = 0; j1 < Shape0[1]; j1++)
      {
       Index4[2] = j1;
       i = i1 + i2*Shape0[1];
       j = j1 + j2*Shape0[1];
       NMatrix(i, j) = NTensor.get(Index4)*IMatrix(i2, j2);
      }
     }
    }
   }
   for (int i2 = 0; i2 < Shape0[2]; i2++)
   {
    Index7[6] = i2;
    for (int i1 = 0; i1 < Shape0[1]; i1++)
    {
     Index7[3] = i1;
     i = i1 + i2*Shape0[1];
     bMatrix(i, 0) = bTensor.get(Index7);
    }
   }
  }
  else if (Direction == "left")
  {
   for (int i2 = 0; i2 < Shape0[2]; i2++)
   {
    for (int i0 = 0; i0 < Shape0[0]; i0++)
    {
     Index4[3] = i0;
     for (int j2 = 0; j2 < Shape0[2]; j2++)
     {
      for (int j0 = 0; j0 < Shape0[0]; j0++)
      {
       Index4[2] = j0;
       i = i0 + i2*Shape0[0];
       j = j0 + j2*Shape0[0];
       NMatrix(i, j) = NTensor.get(Index4)*IMatrix(i2, j2);
      }
     }
    }
   }
   for (int i2 = 0; i2 < Shape0[2]; i2++)
   {
    Index7[6] = i2;
    for (int i0 = 0; i0 < Shape0[0]; i0++)
    {
     Index7[3] = i0;
     i = i0 + i2*Shape0[0];
     bMatrix(i, 0) = bTensor.get(Index7);
    }
   }
  }
// solve system of linear equations:
  if (mode = 0)
  {}
  else if (mode = 1)
  {
// - build closest positive definite matrix to NMatrix:
// -- build closest hermitian NMatrix=(NMatrix+NMatrixH)/2:
   Matrix<T> NMatrixH(NMatrix);
   NMatrixH.adjoint();
   NMatrix.add(NMatrixH);
   NMatrix.multiply(0.5);
// - diagonalize NMatrix and determine its pseudoinverse NMatrixInv:
   NMatrix.setType("hermitian");
   vector<double> W(dim); Matrix<T> Vr(dim, dim);
   NMatrix.eigenDecompose(W, Vr);
   Matrix<T> NMatrixInv(dim, dim); NMatrixInv.fillZeroes();
   for (k = 0; k < dim; k++)
   {
    if (W[k]/abs(W[dim-1]) > cutoff)
    {
     for (j = 0; j < dim; j++)
     {
      for (i = 0; i < dim; i++)
      {
       NMatrixInv(i, j) += Vr(i, k)*MathAuxiliary::complexConjugate(Vr(j, k))/W[k];
      }
     }
    }
   }
// - compute solution:
   NMatrixInv.multiply(bMatrix);
   if (Direction == "right")
   {
    for (i = 0; i < dim; i++)
     this->Tensors[0].set(i, NMatrixInv.get(i));
   }
   else if (Direction == "left")
   {
    for (i = 0; i < dim; i++)
     this->Tensors[this->N-1].set(i, NMatrixInv.get(i));
   }
  }
  if (d2 != 0)
  {
// make solution hermitian:
// - reshape and permute solution Tensor0 to Matrix0:
   Tensor<T> Tensor0;
   if (Direction == "right")
    Tensor0 = this->Tensors[0];
   else if (Direction == "left")
    Tensor0 = this->Tensors[this->N-1];
   vector<unsigned int> Shape6(6);
   Shape6[0] = sqrt(double(Shape0[0])); Shape6[1] = Shape6[0]; Shape6[2] = sqrt(double(Shape0[1]));
   Shape6[3] = Shape6[2]; Shape6[4] = sqrt(double(Shape0[2])); Shape6[5] = Shape6[4];
   Tensor0.reshape(Shape6);
   vector<unsigned int> Order6(6);
   Order6[0] = 0; Order6[1] = 2; Order6[2] = 4; Order6[3] = 1; Order6[4] = 3; Order6[5] = 5;
   Tensor0.permute(Order6);
   unsigned int sqrtDim = sqrt(double(dim));
   vector<unsigned int> Shape2(2);
   Shape2[0] = sqrtDim; Shape2[1] = sqrtDim;
   Tensor0.reshape(Shape2);
   Matrix<T> Matrix0(sqrtDim, sqrtDim);
   vector<unsigned int> Index2(2);
   for (j = 0; j < sqrtDim; j++)
   {
    Index2[1] = j;
    for (i = 0; i < sqrtDim; i++)
    {
     Index2[0] = i;
     Matrix0(i, j) = Tensor0.get(Index2);
    }
   }
// - build closest hermitian matrix to Matrix0:
   Matrix<T> Matrix0H(Matrix0);
   Matrix0H.adjoint();
   double herm = 0.0;
   for (j = 0; j < sqrtDim; j++){
    for (i = 0; i < sqrtDim; i++){
     herm += pow(Matrix0(i, j)-Matrix0H(i, j), 2.0);
    }
   }
   herm = sqrt(abs(herm));
   cout << "hermiticity of boundary solution tensor: " << herm << endl;
   Matrix0.add(Matrix0H);
   Matrix0.multiply(0.5);
// - diagonalize Matrix0 and keep d2 largest eigenvalues and corresponding eigenvectors:
   Matrix0.setType("hermitian");
   vector<double> W(sqrtDim); Matrix<T> Vr(sqrtDim, sqrtDim);
   Matrix0.eigenDecompose(W, Vr);
   cout << "eigenvalues of boundary solution tensor:" << endl;
   for (i = 0; i < sqrtDim; i++){
    cout << "W[" << i << "] = " << W[i] << endl;
   }
   cout << endl;
   d2 = min(d2, sqrtDim);
   Matrix0 = Matrix<T>(sqrtDim, sqrtDim);
   Matrix0.fillZeroes();
   for (j = 0; j < sqrtDim; j++)
   {
    for (i = 0; i < sqrtDim; i++)
    {
     for (k = sqrtDim-d2; k < sqrtDim; k++)
     {
      Matrix0(i, j) += Vr(i, k)*W[k]*MathAuxiliary::complexConjugate(Vr(j, k));
     }
    }
   }
// reshape and permute Matrix0 to solution Tensor0:
   for (j = 0; j < sqrtDim; j++)
   {
    Index2[1] = j;
    for (i = 0; i < sqrtDim; i++)
    {
     Index2[0] = i;
     Tensor0.set(Index2, Matrix0(i, j));
    }
   }
   Shape6[0] = sqrt(double(Shape0[0])); Shape6[1] = sqrt(double(Shape0[1])); Shape6[2] = sqrt(double(Shape0[2]));
   Shape6[3] = Shape6[0]; Shape6[4] = Shape6[1]; Shape6[5] = Shape6[2];
   Tensor0.reshape(Shape6);
   Order6[0] = 0; Order6[1] = 3; Order6[2] = 1; Order6[3] = 4; Order6[4] = 2; Order6[5] = 5;
   Tensor0.permute(Order6);
   Tensor0.reshape(Shape0);
   if (Direction == "right")
    this->Tensors[0] = Tensor0;
   else if (Direction == "left")
    this->Tensors[this->N-1] = Tensor0;
  }
 }
}

template<class T> void MPS<T>::updateTensor(unsigned int position, const string& Direction,
                                            const Tensor<T>& NormTensorLeft,
                                            const Tensor<T>& NormTensorRight, const Tensor<T>& MTensor,
                                            double cutoff, unsigned int mode, unsigned int d2,
                                            double alphaStart, double x, double precision, unsigned int maxNumIter)
{
 Tensor<T> NTensorLeft(NormTensorLeft), NTensorRight(NormTensorRight), bTensor(MTensor);
 vector<unsigned int> Shape0;
 this->Tensors[position].getShape(Shape0);
#ifdef DEBUG
 if ((this->N == 0) || (position > this->N-1) || ((Direction != "right") && (Direction != "left")) ||
     ((mode != 0) && (mode != 1) && (mode != 2)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "updateTensor(unsigned int position, const string& Direction, " <<
                       "const Tensor<T>& NormTensorLeft, const Tensor<T>& NormTensorRight, " <<
                       "const Tensor<T>& MTensor, double cutoff, unsigned int mode, unsigned int d2, " <<
                       "double alphaStart, double x, double precision, unsigned int maxNumIter): " <<
          "((this->N == 0) || (position > this->N-1) || ((Direction != right) && (Direction != left)) || " <<
           "((mode != 0) && (mode != 1) && (mode != 2)))." << endl;
  exit(1);
 }
 vector<unsigned int> NShapeLeft, NShapeRight, bShape;
 NTensorLeft.getShape(NShapeLeft); NTensorRight.getShape(NShapeRight); bTensor.getShape(bShape);
 if ((NTensorLeft.getRank() != 4) || (NShapeLeft[0] != 1) || (NShapeLeft[1] != 1) ||
     (NShapeLeft[2] != Shape0[0]) || (NShapeLeft[3] != Shape0[0]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "updateTensor(unsigned int position, const string& Direction, " <<
                       "const Tensor<T>& NormTensorLeft, const Tensor<T>& NormTensorRight, " <<
                       "const Tensor<T>& MTensor, double cutoff, unsigned int mode, unsigned int d2, " <<
                       "double alphaStart, double x, double precision, unsigned int maxNumIter): " <<
          "NormTensorLeft has wrong Shape." << endl;
  exit(1);
 }
 if ((NTensorRight.getRank() != 4) || (NShapeRight[0] != 1) || (NShapeRight[1] != 1) ||
     (NShapeRight[2] != Shape0[1]) || (NShapeRight[3] != Shape0[1]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "updateTensor(unsigned int position, const string& Direction, " <<
                       "const Tensor<T>& NormTensorLeft, const Tensor<T>& NormTensorRight, " <<
                       "const Tensor<T>& MTensor, double cutoff, unsigned int mode, unsigned int d2, " <<
                       "double alphaStart, double x, double precision, unsigned int maxNumIter): " <<
          "NormTensorRight has wrong Shape." << endl;
  exit(1);
 }
 if ((bTensor.getRank() != 9) || (bShape[0] != 1) || (bShape[1] != 1) || (bShape[2] != 1) ||
     ((bShape[3] != Shape0[0]) && (Direction == "right")) || ((bShape[3] != Shape0[1]) && (Direction == "left")) ||
     (bShape[4] != Shape0[2]) || (bShape[5] != 1) || (bShape[6] != 1) || (bShape[7] != 1) ||
     ((bShape[8] != Shape0[1]) && (Direction == "right")) || ((bShape[8] != Shape0[0]) && (Direction == "left")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void MPS<T>::" <<
          "updateTensor(unsigned int position, const string& Direction, " <<
                       "const Tensor<T>& NormTensorLeft, const Tensor<T>& NormTensorRight, " <<
                       "const Tensor<T>& MTensor, double cutoff, unsigned int mode, unsigned int d2, " <<
                       "double alphaStart, double x, double precision, unsigned int maxNumIter): " <<
          "MTensor has wrong Shape." << endl;
  exit(1);
 }
#endif
 if (mode == 2)
 {
  double alpha, distance0, distance1;
  unsigned int dim = Shape0[0]*Shape0[1]*Shape0[2], sqrtDim = sqrt(double(dim)), i, j, k;
  vector<unsigned int> Indices1(1), Indices12(2), Indices13(3), Indices2(1), Indices22(2), Indices23(3), Shape4(4), Shape6(6), Order6(6);
  vector<unsigned int> Index4(4), Index6(6), Shape2(2), Shape5(5), Index2(2), Indices14(4), Indices24(4);
  vector<double> W;
  Tensor<T> Tensor0, Tensor1, NTensor1, NTensor2, bTensor1, bTensor2;
  Matrix<T> NMatrix, NMatrixH, NMatrixInv, bMatrix;
  Matrix<T> Matrix1(sqrtDim, sqrtDim), Matrix1H, Vr;
// proper d2:
  d2 = min(d2, sqrtDim);
// SqrtShape0 = sqrt(Shape0):
  vector<unsigned int> SqrtShape0(3);
  SqrtShape0[0] = sqrt(double(Shape0[0])); SqrtShape0[1] = sqrt(double(Shape0[1]));
  SqrtShape0[2] = sqrt(double(Shape0[2]));
// define IMatrix:=1 on physical bond:
  Matrix<T> IMatrix(SqrtShape0[2], SqrtShape0[2]);
  IMatrix.fillZeroes();
  for (i = 0; i < SqrtShape0[2]; i++)
   IMatrix(i, i) = 1.0;
// declare purification-Tensor=:ATensor:
  vector<unsigned int> AShape(4);
  AShape[0] = SqrtShape0[0]; AShape[1] = SqrtShape0[1]; AShape[2] = SqrtShape0[2]; AShape[3] = d2;
  Shape5[0] = AShape[0]; Shape5[1] = AShape[1]; Shape5[2] = AShape[2]; Shape5[3] = AShape[3]; Shape5[4] = 1;
  Tensor<T> ATensor(AShape), BTensor;
// reshape NTensorLeft:
  Shape4[0] = SqrtShape0[0]; Shape4[1] = Shape4[0]; Shape4[2] = SqrtShape0[0]; Shape4[3] = Shape4[2];
  NTensorLeft.reshape(Shape4);
// reshape NTensorRight:
  Shape4[0] = SqrtShape0[1]; Shape4[1] = Shape4[0]; Shape4[2] = SqrtShape0[1]; Shape4[3] = Shape4[2];
  NTensorRight.reshape(Shape4);
// reshape bTensor:
  if (Direction == "right"){
   Shape6[0] = SqrtShape0[0]; Shape6[1] = Shape6[0]; Shape6[2] = SqrtShape0[2]; Shape6[3] = Shape6[2];
   Shape6[4] = SqrtShape0[1]; Shape6[5] = Shape6[4];
  }
  else if (Direction == "left"){
   Shape6[0] = SqrtShape0[1]; Shape6[1] = Shape6[0]; Shape6[2] = SqrtShape0[2]; Shape6[3] = Shape6[2];
   Shape6[4] = SqrtShape0[0]; Shape6[5] = Shape6[4];
  }
  bTensor.reshape(Shape6);
// initialize iteration and obtain ATensor:
  Tensor0 = this->Tensors[position];
  Tensor1 = Tensor0;
  Shape6[0] = SqrtShape0[0]; Shape6[1] = Shape6[0]; Shape6[2] = SqrtShape0[1]; Shape6[3] = Shape6[2];
  Shape6[4] = SqrtShape0[2]; Shape6[5] = Shape6[4];
  Tensor1.reshape(Shape6);
  Order6[0] = 0; Order6[1] = 2; Order6[2] = 4; Order6[3] = 1; Order6[4] = 3; Order6[5] = 5;
  Tensor1.permute(Order6);
  Shape2[0] = sqrtDim; Shape2[1] = sqrtDim;
  Tensor1.reshape(Shape2);
  for (j = 0; j < sqrtDim; j++){
   Index2[1] = j;
   for (i = 0; i < sqrtDim; i++){
    Index2[0] = i;
    Matrix1(i, j) = Tensor1.get(Index2);
   }
  }
  Matrix1H = Matrix1;
  Matrix1H.adjoint();
  Matrix1.add(Matrix1H);
  Matrix1.multiply(0.5);
  Matrix1.setType("hermitian");
  W = vector<double>(sqrtDim); Vr = Matrix<T>(sqrtDim, sqrtDim);
  Matrix1.eigenDecompose(W, Vr);
  for (i = 0; i < sqrtDim; i++)
  {
   if (W[i]/abs(W[sqrtDim-1]) < cutoff)
    W[i] = 0.0;
  }
  for (k = sqrtDim-d2; k < sqrtDim; k++){
   for (i = 0; i < sqrtDim; i++){
    ATensor.set(i+(k-sqrtDim+d2)*sqrtDim, Vr(i, k)*sqrt(W[k]));
   }
  }
  Tensor<T> ATensor0(ATensor);
  T NVal, bVal;
  distance0 = 0.0; distance1 = 1.0e10;
  alpha = alphaStart; unsigned int iter = 0;
// iterate:
  while ((abs(distance0-distance1)/abs(distance1) > precision) && (iter < maxNumIter))
  {
// construct NMatrix:
   NMatrix = Matrix<T>(sqrtDim*d2, sqrtDim*d2);
// - contract NTensor:
   NTensor1 = NTensorLeft;
   BTensor = ATensor;
   BTensor.complexConjugate();
   Indices1[0] = 1; Indices2[0] = 0;
   NTensor1.contract(Indices1, BTensor, Indices2);
   BTensor = ATensor;
   Indices12[0] = 2; Indices12[1] = 4; Indices22[0] = 0; Indices22[1] = 2;
   NTensor1.contract(Indices12, BTensor, Indices22);
   Tensor0 = NTensorRight;
   Indices12[0] = 2; Indices12[1] = 4; Indices22[0] = 1; Indices22[1] = 3;
   NTensor1.contract(Indices12, Tensor0, Indices22);
   NTensor2 = NTensor1;
   BTensor = ATensor;
   Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 4; Indices23[0] = 0; Indices23[1] = 3; Indices23[2] = 1;
   NTensor2.contract(Indices13, BTensor, Indices23);
   BTensor = ATensor;
   BTensor.complexConjugate();
   BTensor.reshape(Shape5);
   Indices14[0] = 0; Indices14[1] = 1; Indices14[2] = 2; Indices14[3] = 3; Indices24[0] = 0; Indices24[1] = 3; Indices24[2] = 1; Indices24[3] = 2;
   NTensor2.contract(Indices14, BTensor, Indices24);
   NVal = NTensor2.get(0);
// - build NMatrix:
   for (int j3 = 0; j3 < d2; j3++){
    Index6[2] = j3;
   for (int j2 = 0; j2 < SqrtShape0[2]; j2++){
   for (int j1 = 0; j1 < SqrtShape0[1]; j1++){
    Index6[4] = j1;
   for (int j0 = 0; j0 < SqrtShape0[0]; j0++){
    Index6[0] = j0;
    j = j0 + j1*SqrtShape0[0] + j2*SqrtShape0[0]*SqrtShape0[1] + j3*SqrtShape0[0]*SqrtShape0[1]*SqrtShape0[2];
    for (int i3 = 0; i3 < d2; i3++){
     Index6[3] = i3;
    for (int i2 = 0; i2 < SqrtShape0[2]; i2++){
    for (int i1 = 0; i1 < SqrtShape0[1]; i1++){
     Index6[5] = i1;
    for (int i0 = 0; i0 < SqrtShape0[0]; i0++){
     Index6[1] = i0;
     i = i0 + i1*SqrtShape0[0] + i2*SqrtShape0[0]*SqrtShape0[1] + i3*SqrtShape0[0]*SqrtShape0[1]*SqrtShape0[2];
     NMatrix(i, j) = NTensor1.get(Index6)*IMatrix(i2, j2);
    }}}}
   }}}}
// construct bMatrix:
   bMatrix = Matrix<T>(sqrtDim*d2, 1);
   bTensor1 = bTensor;
   if (Direction == "right"){
// - contract bTensor:
    BTensor = ATensor;
    Indices13[0] = 1; Indices13[1] = 3; Indices13[2] = 5; Indices23[0] = 0; Indices23[1] = 2; Indices23[2] = 1;
    bTensor1.contract(Indices13, BTensor, Indices23);
    bTensor2 = bTensor1;
    BTensor = ATensor;
    BTensor.complexConjugate();
    BTensor.reshape(Shape5);
    Indices14[0] = 0; Indices14[1] = 1; Indices14[2] = 2; Indices14[3] = 3; Indices24[0] = 0; Indices24[1] = 2; Indices24[2] = 1; Indices24[3] = 3;
    bTensor2.contract(Indices14, BTensor, Indices24);
    bVal = bTensor2.get(0);
// - build bMatrix:
    for (int i3 = 0; i3 < d2; i3++){
     Index4[3] = i3;
    for (int i2 = 0; i2 < SqrtShape0[2]; i2++){
     Index4[1] = i2;
    for (int i1 = 0; i1 < SqrtShape0[1]; i1++){
     Index4[2] = i1;
    for (int i0 = 0; i0 < SqrtShape0[0]; i0++){
     Index4[0] = i0;
     i = i0 + i1*SqrtShape0[0] + i2*SqrtShape0[0]*SqrtShape0[1] + i3*SqrtShape0[0]*SqrtShape0[1]*SqrtShape0[2];
     bMatrix(i, 0) = bTensor1.get(Index4);
    }}}}
   }
   else if (Direction == "left"){
// - contract bTensor:
    BTensor = ATensor;
    Indices13[0] = 1; Indices13[1] = 3; Indices13[2] = 5; Indices23[0] = 1; Indices23[1] = 2; Indices23[2] = 0;
    bTensor1.contract(Indices13, BTensor, Indices23);
    bTensor2 = bTensor1;
    BTensor = ATensor;
    BTensor.complexConjugate();
    BTensor.reshape(Shape5);
    Indices14[0] = 0; Indices14[1] = 1; Indices14[2] = 2; Indices14[3] = 3; Indices24[0] = 1; Indices24[1] = 2; Indices24[2] = 0; Indices24[3] = 3;
    bTensor2.contract(Indices14, BTensor, Indices24);
    bVal = bTensor2.get(0);
// - build bMatrix:
    for (int i3 = 0; i3 < d2; i3++){
     Index4[3] = i3;
    for (int i2 = 0; i2 < SqrtShape0[2]; i2++){
     Index4[1] = i2;
    for (int i1 = 0; i1 < SqrtShape0[1]; i1++){
     Index4[0] = i1;
    for (int i0 = 0; i0 < SqrtShape0[0]; i0++){
     Index4[2] = i0;
     i = i0 + i1*SqrtShape0[0] + i2*SqrtShape0[0]*SqrtShape0[1] + i3*SqrtShape0[0]*SqrtShape0[1]*SqrtShape0[2];
     bMatrix(i, 0) = bTensor1.get(Index4);
    }}}}
   }
// solve system of linear equations via pseudoinverse:
   NMatrixH = NMatrix;
   NMatrixH.adjoint();
   NMatrix.add(NMatrixH);
   NMatrix.multiply(0.5);
   NMatrix.setType("hermitian");
   W = vector<double>(sqrtDim*d2); Vr = Matrix<T>(sqrtDim*d2, sqrtDim*d2);
   NMatrix.eigenDecompose(W, Vr);
   NMatrixInv = Matrix<T>(sqrtDim*d2, sqrtDim*d2);
   NMatrixInv.fillZeroes();
   for (k = 0; k < sqrtDim*d2; k++){
    if (W[k]/abs(W[sqrtDim*d2-1]) > cutoff){
     for (j = 0; j < sqrtDim*d2; j++){
      for (i = 0; i < sqrtDim*d2; i++){
       NMatrixInv(i, j) += Vr(i, k)*MathAuxiliary::complexConjugate(Vr(j, k))/W[k];
      }
     }
    }
   }
   NMatrixInv.multiply(bMatrix);
   for (i = 0; i < sqrtDim*d2; i++)
    ATensor.set(i, NMatrixInv.get(i));
   distance0 = distance1;
   distance1 = NVal-2.0*MathAuxiliary::convertToDouble(bVal);
   if (distance1 > distance0)
    alpha *= x;
   ATensor.multiply(alpha);
   ATensor0.multiply(1.0-alpha);
   ATensor.add(ATensor0);
   ATensor0 = ATensor;
   iter++;
  }
  BTensor = ATensor;
  BTensor.complexConjugate();
  Indices1[0] = 3; Indices2[0] = 3;
  ATensor.contract(Indices1, BTensor, Indices2);
  Order6[0] = 0; Order6[1] = 3; Order6[2] = 1; Order6[3] = 4; Order6[4] = 2; Order6[5] = 5;
  ATensor.permute(Order6);
  ATensor.reshape(Shape0);
  this->Tensors[position] = ATensor;
 }
 else
 {
// construct NMatrix and bMatrix:
  unsigned int dim = Shape0[0]*Shape0[1]*Shape0[2], i, j, k;
  Matrix<T> IMatrix(Shape0[2], Shape0[2]), NMatrix(dim, dim), bMatrix(dim, 1);
// - define IMatrix:=1 on physical indices:
  IMatrix.fillZeroes();
  for (i = 0; i < Shape0[2]; i++)
   IMatrix(i, i) = 1.0;
// - fill NMatrix:
  NMatrix.fillZeroes();
  vector<unsigned int> Index4Left(4), Index4Right(4);
  Index4Left[0] = 0; Index4Left[1] = 0; Index4Right[0] = 0; Index4Right[1] = 0;
  for (int i2 = 0; i2 < Shape0[2]; i2++)
  {
  for (int i1 = 0; i1 < Shape0[1]; i1++)
  {
   Index4Right[3] = i1;
  for (int i0 = 0; i0 < Shape0[0]; i0++)
  {
   Index4Left[3] = i0;
   for (int j2 = 0; j2 < Shape0[2]; j2++)
   {
   for (int j1 = 0; j1 < Shape0[1]; j1++)
   {
    Index4Right[2] = j1;
   for (int j0 = 0; j0 < Shape0[0]; j0++)
   {
    Index4Left[2] = j0;
    i = i0 + i1*Shape0[0] + i2*Shape0[0]*Shape0[1];
    j = j0 + j1*Shape0[0] + j2*Shape0[0]*Shape0[1];
    NMatrix(i, j) = NTensorLeft.get(Index4Left)*IMatrix(i2, j2)*NTensorRight.get(Index4Right);
   }
   }
   }
  }
  }
  }
// - fill bMatrix:
  bMatrix.fillZeroes();
  vector<unsigned int> Index9(9);
  Index9[0] = 0; Index9[1] = 0; Index9[2] = 0; Index9[5] = 0; Index9[6] = 0; Index9[7] = 0;
  if (Direction == "right")
  {
   for (int i1 = 0; i1 < Shape0[1]; i1++)
   {
    Index9[8] = i1;
   for (int i2 = 0; i2 < Shape0[2]; i2++)
   {
    Index9[4] = i2;
   for (int i0 = 0; i0 < Shape0[0]; i0++)
   {
    Index9[3] = i0;
    i = i0 + i1*Shape0[0] + i2*Shape0[0]*Shape0[1];
    bMatrix(i, 0) = bTensor.get(Index9);
   }
   }
   }
  }
  else if (Direction == "left")
  {
   for (int i0 = 0; i0 < Shape0[0]; i0++)
   {
    Index9[8] = i0;
   for (int i2 = 0; i2 < Shape0[2]; i2++)
   {
    Index9[4] = i2;
   for (int i1 = 0; i1 < Shape0[1]; i1++)
   {
    Index9[3] = i1;
    i = i0 + i1*Shape0[0] + i2*Shape0[0]*Shape0[1];
    bMatrix(i, 0) = bTensor.get(Index9);
   }
   }
   }
  }
// solve system of linear equations:
  if (mode = 0)
  {}
  else if (mode = 1)
  {
// - build closest positive definite matrix to NMatrix:
// -- build closest hermitian NMatrix=(NMatrix+NMatrixH)/2:
   Matrix<T> NMatrixH(NMatrix);
   NMatrixH.adjoint();
   NMatrix.add(NMatrixH);
   NMatrix.multiply(0.5);
// - diagonalize NMatrix and determine its pseudoinverse NMatrixInv:
   NMatrix.setType("hermitian");
   vector<double> W(dim); Matrix<T> Vr(dim, dim);
   NMatrix.eigenDecompose(W, Vr);
   Matrix<T> NMatrixInv(dim, dim); NMatrixInv.fillZeroes();
   for (k = 0; k < dim; k++)
   {
    if (W[k]/abs(W[dim-1]) > cutoff)
    {
     for (j = 0; j < dim; j++)
     {
      for (i = 0; i < dim; i++)
      {
       NMatrixInv(i, j) += Vr(i, k)*MathAuxiliary::complexConjugate(Vr(j, k))/W[k];
      }
     }
    }
   }
// - compute solution:
   NMatrixInv.multiply(bMatrix);
   for (i = 0; i < dim; i++)
    this->Tensors[position].set(i, NMatrixInv.get(i));
  }
  if (d2 != 0)
  {
// make solution hermitian:
// - reshape and permute solution Tensor0 to Matrix0:
   Tensor<T> Tensor0(this->Tensors[position]);
   vector<unsigned int> Shape6(6);
   Shape6[0] = sqrt(double(Shape0[0])); Shape6[1] = Shape6[0]; Shape6[2] = sqrt(double(Shape0[1]));
   Shape6[3] = Shape6[2]; Shape6[4] = sqrt(double(Shape0[2])); Shape6[5] = Shape6[4];
   Tensor0.reshape(Shape6);
   vector<unsigned int> Order6(6);
   Order6[0] = 0; Order6[1] = 2; Order6[2] = 4; Order6[3] = 1; Order6[4] = 3; Order6[5] = 5;
   Tensor0.permute(Order6);
   unsigned int sqrtDim = sqrt(double(dim));
   vector<unsigned int> Shape2(2);
   Shape2[0] = sqrtDim; Shape2[1] = sqrtDim;
   Tensor0.reshape(Shape2);
   Matrix<T> Matrix0(sqrtDim, sqrtDim);
   vector<unsigned int> Index2(2);
   for (j = 0; j < sqrtDim; j++)
   {
    Index2[1] = j;
    for (i = 0; i < sqrtDim; i++)
    {
     Index2[0] = i;
     Matrix0(i, j) = Tensor0.get(Index2);
    }
   }
// - build closest hermitian matrix to Matrix0:
   Matrix<T> Matrix0H(Matrix0);
   Matrix0H.adjoint();
   double herm = 0.0;
   for (j = 0; j < sqrtDim; j++){
    for (i = 0; i < sqrtDim; i++){
     herm += pow(Matrix0(i, j)-Matrix0H(i, j), 2.0);
    }
   }
   herm = sqrt(abs(herm));
   cout << "hermiticity of bulk solution tensor at position " << position << ": " << herm << endl;
   Matrix0.add(Matrix0H);
   Matrix0.multiply(0.5);
// - diagonalize Matrix0 and keep d2 largest eigenvalues and corresponding eigenvectors:
   Matrix0.setType("hermitian");
   vector<double> W(sqrtDim); Matrix<T> Vr(sqrtDim, sqrtDim);
   Matrix0.eigenDecompose(W, Vr);
   cout << "eigenvalues of bulk solution tensor at position " << position << ":" << endl;
   for (i = 0; i < sqrtDim; i++){
    cout << "W[" << i << "] = " << W[i] << endl;
   }
   cout << endl;
   d2 = min(d2, sqrtDim);
   Matrix0 = Matrix<T>(sqrtDim, sqrtDim);
   Matrix0.fillZeroes();
   for (j = 0; j < sqrtDim; j++)
   {
    for (i = 0; i < sqrtDim; i++)
    {
     for (k = sqrtDim-d2; k < sqrtDim; k++)
     {
      Matrix0(i, j) += Vr(i, k)*W[k]*MathAuxiliary::complexConjugate(Vr(j, k));
     }
    }
   }
// reshape and permute Matrix0 to solution Tensor0:
   for (j = 0; j < sqrtDim; j++)
   {
    Index2[1] = j;
    for (i = 0; i < sqrtDim; i++)
    {
     Index2[0] = i;
     Tensor0.set(Index2, Matrix0(i, j));
    }
   }
   Shape6[0] = sqrt(double(Shape0[0])); Shape6[1] = sqrt(double(Shape0[1])); Shape6[2] = sqrt(double(Shape0[2]));
   Shape6[3] = Shape6[0]; Shape6[4] = Shape6[1]; Shape6[5] = Shape6[2];
   Tensor0.reshape(Shape6);
   Order6[0] = 0; Order6[1] = 3; Order6[2] = 1; Order6[3] = 4; Order6[4] = 2; Order6[5] = 5;
   Tensor0.permute(Order6);
   Tensor0.reshape(Shape0);
   this->Tensors[position] = Tensor0;
  }
 }
}

template<class T> void multiplyEvenOddMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, MPS<T>& MPS1)
{
#ifdef DEBUG
 if ((MPO0.BC != MPS0.BC) || (MPO0.BC != MPS1.BC) || (MPO0.N != MPS0.N) || (MPO0.N != MPS1.N) ||
     (MPO0.N == 0) || (MPO0.d != MPS0.d) || (MPO0.d != MPS1.d) || (MPS0.D != MPS1.D))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> void " <<
          "multiplyEvenOddMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, MPS<T>& MPS1): " <<
          "MPO0, MPS0 and MPS1 are not of the same form or (MPO0.N == 0)." << endl;
  exit(1);
 }
#endif
 if (MPO0.BC == "periodic")
 {
  cout << "The following function is not implemented yet for periodic boundary conditions: " <<
          "template<class T> void " <<
          "multiplyEvenOddMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, MPS<T>& MPS1)." << endl;
  return;
 }
 Tensor<T> Tensor0, Tensor1;
 Matrix<T> TwoBodyMatrix, U, Sigma, Vt;
 vector<unsigned int> Index(6), Index03(3), Index13(3);
 vector<unsigned int> ShapeMPSLeft(3), ShapeMPSRight(3), ShapeMPOLeft(4), ShapeMPORight(4);
 vector<unsigned int> Indices0(1), Indices1(1), Indices02(2), Indices12(2);
 unsigned int dim0, dim1;
// We check if it is the even decomposition {0-1, 2-3, 4-5, ...} or the odd {1-2, 3-4, ...}:
 MPO0.Tensors[0].getShape(ShapeMPOLeft); MPO0.Tensors[1].getShape(ShapeMPORight);
// if ((ShapeMPOLeft[0] == 1) && (ShapeMPOLeft[1] != 1)) it is the even decomposition, else if
// ((ShapeMPORight[0] == 1) && (ShapeMPORight[1] != 1)) it is the odd decomposition:
 if ((ShapeMPOLeft[0] == 1) && (ShapeMPOLeft[1] != 1))
 {
  for (int position = 0; position < MPO0.N-1; position += 2)
  {
   MPS0.Tensors[position].getShape(ShapeMPSLeft);
   MPS0.Tensors[position+1].getShape(ShapeMPSRight);
   MPO0.Tensors[position].getShape(ShapeMPOLeft);
   MPO0.Tensors[position+1].getShape(ShapeMPORight);
   Tensor0 = MPS0.Tensors[position]; Tensor1 = MPO0.Tensors[position];
   Indices0[0] = 2; Indices1[0] = 2;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = MPS0.Tensors[position+1];
   Indices0[0] = 1; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = MPO0.Tensors[position+1];
   Indices02[0] = 2; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 2;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   dim0 = ShapeMPSLeft[0]*ShapeMPOLeft[3];
   dim1 = ShapeMPSRight[1]*ShapeMPORight[3];
   TwoBodyMatrix = Matrix<T>(dim0, dim1); U = Matrix<T>(dim0, dim0); Sigma = Matrix<T>(dim0, dim1);
   Vt = Matrix<T>(dim1, dim1);
   Index[1] = 0; Index[4] = 0;
   for (int i = 0; i < dim0; i++)
   {
    Index[0] = i % ShapeMPSLeft[0]; Index[2] = i / ShapeMPSLeft[0];
    for (int j = 0; j < dim1; j++)
    {
     Index[3] = j % ShapeMPSRight[1]; Index[5] = j / ShapeMPSRight[1];
     TwoBodyMatrix(i, j) = Tensor0.get(Index);
    }
   }
   TwoBodyMatrix.singularValueDecompose(U, Sigma, Vt);
   Tensor0 = MPS1.Tensors[position]; Tensor1 = MPS1.Tensors[position+1];
   for (int i = 0; i < dim0; i++)
   {
    Index03[0] = i % ShapeMPSLeft[0]; Index03[2] = i / ShapeMPSLeft[0];
    for (int j = 0; j < dim1; j++)
    {
     Index13[1] = j % ShapeMPSRight[1]; Index13[2] = j / ShapeMPSRight[1];
     for (int k = 0; k < ShapeMPSLeft[1]; k++)
     {
      Index03[1] = k; Index13[0] = k;
      Tensor0.set(Index03, U(i, k)*sqrt(Sigma(k, k)));
      Tensor1.set(Index13, sqrt(Sigma(k, k))*Vt(k, j));
     }
    }
   }
   MPS1.set(position, Tensor0);
   MPS1.set(position+1, Tensor1);
  }
 }
 else if ((ShapeMPORight[0] == 1) && (ShapeMPORight[1] != 1))
 {
  for (int position = 1; position < MPO0.N-1; position += 2)
  {
   MPS0.Tensors[position].getShape(ShapeMPSLeft);
   MPS0.Tensors[position+1].getShape(ShapeMPSRight);
   MPO0.Tensors[position].getShape(ShapeMPOLeft);
   MPO0.Tensors[position+1].getShape(ShapeMPORight);
   Tensor0 = MPS0.Tensors[position]; Tensor1 = MPO0.Tensors[position];
   Indices0[0] = 2; Indices1[0] = 2;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = MPS0.Tensors[position+1];
   Indices0[0] = 1; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = MPO0.Tensors[position+1];
   Indices02[0] = 2; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 2;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   dim0 = ShapeMPSLeft[0]*ShapeMPOLeft[3];
   dim1 = ShapeMPSRight[1]*ShapeMPORight[3];
   TwoBodyMatrix = Matrix<T>(dim0, dim1); U = Matrix<T>(dim0, dim0); Sigma = Matrix<T>(dim0, dim1);
   Vt = Matrix<T>(dim1, dim1);
   Index[1] = 0; Index[4] = 0;
   for (int i = 0; i < dim0; i++)
   {
    Index[0] = i % ShapeMPSLeft[0]; Index[2] = i / ShapeMPSLeft[0];
    for (int j = 0; j < dim1; j++)
    {
     Index[3] = j % ShapeMPSRight[1]; Index[5] = j / ShapeMPSRight[1];
     TwoBodyMatrix(i, j) = Tensor0.get(Index);
    }
   }
   TwoBodyMatrix.singularValueDecompose(U, Sigma, Vt);
   Tensor0 = MPS1.Tensors[position]; Tensor1 = MPS1.Tensors[position+1];
   for (int i = 0; i < dim0; i++)
   {
    Index03[0] = i % ShapeMPSLeft[0]; Index03[2] = i / ShapeMPSLeft[0];
    for (int j = 0; j < dim1; j++)
    {
     Index13[1] = j % ShapeMPSRight[1]; Index13[2] = j / ShapeMPSRight[1];
     for (int k = 0; k < ShapeMPSLeft[1]; k++)
     {
      Index03[1] = k; Index13[0] = k;
      Tensor0.set(Index03, U(i, k)*sqrt(Sigma(k, k)));
      Tensor1.set(Index13, sqrt(Sigma(k, k))*Vt(k, j));
     }
    }
   }
   MPS1.set(position, Tensor0);
   MPS1.set(position+1, Tensor1);
  }
 }
}

template<class T> void getTensor(const MPS<T>& MPS0, Tensor<T>& Tensor0)
{
#ifdef DEBUG
 if (MPS0.N == 0)
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> void " <<
          "getTensor(const MPS<T>& MPS0, Tensor<T>& Tensor0): " <<
          "(MPS0.N == 0)." << endl;
  exit(1);
 }
#endif
 vector<unsigned int> Index0(1), Index1(1);
 Index1[0] = 0;
 Tensor0 = MPS0.Tensors[0];
 Tensor<T> Tensor1;
 for (int l = 1; l < MPS0.N; l++)
 {
  Tensor1 = MPS0.Tensors[l];
  Index0[0] = l;
  Tensor0.contract(Index0, Tensor1, Index1);
 }
 if (MPS0.BC == "open")
 {
  unsigned int rank0 = Tensor0.rank-2;
  vector<unsigned int> Shape0(rank0);
  for (int l = 0; l < rank0-1; l++)
   Shape0[l] = Tensor0.Shape[l+1];
  Shape0[rank0-1] = Tensor0.Shape[rank0+1];
  Tensor0.rank = rank0;
  Tensor0.Shape = Shape0;
 }
 else if (MPS0.BC == "periodic")
 {
  unsigned int index0 = 0, index1 = Tensor0.rank-2;
  Tensor0.contract(index0, index1);
 }
}

template<class T> inline void MPS<T>::getOpenBCShape(unsigned int position,
                                                     vector<unsigned int>& Shape0) const
{
#ifdef DEBUG
 if ((this->BC != "open") || (this->N == 0) || (position >= this->N) || (Shape0.size() != 3))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline void MPS<T>::" <<
          "getOpenBCShape(unsigned int position, vector<unsigned int>& Shape0) const: " <<
          "((this->BC != open) || (this->N == 0) || (position >= this->N) || " <<
          "(Shape0.size() != 3))." << endl;
  exit(1);
 }
#endif
 double logdSumRight = 0.0; double logD = log((double)this->D); double logdSumLeft = 0.0;
 for (int i = 0; i < (int)position; i++)
 {
  logdSumRight += log((double)this->d[i]);
 }
 for (int i = this->N-1; i >= (int)position; i--)
 {
  logdSumLeft += log((double)this->d[i]);
 }
 if (logdSumRight <= min(logD, logdSumLeft))
 {
  Shape0[0] = 1;
  for (int i = 0; i < (int)position; i++)
  {
   Shape0[0] *= this->d[i];
  }
 }
 else if (logD <= min(logdSumRight, logdSumLeft))
 {
  Shape0[0] = this->D;
 }
 else if (logdSumLeft < min(logdSumRight, logD))
 {
  Shape0[0] = 1;
  for (int i = this->N-1; i >= (int)position; i--)
  {
   Shape0[0] *= this->d[i];
  }
 }
 logdSumRight += log((double)this->d[position]);
 logdSumLeft -= log((double)this->d[position]);
 if (logdSumRight <= min(logD, logdSumLeft))
 {
  Shape0[1] = 1;
  for (int i = 0; i <= (int)position; i++)
  {
   Shape0[1] *= this->d[i];
  }
 }
 else if (logD <= min(logdSumRight, logdSumLeft))
 {
  Shape0[1] = this->D;
 }
 else if (logdSumLeft < min(logdSumRight, logD))
 {
  Shape0[1] = 1;
  for (int i = this->N-1; i > (int)position; i--)
  {
   Shape0[1] *= this->d[i];
  }
 }
 Shape0[2] = this->d[position];
}
