/// Template class PEPS implements Projected Entangled-Pair States.
/** The template class PEPS implements Projected Entangled-Pair States in 2D.
    PEPS consist of Tensors in a matrix.
    \param BC string, the boundary conditions of the PEPS, is "open" or "periodic"
    \param N vector<unsigned int>(2), the number of Tensors of the PEPS, fulfills N.size()==2
    \param d Matrix<unsigned int>(N[0], N[1]), the physical dimensions of the PEPS, fulfills
             d.getDim0()==N[0] and d.getDim1()==N[1]
    \param D unsigned int, the maximal virtual bond dimension of the PEPS
    \param Tensors Tensor<T>*, the tensors of the PEPS, stored in Fortran's column-major order
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

// Declaration of template classes:
template<class T> class PEPS;
template<class T> class PEPO;

// Declaration of friend functions:
template<class T> double distancePEPOPEPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const PEPS<T>& PEPS1,
                                          const string& Direction, const vector<unsigned int>& D2s,
                                          const vector<double>& Epss,
                                          const vector<unsigned int>& MaxNumsSweeps, vector<double>& ErrorsAchieved,
                                          vector<unsigned int>& NumsSweepsDone);
template<class T> void multiplyPEPOPEPSExact(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, PEPS<T>& PEPS1);
template<class T> void multiplyPEPOPEPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0,
                                        double eps, unsigned int maxNumSweeps,
                                        unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv,
                                        const string& UpdateTensor, const string& UpdateMode, double cutoff,
                                        double& epsAchieved, unsigned int& numSweepsDone,
                                        PEPS<T>& PEPS1);

template<class T> class PEPS
{
 public:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0.
    \sa PEPS(const string& BC0, unsigned int Nrows, unsigned int Ncols, const Matrix<unsigned int>& d0,
             unsigned int D0)
    \sa PEPS(const string& BC0, unsigned int Nrows, unsigned int Ncols, unsigned int d0,
             unsigned int D0) */
  PEPS();

/// Constructor for PEPS with specific BC, N, d and D.
/** This constructor initializes a PEPS of a specific form with Nrows rows and Ncols columns.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param Nrows input: unsigned int, the number of rows, must be > 1
    \param Ncols input: unsigned int, the number of columns, must be > 1
    \param d0 input: const Matrix<unsigned int>&, the physical dimensions, must fulfill
                     d0.getDim0()==Nrows and d0.getDim1()==Ncols, and all entries must be > 1
    \param D0 input: unsigned int, the maximal virtual bond dimension, must be > 0
    \sa PEPS()
    \sa PEPS(const string& BC0, unsigned int Nrows, unsigned int Ncols, unsigned int d0,
             unsigned int D0) */
  PEPS(const string& BC0, unsigned int Nrows, unsigned int Ncols, const Matrix<unsigned int>& d0,
       unsigned int D0);

/// Constructor for PEPS with specific BC, N, d and D.
/** This constructor initializes a PEPS of a specific form with Nrows rows, Ncols columns and all
    physical dimensions equal to d0.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param Nrows input: unsigned int, the number of rows, must be > 1
    \param Ncols input: unsigned int, the number of columns, must be > 1
    \param d0 input: unsigned int, the physical dimension, must be > 1
    \param D0 input: unsigned int, the maximal virtual bond dimension, must be > 0
    \sa PEPS()
    \sa PEPS(const string& BC0, unsigned int Nrows, unsigned int Ncols, const Matrix<unsigned int>& d0,
             unsigned int D0) */
  PEPS(const string& BC0, unsigned int Nrows, unsigned int Ncols, unsigned int d0, unsigned int D0);

/// Standard copy constructor.
/** The standard copy constructor copies the input PEPS into this.
    \param PEPS0 input: const PEPS<T>&, to be copied into this
    \sa PEPS<T>& operator=(const PEPS<T>& PEPS0) */
  PEPS(const PEPS<T>& PEPS0);

/// Standard destructor.
/** The standard destructor deletes the elements of the PEPS. */
  ~PEPS();

/// Assigns PEPS to this.
/** The operator= allows to assign a PEPS to this. Hereby this is destroyed and newly constructed to be
    a copy of the right-hand side PEPS.
    \param PEPS0 input: const PEPS<T>&, to be copied into this
    \return PEPS<T>&, a reference to the new this
    \sa PEPS(const PEPS<T>& PEPS0) */
  PEPS<T>& operator=(const PEPS<T>& PEPS0);

/// Returns boundary conditions.
/** The returned boundary conditions are either "open" or "periodic".
    \param BC0 output: string&, the boundary conditions of this PEPS */
  void getBC(string& BC0) const { BC0 = this->BC; }

/// Returns number of rows Nrows.
/** This function returns the number of rows of this PEPS.
    \return unsigned int, the number of rows of this PEPS */
  unsigned int getNrows() const { return this->N[0]; }

/// Returns number of columns Ncols.
/** This function returns the number of columns of this PEPS.
    \return unsigned int, the number of columns of this PEPS */
  unsigned int getNcols() const { return this->N[1]; }

/// Returns physical dimensions d.
/** This function returns the physical dimensions of this PEPS.
    \param d0 output: Matrix<unsigned int>&, the physical dimensions of this PEPS */
  void getd(Matrix<unsigned int>& d0) const { d0 = this->d; }

/// Returns physical dimension d.
/** This function returns the physical dimension of this PEPS at site (0, 0). It is useful, if this
    PEPS has all physical dimensions equal to each other.
    \return unsigned int, the physical dimension of this PEPS at site (0, 0) */
  unsigned int getd() const { return this->d(0, 0); }

/// Sets maximal virtual bond dimension.
/** This function changes the maximal virtual bond dimension of this PEPS.
    If the maximal virtual bond dimension is decreased, D0 < this->D, then the dispensable elements are
    simply discarded.
    If it is increased, D0 > this->D, then the new elements are set as random numbers multiplied by element.
    Element has the default value 0.0, in which case the new elements are zeroes.
    The random numbers for a tensor at position (row, col) are seeded with time(0)+(row+this->N[0]*col)*13.
    \param D0 input: unsigned int, the new maximal virtual bond dimension, must be > 0
    \param element optional input: T, multiplies all random numbers */
  void setD(unsigned int D0, T element = 0.0);

/// Returns maximal virtual bond dimension D.
/** This function returns the maximal virtual bond dimension of this PEPS.
    \return unsigned int, the maximal virtual bond dimension of this PEPS */
  unsigned int getD() const { return this->D; }

/// Sets Tensor at position.
/** This function sets Tensor0 at a given row number positionRow and column number positionCol in this
    PEPS.
    \param positionRow input: unsigned int, the row position
    \param positionCol input: unsigned int, the column position
    \param Tensor0 input: const Tensor<T>&, to be written at position, must have the correct Shape */
  void set(unsigned int positionRow, unsigned int positionCol, const Tensor<T>& Tensor0);

/// Returns Tensor at position.
/** This function returns as Tensor0 a copy of the tensor at a given row number positionRow and column
    number positionCol in this PEPS.
    \param positionRow input: unsigned int, the row position
    \param positionCol input: unsigned int, the column position
    \param Tensor0 output: Tensor<T>&, a copy of the tensor at position */
  void get(unsigned int positionRow, unsigned int positionCol, Tensor<T>& Tensor0) const;

/// Writes this PEPS to binary file.
/** Given a file name FileName, a new binary file is constructed into which this PEPS is written.
    A PEPS is represented in a binary file by:
    {BCSize, BC[0], ..., BC[BCSize-1], N[0], N[1], d[0], ..., d[N[0]*N[1]-1], D, Tensors[0], ...,
     Tensors[N[0]*N[1]-1]}
    where each Tensor is represented by:
    {rank, Shape[0], ..., Shape[rank-1], size, Elements[0], ..., Elements[size-1]}.
    Note that, as always, d is stored in Fortran's column-major order.
    \param FileName input: const string&, the name for the new binary file to which this PEPS is written
    \sa void read(const string& FileName) */
  void write(const string& FileName) const;

/// Reads PEPS from binary file.
/** Given a binary file called FileName, this PEPS is replaced by the PEPS in FileName.
    A PEPS is represented in a binary file by:
    {BCSize, BC[0], ..., BC[BCSize-1], N[0], N[1], d[0], ..., d[N[0]*N[1]-1], D, Tensors[0], ...,
     Tensors[N[0]*N[1]-1]}
    where each Tensor is represented by:
    {rank, Shape[0], ..., Shape[rank-1], size, Elements[0], ..., Elements[size-1]}.
    Note that, as always, d is stored in Fortran's column-major order.
    \param FileName input: const string&, the binary file from which this PEPS is read
    \sa void write(const string& FileName) const */
  void read(const string& FileName);

/// Reads ValentinPEPS from text file.
/** Given a text file called FileName, this PEPS is replaced by the ValentinPEPS in FileName.
    The ValentinPEPS in FileName is initially generated by Matlab with writeValentinPEPS.m. The latter
    stores all tensors in Fortran's column-major order and all tensor elements in consecutive order.
    BC, N, d, and D of the ValentinPEPS in FileName and of PEPS0 must be the same.
    \param FileName input: const string&, the text file from which the ValentinPEPS is read */
  void readValentinPEPS(const string& FileName);

/// Converts float PEPS to complex<float>.
/** This function converts this PEPS of type float to complex<float>. This PEPS must fulfill T==float.
    \param PEPS0 output: PEPS< complex<float> >&, the complex equivalent to this PEPS */
  void convertToComplex(PEPS< complex<float> >& PEPS0) const;

/// Converts double PEPS to complex<double>.
/** This function converts this PEPS of type double to complex<double>. This PEPS must fulfill T==double.
    \param PEPS0 output: PEPS< complex<double> >&, the complex equivalent to this PEPS */
  void convertToComplex(PEPS< complex<double> >& PEPS0) const;

/// Returns this PEPS as vector.
/** This PEPS is returned as a ket vector in the standard basis.
    \param Vector0 output: vector<T>&, this PEPS as a vector, must fulfill Vector0.size()==d^(N[0]*N[1]) */
  void getVector(vector<T>& Vector0) const;

/// This PEPS is filled with random entries.
/** This PEPS is filled with uniformly distributed random entries from [-1,1] which are multiplied by element.
    LAPACK's XLARNV is used.
    \param Seed input: const Matrix<unsigned int>&, contains the seed values for each tensor,
                       must fulfill ((Seed.getDim0() == this->N[0]) && (Seed.getDim1() == this->N[1]))
    \param element optional input: T, multiplies all random numbers
    \sa void PEPS<T>::fillZeroes() */
  void fillRandomly(const Matrix<unsigned int>& Seed, T element = 1.0);

/// Fills this PEPS with zeroes.
/** This PEPS is filled with zeroes.
    \sa void PEPS<T>::fillRandomly(const Matrix<unsigned int>& Seed) */
  void fillZeroes();

/// This PEPS is initialized as separable many-qubit state.
/** This PEPS is initialized as a separable many-qubit state |psi>=|phi>^{\otimes N} where
    |phi>=Coefficients[0]|0>+Coefficients[1]|1>+...+Coefficients[d-1]|d-1>.
    The optional argument element allows to replace all zeroes in this PEPS by uniformly distributed
    random numbers from [-1,1] which are multiplied by element.
    \param Coefficients input: const vector<T>&, the coefficients in the standard basis
    \param element optional input: T, if given all zeroes in this PEPS are replaced by uniformly distributed
                                   random numbers from [-1,1] multiplied by element */
  void setSeparable(const vector<T>& Coefficients, T element = 0.0);

/// This PEPS is initialized as Concatenated PEPS.
/** This PEPS is initialized as a Concatenated PEPS. Each tensor is replaced by a concatenation of tensors,
    consisting of the initial tensor and M^{2}-1 auxiliary tensors having physical dimension 1. Each tensor
    is replaced by a tensor network resembling a PEPS.
    M denotes the concatenation level, such that M==1 denotes this PEPS with no auxiliary tensors, M==2
    denotes this PEPS with 3 auxiliary tensors per physical tensor, M==3 denotes this PEPS with 8 auxiliary
    tensors per physical tensor, and so on. The auxiliary tensors are chosen as Deltas to guarantee that
    the resulting Concatenated PEPS is equivalent to this PEPS for all concatenation levels.
    The optional argument element allows to replace all zeroes in the auxiliary tensors by uniformly distributed
    random numbers from [-1,1] which are multiplied by element.
    \param M input: unsigned int, the concatenation level, must be > 0
    \param element optional input: T, if given all zeroes in the auxiliary tensors are replaced by uniformly
                                   distributed random numbers from [-1,1] multiplied by element */
  void setConcatenated(unsigned int M, T element = 0.0);

/// Multiplies this PEPS with element.
/** This PEPS is multiplied with element, by multiplying the first tensor with element.
    \param element input: T, the scalar with which this PEPS is multiplied */
  void multiply(T element);

/// Decanonicalizes this PEPS.
/** This PEPS in canonical form with vertical lambda matrices LambdasV and horizontal lambda matrices LambdasH is
    decanonicalized. The standard form is obtained by contracting the square root of each lambda matrix with the
    respective neighbouring two tensors of this PEPS.
    \param LambdasV input: const Matrix< Matrix<T> >&, the lambda matrices on the vertical bonds,
                           must have the correct form
    \param LambdasH input: const Matrix< Matrix<T> >&, the lambda matrices on the horizontal bonds,
                           must have the correct form */
  void decanonicalize(const Matrix< Matrix<T> >& LambdasV, const Matrix< Matrix<T> >& LambdasH);

/// Normalizes this PEPS.
/** This function normalizes this PEPS by normalizing each tensor with
       T Tensor<T>::normalize()   .
    \return T, the product of the absolute values of the largest elements of all tensors
    \sa T Tensor<T>::normalize() */
  T normalize();

/// Returns MPS of PEPS sandwich.
/** This function returns the MPS MPS0 representing a boundary of a PEPS sandwich with open boundary
    conditions. Hereby the PEPS sandwich results from the scalar product of this PEPS as bra with itself
    as ket and thus it represents the norm of this PEPS. Direction refers to the direction in which the
    physical indices of the resulting MPS point and hence specifies which boundary of the PEPS sandwich
    is considered. The tensors are taken clockwise.
    \param Direction input: const string&, the direction of the physical bond of MPS0, must be "right",
                            "down", "left", or "up"
    \param MPS0 output: MPS<T>&, the resulting MPS representing a boundary of the PEPS sandwich
    \sa void getMPO(const string& Direction, unsigned int position, MPO<T>& MPO0) const */
  void getMPS(const string& Direction, MPS<T>& MPS0) const;

/// Returns boundary-MPS of PEPS sandwich.
/** This function approximates a boundary of a PEPS sandwich with open boundary conditions by an MPS MPS0
    with a given bond dimension. Hereby the PEPS sandwich results from the scalar product of this PEPS as
    bra with itself as ket and thus it represents the norm of this PEPS.
    Boundary specifies which boundary of the PEPS sandwich is considered. The tensors are taken
    clockwise.
    The simplified approximation error
       error := <MPS0|MPS0>-2Re(<MPS0|ExactBoundaryMPS>)
    is evaluated after each tensor-update, and the function stops if the relative change in the error is
    below eps, i.e. in tensor-update n: eps > abs(error(n)-error(n-1))/abs(error(n)).
    The function always stops if the number of sweeps exceeds maxNumSweeps.
    On input, if (MPS0.D >= this->D^2) then the exact boundary-MPS with bond dimension MPS0.D
    is returned by filling in zeroes.
    On output, epsAchieved is the achieved relative approximation error after the last tensor-update,
    numSweepsDone is the number of complete sweeps done excluding the last partial sweep, and MPS0 is the
    approximating MPS.
    If MPS0.D==1 then we normalize each tensor after its update with
       T Tensor<T>::frobeniusNormalize()   ,
    which guarantees that an initial positive MPS0 remains positive during the updates.
    If MPS0.D!=1 then each tensor is QR/LQ-decomposed after its update.
    \param Boundary input: const string&, the boundary, must be "left", "top", "right", or "bottom"
    \param eps input: double, the convergence precision, i.e. convergence in tensor-update n if
                      eps > abs(error(n)-error(n-1))/abs(error(n))
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param epsAchieved output: double&, the achieved convergence precision
    \param numSweepsDone output: unsigned int&, the number of complete sweeps done excluding the last partial sweep
    \param MPS0 input/output: MPS<T>&, on input the initial approximating MPS, on output the resulting
                              MPS approximating a boundary of the PEPS sandwich,
                              must have the correct form
    \sa void getMPS(const string& Direction, MPS<T>& MPS0) const */
  void getBoundaryMPS(const string& Boundary, double eps, unsigned int maxNumSweeps,
                      double& epsAchieved, unsigned int& numSweepsDone,
                      MPS<T>& MPS0) const;

/// Returns boundary MPO of PEPS in single-layer picture.
/** This function returns the MPO MPO0 representing a boundary of this PEPS in the single-layer picture.
    Direction refers to the direction in which the physical indices of the resulting MPO point and hence specifies
    which boundary of this PEPS is considered. The tensors are taken clockwise.
    \param Direction input: const string&, the direction of the physical bond of MPO0, must be "right",
                            "down", "left", or "up"
    \param MPO0 output: MPO<T>&, the resulting MPO representing a boundary of this PEPS
    \sa void getMPOSL(const string& Direction, unsigned int position, MPO<T>& MPO0) const */
  void getMPOSL(const string& Direction, MPO<T>& MPO0) const;

/// Returns MPO of PEPS sandwich.
/** This function returns the MPO MPO0 representing a row or column of a PEPS sandwich. Hereby the PEPS
    sandwich results from the scalar product of this PEPS as bra with itself as ket and thus it
    represents the norm of this PEPS. Direction refers to the direction in which the physical indices
    of the resulting MPO point and position specifies the row or column number of the PEPS sandwich. The
    tensors are taken clockwise.
    \param Direction input: const string&, the direction of the physical bonds of MPO0, must be "right",
                            "down", "left", or "up"
    \param position input: unsigned int, the row or column number of the PEPS sandwich
    \param MPO0 output: MPO<T>&, the resulting MPO representing a row or column of the PEPS sandwich
    \sa void getMPS(const string& Direction, MPS<T>& MPS0) const */
  void getMPO(const string& Direction, unsigned int position, MPO<T>& MPO0) const;

/// Returns MPO of PEPS in single-layer picture.
/** This function returns the MPO MPO0 representing a row or column of this PEPS in the single-layer picture.
    Direction refers to the direction in which the physical indices of the resulting MPO point and position
    specifies the row or column number of this PEPS. The tensors are taken clockwise.
    \param Direction input: const string&, the direction of the physical bonds of MPO0, must be "right",
                            "down", "left", or "up"
    \param position input: unsigned int, the row or column number of this PEPS
    \param MPO0 output: MPO<T>&, the resulting MPO representing a row or column of this PEPS
    \sa void getMPOSL(const string& Direction, MPO<T>& MPO0) const */
  void getMPOSL(const string& Direction, unsigned int position, MPO<T>& MPO0) const;

/// Multiplies MPS with bulk-MPO of PEPS sandwich and approximates result by MPS.
/** This function approximates the product of the input MPS MPS0 with a bulk-MPO of the PEPS sandwich by
    a MPS MPS1 with a given bond dimension.
    Hereby the PEPS sandwich results from the scalar product of this PEPS as bra with itself as ket and
    thus it represents the norm of this PEPS. Direction refers to the direction in which the physical indices
    of the resulting MPO point and position specifies the row or column number of the PEPS sandwich. The
    tensors are taken clockwise.
    The simplified approximation error
       error := <MPS1|MPS1>-2Re(<MPS1|BulkMPO|MPS0>)
    is evaluated after each tensor-update, and the function stops if the relative change in the error is
    below eps, i.e. in tensor-update n: eps > abs(error(n)-error(n-1))/abs(error(n)).
    The function always stops if the number of sweeps exceeds maxNumSweeps.
    On input, if (MPS1.D >= this->D^2*MPS0.D) then the exact product, as a MPS in normal shape, with bond
    dimension MPS1.D is returned by filling in zeroes.
    On output, epsAchieved is the achieved relative approximation error after the last tensor-update,
    numSweepsDone is the number of complete sweeps done excluding the last partial sweep, and MPS1 is the
    approximating MPS.
    If MPS1.D==1 then we normalize each tensor after its update with
       T Tensor<T>::frobeniusNormalize()   ,
    which guarantees that an initial positive MPS1 remains positive during the updates if MPS0 is positive.
    If MPS1.D!=1 then each tensor is QR/LQ-decomposed after its update.
    \param Direction input: const string&, the direction of the physical bonds of the bulk-MPO,
                            must be "right", "down", "left", or "up"
    \param position input: unsigned int, the row or column number of the PEPS sandwich
    \param MPS0 input: const MPS<T>&, the MPS on which the bulk-MPO acts,
                       must have the correct form
    \param eps input: double, the convergence precision, i.e. convergence in tensor-update n if
                      eps > abs(error(n)-error(n-1))/abs(error(n))
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param epsAchieved output: double&, the achieved convergence precision
    \param numSweepsDone output: unsigned int&, the number of complete sweeps done excluding the last partial sweep
    \param MPS1 input/output: MPS<T>&, on input the initial approximating MPS, on output the resulting MPS,
                              must have the correct form
    \sa void getBoundaryMPS(const string& Boundary, double eps, unsigned int maxNumSweeps, double& epsAchieved,
                            unsigned int& numSweepsDone, MPS<T>& MPS0) const */
  void multiplyBulkMPOMPS(const string& Direction, unsigned int position, const MPS<T>& MPS0,
                          double eps, unsigned int maxNumSweeps,
                          double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const;

/// Contracts MPS in reverse order with bulk-MPO applied to another MPS.
/** This function computes the scalar product of MPS0 with a bulk-MPO applied to MPS1, taking the tensors
    of MPS0 in reverse order and not complex conjugating them.
    The bulk-MPO is taken from the norm-sandwich <thisPEPS|thisPEPS>.
    If Direction=="horizontal" then the bulk-MPO is taken from column position, and MPS0 is the left MPS with
    physical indices pointing right and MPS1 is the right MPS with physical indices pointing left.
    If Direction=="vertical" then the bulk-MPO is taken from row position, and MPS0 is the lower MPS with
    physical indices pointing up and MPS1 is the upper MPS with physical indices pointing down.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their
    physical indices.
    \param Direction input: const string&, the direction of the physical bonds of the bulk-MPO,
                            must be "horizontal" or "vertical"
    \param MPS0 input: const MPS<T>&, the bra MPS that is taken in reverse order and
                       that is not complex conjugated
    \param position input: unsigned int, the column-position if Direction=="horizontal" or
                           the row-position if Direction=="vertical"
    \param MPS1 input: const MPS<T>&, the ket MPS
    \return T, the resulting contraction-value */
  T contractReverseBulkMPOMPS(const string& Direction,
                              const MPS<T>& MPS0, unsigned int position, const MPS<T>& MPS1) const;

/// Returns MPS of PEPS sandwich.
/** This function returns the MPS MPS0 representing a boundary of a PEPS sandwich with open boundary
    conditions. Hereby the PEPS sandwich results from the scalar product of this PEPS as bra with PEPS0
    as ket. Direction refers to the direction in which the physical indices of the resulting MPS point
    and hence specifies which boundary of the PEPS sandwich is considered. The tensors are taken clockwise.
    \param PEPS0 input: const PEPS<T>&, the ket PEPS of the PEPS sandwich
    \param Direction input: const string&, the direction of the physical bond of MPS0, must be "right",
                            "down", "left", or "up"
    \param MPS0 output: MPS<T>&, the resulting MPS representing a boundary of the PEPS sandwich
    \sa void getMPO(const PEPS<T>& PEPS0, const string& Direction, unsigned int position, MPO<T>& MPO0) const */
  void getMPS(const PEPS<T>& PEPS0, const string& Direction, MPS<T>& MPS0) const;

/// Returns boundary-MPS of PEPS sandwich.
/** This function approximates a boundary of a PEPS sandwich with open boundary conditions by an MPS MPS0
    with a given bond dimension. Hereby the PEPS sandwich results from the scalar product of this PEPS as
    bra with PEPS0 as ket.
    Boundary specifies which boundary of the PEPS sandwich is considered. The tensors are taken
    clockwise.
    The simplified approximation error
       error := <MPS0|MPS0>-2Re(<MPS0|ExactBoundaryMPS>)
    is evaluated after each tensor-update, and the function stops if the relative change in the error is
    below eps, i.e. in tensor-update n: eps > abs(error(n)-error(n-1))/abs(error(n)).
    The function always stops if the number of sweeps exceeds maxNumSweeps.
    On input, if (MPS0.D >= PEPS0.D*this->D) then the exact boundary-MPS with bond dimension MPS0.D
    is returned by filling in zeroes.
    On output, epsAchieved is the achieved relative approximation error after the last tensor-update,
    numSweepsDone is the number of complete sweeps done excluding the last partial sweep, and MPS0 is the
    approximating MPS.
    If MPS0.D==1 then we normalize each tensor after its update with
       T Tensor<T>::frobeniusNormalize()   ,
    while if MPS0.D!=1 then each tensor is QR/LQ-decomposed after its update.
    \param PEPS0 input: const PEPS<T>&, the ket PEPS of the PEPS sandwich
    \param Boundary input: const string&, the boundary, must be "left", "top", "right", or "bottom"
    \param eps input: double, the convergence precision, i.e. convergence in tensor-update n if
                      eps > abs(error(n)-error(n-1))/abs(error(n))
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param epsAchieved output: double&, the achieved convergence precision
    \param numSweepsDone output: unsigned int&, the number of complete sweeps done excluding the last partial sweep
    \param MPS0 input/output: MPS<T>&, on input the initial approximating MPS, on output the resulting
                              MPS approximating a boundary of the PEPS sandwich,
                              must have the correct form
    \sa void getMPS(const PEPS<T>& PEPS0, const string& Direction, MPS<T>& MPS0) const */
  void getBoundaryMPS(const PEPS<T>& PEPS0, const string& Boundary, double eps, unsigned int maxNumSweeps,
                      double& epsAchieved, unsigned int& numSweepsDone,
                      MPS<T>& MPS0) const;

/// Returns MPO of PEPS sandwich.
/** This function returns the MPO MPO0 representing a row or column of a PEPS sandwich. Hereby the PEPS
    sandwich results from the scalar product of this PEPS as bra with PEPS0 as ket. Direction refers to
    the direction in which the physical indices of the resulting MPO point and position specifies the
    row or column number of the PEPS sandwich. The tensors are taken clockwise.
    \param PEPS0 input: const PEPS<T>&, the ket PEPS of the PEPS sandwich
    \param Direction input: const string&, the direction of the physical bonds of MPO0, must be "right",
                            "down", "left", or "up"
    \param position input: unsigned int, the row or column number of the PEPS sandwich
    \param MPO0 output: MPO<T>&, the resulting MPO representing a row or column of the PEPS sandwich
    \sa void getMPS(const PEPS<T>& PEPS0, const string& Direction, MPS<T>& MPS0) const */
  void getMPO(const PEPS<T>& PEPS0, const string& Direction, unsigned int position, MPO<T>& MPO0) const;

/// Multiplies MPS with bulk-MPO of PEPS sandwich and approximates result by MPS.
/** This function approximates the product of the input MPS MPS0 with a bulk-MPO of the PEPS sandwich by
    a MPS MPS1 with a given bond dimension. Hereby the PEPS sandwich results from the scalar product of
    this PEPS as bra with PEPS0 as ket.
    Direction refers to the direction in which the physical indices of the resulting MPO point and position
    specifies the row or column number of the PEPS sandwich. The tensors are taken clockwise.
    The simplified approximation error
       error := <MPS1|MPS1>-2Re(<MPS1|BulkMPO|MPS0>)
    is evaluated after each tensor-update, and the function stops if the relative change in the error is
    below eps, i.e. in tensor-update n: eps > abs(error(n)-error(n-1))/abs(error(n)).
    The function always stops if the number of sweeps exceeds maxNumSweeps.
    On input, if (MPS1.D >= PEPS0.D*this->D*MPS0.D) then the exact product, as a MPS in normal shape, with bond
    dimension MPS1.D is returned by filling in zeroes.
    On output, epsAchieved is the achieved relative approximation error after the last tensor-update,
    numSweepsDone is the number of complete sweeps done excluding the last partial sweep, and MPS1 is the
    approximating MPS.
    If MPS1.D==1 then we normalize each tensor after its update with
       T Tensor<T>::frobeniusNormalize()   ,
    while if MPS1.D!=1 then each tensor is QR/LQ-decomposed after its update.
    \param PEPS0 input: const PEPS<T>&, the ket PEPS of the PEPS sandwich
    \param Direction input: const string&, the direction of the physical bonds of the bulk-MPO,
                            must be "right", "down", "left", or "up"
    \param position input: unsigned int, the row or column number of the PEPS sandwich
    \param MPS0 input: const MPS<T>&, the MPS on which the bulk-MPO acts,
                       must have the correct form
    \param eps input: double, the convergence precision, i.e. convergence in tensor-update n if
                      eps > abs(error(n)-error(n-1))/abs(error(n))
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param epsAchieved output: double&, the achieved convergence precision
    \param numSweepsDone output: unsigned int&, the number of complete sweeps done excluding the last partial sweep
    \param MPS1 input/output: MPS<T>&, on input the initial approximating MPS, on output the resulting MPS,
                              must have the correct form
    \sa void getBoundaryMPS(const PEPS<T>& PEPS0, const string& Boundary, double eps, unsigned int maxNumSweeps,
                            double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const */
  void multiplyBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, unsigned int position, const MPS<T>& MPS0,
                          double eps, unsigned int maxNumSweeps,
                          double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const;

/// Contracts MPS in reverse order with bulk-MPO applied to another MPS.
/** This function computes the scalar product of MPS0 with a bulk-MPO applied to MPS1, taking the tensors
    of MPS0 in reverse order and not complex conjugating them.
    The bulk-MPO is taken from the sandwich <thisPEPS|PEPS0>.
    If Direction=="horizontal" then the bulk-MPO is taken from column position, and MPS0 is the left MPS with
    physical indices pointing right and MPS1 is the right MPS with physical indices pointing left.
    If Direction=="vertical" then the bulk-MPO is taken from row position, and MPS0 is the lower MPS with
    physical indices pointing up and MPS1 is the upper MPS with physical indices pointing down.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their
    physical indices.
    \param PEPS0 input: const PEPS<T>&, the ket PEPS of the PEPS sandwich
    \param Direction input: const string&, the direction of the physical bonds of the bulk-MPO,
                            must be "horizontal" or "vertical"
    \param MPS0 input: const MPS<T>&, the bra MPS that is taken in reverse order and
                       that is not complex conjugated
    \param position input: unsigned int, the column-position if Direction=="horizontal" or
                           the row-position if Direction=="vertical"
    \param MPS1 input: const MPS<T>&, the ket MPS
    \return T, the resulting contraction-value */
  T contractReverseBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction,
                              const MPS<T>& MPS0, unsigned int position, const MPS<T>& MPS1) const;

/// Returns MPS of PEPS sandwich.
/** This function returns the MPS MPS0 representing a boundary of a PEPS sandwich with open boundary
    conditions. Hereby the PEPS sandwich results from the expectation value of PEPO0 with this PEPS.
    Direction refers to the direction in which the physical indices of the resulting MPS point
    and hence specifies which boundary of the PEPS sandwich is considered. The tensors are taken clockwise.
    \param PEPO0 input: const PEPO<T>&, the PEPO of the PEPS sandwich
    \param Direction input: const string&, the direction of the physical bond of MPS0, must be "right",
                            "down", "left", or "up"
    \param MPS0 output: MPS<T>&, the resulting MPS representing a boundary of the PEPS sandwich
    \sa void getMPO(const PEPO<T>& PEPO0, const string& Direction, unsigned int position, MPO<T>& MPO0) const */
  void getMPS(const PEPO<T>& PEPO0, const string& Direction, MPS<T>& MPS0) const;

/// Returns boundary-MPS of PEPS sandwich.
/** This function approximates a boundary of a PEPS sandwich with open boundary conditions by an MPS MPS0
    with a given bond dimension. Hereby the PEPS sandwich results from the expectation value of PEPO0 with this PEPS.
    Boundary specifies which boundary of the PEPS sandwich is considered. The tensors are taken
    clockwise.
    The simplified approximation error
       error := <MPS0|MPS0>-2Re(<MPS0|ExactBoundaryMPS>)
    is evaluated after each tensor-update, and the function stops if the relative change in the error is
    below eps, i.e. in tensor-update n: eps > abs(error(n)-error(n-1))/abs(error(n)).
    The function always stops if the number of sweeps exceeds maxNumSweeps.
    On input, if (MPS0.D >= this->D*PEPO0.D*this->D) then the exact boundary-MPS with bond dimension MPS0.D
    is returned by filling in zeroes.
    On output, epsAchieved is the achieved relative approximation error after the last tensor-update,
    numSweepsDone is the number of complete sweeps done excluding the last partial sweep, and MPS0 is the
    approximating MPS.
    If MPS0.D==1 then we normalize each tensor after its update with
       T Tensor<T>::frobeniusNormalize()   ,
    while if MPS0.D!=1 then each tensor is QR/LQ-decomposed after its update.
    \param PEPO0 input: const PEPO<T>&, the PEPO of the PEPS sandwich
    \param Boundary input: const string&, the boundary, must be "left", "top", "right", or "bottom"
    \param eps input: double, the convergence precision, i.e. convergence in tensor-update n if
                      eps > abs(error(n)-error(n-1))/abs(error(n))
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param epsAchieved output: double&, the achieved convergence precision
    \param numSweepsDone output: unsigned int&, the number of complete sweeps done excluding the last partial sweep
    \param MPS0 input/output: MPS<T>&, on input the initial approximating MPS, on output the resulting
                              MPS approximating a boundary of the PEPS sandwich,
                              must have the correct form */
  void getBoundaryMPS(const PEPO<T>& PEPO0, const string& Boundary,
                      double eps, unsigned int maxNumSweeps,
                      double& epsAchieved, unsigned int& numSweepsDone,
                      MPS<T>& MPS0) const;

/// Returns MPO of PEPS sandwich.
/** This function returns the MPO MPO0 representing a row or column of a PEPS sandwich. Hereby the PEPS
    sandwich results from the expectation value of PEPO0 with this PEPS. Direction refers to
    the direction in which the physical indices of the resulting MPO point and position specifies the
    row or column number of the PEPS sandwich. The tensors are taken clockwise.
    \param PEPO0 input: const PEPO<T>&, the PEPO of the PEPS sandwich
    \param Direction input: const string&, the direction of the physical bonds of MPO0, must be "right",
                            "down", "left", or "up"
    \param position input: unsigned int, the row or column number of the PEPS sandwich
    \param MPO0 output: MPO<T>&, the resulting MPO representing a row or column of the PEPS sandwich
    \sa void getMPS(const PEPO<T>& PEPO0, const string& Direction, MPS<T>& MPS0) const */
  void getMPO(const PEPO<T>& PEPO0, const string& Direction, unsigned int position, MPO<T>& MPO0) const;

/// Multiplies MPS with bulk-MPO of PEPS sandwich and approximates result by MPS.
/** This function approximates the product of the input MPS MPS0 with a bulk-MPO of the PEPS sandwich by
    a MPS MPS1 with a given bond dimension. Hereby the PEPS sandwich results from the expectation value of
    PEPO0 with this PEPS.
    Direction refers to the direction in which the physical indices of the resulting MPO point and position
    specifies the row or column number of the PEPS sandwich. The tensors are taken clockwise.
    The simplified approximation error
       error := <MPS1|MPS1>-2Re(<MPS1|BulkMPO|MPS0>)
    is evaluated after each tensor-update, and the function stops if the relative change in the error is
    below eps, i.e. in tensor-update n: eps > abs(error(n)-error(n-1))/abs(error(n)).
    The function always stops if the number of sweeps exceeds maxNumSweeps.
    On input, if (MPS1.D >= this->D*PEPO0.D*this->D*MPS0.D) then the exact product, as a MPS in normal shape,
    with bond dimension MPS1.D is returned by filling in zeroes.
    On output, epsAchieved is the achieved relative approximation error after the last tensor-update,
    numSweepsDone is the number of complete sweeps done excluding the last partial sweep, and MPS1 is the
    approximating MPS.
    If MPS1.D==1 then we normalize each tensor after its update with
       T Tensor<T>::frobeniusNormalize()   ,
    while if MPS1.D!=1 then each tensor is QR/LQ-decomposed after its update.
    \param PEPO0 input: const PEPO<T>&, the PEPO of the PEPS sandwich
    \param Direction input: const string&, the direction of the physical bonds of the bulk-MPO,
                            must be "right", "down", "left", or "up"
    \param position input: unsigned int, the row or column number of the PEPS sandwich
    \param MPS0 input: const MPS<T>&, the MPS on which the bulk-MPO acts,
                       must have the correct form
    \param eps input: double, the convergence precision, i.e. convergence in tensor-update n if
                      eps > abs(error(n)-error(n-1))/abs(error(n))
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param epsAchieved output: double&, the achieved convergence precision
    \param numSweepsDone output: unsigned int&, the number of complete sweeps done excluding the last partial sweep
    \param MPS1 input/output: MPS<T>&, on input the initial approximating MPS, on output the resulting MPS,
                              must have the correct form
    \sa void getBoundaryMPS(const PEPO<T>& PEPO0, const string& Boundary,
                            double eps, unsigned int maxNumSweeps,
                            double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const */
  void multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const string& Direction, unsigned int position, const MPS<T>& MPS0,
                          double eps, unsigned int maxNumSweeps,
                          double& epsAchieved, unsigned int& numSweepsDone,
                          MPS<T>& MPS1) const;

/// Returns boundary-MPS of PEPS sandwich.
/** This function approximates a boundary of a PEPS sandwich with open boundary conditions by an MPS MPS0
    with a given bond dimension. Hereby the PEPS sandwich results from the scalar product of this PEPS as
    bra with PEPO0 and PEPS0 as ket.
    Boundary specifies which boundary of the PEPS sandwich is considered. The tensors are taken
    clockwise.
    The simplified approximation error
       error := <MPS0|MPS0>-2Re(<MPS0|ExactBoundaryMPS>)
    is evaluated after each tensor-update, and the function stops if the relative change in the error is
    below eps, i.e. in tensor-update n: eps > abs(error(n)-error(n-1))/abs(error(n)).
    The function always stops if the number of sweeps exceeds maxNumSweeps.
    On input, if (MPS0.D >= PEPS0.D*PEPO0.D*this->D) then the exact boundary-MPS with bond dimension MPS0.D
    is returned by filling in zeroes.
    On output, epsAchieved is the achieved relative approximation error after the last tensor-update,
    numSweepsDone is the number of complete sweeps done excluding the last partial sweep, and MPS0 is the
    approximating MPS.
    If MPS0.D==1 then we normalize each tensor after its update with
       T Tensor<T>::frobeniusNormalize()   ,
    while if MPS0.D!=1 then each tensor is QR/LQ-decomposed after its update.
    \param PEPO0 input: const PEPO<T>&, the PEPO of the PEPS sandwich
    \param PEPS0 input: const PEPS<T>&, the ket PEPS of the PEPS sandwich
    \param Boundary input: const string&, the boundary, must be "left", "top", "right", or "bottom"
    \param eps input: double, the convergence precision, i.e. convergence in tensor-update n if
                      eps > abs(error(n)-error(n-1))/abs(error(n))
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param epsAchieved output: double&, the achieved convergence precision
    \param numSweepsDone output: unsigned int&, the number of complete sweeps done excluding the last partial sweep
    \param MPS0 input/output: MPS<T>&, on input the initial approximating MPS, on output the resulting
                              MPS approximating a boundary of the PEPS sandwich,
                              must have the correct form */
  void getBoundaryMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary,
                      double eps, unsigned int maxNumSweeps,
                      double& epsAchieved, unsigned int& numSweepsDone,
                      MPS<T>& MPS0) const;

/// Multiplies MPS with bulk-MPO of PEPS sandwich and approximates result by MPS.
/** This function approximates the product of the input MPS MPS0 with a bulk-MPO of the PEPS sandwich by
    a MPS MPS1 with a given bond dimension. Hereby the PEPS sandwich results from the scalar product of
    this PEPS as bra with PEPO0 and PEPS0 as ket.
    Direction refers to the direction in which the physical indices of the resulting MPO point and position
    specifies the row or column number of the PEPS sandwich. The tensors are taken clockwise.
    The simplified approximation error
       error := <MPS1|MPS1>-2Re(<MPS1|BulkMPO|MPS0>)
    is evaluated after each tensor-update, and the function stops if the relative change in the error is
    below eps, i.e. in tensor-update n: eps > abs(error(n)-error(n-1))/abs(error(n)).
    The function always stops if the number of sweeps exceeds maxNumSweeps.
    On input, if (MPS1.D >= PEPS0.D*PEPO0.D*this->D*MPS0.D) then the exact product, as a MPS in normal shape,
    with bond dimension MPS1.D is returned by filling in zeroes.
    On output, epsAchieved is the achieved relative approximation error after the last tensor-update,
    numSweepsDone is the number of complete sweeps done excluding the last partial sweep, and MPS1 is the
    approximating MPS.
    If MPS1.D==1 then we normalize each tensor after its update with
       T Tensor<T>::frobeniusNormalize()   ,
    while if MPS1.D!=1 then each tensor is QR/LQ-decomposed after its update.
    \param PEPO0 input: const PEPO<T>&, the PEPO of the PEPS sandwich
    \param PEPS0 input: const PEPS<T>&, the ket PEPS of the PEPS sandwich
    \param Direction input: const string&, the direction of the physical bonds of the bulk-MPO,
                            must be "right", "down", "left", or "up"
    \param position input: unsigned int, the row or column number of the PEPS sandwich
    \param MPS0 input: const MPS<T>&, the MPS on which the bulk-MPO acts,
                       must have the correct form
    \param eps input: double, the convergence precision, i.e. convergence in tensor-update n if
                      eps > abs(error(n)-error(n-1))/abs(error(n))
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param epsAchieved output: double&, the achieved convergence precision
    \param numSweepsDone output: unsigned int&, the number of complete sweeps done excluding the last partial sweep
    \param MPS1 input/output: MPS<T>&, on input the initial approximating MPS, on output the resulting MPS,
                              must have the correct form
    \sa void getBoundaryMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary,
                            double eps, unsigned int maxNumSweeps,
                            double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const */
  void multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0,
                          const string& Direction, unsigned int position, const MPS<T>& MPS0,
                          double eps, unsigned int maxNumSweeps,
                          double& epsAchieved, unsigned int& numSweepsDone,
                          MPS<T>& MPS1) const;

/// Takes scalar product of this PEPS with itself.
/** This function computes the scalar product of this PEPS with itself.
    If the contraction Direction is set to "horizontal", two neighbouring columns are approximated by one
    new column; whereas if Direction is "vertical", two neigbouring rows are approximated by one row. In
    the case of open BC the MPO-MPS approximation is done with the friend function
       friend void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                  unsigned int maxNumSweeps, double& errorAchieved,
                                  unsigned int& numSweepsDone, MPS<T>& MPS1)   .
    The vector D2s specifies the maximal virtual bond dimensions D2 for the approximating MPSs,
    enumerated from left to right if Direction is "horizontal" or up to down if Direction is "vertical".
    Epss defines the convergence precisions and MaxNumsSweeps the maximal numbers of sweeps allowed.
    The vectors ErrorsAchieved and NumsSweepsDone capture the respective resulting values.
    The two halves are contracted independently, such that this->N[1]-3 column-column approximations are
    done if Direction is "horizontal", or this->N[0]-3 row-row approximations are done if Direction is
    "vertical".
    \param Direction input: const string&, the contraction direction, must be "horizontal" or "vertical"
    \param D2s input: const vector<unsigned int>&, the maximal virtual bond dimensions,
                      must fulfill D2s.size()==this->N[1]-3 if Direction=="horizontal" or
                      D2s.size()==this->N[0]-3 if Direction=="vertical"
    \param Epss input: const vector<double>&, the convergence precisions,
                       must fulfill Epss.size()==this->N[1]-3 if Direction=="horizontal" or
                       Epss.size()==this->N[0]-3 if Direction=="vertical"
    \param MaxNumsSweeps input: const vector<unsigned int>&, the maximal numbers of sweeps allowed,
                                must fulfill MaxNumsSweeps.size()==this->N[1]-3 if Direction=="horizontal"
                                or MaxNumsSweeps.size()==this->N[0]-3 if Direction=="vertical"
    \param ErrorsAchieved output: vector<double>&, the achieved approximation errors,
                                  must fulfill ErrorsAchieved.size()==this->N[1]-3 if
                                  Direction=="horizontal" or
                                  ErrorsAchieved.size()==this->N[0]-3 if Direction=="vertical"
    \param NumsSweepsDone output: vector<unsigned int>&, the numbers of sweeps done,
                                  must fulfill NumsSweepsDone.size()==this->N[1]-3 if
                                  Direction=="horizontal" or
                                  NumsSweepsDone.size()==this->N[0]-3 if Direction=="vertical"
    \return T, the resulting value of the scalar product
    \sa template<class T> friend void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                                     unsigned int maxNumSweeps, double& errorAchieved,
                                                     unsigned int& numSweepsDone, MPS<T>& MPS1)
    \sa template<class T> T MPS<T>::contractReverse(const MPS<T>& MPS0) const
    \sa template<class T> T contractReverseMPOMPS(const MPS<T>& MPS0, const MPO<T>& MPO0,
                                                  const MPS<T>& MPS1) */
  T scalarProduct(const string& Direction, const vector<unsigned int>& D2s, const vector<double>& Epss,
                  const vector<unsigned int>& MaxNumsSweeps, vector<double>& ErrorsAchieved,
                  vector<unsigned int>& NumsSweepsDone) const;

/// Takes scalar product of this PEPS with itself.
/** This function computes the scalar product of this PEPS with itself, after normalizing the columns or
    rows of this PEPS sandwich.
    If the contraction Direction is set to "horizontal", two neighbouring columns are approximated by one
    new column; whereas if Direction is "vertical", two neigbouring rows are approximated by one row. In
    the case of open BC the MPO-MPS approximation is done with the friend function
       friend void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                  unsigned int maxNumSweeps, double& errorAchieved,
                                  unsigned int& numSweepsDone, MPS<T>& MPS1)   .
    The vector D2s specifies the maximal virtual bond dimensions D2 for the approximating MPSs,
    enumerated from left to right if Direction is "horizontal" or up to down if Direction is "vertical".
    Epss defines the convergence precisions and MaxNumsSweeps the maximal numbers of sweeps allowed.
    The vectors ErrorsAchieved and NumsSweepsDone capture the respective resulting values.
    Before the contraction of this PEPS sandwich, its columns or rows are normalized with
       T MPS<T>::simplifiedNormalize()   and   T MPO<T>::simplifiedNormalize()
    and their norms are returned in the vector Norms.
    The two halves are contracted independently, such that this->N[1]-3 column-column approximations are
    done if Direction is "horizontal", or this->N[0]-3 row-row approximations are done if Direction is
    "vertical".
    \param Direction input: const string&, the contraction direction, must be "horizontal" or "vertical"
    \param D2s input: const vector<unsigned int>&, the maximal virtual bond dimensions,
                      must fulfill D2s.size()==this->N[1]-3 if Direction=="horizontal" or
                      D2s.size()==this->N[0]-3 if Direction=="vertical"
    \param Epss input: const vector<double>&, the convergence precisions,
                       must fulfill Epss.size()==this->N[1]-3 if Direction=="horizontal" or
                       Epss.size()==this->N[0]-3 if Direction=="vertical"
    \param MaxNumsSweeps input: const vector<unsigned int>&, the maximal numbers of sweeps allowed,
                                must fulfill MaxNumsSweeps.size()==this->N[1]-3 if Direction=="horizontal"
                                or MaxNumsSweeps.size()==this->N[0]-3 if Direction=="vertical"
    \param ErrorsAchieved output: vector<double>&, the achieved approximation errors,
                                  must fulfill ErrorsAchieved.size()==this->N[1]-3 if
                                  Direction=="horizontal" or
                                  ErrorsAchieved.size()==this->N[0]-3 if Direction=="vertical"
    \param NumsSweepsDone output: vector<unsigned int>&, the numbers of sweeps done,
                                  must fulfill NumsSweepsDone.size()==this->N[1]-3 if
                                  Direction=="horizontal" or
                                  NumsSweepsDone.size()==this->N[0]-3 if Direction=="vertical"
    \param Norms output: vector<T>&, the norms of the columns or rows of this PEPS sandwich,
                         must fulfill Norms.size()==this->N[1] if Direction=="horizontal" or
                         Norms.size()==this->N[0] if Direction=="vertical"
    \return T, the resulting value of the scalar product excluding the norms
    \sa template<class T> T MPS<T>::simplifiedNormalize()
    \sa template<class T> T MPO<T>::simplifiedNormalize()
    \sa template<class T> friend void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                                     unsigned int maxNumSweeps, double& errorAchieved,
                                                     unsigned int& numSweepsDone, MPS<T>& MPS1)
    \sa template<class T> T MPS<T>::contractReverse(const MPS<T>& MPS0) const
    \sa template<class T> T contractReverseMPOMPS(const MPS<T>& MPS0, const MPO<T>& MPO0,
                                                  const MPS<T>& MPS1) */
  T scalarProduct(const string& Direction, const vector<unsigned int>& D2s, const vector<double>& Epss,
                  const vector<unsigned int>& MaxNumsSweeps, vector<double>& ErrorsAchieved,
                  vector<unsigned int>& NumsSweepsDone, vector<T>& Norms) const;

/// Takes scalar product of this PEPS with itself.
/** This function computes the scalar product of this PEPS with itself, i.e. the PEPS sandwich
       <thisPEPS|thisPEPS>   .
    If the contraction Direction is set to "horizontal" then two neighbouring columns are approximated by one
    new column, whereas if Direction is "vertical" then two neigbouring rows are approximated by one row.
    D2 specifies the maximal virtual bond dimension for the approximating boundary-MPSs, eps defines the
    convergence precision and maxNumSweeps the maximal number of sweeps allowed.
    The vectors EpssAchieved and NumsSweepsDone capture the respective resulting values.
    All boundary-MPSs are normalized and their norms are returned in the vector Norms.
    \param Direction input: const string&, the contraction direction, must be "horizontal" or "vertical"
    \param D2 input: unsigned int, the maximal virtual bond dimension, must be > 0
    \param eps input: double, the convergence precision
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed, must be > 0
    \param EpssAchieved output: vector<double>&, the achieved convergence precisions,
                                must fulfill EpssAchieved.size()==this->N[1]-1 if
                                Direction=="horizontal" or
                                EpssAchieved.size()==this->N[0]-1 if Direction=="vertical"
    \param NumsSweepsDone output: vector<unsigned int>&, the numbers of sweeps done,
                                  must fulfill NumsSweepsDone.size()==this->N[1]-1 if
                                  Direction=="horizontal" or
                                  NumsSweepsDone.size()==this->N[0]-1 if Direction=="vertical"
    \param Norms output: vector<T>&, the norms of the boundary-MPSs of the PEPS sandwich,
                         must fulfill Norms.size()==this->N[1]-1 if Direction=="horizontal" or
                         Norms.size()==this->N[0]-1 if Direction=="vertical"
    \return T, the resulting value of the scalar product excluding the norms,
            has to be multiplied by the product of all values in Norms to get the true scalar product
    \sa template<class T> T MPS<T>::normalize()
    \sa template<class T> void PEPS<T>::getBoundaryMPS(const string& Boundary,
                                                       double eps, unsigned int maxNumSweeps,
                                                       double& epsAchieved, unsigned int& numSweepsDone,
                                                       MPS<T>& MPS0) const
    \sa template<class T> void PEPS<T>::multiplyBulkMPOMPS(const string& Direction,
                                                           unsigned int position, const MPS<T>& MPS0,
                                                           double eps, unsigned int maxNumSweeps,
                                                           double& epsAchieved, unsigned int& numSweepsDone,
                                                           MPS<T>& MPS1) const
    \sa template<class T> T PEPS<T>::normalize(const string& Direction,
                                               unsigned int D2, double eps, unsigned int maxNumSweeps)
    \sa template<class T> T PEPS<T>::expectationValue(const vector< vector<unsigned int> >& PositionsRow,
                                                      const vector< vector<unsigned int> >& PositionsCol,
                                                      const vector< vector< Matrix<T> > >& Interactions,
                                                      const string& Direction,
                                                      unsigned int D2, double eps, unsigned int maxNumSweeps) const */
  T scalarProduct(const string& Direction,
                  unsigned int D2, double eps, unsigned int maxNumSweeps,
                  vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone,
                  vector<T>& Norms) const;

/// Takes scalar product of this PEPS with itself in single-layer picture.
/** This function computes the scalar product of this PEPS with itself in the single-layer picture.
    If the contraction Direction is set to "horizontal", two neighbouring columns are approximated by one
    new column; whereas if Direction is "vertical", two neigbouring rows are approximated by one row.
    The MPO-MPO approximation is done with the function
       void MPO<T>::multiply(const MPO<T>& MPO0, double eps, unsigned int maxNumSweeps,
                             double& errorAchieved, unsigned int& numSweepsDone, MPO<T>& MPO1) const   .
    The vector D2s specifies the maximal virtual bond dimensions D2 for the approximating MPOs,
    enumerated from left to right if Direction is "horizontal" or up to down if Direction is "vertical".
    Epss defines the convergence precisions and MaxNumsSweeps the maximal numbers of sweeps allowed.
    The vectors ErrorsAchieved and NumsSweepsDone capture the respective resulting values.
    d2s defines the maximal physical bond dimensions enforced with
       void MPO<T>::setd2(unsigned int d2, double eps, vector< vector<double> >& Eigenvalues)   .
    The two halves are contracted independently, such that this->N[1]-2 column-column approximations are
    done if Direction is "horizontal", or this->N[0]-2 row-row approximations are done if Direction is
    "vertical".
    \param Direction input: const string&, the contraction direction, must be "horizontal" or "vertical"
    \param D2s input: const vector<unsigned int>&, the maximal virtual bond dimensions,
                      must fulfill D2s.size()==this->N[1]-2 if Direction=="horizontal" or
                      D2s.size()==this->N[0]-2 if Direction=="vertical"
    \param Epss input: const vector<double>&, the convergence precisions,
                       must fulfill Epss.size()==this->N[1]-2 if Direction=="horizontal" or
                       Epss.size()==this->N[0]-2 if Direction=="vertical"
    \param MaxNumsSweeps input: const vector<unsigned int>&, the maximal numbers of sweeps allowed,
                                must fulfill MaxNumsSweeps.size()==this->N[1]-2 if Direction=="horizontal"
                                or MaxNumsSweeps.size()==this->N[0]-2 if Direction=="vertical"
    \param ErrorsAchieved output: vector<double>&, the achieved approximation errors,
                                  must fulfill ErrorsAchieved.size()==this->N[1]-2 if
                                  Direction=="horizontal" or
                                  ErrorsAchieved.size()==this->N[0]-2 if Direction=="vertical"
    \param NumsSweepsDone output: vector<unsigned int>&, the numbers of sweeps done,
                                  must fulfill NumsSweepsDone.size()==this->N[1]-2 if
                                  Direction=="horizontal" or
                                  NumsSweepsDone.size()==this->N[0]-2 if Direction=="vertical"
    \param d2s input: const vector<unsigned int>&, the maximal physical bond dimensions,
                      must fulfill d2s.size()==this->N[1]-2 if Direction=="horizontal" or
                      d2s.size()==this->N[0]-2 if Direction=="vertical"
    \param eps input: double, the cutoff for the pseudoinverse, elements <= eps are assumed to be zero
    \param Eigenvalues output: Matrix< vector<double> >&, the eigenvalues returned from MPO<T>::setd2,
                               must fulfill (Eigenvalues.getDim0()==this->N[0] &&
                               Eigenvalues.getDim1()==this->N[1]-2) if Direction=="horizontal" or
                               (Eigenvalues.getDim0()==this->N[0]-2 && Eigenvalues.getDim1()==this->N[1])
                               if Direction=="vertical"
    \return T, the resulting value of the scalar product
    \sa template<class T> void MPO<T>::multiply(const MPO<T>& MPO0, double eps, unsigned int maxNumSweeps,
                                                double& errorAchieved, unsigned int& numSweepsDone,
                                                MPO<T>& MPO1) const
    \sa template<class T> void MPO<T>::setd2(unsigned int d2, double eps,
                                             vector< vector<double> >& Eigenvalues)
    \sa template<class T> T MPS<T>::contractReverse(const MPS<T>& MPS0) const */
  T scalarProduct(const string& Direction, const vector<unsigned int>& D2s, const vector<double>& Epss,
                  const vector<unsigned int>& MaxNumsSweeps, vector<double>& ErrorsAchieved,
                  vector<unsigned int>& NumsSweepsDone, const vector<unsigned int>& d2s,
                  double eps, Matrix< vector<double> >& Eigenvalues) const;

/// Takes scalar product of this PEPS with itself in single-layer picture.
/** This function computes the scalar product of this PEPS with itself in the single-layer picture.
    If the contraction Direction is set to "horizontal", two neighbouring columns are approximated by one
    new column; whereas if Direction is "vertical", two neigbouring rows are approximated by one row.
    The MPO-MPO approximation is done with the function
       void MPO<T>::multiply(const MPO<T>& MPO0, double eps, unsigned int maxNumSweeps,
                             double& errorAchieved, unsigned int& numSweepsDone, MPO<T>& MPO1) const
    and the purification bond is set with
       void MPO<T>::setd2(unsigned int d2, double cutoff, double eps, unsigned int maxNumSweeps,
                          double& errorAchieved, unsigned int& numSweepsDone)   .
    D2 specifies the maximal virtual bond dimension and d2 the maximal purification bond dimension for the
    approximating MPOs. eps defines the convergence precision, maxNumSweeps the maximal number of sweeps allowed,
    and cutoff the cutoff for the pseudoinverse.
    The vectors ErrorsAchieved1 and NumsSweepsDone1 capture the resulting values from the MPO-MPO approximation,
    enumerated from left to right if Direction is "horizontal" or up to down if Direction is "vertical".
    ErrorsAchieved2 and NumsSweepsDone2 capture the resulting values from setd2.
    The two halves are contracted independently, such that this->N[1]-2 column-column approximations are
    done if Direction is "horizontal", or this->N[0]-2 row-row approximations are done if Direction is
    "vertical".
    \param Direction input: const string&, the contraction direction, must be "horizontal" or "vertical"
    \param D2 input: unsigned int, the maximal virtual bond dimension
    \param d2 input: unsigned int, the maximal purification bond dimension
    \param eps input: double, the convergence precision
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param cutoff input: double, the cutoff for the pseudoinverse
    \param ErrorsAchieved1 output: vector<double>&, the achieved approximation errors from the MPO-MPO approximation,
                                   must fulfill ErrorsAchieved1.size()==this->N[1]-2 if Direction=="horizontal" or
                                   ErrorsAchieved1.size()==this->N[0]-2 if Direction=="vertical"
    \param NumsSweepsDone1 output: vector<unsigned int>&, the numbers of sweeps done from the MPO-MPO approximation,
                                   must fulfill NumsSweepsDone1.size()==this->N[1]-2 if Direction=="horizontal" or
                                   NumsSweepsDone1.size()==this->N[0]-2 if Direction=="vertical"
    \param ErrorsAchieved2 output: vector<double>&, the achieved approximation errors from setd2,
                                   must fulfill ErrorsAchieved2.size()==this->N[1]-2 if Direction=="horizontal" or
                                   ErrorsAchieved2.size()==this->N[0]-2 if Direction=="vertical"
    \param NumsSweepsDone2 output: vector<unsigned int>&, the numbers of sweeps done from setd2,
                                   must fulfill NumsSweepsDone2.size()==this->N[1]-2 if Direction=="horizontal" or
                                   NumsSweepsDone2.size()==this->N[0]-2 if Direction=="vertical"
    \return T, the resulting value of the scalar product
    \sa template<class T> void MPO<T>::multiply(const MPO<T>& MPO0, double eps, unsigned int maxNumSweeps,
                                                double& errorAchieved, unsigned int& numSweepsDone,
                                                MPO<T>& MPO1) const
    \sa template<class T> void MPO<T>::setd2(unsigned int d2, double cutoff, double eps, unsigned int maxNumSweeps,
                                             double& errorAchieved, unsigned int& numSweepsDone)
    \sa template<class T> T MPS<T>::contractReverse(const MPS<T>& MPS0) const */
  T scalarProduct(const string& Direction, unsigned int D2, unsigned int d2,
                  double eps, unsigned int maxNumSweeps, double cutoff,
                  vector<double>& ErrorsAchieved1, vector<unsigned int>& NumsSweepsDone1,
                  vector<double>& ErrorsAchieved2, vector<unsigned int>& NumsSweepsDone2) const;

/// Takes scalar product of this PEPS with itself using MPDOs.
/** This function computes the scalar product of this PEPS with itself using MPDOs.
    If the contraction Direction is set to "horizontal", two neighbouring columns are approximated by one
    new column; whereas if Direction is "vertical", two neigbouring rows are approximated by one row. In
    the case of open BC the MPO-MPS approximation is done with the friend function
       friend void multiplyMPOMPS<>(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                    unsigned int maxNumSweeps, double& errorAchieved,
                                    unsigned int& numSweepsDone, MPS<T>& MPS1,
                                    double cutoff, unsigned int mode, unsigned int d2)   .
    The vector D2s specifies the maximal virtual bond dimensions D2 for the approximating MPSs,
    enumerated from left to right if Direction is "horizontal" or up to down if Direction is "vertical".
    d2s captures the maximal ranks of the hermitian tensors of the approximating MPSs.
    Epss defines the convergence precisions and MaxNumsSweeps the maximal numbers of sweeps allowed.
    The vectors ErrorsAchieved and NumsSweepsDone capture the respective resulting values.
    The two halves are contracted independently, such that this->N[1]-3 column-column approximations are
    done if Direction is "horizontal", or this->N[0]-3 row-row approximations are done if Direction is
    "vertical".
    \param Direction input: const string&, the contraction direction, must be "horizontal" or "vertical"
    \param D2s input: const vector<unsigned int>&, the maximal virtual bond dimensions,
                      must fulfill D2s.size()==this->N[1]-3 if Direction=="horizontal" or
                      D2s.size()==this->N[0]-3 if Direction=="vertical"
    \param d2s input: const vector<unsigned int>&, the maximal tensor ranks,
                      must fulfill d2s.size()==this->N[1]-3 if Direction=="horizontal" or
                      d2s.size()==this->N[0]-3 if Direction=="vertical"
    \param Epss input: const vector<double>&, the convergence precisions,
                       must fulfill Epss.size()==this->N[1]-3 if Direction=="horizontal" or
                       Epss.size()==this->N[0]-3 if Direction=="vertical"
    \param MaxNumsSweeps input: const vector<unsigned int>&, the maximal numbers of sweeps allowed,
                                must fulfill MaxNumsSweeps.size()==this->N[1]-3 if Direction=="horizontal"
                                or MaxNumsSweeps.size()==this->N[0]-3 if Direction=="vertical"
    \param ErrorsAchieved output: vector<double>&, the achieved approximation errors,
                                  must fulfill ErrorsAchieved.size()==this->N[1]-3 if
                                  Direction=="horizontal" or
                                  ErrorsAchieved.size()==this->N[0]-3 if Direction=="vertical"
    \param NumsSweepsDone output: vector<unsigned int>&, the numbers of sweeps done,
                                  must fulfill NumsSweepsDone.size()==this->N[1]-3 if
                                  Direction=="horizontal" or
                                  NumsSweepsDone.size()==this->N[0]-3 if Direction=="vertical"
    \return T, the resulting value of the scalar product
    \sa template<class T> void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                              unsigned int maxNumSweeps, double& errorAchieved,
                                              unsigned int& numSweepsDone, MPS<T>& MPS1,
                                              double cutoff, unsigned int mode, unsigned int d2)
    \sa template<class T> T MPS<T>::contractReverse(const MPS<T>& MPS0) const
    \sa template<class T> T contractReverseMPOMPS(const MPS<T>& MPS0, const MPO<T>& MPO0,
                                                  const MPS<T>& MPS1) */
  T scalarProduct(const string& Direction, const vector<unsigned int>& D2s, const vector<unsigned int>& d2s,
                  const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps,
                  vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const;

/// Normalizes this PEPS.
/** This function normalizes this PEPS.
    This PEPS is normalized such that the PEPS sandwich <thisPEPS|thisPEPS> is 1, and each tensor is normalized
    such that all tensors have the same largest element in absolute value.
    If the contraction Direction is set to "horizontal" then two neighbouring columns are approximated by one
    new column, whereas if Direction is "vertical" then two neigbouring rows are approximated by one row.
    D2 specifies the maximal virtual bond dimension for the approximating boundary-MPSs, eps defines the
    convergence precision and maxNumSweeps the maximal number of sweeps allowed.
    This function first uses
       T PEPS<T>::normalize()
    and then
       T PEPS<T>::scalarProduct(const string& Direction,
                                unsigned int D2, double eps, unsigned int maxNumSweeps,
                                vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone,
                                vector<T>& Norms) const   .
    \param Direction input: const string&, the contraction direction, must be "horizontal" or "vertical"
    \param D2 input: unsigned int, the maximal virtual bond dimension, must be > 0
    \param eps input: double, the convergence precision
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed, must be > 0
    \return T, the absolute value of the new largest element in absolute value of all tensors
    \sa template<class T> T PEPS<T>::normalize()
    \sa template<class T> T PEPS<T>::scalarProduct(const string& Direction,
                                                   unsigned int D2, double eps, unsigned int maxNumSweeps,
                                                   vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone,
                                                   vector<T>& Norms) const */
  T normalize(const string& Direction, unsigned int D2, double eps, unsigned int maxNumSweeps);

/// Takes scalar product of this PEPS with another PEPS.
/** This function computes the scalar product of this PEPS as bra with PEPS0 as ket.
    If the contraction Direction is set to "horizontal", two neighbouring columns are approximated by one
    new column; whereas if Direction is "vertical", two neigbouring rows are approximated by one row. In
    the case of open BC the MPO-MPS approximation is done with the friend function
       friend void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                  unsigned int maxNumSweeps, double& errorAchieved,
                                  unsigned int& numSweepsDone, MPS<T>& MPS1)   .
    The vector D2s specifies the maximal virtual bond dimensions D2 for the approximating MPSs,
    enumerated from left to right if Direction is "horizontal" or up to down if Direction is "vertical".
    Epss defines the convergence precisions and MaxNumsSweeps the maximal numbers of sweeps allowed.
    The vectors ErrorsAchieved and NumsSweepsDone capture the respective resulting values.
    The two halves are contracted independently, such that this->N[1]-3 column-column approximations are
    done if Direction is "horizontal", or this->N[0]-3 row-row approximations are done if Direction is
    "vertical".
    \param PEPS0 input: const PEPS<T>&, the ket PEPS
    \param Direction input: const string&, the contraction direction, must be "horizontal" or "vertical"
    \param D2s input: const vector<unsigned int>&, the maximal virtual bond dimensions,
                      must fulfill D2s.size()==this->N[1]-3 if Direction=="horizontal" or
                      D2s.size()==this->N[0]-3 if Direction=="vertical"
    \param Epss input: const vector<double>&, the convergence precisions,
                       must fulfill Epss.size()==this->N[1]-3 if Direction=="horizontal" or
                       Epss.size()==this->N[0]-3 if Direction=="vertical"
    \param MaxNumsSweeps input: const vector<unsigned int>&, the maximal numbers of sweeps allowed,
                                must fulfill MaxNumsSweeps.size()==this->N[1]-3 if Direction=="horizontal"
                                or MaxNumsSweeps.size()==this->N[0]-3 if Direction=="vertical"
    \param ErrorsAchieved output: vector<double>&, the achieved approximation errors,
                                  must fulfill ErrorsAchieved.size()==this->N[1]-3 if
                                  Direction=="horizontal" or
                                  ErrorsAchieved.size()==this->N[0]-3 if Direction=="vertical"
    \param NumsSweepsDone output: vector<unsigned int>&, the numbers of sweeps done,
                                  must fulfill NumsSweepsDone.size()==this->N[1]-3 if
                                  Direction=="horizontal" or
                                  NumsSweepsDone.size()==this->N[0]-3 if Direction=="vertical"
    \return T, the resulting value of the scalar product
    \sa template<class T> friend void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                                     unsigned int maxNumSweeps, double& errorAchieved,
                                                     unsigned int& numSweepsDone, MPS<T>& MPS1)
    \sa template<class T> T MPS<T>::contractReverse(const MPS<T>& MPS0) const
    \sa template<class T> T contractReverseMPOMPS(const MPS<T>& MPS0, const MPO<T>& MPO0,
                                                  const MPS<T>& MPS1) */
  T scalarProduct(const PEPS<T>& PEPS0, const string& Direction, const vector<unsigned int>& D2s,
                  const vector<double>& Epss,
                  const vector<unsigned int>& MaxNumsSweeps, vector<double>& ErrorsAchieved,
                  vector<unsigned int>& NumsSweepsDone) const;

/// Takes scalar product of this PEPS with another PEPS.
/** This function computes the scalar product of this PEPS as bra with PEPS0 as ket, after normalizing the
    columns or rows of this PEPS sandwich.
    If the contraction Direction is set to "horizontal", two neighbouring columns are approximated by one
    new column; whereas if Direction is "vertical", two neigbouring rows are approximated by one row. In
    the case of open BC the MPO-MPS approximation is done with the friend function
       friend void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                  unsigned int maxNumSweeps, double& errorAchieved,
                                  unsigned int& numSweepsDone, MPS<T>& MPS1)   .
    The vector D2s specifies the maximal virtual bond dimensions D2 for the approximating MPSs,
    enumerated from left to right if Direction is "horizontal" or up to down if Direction is "vertical".
    Epss defines the convergence precisions and MaxNumsSweeps the maximal numbers of sweeps allowed.
    The vectors ErrorsAchieved and NumsSweepsDone capture the respective resulting values.
    Before the contraction of this PEPS sandwich, its columns or rows are normalized with
       T MPS<T>::simplifiedNormalize()   and   T MPO<T>::simplifiedNormalize()
    and their norms are returned in the vector Norms.
    The two halves are contracted independently, such that this->N[1]-3 column-column approximations are
    done if Direction is "horizontal", or this->N[0]-3 row-row approximations are done if Direction is
    "vertical".
    \param PEPS0 input: const PEPS<T>&, the ket PEPS
    \param Direction input: const string&, the contraction direction, must be "horizontal" or "vertical"
    \param D2s input: const vector<unsigned int>&, the maximal virtual bond dimensions,
                      must fulfill D2s.size()==this->N[1]-3 if Direction=="horizontal" or
                      D2s.size()==this->N[0]-3 if Direction=="vertical"
    \param Epss input: const vector<double>&, the convergence precisions,
                       must fulfill Epss.size()==this->N[1]-3 if Direction=="horizontal" or
                       Epss.size()==this->N[0]-3 if Direction=="vertical"
    \param MaxNumsSweeps input: const vector<unsigned int>&, the maximal numbers of sweeps allowed,
                                must fulfill MaxNumsSweeps.size()==this->N[1]-3 if Direction=="horizontal"
                                or MaxNumsSweeps.size()==this->N[0]-3 if Direction=="vertical"
    \param ErrorsAchieved output: vector<double>&, the achieved approximation errors,
                                  must fulfill ErrorsAchieved.size()==this->N[1]-3 if
                                  Direction=="horizontal" or
                                  ErrorsAchieved.size()==this->N[0]-3 if Direction=="vertical"
    \param NumsSweepsDone output: vector<unsigned int>&, the numbers of sweeps done,
                                  must fulfill NumsSweepsDone.size()==this->N[1]-3 if
                                  Direction=="horizontal" or
                                  NumsSweepsDone.size()==this->N[0]-3 if Direction=="vertical"
    \param Norms output: vector<T>&, the norms of the columns or rows of this PEPS sandwich,
                         must fulfill Norms.size()==this->N[1] if Direction=="horizontal" or
                         Norms.size()==this->N[0] if Direction=="vertical"
    \return T, the resulting value of the scalar product excluding the norms
    \sa template<class T> T MPS<T>::simplifiedNormalize()
    \sa template<class T> T MPO<T>::simplifiedNormalize()
    \sa template<class T> friend void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                                     unsigned int maxNumSweeps, double& errorAchieved,
                                                     unsigned int& numSweepsDone, MPS<T>& MPS1)
    \sa template<class T> T MPS<T>::contractReverse(const MPS<T>& MPS0) const
    \sa template<class T> T contractReverseMPOMPS(const MPS<T>& MPS0, const MPO<T>& MPO0,
                                                  const MPS<T>& MPS1) */
  T scalarProduct(const PEPS<T>& PEPS0, const string& Direction, const vector<unsigned int>& D2s,
                  const vector<double>& Epss,
                  const vector<unsigned int>& MaxNumsSweeps, vector<double>& ErrorsAchieved,
                  vector<unsigned int>& NumsSweepsDone, vector<T>& Norms) const;

/// Takes scalar product of this PEPS with another PEPS.
/** This function computes the scalar product of this PEPS as bra with PEPS0 as ket, i.e. the PEPS sandwich
       <thisPEPS|PEPS0>   .
    If the contraction Direction is set to "horizontal" then two neighbouring columns are approximated by one
    new column, whereas if Direction is "vertical" then two neigbouring rows are approximated by one row.
    D2 specifies the maximal virtual bond dimension for the approximating boundary-MPSs, eps defines the
    convergence precision and maxNumSweeps the maximal number of sweeps allowed.
    The vectors EpssAchieved and NumsSweepsDone capture the respective resulting values.
    All boundary-MPSs are normalized and their norms are returned in the vector Norms.
    \param PEPS0 input: const PEPS<T>&, the ket PEPS
    \param Direction input: const string&, the contraction direction, must be "horizontal" or "vertical"
    \param D2 input: unsigned int, the maximal virtual bond dimension, must be > 0
    \param eps input: double, the convergence precision
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed, must be > 0
    \param EpssAchieved output: vector<double>&, the achieved convergence precisions,
                                must fulfill EpssAchieved.size()==this->N[1]-1 if
                                Direction=="horizontal" or
                                EpssAchieved.size()==this->N[0]-1 if Direction=="vertical"
    \param NumsSweepsDone output: vector<unsigned int>&, the numbers of sweeps done,
                                  must fulfill NumsSweepsDone.size()==this->N[1]-1 if
                                  Direction=="horizontal" or
                                  NumsSweepsDone.size()==this->N[0]-1 if Direction=="vertical"
    \param Norms output: vector<T>&, the norms of the boundary-MPSs of the PEPS sandwich,
                         must fulfill Norms.size()==this->N[1]-1 if Direction=="horizontal" or
                         Norms.size()==this->N[0]-1 if Direction=="vertical"
    \return T, the resulting value of the scalar product excluding the norms,
            has to be multiplied by the product of all values in Norms to get the true scalar product
    \sa template<class T> T MPS<T>::normalize()
    \sa template<class T> void PEPS<T>::void getBoundaryMPS(const PEPS<T>& PEPS0, const string& Boundary,
                                                            double eps, unsigned int maxNumSweeps,
                                                            double& epsAchieved, unsigned int& numSweepsDone,
                                                            MPS<T>& MPS0) const
    \sa template<class T> void PEPS<T>::multiplyBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction,
                                                           unsigned int position, const MPS<T>& MPS0,
                                                           double eps, unsigned int maxNumSweeps,
                                                           double& epsAchieved, unsigned int& numSweepsDone,
                                                           MPS<T>& MPS1) const
    \sa template<class T> T PEPS<T>::expectationValue(const vector< vector<unsigned int> >& PositionsRow,
                                                      const vector< vector<unsigned int> >& PositionsCol,
                                                      const vector< vector< Matrix<T> > >& Interactions,
                                                      const string& Direction,
                                                      unsigned int D2, double eps, unsigned int maxNumSweeps) const */
  T scalarProduct(const PEPS<T>& PEPS0, const string& Direction,
                  unsigned int D2, double eps, unsigned int maxNumSweeps,
                  vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone,
                  vector<T>& Norms) const;

/// Computes expectation value of PEPO.
/** This function computes the expectation value of PEPO0 with this PEPS.
    If the contraction Direction is set to "horizontal", two neighbouring columns are approximated by one
    new column; whereas if Direction is "vertical", two neigbouring rows are approximated by one row. In
    the case of open BC the MPO-MPS approximation is done with the friend function
       friend void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                  unsigned int maxNumSweeps, double& errorAchieved,
                                  unsigned int& numSweepsDone, MPS<T>& MPS1)   .
    The vector D2s specifies the maximal virtual bond dimensions D2 for the approximating MPSs,
    enumerated from left to right if Direction is "horizontal" or up to down if Direction is "vertical".
    Epss defines the convergence precisions and MaxNumsSweeps the maximal numbers of sweeps allowed.
    The vectors ErrorsAchieved and NumsSweepsDone capture the respective resulting values.
    The two halves are contracted independently, such that this->N[1]-3 column-column approximations are
    done if Direction is "horizontal", or this->N[0]-3 row-row approximations are done if Direction is
    "vertical".
    \param PEPO0 input: const PEPO<T>&, the PEPO
    \param Direction input: const string&, the contraction direction, must be "horizontal" or "vertical"
    \param D2s input: const vector<unsigned int>&, the maximal virtual bond dimensions,
                      must fulfill D2s.size()==this->N[1]-3 if Direction=="horizontal" or
                      D2s.size()==this->N[0]-3 if Direction=="vertical"
    \param Epss input: const vector<double>&, the convergence precisions,
                       must fulfill Epss.size()==this->N[1]-3 if Direction=="horizontal" or
                       Epss.size()==this->N[0]-3 if Direction=="vertical"
    \param MaxNumsSweeps input: const vector<unsigned int>&, the maximal numbers of sweeps allowed,
                                must fulfill MaxNumsSweeps.size()==this->N[1]-3 if Direction=="horizontal"
                                or MaxNumsSweeps.size()==this->N[0]-3 if Direction=="vertical"
    \param ErrorsAchieved output: vector<double>&, the achieved approximation errors,
                                  must fulfill ErrorsAchieved.size()==this->N[1]-3 if
                                  Direction=="horizontal" or
                                  ErrorsAchieved.size()==this->N[0]-3 if Direction=="vertical"
    \param NumsSweepsDone output: vector<unsigned int>&, the numbers of sweeps done,
                                  must fulfill NumsSweepsDone.size()==this->N[1]-3 if
                                  Direction=="horizontal" or
                                  NumsSweepsDone.size()==this->N[0]-3 if Direction=="vertical"
    \return T, the resulting expectation value
    \sa template<class T> friend void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                                     unsigned int maxNumSweeps, double& errorAchieved,
                                                     unsigned int& numSweepsDone, MPS<T>& MPS1)
    \sa template<class T> T MPS<T>::contractReverse(const MPS<T>& MPS0) const
    \sa template<class T> T contractReverseMPOMPS(const MPS<T>& MPS0, const MPO<T>& MPO0,
                                                  const MPS<T>& MPS1) */
  T expectationValue(const PEPO<T>& PEPO0, const string& Direction,
                     const vector<unsigned int>& D2s, const vector<double>& Epss,
                     const vector<unsigned int>& MaxNumsSweeps, vector<double>& ErrorsAchieved,
                     vector<unsigned int>& NumsSweepsDone) const;

/// Computes expectation value of Interactions.
/** This function computes the expectation value of Interactions with this PEPS divided by the norm of this PEPS.
    The Interactions are implemented as a vector of row positions, a vector of column positions, and a vector of
    matrices, where each vector entry corresponds to one Interaction, i.e. one term of the sum making up the
    operator.
    This function uses
       T PEPS<T>::scalarProduct(const string& Direction, const vector<unsigned int>& D2s,
                                const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps,
                                vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone,
                                vector<T>& Norms) const
    for computing the norm of this PEPS, and
       T PEPS<T>::scalarProduct(const PEPS<T>& PEPS0, const string& Direction, const vector<unsigned int>& D2s,
                                const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps,
                                vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone,
                                vector<T>& Norms) const
    for computing the expectation values of the individual Interactions.
    If the contraction Direction is set to "horizontal", two neighbouring columns are approximated by one
    new column; whereas if Direction is "vertical", two neigbouring rows are approximated by one row. In
    the case of open BC the MPO-MPS approximation is done with the friend function
       friend void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                  unsigned int maxNumSweeps, double& errorAchieved,
                                  unsigned int& numSweepsDone, MPS<T>& MPS1)   .
    The vector D2s specifies the maximal virtual bond dimensions D2 for the approximating MPSs,
    enumerated from left to right if Direction is "horizontal" or up to down if Direction is "vertical".
    Epss defines the convergence precisions and MaxNumsSweeps the maximal numbers of sweeps allowed.
    The vectors ErrorsAchieved and NumsSweepsDone capture the respective resulting values for the Interactions
    contractions.
    Before the contraction of each PEPS sandwich, its columns or rows are normalized with
       T MPS<T>::simplifiedNormalize()   and   T MPO<T>::simplifiedNormalize()   .
    In each sandwich contraction, the two halves are contracted independently, such that this->N[1]-3 column-column
    approximations are done if Direction is "horizontal", or this->N[0]-3 row-row approximations are done if
    Direction is "vertical".
    \param PositionsRow input: const vector< vector<unsigned int> >&, the row positions of the Interactions
    \param PositionsCol input: const vector< vector<unsigned int> >&, the column positions of the Interactions
    \param Interactions input: const vector< vector< Matrix<T> > >&, the Interactions
    \param Direction input: const string&, the contraction direction, must be "horizontal" or "vertical"
    \param D2s input: const vector<unsigned int>&, the maximal virtual bond dimensions,
                      must fulfill D2s.size()==this->N[1]-3 if Direction=="horizontal" or
                      D2s.size()==this->N[0]-3 if Direction=="vertical"
    \param Epss input: const vector<double>&, the convergence precisions,
                       must fulfill Epss.size()==this->N[1]-3 if Direction=="horizontal" or
                       Epss.size()==this->N[0]-3 if Direction=="vertical"
    \param MaxNumsSweeps input: const vector<unsigned int>&, the maximal numbers of sweeps allowed,
                                must fulfill MaxNumsSweeps.size()==this->N[1]-3 if Direction=="horizontal"
                                or MaxNumsSweeps.size()==this->N[0]-3 if Direction=="vertical"
    \param ErrorsAchieved output: vector<double>&, the achieved approximation errors in the Interactions
                                  contractions,
                                  must fulfill ErrorsAchieved.size()==(this->N[1]-3)*Interactions.size() if
                                  Direction=="horizontal" or
                                  ErrorsAchieved.size()==(this->N[0]-3)*Interactions.size() if Direction=="vertical"
    \param NumsSweepsDone output: vector<unsigned int>&, the numbers of sweeps done in the Interactions contractions,
                                  must fulfill NumsSweepsDone.size()==(this->N[1]-3)*Interactions.size() if
                                  Direction=="horizontal" or
                                  NumsSweepsDone.size()==(this->N[0]-3)*Interactions.size() if Direction=="vertical"
    \return T, the resulting expectation value divided by the norm
    \sa template<class T> T PEPS<T>::scalarProduct(const string& Direction, const vector<unsigned int>& D2s,
                                                   const vector<double>& Epss,
                                                   const vector<unsigned int>& MaxNumsSweeps,
                                                   vector<double>& ErrorsAchieved,
                                                   vector<unsigned int>& NumsSweepsDone, vector<T>& Norms) const
    \sa template<class T> T PEPS<T>::scalarProduct(const PEPS<T>& PEPS0, const string& Direction,
                                                   const vector<unsigned int>& D2s,
                                                   const vector<double>& Epss,
                                                   const vector<unsigned int>& MaxNumsSweeps,
                                                   vector<double>& ErrorsAchieved,
                                                   vector<unsigned int>& NumsSweepsDone, vector<T>& Norms) const
    \sa template<class T> T MPS<T>::simplifiedNormalize()
    \sa template<class T> T MPO<T>::simplifiedNormalize()
    \sa template<class T> friend void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                                     unsigned int maxNumSweeps, double& errorAchieved,
                                                     unsigned int& numSweepsDone, MPS<T>& MPS1)
    \sa template<class T> T MPS<T>::contractReverse(const MPS<T>& MPS0) const
    \sa template<class T> T contractReverseMPOMPS(const MPS<T>& MPS0, const MPO<T>& MPO0,
                                                  const MPS<T>& MPS1) */
  T expectationValue(const vector< vector<unsigned int> >& PositionsRow,
                     const vector< vector<unsigned int> >& PositionsCol,
                     const vector< vector< Matrix<T> > >& Interactions,
                     const string& Direction, const vector<unsigned int>& D2s,
                     const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps,
                     vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const;

/// Computes expectation value of Interactions.
/** This function computes the expectation value of Interactions with this PEPS divided by the norm of this PEPS.
    The Interactions are implemented as a vector of row-positions, a vector of column-positions, and a vector of
    matrices, where each vector entry corresponds to one Interaction, i.e. one term of the sum making up the
    operator.
    If the contraction Direction is set to "horizontal" then two neighbouring columns are approximated by one
    new column, whereas if Direction is "vertical" then two neigbouring rows are approximated by one row.
    D2 specifies the maximal virtual bond dimension for the approximating boundary-MPSs, eps defines the
    convergence precision and maxNumSweeps the maximal number of sweeps allowed.
    This PEPS can be a Concatenated PEPS of concatenation level M.
    In each sandwich contraction, the two halves are contracted independently, such that this->N[1]-1 column-column
    approximations are done if Direction is "horizontal", or this->N[0]-1 row-row approximations are done if
    Direction is "vertical".
    This function uses
       T PEPS<T>::scalarProduct(const string& Direction,
                                unsigned int D2, double eps, unsigned int maxNumSweeps,
                                vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone,
                                vector<T>& Norms) const
    for computing the norm of this PEPS, and
       T PEPS<T>::scalarProduct(const PEPS<T>& PEPS0, const string& Direction,
                                unsigned int D2, double eps, unsigned int maxNumSweeps,
                                vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone,
                                vector<T>& Norms) const
    for computing the expectation values of the individual Interactions.
    \param PositionsRow input: const vector< vector<unsigned int> >&, the row-positions of the Interactions
    \param PositionsCol input: const vector< vector<unsigned int> >&, the column-positions of the Interactions
    \param Interactions input: const vector< vector< Matrix<T> > >&, the Interactions
    \param Direction input: const string&, the contraction direction, must be "horizontal" or "vertical"
    \param D2 input: unsigned int, the maximal virtual bond dimension, must be > 0
    \param eps input: double, the convergence precision
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed, must be > 0
    \param M optional input: unsigned int, if given the concatenation level, must be > 0
    \return T, the resulting expectation value divided by the norm
    \sa template<class T> T PEPS<T>::scalarProduct(const string& Direction,
                                                   unsigned int D2, double eps, unsigned int maxNumSweeps,
                                                   vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone,
                                                   vector<T>& Norms) const
    \sa template<class T> T PEPS<T>::scalarProduct(const PEPS<T>& PEPS0, const string& Direction,
                                                   unsigned int D2, double eps, unsigned int maxNumSweeps,
                                                   vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone,
                                                   vector<T>& Norms) const */
  T expectationValue(const vector< vector<unsigned int> >& PositionsRow,
                     const vector< vector<unsigned int> >& PositionsCol,
                     const vector< vector< Matrix<T> > >& Interactions,
                     const string& Direction,
                     unsigned int D2, double eps, unsigned int maxNumSweeps,
                     unsigned int M = 1) const;

/// Computes distance between product PEPO*PEPS and another PEPS.
/** This friend function computes the distance
       ||PEPO0*|PEPS0>-|PEPS1>||
    between PEPO0*PEPS0 and PEPS1.
    If the contraction Direction is set to "horizontal", two neighbouring columns are approximated by one
    new column; whereas if Direction is "vertical", two neigbouring rows are approximated by one row. In
    the case of open BC the MPO-MPS approximation is done with the friend function
       friend void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps,
                                  unsigned int maxNumSweeps, double& errorAchieved,
                                  unsigned int& numSweepsDone, MPS<T>& MPS1)   .
    The vector D2s specifies the maximal virtual bond dimensions D2 for the approximating MPSs,
    enumerated from left to right if Direction is "horizontal" or up to down if Direction is "vertical".
    Epss defines the convergence precisions and MaxNumsSweeps the maximal numbers of sweeps allowed.
    The vectors ErrorsAchieved and NumsSweepsDone capture the respective resulting values.
    All vectors contain values for the contraction of the three PEPS sandwiches
    <PEPS0|PEPO0^{+}PEPO0|PEPS0>, <PEPS1|PEPO0|PEPS0>, and <PEPS1|PEPS1>, in this order.
    For each PEPS sandwich, the two halves are contracted independently, such that this->N[1]-3 column-column
    approximations are done if Direction is "horizontal", or this->N[0]-3 row-row approximations are done if
    Direction is "vertical".
    \param PEPO0 input: const PEPO<T>&, the PEPO
    \param PEPS0 input: const PEPS<T>&, the PEPS on which the PEPO acts
    \param PEPS1 input: const PEPS<T>&, the approximating PEPS
    \param Direction input: const string&, the contraction direction, must be "horizontal" or "vertical"
    \param D2s input: const vector<unsigned int>&, the maximal virtual bond dimensions,
                      must fulfill D2s.size()==3*(this->N[1]-3) if Direction=="horizontal" or
                      D2s.size()==3*(this->N[0]-3) if Direction=="vertical"
    \param Epss input: const vector<double>&, the convergence precisions,
                       must fulfill Epss.size()==3*(this->N[1]-3) if Direction=="horizontal" or
                       Epss.size()==3*(this->N[0]-3) if Direction=="vertical"
    \param MaxNumsSweeps input: const vector<unsigned int>&, the maximal numbers of sweeps allowed,
                                must fulfill MaxNumsSweeps.size()==3*(this->N[1]-3) if Direction=="horizontal"
                                or MaxNumsSweeps.size()==3*(this->N[0]-3) if Direction=="vertical"
    \param ErrorsAchieved output: vector<double>&, the achieved approximation errors,
                                  must fulfill ErrorsAchieved.size()==3*(this->N[1]-3) if
                                  Direction=="horizontal" or
                                  ErrorsAchieved.size()==3*(this->N[0]-3) if Direction=="vertical"
    \param NumsSweepsDone output: vector<unsigned int>&, the numbers of sweeps done,
                                  must fulfill NumsSweepsDone.size()==3*(this->N[1]-3) if
                                  Direction=="horizontal" or
                                  NumsSweepsDone.size()==3*(this->N[0]-3) if Direction=="vertical"
    \return double, the resulting distance */
  friend double distancePEPOPEPS<>(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const PEPS<T>& PEPS1,
                                   const string& Direction, const vector<unsigned int>& D2s, const vector<double>& Epss,
                                   const vector<unsigned int>& MaxNumsSweeps, vector<double>& ErrorsAchieved,
                                   vector<unsigned int>& NumsSweepsDone);

/// Multiplies PEPO with PEPS exactly.
/** This friend function multiplies PEPO0 with PEPS0 and stores the result in PEPS1.
    \param PEPO0 input: const PEPO<T>&, the PEPO
    \param PEPS0 input: const PEPS<T>&, the PEPS on which the PEPO acts
    \param PEPS1 output: PEPS<T>&, the resulting product of PEPO0 and PEPS0,
                         must have the correct form */
  friend void multiplyPEPOPEPSExact<>(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, PEPS<T>& PEPS1);

/// Multiplies PEPO with PEPS and approximates result by PEPS.
/** This friend function approximates the product of PEPO0 with PEPS0 by PEPS1 with a given bond dimension.
    The cost function
       |||PEPS1>-PEPO0*|PEPS0>||^{2}
    is minimized by sweeping over the tensors and updating them column by column.
    The simplified approximation error
       error := <PEPS1|PEPS1>-2Re(<PEPS1|PEPO0|PEPS0>)
    is evaluated after each tensor-update, and the function stops if the relative change in the error is
    below eps, i.e. in tensor-update n: eps > abs(error(n)-error(n-1))/abs(error(n)).
    The function always stops if the number of sweeps exceeds maxNumSweeps.
    Environment boundary-MPSs are computed with D2Env, epsEnv, and maxNumSweepsEnv.
    The local tensor-update can be reduced from the full tensor to the smallest subtensor including the physical bond.
    Either the pseudoinverse or the closest positive norm-matrix is used in the local update algorithm.
    On input, if (PEPS1.D >= PEPO0.D*PEPS0.D) then the exact product with bond dimension PEPS1.D is returned by filling in zeroes.
    On output, epsAchieved is the achieved relative approximation error after the last tensor-update,
    numSweepsDone is the number of complete sweeps done excluding the last partial sweep, and PEPS1 is the
    approximating PEPS.
    \param PEPO0 input: const PEPO<T>&, the PEPO
    \param PEPS0 input: const PEPS<T>&, the PEPS on which the PEPO acts
    \param eps input: double, the convergence precision, i.e. convergence in tensor-update n if
                      eps > abs(error(n)-error(n-1))/abs(error(n))
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param D2Env input: unsigned int, the maximal virtual bond dimension in the environment approximation
    \param epsEnv input: double, the convergence precision in the environment approximation
    \param maxNumSweepsEnv input: unsigned int, the maximal number of sweeps allowed in the environment approximation
    \param UpdateTensor input: const string&, the tensor in the local tensor update, must be "full" or "reduced"
    \param UpdateMode input: const string&, the local update algorithm, must be "Pseudoinverse" for pseudoinverse or
                             "PositiveNormMatrix" for the closest positive norm-matrix
    \param cutoff input: double, the cutoff, singular values or eigenvalues x_{i} with (x_{i}/x_{max} <= cutoff) are put to zero
    \param epsAchieved output: double&, the achieved convergence precision
    \param numSweepsDone output: unsigned int&, the number of complete sweeps done excluding the last partial sweep
    \param PEPS1 input/output: PEPS<T>&, on input the initial approximating PEPS, on output the resulting PEPS,
                               must have the correct form */
  friend void multiplyPEPOPEPS<>(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0,
                                 double eps, unsigned int maxNumSweeps,
                                 unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv,
                                 const string& UpdateTensor, const string& UpdateMode, double cutoff,
                                 double& epsAchieved, unsigned int& numSweepsDone,
                                 PEPS<T>& PEPS1);

  void getMPS(const string& Direction, const vector<unsigned int>& dIndices, MPS<T>& MPS0) const;

  void getMPO(const string& Direction, unsigned int position, const vector<unsigned int>& dIndices, MPO<T>& MPO0)
              const;

  void sampleExpectationValue(unsigned int positionRow, unsigned int positionCol, const Matrix<T>& Observable,
                              unsigned int D1, const string& Filling, double eps, unsigned int maxNumSweeps,
                              unsigned int MEqu, unsigned int M, unsigned int MCSteps,
                              const string& FileName) const;

  void sampleExpectationValue(const Matrix<T>& Observable,
                              unsigned int D1, const string& Filling, double eps, unsigned int maxNumSweeps,
                              unsigned int MEqu, unsigned int M, unsigned int MCSteps,
                              const string& FileName) const;

/// Returns Shape to position for open boundary conditions.
/** Given a row position positionRow and a column position positionCol in this PEPS with open boundary
    conditions, the correct Shape of the corresponding Tensor is returned.
    This PEPS is not changed.
    \param positionRow input: unsigned int, the row position
    \param positionCol input: unsigned int, the column position
    \param Shape output: vector<unsigned int>&, the resulting Shape, on input it must fulfill
                          Shape.size()==5, on output it will contain {Dleft, Dup, Dright, Ddown, d} */
  void getOpenBCShape(unsigned int positionRow, unsigned int positionCol,
                      vector<unsigned int>& Shape);

 protected:

/// Boundary conditions BC.
/** The boundary conditions are either "open" or "periodic". */
  string BC;

/// Number of tensors N.
/** N[0]=Nrows is the number of rows and N[1]=Ncols is the number of columns of this PEPS. */
  vector<unsigned int> N;

/// Physical dimensions d.
/** The Matrix d has dimensions d.getDim0()==N[0] and d.getDim1()==N[1]. */
  Matrix<unsigned int> d;

/// Maximal virtual bond dimension D.
  unsigned int D;

/// Tensors.
/** The Tensors of this PEPS are stored in Fortran's column-major order, such that the Tensor at matrix
    position (row, col) is located in Tensors at linear position (row + Nrows*col). Each PEPS tensor has the
    Shape {Dleft, Dup, Dright, Ddown, d}. */
  Tensor<T>* Tensors;

 private:

/// Returns initial column norm-MPSs.
/** This function computes the initial boundary-MPSs for the columns of the norm-sandwich <thisPEPS|thisPEPS> of this PEPS, before the columns get updated.
    ColumnNormMPSs will have size this->N[1] and comprise all norm-MPSs for column-positions > 0.
    The norm-MPS for column position will be ColumnNormMPSs[position] and its physical indices will be pointing left.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices. */
  void getInitialColumnNormMPSs(unsigned int D2, double eps, unsigned int maxNumSweeps,
                                vector< MPS<T> >& ColumnNormMPSs) const;

/// Returns initial column b-MPSs.
/** This function computes the initial boundary-MPSs for the columns of the sandwich <thisPEPS|PEPO0|PEPS0>, before the columns of this PEPS get updated.
    ColumnbMPSs will have size this->N[1] and comprise all b-MPSs for column-positions > 0.
    The b-MPS for column position will be ColumnbMPSs[position] and its physical indices will be pointing left.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices. */
  void getInitialColumnbMPSs(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0,
                             unsigned int D2, double eps, unsigned int maxNumSweeps,
                             vector< MPS<T> >& ColumnbMPSs) const;

/// Updates column norm-MPSs on boundary.
/** This function updates the boundary-MPSs for the columns of the norm-sandwich <thisPEPS|thisPEPS> of this PEPS, after the update of a boundary column.
    If (Boundary=="left") then the boundary-MPS for the left boundary will be written into ColumnNormMPSs[0] with physical indices pointing right.
    If (Boundary=="right") then the boundary-MPS for the right boundary will be written into ColumnNormMPSs[this->N[1]-1] with physical indices
    pointing left.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices. */
  void updateColumnNormMPSsBoundary(const string& Boundary,
                                    unsigned int D2, double eps, unsigned int maxNumSweeps,
                                    vector< MPS<T> >& ColumnNormMPSs) const;

/// Updates column b-MPSs on boundary.
/** This function updates the boundary-MPSs for the columns of the sandwich <thisPEPS|PEPO0|PEPS0>, after the update of a boundary column of this PEPS.
    If (Boundary=="left") then the boundary-MPS for the left boundary will be written into ColumnbMPSs[0] with physical indices pointing right.
    If (Boundary=="right") then the boundary-MPS for the right boundary will be written into ColumnbMPSs[this->N[1]-1] with physical indices
    pointing left.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices. */
  void updateColumnbMPSsBoundary(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary,
                                 unsigned int D2, double eps, unsigned int maxNumSweeps,
                                 vector< MPS<T> >& ColumnbMPSs) const;

/// Updates column norm-MPSs in bulk.
/** This function updates the boundary-MPSs for the columns of the norm-sandwich <thisPEPS|thisPEPS> of this PEPS, after the update of a bulk column.
    If (Direction=="right") then the boundary-MPS for column position will be written into ColumnNormMPSs[position] with physical
    indices pointing right.
    If (Direction=="left") then the boundary-MPS for column position will be written into ColumnNormMPSs[position] with physical
    indices pointing left.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices. */
  void updateColumnNormMPSsBulk(unsigned int position, const string& Direction,
                                unsigned int D2, double eps, unsigned int maxNumSweeps,
                                vector< MPS<T> >& ColumnNormMPSs) const;

/// Updates column b-MPSs in bulk.
/** This function updates the boundary-MPSs for the columns of the sandwich <thisPEPS|PEPO0|PEPS0>, after the update of a bulk column of this PEPS.
    If (Direction=="right") then the boundary-MPS for column position will be written into ColumnbMPSs[position] with physical
    indices pointing right.
    If (Direction=="left") then the boundary-MPS for column position will be written into ColumnbMPSs[position] with physical
    indices pointing left.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices. */
  void updateColumnbMPSsBulk(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, unsigned int position, const string& Direction,
                             unsigned int D2, double eps, unsigned int maxNumSweeps,
                             vector< MPS<T> >& ColumnbMPSs) const;

/// Returns initial column norm-tensors for boundary.
/** This function computes the initial norm-tensors for a boundary-column of the sandwich <thisPEPS|thisPEPS>, by contracting that boundary with
    its neighbouring column norm-MPS. The resulting tensors are the initial norm-tensors for the update of that boundary.
    They lie above and below the updated tensor, and are enumerated like the rows in this PEPS.
    Each norm-tensor has the shape required by the norm-MPO. */
  void getInitialColumnNormTensorsBoundary(const string& Boundary, const MPS<T>& ColumnNormMPS,
                                           vector< Tensor<T> >& ColumnNormTensors) const;

/// Returns initial column b-tensors for boundary.
/** This function computes the initial b-tensors for a boundary-column of the sandwich <thisPEPS|PEPO0|PEPS0>, by contracting that boundary with
    its neighbouring column b-MPS. The resulting tensors are the initial b-tensors for the update of that boundary.
    They lie above and below the updated tensor, and are enumerated like the rows in this PEPS.
    Each b-tensor still has to be reshaped for the b-MPO. */
  void getInitialColumnbTensorsBoundary(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary,
                                        const MPS<T>& ColumnbMPS,
                                        vector< Tensor<T> >& ColumnbTensors) const;

/// Returns initial column norm-tensors for bulk.
/** This function computes the initial norm-tensors for a bulk-column of the sandwich <thisPEPS|thisPEPS>, by contracting that column with
    its neighbouring column norm-MPSs. The resulting tensors are the initial norm-tensors for the update of that column.
    ColumnNormMPSL denotes the left and ColumnNormMPSR the right norm-MPS to column position.
    The norm-tensors lie above and below the updated tensor, and are enumerated like the rows in this PEPS.
    Each norm-tensor has the shape required by the norm-MPO. */
  void getInitialColumnNormTensorsBulk(unsigned int position, const MPS<T>& ColumnNormMPSL, const MPS<T>& ColumnNormMPSR,
                                       vector< Tensor<T> >& ColumnNormTensors) const;

/// Returns initial column b-tensors for bulk.
/** This function computes the initial b-tensors for a bulk-column of the sandwich <thisPEPS|PEPO0|PEPS0>, by contracting that column with
    its neighbouring column b-MPSs. The resulting tensors are the initial b-tensors for the update of that column.
    ColumnbMPSL denotes the left and ColumnbMPSR the right b-MPS to column position.
    The b-tensors lie above and below the updated tensor, and are enumerated like the rows in this PEPS.
    Each b-tensor still has to be reshaped for the b-MPO. */
  void getInitialColumnbTensorsBulk(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, unsigned int position,
                                    const MPS<T>& ColumnbMPSL, const MPS<T>& ColumnbMPSR,
                                    vector< Tensor<T> >& ColumnbTensors) const;

/// Updates column norm-tensors on boundary.
/** This function updates the norm-tensors of a boundary-column of the sandwich <thisPEPS|thisPEPS>, after the update of a tensor on the boundary.
    If (Boundary=="left") then ColumnNormMPS is the boundary-MPS on column 1 with physical indices pointing left, and tensorPosition denotes the
    row-position of the previously updated tensor.
    If (Boundary=="right") then ColumnNormMPS is the boundary-MPS on column this->N[1]-2 with physical indices pointing right, and tensorPosition denotes the
    row-position of the previously updated tensor.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices.
    Each norm-tensor has the shape required by the norm-MPO. */
  void updateColumnNormTensorsBoundary(const string& Boundary, unsigned int tensorPosition, const MPS<T>& ColumnNormMPS,
                                       vector< Tensor<T> >& ColumnNormTensors) const;

/// Updates column b-tensors on boundary.
/** This function updates the b-tensors of a boundary-column of the sandwich <thisPEPS|PEPO0|PEPS0>, after the update of a tensor on the boundary.
    If (Boundary=="left") then ColumnbMPS is the boundary-MPS on column 1 with physical indices pointing left, and tensorPosition denotes the
    row-position of the previously updated tensor.
    If (Boundary=="right") then ColumnbMPS is the boundary-MPS on column this->N[1]-2 with physical indices pointing right, and tensorPosition denotes the
    row-position of the previously updated tensor.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices.
    Each b-tensor still has to be reshaped for the b-MPO. */
  void updateColumnbTensorsBoundary(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary,
                                    unsigned int tensorPosition, const MPS<T>& ColumnbMPS,
                                    vector< Tensor<T> >& ColumnbTensors) const;

/// Updates column norm-tensors in bulk.
/** This function updates the norm-tensors of a bulk-column of the sandwich <thisPEPS|thisPEPS>, after the update of a tensor in the bulk.
    ColumnNormMPSL denotes the left and ColumnNormMPSR the right boundary-MPS to column position, and tensorPosition denotes the
    row-position of the previously updated tensor.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices.
    Each norm-tensor has the shape required by the norm-MPO. */
  void updateColumnNormTensorsBulk(unsigned int position, unsigned int tensorPosition,
                                   const MPS<T>& ColumnNormMPSL, const MPS<T>& ColumnNormMPSR,
                                   vector< Tensor<T> >& ColumnNormTensors) const;

/// Updates column b-tensors in bulk.
/** This function updates the b-tensors of a bulk-column of the sandwich <thisPEPS|PEPO0|PEPS0>, after the update of a tensor in the bulk.
    ColumnbMPSL denotes the left and ColumnbMPSR the right boundary-MPS to column position, and tensorPosition denotes the
    row-position of the previously updated tensor.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices.
    Each b-tensor still has to be reshaped for the b-MPO. */
  void updateColumnbTensorsBulk(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, unsigned int position,
                                unsigned int tensorPosition, const MPS<T>& ColumnbMPSL, const MPS<T>& ColumnbMPSR,
                                vector< Tensor<T> >& ColumnbTensors) const;

/// Returns column norm-MPO for boundary.
/** This function computes the norm-MPO for the update of a tensor in a boundary-column. The norm-MPO is the periodic MPO consisting of the 4 tensors
    surrounding the tensor in the norm-sandwich <thisPEPS|thisPEPS>.
    If (Boundary=="left") then ColumnNormMPS denotes the boundary-MPS on column 1 with physical indices pointing left, and tensorPosition denotes
    the row-position of the tensor to be updated.
    If (Boundary=="right") then ColumnNormMPS denotes the boundary-MPS on column this->N[1]-2 with physical indices pointing right, and tensorPosition denotes
    the row-position of the tensor to be updated.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices.
    ColumnNormTensors stores the norm-tensors with required shape for ColumnNormMPO.
    ColumnNormMPO is enumerated clockwise around the tensor in the sequence of its virtual indices, starting from the left virtual index 0. */
  void getColumnNormMPOBoundary(const string& Boundary, unsigned int tensorPosition,
                                const MPS<T>& ColumnNormMPS, const vector< Tensor<T> >& ColumnNormTensors,
                                MPO<T>& ColumnNormMPO) const;

/// Returns column b-MPO for boundary.
/** This function computes the b-MPO for the update of a tensor in a boundary-column. The b-MPO is the periodic MPO consisting of the 4 tensors
    surrounding the tensor in the sandwich <thisPEPS|PEPO0|PEPS0>. Its physical index 2 comprises the virtual indices of PEPS0 and PEPO0, and
    its physical index 3 the virtual index of this PEPS.
    If (Boundary=="left") then ColumnbMPS denotes the boundary-MPS on column 1 with physical indices pointing left, and tensorPosition denotes
    the row-position of the tensor to be updated.
    If (Boundary=="right") then ColumnbMPS denotes the boundary-MPS on column this->N[1]-2 with physical indices pointing right, and tensorPosition denotes
    the row-position of the tensor to be updated.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices.
    ColumnbTensors stores the b-tensors, which still have to be reshaped for ColumnbMPO.
    ColumnbMPO is enumerated clockwise around the tensor in the sequence of its virtual indices, starting from the left virtual index 0. */
  void getColumnbMPOBoundary(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary, unsigned int tensorPosition,
                             const MPS<T>& ColumnbMPS, const vector< Tensor<T> >& ColumnbTensors,
                             MPO<T>& ColumnbMPO) const;

/// Returns column norm-MPO for bulk.
/** This function computes the norm-MPO for the update of a tensor in a bulk-column. The norm-MPO is the periodic MPO consisting of the 4 tensors
    surrounding the tensor in the norm-sandwich <thisPEPS|thisPEPS>.
    ColumnNormMPSL denotes the left and ColumnNormMPSR the right boundary-MPS to column position, and tensorPosition denotes the
    row-position of the tensor to be updated.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices.
    ColumnNormTensors stores the norm-tensors with required shape for ColumnNormMPO.
    ColumnNormMPO is enumerated clockwise around the tensor in the sequence of its virtual indices, starting from the left virtual index 0. */
  void getColumnNormMPOBulk(unsigned int position, unsigned int tensorPosition,
                            const MPS<T>& ColumnNormMPSL, const MPS<T>& ColumnNormMPSR, const vector< Tensor<T> >& ColumnNormTensors,
                            MPO<T>& ColumnNormMPO) const;

/// Returns column b-MPO for bulk.
/** This function computes the b-MPO for the update of a tensor in a bulk-column. The b-MPO is the periodic MPO consisting of the 4 tensors
    surrounding the tensor in the sandwich <thisPEPS|PEPO0|PEPS0>. Its physical index 2 comprises the virtual indices of PEPS0 and PEPO0, and
    its physical index 3 the virtual index of this PEPS.
    ColumnbMPSL denotes the left and ColumnbMPSR the right boundary-MPS to column position, and tensorPosition denotes the
    row-position of the tensor to be updated.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices.
    ColumnbTensors stores the b-tensors, which still have to be reshaped for ColumnbMPO.
    ColumnbMPO is enumerated clockwise around the tensor in the sequence of its virtual indices, starting from the left virtual index 0. */
  void getColumnbMPOBulk(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, unsigned int position, unsigned int tensorPosition,
                         const MPS<T>& ColumnbMPSL, const MPS<T>& ColumnbMPSR, const vector< Tensor<T> >& ColumnbTensors,
                         MPO<T>& ColumnbMPO) const;

/// Updates tensor.
/** This function updates a tensor of this PEPS for multiplyPEPOPEPS.
    The cost function
       |||thisPEPS>-PEPO0*|PEPS0>||^{2}
    is minimized for the tensor at row positionRow and column positionCol.
    NormMPO denotes the environment from the sandwich <thisPEPS|thisPEPS> and bMPO from <thisPEPS|PEPO0|PEPS0>.
    The local tensor-update can be reduced from the full tensor to the smallest subtensor including the physical bond.
    Either the pseudoinverse or the closest positive norm-matrix is used in the local update algorithm.
    The simplified approximation error
       error := <thisPEPS|thisPEPS>-2Re(<thisPEPS|PEPO0|PEPS0>)
    is evaluated before and after the tensor-update, and the relative change
       epsAchieved := abs(error_{after}-error_{before})/abs(error_{after})
    is returned.
    On output, this PEPS stores the updated tensor.
    \param PEPO0 input: const PEPO<T>&, the PEPO
    \param PEPS0 input: const PEPS<T>&, the PEPS on which the PEPO acts
    \param positionRow input: unsigned int, the row-position
    \param positionCol input: unsigned int, the column-position
    \param NormMPO input: const MPO<T>&, the periodic norm-MPO,
                          consisting of the 4 tensors surrounding the tensor in <thisPEPS|thisPEPS>,
                          enumerated clockwise around the tensor starting from the left virtual index 0
    \param bMPO input: const MPO<T>&, the periodic b-MPO,
                       consisting of the 4 tensors surrounding the tensor in <thisPEPS|PEPO0|PEPS0>,
                       enumerated clockwise around the tensor starting from the left virtual index 0
    \param UpdateTensor input: const string&, the tensor in the local tensor-update, must be "full" or "reduced"
    \param UpdateMode input: const string&, the local update algorithm, must be "Pseudoinverse" for pseudoinverse or
                             "PositiveNormMatrix" for the closest positive norm-matrix
    \param cutoff input: double, the cutoff, singular values or eigenvalues x_{i} with (x_{i}/x_{max} <= cutoff) are put to zero
    \param epsAchieved output: double&, the achieved convergence precision */
  void updateTensor(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, unsigned int positionRow, unsigned int positionCol,
                    const MPO<T>& NormMPO, const MPO<T>& bMPO,
                    const string& UpdateTensor, const string& UpdateMode, double cutoff,
                    double& epsAchieved);

/// Updates boundary-column.
/** This function updates a boundary-column of this PEPS for multiplyPEPOPEPS.
    The cost function
       |||thisPEPS>-PEPO0*|PEPS0>||^{2}
    is minimized by sweeping once over all tensors in the boundary-column from top to bottom and updating them one after the other.
    The simplified approximation error
       error := <thisPEPS|thisPEPS>-2Re(<thisPEPS|PEPO0|PEPS0>)
    is evaluated after each tensor-update, and the function returns if the relative change in the error is
    below eps, i.e. in tensor-update n: eps > abs(error(n)-error(n-1))/abs(error(n)).
    Environment boundary-MPSs for <thisPEPS|thisPEPS> are stored in ColumnNormMPSs and for <thisPEPS|PEPO0|PEPS0> in ColumnbMPSs.
    The local tensor-update can be reduced from the full tensor to the smallest subtensor including the physical bond.
    Either the pseudoinverse or the closest positive norm-matrix is used in the local update algorithm.
    On output, epsAchieved is the achieved relative approximation error after the last tensor-update, and this PEPS stores the
    updated boundary-column.
    \param PEPO0 input: const PEPO<T>&, the PEPO
    \param PEPS0 input: const PEPS<T>&, the PEPS on which the PEPO acts
    \param Boundary input: const string&, the boundary, must be "left" or "right"
    \param eps input: double, the convergence precision, i.e. convergence in tensor-update n if
                      eps > abs(error(n)-error(n-1))/abs(error(n))
    \param ColumnNormMPSs input: const vector< MPS<T> >&, the column norm-MPSs for <thisPEPS|thisPEPS>
    \param ColumnbMPSs input: const vector< MPS<T> >&, the column b-MPSs for <thisPEPS|PEPO0|PEPS0>
    \param UpdateTensor input: const string&, the tensor in the local tensor-update, must be "full" or "reduced"
    \param UpdateMode input: const string&, the local update algorithm, must be "Pseudoinverse" for pseudoinverse or
                             "PositiveNormMatrix" for the closest positive norm-matrix
    \param cutoff input: double, the cutoff, singular values or eigenvalues x_{i} with (x_{i}/x_{max} <= cutoff) are put to zero
    \param epsAchieved output: double&, the achieved convergence precision */
  void updateBoundaryColumn(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary,
                            double eps,
                            const vector< MPS<T> >& ColumnNormMPSs, const vector< MPS<T> >& ColumnbMPSs,
                            const string& UpdateTensor, const string& UpdateMode, double cutoff,
                            double& epsAchieved);

/// Updates bulk-column.
/** This function updates a bulk-column of this PEPS for multiplyPEPOPEPS.
    The cost function
       |||thisPEPS>-PEPO0*|PEPS0>||^{2}
    is minimized by sweeping once over all tensors in the bulk-column from top to bottom and updating them one after the other.
    The simplified approximation error
       error := <thisPEPS|thisPEPS>-2Re(<thisPEPS|PEPO0|PEPS0>)
    is evaluated after each tensor-update, and the function returns if the relative change in the error is
    below eps, i.e. in tensor-update n: eps > abs(error(n)-error(n-1))/abs(error(n)).
    Environment boundary-MPSs for <thisPEPS|thisPEPS> are stored in ColumnNormMPSs and for <thisPEPS|PEPO0|PEPS0> in ColumnbMPSs.
    As always, the boundary-MPS on the left of column position has physical indices pointing right and the one on the right has
    physical indices pointing left, and the tensors in an MPS are enumerated in the natural way, given by the direction of their
    physical indices.
    The local tensor-update can be reduced from the full tensor to the smallest subtensor including the physical bond.
    Either the pseudoinverse or the closest positive norm-matrix is used in the local update algorithm.
    On output, epsAchieved is the achieved relative approximation error after the last tensor-update, and this PEPS stores the
    updated bulk-column.
    \param PEPO0 input: const PEPO<T>&, the PEPO
    \param PEPS0 input: const PEPS<T>&, the PEPS on which the PEPO acts
    \param position input: unsigned int, the column-position
    \param eps input: double, the convergence precision, i.e. convergence in tensor-update n if
                      eps > abs(error(n)-error(n-1))/abs(error(n))
    \param ColumnNormMPSs input: const vector< MPS<T> >&, the column norm-MPSs for <thisPEPS|thisPEPS>
    \param ColumnbMPSs input: const vector< MPS<T> >&, the column b-MPSs for <thisPEPS|PEPO0|PEPS0>
    \param UpdateTensor input: const string&, the tensor in the local tensor-update, must be "full" or "reduced"
    \param UpdateMode input: const string&, the local update algorithm, must be "Pseudoinverse" for pseudoinverse or
                             "PositiveNormMatrix" for the closest positive norm-matrix
    \param cutoff input: double, the cutoff, singular values or eigenvalues x_{i} with (x_{i}/x_{max} <= cutoff) are put to zero
    \param epsAchieved output: double&, the achieved convergence precision */
  void updateBulkColumn(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, unsigned int position,
                        double eps,
                        const vector< MPS<T> >& ColumnNormMPSs, const vector< MPS<T> >& ColumnbMPSs,
                        const string& UpdateTensor, const string& UpdateMode, double cutoff,
                        double& epsAchieved);

};

template<class T> PEPS<T>::PEPS()
{
 this->N = vector<unsigned int>(2);
 this->N[0] = 0; this->N[1] = 0;
 this->D = 0;
 this->Tensors = 0;
}

template<class T> PEPS<T>::PEPS(const string& BC0, unsigned int Nrows, unsigned int Ncols,
                                const Matrix<unsigned int>& d0, unsigned int D0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (Nrows < 2) || (Ncols < 2) ||
     (d0.getDim0() != Nrows) || (d0.getDim1() != Ncols) || (D0 == 0))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> PEPS<T>::" <<
          "PEPS(const string& BC0, unsigned int Nrows, unsigned int Ncols, " <<
               "const Matrix<unsigned int>& d0, unsigned int D0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (Nrows < 2) || (Ncols < 2) || " <<
           "(d0.getDim0() != Nrows) || (d0.getDim1() != Ncols) || (D0 == 0))." << endl;
  exit(1);
 }
 for (int j = 0; j < d0.getDim1(); j++)
 {
  for (int i = 0; i < d0.getDim0(); i++)
  {
   if (d0(i, j) < 2)
   {
    cerr << "Program terminated because of error in constructor " <<
            "template<class T> PEPS<T>::" <<
            "PEPS(const string& BC0, unsigned int Nrows, unsigned int Ncols, " <<
                 "const Matrix<unsigned int>& d0, unsigned int D0): " <<
            "(d0(" << i << ", " << j << ") < 2)." << endl;
    exit(1);
   }
  }
 }
#endif
 this->BC = BC0;
 this->N = vector<unsigned int>(2);
 this->N[0] = Nrows; this->N[1] = Ncols;
 this->d = d0;
 this->D = D0;
 unsigned int numTensors = this->N[0] * this->N[1];
 this->Tensors = new Tensor<T>[numTensors];
 vector<unsigned int> Shape(5);
 unsigned int position;
 if (this->BC == "open")
 {
  for (int j = 0; j < this->N[1]; j++)
  {
   for (int i = 0; i < this->N[0]; i++)
   {
    this->getOpenBCShape(i, j, Shape);
    position = i + j*this->N[0];
    this->Tensors[position] = Tensor<T>(Shape);
   }
  }
 }
 else if (this->BC == "periodic")
 {
// we put the same virtual shape everywhere, only the physical indices are site dependent:
  Shape[0] = this->D; Shape[1] = this->D; Shape[2] = this->D; Shape[3] = this->D;
  for (int j = 0; j < this->N[1]; j++)
  {
   for (int i = 0; i < this->N[0]; i++)
   {
    Shape[4] = this->d(i, j);
    position = i + j*this->N[0];
    this->Tensors[position] = Tensor<T>(Shape);
   }
  }
 }
}

template<class T> PEPS<T>::PEPS(const string& BC0, unsigned int Nrows, unsigned int Ncols,
                                unsigned int d0, unsigned int D0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (Nrows < 2) || (Ncols < 2) || (d0 < 2) || (D0 == 0))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> PEPS<T>::" <<
          "PEPS(const string& BC0, unsigned int Nrows, unsigned int Ncols, unsigned int d0, " <<
               "unsigned int D0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (Nrows < 2) || (Ncols < 2) || (d0 < 2) || " <<
           "(D0 == 0))." << endl;
  exit(1);
 }
#endif
 this->BC = BC0;
 this->N = vector<unsigned int>(2);
 this->N[0] = Nrows; this->N[1] = Ncols;
 this->d = Matrix<unsigned int>(this->N[0], this->N[1]);
 for (int j = 0; j < this->N[1]; j++)
 {
  for (int i = 0; i < this->N[0]; i++)
  {
   this->d(i, j) = d0;
  }
 }
 this->D = D0;
 unsigned int numTensors = this->N[0] * this->N[1];
 this->Tensors = new Tensor<T>[numTensors];
 vector<unsigned int> Shape(5);
 unsigned int position;
 if (this->BC == "open")
 {
  for (int j = 0; j < this->N[1]; j++)
  {
   for (int i = 0; i < this->N[0]; i++)
   {
    this->getOpenBCShape(i, j, Shape);
    position = i + j*this->N[0];
    this->Tensors[position] = Tensor<T>(Shape);
   }
  }
 }
 else if (this->BC == "periodic")
 {
// we put the same shape everywhere:
  Shape[0] = this->D; Shape[1] = this->D; Shape[2] = this->D; Shape[3] = this->D;
  Shape[4] = this->d(0, 0);
  for (int j = 0; j < this->N[1]; j++)
  {
   for (int i = 0; i < this->N[0]; i++)
   {
    position = i + j*this->N[0];
    this->Tensors[position] = Tensor<T>(Shape);
   }
  }
 }
}

template<class T> PEPS<T>::PEPS(const PEPS<T>& PEPS0)
{
 this->BC = PEPS0.BC;
 this->N = PEPS0.N;
 this->d = PEPS0.d;
 this->D = PEPS0.D;
 unsigned int numTensors = this->N[0] * this->N[1];
 this->Tensors = new Tensor<T>[numTensors];
 for (int i = 0; i < numTensors; i++)
 {
  this->Tensors[i] = PEPS0.Tensors[i];
 }
}

template<class T> PEPS<T>::~PEPS()
{
 delete[] this->Tensors;
}

template<class T> PEPS<T>& PEPS<T>::operator=(const PEPS<T>& PEPS0)
{
 if (this != &PEPS0)
 {
  this->BC = PEPS0.BC;
  this->N = PEPS0.N;
  this->d = PEPS0.d;
  this->D = PEPS0.D;
  delete[] this->Tensors;
  unsigned int numTensors = this->N[0] * this->N[1];
  this->Tensors = new Tensor<T>[numTensors];
  for (int i = 0; i < numTensors; i++)
  {
   this->Tensors[i] = PEPS0.Tensors[i];
  }
 }
 return *this;
}

template<class T> void PEPS<T>::setD(unsigned int D0, T element)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (D0 == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "setD(unsigned int D0, T element = 0.0): " <<
          "((this->N[0] == 0) || (D0 == 0))." << endl;
  exit(1);
 }
#endif
 if (D0 == this->D)
 {
  return;
 }
 unsigned int D1 = this->D;
 T element0;
 vector<unsigned int> Index(5), Shape0(5), Shape1(5);
 Tensor<T> Tensor0, Tensor1;
 this->D = D0;
 if (this->BC == "open")
 {
  for (int col = 0; col < this->N[1]; col++)
  {
   for (int row = 0; row < this->N[0]; row++)
   {
    this->getOpenBCShape(row, col, Shape0);
    Tensor0 = Tensor<T>(Shape0);
    if (element == 0.0)
    {
     Tensor0.fillZeroes();
    }
    else if (element != 0.0)
    {
     Tensor0.fillRandomly(time(0)+(row+this->N[0]*col)*13, element);
    }
    Tensor1 = this->Tensors[row+this->N[0]*col];
    Tensor1.getShape(Shape1);
    for (int j = 0; j < this->d(row, col); j++){
     Index[4] = j;
    for (int i3 = 0; i3 < min(Shape0[3], Shape1[3]); i3++){
     Index[3] = i3;
    for (int i2 = 0; i2 < min(Shape0[2], Shape1[2]); i2++){
     Index[2] = i2;
    for (int i1 = 0; i1 < min(Shape0[1], Shape1[1]); i1++){
     Index[1] = i1;
    for (int i0 = 0; i0 < min(Shape0[0], Shape1[0]); i0++){
     Index[0] = i0;
     element0 = Tensor1.get(Index);
     Tensor0.set(Index, element0);
    }}}}}
    this->Tensors[row+this->N[0]*col] = Tensor0;
   }
  }
 }
 else if (this->BC == "periodic")
 {
  Shape0[0] = D0; Shape0[1] = D0; Shape0[2] = D0; Shape0[3] = D0;
  for (int col = 0; col < this->N[1]; col++)
  {
   for (int row = 0; row < this->N[0]; row++)
   {
    Shape0[4] = this->d(row, col);
    Tensor0 = Tensor<T>(Shape0);
    if (element == 0.0)
    {
     Tensor0.fillZeroes();
    }
    else if (element != 0.0)
    {
     Tensor0.fillRandomly(time(0)+(row+this->N[0]*col)*13, element);
    }
    Tensor1 = this->Tensors[row+this->N[0]*col];
    Tensor1.getShape(Shape1);
    for (int j = 0; j < this->d(row, col); j++){
     Index[4] = j;
    for (int i3 = 0; i3 < min(Shape0[3], Shape1[3]); i3++){
     Index[3] = i3;
    for (int i2 = 0; i2 < min(Shape0[2], Shape1[2]); i2++){
     Index[2] = i2;
    for (int i1 = 0; i1 < min(Shape0[1], Shape1[1]); i1++){
     Index[1] = i1;
    for (int i0 = 0; i0 < min(Shape0[0], Shape1[0]); i0++){
     Index[0] = i0;
     element0 = Tensor1.get(Index);
     Tensor0.set(Index, element0);
    }}}}}
    this->Tensors[row+this->N[0]*col] = Tensor0;
   }
  }
 }
}

template<class T> inline void PEPS<T>::set(unsigned int positionRow, unsigned int positionCol,
                                           const Tensor<T>& Tensor0)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (positionRow >= this->N[0]) || (positionCol >= this->N[1]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> inline void PEPS<T>::" <<
          "set(unsigned int positionRow, unsigned int positionCol, const Tensor<T>& Tensor0): " <<
          "((this->N[0] == 0) || (positionRow >= this->N[0]) || (positionCol >= this->N[1]))." << endl;
  exit(1);
 }
 vector<unsigned int> Shape0;
 Tensor0.getShape(Shape0);
 if (this->BC == "open")
 {
  vector<unsigned int> Shape1(5);
  this->getOpenBCShape(positionRow, positionCol, Shape1);
  if (Shape0 != Shape1)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> inline void PEPS<T>::" <<
           "set(unsigned int positionRow, unsigned int positionCol, const Tensor<T>& Tensor0): " <<
           "Tensor0 has an incorrect shape." << endl;
   exit(1);
  }
 }
 else if (this->BC == "periodic")
 {
  if ((Shape0[0] != this->D) || (Shape0[1] != this->D) || (Shape0[2] != this->D) ||
      (Shape0[3] != this->D) || (Shape0[4] != this->d(positionRow, positionCol)))
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> inline void PEPS<T>::" <<
           "set(unsigned int positionRow, unsigned int positionCol, const Tensor<T>& Tensor0): " <<
           "Tensor0 has an incorrect shape." << endl;
   exit(1);
  }
 }
#endif
 this->Tensors[positionRow + this->N[0]*positionCol] = Tensor0;
}

template<class T> inline void PEPS<T>::get(unsigned int positionRow, unsigned int positionCol,
                                           Tensor<T>& Tensor0) const
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (positionRow >= this->N[0]) || (positionCol >= this->N[1]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> inline void PEPS<T>::" <<
          "get(unsigned int positionRow, unsigned int positionCol, Tensor<T>& Tensor0) const: " <<
          "((this->N[0] == 0) || (positionRow >= this->N[0]) || (positionCol >= this->N[1]))." << endl;
  exit(1);
 }
#endif
 Tensor0 = this->Tensors[positionRow + this->N[0]*positionCol];
}

template<class T> void PEPS<T>::write(const string& FileName) const
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
  File.write((char*)&(this->N[0]), sizeof(this->N[0]));
  File.write((char*)&(this->N[1]), sizeof(this->N[1]));
  unsigned int size = this->N[0]*this->N[1];
  unsigned int elementd;
  for (int i = 0; i < size; i++)
  {
   elementd = this->d.get(i);
   File.write((char*)&elementd, sizeof(elementd));
  }
  File.write((char*)&(this->D), sizeof(this->D));
  Tensor<T> Tensor0;
  unsigned int rank0, size0;
  vector<unsigned int> Shape0;
  T element0;
  for (int i = 0; i < size; i++)
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
          "template<class T> void PEPS<T>::" <<
          "write(const string& FileName) const: " <<
          "Binary output file " << FileName << " could not be opened." << endl;
  exit(1);
 }
}

template<class T> void PEPS<T>::read(const string& FileName)
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
  File.read((char*)&(this->N[0]), sizeof(this->N[0]));
  File.read((char*)&(this->N[1]), sizeof(this->N[1]));
  unsigned int size = this->N[0]*this->N[1];
  this->d = Matrix<unsigned int>(this->N[0], this->N[1]);
  unsigned int elementd;
  for (int i = 0; i < size; i++)
  {
   File.read((char*)&elementd, sizeof(elementd));
   d.set(i, elementd);
  }
  File.read((char*)&(this->D), sizeof(this->D));
  delete[] this->Tensors;
  this->Tensors = new Tensor<T>[size];
  unsigned int rank0, size0;
  vector<unsigned int> Shape0;
  Tensor<T> Tensor0;
  T element0;
  for (int i = 0; i < size; i++)
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
          "template<class T> void PEPS<T>::" <<
          "read(const string& FileName): " <<
          "Binary input file " << FileName << " could not be opened." << endl;
  exit(1);
 }
}

template<class T> void PEPS<T>::readValentinPEPS(const string& FileName)
{
 ifstream File(FileName.c_str(), ios::in);
 if (File.is_open())
 {
  vector<unsigned int> Shape0(5), Index0(5);
  Tensor<T> Tensor0;
  T element;
  for (int col = 0; col < this->N[1]; col++)
  {
   for (int row = 0; row < this->N[0]; row++)
   {
    this->getOpenBCShape(row, col, Shape0);
    Tensor0 = Tensor<T>(Shape0);
    for (int i0 = 0; i0 < Shape0[0]; i0++)
    {
     Index0[0] = i0;
     for (int i1 = 0; i1 < Shape0[1]; i1++)
     {
      Index0[1] = i1;
      for (int i2 = 0; i2 < Shape0[2]; i2++)
      {
       Index0[2] = i2;
       for (int i3 = 0; i3 < Shape0[3]; i3++)
       {
        Index0[3] = i3;
        for (int j = 0; j < Shape0[4]; j++)
        {
         Index0[4] = j;
         File >> element;
         Tensor0.set(Index0, element);
        }
       }
      }
     }
    }
    this->Tensors[row+this->N[0]*col] = Tensor0;
   }
  }
  File.close();
 }
 else
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "readValentinPEPS(const string& FileName): " <<
          "Text input file " << FileName << " could not be opened." << endl;
  exit(1);
 }
}

template<class T> void PEPS<T>::convertToComplex(PEPS< complex<float> >& PEPS0) const
{
#ifdef DEBUG
 if (typeid(T) != typeid(float))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "convertToComplex(PEPS< complex<float> >& PEPS0) const: " <<
          "(typeid(T) != typeid(float))." << endl;
  exit(1);
 }
#endif
 PEPS0 = PEPS< complex<float> >(this->BC, this->N[0], this->N[1], this->d, this->D);
 Tensor< complex<float> > Tensor0;
 for (int col = 0; col < this->N[1]; col++)
 {
  for (int row = 0; row < this->N[0]; row++)
  {
   this->Tensors[row+this->N[0]*col].convertToComplex(Tensor0);
   PEPS0.set(row, col, Tensor0);
  }
 }
}

template<class T> void PEPS<T>::convertToComplex(PEPS< complex<double> >& PEPS0) const
{
#ifdef DEBUG
 if (typeid(T) != typeid(double))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "convertToComplex(PEPS< complex<double> >& PEPS0) const: " <<
          "(typeid(T) != typeid(double))." << endl;
  exit(1);
 }
#endif
 PEPS0 = PEPS< complex<double> >(this->BC, this->N[0], this->N[1], this->d, this->D);
 Tensor< complex<double> > Tensor0;
 for (int col = 0; col < this->N[1]; col++)
 {
  for (int row = 0; row < this->N[0]; row++)
  {
   this->Tensors[row+this->N[0]*col].convertToComplex(Tensor0);
   PEPS0.set(row, col, Tensor0);
  }
 }
}

template<class T> void PEPS<T>::getVector(vector<T>& Vector0) const
{}

template<class T> void PEPS<T>::fillRandomly(const Matrix<unsigned int>& Seed, T element)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (Seed.getDim0() != this->N[0]) || (Seed.getDim1() != this->N[1]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "fillRandomly(const Matrix<unsigned int>& Seed, T element): " <<
          "((this->N[0] == 0) || (Seed.getDim0() != this->N[0]) || " <<
           "(Seed.getDim1() != this->N[1]))." << endl;
  exit(1);
 }
#endif
 for (int col = 0; col < this->N[1]; col++)
 {
  for (int row = 0; row < this->N[0]; row++)
  {
   this->Tensors[row+this->N[0]*col].fillRandomly(Seed(row, col), element);
  }
 }
}

template<class T> void PEPS<T>::fillZeroes()
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "fillZeroes(): " <<
          "(this->N[0] == 0)." << endl;
  exit(1);
 }
#endif
 unsigned int numTensors = this->N[0] * this->N[1];
 for (int i = 0; i < numTensors; i++)
 {
  this->Tensors[i].fillZeroes();
 }
}

template<class T> void PEPS<T>::setSeparable(const vector<T>& Coefficients, T element)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (Coefficients.size() == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "setSeparable(const vector<T>& Coefficients, T element): " <<
          "((this->N[0] == 0) || (Coefficients.size() == 0))." << endl;
  exit(1);
 }
#endif
 if (element == 0.0)
 {
  this->fillZeroes();
 }
 else
 {
  Matrix<unsigned int> Seed(this->N[0], this->N[1]);
  for (int col = 0; col < this->N[1]; col++)
  {
   for (int row = 0; row < this->N[0]; row++)
   {
    Seed(row, col) = 13+(row+this->N[0]*col)*2;
   }
  }
  this->fillRandomly(Seed, element);
 }
 vector<unsigned int> Index(5);
 Index[0] = 0; Index[1] = 0; Index[2] = 0; Index[3] = 0;
 for (int col = 0; col < this->N[1]; col++)
 {
  for (int row = 0; row < this->N[0]; row++)
  {
   for (int j = 0; j < this->d(row, col); j++)
   {
    if (j < Coefficients.size())
    {
     Index[4] = j;
     this->Tensors[row+this->N[0]*col].set(Index, Coefficients[j]);
    }
    else
    {
     break;
    }
   }
  }
 }
}

template<class T> void PEPS<T>::setConcatenated(unsigned int M, T element)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (M == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "setConcatenated(unsigned int M, T element): " <<
          "((this->N[0] == 0) || (M == 0))." << endl;
  exit(1);
 }
#endif
 if (M == 1)
  return;
 unsigned int positionCol, positionRow;
 vector<unsigned int> Index(5), N0, Shape(5);
 Tensor<T> Tensor0;
 Tensor<T>* TensorsConc;
 Matrix<unsigned int> d0;
// 1. get N and d of Concatenated PEPS:
 N0 = this->N;
 this->N[0] = N0[0]*M;
 this->N[1] = N0[1]*M;
 d0 = this->d;
 this->d = Matrix<unsigned int>(this->N[0], this->N[1]);
 for (int col0 = 0; col0 < N0[1]; col0++){
  for (int col1 = 0; col1 < M; col1++){
   positionCol = col0*M+col1;
   for (int row0 = 0; row0 < N0[0]; row0++){
    for (int row1 = 0; row1 < M; row1++){
     positionRow = row0*M+row1;
// 1.a) ((row1 == 0) && (col1 == 0)): physical tensors with this physical dimension d
     if ((row1 == 0) && (col1 == 0))
      this->d(positionRow, positionCol) = d0(row0, col0);
// 1.b) !((row1 == 0) && (col1 == 0)): auxiliary tensors with physical dimension 1
     else
      this->d(positionRow, positionCol) = 1;
 }}}}
// 2. get tensors of Concatenated PEPS:
 TensorsConc = new Tensor<T>[this->N[0]*this->N[1]];
 Shape[0] = this->D; Shape[1] = this->D; Shape[2] = this->D; Shape[3] = this->D;
 for (int col0 = 0; col0 < N0[1]; col0++){
  for (int col1 = 0; col1 < M; col1++){
   positionCol = col0*M+col1;
   for (int row0 = 0; row0 < N0[0]; row0++){
    for (int row1 = 0; row1 < M; row1++){
     positionRow = row0*M+row1;
     if (this->BC == "open")
      this->getOpenBCShape(positionRow, positionCol, Shape);
     else if (this->BC == "periodic")
      Shape[4] = this->d(positionRow, positionCol);
// 2.a) ((row1 == 0) && (col1 == 0)): physical tensors taken from this PEPS
     if ((row1 == 0) && (col1 == 0))
     {
      TensorsConc[positionRow+this->N[0]*positionCol] = this->Tensors[row0+N0[0]*col0];
      TensorsConc[positionRow+this->N[0]*positionCol].setShape(Shape, element, time(0)+(positionRow+this->N[0]*positionCol)*13);
     }
// 2.b) !((row1 == 0) && (col1 == 0)): auxiliary tensors set as Deltas
     else
     {
      Tensor0 = Tensor<T>(Shape);
      if (element == 0.0)
       Tensor0.fillZeroes();
      else
       Tensor0.fillRandomly(time(0)+(positionRow+this->N[0]*positionCol)*13, element);
      Index[4] = 0;
// 2.b)a) ((row1 == 0) && (col1 != 0)): auxiliary tensor set as Delta in horizontal direction
      if (row1 == 0)
      {
       Index[1] = 0; Index[3] = 0;
       for (int i = 0; i < min(Shape[0], Shape[2]); i++){
        Index[0] = i; Index[2] = i;
        Tensor0.set(Index, 1.0);
       }
      }
// 2.b)b) ((row1 != 0) && (col1 == 0)): auxiliary tensor set as Delta in vertical direction
      else if (col1 == 0)
      {
       Index[0] = 0; Index[2] = 0;
       for (int i = 0; i < min(Shape[1], Shape[3]); i++){
        Index[1] = i; Index[3] = i;
        Tensor0.set(Index, 1.0);
       }
      }
// 2.b)c) ((row1 != 0) && (col1 != 0)): auxiliary tensor set as scalar with value 1
      else
      {
       Index[0] = 0; Index[1] = 0; Index[2] = 0; Index[3] = 0;
       Tensor0.set(Index, 1.0);
      }
      TensorsConc[positionRow+this->N[0]*positionCol] = Tensor0;
     }
 }}}}
// 3. write Concatenated PEPS into this PEPS:
 delete[] this->Tensors;
 this->Tensors = TensorsConc;
}

template<class T> void PEPS<T>::multiply(T element)
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiply(T element): " <<
          "(this->N[0] == 0)." << endl;
  exit(1);
 }
#endif
 this->Tensors[0].multiply(element);
}

template<class T> void PEPS<T>::decanonicalize(const Matrix< Matrix<T> >& LambdasV,
                                               const Matrix< Matrix<T> >& LambdasH)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (LambdasV.getDim0() != this->N[0]-1) || (LambdasV.getDim1() != this->N[1]) ||
     (LambdasH.getDim0() != this->N[0]) || (LambdasH.getDim1() != this->N[1]-1))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "decanonicalize(const Matrix< Matrix<T> >& LambdasV, " <<
                         "const Matrix< Matrix<T> >& LambdasH): " <<
          "((this->N[0] == 0) || (LambdasV.getDim0() != this->N[0]-1) || (LambdasV.getDim1() != this->N[1]) || " <<
           "(LambdasH.getDim0() != this->N[0]) || (LambdasH.getDim1() != this->N[1]-1))." << endl;
  exit(1);
 }
 for (int col = 0; col < this->N[1]; col++)
 {
  for (int row = 0; row < this->N[0]-1; row++)
  {
   if ((LambdasV(row, col).getDim0() != this->D) || (LambdasV(row, col).getDim1() != this->D))
   {
    cerr << "Program terminated because of error in function " <<
            "template<class T> void PEPS<T>::" <<
            "decanonicalize(const Matrix< Matrix<T> >& LambdasV, " <<
                           "const Matrix< Matrix<T> >& LambdasH): " <<
            "LambdasV(" << row << ", " << col << ") has incorrect shape." << endl;
    exit(1);
   }
  }
 }
 for (int col = 0; col < this->N[1]-1; col++)
 {
  for (int row = 0; row < this->N[0]; row++)
  {
   if ((LambdasH(row, col).getDim0() != this->D) || (LambdasH(row, col).getDim1() != this->D))
   {
    cerr << "Program terminated because of error in function " <<
            "template<class T> void PEPS<T>::" <<
            "decanonicalize(const Matrix< Matrix<T> >& LambdasV, " <<
                           "const Matrix< Matrix<T> >& LambdasH): " <<
            "LambdasH(" << row << ", " << col << ") has incorrect shape." << endl;
    exit(1);
   }
  }
 }
#endif
 vector<unsigned int> Indices0(1), Indices1(1), Order(5);
 Matrix<T> SquareRootLambda(this->D, this->D);
 for (int col = 0; col < this->N[1]; col++)
 {
  for (int row = 0; row < this->N[0]; row++)
  {
   if (col != 0)
   {
    Indices0[0] = 0; Indices1[0] = 1;
    SquareRootLambda.fillZeroes();
    for (int i = 0; i < this->D; i++)
     SquareRootLambda(i, i) = sqrt(LambdasH(row, col-1)(i, i));
    this->Tensors[row+this->N[0]*col].contract(Indices0, SquareRootLambda, Indices1);
    Order[0] = 4; Order[1] = 0; Order[2] = 1; Order[3] = 2; Order[4] = 3;
    this->Tensors[row+this->N[0]*col].permute(Order);
   }
   if (row != 0)
   {
    Indices0[0] = 1; Indices1[0] = 1;
    SquareRootLambda.fillZeroes();
    for (int i = 0; i < this->D; i++)
     SquareRootLambda(i, i) = sqrt(LambdasV(row-1, col)(i, i));
    this->Tensors[row+this->N[0]*col].contract(Indices0, SquareRootLambda, Indices1);
    Order[0] = 0; Order[1] = 4; Order[2] = 1; Order[3] = 2; Order[4] = 3;
    this->Tensors[row+this->N[0]*col].permute(Order);
   }
   if (col != this->N[1]-1)
   {
    Indices0[0] = 2; Indices1[0] = 0;
    SquareRootLambda.fillZeroes();
    for (int i = 0; i < this->D; i++)
     SquareRootLambda(i, i) = sqrt(LambdasH(row, col)(i, i));
    this->Tensors[row+this->N[0]*col].contract(Indices0, SquareRootLambda, Indices1);
    Order[0] = 0; Order[1] = 1; Order[2] = 4; Order[3] = 2; Order[4] = 3;
    this->Tensors[row+this->N[0]*col].permute(Order);
   }
   if (row != this->N[0]-1)
   {
    Indices0[0] = 3; Indices1[0] = 0;
    SquareRootLambda.fillZeroes();
    for (int i = 0; i < this->D; i++)
     SquareRootLambda(i, i) = sqrt(LambdasV(row, col)(i, i));
    this->Tensors[row+this->N[0]*col].contract(Indices0, SquareRootLambda, Indices1);
    Order[0] = 0; Order[1] = 1; Order[2] = 2; Order[3] = 4; Order[4] = 3;
    this->Tensors[row+this->N[0]*col].permute(Order);
   }
  }
 }
}

template<class T> T PEPS<T>::normalize()
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "normalize(): " <<
          "(this->N[0] == 0)." << endl;
  exit(1);
 }
#endif
 T result = 1.0;
 for (int col = 0; col < this->N[1]; col++)
 {
  for (int row = 0; row < this->N[0]; row++)
  {
   result *= this->Tensors[row+this->N[0]*col].normalize();
  }
 }
 return result;
}

template<class T> void PEPS<T>::getMPS(const string& Direction, MPS<T>& MPS0) const
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (this->BC != "open") ||
     ((Direction != "right") && (Direction != "down") && (Direction != "left") && (Direction != "up")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getMPS(const string& Direction, MPS<T>& MPS0) const: " <<
          "((this->N[0] == 0) || (this->BC != open) || " <<
           "((Direction != right) && (Direction != down) && (Direction != left) && (Direction != up)))." << endl;
  exit(1);
 }
#endif
 string BCMPS0 = "open";
 unsigned int NMPS0;
 unsigned int dMPS0 = this->D*this->D, DMPS0 = this->D*this->D;
 vector<unsigned int> Indices0(1), Indices1(1);
 Indices0[0] = 4; Indices1[0] = 4;
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Shape0(3), Order0(8);
 T element;
 if (Direction == "right")
 {
  NMPS0 = this->N[0];
  MPS0 = MPS<T>(BCMPS0, NMPS0, dMPS0, DMPS0);
  Order0[0] = 3; Order0[1] = 7; Order0[2] = 0; Order0[3] = 4;
  Order0[4] = 1; Order0[5] = 5; Order0[6] = 2; Order0[7] = 6;
// the left end of MPS0:
  Shape0[0] = 1; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  Tensor0 = this->Tensors[this->N[0]-1];
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(0, Tensor0);
// the intermediate part of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  for (int i = 1; i < NMPS0-1; i++)
  {
   Tensor0 = this->Tensors[this->N[0]-1-i];
   Tensor0.complexConjugate(Tensor1);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPS0.set(i, Tensor0);
  }
// the right end of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = 1; Shape0[2] = dMPS0;
  Tensor0 = this->Tensors[0];
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(NMPS0-1, Tensor0);
 }
 else if (Direction == "down")
 {
  NMPS0 = this->N[1];
  MPS0 = MPS<T>(BCMPS0, NMPS0, dMPS0, DMPS0);
  Order0[0] = 0; Order0[1] = 4; Order0[2] = 1; Order0[3] = 5;
  Order0[4] = 2; Order0[5] = 6; Order0[6] = 3; Order0[7] = 7;
// the left end of MPS0:
  Shape0[0] = 1; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  Tensor0 = this->Tensors[0];
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(0, Tensor0);
// the intermediate part of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  for (int i = 1; i < NMPS0-1; i++)
  {
   Tensor0 = this->Tensors[this->N[0]*i];
   Tensor0.complexConjugate(Tensor1);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPS0.set(i, Tensor0);
  }
// the right end of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = 1; Shape0[2] = dMPS0;
  Tensor0 = this->Tensors[this->N[0]*(this->N[1]-1)];
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(NMPS0-1, Tensor0);
 }
 else if (Direction == "left")
 {
  NMPS0 = this->N[0];
  MPS0 = MPS<T>(BCMPS0, NMPS0, dMPS0, DMPS0);
  Order0[0] = 1; Order0[1] = 5; Order0[2] = 2; Order0[3] = 6;
  Order0[4] = 3; Order0[5] = 7; Order0[6] = 0; Order0[7] = 4;
// the left end of MPS0:
  Shape0[0] = 1; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  Tensor0 = this->Tensors[this->N[0]*(this->N[1]-1)];
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(0, Tensor0);
// the intermediate part of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  for (int i = 1; i < NMPS0-1; i++)
  {
   Tensor0 = this->Tensors[i+this->N[0]*(this->N[1]-1)];
   Tensor0.complexConjugate(Tensor1);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPS0.set(i, Tensor0);
  }
// the right end of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = 1; Shape0[2] = dMPS0;
  Tensor0 = this->Tensors[this->N[0]*this->N[1]-1];
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(NMPS0-1, Tensor0);
 }
 else if (Direction == "up")
 {
  NMPS0 = this->N[1];
  MPS0 = MPS<T>(BCMPS0, NMPS0, dMPS0, DMPS0);
  Order0[0] = 2; Order0[1] = 6; Order0[2] = 3; Order0[3] = 7;
  Order0[4] = 0; Order0[5] = 4; Order0[6] = 1; Order0[7] = 5;
// the left end of MPS0:
  Shape0[0] = 1; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  Tensor0 = this->Tensors[this->N[0]*this->N[1]-1];
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(0, Tensor0);
// the intermediate part of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  for (int i = 1; i < NMPS0-1; i++)
  {
   Tensor0 = this->Tensors[this->N[0]*(this->N[1]-i)-1];
   Tensor0.complexConjugate(Tensor1);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPS0.set(i, Tensor0);
  }
// the right end of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = 1; Shape0[2] = dMPS0;
  Tensor0 = this->Tensors[this->N[0]-1];
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(NMPS0-1, Tensor0);
 }
}

template<class T> void PEPS<T>::getBoundaryMPS(const string& Boundary,
                                               double eps, unsigned int maxNumSweeps,
                                               double& epsAchieved, unsigned int& numSweepsDone,
                                               MPS<T>& MPS0) const
{
#ifdef DEBUG
 if ((this->BC != "open") || (this->N[0] == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getBoundaryMPS(const string& Boundary, double eps, unsigned int maxNumSweeps, " <<
                         "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "((this->BC != open) || (this->N[0] == 0))." << endl;
  exit(1);
 }
 if ((Boundary != "left") && (Boundary != "top") && (Boundary != "right") && (Boundary != "bottom"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getBoundaryMPS(const string& Boundary, double eps, unsigned int maxNumSweeps, " <<
                         "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "((Boundary != left) && (Boundary != top) && (Boundary != right) && (Boundary != bottom))." << endl;
  exit(1);
 }
 string BCMPS0; MPS0.getBC(BCMPS0);
 if ((BCMPS0 != "open") ||
     (((Boundary == "left") || (Boundary == "right")) && (MPS0.getN() != this->N[0])) ||
     (((Boundary == "top") || (Boundary == "bottom")) && (MPS0.getN() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getBoundaryMPS(const string& Boundary, double eps, unsigned int maxNumSweeps, " <<
                         "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "((MPS0.BC != open) || (((Boundary == left) || (Boundary == right)) && (MPS0.N != this->N[0])) || " <<
           "(((Boundary == top) || (Boundary == bottom)) && (MPS0.N != this->N[1])))." << endl;
  exit(1);
 }
#endif
 unsigned int DMPS0, L;
 double error0, error1;
 T bValue, NValue, tensorFrobeniusNorm, tensorFrobeniusNormProduct;
 string SweepDirection;
 vector<unsigned int> dMPS0, Indices0(1), Indices02(2), Indices03(3), Indices1(1), Indices12(2), Indices13(3), Indices2(2);
 vector<unsigned int> Order3(3), Order5(5), Order6(6), Shape3(3), Shape4(4), Shape5(5), Shape6(6);
 Tensor<T> bTensor, Tensor0, Tensor1, TensorL, TensorR;
 vector< Tensor<T> > BoundaryTensors, bTensors;
// 0. initialize:
// 0. - set boundary-length L:
 if ((Boundary == "left") || (Boundary == "right"))
  L = this->N[0];
 else if ((Boundary == "top") || (Boundary == "bottom"))
  L = this->N[1];
// 0. - get boundary-tensors clockwise, permute them as "top" boundary, and write them in BoundaryTensors:
 BoundaryTensors = vector< Tensor<T> >(L);
 if (Boundary == "left")
 {
  Order5[0] = 3; Order5[1] = 0; Order5[2] = 1; Order5[3] = 2; Order5[4] = 4;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[this->N[0]-1-pos];
   Tensor0.permute(Order5);
   BoundaryTensors[pos] = Tensor0;
  }
 }
 else if (Boundary == "top")
 {
  for (int pos = 0; pos < L; pos++)
  {
   BoundaryTensors[pos] = this->Tensors[this->N[0]*pos];
  }
 }
 else if (Boundary == "right")
 {
  Order5[0] = 1; Order5[1] = 2; Order5[2] = 3; Order5[3] = 0; Order5[4] = 4;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[pos+this->N[0]*(this->N[1]-1)];
   Tensor0.permute(Order5);
   BoundaryTensors[pos] = Tensor0;
  }
 }
 else if (Boundary == "bottom")
 {
  Order5[0] = 2; Order5[1] = 3; Order5[2] = 0; Order5[3] = 1; Order5[4] = 4;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[this->N[0]-1+this->N[0]*(this->N[1]-1-pos)];
   Tensor0.permute(Order5);
   BoundaryTensors[pos] = Tensor0;
  }
 }
// 0. - reshape boundary-tensors eliminating index 1:
 for (int pos = 0; pos < L; pos++)
 {
  BoundaryTensors[pos].getShape(Shape5);
  Shape4[0] = Shape5[0]; Shape4[1] = Shape5[2]; Shape4[2] = Shape5[3]; Shape4[3] = Shape5[4];
  BoundaryTensors[pos].reshape(Shape4);
 }
// 0. - check MPS0.d:
 MPS0.getd(dMPS0);
#ifdef DEBUG
 for (int pos = 0; pos < L; pos++)
 {
  BoundaryTensors[pos].getShape(Shape4);
  if (dMPS0[pos] != Shape4[2]*Shape4[2])
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void PEPS<T>::" <<
           "getBoundaryMPS(const string& Boundary, double eps, unsigned int maxNumSweeps, " <<
                          "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
           "(MPS0.d[" << pos << "] != this->D[" << pos << "]*this->D[" << pos << "])." << endl;
   exit(1);
  }
 }
#endif
// 1. MPS0.D >= this->D*this->D:
 DMPS0 = MPS0.getD();
 if (DMPS0 >= this->D*this->D)
 {
  MPS0 = MPS<T>("open", L, dMPS0, this->D*this->D);
  Indices0[0] = 3; Indices1[0] = 3;
  Order6[0] = 0; Order6[1] = 3; Order6[2] = 1; Order6[3] = 4; Order6[4] = 2; Order6[5] = 5;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = BoundaryTensors[pos];
   Tensor0.complexConjugate(Tensor1);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order6);
   Tensor0.getShape(Shape6);
   Shape3[0] = Shape6[0]*Shape6[1]; Shape3[1] = Shape6[2]*Shape6[3]; Shape3[2] = Shape6[4]*Shape6[5];
   Tensor0.reshape(Shape3);
   MPS0.set(pos, Tensor0);
  }
  MPS0.bringIntoNormalShape();
  MPS0.setD(DMPS0);
  epsAchieved = 0.0;
  numSweepsDone = 0;
  return;
 }
// 2. bring MPS0 into normal form from right to left:
// 2.a) MPS0.D == 1:
 if (DMPS0 == 1)
 {
  tensorFrobeniusNormProduct = 1.0;
  for (int pos = L-1; pos > 0; pos--)
  {
   MPS0.get(pos, Tensor0);
   tensorFrobeniusNormProduct *= Tensor0.frobeniusNormalize();
   MPS0.set(pos, Tensor0);
  }
  MPS0.get(0, Tensor0);
  Tensor0.multiply(tensorFrobeniusNormProduct);
  MPS0.set(0, Tensor0);
 }
// 2.b) MPS0.D != 1:
 else
 {
  Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
  Indices0[0] = 1;
  Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
  MPS0.get(L-1, Tensor0);
  for (int pos = L-1; pos > 0; pos--)
  {
   Tensor0.LQDecompose(Indices1, Indices2, TensorL);
   MPS0.set(pos, Tensor0);
   MPS0.get(pos-1, Tensor0);
   Tensor0.contract(Indices0, TensorL, Indices1);
   Tensor0.permute(Order3);
  }
  MPS0.set(0, Tensor0);
 }
// 3. construct right half of b for each site:
 bTensors = vector< Tensor<T> >(L);
 Tensor0 = BoundaryTensors[L-1];
 Tensor0.complexConjugate(Tensor1);
 Indices0[0] = 3; Indices1[0] = 3;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor0.getShape(Shape6);
 MPS0.get(L-1, Tensor1);
 Tensor1.getShape(Shape3);
 Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape6[2]; Shape4[3] = Shape6[5];
 Tensor1.reshape(Shape4);
 Tensor1.complexConjugate();
 Indices02[0] = 2; Indices02[1] = 5; Indices12[0] = 2; Indices12[1] = 3;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 Tensor0.getShape(Shape6);
 Shape3[0] = Shape6[0]; Shape3[1] = Shape6[2]; Shape3[2] = Shape6[4];
 Tensor0.reshape(Shape3);
 bTensors[L-1] = Tensor0;
 for (int pos = L-2; pos > 0; pos--)
 {
  Tensor1 = BoundaryTensors[pos];
  Indices0[0] = 0; Indices1[0] = 1;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  BoundaryTensors[pos].complexConjugate(Tensor1);
  Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 3;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape5);
  MPS0.get(pos, Tensor1);
  Tensor1.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[4];
  Tensor1.reshape(Shape4);
  Tensor1.complexConjugate();
  Indices03[0] = 0; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  bTensors[pos] = Tensor0;
 }
// 4. compute initial error:
 Tensor1 = BoundaryTensors[0];
 Indices0[0] = 0; Indices1[0] = 1;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 BoundaryTensors[0].complexConjugate(Tensor1);
 Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 3;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 Tensor0.getShape(Shape5);
 MPS0.get(0, Tensor1);
 Tensor1.getShape(Shape3);
 Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[4];
 Tensor1.reshape(Shape4);
 Tensor1.complexConjugate();
 Indices03[0] = 0; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
 Tensor0.contract(Indices03, Tensor1, Indices13);
 bValue = Tensor0.get(0);
 MPS0.get(0, Tensor0);
 NValue = Tensor0.scalarProduct();
 error0 = MathAuxiliary::convertToDouble(NValue)-2.0*MathAuxiliary::convertToDouble(bValue);
// 5. sweep:
 numSweepsDone = 0;
 SweepDirection = "right";
 while (numSweepsDone < maxNumSweeps)
 {
// 5.a) sweep to the right:
  if (SweepDirection == "right")
  {
// 5.a)1. position 0:
// 5.a)1. - compute bTensor:
   Tensor0 = bTensors[1];
   Tensor1 = BoundaryTensors[0];
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   BoundaryTensors[0].complexConjugate(Tensor1);
   Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
// 5.a)1. - update tensor:
   Tensor0.getShape(Shape5);
   Shape3[0] = 1; Shape3[1] = Shape5[0]; Shape3[2] = Shape5[2]*Shape5[4];
   Tensor0.reshape(Shape3);
   MPS0.set(0, Tensor0);
// 5.a)1. - compute error:
   NValue = Tensor0.scalarProduct();
   error1 = -MathAuxiliary::convertToDouble(NValue);
   epsAchieved = abs((error1-error0)/error1);
   if (epsAchieved < eps)
    return;
   error0 = error1;
// 5.a)1. - normalize tensor:
   if (DMPS0 == 1)
    tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
   else
   {
    Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
    Tensor0.QRDecompose(Indices2, Indices1, TensorR);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
   }
   MPS0.set(0, Tensor0);
// 5.a)1. - set bTensors:
   Tensor0 = BoundaryTensors[0];
   Tensor0.complexConjugate(Tensor1);
   Indices0[0] = 3; Indices1[0] = 3;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.getShape(Shape6);
   MPS0.get(0, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape6[2]; Shape4[3] = Shape6[5];
   Tensor1.reshape(Shape4);
   Tensor1.complexConjugate();
   Indices02[0] = 2; Indices02[1] = 5; Indices12[0] = 2; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape6);
   Shape3[0] = Shape6[1]; Shape3[1] = Shape6[3]; Shape3[2] = Shape6[5];
   Tensor0.reshape(Shape3);
   bTensors[0] = Tensor0;
// 5.a)2. 0 < position < L-1:
   for (int pos = 1; pos < L-1; pos++)
   {
// 5.a)2. - compute bTensor:
    Tensor0 = bTensors[pos-1];
    Tensor1 = BoundaryTensors[pos];
    Indices0[0] = 0; Indices1[0] = 0;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    BoundaryTensors[pos].complexConjugate(Tensor1);
    Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 3;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    bTensor = Tensor0;
// 5.a)2. - update tensor:
    Tensor1 = bTensors[pos+1];
    Indices02[0] = 1; Indices02[1] = 3; Indices12[0] = 0; Indices12[1] = 1;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    Tensor0.getShape(Shape4);
    Shape3[0] = Shape4[0]; Shape3[1] = Shape4[1]*Shape4[2]; Shape3[2] = Shape4[3];
    Tensor0.reshape(Shape3);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
    MPS0.set(pos, Tensor0);
// 5.a)2. - compute error:
    NValue = Tensor0.scalarProduct();
    error1 = -MathAuxiliary::convertToDouble(NValue);
    epsAchieved = abs((error1-error0)/error1);
    if (epsAchieved < eps)
     return;
    error0 = error1;
// 5.a)2. - normalize tensor:
    if (DMPS0 == 1)
     tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
    else
    {
     Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
     Tensor0.QRDecompose(Indices2, Indices1, TensorR);
     Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
     Tensor0.permute(Order3);
    }
    MPS0.set(pos, Tensor0);
// 5.a)2. - set bTensors:
    bTensor.getShape(Shape5);
    Tensor0.getShape(Shape3);
    Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[4];
    Tensor0.reshape(Shape4);
    Tensor0.complexConjugate();
    Indices03[0] = 0; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
    bTensor.contract(Indices03, Tensor0, Indices13);
    bTensors[pos] = bTensor;
   }
// 5.a)3. position L-1:
   MPS0.get(L-1, Tensor0);
   if (DMPS0 == 1)
    Tensor0.multiply(tensorFrobeniusNorm);
   else
   {
    Indices0[0] = 1; Indices1[0] = 0;
    TensorR.contract(Indices0, Tensor0, Indices1);
    Tensor0 = TensorR;
   }
   MPS0.set(L-1, Tensor0);
   SweepDirection = "left";
  }
// 5.b) sweep to the left:
  else if (SweepDirection == "left")
  {
// 5.b)1. position L-1:
// 5.b)1. - compute bTensor:
   Tensor0 = bTensors[L-2];
   Tensor1 = BoundaryTensors[L-1];
   Indices0[0] = 0; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   BoundaryTensors[L-1].complexConjugate(Tensor1);
   Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
// 5.b)1. - update tensor:
   Tensor0.getShape(Shape5);
   Shape3[0] = Shape5[0]; Shape3[1] = 1; Shape3[2] = Shape5[2]*Shape5[4];
   Tensor0.reshape(Shape3);
   MPS0.set(L-1, Tensor0);
// 5.b)1. - compute error:
   NValue = Tensor0.scalarProduct();
   error1 = -MathAuxiliary::convertToDouble(NValue);
   epsAchieved = abs((error1-error0)/error1);
   if (epsAchieved < eps)
    return;
   error0 = error1;
// 5.b)1. - normalize tensor:
   if (DMPS0 == 1)
    tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
   else
   {
    Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
    Tensor0.LQDecompose(Indices1, Indices2, TensorL);
   }
   MPS0.set(L-1, Tensor0);
// 5.b)1. - set bTensors:
   Tensor0 = BoundaryTensors[L-1];
   Tensor0.complexConjugate(Tensor1);
   Indices0[0] = 3; Indices1[0] = 3;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.getShape(Shape6);
   MPS0.get(L-1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape6[2]; Shape4[3] = Shape6[5];
   Tensor1.reshape(Shape4);
   Tensor1.complexConjugate();
   Indices02[0] = 2; Indices02[1] = 5; Indices12[0] = 2; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape6);
   Shape3[0] = Shape6[0]; Shape3[1] = Shape6[2]; Shape3[2] = Shape6[4];
   Tensor0.reshape(Shape3);
   bTensors[L-1] = Tensor0;
// 5.b)2. 0 < position < L-1:
   for (int pos = L-2; pos > 0; pos--)
   {
// 5.b)2. - compute bTensor:
    Tensor0 = bTensors[pos+1];
    Tensor1 = BoundaryTensors[pos];
    Indices0[0] = 0; Indices1[0] = 1;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    BoundaryTensors[pos].complexConjugate(Tensor1);
    Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 3;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    bTensor = Tensor0;
// 5.b)2. - update tensor:
    Tensor1 = bTensors[pos-1];
    Indices02[0] = 1; Indices02[1] = 3; Indices12[0] = 0; Indices12[1] = 1;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    Tensor0.getShape(Shape4);
    Shape3[0] = Shape4[0]; Shape3[1] = Shape4[1]*Shape4[2]; Shape3[2] = Shape4[3];
    Tensor0.reshape(Shape3);
    Order3[0] = 2; Order3[1] = 0; Order3[2] = 1;
    Tensor0.permute(Order3);
    MPS0.set(pos, Tensor0);
// 5.b)2. - compute error:
    NValue = Tensor0.scalarProduct();
    error1 = -MathAuxiliary::convertToDouble(NValue);
    epsAchieved = abs((error1-error0)/error1);
    if (epsAchieved < eps)
     return;
    error0 = error1;
// 5.b)2. - normalize tensor:
    if (DMPS0 == 1)
     tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
    else
    {
     Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
     Tensor0.LQDecompose(Indices1, Indices2, TensorL);
    }
    MPS0.set(pos, Tensor0);
// 5.b)2. - set bTensors:
    bTensor.getShape(Shape5);
    Tensor0.getShape(Shape3);
    Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[4];
    Tensor0.reshape(Shape4);
    Tensor0.complexConjugate();
    Indices03[0] = 0; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
    bTensor.contract(Indices03, Tensor0, Indices13);
    bTensors[pos] = bTensor;
   }
// 5.b)3. position 0:
   MPS0.get(0, Tensor0);
   if (DMPS0 == 1)
    Tensor0.multiply(tensorFrobeniusNorm);
   else
   {
    Indices0[0] = 1; Indices1[0] = 0;
    Tensor0.contract(Indices0, TensorL, Indices1);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
   }
   MPS0.set(0, Tensor0);
   SweepDirection = "right";
  }
  numSweepsDone++;
 }
}

template<class T> void PEPS<T>::getMPOSL(const string& Direction, MPO<T>& MPO0) const
{
#ifdef DEBUG
 if ((this->N[0] == 0) ||
     ((Direction != "right") && (Direction != "down") && (Direction != "left") && (Direction != "up")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getMPOSL(const string& Direction, MPO<T>& MPO0) const: " <<
          "((this->N[0] == 0) || " <<
           "((Direction != right) && (Direction != down) && (Direction != left) && (Direction != up)))." << endl;
  exit(1);
 }
#endif
 string BCMPO0 = this->BC;
 unsigned int NMPO0;
 unsigned int dMPO0 = this->D, DMPO0 = this->D;
 Tensor<T> Tensor0;
 vector<unsigned int> Shape4(4), Shape5(5), Order(5);
 unsigned int position;
 if (Direction == "right")
 {
  NMPO0 = this->N[0];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order[0] = 3; Order[1] = 1; Order[2] = 0; Order[3] = 4; Order[4] = 2;
  for (position = 0; position < NMPO0; position++)
  {
   this->get(this->N[0]-1-position, 0, Tensor0);
   Tensor0.permute(Order);
   Tensor0.getShape(Shape5);
   Shape4[0] = Shape5[0]; Shape4[1] = Shape5[1]; Shape4[2] = Shape5[2]*Shape5[3]; Shape4[3] = Shape5[4];
   Tensor0.reshape(Shape4);
   MPO0.set(position, Tensor0);
  }
 }
 else if (Direction == "down")
 {
  NMPO0 = this->N[1];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order[0] = 0; Order[1] = 2; Order[2] = 1; Order[3] = 4; Order[4] = 3;
  for (position = 0; position < NMPO0; position++)
  {
   this->get(0, position, Tensor0);
   Tensor0.permute(Order);
   Tensor0.getShape(Shape5);
   Shape4[0] = Shape5[0]; Shape4[1] = Shape5[1]; Shape4[2] = Shape5[2]*Shape5[3]; Shape4[3] = Shape5[4];
   Tensor0.reshape(Shape4);
   MPO0.set(position, Tensor0);
  }
 }
 else if (Direction == "left")
 {
  NMPO0 = this->N[0];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order[0] = 1; Order[1] = 3; Order[2] = 2; Order[3] = 4; Order[4] = 0;
  for (position = 0; position < NMPO0; position++)
  {
   this->get(position, this->N[1]-1, Tensor0);
   Tensor0.permute(Order);
   Tensor0.getShape(Shape5);
   Shape4[0] = Shape5[0]; Shape4[1] = Shape5[1]; Shape4[2] = Shape5[2]*Shape5[3]; Shape4[3] = Shape5[4];
   Tensor0.reshape(Shape4);
   MPO0.set(position, Tensor0);
  }
 }
 else if (Direction == "up")
 {
  NMPO0 = this->N[1];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order[0] = 2; Order[1] = 0; Order[2] = 3; Order[3] = 4; Order[4] = 1;
  for (position = 0; position < NMPO0; position++)
  {
   this->get(this->N[0]-1, this->N[1]-1-position, Tensor0);
   Tensor0.permute(Order);
   Tensor0.getShape(Shape5);
   Shape4[0] = Shape5[0]; Shape4[1] = Shape5[1]; Shape4[2] = Shape5[2]*Shape5[3]; Shape4[3] = Shape5[4];
   Tensor0.reshape(Shape4);
   MPO0.set(position, Tensor0);
  }
 }
}

template<class T> void PEPS<T>::getMPO(const string& Direction, unsigned int position,
                                       MPO<T>& MPO0) const
{
#ifdef DEBUG
 if ((this->N[0] == 0) ||
     ((Direction != "right") && (Direction != "down") && (Direction != "left") &&
      (Direction != "up")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getMPO(const string& Direction, unsigned int position, MPO<T>& MPO0) const: " <<
          "((this->N[0] == 0) || " <<
           "((Direction != right) && (Direction != down) && (Direction != left) && " <<
            "(Direction != up)))." << endl;
  exit(1);
 }
 else if (((position > this->N[1]-1) && ((Direction == "right") || (Direction == "left"))) ||
          ((position > this->N[0]-1) && ((Direction == "down") || (Direction == "up"))))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getMPO(const string& Direction, unsigned int position, MPO<T>& MPO0) const: " <<
          "(((position > this->N[1]-1) && ((Direction == right) || (Direction == left))) || " <<
           "((position > this->N[0]-1) && ((Direction == down) || (Direction == up))))." << endl;
  exit(1);
 }
#endif
 string BCMPO0 = this->BC;
 unsigned int NMPO0;
 unsigned int dMPO0 = this->D*this->D, DMPO0 = this->D*this->D;
 vector<unsigned int> Indices0(1), Indices1(1);
 Indices0[0] = 4; Indices1[0] = 4;
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Shape0(4), Order0(8);
 T element;
 if (Direction == "right")
 {
  NMPO0 = this->N[0];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order0[0] = 3; Order0[1] = 7; Order0[2] = 1; Order0[3] = 5;
  Order0[4] = 0; Order0[5] = 4; Order0[6] = 2; Order0[7] = 6;
// the left end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = 1; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = this->Tensors[this->N[0]*(position+1)-1];
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(0, Tensor0);
// the intermediate part of MPO0:
  Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  for (int i = 1; i < NMPO0-1; i++)
  {
   Tensor0 = this->Tensors[this->N[0]*(position+1)-1-i];
   Tensor0.complexConjugate(Tensor1);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPO0.set(i, Tensor0);
  }
// the right end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = DMPO0; Shape0[1] = 1; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = this->Tensors[this->N[0]*position];
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(NMPO0-1, Tensor0);
 }
 else if (Direction == "down")
 {
  NMPO0 = this->N[1];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order0[0] = 0; Order0[1] = 4; Order0[2] = 2; Order0[3] = 6;
  Order0[4] = 1; Order0[5] = 5; Order0[6] = 3; Order0[7] = 7;
// the left end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = 1; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = this->Tensors[position];
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(0, Tensor0);
// the intermediate part of MPO0:
  Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  for (int i = 1; i < NMPO0-1; i++)
  {
   Tensor0 = this->Tensors[position+this->N[0]*i];
   Tensor0.complexConjugate(Tensor1);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPO0.set(i, Tensor0);
  }
// the right end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = DMPO0; Shape0[1] = 1; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = this->Tensors[position+this->N[0]*(this->N[1]-1)];
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(NMPO0-1, Tensor0);
 }
 else if (Direction == "left")
 {
  NMPO0 = this->N[0];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order0[0] = 1; Order0[1] = 5; Order0[2] = 3; Order0[3] = 7;
  Order0[4] = 2; Order0[5] = 6; Order0[6] = 0; Order0[7] = 4;
// the left end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = 1; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = this->Tensors[this->N[0]*position];
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(0, Tensor0);
// the intermediate part of MPO0:
  Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  for (int i = 1; i < NMPO0-1; i++)
  {
   Tensor0 = this->Tensors[i+this->N[0]*position];
   Tensor0.complexConjugate(Tensor1);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPO0.set(i, Tensor0);
  }
// the right end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = DMPO0; Shape0[1] = 1; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = this->Tensors[this->N[0]*(position+1)-1];
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(NMPO0-1, Tensor0);
 }
 else if (Direction == "up")
 {
  NMPO0 = this->N[1];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order0[0] = 2; Order0[1] = 6; Order0[2] = 0; Order0[3] = 4;
  Order0[4] = 3; Order0[5] = 7; Order0[6] = 1; Order0[7] = 5;
// the left end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = 1; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = this->Tensors[position+this->N[0]*(this->N[1]-1)];
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(0, Tensor0);
// the intermediate part of MPO0:
  Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  for (int i = 1; i < NMPO0-1; i++)
  {
   Tensor0 = this->Tensors[position+this->N[0]*(this->N[1]-1-i)];
   Tensor0.complexConjugate(Tensor1);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPO0.set(i, Tensor0);
  }
// the right end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = DMPO0; Shape0[1] = 1; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = this->Tensors[position];
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(NMPO0-1, Tensor0);
 }
}

template<class T> void PEPS<T>::getMPOSL(const string& Direction, unsigned int position, MPO<T>& MPO0) const
{
#ifdef DEBUG
 if ((this->N[0] == 0) ||
     ((Direction != "right") && (Direction != "down") && (Direction != "left") && (Direction != "up")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getMPOSL(const string& Direction, unsigned int position, MPO<T>& MPO0) const: " <<
          "((this->N[0] == 0) || " <<
           "((Direction != right) && (Direction != down) && (Direction != left) && (Direction != up)))." << endl;
  exit(1);
 }
 else if (((position > this->N[1]-1) && ((Direction == "right") || (Direction == "left"))) ||
          ((position > this->N[0]-1) && ((Direction == "down") || (Direction == "up"))))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getMPOSL(const string& Direction, unsigned int position, MPO<T>& MPO0) const: " <<
          "(((position > this->N[1]-1) && ((Direction == right) || (Direction == left))) || " <<
           "((position > this->N[0]-1) && ((Direction == down) || (Direction == up))))." << endl;
  exit(1);
 }
#endif
 string BCMPO0 = this->BC;
 unsigned int NMPO0;
 unsigned int dMPO0 = this->D, DMPO0 = this->D;
 Tensor<T> Tensor0;
 vector<unsigned int> Shape4(4), Shape5(5), Order(5);
 unsigned int position0;
 if (Direction == "right")
 {
  NMPO0 = this->N[0];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order[0] = 3; Order[1] = 1; Order[2] = 0; Order[3] = 2; Order[4] = 4;
  for (position0 = 0; position0 < NMPO0; position0++)
  {
   this->get(this->N[0]-1-position0, position, Tensor0);
   Tensor0.permute(Order);
   Tensor0.getShape(Shape5);
   Shape4[0] = Shape5[0]; Shape4[1] = Shape5[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[3]*Shape5[4];
   Tensor0.reshape(Shape4);
   MPO0.set(position0, Tensor0);
  }
 }
 else if (Direction == "down")
 {
  NMPO0 = this->N[1];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order[0] = 0; Order[1] = 2; Order[2] = 1; Order[3] = 3; Order[4] = 4;
  for (position0 = 0; position0 < NMPO0; position0++)
  {
   this->get(position, position0, Tensor0);
   Tensor0.permute(Order);
   Tensor0.getShape(Shape5);
   Shape4[0] = Shape5[0]; Shape4[1] = Shape5[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[3]*Shape5[4];
   Tensor0.reshape(Shape4);
   MPO0.set(position0, Tensor0);
  }
 }
 else if (Direction == "left")
 {
  NMPO0 = this->N[0];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order[0] = 1; Order[1] = 3; Order[2] = 2; Order[3] = 0; Order[4] = 4;
  for (position0 = 0; position0 < NMPO0; position0++)
  {
   this->get(position0, position, Tensor0);
   Tensor0.permute(Order);
   Tensor0.getShape(Shape5);
   Shape4[0] = Shape5[0]; Shape4[1] = Shape5[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[3]*Shape5[4];
   Tensor0.reshape(Shape4);
   MPO0.set(position0, Tensor0);
  }
 }
 else if (Direction == "up")
 {
  NMPO0 = this->N[1];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order[0] = 2; Order[1] = 0; Order[2] = 3; Order[3] = 1; Order[4] = 4;
  for (position0 = 0; position0 < NMPO0; position0++)
  {
   this->get(position, this->N[1]-1-position0, Tensor0);
   Tensor0.permute(Order);
   Tensor0.getShape(Shape5);
   Shape4[0] = Shape5[0]; Shape4[1] = Shape5[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[3]*Shape5[4];
   Tensor0.reshape(Shape4);
   MPO0.set(position0, Tensor0);
  }
 }
}

template<class T> void PEPS<T>::multiplyBulkMPOMPS(const string& Direction, unsigned int position, const MPS<T>& MPS0,
                                                   double eps, unsigned int maxNumSweeps,
                                                   double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const string& Direction, unsigned int position, const MPS<T>& MPS0, " <<
                             "double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "(this->N[0] == 0)." << endl;
  exit(1);
 }
 if ((Direction != "right") && (Direction != "down") && (Direction != "left") && (Direction != "up"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const string& Direction, unsigned int position, const MPS<T>& MPS0, " <<
                             "double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((Direction != right) && (Direction != down) && (Direction != left) && (Direction != up))." << endl;
  exit(1);
 }
 if ((((Direction == "right") || (Direction == "left")) && (position > this->N[1]-1)) ||
     (((Direction == "down") || (Direction == "up")) && (position > this->N[0]-1)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const string& Direction, unsigned int position, const MPS<T>& MPS0, " <<
                             "double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((((Direction == right) || (Direction == left)) && (position > this->N[1]-1)) || " <<
           "(((Direction == down) || (Direction == up)) && (position > this->N[0]-1)))." << endl;
  exit(1);
 }
 string BCMPS;
 MPS0.getBC(BCMPS);
 if ((BCMPS != "open") ||
     (((Direction == "right") || (Direction == "left")) && (MPS0.getN() != this->N[0])) ||
     (((Direction == "down") || (Direction == "up")) && (MPS0.getN() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const string& Direction, unsigned int position, const MPS<T>& MPS0, " <<
                             "double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((MPS0.BC != open) || (((Direction == right) || (Direction == left)) && (MPS0.N != this->N[0])) || " <<
           "(((Direction == down) || (Direction == up)) && (MPS0.N != this->N[1])))." << endl;
  exit(1);
 }
 MPS1.getBC(BCMPS);
 if ((BCMPS != "open") ||
     (((Direction == "right") || (Direction == "left")) && (MPS1.getN() != this->N[0])) ||
     (((Direction == "down") || (Direction == "up")) && (MPS1.getN() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const string& Direction, unsigned int position, const MPS<T>& MPS0, " <<
                             "double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((MPS1.BC != open) || (((Direction == right) || (Direction == left)) && (MPS1.N != this->N[0])) || " <<
           "(((Direction == down) || (Direction == up)) && (MPS1.N != this->N[1])))." << endl;
  exit(1);
 }
#endif
 unsigned int DMPS1, L;
 double error0, error1;
 T bValue, NValue, tensorFrobeniusNorm, tensorFrobeniusNormProduct;
 string SweepDirection;
 vector<unsigned int> dMPS, Indices0(1), Indices02(2), Indices03(3), Indices1(1), Indices12(2), Indices13(3), Indices2(2);
 vector<unsigned int> Order3(3), Order5(5), Order8(8), Shape3(3), Shape4(4), Shape5(5), Shape6(6), Shape8(8);
 Tensor<T> bTensor, Tensor0, Tensor1, TensorL, TensorR;
 vector< Tensor<T> > BulkTensors, bTensors;
// 0. initialize:
// 0. - set length L:
 if ((Direction == "right") || (Direction == "left"))
  L = this->N[0];
 else if ((Direction == "down") || (Direction == "up"))
  L = this->N[1];
// 0. - get bulk-tensors clockwise, permute them as Direction "down", and write them in BulkTensors:
 BulkTensors = vector< Tensor<T> >(L);
 if (Direction == "right")
 {
  Order5[0] = 3; Order5[1] = 0; Order5[2] = 1; Order5[3] = 2; Order5[4] = 4;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[this->N[0]-1-pos+this->N[0]*position];
   Tensor0.permute(Order5);
   BulkTensors[pos] = Tensor0;
  }
 }
 else if (Direction == "down")
 {
  for (int pos = 0; pos < L; pos++)
  {
   BulkTensors[pos] = this->Tensors[position+this->N[0]*pos];
  }
 }
 else if (Direction == "left")
 {
  Order5[0] = 1; Order5[1] = 2; Order5[2] = 3; Order5[3] = 0; Order5[4] = 4;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[pos+this->N[0]*position];
   Tensor0.permute(Order5);
   BulkTensors[pos] = Tensor0;
  }
 }
 else if (Direction == "up")
 {
  Order5[0] = 2; Order5[1] = 3; Order5[2] = 0; Order5[3] = 1; Order5[4] = 4;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[position+this->N[0]*(this->N[1]-1-pos)];
   Tensor0.permute(Order5);
   BulkTensors[pos] = Tensor0;
  }
 }
// 0. - check MPS0.d:
#ifdef DEBUG
 MPS0.getd(dMPS);
 for (int pos = 0; pos < L; pos++)
 {
  BulkTensors[pos].getShape(Shape5);
  if (dMPS[pos] != Shape5[1]*Shape5[1])
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void PEPS<T>::" <<
           "multiplyBulkMPOMPS(const string& Direction, unsigned int position, const MPS<T>& MPS0, " <<
                              "double eps, unsigned int maxNumSweeps, " <<
                              "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
           "(MPS0.d[" << pos << "] != this->D[" << pos << "]*this->D[" << pos << "])." << endl;
   exit(1);
  }
 }
#endif
// 0. - check MPS1.d:
 MPS1.getd(dMPS);
#ifdef DEBUG
 for (int pos = 0; pos < L; pos++)
 {
  BulkTensors[pos].getShape(Shape5);
  if (dMPS[pos] != Shape5[3]*Shape5[3])
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void PEPS<T>::" <<
           "multiplyBulkMPOMPS(const string& Direction, unsigned int position, const MPS<T>& MPS0, " <<
                              "double eps, unsigned int maxNumSweeps, " <<
                              "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
           "(MPS1.d[" << pos << "] != this->D[" << pos << "]*this->D[" << pos << "])." << endl;
   exit(1);
  }
 }
#endif
// 1. MPS1.D >= this->D^2*MPS0.D:
 DMPS1 = MPS1.getD();
 if (DMPS1 >= this->D*this->D*MPS0.getD())
 {
  MPS1 = MPS<T>("open", L, dMPS, this->D*this->D*MPS0.getD());
  Indices0[0] = 2; Indices1[0] = 1;
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 4;
  Order8[0] = 0; Order8[1] = 2; Order8[2] = 5; Order8[3] = 1; Order8[4] = 3; Order8[5] = 6; Order8[6] = 4; Order8[7] = 7;
  for (int pos = 0; pos < L; pos++)
  {
   MPS0.get(pos, Tensor0);
   Tensor0.getShape(Shape3);
   BulkTensors[pos].getShape(Shape5);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[1]; Shape4[3] = Shape5[1];
   Tensor0.reshape(Shape4);
   Tensor1 = BulkTensors[pos];
   Tensor0.contract(Indices0, Tensor1, Indices1);
   BulkTensors[pos].complexConjugate(Tensor1);
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.permute(Order8);
   Tensor0.getShape(Shape8);
   Shape3[0] = Shape8[0]*Shape8[1]*Shape8[2]; Shape3[1] = Shape8[3]*Shape8[4]*Shape8[5]; Shape3[2] = Shape8[6]*Shape8[7];
   Tensor0.reshape(Shape3);
   MPS1.set(pos, Tensor0);
  }
  MPS1.bringIntoNormalShape();
  MPS1.setD(DMPS1);
  epsAchieved = 0.0;
  numSweepsDone = 0;
  return;
 }
// 2. bring MPS1 into normal form from right to left:
// 2.a) MPS1.D == 1:
 if (DMPS1 == 1)
 {
  tensorFrobeniusNormProduct = 1.0;
  for (int pos = L-1; pos > 0; pos--)
  {
   MPS1.get(pos, Tensor0);
   tensorFrobeniusNormProduct *= Tensor0.frobeniusNormalize();
   MPS1.set(pos, Tensor0);
  }
  MPS1.get(0, Tensor0);
  Tensor0.multiply(tensorFrobeniusNormProduct);
  MPS1.set(0, Tensor0);
 }
// 2.b) MPS1.D != 1:
 else
 {
  Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
  Indices0[0] = 1;
  Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
  MPS1.get(L-1, Tensor0);
  for (int pos = L-1; pos > 0; pos--)
  {
   Tensor0.LQDecompose(Indices1, Indices2, TensorL);
   MPS1.set(pos, Tensor0);
   MPS1.get(pos-1, Tensor0);
   Tensor0.contract(Indices0, TensorL, Indices1);
   Tensor0.permute(Order3);
  }
  MPS1.set(0, Tensor0);
 }
// 3. construct right half of b for each site:
 bTensors = vector< Tensor<T> >(L);
 MPS0.get(L-1, Tensor0);
 Tensor0.getShape(Shape3);
 BulkTensors[L-1].getShape(Shape5);
 Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[1]; Shape4[3] = Shape5[1];
 Tensor0.reshape(Shape4);
 Tensor1 = BulkTensors[L-1];
 Indices0[0] = 2; Indices1[0] = 1;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 BulkTensors[L-1].complexConjugate(Tensor1);
 Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 4;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 MPS1.get(L-1, Tensor1);
 Tensor1.getShape(Shape3);
 Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[3]; Shape4[3] = Shape5[3];
 Tensor1.reshape(Shape4);
 Tensor1.complexConjugate();
 Indices02[0] = 4; Indices02[1] = 7; Indices12[0] = 2; Indices12[1] = 3;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 Tensor0.getShape(Shape8);
 Shape4[0] = Shape8[0]; Shape4[1] = Shape8[2]; Shape4[2] = Shape8[4]; Shape4[3] = Shape8[6];
 Tensor0.reshape(Shape4);
 bTensors[L-1] = Tensor0;
 Indices0[0] = 0; Indices1[0] = 1;
 Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 2; Indices12[1] = 1;
 Indices03[0] = 0; Indices03[1] = 3;
 for (int pos = L-2; pos > 0; pos--)
 {
  MPS0.get(pos, Tensor1);
  Tensor1.getShape(Shape3);
  BulkTensors[pos].getShape(Shape5);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[1]; Shape4[3] = Shape5[1];
  Tensor1.reshape(Shape4);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor1 = BulkTensors[pos];
  Tensor0.contract(Indices02, Tensor1, Indices12);
  BulkTensors[pos].complexConjugate(Tensor1);
  Indices03[2] = 6; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  MPS1.get(pos, Tensor1);
  Tensor1.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[3]; Shape4[3] = Shape5[3];
  Tensor1.reshape(Shape4);
  Tensor1.complexConjugate();
  Indices03[2] = 5; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  bTensors[pos] = Tensor0;
 }
// 4. compute initial error:
 MPS0.get(0, Tensor1);
 Tensor1.getShape(Shape3);
 BulkTensors[0].getShape(Shape5);
 Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[1]; Shape4[3] = Shape5[1];
 Tensor1.reshape(Shape4);
 Indices0[0] = 0; Indices1[0] = 1;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = BulkTensors[0];
 Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 2; Indices12[1] = 1;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 BulkTensors[0].complexConjugate(Tensor1);
 Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
 Tensor0.contract(Indices03, Tensor1, Indices13);
 MPS1.get(0, Tensor1);
 Tensor1.getShape(Shape3);
 Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[3]; Shape4[3] = Shape5[3];
 Tensor1.reshape(Shape4);
 Tensor1.complexConjugate();
 Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 5; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
 Tensor0.contract(Indices03, Tensor1, Indices13);
 bValue = Tensor0.get(0);
 MPS1.get(0, Tensor0);
 NValue = Tensor0.scalarProduct();
 error0 = MathAuxiliary::convertToDouble(NValue)-2.0*MathAuxiliary::convertToDouble(bValue);
// 5. sweep:
 numSweepsDone = 0;
 SweepDirection = "right";
 while (numSweepsDone < maxNumSweeps)
 {
// 5.a) sweep to the right:
  if (SweepDirection == "right")
  {
// 5.a)1. position 0:
// 5.a)1. - compute bTensor:
   Tensor0 = bTensors[1];
   MPS0.get(0, Tensor1);
   Tensor1.getShape(Shape3);
   BulkTensors[0].getShape(Shape5);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[1]; Shape4[3] = Shape5[1];
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BulkTensors[0];
   Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 2; Indices12[1] = 1;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   BulkTensors[0].complexConjugate(Tensor1);
   Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
// 5.a)1. - update tensor:
   Tensor0.getShape(Shape6);
   Shape3[0] = 1; Shape3[1] = Shape6[0]; Shape3[2] = Shape6[3]*Shape6[5];
   Tensor0.reshape(Shape3);
   MPS1.set(0, Tensor0);
// 5.a)1. - compute error:
   NValue = Tensor0.scalarProduct();
   error1 = -MathAuxiliary::convertToDouble(NValue);
   epsAchieved = abs((error1-error0)/error1);
   if (epsAchieved < eps)
    return;
   error0 = error1;
// 5.a)1. - normalize tensor:
   if (DMPS1 == 1)
    tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
   else
   {
    Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
    Tensor0.QRDecompose(Indices2, Indices1, TensorR);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
   }
   MPS1.set(0, Tensor0);
// 5.a)1. - set bTensors:
   MPS0.get(0, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[1]; Shape4[3] = Shape5[1];
   Tensor0.reshape(Shape4);
   Tensor1 = BulkTensors[0];
   Indices0[0] = 2; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   BulkTensors[0].complexConjugate(Tensor1);
   Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 4;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   MPS1.get(0, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[3]; Shape4[3] = Shape5[3];
   Tensor1.reshape(Shape4);
   Tensor1.complexConjugate();
   Indices02[0] = 4; Indices02[1] = 7; Indices12[0] = 2; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape8);
   Shape4[0] = Shape8[1]; Shape4[1] = Shape8[3]; Shape4[2] = Shape8[5]; Shape4[3] = Shape8[7];
   Tensor0.reshape(Shape4);
   bTensors[0] = Tensor0;
// 5.a)2. 0 < pos < L-1:
   for (int pos = 1; pos < L-1; pos++)
   {
// 5.a)2. - compute bTensor:
    Tensor0 = bTensors[pos-1];
    MPS0.get(pos, Tensor1);
    Tensor1.getShape(Shape3);
    BulkTensors[pos].getShape(Shape5);
    Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[1]; Shape4[3] = Shape5[1];
    Tensor1.reshape(Shape4);
    Indices0[0] = 0; Indices1[0] = 0;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = BulkTensors[pos];
    Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 1;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    BulkTensors[pos].complexConjugate(Tensor1);
    Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 4;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    bTensor = Tensor0;
// 5.a)2. - update tensor:
    Tensor1 = bTensors[pos+1];
    Indices03[0] = 1; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 2;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    Tensor0.getShape(Shape4);
    Shape3[0] = Shape4[0]; Shape3[1] = Shape4[1]*Shape4[2]; Shape3[2] = Shape4[3];
    Tensor0.reshape(Shape3);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
    MPS1.set(pos, Tensor0);
// 5.a)2. - compute error:
    NValue = Tensor0.scalarProduct();
    error1 = -MathAuxiliary::convertToDouble(NValue);
    epsAchieved = abs((error1-error0)/error1);
    if (epsAchieved < eps)
     return;
    error0 = error1;
// 5.a)2. - normalize tensor:
    if (DMPS1 == 1)
     tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
    else
    {
     Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
     Tensor0.QRDecompose(Indices2, Indices1, TensorR);
     Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
     Tensor0.permute(Order3);
    }
    MPS1.set(pos, Tensor0);
// 5.a)2. - set bTensors:
    Tensor0.getShape(Shape3);
    Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[3]; Shape4[3] = Shape5[3];
    Tensor0.reshape(Shape4);
    Tensor0.complexConjugate();
    Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 5; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
    bTensor.contract(Indices03, Tensor0, Indices13);
    bTensors[pos] = bTensor;
   }
// 5.a)3. position L-1:
   MPS1.get(L-1, Tensor0);
   if (DMPS1 == 1)
    Tensor0.multiply(tensorFrobeniusNorm);
   else
   {
    Indices0[0] = 1; Indices1[0] = 0;
    TensorR.contract(Indices0, Tensor0, Indices1);
    Tensor0 = TensorR;
   }
   MPS1.set(L-1, Tensor0);
   SweepDirection = "left";
  }
// 5.b) sweep to the left:
  else if (SweepDirection == "left")
  {
// 5.b)1. position L-1:
// 5.b)1. - compute bTensor:
   Tensor0 = bTensors[L-2];
   MPS0.get(L-1, Tensor1);
   Tensor1.getShape(Shape3);
   BulkTensors[L-1].getShape(Shape5);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape5[1]; Shape4[3] = Shape5[1];
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BulkTensors[L-1];
   Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 1;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   BulkTensors[L-1].complexConjugate(Tensor1);
   Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
// 5.b)1. - update tensor:
   Tensor0.getShape(Shape6);
   Shape3[0] = Shape6[0]; Shape3[1] = 1; Shape3[2] = Shape6[3]*Shape6[5];
   Tensor0.reshape(Shape3);
   MPS1.set(L-1, Tensor0);
// 5.b)1. - compute error:
   NValue = Tensor0.scalarProduct();
   error1 = -MathAuxiliary::convertToDouble(NValue);
   epsAchieved = abs((error1-error0)/error1);
   if (epsAchieved < eps)
    return;
   error0 = error1;
// 5.b)1. - normalize tensor:
   if (DMPS1 == 1)
    tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
   else
   {
    Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
    Tensor0.LQDecompose(Indices1, Indices2, TensorL);
   }
   MPS1.set(L-1, Tensor0);
// 5.b)1. - set bTensors:
   MPS0.get(L-1, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape5[1]; Shape4[3] = Shape5[1];
   Tensor0.reshape(Shape4);
   Tensor1 = BulkTensors[L-1];
   Indices0[0] = 2; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   BulkTensors[L-1].complexConjugate(Tensor1);
   Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 4;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   MPS1.get(L-1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape5[3]; Shape4[3] = Shape5[3];
   Tensor1.reshape(Shape4);
   Tensor1.complexConjugate();
   Indices02[0] = 4; Indices02[1] = 7; Indices12[0] = 2; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape8);
   Shape4[0] = Shape8[0]; Shape4[1] = Shape8[2]; Shape4[2] = Shape8[4]; Shape4[3] = Shape8[6];
   Tensor0.reshape(Shape4);
   bTensors[L-1] = Tensor0;
// 5.b)2. 0 < pos < L-1:
   for (int pos = L-2; pos > 0; pos--)
   {
// 5.b)2. - compute bTensor:
    Tensor0 = bTensors[pos+1];
    MPS0.get(pos, Tensor1);
    Tensor1.getShape(Shape3);
    BulkTensors[pos].getShape(Shape5);
    Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[1]; Shape4[3] = Shape5[1];
    Tensor1.reshape(Shape4);
    Indices0[0] = 0; Indices1[0] = 1;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = BulkTensors[pos];
    Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 2; Indices12[1] = 1;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    BulkTensors[pos].complexConjugate(Tensor1);
    Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    bTensor = Tensor0;
// 5.b)2. - update tensor:
    Tensor1 = bTensors[pos-1];
    Indices03[0] = 1; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 2;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    Tensor0.getShape(Shape4);
    Shape3[0] = Shape4[0]; Shape3[1] = Shape4[1]*Shape4[2]; Shape3[2] = Shape4[3];
    Tensor0.reshape(Shape3);
    Order3[0] = 2; Order3[1] = 0; Order3[2] = 1;
    Tensor0.permute(Order3);
    MPS1.set(pos, Tensor0);
// 5.b)2. - compute error:
    NValue = Tensor0.scalarProduct();
    error1 = -MathAuxiliary::convertToDouble(NValue);
    epsAchieved = abs((error1-error0)/error1);
    if (epsAchieved < eps)
     return;
    error0 = error1;
// 5.b)2. - normalize tensor:
    if (DMPS1 == 1)
     tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
    else
    {
     Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
     Tensor0.LQDecompose(Indices1, Indices2, TensorL);
    }
    MPS1.set(pos, Tensor0);
// 5.b)2. - set bTensors:
    Tensor0.getShape(Shape3);
    Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[3]; Shape4[3] = Shape5[3];
    Tensor0.reshape(Shape4);
    Tensor0.complexConjugate();
    Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 5; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
    bTensor.contract(Indices03, Tensor0, Indices13);
    bTensors[pos] = bTensor;
   }
// 5.b)3. position 0:
   MPS1.get(0, Tensor0);
   if (DMPS1 == 1)
    Tensor0.multiply(tensorFrobeniusNorm);
   else
   {
    Indices0[0] = 1; Indices1[0] = 0;
    Tensor0.contract(Indices0, TensorL, Indices1);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
   }
   MPS1.set(0, Tensor0);
   SweepDirection = "right";
  }
  numSweepsDone++;
 }
}

template<class T> T PEPS<T>::contractReverseBulkMPOMPS(const string& Direction,
                                                       const MPS<T>& MPS0, unsigned int position, const MPS<T>& MPS1) const
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "contractReverseBulkMPOMPS(const string& Direction, const MPS<T>& MPS0, " <<
                                    "unsigned int position, const MPS<T>& MPS1) const: " <<
          "(this->N[0] == 0)." << endl;
  exit(1);
 }
 if ((Direction != "horizontal") && (Direction != "vertical"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "contractReverseBulkMPOMPS(const string& Direction, const MPS<T>& MPS0, " <<
                                    "unsigned int position, const MPS<T>& MPS1) const: " <<
          "((Direction != horizontal) && (Direction != vertical))." << endl;
  exit(1);
 }
 string BCMPS;
 MPS0.getBC(BCMPS);
 if ((BCMPS != "open") ||
     ((Direction == "horizontal") && (MPS0.getN() != this->N[0])) ||
     ((Direction == "vertical") && (MPS0.getN() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "contractReverseBulkMPOMPS(const string& Direction, const MPS<T>& MPS0, " <<
                                    "unsigned int position, const MPS<T>& MPS1) const: " <<
          "((MPS0.BC != open) || ((Direction == horizontal) && (MPS0.N != this->N[0])) || " <<
           "((Direction == vertical) && (MPS0.N != this->N[1])))." << endl;
  exit(1);
 }
 if (((Direction == "horizontal") && (position > this->N[1]-1)) ||
     ((Direction == "vertical") && (position > this->N[0]-1)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "contractReverseBulkMPOMPS(const string& Direction, const MPS<T>& MPS0, " <<
                                    "unsigned int position, const MPS<T>& MPS1) const: " <<
          "(((Direction == horizontal) && (position > this->N[1]-1)) || " <<
           "((Direction == vertical) && (position > this->N[0]-1)))." << endl;
  exit(1);
 }
 MPS1.getBC(BCMPS);
 if ((BCMPS != "open") ||
     ((Direction == "horizontal") && (MPS1.getN() != this->N[0])) ||
     ((Direction == "vertical") && (MPS1.getN() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "contractReverseBulkMPOMPS(const string& Direction, const MPS<T>& MPS0, " <<
                                    "unsigned int position, const MPS<T>& MPS1) const: " <<
          "((MPS1.BC != open) || ((Direction == horizontal) && (MPS1.N != this->N[0])) || " <<
           "((Direction == vertical) && (MPS1.N != this->N[1])))." << endl;
  exit(1);
 }
#endif
 unsigned int L;
 vector<unsigned int> dMPS, Indices0(1), Indices02(2), Indices03(3), Indices1(1), Indices12(2), Indices13(3);
 vector<unsigned int> Order5(5), Shape3(3), Shape4(4), Shape5(5);
 Tensor<T> Tensor0, Tensor1;
 vector< Tensor<T> > BulkTensors;
// 0. initialize:
// 0. - set length L:
 if (Direction == "horizontal")
  L = this->N[0];
 else if (Direction == "vertical")
  L = this->N[1];
// 0. - get bulk-tensors, permute them as Direction "vertical", and write them in BulkTensors:
 BulkTensors = vector< Tensor<T> >(L);
 if (Direction == "horizontal")
 {
  Order5[0] = 1; Order5[1] = 2; Order5[2] = 3; Order5[3] = 0; Order5[4] = 4;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[pos+this->N[0]*position];
   Tensor0.permute(Order5);
   BulkTensors[pos] = Tensor0;
  }
 }
 else if (Direction == "vertical")
 {
  for (int pos = 0; pos < L; pos++)
  {
   BulkTensors[pos] = this->Tensors[position+this->N[0]*pos];
  }
 }
// 0. - check MPS0.d and MPS1.d:
#ifdef DEBUG
 MPS0.getd(dMPS);
 for (int pos = 0; pos < L; pos++)
 {
  BulkTensors[L-1-pos].getShape(Shape5);
  if (dMPS[pos] != Shape5[3]*Shape5[3])
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> T PEPS<T>::" <<
           "contractReverseBulkMPOMPS(const string& Direction, const MPS<T>& MPS0, " <<
                                     "unsigned int position, const MPS<T>& MPS1) const: " <<
           "(MPS0.d[" << pos << "] != this->D[" << L-1-pos << "]*this->D[" << L-1-pos << "])." << endl;
   exit(1);
  }
 }
 MPS1.getd(dMPS);
 for (int pos = 0; pos < L; pos++)
 {
  BulkTensors[pos].getShape(Shape5);
  if (dMPS[pos] != Shape5[1]*Shape5[1])
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> T PEPS<T>::" <<
           "contractReverseBulkMPOMPS(const string& Direction, const MPS<T>& MPS0, " <<
                                     "unsigned int position, const MPS<T>& MPS1) const: " <<
           "(MPS1.d[" << pos << "] != this->D[" << pos << "]*this->D[" << pos << "])." << endl;
   exit(1);
  }
 }
#endif
// 1. position 0:
 MPS1.get(0, Tensor0);
 Tensor0.getShape(Shape3);
 BulkTensors[0].getShape(Shape5);
 Shape3[0] = Shape3[1]; Shape3[1] = Shape5[1]; Shape3[2] = Shape5[1];
 Tensor0.reshape(Shape3);
 Shape4[0] = Shape5[1]; Shape4[1] = Shape5[2]; Shape4[2] = Shape5[3]; Shape4[3] = Shape5[4];
 BulkTensors[0].reshape(Shape4);
 Tensor1 = BulkTensors[0];
 Indices0[0] = 1; Indices1[0] = 0;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 BulkTensors[0].complexConjugate();
 Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 3;
 Tensor0.contract(Indices02, BulkTensors[0], Indices12);
 MPS0.get(L-1, Tensor1);
 Tensor1.getShape(Shape3);
 Shape3[1] = Shape5[3]; Shape3[2] = Shape5[3];
 Tensor1.reshape(Shape3);
 Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
 Tensor0.contract(Indices02, Tensor1, Indices12);
// 2. 0 < position <= L-1:
 Indices0[0] = 0; Indices1[0] = 0;
 Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 1;
 Indices03[0] = 0; Indices03[1] = 3;
 for (int pos = 1; pos < L; pos++)
 {
  MPS1.get(pos, Tensor1);
  Tensor1.getShape(Shape3);
  BulkTensors[pos].getShape(Shape5);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[1]; Shape4[3] = Shape5[1];
  Tensor1.reshape(Shape4);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor1 = BulkTensors[pos];
  Tensor0.contract(Indices02, Tensor1, Indices12);
  BulkTensors[pos].complexConjugate();
  Indices03[2] = 6; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 4;
  Tensor0.contract(Indices03, BulkTensors[pos], Indices13);
  MPS0.get(L-1-pos, Tensor1);
  Tensor1.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[3]; Shape4[3] = Shape5[3];
  Tensor1.reshape(Shape4);
  Indices03[2] = 5; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
  Tensor0.contract(Indices03, Tensor1, Indices13);
 }
// 3. return contraction-value:
 return Tensor0.get(0);
}

template<class T> void PEPS<T>::getMPS(const PEPS<T>& PEPS0, const string& Direction,
                                       MPS<T>& MPS0) const
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (this->BC != "open") || (this->BC != PEPS0.BC) ||
     (this->N != PEPS0.N) || (this->d != PEPS0.d) ||
     ((Direction != "right") && (Direction != "down") && (Direction != "left") &&
      (Direction != "up")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getMPS(const PEPS<T>& PEPS0, const string& Direction, " <<
                 "MPS<T>& MPS0) const: " <<
          "((this->N[0] == 0) || (this->BC != open) || (this->BC != PEPS0.BC) || " <<
           "(this->N != PEPS0.N) || (this->d != PEPS0.d) || " <<
           "((Direction != right) && (Direction != down) && (Direction != left) && " <<
            "(Direction != up)))." << endl;
  exit(1);
 }
#endif
 string BCMPS0 = "open";
 unsigned int NMPS0;
 unsigned int dMPS0 = PEPS0.D*this->D, DMPS0 = PEPS0.D*this->D;
 vector<unsigned int> Indices0(1), Indices1(1);
 Indices0[0] = 4; Indices1[0] = 4;
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Shape0(3), Order0(8);
 T element;
 if (Direction == "right")
 {
  NMPS0 = this->N[0];
  MPS0 = MPS<T>(BCMPS0, NMPS0, dMPS0, DMPS0);
  Order0[0] = 3; Order0[1] = 7; Order0[2] = 0; Order0[3] = 4;
  Order0[4] = 1; Order0[5] = 5; Order0[6] = 2; Order0[7] = 6;
// the left end of MPS0:
  Shape0[0] = 1; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  Tensor0 = PEPS0.Tensors[this->N[0]-1];
  Tensor1 = this->Tensors[this->N[0]-1];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(0, Tensor0);
// the intermediate part of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  for (int i = 1; i < NMPS0-1; i++)
  {
   Tensor0 = PEPS0.Tensors[this->N[0]-1-i];
   Tensor1 = this->Tensors[this->N[0]-1-i];
   Tensor1.complexConjugate();
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPS0.set(i, Tensor0);
  }
// the right end of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = 1; Shape0[2] = dMPS0;
  Tensor0 = PEPS0.Tensors[0];
  Tensor1 = this->Tensors[0];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(NMPS0-1, Tensor0);
 }
 else if (Direction == "down")
 {
  NMPS0 = this->N[1];
  MPS0 = MPS<T>(BCMPS0, NMPS0, dMPS0, DMPS0);
  Order0[0] = 0; Order0[1] = 4; Order0[2] = 1; Order0[3] = 5;
  Order0[4] = 2; Order0[5] = 6; Order0[6] = 3; Order0[7] = 7;
// the left end of MPS0:
  Shape0[0] = 1; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  Tensor0 = PEPS0.Tensors[0];
  Tensor1 = this->Tensors[0];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(0, Tensor0);
// the intermediate part of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  for (int i = 1; i < NMPS0-1; i++)
  {
   Tensor0 = PEPS0.Tensors[this->N[0]*i];
   Tensor1 = this->Tensors[this->N[0]*i];
   Tensor1.complexConjugate();
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPS0.set(i, Tensor0);
  }
// the right end of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = 1; Shape0[2] = dMPS0;
  Tensor0 = PEPS0.Tensors[this->N[0]*(this->N[1]-1)];
  Tensor1 = this->Tensors[this->N[0]*(this->N[1]-1)];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(NMPS0-1, Tensor0);
 }
 else if (Direction == "left")
 {
  NMPS0 = this->N[0];
  MPS0 = MPS<T>(BCMPS0, NMPS0, dMPS0, DMPS0);
  Order0[0] = 1; Order0[1] = 5; Order0[2] = 2; Order0[3] = 6;
  Order0[4] = 3; Order0[5] = 7; Order0[6] = 0; Order0[7] = 4;
// the left end of MPS0:
  Shape0[0] = 1; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  Tensor0 = PEPS0.Tensors[this->N[0]*(this->N[1]-1)];
  Tensor1 = this->Tensors[this->N[0]*(this->N[1]-1)];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(0, Tensor0);
// the intermediate part of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  for (int i = 1; i < NMPS0-1; i++)
  {
   Tensor0 = PEPS0.Tensors[i+this->N[0]*(this->N[1]-1)];
   Tensor1 = this->Tensors[i+this->N[0]*(this->N[1]-1)];
   Tensor1.complexConjugate();
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPS0.set(i, Tensor0);
  }
// the right end of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = 1; Shape0[2] = dMPS0;
  Tensor0 = PEPS0.Tensors[this->N[0]*this->N[1]-1];
  Tensor1 = this->Tensors[this->N[0]*this->N[1]-1];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(NMPS0-1, Tensor0);
 }
 else if (Direction == "up")
 {
  NMPS0 = this->N[1];
  MPS0 = MPS<T>(BCMPS0, NMPS0, dMPS0, DMPS0);
  Order0[0] = 2; Order0[1] = 6; Order0[2] = 3; Order0[3] = 7;
  Order0[4] = 0; Order0[5] = 4; Order0[6] = 1; Order0[7] = 5;
// the left end of MPS0:
  Shape0[0] = 1; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  Tensor0 = PEPS0.Tensors[this->N[0]*this->N[1]-1];
  Tensor1 = this->Tensors[this->N[0]*this->N[1]-1];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(0, Tensor0);
// the intermediate part of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  for (int i = 1; i < NMPS0-1; i++)
  {
   Tensor0 = PEPS0.Tensors[this->N[0]*(this->N[1]-i)-1];
   Tensor1 = this->Tensors[this->N[0]*(this->N[1]-i)-1];
   Tensor1.complexConjugate();
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPS0.set(i, Tensor0);
  }
// the right end of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = 1; Shape0[2] = dMPS0;
  Tensor0 = PEPS0.Tensors[this->N[0]-1];
  Tensor1 = this->Tensors[this->N[0]-1];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(NMPS0-1, Tensor0);
 }
}

template<class T> void PEPS<T>::getBoundaryMPS(const PEPS<T>& PEPS0, const string& Boundary,
                                               double eps, unsigned int maxNumSweeps,
                                               double& epsAchieved, unsigned int& numSweepsDone,
                                               MPS<T>& MPS0) const
{
#ifdef DEBUG
 if ((this->BC != "open") || (this->N[0] == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getBoundaryMPS(const PEPS<T>& PEPS0, const string& Boundary, double eps, unsigned int maxNumSweeps, " <<
                         "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "((this->BC != open) || (this->N[0] == 0))." << endl;
  exit(1);
 }
 if ((this->BC != PEPS0.BC) || (this->N != PEPS0.N) || (this->d != PEPS0.d))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getBoundaryMPS(const PEPS<T>& PEPS0, const string& Boundary, double eps, unsigned int maxNumSweeps, " <<
                         "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "((this->BC != PEPS0.BC) || (this->N != PEPS0.N) || (this->d != PEPS0.d))." << endl;
  exit(1);
 }
 if ((Boundary != "left") && (Boundary != "top") && (Boundary != "right") && (Boundary != "bottom"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getBoundaryMPS(const PEPS<T>& PEPS0, const string& Boundary, double eps, unsigned int maxNumSweeps, " <<
                         "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "((Boundary != left) && (Boundary != top) && (Boundary != right) && (Boundary != bottom))." << endl;
  exit(1);
 }
 string BCMPS0; MPS0.getBC(BCMPS0);
 if ((BCMPS0 != "open") ||
     (((Boundary == "left") || (Boundary == "right")) && (MPS0.getN() != this->N[0])) ||
     (((Boundary == "top") || (Boundary == "bottom")) && (MPS0.getN() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getBoundaryMPS(const PEPS<T>& PEPS0, const string& Boundary, double eps, unsigned int maxNumSweeps, " <<
                         "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "((MPS0.BC != open) || (((Boundary == left) || (Boundary == right)) && (MPS0.N != this->N[0])) || " <<
           "(((Boundary == top) || (Boundary == bottom)) && (MPS0.N != this->N[1])))." << endl;
  exit(1);
 }
#endif
 unsigned int DMPS0, L;
 double error0, error1;
 T bValue, NValue, tensorFrobeniusNorm, tensorFrobeniusNormProduct;
 string SweepDirection;
 vector<unsigned int> dMPS0, Indices0(1), Indices02(2), Indices03(3), Indices1(1), Indices12(2), Indices13(3), Indices2(2);
 vector<unsigned int> Order3(3), Order5(5), Order6(6), Shape3(3), Shape4(4), Shape04(4), Shape5(5), Shape6(6);
 Tensor<T> bTensor, Tensor0, Tensor1, TensorL, TensorR;
 vector< Tensor<T> > BoundaryTensors, BoundaryTensors0, bTensors;
// 0. initialize:
// 0. - set boundary-length L:
 if ((Boundary == "left") || (Boundary == "right"))
  L = this->N[0];
 else if ((Boundary == "top") || (Boundary == "bottom"))
  L = this->N[1];
// 0. - get boundary-tensors clockwise, permute them as "top" boundary, and write the ones from this PEPS in BoundaryTensors and the ones from PEPS0 in BoundaryTensors0:
 BoundaryTensors = vector< Tensor<T> >(L);
 BoundaryTensors0 = vector< Tensor<T> >(L);
 if (Boundary == "left")
 {
  Order5[0] = 3; Order5[1] = 0; Order5[2] = 1; Order5[3] = 2; Order5[4] = 4;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[this->N[0]-1-pos];
   Tensor0.permute(Order5);
   BoundaryTensors[pos] = Tensor0;
   Tensor0 = PEPS0.Tensors[this->N[0]-1-pos];
   Tensor0.permute(Order5);
   BoundaryTensors0[pos] = Tensor0;
  }
 }
 else if (Boundary == "top")
 {
  for (int pos = 0; pos < L; pos++)
  {
   BoundaryTensors[pos] = this->Tensors[this->N[0]*pos];
   BoundaryTensors0[pos] = PEPS0.Tensors[this->N[0]*pos];
  }
 }
 else if (Boundary == "right")
 {
  Order5[0] = 1; Order5[1] = 2; Order5[2] = 3; Order5[3] = 0; Order5[4] = 4;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[pos+this->N[0]*(this->N[1]-1)];
   Tensor0.permute(Order5);
   BoundaryTensors[pos] = Tensor0;
   Tensor0 = PEPS0.Tensors[pos+this->N[0]*(this->N[1]-1)];
   Tensor0.permute(Order5);
   BoundaryTensors0[pos] = Tensor0;
  }
 }
 else if (Boundary == "bottom")
 {
  Order5[0] = 2; Order5[1] = 3; Order5[2] = 0; Order5[3] = 1; Order5[4] = 4;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[this->N[0]-1+this->N[0]*(this->N[1]-1-pos)];
   Tensor0.permute(Order5);
   BoundaryTensors[pos] = Tensor0;
   Tensor0 = PEPS0.Tensors[this->N[0]-1+this->N[0]*(this->N[1]-1-pos)];
   Tensor0.permute(Order5);
   BoundaryTensors0[pos] = Tensor0;
  }
 }
// 0. - reshape boundary-tensors eliminating index 1:
 for (int pos = 0; pos < L; pos++)
 {
  BoundaryTensors[pos].getShape(Shape5);
  Shape4[0] = Shape5[0]; Shape4[1] = Shape5[2]; Shape4[2] = Shape5[3]; Shape4[3] = Shape5[4];
  BoundaryTensors[pos].reshape(Shape4);
  BoundaryTensors0[pos].getShape(Shape5);
  Shape4[0] = Shape5[0]; Shape4[1] = Shape5[2]; Shape4[2] = Shape5[3]; Shape4[3] = Shape5[4];
  BoundaryTensors0[pos].reshape(Shape4);
 }
// 0. - check MPS0.d:
 MPS0.getd(dMPS0);
#ifdef DEBUG
 for (int pos = 0; pos < L; pos++)
 {
  BoundaryTensors0[pos].getShape(Shape04);
  BoundaryTensors[pos].getShape(Shape4);
  if (dMPS0[pos] != Shape04[2]*Shape4[2])
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void PEPS<T>::" <<
           "getBoundaryMPS(const PEPS<T>& PEPS0, const string& Boundary, double eps, unsigned int maxNumSweeps, " <<
                          "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
           "(MPS0.d[" << pos << "] != PEPS0.D[" << pos << "]*this->D[" << pos << "])." << endl;
   exit(1);
  }
 }
#endif
// 1. MPS0.D >= PEPS0.D*this->D:
 DMPS0 = MPS0.getD();
 if (DMPS0 >= PEPS0.D*this->D)
 {
  MPS0 = MPS<T>("open", L, dMPS0, PEPS0.D*this->D);
  Indices0[0] = 3; Indices1[0] = 3;
  Order6[0] = 0; Order6[1] = 3; Order6[2] = 1; Order6[3] = 4; Order6[4] = 2; Order6[5] = 5;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = BoundaryTensors0[pos];
   BoundaryTensors[pos].complexConjugate(Tensor1);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order6);
   Tensor0.getShape(Shape6);
   Shape3[0] = Shape6[0]*Shape6[1]; Shape3[1] = Shape6[2]*Shape6[3]; Shape3[2] = Shape6[4]*Shape6[5];
   Tensor0.reshape(Shape3);
   MPS0.set(pos, Tensor0);
  }
  MPS0.bringIntoNormalShape();
  MPS0.setD(DMPS0);
  epsAchieved = 0.0;
  numSweepsDone = 0;
  return;
 }
// 2. bring MPS0 into normal form from right to left:
// 2.a) MPS0.D == 1:
 if (DMPS0 == 1)
 {
  tensorFrobeniusNormProduct = 1.0;
  for (int pos = L-1; pos > 0; pos--)
  {
   MPS0.get(pos, Tensor0);
   tensorFrobeniusNormProduct *= Tensor0.frobeniusNormalize();
   MPS0.set(pos, Tensor0);
  }
  MPS0.get(0, Tensor0);
  Tensor0.multiply(tensorFrobeniusNormProduct);
  MPS0.set(0, Tensor0);
 }
// 2.b) MPS0.D != 1:
 else
 {
  Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
  Indices0[0] = 1;
  Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
  MPS0.get(L-1, Tensor0);
  for (int pos = L-1; pos > 0; pos--)
  {
   Tensor0.LQDecompose(Indices1, Indices2, TensorL);
   MPS0.set(pos, Tensor0);
   MPS0.get(pos-1, Tensor0);
   Tensor0.contract(Indices0, TensorL, Indices1);
   Tensor0.permute(Order3);
  }
  MPS0.set(0, Tensor0);
 }
// 3. construct right half of b for each site:
 bTensors = vector< Tensor<T> >(L);
 Tensor0 = BoundaryTensors0[L-1];
 BoundaryTensors[L-1].complexConjugate(Tensor1);
 Indices0[0] = 3; Indices1[0] = 3;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor0.getShape(Shape6);
 MPS0.get(L-1, Tensor1);
 Tensor1.getShape(Shape3);
 Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape6[2]; Shape4[3] = Shape6[5];
 Tensor1.reshape(Shape4);
 Tensor1.complexConjugate();
 Indices02[0] = 2; Indices02[1] = 5; Indices12[0] = 2; Indices12[1] = 3;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 Tensor0.getShape(Shape6);
 Shape3[0] = Shape6[0]; Shape3[1] = Shape6[2]; Shape3[2] = Shape6[4];
 Tensor0.reshape(Shape3);
 bTensors[L-1] = Tensor0;
 for (int pos = L-2; pos > 0; pos--)
 {
  Tensor1 = BoundaryTensors0[pos];
  Indices0[0] = 0; Indices1[0] = 1;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  BoundaryTensors[pos].complexConjugate(Tensor1);
  Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 3;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape5);
  MPS0.get(pos, Tensor1);
  Tensor1.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[4];
  Tensor1.reshape(Shape4);
  Tensor1.complexConjugate();
  Indices03[0] = 0; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  bTensors[pos] = Tensor0;
 }
// 4. compute initial error:
 Tensor1 = BoundaryTensors0[0];
 Indices0[0] = 0; Indices1[0] = 1;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 BoundaryTensors[0].complexConjugate(Tensor1);
 Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 3;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 Tensor0.getShape(Shape5);
 MPS0.get(0, Tensor1);
 Tensor1.getShape(Shape3);
 Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[4];
 Tensor1.reshape(Shape4);
 Tensor1.complexConjugate();
 Indices03[0] = 0; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
 Tensor0.contract(Indices03, Tensor1, Indices13);
 bValue = Tensor0.get(0);
 MPS0.get(0, Tensor0);
 NValue = Tensor0.scalarProduct();
 error0 = MathAuxiliary::convertToDouble(NValue)-2.0*MathAuxiliary::convertToDouble(bValue);
// 5. sweep:
 numSweepsDone = 0;
 SweepDirection = "right";
 while (numSweepsDone < maxNumSweeps)
 {
// 5.a) sweep to the right:
  if (SweepDirection == "right")
  {
// 5.a)1. position 0:
// 5.a)1. - compute bTensor:
   Tensor0 = bTensors[1];
   Tensor1 = BoundaryTensors0[0];
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   BoundaryTensors[0].complexConjugate(Tensor1);
   Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
// 5.a)1. - update tensor:
   Tensor0.getShape(Shape5);
   Shape3[0] = 1; Shape3[1] = Shape5[0]; Shape3[2] = Shape5[2]*Shape5[4];
   Tensor0.reshape(Shape3);
   MPS0.set(0, Tensor0);
// 5.a)1. - compute error:
   NValue = Tensor0.scalarProduct();
   error1 = -MathAuxiliary::convertToDouble(NValue);
   epsAchieved = abs((error1-error0)/error1);
   if (epsAchieved < eps)
    return;
   error0 = error1;
// 5.a)1. - normalize tensor:
   if (DMPS0 == 1)
    tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
   else
   {
    Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
    Tensor0.QRDecompose(Indices2, Indices1, TensorR);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
   }
   MPS0.set(0, Tensor0);
// 5.a)1. - set bTensors:
   Tensor0 = BoundaryTensors0[0];
   BoundaryTensors[0].complexConjugate(Tensor1);
   Indices0[0] = 3; Indices1[0] = 3;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.getShape(Shape6);
   MPS0.get(0, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape6[2]; Shape4[3] = Shape6[5];
   Tensor1.reshape(Shape4);
   Tensor1.complexConjugate();
   Indices02[0] = 2; Indices02[1] = 5; Indices12[0] = 2; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape6);
   Shape3[0] = Shape6[1]; Shape3[1] = Shape6[3]; Shape3[2] = Shape6[5];
   Tensor0.reshape(Shape3);
   bTensors[0] = Tensor0;
// 5.a)2. 0 < position < L-1:
   for (int pos = 1; pos < L-1; pos++)
   {
// 5.a)2. - compute bTensor:
    Tensor0 = bTensors[pos-1];
    Tensor1 = BoundaryTensors0[pos];
    Indices0[0] = 0; Indices1[0] = 0;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    BoundaryTensors[pos].complexConjugate(Tensor1);
    Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 3;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    bTensor = Tensor0;
// 5.a)2. - update tensor:
    Tensor1 = bTensors[pos+1];
    Indices02[0] = 1; Indices02[1] = 3; Indices12[0] = 0; Indices12[1] = 1;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    Tensor0.getShape(Shape4);
    Shape3[0] = Shape4[0]; Shape3[1] = Shape4[1]*Shape4[2]; Shape3[2] = Shape4[3];
    Tensor0.reshape(Shape3);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
    MPS0.set(pos, Tensor0);
// 5.a)2. - compute error:
    NValue = Tensor0.scalarProduct();
    error1 = -MathAuxiliary::convertToDouble(NValue);
    epsAchieved = abs((error1-error0)/error1);
    if (epsAchieved < eps)
     return;
    error0 = error1;
// 5.a)2. - normalize tensor:
    if (DMPS0 == 1)
     tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
    else
    {
     Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
     Tensor0.QRDecompose(Indices2, Indices1, TensorR);
     Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
     Tensor0.permute(Order3);
    }
    MPS0.set(pos, Tensor0);
// 5.a)2. - set bTensors:
    bTensor.getShape(Shape5);
    Tensor0.getShape(Shape3);
    Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[4];
    Tensor0.reshape(Shape4);
    Tensor0.complexConjugate();
    Indices03[0] = 0; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
    bTensor.contract(Indices03, Tensor0, Indices13);
    bTensors[pos] = bTensor;
   }
// 5.a)3. position L-1:
   MPS0.get(L-1, Tensor0);
   if (DMPS0 == 1)
    Tensor0.multiply(tensorFrobeniusNorm);
   else
   {
    Indices0[0] = 1; Indices1[0] = 0;
    TensorR.contract(Indices0, Tensor0, Indices1);
    Tensor0 = TensorR;
   }
   MPS0.set(L-1, Tensor0);
   SweepDirection = "left";
  }
// 5.b) sweep to the left:
  else if (SweepDirection == "left")
  {
// 5.b)1. position L-1:
// 5.b)1. - compute bTensor:
   Tensor0 = bTensors[L-2];
   Tensor1 = BoundaryTensors0[L-1];
   Indices0[0] = 0; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   BoundaryTensors[L-1].complexConjugate(Tensor1);
   Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
// 5.b)1. - update tensor:
   Tensor0.getShape(Shape5);
   Shape3[0] = Shape5[0]; Shape3[1] = 1; Shape3[2] = Shape5[2]*Shape5[4];
   Tensor0.reshape(Shape3);
   MPS0.set(L-1, Tensor0);
// 5.b)1. - compute error:
   NValue = Tensor0.scalarProduct();
   error1 = -MathAuxiliary::convertToDouble(NValue);
   epsAchieved = abs((error1-error0)/error1);
   if (epsAchieved < eps)
    return;
   error0 = error1;
// 5.b)1. - normalize tensor:
   if (DMPS0 == 1)
    tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
   else
   {
    Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
    Tensor0.LQDecompose(Indices1, Indices2, TensorL);
   }
   MPS0.set(L-1, Tensor0);
// 5.b)1. - set bTensors:
   Tensor0 = BoundaryTensors0[L-1];
   BoundaryTensors[L-1].complexConjugate(Tensor1);
   Indices0[0] = 3; Indices1[0] = 3;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.getShape(Shape6);
   MPS0.get(L-1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape6[2]; Shape4[3] = Shape6[5];
   Tensor1.reshape(Shape4);
   Tensor1.complexConjugate();
   Indices02[0] = 2; Indices02[1] = 5; Indices12[0] = 2; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape6);
   Shape3[0] = Shape6[0]; Shape3[1] = Shape6[2]; Shape3[2] = Shape6[4];
   Tensor0.reshape(Shape3);
   bTensors[L-1] = Tensor0;
// 5.b)2. 0 < position < L-1:
   for (int pos = L-2; pos > 0; pos--)
   {
// 5.b)2. - compute bTensor:
    Tensor0 = bTensors[pos+1];
    Tensor1 = BoundaryTensors0[pos];
    Indices0[0] = 0; Indices1[0] = 1;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    BoundaryTensors[pos].complexConjugate(Tensor1);
    Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 3;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    bTensor = Tensor0;
// 5.b)2. - update tensor:
    Tensor1 = bTensors[pos-1];
    Indices02[0] = 1; Indices02[1] = 3; Indices12[0] = 0; Indices12[1] = 1;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    Tensor0.getShape(Shape4);
    Shape3[0] = Shape4[0]; Shape3[1] = Shape4[1]*Shape4[2]; Shape3[2] = Shape4[3];
    Tensor0.reshape(Shape3);
    Order3[0] = 2; Order3[1] = 0; Order3[2] = 1;
    Tensor0.permute(Order3);
    MPS0.set(pos, Tensor0);
// 5.b)2. - compute error:
    NValue = Tensor0.scalarProduct();
    error1 = -MathAuxiliary::convertToDouble(NValue);
    epsAchieved = abs((error1-error0)/error1);
    if (epsAchieved < eps)
     return;
    error0 = error1;
// 5.b)2. - normalize tensor:
    if (DMPS0 == 1)
     tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
    else
    {
     Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
     Tensor0.LQDecompose(Indices1, Indices2, TensorL);
    }
    MPS0.set(pos, Tensor0);
// 5.b)2. - set bTensors:
    bTensor.getShape(Shape5);
    Tensor0.getShape(Shape3);
    Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[4];
    Tensor0.reshape(Shape4);
    Tensor0.complexConjugate();
    Indices03[0] = 0; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
    bTensor.contract(Indices03, Tensor0, Indices13);
    bTensors[pos] = bTensor;
   }
// 5.b)3. position 0:
   MPS0.get(0, Tensor0);
   if (DMPS0 == 1)
    Tensor0.multiply(tensorFrobeniusNorm);
   else
   {
    Indices0[0] = 1; Indices1[0] = 0;
    Tensor0.contract(Indices0, TensorL, Indices1);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
   }
   MPS0.set(0, Tensor0);
   SweepDirection = "right";
  }
  numSweepsDone++;
 }
}

template<class T> void PEPS<T>::getMPO(const PEPS<T>& PEPS0, const string& Direction,
                                       unsigned int position, MPO<T>& MPO0) const
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (this->BC != PEPS0.BC) || (this->N != PEPS0.N) || (this->d != PEPS0.d) ||
     ((Direction != "right") && (Direction != "down") && (Direction != "left") &&
      (Direction != "up")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getMPO(const PEPS<T>& PEPS0, const string& Direction, " <<
                 "unsigned int position, MPO<T>& MPO0) const: " <<
          "((this->N[0] == 0) || (this->BC != PEPS0.BC) || (this->N != PEPS0.N) || " <<
           "(this->d != PEPS0.d) || " <<
           "((Direction != right) && (Direction != down) && (Direction != left) && " <<
            "(Direction != up)))." << endl;
  exit(1);
 }
 else if (((position > this->N[1]-1) && ((Direction == "right") || (Direction == "left"))) ||
          ((position > this->N[0]-1) && ((Direction == "down") || (Direction == "up"))))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getMPO(const PEPS<T>& PEPS0, const string& Direction, " <<
                 "unsigned int position, MPO<T>& MPO0) const: " <<
          "(((position > this->N[1]-1) && ((Direction == right) || (Direction == left))) || " <<
           "((position > this->N[0]-1) && ((Direction == down) || (Direction == up))))." << endl;
  exit(1);
 }
#endif
 string BCMPO0 = this->BC;
 unsigned int NMPO0;
 unsigned int dMPO0 = PEPS0.D*this->D, DMPO0 = PEPS0.D*this->D;
 vector<unsigned int> Indices0(1), Indices1(1);
 Indices0[0] = 4; Indices1[0] = 4;
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Shape0(4), Order0(8);
 T element;
 if (Direction == "right")
 {
  NMPO0 = this->N[0];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order0[0] = 3; Order0[1] = 7; Order0[2] = 1; Order0[3] = 5;
  Order0[4] = 0; Order0[5] = 4; Order0[6] = 2; Order0[7] = 6;
// the left end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = 1; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = PEPS0.Tensors[this->N[0]*(position+1)-1];
  Tensor1 = this->Tensors[this->N[0]*(position+1)-1];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(0, Tensor0);
// the intermediate part of MPO0:
  Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  for (int i = 1; i < NMPO0-1; i++)
  {
   Tensor0 = PEPS0.Tensors[this->N[0]*(position+1)-1-i];
   Tensor1 = this->Tensors[this->N[0]*(position+1)-1-i];
   Tensor1.complexConjugate();
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPO0.set(i, Tensor0);
  }
// the right end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = DMPO0; Shape0[1] = 1; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = PEPS0.Tensors[this->N[0]*position];
  Tensor1 = this->Tensors[this->N[0]*position];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(NMPO0-1, Tensor0);
 }
 else if (Direction == "down")
 {
  NMPO0 = this->N[1];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order0[0] = 0; Order0[1] = 4; Order0[2] = 2; Order0[3] = 6;
  Order0[4] = 1; Order0[5] = 5; Order0[6] = 3; Order0[7] = 7;
// the left end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = 1; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = PEPS0.Tensors[position];
  Tensor1 = this->Tensors[position];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(0, Tensor0);
// the intermediate part of MPO0:
  Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  for (int i = 1; i < NMPO0-1; i++)
  {
   Tensor0 = PEPS0.Tensors[position+this->N[0]*i];
   Tensor1 = this->Tensors[position+this->N[0]*i];
   Tensor1.complexConjugate();
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPO0.set(i, Tensor0);
  }
// the right end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = DMPO0; Shape0[1] = 1; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = PEPS0.Tensors[position+this->N[0]*(this->N[1]-1)];
  Tensor1 = this->Tensors[position+this->N[0]*(this->N[1]-1)];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(NMPO0-1, Tensor0);
 }
 else if (Direction == "left")
 {
  NMPO0 = this->N[0];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order0[0] = 1; Order0[1] = 5; Order0[2] = 3; Order0[3] = 7;
  Order0[4] = 2; Order0[5] = 6; Order0[6] = 0; Order0[7] = 4;
// the left end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = 1; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = PEPS0.Tensors[this->N[0]*position];
  Tensor1 = this->Tensors[this->N[0]*position];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(0, Tensor0);
// the intermediate part of MPO0:
  Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  for (int i = 1; i < NMPO0-1; i++)
  {
   Tensor0 = PEPS0.Tensors[i+this->N[0]*position];
   Tensor1 = this->Tensors[i+this->N[0]*position];
   Tensor1.complexConjugate();
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPO0.set(i, Tensor0);
  }
// the right end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = DMPO0; Shape0[1] = 1; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = PEPS0.Tensors[this->N[0]*(position+1)-1];
  Tensor1 = this->Tensors[this->N[0]*(position+1)-1];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(NMPO0-1, Tensor0);
 }
 else if (Direction == "up")
 {
  NMPO0 = this->N[1];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order0[0] = 2; Order0[1] = 6; Order0[2] = 0; Order0[3] = 4;
  Order0[4] = 3; Order0[5] = 7; Order0[6] = 1; Order0[7] = 5;
// the left end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = 1; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = PEPS0.Tensors[position+this->N[0]*(this->N[1]-1)];
  Tensor1 = this->Tensors[position+this->N[0]*(this->N[1]-1)];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(0, Tensor0);
// the intermediate part of MPO0:
  Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  for (int i = 1; i < NMPO0-1; i++)
  {
   Tensor0 = PEPS0.Tensors[position+this->N[0]*(this->N[1]-1-i)];
   Tensor1 = this->Tensors[position+this->N[0]*(this->N[1]-1-i)];
   Tensor1.complexConjugate();
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPO0.set(i, Tensor0);
  }
// the right end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = DMPO0; Shape0[1] = 1; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = PEPS0.Tensors[position];
  Tensor1 = this->Tensors[position];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(NMPO0-1, Tensor0);
 }
}

template<class T> void PEPS<T>::multiplyBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, unsigned int position,
                                                   const MPS<T>& MPS0,
                                                   double eps, unsigned int maxNumSweeps,
                                                   double& epsAchieved, unsigned int& numSweepsDone,
                                                   MPS<T>& MPS1) const
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, unsigned int position, const MPS<T>& MPS0, " <<
                             "double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "(this->N[0] == 0)." << endl;
  exit(1);
 }
 if ((this->BC != PEPS0.BC) || (this->N != PEPS0.N) || (this->d != PEPS0.d))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, unsigned int position, const MPS<T>& MPS0, " <<
                             "double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((this->BC != PEPS0.BC) || (this->N != PEPS0.N) || (this->d != PEPS0.d))." << endl;
  exit(1);
 }
 if ((Direction != "right") && (Direction != "down") && (Direction != "left") && (Direction != "up"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, unsigned int position, const MPS<T>& MPS0, " <<
                             "double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((Direction != right) && (Direction != down) && (Direction != left) && (Direction != up))." << endl;
  exit(1);
 }
 if ((((Direction == "right") || (Direction == "left")) && (position > this->N[1]-1)) ||
     (((Direction == "down") || (Direction == "up")) && (position > this->N[0]-1)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, unsigned int position, const MPS<T>& MPS0, " <<
                             "double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((((Direction == right) || (Direction == left)) && (position > this->N[1]-1)) || " <<
           "(((Direction == down) || (Direction == up)) && (position > this->N[0]-1)))." << endl;
  exit(1);
 }
 string BCMPS;
 MPS0.getBC(BCMPS);
 if ((BCMPS != "open") ||
     (((Direction == "right") || (Direction == "left")) && (MPS0.getN() != this->N[0])) ||
     (((Direction == "down") || (Direction == "up")) && (MPS0.getN() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, unsigned int position, const MPS<T>& MPS0, " <<
                             "double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((MPS0.BC != open) || (((Direction == right) || (Direction == left)) && (MPS0.N != this->N[0])) || " <<
           "(((Direction == down) || (Direction == up)) && (MPS0.N != this->N[1])))." << endl;
  exit(1);
 }
 MPS1.getBC(BCMPS);
 if ((BCMPS != "open") ||
     (((Direction == "right") || (Direction == "left")) && (MPS1.getN() != this->N[0])) ||
     (((Direction == "down") || (Direction == "up")) && (MPS1.getN() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, unsigned int position, const MPS<T>& MPS0, " <<
                             "double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((MPS1.BC != open) || (((Direction == right) || (Direction == left)) && (MPS1.N != this->N[0])) || " <<
           "(((Direction == down) || (Direction == up)) && (MPS1.N != this->N[1])))." << endl;
  exit(1);
 }
#endif
 unsigned int DMPS1, L;
 double error0, error1;
 T bValue, NValue, tensorFrobeniusNorm, tensorFrobeniusNormProduct;
 string SweepDirection;
 vector<unsigned int> dMPS, Indices0(1), Indices02(2), Indices03(3), Indices1(1), Indices12(2), Indices13(3), Indices2(2);
 vector<unsigned int> Order3(3), Order5(5), Order8(8), Shape3(3), Shape4(4), Shape5(5), Shape05(5), Shape6(6), Shape8(8);
 Tensor<T> bTensor, Tensor0, Tensor1, TensorL, TensorR;
 vector< Tensor<T> > BulkTensors, BulkTensors0, bTensors;
// 0. initialize:
// 0. - set length L:
 if ((Direction == "right") || (Direction == "left"))
  L = this->N[0];
 else if ((Direction == "down") || (Direction == "up"))
  L = this->N[1];
// 0. - get bulk-tensors clockwise, permute them as Direction "down", and write the ones from this PEPS in BulkTensors and
//      the ones from PEPS0 in BulkTensors0:
 BulkTensors = vector< Tensor<T> >(L);
 BulkTensors0 = vector< Tensor<T> >(L);
 if (Direction == "right")
 {
  Order5[0] = 3; Order5[1] = 0; Order5[2] = 1; Order5[3] = 2; Order5[4] = 4;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[this->N[0]-1-pos+this->N[0]*position];
   Tensor0.permute(Order5);
   BulkTensors[pos] = Tensor0;
   Tensor0 = PEPS0.Tensors[this->N[0]-1-pos+this->N[0]*position];
   Tensor0.permute(Order5);
   BulkTensors0[pos] = Tensor0;
  }
 }
 else if (Direction == "down")
 {
  for (int pos = 0; pos < L; pos++)
  {
   BulkTensors[pos] = this->Tensors[position+this->N[0]*pos];
   BulkTensors0[pos] = PEPS0.Tensors[position+this->N[0]*pos];
  }
 }
 else if (Direction == "left")
 {
  Order5[0] = 1; Order5[1] = 2; Order5[2] = 3; Order5[3] = 0; Order5[4] = 4;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[pos+this->N[0]*position];
   Tensor0.permute(Order5);
   BulkTensors[pos] = Tensor0;
   Tensor0 = PEPS0.Tensors[pos+this->N[0]*position];
   Tensor0.permute(Order5);
   BulkTensors0[pos] = Tensor0;
  }
 }
 else if (Direction == "up")
 {
  Order5[0] = 2; Order5[1] = 3; Order5[2] = 0; Order5[3] = 1; Order5[4] = 4;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[position+this->N[0]*(this->N[1]-1-pos)];
   Tensor0.permute(Order5);
   BulkTensors[pos] = Tensor0;
   Tensor0 = PEPS0.Tensors[position+this->N[0]*(this->N[1]-1-pos)];
   Tensor0.permute(Order5);
   BulkTensors0[pos] = Tensor0;
  }
 }
// 0. - check MPS0.d:
#ifdef DEBUG
 MPS0.getd(dMPS);
 for (int pos = 0; pos < L; pos++)
 {
  BulkTensors0[pos].getShape(Shape05);
  BulkTensors[pos].getShape(Shape5);
  if (dMPS[pos] != Shape05[1]*Shape5[1])
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void PEPS<T>::" <<
           "multiplyBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, unsigned int position, const MPS<T>& MPS0, " <<
                              "double eps, unsigned int maxNumSweeps, " <<
                              "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
           "(MPS0.d[" << pos << "] != PEPS0.D[" << pos << "]*this->D[" << pos << "])." << endl;
   exit(1);
  }
 }
#endif
// 0. - check MPS1.d:
 MPS1.getd(dMPS);
#ifdef DEBUG
 for (int pos = 0; pos < L; pos++)
 {
  BulkTensors0[pos].getShape(Shape05);
  BulkTensors[pos].getShape(Shape5);
  if (dMPS[pos] != Shape05[3]*Shape5[3])
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void PEPS<T>::" <<
           "multiplyBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, unsigned int position, const MPS<T>& MPS0, " <<
                              "double eps, unsigned int maxNumSweeps, " <<
                              "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
           "(MPS1.d[" << pos << "] != PEPS0.D[" << pos << "]*this->D[" << pos << "])." << endl;
   exit(1);
  }
 }
#endif
// 1. MPS1.D >= PEPS0.D*this->D*MPS0.D:
 DMPS1 = MPS1.getD();
 if (DMPS1 >= PEPS0.D*this->D*MPS0.getD())
 {
  MPS1 = MPS<T>("open", L, dMPS, PEPS0.D*this->D*MPS0.getD());
  Indices0[0] = 2; Indices1[0] = 1;
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 4;
  Order8[0] = 0; Order8[1] = 2; Order8[2] = 5; Order8[3] = 1; Order8[4] = 3; Order8[5] = 6; Order8[6] = 4; Order8[7] = 7;
  for (int pos = 0; pos < L; pos++)
  {
   MPS0.get(pos, Tensor0);
   Tensor0.getShape(Shape3);
   BulkTensors0[pos].getShape(Shape05);
   BulkTensors[pos].getShape(Shape5);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[1]; Shape4[3] = Shape5[1];
   Tensor0.reshape(Shape4);
   Tensor1 = BulkTensors0[pos];
   Tensor0.contract(Indices0, Tensor1, Indices1);
   BulkTensors[pos].complexConjugate(Tensor1);
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.permute(Order8);
   Tensor0.getShape(Shape8);
   Shape3[0] = Shape8[0]*Shape8[1]*Shape8[2]; Shape3[1] = Shape8[3]*Shape8[4]*Shape8[5]; Shape3[2] = Shape8[6]*Shape8[7];
   Tensor0.reshape(Shape3);
   MPS1.set(pos, Tensor0);
  }
  MPS1.bringIntoNormalShape();
  MPS1.setD(DMPS1);
  epsAchieved = 0.0;
  numSweepsDone = 0;
  return;
 }
// 2. bring MPS1 into normal form from right to left:
// 2.a) MPS1.D == 1:
 if (DMPS1 == 1)
 {
  tensorFrobeniusNormProduct = 1.0;
  for (int pos = L-1; pos > 0; pos--)
  {
   MPS1.get(pos, Tensor0);
   tensorFrobeniusNormProduct *= Tensor0.frobeniusNormalize();
   MPS1.set(pos, Tensor0);
  }
  MPS1.get(0, Tensor0);
  Tensor0.multiply(tensorFrobeniusNormProduct);
  MPS1.set(0, Tensor0);
 }
// 2.b) MPS1.D != 1:
 else
 {
  Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
  Indices0[0] = 1;
  Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
  MPS1.get(L-1, Tensor0);
  for (int pos = L-1; pos > 0; pos--)
  {
   Tensor0.LQDecompose(Indices1, Indices2, TensorL);
   MPS1.set(pos, Tensor0);
   MPS1.get(pos-1, Tensor0);
   Tensor0.contract(Indices0, TensorL, Indices1);
   Tensor0.permute(Order3);
  }
  MPS1.set(0, Tensor0);
 }
// 3. construct right half of b for each site:
 bTensors = vector< Tensor<T> >(L);
 MPS0.get(L-1, Tensor0);
 Tensor0.getShape(Shape3);
 BulkTensors0[L-1].getShape(Shape05);
 BulkTensors[L-1].getShape(Shape5);
 Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[1]; Shape4[3] = Shape5[1];
 Tensor0.reshape(Shape4);
 Tensor1 = BulkTensors0[L-1];
 Indices0[0] = 2; Indices1[0] = 1;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 BulkTensors[L-1].complexConjugate(Tensor1);
 Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 4;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 MPS1.get(L-1, Tensor1);
 Tensor1.getShape(Shape3);
 Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[3]; Shape4[3] = Shape5[3];
 Tensor1.reshape(Shape4);
 Tensor1.complexConjugate();
 Indices02[0] = 4; Indices02[1] = 7; Indices12[0] = 2; Indices12[1] = 3;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 Tensor0.getShape(Shape8);
 Shape4[0] = Shape8[0]; Shape4[1] = Shape8[2]; Shape4[2] = Shape8[4]; Shape4[3] = Shape8[6];
 Tensor0.reshape(Shape4);
 bTensors[L-1] = Tensor0;
 Indices0[0] = 0; Indices1[0] = 1;
 Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 2; Indices12[1] = 1;
 Indices03[0] = 0; Indices03[1] = 3;
 for (int pos = L-2; pos > 0; pos--)
 {
  MPS0.get(pos, Tensor1);
  Tensor1.getShape(Shape3);
  BulkTensors0[pos].getShape(Shape05);
  BulkTensors[pos].getShape(Shape5);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[1]; Shape4[3] = Shape5[1];
  Tensor1.reshape(Shape4);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor1 = BulkTensors0[pos];
  Tensor0.contract(Indices02, Tensor1, Indices12);
  BulkTensors[pos].complexConjugate(Tensor1);
  Indices03[2] = 6; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  MPS1.get(pos, Tensor1);
  Tensor1.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[3]; Shape4[3] = Shape5[3];
  Tensor1.reshape(Shape4);
  Tensor1.complexConjugate();
  Indices03[2] = 5; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  bTensors[pos] = Tensor0;
 }
// 4. compute initial error:
 MPS0.get(0, Tensor1);
 Tensor1.getShape(Shape3);
 BulkTensors0[0].getShape(Shape05);
 BulkTensors[0].getShape(Shape5);
 Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[1]; Shape4[3] = Shape5[1];
 Tensor1.reshape(Shape4);
 Indices0[0] = 0; Indices1[0] = 1;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = BulkTensors0[0];
 Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 2; Indices12[1] = 1;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 BulkTensors[0].complexConjugate(Tensor1);
 Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
 Tensor0.contract(Indices03, Tensor1, Indices13);
 MPS1.get(0, Tensor1);
 Tensor1.getShape(Shape3);
 Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[3]; Shape4[3] = Shape5[3];
 Tensor1.reshape(Shape4);
 Tensor1.complexConjugate();
 Indices03[2] = 5; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
 Tensor0.contract(Indices03, Tensor1, Indices13);
 bValue = Tensor0.get(0);
 MPS1.get(0, Tensor0);
 NValue = Tensor0.scalarProduct();
 error0 = MathAuxiliary::convertToDouble(NValue)-2.0*MathAuxiliary::convertToDouble(bValue);
// 5. sweep:
 numSweepsDone = 0;
 SweepDirection = "right";
 while (numSweepsDone < maxNumSweeps)
 {
// 5.a) sweep to the right:
  if (SweepDirection == "right")
  {
// 5.a)1. position 0:
// 5.a)1. - compute bTensor:
   Tensor0 = bTensors[1];
   MPS0.get(0, Tensor1);
   Tensor1.getShape(Shape3);
   BulkTensors0[0].getShape(Shape05);
   BulkTensors[0].getShape(Shape5);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[1]; Shape4[3] = Shape5[1];
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BulkTensors0[0];
   Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 2; Indices12[1] = 1;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   BulkTensors[0].complexConjugate(Tensor1);
   Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
// 5.a)1. - update tensor:
   Tensor0.getShape(Shape6);
   Shape3[0] = 1; Shape3[1] = Shape6[0]; Shape3[2] = Shape6[3]*Shape6[5];
   Tensor0.reshape(Shape3);
   MPS1.set(0, Tensor0);
// 5.a)1. - compute error:
   NValue = Tensor0.scalarProduct();
   error1 = -MathAuxiliary::convertToDouble(NValue);
   epsAchieved = abs((error1-error0)/error1);
   if (epsAchieved < eps)
    return;
   error0 = error1;
// 5.a)1. - normalize tensor:
   if (DMPS1 == 1)
    tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
   else
   {
    Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
    Tensor0.QRDecompose(Indices2, Indices1, TensorR);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
   }
   MPS1.set(0, Tensor0);
// 5.a)1. - set bTensors:
   MPS0.get(0, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[1]; Shape4[3] = Shape5[1];
   Tensor0.reshape(Shape4);
   Tensor1 = BulkTensors0[0];
   Indices0[0] = 2; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   BulkTensors[0].complexConjugate(Tensor1);
   Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 4;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   MPS1.get(0, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[3]; Shape4[3] = Shape5[3];
   Tensor1.reshape(Shape4);
   Tensor1.complexConjugate();
   Indices02[0] = 4; Indices02[1] = 7; Indices12[0] = 2; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape8);
   Shape4[0] = Shape8[1]; Shape4[1] = Shape8[3]; Shape4[2] = Shape8[5]; Shape4[3] = Shape8[7];
   Tensor0.reshape(Shape4);
   bTensors[0] = Tensor0;
// 5.a)2. 0 < pos < L-1:
   for (int pos = 1; pos < L-1; pos++)
   {
// 5.a)2. - compute bTensor:
    Tensor0 = bTensors[pos-1];
    MPS0.get(pos, Tensor1);
    Tensor1.getShape(Shape3);
    BulkTensors0[pos].getShape(Shape05);
    BulkTensors[pos].getShape(Shape5);
    Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[1]; Shape4[3] = Shape5[1];
    Tensor1.reshape(Shape4);
    Indices0[0] = 0; Indices1[0] = 0;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = BulkTensors0[pos];
    Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 1;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    BulkTensors[pos].complexConjugate(Tensor1);
    Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 4;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    bTensor = Tensor0;
// 5.a)2. - update tensor:
    Tensor1 = bTensors[pos+1];
    Indices03[0] = 1; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 2;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    Tensor0.getShape(Shape4);
    Shape3[0] = Shape4[0]; Shape3[1] = Shape4[1]*Shape4[2]; Shape3[2] = Shape4[3];
    Tensor0.reshape(Shape3);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
    MPS1.set(pos, Tensor0);
// 5.a)2. - compute error:
    NValue = Tensor0.scalarProduct();
    error1 = -MathAuxiliary::convertToDouble(NValue);
    epsAchieved = abs((error1-error0)/error1);
    if (epsAchieved < eps)
     return;
    error0 = error1;
// 5.a)2. - normalize tensor:
    if (DMPS1 == 1)
     tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
    else
    {
     Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
     Tensor0.QRDecompose(Indices2, Indices1, TensorR);
     Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
     Tensor0.permute(Order3);
    }
    MPS1.set(pos, Tensor0);
// 5.a)2. - set bTensors:
    Tensor0.getShape(Shape3);
    Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[3]; Shape4[3] = Shape5[3];
    Tensor0.reshape(Shape4);
    Tensor0.complexConjugate();
    Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 5; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
    bTensor.contract(Indices03, Tensor0, Indices13);
    bTensors[pos] = bTensor;
   }
// 5.a)3. position L-1:
   MPS1.get(L-1, Tensor0);
   if (DMPS1 == 1)
    Tensor0.multiply(tensorFrobeniusNorm);
   else
   {
    Indices0[0] = 1; Indices1[0] = 0;
    TensorR.contract(Indices0, Tensor0, Indices1);
    Tensor0 = TensorR;
   }
   MPS1.set(L-1, Tensor0);
   SweepDirection = "left";
  }
// 5.b) sweep to the left:
  else if (SweepDirection == "left")
  {
// 5.b)1. position L-1:
// 5.b)1. - compute bTensor:
   Tensor0 = bTensors[L-2];
   MPS0.get(L-1, Tensor1);
   Tensor1.getShape(Shape3);
   BulkTensors0[L-1].getShape(Shape05);
   BulkTensors[L-1].getShape(Shape5);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape05[1]; Shape4[3] = Shape5[1];
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BulkTensors0[L-1];
   Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 1;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   BulkTensors[L-1].complexConjugate(Tensor1);
   Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
// 5.b)1. - update tensor:
   Tensor0.getShape(Shape6);
   Shape3[0] = Shape6[0]; Shape3[1] = 1; Shape3[2] = Shape6[3]*Shape6[5];
   Tensor0.reshape(Shape3);
   MPS1.set(L-1, Tensor0);
// 5.b)1. - compute error:
   NValue = Tensor0.scalarProduct();
   error1 = -MathAuxiliary::convertToDouble(NValue);
   epsAchieved = abs((error1-error0)/error1);
   if (epsAchieved < eps)
    return;
   error0 = error1;
// 5.b)1. - normalize tensor:
   if (DMPS1 == 1)
    tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
   else
   {
    Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
    Tensor0.LQDecompose(Indices1, Indices2, TensorL);
   }
   MPS1.set(L-1, Tensor0);
// 5.b)1. - set bTensors:
   MPS0.get(L-1, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape05[1]; Shape4[3] = Shape5[1];
   Tensor0.reshape(Shape4);
   Tensor1 = BulkTensors0[L-1];
   Indices0[0] = 2; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   BulkTensors[L-1].complexConjugate(Tensor1);
   Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 4;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   MPS1.get(L-1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape05[3]; Shape4[3] = Shape5[3];
   Tensor1.reshape(Shape4);
   Tensor1.complexConjugate();
   Indices02[0] = 4; Indices02[1] = 7; Indices12[0] = 2; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape8);
   Shape4[0] = Shape8[0]; Shape4[1] = Shape8[2]; Shape4[2] = Shape8[4]; Shape4[3] = Shape8[6];
   Tensor0.reshape(Shape4);
   bTensors[L-1] = Tensor0;
// 5.b)2. 0 < pos < L-1:
   for (int pos = L-2; pos > 0; pos--)
   {
// 5.b)2. - compute bTensor:
    Tensor0 = bTensors[pos+1];
    MPS0.get(pos, Tensor1);
    Tensor1.getShape(Shape3);
    BulkTensors0[pos].getShape(Shape05);
    BulkTensors[pos].getShape(Shape5);
    Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[1]; Shape4[3] = Shape5[1];
    Tensor1.reshape(Shape4);
    Indices0[0] = 0; Indices1[0] = 1;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = BulkTensors0[pos];
    Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 2; Indices12[1] = 1;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    BulkTensors[pos].complexConjugate(Tensor1);
    Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    bTensor = Tensor0;
// 5.b)2. - update tensor:
    Tensor1 = bTensors[pos-1];
    Indices03[0] = 1; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 2;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    Tensor0.getShape(Shape4);
    Shape3[0] = Shape4[0]; Shape3[1] = Shape4[1]*Shape4[2]; Shape3[2] = Shape4[3];
    Tensor0.reshape(Shape3);
    Order3[0] = 2; Order3[1] = 0; Order3[2] = 1;
    Tensor0.permute(Order3);
    MPS1.set(pos, Tensor0);
// 5.b)2. - compute error:
    NValue = Tensor0.scalarProduct();
    error1 = -MathAuxiliary::convertToDouble(NValue);
    epsAchieved = abs((error1-error0)/error1);
    if (epsAchieved < eps)
     return;
    error0 = error1;
// 5.b)2. - normalize tensor:
    if (DMPS1 == 1)
     tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
    else
    {
     Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
     Tensor0.LQDecompose(Indices1, Indices2, TensorL);
    }
    MPS1.set(pos, Tensor0);
// 5.b)2. - set bTensors:
    Tensor0.getShape(Shape3);
    Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[3]; Shape4[3] = Shape5[3];
    Tensor0.reshape(Shape4);
    Tensor0.complexConjugate();
    Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 5; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
    bTensor.contract(Indices03, Tensor0, Indices13);
    bTensors[pos] = bTensor;
   }
// 5.b)3. position 0:
   MPS1.get(0, Tensor0);
   if (DMPS1 == 1)
    Tensor0.multiply(tensorFrobeniusNorm);
   else
   {
    Indices0[0] = 1; Indices1[0] = 0;
    Tensor0.contract(Indices0, TensorL, Indices1);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
   }
   MPS1.set(0, Tensor0);
   SweepDirection = "right";
  }
  numSweepsDone++;
 }
}

template<class T> T PEPS<T>::contractReverseBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction,
                                                       const MPS<T>& MPS0, unsigned int position, const MPS<T>& MPS1) const
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (PEPS0.N[0] == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "contractReverseBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, const MPS<T>& MPS0, " <<
                                    "unsigned int position, const MPS<T>& MPS1) const: " <<
          "((this->N[0] == 0) || (PEPS0.N[0] == 0))." << endl;
  exit(1);
 }
 if ((Direction != "horizontal") && (Direction != "vertical"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "contractReverseBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, const MPS<T>& MPS0, " <<
                                    "unsigned int position, const MPS<T>& MPS1) const: " <<
          "((Direction != horizontal) && (Direction != vertical))." << endl;
  exit(1);
 }
 string BCMPS;
 MPS0.getBC(BCMPS);
 if ((BCMPS != "open") ||
     ((Direction == "horizontal") && (MPS0.getN() != this->N[0])) ||
     ((Direction == "vertical") && (MPS0.getN() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "contractReverseBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, const MPS<T>& MPS0, " <<
                                    "unsigned int position, const MPS<T>& MPS1) const: " <<
          "((MPS0.BC != open) || ((Direction == horizontal) && (MPS0.N != this->N[0])) || " <<
           "((Direction == vertical) && (MPS0.N != this->N[1])))." << endl;
  exit(1);
 }
 if (((Direction == "horizontal") && (position > this->N[1]-1)) ||
     ((Direction == "vertical") && (position > this->N[0]-1)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "contractReverseBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, const MPS<T>& MPS0, " <<
                                    "unsigned int position, const MPS<T>& MPS1) const: " <<
          "(((Direction == horizontal) && (position > this->N[1]-1)) || " <<
           "((Direction == vertical) && (position > this->N[0]-1)))." << endl;
  exit(1);
 }
 MPS1.getBC(BCMPS);
 if ((BCMPS != "open") ||
     ((Direction == "horizontal") && (MPS1.getN() != this->N[0])) ||
     ((Direction == "vertical") && (MPS1.getN() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "contractReverseBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, const MPS<T>& MPS0, " <<
                                    "unsigned int position, const MPS<T>& MPS1) const: " <<
          "((MPS1.BC != open) || ((Direction == horizontal) && (MPS1.N != this->N[0])) || " <<
           "((Direction == vertical) && (MPS1.N != this->N[1])))." << endl;
  exit(1);
 }
#endif
 unsigned int L;
 vector<unsigned int> dMPS, Indices0(1), Indices02(2), Indices03(3), Indices1(1), Indices12(2), Indices13(3);
 vector<unsigned int> Order5(5), Shape3(3), Shape4(4), Shape5(5), Shape05(5);
 Tensor<T> Tensor0, Tensor1;
 vector< Tensor<T> > BulkTensors, BulkTensors0;
// 0. initialize:
// 0. - set length L:
 if (Direction == "horizontal")
  L = this->N[0];
 else if (Direction == "vertical")
  L = this->N[1];
// 0. - get bulk-tensors, permute them as Direction "vertical", and write the ones from this PEPS in BulkTensors and
//      the ones from PEPS0 in BulkTensors0:
 BulkTensors = vector< Tensor<T> >(L);
 BulkTensors0 = vector< Tensor<T> >(L);
 if (Direction == "horizontal")
 {
  Order5[0] = 1; Order5[1] = 2; Order5[2] = 3; Order5[3] = 0; Order5[4] = 4;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[pos+this->N[0]*position];
   Tensor0.permute(Order5);
   BulkTensors[pos] = Tensor0;
   Tensor0 = PEPS0.Tensors[pos+this->N[0]*position];
   Tensor0.permute(Order5);
   BulkTensors0[pos] = Tensor0;
  }
 }
 else if (Direction == "vertical")
 {
  for (int pos = 0; pos < L; pos++)
  {
   BulkTensors[pos] = this->Tensors[position+this->N[0]*pos];
   BulkTensors0[pos] = PEPS0.Tensors[position+this->N[0]*pos];
  }
 }
// 0. - check MPS0.d and MPS1.d:
#ifdef DEBUG
 MPS0.getd(dMPS);
 for (int pos = 0; pos < L; pos++)
 {
  BulkTensors0[L-1-pos].getShape(Shape05);
  BulkTensors[L-1-pos].getShape(Shape5);
  if (dMPS[pos] != Shape05[3]*Shape5[3])
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> T PEPS<T>::" <<
           "contractReverseBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, const MPS<T>& MPS0, " <<
                                     "unsigned int position, const MPS<T>& MPS1) const: " <<
           "(MPS0.d[" << pos << "] != PEPS0.D[" << L-1-pos << "]*this->D[" << L-1-pos << "])." << endl;
   exit(1);
  }
 }
 MPS1.getd(dMPS);
 for (int pos = 0; pos < L; pos++)
 {
  BulkTensors0[pos].getShape(Shape05);
  BulkTensors[pos].getShape(Shape5);
  if (dMPS[pos] != Shape05[1]*Shape5[1])
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> T PEPS<T>::" <<
           "contractReverseBulkMPOMPS(const PEPS<T>& PEPS0, const string& Direction, const MPS<T>& MPS0, " <<
                                     "unsigned int position, const MPS<T>& MPS1) const: " <<
           "(MPS1.d[" << pos << "] != PEPS0.D[" << pos << "]*this->D[" << pos << "])." << endl;
   exit(1);
  }
 }
#endif
// 1. position 0:
 MPS1.get(0, Tensor0);
 Tensor0.getShape(Shape3);
 BulkTensors0[0].getShape(Shape05);
 BulkTensors[0].getShape(Shape5);
 Shape3[0] = Shape3[1]; Shape3[1] = Shape05[1]; Shape3[2] = Shape5[1];
 Tensor0.reshape(Shape3);
 Shape4[0] = Shape05[1]; Shape4[1] = Shape05[2]; Shape4[2] = Shape05[3]; Shape4[3] = Shape05[4];
 BulkTensors0[0].reshape(Shape4);
 Indices0[0] = 1; Indices1[0] = 0;
 Tensor0.contract(Indices0, BulkTensors0[0], Indices1);
 Shape4[0] = Shape5[1]; Shape4[1] = Shape5[2]; Shape4[2] = Shape5[3]; Shape4[3] = Shape5[4];
 BulkTensors[0].reshape(Shape4);
 BulkTensors[0].complexConjugate();
 Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 3;
 Tensor0.contract(Indices02, BulkTensors[0], Indices12);
 MPS0.get(L-1, Tensor1);
 Tensor1.getShape(Shape3);
 Shape3[1] = Shape05[3]; Shape3[2] = Shape5[3];
 Tensor1.reshape(Shape3);
 Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
 Tensor0.contract(Indices02, Tensor1, Indices12);
// 2. 0 < position <= L-1:
 Indices0[0] = 0; Indices1[0] = 0;
 Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 1;
 Indices03[0] = 0; Indices03[1] = 3;
 for (int pos = 1; pos < L; pos++)
 {
  MPS1.get(pos, Tensor1);
  Tensor1.getShape(Shape3);
  BulkTensors0[pos].getShape(Shape05);
  BulkTensors[pos].getShape(Shape5);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[1]; Shape4[3] = Shape5[1];
  Tensor1.reshape(Shape4);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.contract(Indices02, BulkTensors0[pos], Indices12);
  BulkTensors[pos].complexConjugate();
  Indices03[2] = 6; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 4;
  Tensor0.contract(Indices03, BulkTensors[pos], Indices13);
  MPS0.get(L-1-pos, Tensor1);
  Tensor1.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[3]; Shape4[3] = Shape5[3];
  Tensor1.reshape(Shape4);
  Indices03[2] = 5; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
  Tensor0.contract(Indices03, Tensor1, Indices13);
 }
// 3. return contraction-value:
 return Tensor0.get(0);
}

template<class T> void PEPS<T>::getMPS(const PEPO<T>& PEPO0, const string& Direction,
                                       MPS<T>& MPS0) const
{
#ifdef DEBUG
 string BCPEPO0; PEPO0.getBC(BCPEPO0);
 vector<unsigned int> NPEPO0(2); NPEPO0[0] = PEPO0.getNrows(); NPEPO0[1] = PEPO0.getNcols();
 Matrix<unsigned int> dPEPO0; PEPO0.getd(dPEPO0);
 if ((this->N[0] == 0) || (this->BC != "open") || (this->BC != BCPEPO0) ||
     (this->N != NPEPO0) || (this->d != dPEPO0) ||
     ((Direction != "right") && (Direction != "down") && (Direction != "left") &&
      (Direction != "up")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getMPS(const PEPO<T>& PEPO0, const string& Direction, " <<
                 "MPS<T>& MPS0) const: " <<
          "((this->N[0] == 0) || (this->BC != open) || (this->BC != BCPEPO0) || " <<
           "(this->N != NPEPO0) || (this->d != dPEPO0) || " <<
           "((Direction != right) && (Direction != down) && (Direction != left) && " <<
            "(Direction != up)))." << endl;
  exit(1);
 }
#endif
// compute |PEPS0>=PEPO0*|this>:
 unsigned D0 = this->D*PEPO0.getD();
 PEPS<T> PEPS0(this->BC, this->N[0], this->N[1], this->d, D0);
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Indices0(1), Indices1(1);
 Indices0[0] = 4; Indices1[0] = 4;
 vector<unsigned int> Order(9), Shape(5);
 Order[0] = 0; Order[1] = 4; Order[2] = 1; Order[3] = 5; Order[4] = 2; Order[5] = 6; Order[6] = 3; Order[7] = 7; Order[8] = 8;
 Shape[0] = D0; Shape[1] = D0; Shape[2] = D0; Shape[3] = D0;
 for (int col = 0; col < this->N[1]; col++)
 {
  for (int row = 0; row < this->N[0]; row++)
  {
   this->get(row, col, Tensor0); PEPO0.get(row, col, Tensor1);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order);
   Shape[4] = this->d(row, col);
   if (this->BC == "open")
    PEPS0.getOpenBCShape(row, col, Shape);
   Tensor0.reshape(Shape);
   PEPS0.set(row, col, Tensor0);
  }
 }
// compute MPS0:
 string BCMPS0 = "open";
 unsigned int NMPS0;
 unsigned int dMPS0 = D0*this->D, DMPS0 = D0*this->D;
 vector<unsigned int> Shape0(3), Order0(8);
 T element;
 if (Direction == "right")
 {
  NMPS0 = this->N[0];
  MPS0 = MPS<T>(BCMPS0, NMPS0, dMPS0, DMPS0);
  Order0[0] = 3; Order0[1] = 7; Order0[2] = 0; Order0[3] = 4;
  Order0[4] = 1; Order0[5] = 5; Order0[6] = 2; Order0[7] = 6;
// the left end of MPS0:
  Shape0[0] = 1; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  Tensor0 = PEPS0.Tensors[this->N[0]-1];
  Tensor1 = this->Tensors[this->N[0]-1];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(0, Tensor0);
// the intermediate part of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  for (int i = 1; i < NMPS0-1; i++)
  {
   Tensor0 = PEPS0.Tensors[this->N[0]-1-i];
   Tensor1 = this->Tensors[this->N[0]-1-i];
   Tensor1.complexConjugate();
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPS0.set(i, Tensor0);
  }
// the right end of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = 1; Shape0[2] = dMPS0;
  Tensor0 = PEPS0.Tensors[0];
  Tensor1 = this->Tensors[0];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(NMPS0-1, Tensor0);
 }
 else if (Direction == "down")
 {
  NMPS0 = this->N[1];
  MPS0 = MPS<T>(BCMPS0, NMPS0, dMPS0, DMPS0);
  Order0[0] = 0; Order0[1] = 4; Order0[2] = 1; Order0[3] = 5;
  Order0[4] = 2; Order0[5] = 6; Order0[6] = 3; Order0[7] = 7;
// the left end of MPS0:
  Shape0[0] = 1; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  Tensor0 = PEPS0.Tensors[0];
  Tensor1 = this->Tensors[0];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(0, Tensor0);
// the intermediate part of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  for (int i = 1; i < NMPS0-1; i++)
  {
   Tensor0 = PEPS0.Tensors[this->N[0]*i];
   Tensor1 = this->Tensors[this->N[0]*i];
   Tensor1.complexConjugate();
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPS0.set(i, Tensor0);
  }
// the right end of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = 1; Shape0[2] = dMPS0;
  Tensor0 = PEPS0.Tensors[this->N[0]*(this->N[1]-1)];
  Tensor1 = this->Tensors[this->N[0]*(this->N[1]-1)];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(NMPS0-1, Tensor0);
 }
 else if (Direction == "left")
 {
  NMPS0 = this->N[0];
  MPS0 = MPS<T>(BCMPS0, NMPS0, dMPS0, DMPS0);
  Order0[0] = 1; Order0[1] = 5; Order0[2] = 2; Order0[3] = 6;
  Order0[4] = 3; Order0[5] = 7; Order0[6] = 0; Order0[7] = 4;
// the left end of MPS0:
  Shape0[0] = 1; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  Tensor0 = PEPS0.Tensors[this->N[0]*(this->N[1]-1)];
  Tensor1 = this->Tensors[this->N[0]*(this->N[1]-1)];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(0, Tensor0);
// the intermediate part of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  for (int i = 1; i < NMPS0-1; i++)
  {
   Tensor0 = PEPS0.Tensors[i+this->N[0]*(this->N[1]-1)];
   Tensor1 = this->Tensors[i+this->N[0]*(this->N[1]-1)];
   Tensor1.complexConjugate();
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPS0.set(i, Tensor0);
  }
// the right end of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = 1; Shape0[2] = dMPS0;
  Tensor0 = PEPS0.Tensors[this->N[0]*this->N[1]-1];
  Tensor1 = this->Tensors[this->N[0]*this->N[1]-1];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(NMPS0-1, Tensor0);
 }
 else if (Direction == "up")
 {
  NMPS0 = this->N[1];
  MPS0 = MPS<T>(BCMPS0, NMPS0, dMPS0, DMPS0);
  Order0[0] = 2; Order0[1] = 6; Order0[2] = 3; Order0[3] = 7;
  Order0[4] = 0; Order0[5] = 4; Order0[6] = 1; Order0[7] = 5;
// the left end of MPS0:
  Shape0[0] = 1; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  Tensor0 = PEPS0.Tensors[this->N[0]*this->N[1]-1];
  Tensor1 = this->Tensors[this->N[0]*this->N[1]-1];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(0, Tensor0);
// the intermediate part of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = DMPS0; Shape0[2] = dMPS0;
  for (int i = 1; i < NMPS0-1; i++)
  {
   Tensor0 = PEPS0.Tensors[this->N[0]*(this->N[1]-i)-1];
   Tensor1 = this->Tensors[this->N[0]*(this->N[1]-i)-1];
   Tensor1.complexConjugate();
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPS0.set(i, Tensor0);
  }
// the right end of MPS0:
  Shape0[0] = DMPS0; Shape0[1] = 1; Shape0[2] = dMPS0;
  Tensor0 = PEPS0.Tensors[this->N[0]-1];
  Tensor1 = this->Tensors[this->N[0]-1];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPS0.set(NMPS0-1, Tensor0);
 }
}

template<class T> void PEPS<T>::getBoundaryMPS(const PEPO<T>& PEPO0, const string& Boundary,
                                               double eps, unsigned int maxNumSweeps,
                                               double& epsAchieved, unsigned int& numSweepsDone,
                                               MPS<T>& MPS0) const
{
#ifdef DEBUG
 if ((this->BC != "open") || (this->N[0] == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getBoundaryMPS(const PEPO<T>& PEPO0, const string& Boundary, double eps, unsigned int maxNumSweeps, " <<
                         "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "((this->BC != open) || (this->N[0] == 0))." << endl;
  exit(1);
 }
 string BCPEPO0; PEPO0.getBC(BCPEPO0);
 vector<unsigned int> NPEPO0(2); NPEPO0[0] = PEPO0.getNrows(); NPEPO0[1] = PEPO0.getNcols();
 Matrix<unsigned int> dPEPO0; PEPO0.getd(dPEPO0);
 if ((this->BC != BCPEPO0) || (this->N != NPEPO0) || (this->d != dPEPO0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getBoundaryMPS(const PEPO<T>& PEPO0, const string& Boundary, double eps, unsigned int maxNumSweeps, " <<
                         "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "((this->BC != PEPO0.BC) || (this->N != PEPO0.N) || (this->d != PEPO0.d))." << endl;
  exit(1);
 }
 if ((Boundary != "left") && (Boundary != "top") && (Boundary != "right") && (Boundary != "bottom"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getBoundaryMPS(const PEPO<T>& PEPO0, const string& Boundary, double eps, unsigned int maxNumSweeps, " <<
                         "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "((Boundary != left) && (Boundary != top) && (Boundary != right) && (Boundary != bottom))." << endl;
  exit(1);
 }
 string BCMPS0; MPS0.getBC(BCMPS0);
 if ((BCMPS0 != "open") ||
     (((Boundary == "left") || (Boundary == "right")) && (MPS0.getN() != this->N[0])) ||
     (((Boundary == "top") || (Boundary == "bottom")) && (MPS0.getN() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getBoundaryMPS(const PEPO<T>& PEPO0, const string& Boundary, double eps, unsigned int maxNumSweeps, " <<
                         "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "((MPS0.BC != open) || (((Boundary == left) || (Boundary == right)) && (MPS0.N != this->N[0])) || " <<
           "(((Boundary == top) || (Boundary == bottom)) && (MPS0.N != this->N[1])))." << endl;
  exit(1);
 }
#endif
 unsigned int DMPS0, L;
 double error0, error1;
 T bValue, NValue, tensorFrobeniusNorm, tensorFrobeniusNormProduct;
 string SweepDirection;
 vector<unsigned int> dMPS0, Indices0(1), Indices02(2), Indices03(3), Indices04(4), Indices1(1), Indices12(2), Indices13(3), Indices14(4), Indices2(2);
 vector<unsigned int> Order3(3), Order5(5), Order6(6), Order9(9), Shape3(3), Shape4(4), Shape5(5), Shape6(6), Shape7(7), Shape8(8), Shape9(9);
 Tensor<T> bTensor, Tensor0, Tensor1, TensorL, TensorR;
 vector< Tensor<T> > BoundaryTensors, BoundaryTensorsPEPO, bTensors;
// 0. initialize:
// 0. - set boundary-length L:
 if ((Boundary == "left") || (Boundary == "right"))
  L = this->N[0];
 else if ((Boundary == "top") || (Boundary == "bottom"))
  L = this->N[1];
// 0. - get boundary-tensors clockwise, permute them as "top" boundary, and write the ones from this PEPS in BoundaryTensors,
//      and the ones from PEPO0 in BoundaryTensorsPEPO:
 BoundaryTensors = vector< Tensor<T> >(L);
 BoundaryTensorsPEPO = vector< Tensor<T> >(L);
 if (Boundary == "left")
 {
  Order5[0] = 3; Order5[1] = 0; Order5[2] = 1; Order5[3] = 2; Order5[4] = 4;
  Order6[0] = 3; Order6[1] = 0; Order6[2] = 1; Order6[3] = 2; Order6[4] = 4; Order6[5] = 5;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[this->N[0]-1-pos];
   Tensor0.permute(Order5);
   BoundaryTensors[pos] = Tensor0;
   PEPO0.get(this->N[0]-1-pos, 0, Tensor0);
   Tensor0.permute(Order6);
   BoundaryTensorsPEPO[pos] = Tensor0;
  }
 }
 else if (Boundary == "top")
 {
  for (int pos = 0; pos < L; pos++)
  {
   BoundaryTensors[pos] = this->Tensors[this->N[0]*pos];
   PEPO0.get(0, pos, BoundaryTensorsPEPO[pos]);
  }
 }
 else if (Boundary == "right")
 {
  Order5[0] = 1; Order5[1] = 2; Order5[2] = 3; Order5[3] = 0; Order5[4] = 4;
  Order6[0] = 1; Order6[1] = 2; Order6[2] = 3; Order6[3] = 0; Order6[4] = 4; Order6[5] = 5;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[pos+this->N[0]*(this->N[1]-1)];
   Tensor0.permute(Order5);
   BoundaryTensors[pos] = Tensor0;
   PEPO0.get(pos, this->N[1]-1, Tensor0);
   Tensor0.permute(Order6);
   BoundaryTensorsPEPO[pos] = Tensor0;
  }
 }
 else if (Boundary == "bottom")
 {
  Order5[0] = 2; Order5[1] = 3; Order5[2] = 0; Order5[3] = 1; Order5[4] = 4;
  Order6[0] = 2; Order6[1] = 3; Order6[2] = 0; Order6[3] = 1; Order6[4] = 4; Order6[5] = 5;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[this->N[0]-1+this->N[0]*(this->N[1]-1-pos)];
   Tensor0.permute(Order5);
   BoundaryTensors[pos] = Tensor0;
   PEPO0.get(this->N[0]-1, this->N[1]-1-pos, Tensor0);
   Tensor0.permute(Order6);
   BoundaryTensorsPEPO[pos] = Tensor0;
  }
 }
// 0. - reshape boundary-tensors eliminating index 1:
 for (int pos = 0; pos < L; pos++)
 {
  BoundaryTensors[pos].getShape(Shape5);
  Shape4[0] = Shape5[0]; Shape4[1] = Shape5[2]; Shape4[2] = Shape5[3]; Shape4[3] = Shape5[4];
  BoundaryTensors[pos].reshape(Shape4);
  BoundaryTensorsPEPO[pos].getShape(Shape6);
  Shape5[0] = Shape6[0]; Shape5[1] = Shape6[2]; Shape5[2] = Shape6[3]; Shape5[3] = Shape6[4]; Shape5[4] = Shape6[5];
  BoundaryTensorsPEPO[pos].reshape(Shape5);
 }
// 0. - check MPS0.d:
 MPS0.getd(dMPS0);
#ifdef DEBUG
 for (int pos = 0; pos < L; pos++)
 {
  BoundaryTensors[pos].getShape(Shape4);
  BoundaryTensorsPEPO[pos].getShape(Shape5);
  if (dMPS0[pos] != Shape4[2]*Shape5[2]*Shape4[2])
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void PEPS<T>::" <<
           "getBoundaryMPS(const PEPO<T>& PEPO0, const string& Boundary, double eps, unsigned int maxNumSweeps, " <<
                          "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
           "(MPS0.d[" << pos << "] != this->D[" << pos << "]*PEPO0.D[" << pos << "]*this->D[" << pos << "])." << endl;
   exit(1);
  }
 }
#endif
// 1. MPS0.D >= this->D*PEPO0.D*this->D:
 DMPS0 = MPS0.getD();
 if (DMPS0 >= this->D*PEPO0.getD()*this->D)
 {
  MPS0 = MPS<T>("open", L, dMPS0, this->D*PEPO0.getD()*this->D);
  Indices1[0] = 3;
  Order9[0] = 0; Order9[1] = 3; Order9[2] = 6; Order9[3] = 1; Order9[4] = 4; Order9[5] = 7; Order9[6] = 2; Order9[7] = 5; Order9[8] = 8;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = BoundaryTensors[pos];
   Tensor1 = BoundaryTensorsPEPO[pos];
   Indices0[0] = 3;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   BoundaryTensors[pos].complexConjugate(Tensor1);
   Indices0[0] = 6;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order9);
   Tensor0.getShape(Shape9);
   Shape3[0] = Shape9[0]*Shape9[1]*Shape9[2]; Shape3[1] = Shape9[3]*Shape9[4]*Shape9[5]; Shape3[2] = Shape9[6]*Shape9[7]*Shape9[8];
   Tensor0.reshape(Shape3);
   MPS0.set(pos, Tensor0);
  }
  MPS0.bringIntoNormalShape();
  MPS0.setD(DMPS0);
  epsAchieved = 0.0;
  numSweepsDone = 0;
  return;
 }
// 2. bring MPS0 into normal form from right to left:
// 2.a) MPS0.D == 1:
 if (DMPS0 == 1)
 {
  tensorFrobeniusNormProduct = 1.0;
  for (int pos = L-1; pos > 0; pos--)
  {
   MPS0.get(pos, Tensor0);
   tensorFrobeniusNormProduct *= Tensor0.frobeniusNormalize();
   MPS0.set(pos, Tensor0);
  }
  MPS0.get(0, Tensor0);
  Tensor0.multiply(tensorFrobeniusNormProduct);
  MPS0.set(0, Tensor0);
 }
// 2.b) MPS0.D != 1:
 else
 {
  Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
  Indices0[0] = 1;
  Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
  MPS0.get(L-1, Tensor0);
  for (int pos = L-1; pos > 0; pos--)
  {
   Tensor0.LQDecompose(Indices1, Indices2, TensorL);
   MPS0.set(pos, Tensor0);
   MPS0.get(pos-1, Tensor0);
   Tensor0.contract(Indices0, TensorL, Indices1);
   Tensor0.permute(Order3);
  }
  MPS0.set(0, Tensor0);
 }
// 3. construct right half of b for each site:
 bTensors = vector< Tensor<T> >(L);
 Tensor0 = BoundaryTensors[L-1];
 Tensor1 = BoundaryTensorsPEPO[L-1];
 Indices0[0] = 3; Indices1[0] = 3;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 BoundaryTensors[L-1].complexConjugate(Tensor1);
 Indices0[0] = 6;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor0.getShape(Shape9);
 MPS0.get(L-1, Tensor1);
 Tensor1.getShape(Shape3);
 Shape5[0] = Shape3[0]; Shape5[1] = 1; Shape5[2] = Shape9[2]; Shape5[3] = Shape9[5]; Shape5[4] = Shape9[8];
 Tensor1.reshape(Shape5);
 Tensor1.complexConjugate();
 Indices03[0] = 2; Indices03[1] = 5; Indices03[2] = 8; Indices13[0] = 2; Indices13[1] = 3; Indices13[2] = 4;
 Tensor0.contract(Indices03, Tensor1, Indices13);
 Tensor0.getShape(Shape8);
 Shape4[0] = Shape8[0]; Shape4[1] = Shape8[2]; Shape4[2] = Shape8[4]; Shape4[3] = Shape8[6];
 Tensor0.reshape(Shape4);
 bTensors[L-1] = Tensor0;
 Indices0[0] = 0; Indices1[0] = 1;
 Indices02[0] = 0; Indices12[0] = 1; Indices12[1] = 3;
 Indices04[0] = 0; Indices04[1] = 2; Indices04[2] = 4; Indices04[3] = 6; Indices14[0] = 1; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
 for (int pos = L-2; pos > 0; pos--)
 {
  Tensor1 = BoundaryTensors[pos];
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor1 = BoundaryTensorsPEPO[pos];
  Indices02[1] = 5;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  BoundaryTensors[pos].complexConjugate(Tensor1);
  Indices02[1] = 6;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape7);
  MPS0.get(pos, Tensor1);
  Tensor1.getShape(Shape3);
  Shape5[0] = Shape3[0]; Shape5[1] = Shape3[1]; Shape5[2] = Shape7[2]; Shape5[3] = Shape7[4]; Shape5[4] = Shape7[6];
  Tensor1.reshape(Shape5);
  Tensor1.complexConjugate();
  Tensor0.contract(Indices04, Tensor1, Indices14);
  bTensors[pos] = Tensor0;
 }
// 4. compute initial error:
 Tensor1 = BoundaryTensors[0];
 Indices0[0] = 0; Indices1[0] = 1;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = BoundaryTensorsPEPO[0];
 Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 1; Indices12[1] = 3;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 BoundaryTensors[0].complexConjugate(Tensor1);
 Indices02[1] = 6;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 Tensor0.getShape(Shape7);
 MPS0.get(0, Tensor1);
 Tensor1.getShape(Shape3);
 Shape5[0] = 1; Shape5[1] = Shape3[1]; Shape5[2] = Shape7[2]; Shape5[3] = Shape7[4]; Shape5[4] = Shape7[6];
 Tensor1.reshape(Shape5);
 Tensor1.complexConjugate();
 Indices04[0] = 0; Indices04[1] = 2; Indices04[2] = 4; Indices04[3] = 6; Indices14[0] = 1; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
 Tensor0.contract(Indices04, Tensor1, Indices14);
 bValue = Tensor0.get(0);
 MPS0.get(0, Tensor0);
 NValue = Tensor0.scalarProduct();
 error0 = MathAuxiliary::convertToDouble(NValue)-2.0*MathAuxiliary::convertToDouble(bValue);
// 5. sweep:
 numSweepsDone = 0;
 SweepDirection = "right";
 while (numSweepsDone < maxNumSweeps)
 {
// 5.a) sweep to the right:
  if (SweepDirection == "right")
  {
// 5.a)1. position 0:
// 5.a)1. - compute bTensor:
   Tensor0 = bTensors[1];
   Tensor1 = BoundaryTensors[0];
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BoundaryTensorsPEPO[0];
   Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 1; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   BoundaryTensors[0].complexConjugate(Tensor1);
   Indices02[1] = 6;
   Tensor0.contract(Indices02, Tensor1, Indices12);
// 5.a)1. - update tensor:
   Tensor0.getShape(Shape7);
   Shape3[0] = 1; Shape3[1] = Shape7[0]; Shape3[2] = Shape7[2]*Shape7[4]*Shape7[6];
   Tensor0.reshape(Shape3);
   MPS0.set(0, Tensor0);
// 5.a)1. - compute error:
   NValue = Tensor0.scalarProduct();
   error1 = -MathAuxiliary::convertToDouble(NValue);
   epsAchieved = abs((error1-error0)/error1);
   if (epsAchieved < eps)
    return;
   error0 = error1;
// 5.a)1. - normalize tensor:
   if (DMPS0 == 1)
    tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
   else
   {
    Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
    Tensor0.QRDecompose(Indices2, Indices1, TensorR);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
   }
   MPS0.set(0, Tensor0);
// 5.a)1. - set bTensors:
   Tensor0 = BoundaryTensors[0];
   Tensor1 = BoundaryTensorsPEPO[0];
   Indices0[0] = 3; Indices1[0] = 3;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   BoundaryTensors[0].complexConjugate(Tensor1);
   Indices0[0] = 6;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.getShape(Shape9);
   MPS0.get(0, Tensor1);
   Tensor1.getShape(Shape3);
   Shape5[0] = 1; Shape5[1] = Shape3[1]; Shape5[2] = Shape9[2]; Shape5[3] = Shape9[5]; Shape5[4] = Shape9[8];
   Tensor1.reshape(Shape5);
   Tensor1.complexConjugate();
   Indices03[0] = 2; Indices03[1] = 5; Indices03[2] = 8; Indices13[0] = 2; Indices13[1] = 3; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape8);
   Shape4[0] = Shape8[1]; Shape4[1] = Shape8[3]; Shape4[2] = Shape8[5]; Shape4[3] = Shape8[7];
   Tensor0.reshape(Shape4);
   bTensors[0] = Tensor0;
// 5.a)2. 0 < position < L-1:
   for (int pos = 1; pos < L-1; pos++)
   {
// 5.a)2. - compute bTensor:
    Tensor0 = bTensors[pos-1];
    Tensor1 = BoundaryTensors[pos];
    Indices0[0] = 0; Indices1[0] = 0;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = BoundaryTensorsPEPO[pos];
    Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 3;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    BoundaryTensors[pos].complexConjugate(Tensor1);
    Indices02[1] = 6;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    bTensor = Tensor0;
// 5.a)2. - update tensor:
    Tensor1 = bTensors[pos+1];
    Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 5; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 2;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    Tensor0.getShape(Shape5);
    Shape3[0] = Shape5[0]; Shape3[1] = Shape5[1]*Shape5[2]*Shape5[3]; Shape3[2] = Shape5[4];
    Tensor0.reshape(Shape3);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
    MPS0.set(pos, Tensor0);
// 5.a)2. - compute error:
    NValue = Tensor0.scalarProduct();
    error1 = -MathAuxiliary::convertToDouble(NValue);
    epsAchieved = abs((error1-error0)/error1);
    if (epsAchieved < eps)
     return;
    error0 = error1;
// 5.a)2. - normalize tensor:
    if (DMPS0 == 1)
     tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
    else
    {
     Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
     Tensor0.QRDecompose(Indices2, Indices1, TensorR);
     Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
     Tensor0.permute(Order3);
    }
    MPS0.set(pos, Tensor0);
// 5.a)2. - set bTensors:
    bTensor.getShape(Shape7);
    Tensor0.getShape(Shape3);
    Shape5[0] = Shape3[0]; Shape5[1] = Shape3[1]; Shape5[2] = Shape7[2]; Shape5[3] = Shape7[4]; Shape5[4] = Shape7[6];
    Tensor0.reshape(Shape5);
    Tensor0.complexConjugate();
    Indices04[0] = 0; Indices04[1] = 2; Indices04[2] = 4; Indices04[3] = 6; Indices14[0] = 0; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
    bTensor.contract(Indices04, Tensor0, Indices14);
    bTensors[pos] = bTensor;
   }
// 5.a)3. position L-1:
   MPS0.get(L-1, Tensor0);
   if (DMPS0 == 1)
    Tensor0.multiply(tensorFrobeniusNorm);
   else
   {
    Indices0[0] = 1; Indices1[0] = 0;
    TensorR.contract(Indices0, Tensor0, Indices1);
    Tensor0 = TensorR;
   }
   MPS0.set(L-1, Tensor0);
   SweepDirection = "left";
  }
// 5.b) sweep to the left:
  else if (SweepDirection == "left")
  {
// 5.b)1. position L-1:
// 5.b)1. - compute bTensor:
   Tensor0 = bTensors[L-2];
   Tensor1 = BoundaryTensors[L-1];
   Indices0[0] = 0; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BoundaryTensorsPEPO[L-1];
   Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   BoundaryTensors[L-1].complexConjugate(Tensor1);
   Indices02[1] = 6;
   Tensor0.contract(Indices02, Tensor1, Indices12);
// 5.b)1. - update tensor:
   Tensor0.getShape(Shape7);
   Shape3[0] = Shape7[0]; Shape3[1] = 1; Shape3[2] = Shape7[2]*Shape7[4]*Shape7[6];
   Tensor0.reshape(Shape3);
   MPS0.set(L-1, Tensor0);
// 5.b)1. - compute error:
   NValue = Tensor0.scalarProduct();
   error1 = -MathAuxiliary::convertToDouble(NValue);
   epsAchieved = abs((error1-error0)/error1);
   if (epsAchieved < eps)
    return;
   error0 = error1;
// 5.b)1. - normalize tensor:
   if (DMPS0 == 1)
    tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
   else
   {
    Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
    Tensor0.LQDecompose(Indices1, Indices2, TensorL);
   }
   MPS0.set(L-1, Tensor0);
// 5.b)1. - set bTensors:
   Tensor0 = BoundaryTensors[L-1];
   Tensor1 = BoundaryTensorsPEPO[L-1];
   Indices0[0] = 3; Indices1[0] = 3;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   BoundaryTensors[L-1].complexConjugate(Tensor1);
   Indices0[0] = 6;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.getShape(Shape9);
   MPS0.get(L-1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape5[0] = Shape3[0]; Shape5[1] = 1; Shape5[2] = Shape9[2]; Shape5[3] = Shape9[5]; Shape5[4] = Shape9[8];
   Tensor1.reshape(Shape5);
   Tensor1.complexConjugate();
   Indices03[0] = 2; Indices03[1] = 5; Indices03[2] = 8; Indices13[0] = 2; Indices13[1] = 3; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape8);
   Shape4[0] = Shape8[0]; Shape4[1] = Shape8[2]; Shape4[2] = Shape8[4]; Shape4[3] = Shape8[6];
   Tensor0.reshape(Shape4);
   bTensors[L-1] = Tensor0;
// 5.b)2. 0 < position < L-1:
   for (int pos = L-2; pos > 0; pos--)
   {
// 5.b)2. - compute bTensor:
    Tensor0 = bTensors[pos+1];
    Tensor1 = BoundaryTensors[pos];
    Indices0[0] = 0; Indices1[0] = 1;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = BoundaryTensorsPEPO[pos];
    Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 1; Indices12[1] = 3;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    BoundaryTensors[pos].complexConjugate(Tensor1);
    Indices02[1] = 6;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    bTensor = Tensor0;
// 5.b)2. - update tensor:
    Tensor1 = bTensors[pos-1];
    Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 5; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 2;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    Tensor0.getShape(Shape5);
    Shape3[0] = Shape5[0]; Shape3[1] = Shape5[1]*Shape5[2]*Shape5[3]; Shape3[2] = Shape5[4];
    Tensor0.reshape(Shape3);
    Order3[0] = 2; Order3[1] = 0; Order3[2] = 1;
    Tensor0.permute(Order3);
    MPS0.set(pos, Tensor0);
// 5.b)2. - compute error:
    NValue = Tensor0.scalarProduct();
    error1 = -MathAuxiliary::convertToDouble(NValue);
    epsAchieved = abs((error1-error0)/error1);
    if (epsAchieved < eps)
     return;
    error0 = error1;
// 5.b)2. - normalize tensor:
    if (DMPS0 == 1)
     tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
    else
    {
     Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
     Tensor0.LQDecompose(Indices1, Indices2, TensorL);
    }
    MPS0.set(pos, Tensor0);
// 5.b)2. - set bTensors:
    bTensor.getShape(Shape7);
    Tensor0.getShape(Shape3);
    Shape5[0] = Shape3[0]; Shape5[1] = Shape3[1]; Shape5[2] = Shape7[2]; Shape5[3] = Shape7[4]; Shape5[4] = Shape7[6];
    Tensor0.reshape(Shape5);
    Tensor0.complexConjugate();
    Indices04[0] = 0; Indices04[1] = 2; Indices04[2] = 4; Indices04[3] = 6; Indices14[0] = 1; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
    bTensor.contract(Indices04, Tensor0, Indices14);
    bTensors[pos] = bTensor;
   }
// 5.b)3. position 0:
   MPS0.get(0, Tensor0);
   if (DMPS0 == 1)
    Tensor0.multiply(tensorFrobeniusNorm);
   else
   {
    Indices0[0] = 1; Indices1[0] = 0;
    Tensor0.contract(Indices0, TensorL, Indices1);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
   }
   MPS0.set(0, Tensor0);
   SweepDirection = "right";
  }
  numSweepsDone++;
 } 
}

template<class T> void PEPS<T>::getMPO(const PEPO<T>& PEPO0, const string& Direction,
                                       unsigned int position, MPO<T>& MPO0) const
{
#ifdef DEBUG
 string BCPEPO0; PEPO0.getBC(BCPEPO0);
 vector<unsigned int> NPEPO0(2); NPEPO0[0] = PEPO0.getNrows(); NPEPO0[1] = PEPO0.getNcols();
 Matrix<unsigned int> dPEPO0; PEPO0.getd(dPEPO0);
 if ((this->N[0] == 0) || (this->BC != BCPEPO0) || (this->N != NPEPO0) || (this->d != dPEPO0) ||
     ((Direction != "right") && (Direction != "down") && (Direction != "left") &&
      (Direction != "up")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getMPO(const PEPO<T>& PEPO0, const string& Direction, " <<
                 "unsigned int position, MPO<T>& MPO0) const: " <<
          "((this->N[0] == 0) || (this->BC != BCPEPO0) || (this->N != NPEPO0) || " <<
           "(this->d != dPEPO0) || " <<
           "((Direction != right) && (Direction != down) && (Direction != left) && " <<
            "(Direction != up)))." << endl;
  exit(1);
 }
 else if (((position > this->N[1]-1) && ((Direction == "right") || (Direction == "left"))) ||
          ((position > this->N[0]-1) && ((Direction == "down") || (Direction == "up"))))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getMPO(const PEPO<T>& PEPO0, const string& Direction, " <<
                 "unsigned int position, MPO<T>& MPO0) const: " <<
          "(((position > this->N[1]-1) && ((Direction == right) || (Direction == left))) || " <<
           "((position > this->N[0]-1) && ((Direction == down) || (Direction == up))))." << endl;
  exit(1);
 }
#endif
// compute |PEPS0>=PEPO0*|this>:
 unsigned D0 = this->D*PEPO0.getD();
 PEPS<T> PEPS0(this->BC, this->N[0], this->N[1], this->d, D0);
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Indices0(1), Indices1(1);
 Indices0[0] = 4; Indices1[0] = 4;
 vector<unsigned int> Order(9), Shape(5);
 Order[0] = 0; Order[1] = 4; Order[2] = 1; Order[3] = 5; Order[4] = 2; Order[5] = 6; Order[6] = 3; Order[7] = 7; Order[8] = 8;
 Shape[0] = D0; Shape[1] = D0; Shape[2] = D0; Shape[3] = D0;
 for (int col = 0; col < this->N[1]; col++)
 {
  for (int row = 0; row < this->N[0]; row++)
  {
   this->get(row, col, Tensor0); PEPO0.get(row, col, Tensor1);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order);
   Shape[4] = this->d(row, col);
   if (this->BC == "open")
    PEPS0.getOpenBCShape(row, col, Shape);
   Tensor0.reshape(Shape);
   PEPS0.set(row, col, Tensor0);
  }
 }
// compute MPO0:
 string BCMPO0 = this->BC;
 unsigned int NMPO0;
 unsigned int dMPO0 = D0*this->D, DMPO0 = D0*this->D;
 vector<unsigned int> Shape0(4), Order0(8);
 T element;
 if (Direction == "right")
 {
  NMPO0 = this->N[0];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order0[0] = 3; Order0[1] = 7; Order0[2] = 1; Order0[3] = 5;
  Order0[4] = 0; Order0[5] = 4; Order0[6] = 2; Order0[7] = 6;
// the left end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = 1; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = PEPS0.Tensors[this->N[0]*(position+1)-1];
  Tensor1 = this->Tensors[this->N[0]*(position+1)-1];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(0, Tensor0);
// the intermediate part of MPO0:
  Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  for (int i = 1; i < NMPO0-1; i++)
  {
   Tensor0 = PEPS0.Tensors[this->N[0]*(position+1)-1-i];
   Tensor1 = this->Tensors[this->N[0]*(position+1)-1-i];
   Tensor1.complexConjugate();
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPO0.set(i, Tensor0);
  }
// the right end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = DMPO0; Shape0[1] = 1; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = PEPS0.Tensors[this->N[0]*position];
  Tensor1 = this->Tensors[this->N[0]*position];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(NMPO0-1, Tensor0);
 }
 else if (Direction == "down")
 {
  NMPO0 = this->N[1];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order0[0] = 0; Order0[1] = 4; Order0[2] = 2; Order0[3] = 6;
  Order0[4] = 1; Order0[5] = 5; Order0[6] = 3; Order0[7] = 7;
// the left end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = 1; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = PEPS0.Tensors[position];
  Tensor1 = this->Tensors[position];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(0, Tensor0);
// the intermediate part of MPO0:
  Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  for (int i = 1; i < NMPO0-1; i++)
  {
   Tensor0 = PEPS0.Tensors[position+this->N[0]*i];
   Tensor1 = this->Tensors[position+this->N[0]*i];
   Tensor1.complexConjugate();
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPO0.set(i, Tensor0);
  }
// the right end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = DMPO0; Shape0[1] = 1; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = PEPS0.Tensors[position+this->N[0]*(this->N[1]-1)];
  Tensor1 = this->Tensors[position+this->N[0]*(this->N[1]-1)];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(NMPO0-1, Tensor0);
 }
 else if (Direction == "left")
 {
  NMPO0 = this->N[0];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order0[0] = 1; Order0[1] = 5; Order0[2] = 3; Order0[3] = 7;
  Order0[4] = 2; Order0[5] = 6; Order0[6] = 0; Order0[7] = 4;
// the left end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = 1; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = PEPS0.Tensors[this->N[0]*position];
  Tensor1 = this->Tensors[this->N[0]*position];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(0, Tensor0);
// the intermediate part of MPO0:
  Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  for (int i = 1; i < NMPO0-1; i++)
  {
   Tensor0 = PEPS0.Tensors[i+this->N[0]*position];
   Tensor1 = this->Tensors[i+this->N[0]*position];
   Tensor1.complexConjugate();
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPO0.set(i, Tensor0);
  }
// the right end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = DMPO0; Shape0[1] = 1; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = PEPS0.Tensors[this->N[0]*(position+1)-1];
  Tensor1 = this->Tensors[this->N[0]*(position+1)-1];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(NMPO0-1, Tensor0);
 }
 else if (Direction == "up")
 {
  NMPO0 = this->N[1];
  MPO0 = MPO<T>(BCMPO0, NMPO0, dMPO0, DMPO0);
  Order0[0] = 2; Order0[1] = 6; Order0[2] = 0; Order0[3] = 4;
  Order0[4] = 3; Order0[5] = 7; Order0[6] = 1; Order0[7] = 5;
// the left end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = 1; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = PEPS0.Tensors[position+this->N[0]*(this->N[1]-1)];
  Tensor1 = this->Tensors[position+this->N[0]*(this->N[1]-1)];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(0, Tensor0);
// the intermediate part of MPO0:
  Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  for (int i = 1; i < NMPO0-1; i++)
  {
   Tensor0 = PEPS0.Tensors[position+this->N[0]*(this->N[1]-1-i)];
   Tensor1 = this->Tensors[position+this->N[0]*(this->N[1]-1-i)];
   Tensor1.complexConjugate();
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order0);
   Tensor0.reshape(Shape0);
   MPO0.set(i, Tensor0);
  }
// the right end of MPO0:
  if (BCMPO0 == "open")
  {
   Shape0[0] = DMPO0; Shape0[1] = 1; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  else if (BCMPO0 == "periodic")
  {
   Shape0[0] = DMPO0; Shape0[1] = DMPO0; Shape0[2] = dMPO0; Shape0[3] = dMPO0;
  }
  Tensor0 = PEPS0.Tensors[position];
  Tensor1 = this->Tensors[position];
  Tensor1.complexConjugate();
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.permute(Order0);
  Tensor0.reshape(Shape0);
  MPO0.set(NMPO0-1, Tensor0);
 }
}

template<class T> void PEPS<T>::multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const string& Direction, unsigned int position,
                                                   const MPS<T>& MPS0,
                                                   double eps, unsigned int maxNumSweeps,
                                                   double& epsAchieved, unsigned int& numSweepsDone,
                                                   MPS<T>& MPS1) const
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const string& Direction, unsigned int position, " <<
                             "const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "(this->N[0] == 0)." << endl;
  exit(1);
 }
 string BCPEPO0; PEPO0.getBC(BCPEPO0);
 vector<unsigned int> NPEPO0(2); NPEPO0[0] = PEPO0.getNrows(); NPEPO0[1] = PEPO0.getNcols();
 Matrix<unsigned int> dPEPO0; PEPO0.getd(dPEPO0);
 if ((this->BC != BCPEPO0) || (this->N != NPEPO0) || (this->d != dPEPO0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const string& Direction, unsigned int position, " <<
                             "const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((this->BC != PEPO0.BC) || (this->N != PEPO0.N) || (this->d != PEPO0.d))." << endl;
  exit(1);
 }
 if ((Direction != "right") && (Direction != "down") && (Direction != "left") && (Direction != "up"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const string& Direction, unsigned int position, " <<
                             "const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((Direction != right) && (Direction != down) && (Direction != left) && (Direction != up))." << endl;
  exit(1);
 }
 if ((((Direction == "right") || (Direction == "left")) && (position > this->N[1]-1)) ||
     (((Direction == "down") || (Direction == "up")) && (position > this->N[0]-1)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const string& Direction, unsigned int position, " <<
                             "const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((((Direction == right) || (Direction == left)) && (position > this->N[1]-1)) || " <<
           "(((Direction == down) || (Direction == up)) && (position > this->N[0]-1)))." << endl;
  exit(1);
 }
 string BCMPS;
 MPS0.getBC(BCMPS);
 if ((BCMPS != "open") ||
     (((Direction == "right") || (Direction == "left")) && (MPS0.getN() != this->N[0])) ||
     (((Direction == "down") || (Direction == "up")) && (MPS0.getN() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const string& Direction, unsigned int position, " <<
                             "const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((MPS0.BC != open) || (((Direction == right) || (Direction == left)) && (MPS0.N != this->N[0])) || " <<
           "(((Direction == down) || (Direction == up)) && (MPS0.N != this->N[1])))." << endl;
  exit(1);
 }
 MPS1.getBC(BCMPS);
 if ((BCMPS != "open") ||
     (((Direction == "right") || (Direction == "left")) && (MPS1.getN() != this->N[0])) ||
     (((Direction == "down") || (Direction == "up")) && (MPS1.getN() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const string& Direction, unsigned int position, " <<
                             "const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((MPS1.BC != open) || (((Direction == right) || (Direction == left)) && (MPS1.N != this->N[0])) || " <<
           "(((Direction == down) || (Direction == up)) && (MPS1.N != this->N[1])))." << endl;
  exit(1);
 }
#endif
 unsigned int DMPS1, L;
 double error0, error1;
 T bValue, NValue, tensorFrobeniusNorm, tensorFrobeniusNormProduct;
 string SweepDirection;
 vector<unsigned int> dMPS, Indices0(1), Indices02(2), Indices03(3), Indices04(4), Indices1(1), Indices12(2), Indices13(3), Indices14(4), Indices2(2);
 vector<unsigned int> Order3(3), Order5(5), Order6(6), Order11(11);
 vector<unsigned int> Shape3(3), Shape5(5), Shape15(5), Shape6(6), Shape8(8), Shape10(10), Shape11(11);
 Tensor<T> bTensor, Tensor0, Tensor1, TensorL, TensorR;
 vector< Tensor<T> > BulkTensors, BulkTensorsPEPO, bTensors;
// 0. initialize:
// 0. - set length L:
 if ((Direction == "right") || (Direction == "left"))
  L = this->N[0];
 else if ((Direction == "down") || (Direction == "up"))
  L = this->N[1];
// 0. - get bulk-tensors clockwise, permute them as Direction "down", and write the ones from this PEPS in BulkTensors, and
//      the ones from PEPO0 in BulkTensorsPEPO:
 BulkTensors = vector< Tensor<T> >(L);
 BulkTensorsPEPO = vector< Tensor<T> >(L);
 if (Direction == "right")
 {
  Order5[0] = 3; Order5[1] = 0; Order5[2] = 1; Order5[3] = 2; Order5[4] = 4;
  Order6[0] = 3; Order6[1] = 0; Order6[2] = 1; Order6[3] = 2; Order6[4] = 4; Order6[5] = 5;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[this->N[0]-1-pos+this->N[0]*position];
   Tensor0.permute(Order5);
   BulkTensors[pos] = Tensor0;
   PEPO0.get(this->N[0]-1-pos, position, Tensor0);
   Tensor0.permute(Order6);
   BulkTensorsPEPO[pos] = Tensor0;
  }
 }
 else if (Direction == "down")
 {
  for (int pos = 0; pos < L; pos++)
  {
   BulkTensors[pos] = this->Tensors[position+this->N[0]*pos];
   PEPO0.get(position, pos, BulkTensorsPEPO[pos]);
  }
 }
 else if (Direction == "left")
 {
  Order5[0] = 1; Order5[1] = 2; Order5[2] = 3; Order5[3] = 0; Order5[4] = 4;
  Order6[0] = 1; Order6[1] = 2; Order6[2] = 3; Order6[3] = 0; Order6[4] = 4; Order6[5] = 5;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[pos+this->N[0]*position];
   Tensor0.permute(Order5);
   BulkTensors[pos] = Tensor0;
   PEPO0.get(pos, position, Tensor0);
   Tensor0.permute(Order6);
   BulkTensorsPEPO[pos] = Tensor0;
  }
 }
 else if (Direction == "up")
 {
  Order5[0] = 2; Order5[1] = 3; Order5[2] = 0; Order5[3] = 1; Order5[4] = 4;
  Order6[0] = 2; Order6[1] = 3; Order6[2] = 0; Order6[3] = 1; Order6[4] = 4; Order6[5] = 5;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[position+this->N[0]*(this->N[1]-1-pos)];
   Tensor0.permute(Order5);
   BulkTensors[pos] = Tensor0;
   PEPO0.get(position, this->N[1]-1-pos, Tensor0);
   Tensor0.permute(Order6);
   BulkTensorsPEPO[pos] = Tensor0;
  }
 }
// 0. - check MPS0.d:
#ifdef DEBUG
 MPS0.getd(dMPS);
 for (int pos = 0; pos < L; pos++)
 {
  BulkTensors[pos].getShape(Shape5);
  BulkTensorsPEPO[pos].getShape(Shape6);
  if (dMPS[pos] != Shape5[1]*Shape6[1]*Shape5[1])
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void PEPS<T>::" <<
           "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const string& Direction, unsigned int position, " <<
                              "const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                              "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
           "(MPS0.d[" << pos << "] != this->D[" << pos << "]*PEPO0.D[" << pos << "]*this->D[" << pos << "])." << endl;
   exit(1);
  }
 }
#endif
// 0. - check MPS1.d:
 MPS1.getd(dMPS);
#ifdef DEBUG
 for (int pos = 0; pos < L; pos++)
 {
  BulkTensors[pos].getShape(Shape5);
  BulkTensorsPEPO[pos].getShape(Shape6);
  if (dMPS[pos] != Shape5[3]*Shape6[3]*Shape5[3])
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void PEPS<T>::" <<
           "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const string& Direction, unsigned int position, " <<
                              "const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                              "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
           "(MPS1.d[" << pos << "] != this->D[" << pos << "]*PEPO0.D[" << pos << "]*this->D[" << pos << "])." << endl;
   exit(1);
  }
 }
#endif
// 1. MPS1.D >= this->D*PEPO0.D*this->D*MPS0.D:
 DMPS1 = MPS1.getD();
 if (DMPS1 >= this->D*PEPO0.getD()*this->D*MPS0.getD())
 {
  MPS1 = MPS<T>("open", L, dMPS, this->D*PEPO0.getD()*this->D*MPS0.getD());
  Indices0[0] = 2; Indices1[0] = 1;
  Indices02[0] = 2; Indices12[0] = 1; Indices12[1] = 4;
  Order11[0] = 0; Order11[1] = 2; Order11[2] = 5; Order11[3] = 8; Order11[4] = 1; Order11[5] = 3;
  Order11[6] = 6; Order11[7] = 9; Order11[8] = 4; Order11[9] = 7; Order11[10] = 10;
  for (int pos = 0; pos < L; pos++)
  {
   MPS0.get(pos, Tensor0);
   Tensor0.getShape(Shape3);
   BulkTensors[pos].getShape(Shape5);
   BulkTensorsPEPO[pos].getShape(Shape6);
   Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape5[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
   Tensor0.reshape(Shape15);
   Tensor1 = BulkTensors[pos];
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BulkTensorsPEPO[pos];
   Indices02[1] = 7;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   BulkTensors[pos].complexConjugate(Tensor1);
   Indices02[1] = 9;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.permute(Order11);
   Tensor0.getShape(Shape11);
   Shape3[0] = Shape11[0]*Shape11[1]*Shape11[2]*Shape11[3]; Shape3[1] = Shape11[4]*Shape11[5]*Shape11[6]*Shape11[7];
   Shape3[2] = Shape11[8]*Shape11[9]*Shape11[10];
   Tensor0.reshape(Shape3);
   MPS1.set(pos, Tensor0);
  }
  MPS1.bringIntoNormalShape();
  MPS1.setD(DMPS1);
  epsAchieved = 0.0;
  numSweepsDone = 0;
  return;
 }
// 2. bring MPS1 into normal form from right to left:
// 2.a) MPS1.D == 1:
 if (DMPS1 == 1)
 {
  tensorFrobeniusNormProduct = 1.0;
  for (int pos = L-1; pos > 0; pos--)
  {
   MPS1.get(pos, Tensor0);
   tensorFrobeniusNormProduct *= Tensor0.frobeniusNormalize();
   MPS1.set(pos, Tensor0);
  }
  MPS1.get(0, Tensor0);
  Tensor0.multiply(tensorFrobeniusNormProduct);
  MPS1.set(0, Tensor0);
 }
// 2.b) MPS1.D != 1:
 else
 {
  Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
  Indices0[0] = 1;
  Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
  MPS1.get(L-1, Tensor0);
  for (int pos = L-1; pos > 0; pos--)
  {
   Tensor0.LQDecompose(Indices1, Indices2, TensorL);
   MPS1.set(pos, Tensor0);
   MPS1.get(pos-1, Tensor0);
   Tensor0.contract(Indices0, TensorL, Indices1);
   Tensor0.permute(Order3);
  }
  MPS1.set(0, Tensor0);
 }
// 3. construct right half of b for each site:
 bTensors = vector< Tensor<T> >(L);
 MPS0.get(L-1, Tensor0);
 Tensor0.getShape(Shape3);
 BulkTensors[L-1].getShape(Shape5);
 BulkTensorsPEPO[L-1].getShape(Shape6);
 Shape15[0] = Shape3[0]; Shape15[1] = 1; Shape15[2] = Shape5[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
 Tensor0.reshape(Shape15);
 Tensor1 = BulkTensors[L-1];
 Indices0[0] = 2; Indices1[0] = 1;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = BulkTensorsPEPO[L-1];
 Indices02[0] = 2; Indices02[1] = 7; Indices12[0] = 1; Indices12[1] = 4;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 BulkTensors[L-1].complexConjugate(Tensor1);
 Indices02[1] = 9;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 MPS1.get(L-1, Tensor1);
 Tensor1.getShape(Shape3);
 Shape15[0] = Shape3[0]; Shape15[1] = 1; Shape15[2] = Shape5[3]; Shape15[3] = Shape6[3]; Shape15[4] = Shape5[3];
 Tensor1.reshape(Shape15);
 Tensor1.complexConjugate();
 Indices03[0] = 4; Indices03[1] = 7; Indices03[2] = 10; Indices13[0] = 2; Indices13[1] = 3; Indices13[2] = 4;
 Tensor0.contract(Indices03, Tensor1, Indices13);
 Tensor0.getShape(Shape10);
 Shape5[0] = Shape10[0]; Shape5[1] = Shape10[2]; Shape5[2] = Shape10[4]; Shape5[3] = Shape10[6]; Shape5[4] = Shape10[8];
 Tensor0.reshape(Shape5);
 bTensors[L-1] = Tensor0;
 Indices0[0] = 0; Indices1[0] = 1;
 Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 2; Indices12[1] = 1;
 Indices03[0] = 0; Indices03[2] = 8; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
 Indices04[0] = 0; Indices04[1] = 3; Indices04[2] = 5; Indices04[3] = 7; Indices14[0] = 1; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
 for (int pos = L-2; pos > 0; pos--)
 {
  MPS0.get(pos, Tensor1);
  Tensor1.getShape(Shape3);
  BulkTensors[pos].getShape(Shape5);
  BulkTensorsPEPO[pos].getShape(Shape6);
  Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape5[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
  Tensor1.reshape(Shape15);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor1 = BulkTensors[pos];
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor1 = BulkTensorsPEPO[pos];
  Indices03[1] = 4;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  BulkTensors[pos].complexConjugate(Tensor1);
  Indices03[1] = 3;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  MPS1.get(pos, Tensor1);
  Tensor1.getShape(Shape3);
  Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape5[3]; Shape15[3] = Shape6[3]; Shape15[4] = Shape5[3];
  Tensor1.reshape(Shape15);
  Tensor1.complexConjugate();
  Tensor0.contract(Indices04, Tensor1, Indices14);
  bTensors[pos] = Tensor0;
 }
// 4. compute initial error:
 MPS0.get(0, Tensor1);
 Tensor1.getShape(Shape3);
 BulkTensors[0].getShape(Shape5);
 BulkTensorsPEPO[0].getShape(Shape6);
 Shape15[0] = 1; Shape15[1] = Shape3[1]; Shape15[2] = Shape5[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
 Tensor1.reshape(Shape15);
 Indices0[0] = 0; Indices1[0] = 1;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = BulkTensors[0];
 Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 2; Indices12[1] = 1;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 Tensor1 = BulkTensorsPEPO[0];
 Indices03[0] = 0; Indices03[1] = 4; Indices03[2] = 8; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
 Tensor0.contract(Indices03, Tensor1, Indices13);
 BulkTensors[0].complexConjugate(Tensor1);
 Indices03[1] = 3;
 Tensor0.contract(Indices03, Tensor1, Indices13);
 MPS1.get(0, Tensor1);
 Tensor1.getShape(Shape3);
 Shape15[0] = 1; Shape15[1] = Shape3[1]; Shape15[2] = Shape5[3]; Shape15[3] = Shape6[3]; Shape15[4] = Shape5[3];
 Tensor1.reshape(Shape15);
 Tensor1.complexConjugate();
 Indices04[0] = 0; Indices04[1] = 3; Indices04[2] = 5; Indices04[3] = 7; Indices14[0] = 1; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
 Tensor0.contract(Indices04, Tensor1, Indices14);
 bValue = Tensor0.get(0);
 MPS1.get(0, Tensor0);
 NValue = Tensor0.scalarProduct();
 error0 = MathAuxiliary::convertToDouble(NValue)-2.0*MathAuxiliary::convertToDouble(bValue);
// 5. sweep:
 numSweepsDone = 0;
 SweepDirection = "right";
 while (numSweepsDone < maxNumSweeps)
 {
// 5.a) sweep to the right:
  if (SweepDirection == "right")
  {
// 5.a)1. position 0:
// 5.a)1. - compute bTensor:
   Tensor0 = bTensors[1];
   MPS0.get(0, Tensor1);
   Tensor1.getShape(Shape3);
   BulkTensors[0].getShape(Shape5);
   BulkTensorsPEPO[0].getShape(Shape6);
   Shape15[0] = 1; Shape15[1] = Shape3[1]; Shape15[2] = Shape5[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
   Tensor1.reshape(Shape15);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BulkTensors[0];
   Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 2; Indices12[1] = 1;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor1 = BulkTensorsPEPO[0];
   Indices03[0] = 0; Indices03[1] = 4; Indices03[2] = 8; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   BulkTensors[0].complexConjugate(Tensor1);
   Indices03[1] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
// 5.a)1. - update tensor:
   Tensor0.getShape(Shape8);
   Shape3[0] = 1; Shape3[1] = Shape8[0]; Shape3[2] = Shape8[3]*Shape8[5]*Shape8[7];
   Tensor0.reshape(Shape3);
   MPS1.set(0, Tensor0);
// 5.a)1. - compute error:
   NValue = Tensor0.scalarProduct();
   error1 = -MathAuxiliary::convertToDouble(NValue);
   epsAchieved = abs((error1-error0)/error1);
   if (epsAchieved < eps)
    return;
   error0 = error1;
// 5.a)1. - normalize tensor:
   if (DMPS1 == 1)
    tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
   else
   {
    Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
    Tensor0.QRDecompose(Indices2, Indices1, TensorR);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
   }
   MPS1.set(0, Tensor0);
// 5.a)1. - set bTensors:
   MPS0.get(0, Tensor0);
   Tensor0.getShape(Shape3);
   Shape15[0] = 1; Shape15[1] = Shape3[1]; Shape15[2] = Shape5[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
   Tensor0.reshape(Shape15);
   Tensor1 = BulkTensors[0];
   Indices0[0] = 2; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BulkTensorsPEPO[0];
   Indices02[0] = 2; Indices02[1] = 7; Indices12[0] = 1; Indices12[1] = 4;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   BulkTensors[0].complexConjugate(Tensor1);
   Indices02[1] = 9;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   MPS1.get(0, Tensor1);
   Tensor1.getShape(Shape3);
   Shape15[0] = 1; Shape15[1] = Shape3[1]; Shape15[2] = Shape5[3]; Shape15[3] = Shape6[3]; Shape15[4] = Shape5[3];
   Tensor1.reshape(Shape15);
   Tensor1.complexConjugate();
   Indices03[0] = 4; Indices03[1] = 7; Indices03[2] = 10; Indices13[0] = 2; Indices13[1] = 3; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape10);
   Shape5[0] = Shape10[1]; Shape5[1] = Shape10[3]; Shape5[2] = Shape10[5]; Shape5[3] = Shape10[7]; Shape5[4] = Shape10[9];
   Tensor0.reshape(Shape5);
   bTensors[0] = Tensor0;
// 5.a)2. 0 < pos < L-1:
   for (int pos = 1; pos < L-1; pos++)
   {
// 5.a)2. - compute bTensor:
    Tensor0 = bTensors[pos-1];
    MPS0.get(pos, Tensor1);
    Tensor1.getShape(Shape3);
    BulkTensors[pos].getShape(Shape5);
    BulkTensorsPEPO[pos].getShape(Shape6);
    Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape5[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
    Tensor1.reshape(Shape15);
    Indices0[0] = 0; Indices1[0] = 0;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = BulkTensors[pos];
    Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 1;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    Tensor1 = BulkTensorsPEPO[pos];
    Indices03[0] = 0; Indices03[1] = 4; Indices03[2] = 8; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 4;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    BulkTensors[pos].complexConjugate(Tensor1);
    Indices03[1] = 3;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    bTensor = Tensor0;
// 5.a)2. - update tensor:
    Tensor1 = bTensors[pos+1];
    Indices04[0] = 1; Indices04[1] = 2; Indices04[2] = 4; Indices04[3] = 6; Indices14[0] = 0; Indices14[1] = 1; Indices14[2] = 2; Indices14[3] = 3;
    Tensor0.contract(Indices04, Tensor1, Indices14);
    Tensor0.getShape(Shape15);
    Shape3[0] = Shape15[0]; Shape3[1] = Shape15[1]*Shape15[2]*Shape15[3]; Shape3[2] = Shape15[4];
    Tensor0.reshape(Shape3);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
    MPS1.set(pos, Tensor0);
// 5.a)2. - compute error:
    NValue = Tensor0.scalarProduct();
    error1 = -MathAuxiliary::convertToDouble(NValue);
    epsAchieved = abs((error1-error0)/error1);
    if (epsAchieved < eps)
     return;
    error0 = error1;
// 5.a)2. - normalize tensor:
    if (DMPS1 == 1)
     tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
    else
    {
     Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
     Tensor0.QRDecompose(Indices2, Indices1, TensorR);
     Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
     Tensor0.permute(Order3);
    }
    MPS1.set(pos, Tensor0);
// 5.a)2. - set bTensors:
    Tensor0.getShape(Shape3);
    Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape5[3]; Shape15[3] = Shape6[3]; Shape15[4] = Shape5[3];
    Tensor0.reshape(Shape15);
    Tensor0.complexConjugate();
    Indices04[0] = 0; Indices04[1] = 3; Indices04[2] = 5; Indices04[3] = 7; Indices14[0] = 0; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
    bTensor.contract(Indices04, Tensor0, Indices14);
    bTensors[pos] = bTensor;
   }
// 5.a)3. position L-1:
   MPS1.get(L-1, Tensor0);
   if (DMPS1 == 1)
    Tensor0.multiply(tensorFrobeniusNorm);
   else
   {
    Indices0[0] = 1; Indices1[0] = 0;
    TensorR.contract(Indices0, Tensor0, Indices1);
    Tensor0 = TensorR;
   }
   MPS1.set(L-1, Tensor0);
   SweepDirection = "left";
  }
// 5.b) sweep to the left:
  else if (SweepDirection == "left")
  {
// 5.b)1. position L-1:
// 5.b)1. - compute bTensor:
   Tensor0 = bTensors[L-2];
   MPS0.get(L-1, Tensor1);
   Tensor1.getShape(Shape3);
   BulkTensors[L-1].getShape(Shape5);
   BulkTensorsPEPO[L-1].getShape(Shape6);
   Shape15[0] = Shape3[0]; Shape15[1] = 1; Shape15[2] = Shape5[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
   Tensor1.reshape(Shape15);
   Indices0[0] = 0; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BulkTensors[L-1];
   Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 1;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor1 = BulkTensorsPEPO[L-1];
   Indices03[0] = 0; Indices03[1] = 4; Indices03[2] = 8; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   BulkTensors[L-1].complexConjugate(Tensor1);
   Indices03[1] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
// 5.b)1. - update tensor:
   Tensor0.getShape(Shape8);
   Shape3[0] = Shape8[0]; Shape3[1] = 1; Shape3[2] = Shape8[3]*Shape8[5]*Shape8[7];
   Tensor0.reshape(Shape3);
   MPS1.set(L-1, Tensor0);
// 5.b)1. - compute error:
   NValue = Tensor0.scalarProduct();
   error1 = -MathAuxiliary::convertToDouble(NValue);
   epsAchieved = abs((error1-error0)/error1);
   if (epsAchieved < eps)
    return;
   error0 = error1;
// 5.b)1. - normalize tensor:
   if (DMPS1 == 1)
    tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
   else
   {
    Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
    Tensor0.LQDecompose(Indices1, Indices2, TensorL);
   }
   MPS1.set(L-1, Tensor0);
// 5.b)1. - set bTensors:
   MPS0.get(L-1, Tensor0);
   Tensor0.getShape(Shape3);
   Shape15[0] = Shape3[0]; Shape15[1] = 1; Shape15[2] = Shape5[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
   Tensor0.reshape(Shape15);
   Tensor1 = BulkTensors[L-1];
   Indices0[0] = 2; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BulkTensorsPEPO[L-1];
   Indices02[0] = 2; Indices02[1] = 7; Indices12[0] = 1; Indices12[1] = 4;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   BulkTensors[L-1].complexConjugate(Tensor1);
   Indices02[1] = 9;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   MPS1.get(L-1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape15[0] = Shape3[0]; Shape15[1] = 1; Shape15[2] = Shape5[3]; Shape15[3] = Shape6[3]; Shape15[4] = Shape5[3];
   Tensor1.reshape(Shape15);
   Tensor1.complexConjugate();
   Indices03[0] = 4; Indices03[1] = 7; Indices03[2] = 10; Indices13[0] = 2; Indices13[1] = 3; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape10);
   Shape5[0] = Shape10[0]; Shape5[1] = Shape10[2]; Shape5[2] = Shape10[4]; Shape5[3] = Shape10[6]; Shape5[4] = Shape10[8];
   Tensor0.reshape(Shape5);
   bTensors[L-1] = Tensor0;
// 5.b)2. 0 < pos < L-1:
   for (int pos = L-2; pos > 0; pos--)
   {
// 5.b)2. - compute bTensor:
    Tensor0 = bTensors[pos+1];
    MPS0.get(pos, Tensor1);
    Tensor1.getShape(Shape3);
    BulkTensors[pos].getShape(Shape5);
    BulkTensorsPEPO[pos].getShape(Shape6);
    Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape5[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
    Tensor1.reshape(Shape15);
    Indices0[0] = 0; Indices1[0] = 1;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = BulkTensors[pos];
    Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 2; Indices12[1] = 1;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    Tensor1 = BulkTensorsPEPO[pos];
    Indices03[0] = 0; Indices03[1] = 4; Indices03[2] = 8; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    BulkTensors[pos].complexConjugate(Tensor1);
    Indices03[1] = 3;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    bTensor = Tensor0;
// 5.b)2. - update tensor:
    Tensor1 = bTensors[pos-1];
    Indices04[0] = 1; Indices04[1] = 2; Indices04[2] = 4; Indices04[3] = 6; Indices14[0] = 0; Indices14[1] = 1; Indices14[2] = 2; Indices14[3] = 3;
    Tensor0.contract(Indices04, Tensor1, Indices14);
    Tensor0.getShape(Shape15);
    Shape3[0] = Shape15[0]; Shape3[1] = Shape15[1]*Shape15[2]*Shape15[3]; Shape3[2] = Shape15[4];
    Tensor0.reshape(Shape3);
    Order3[0] = 2; Order3[1] = 0; Order3[2] = 1;
    Tensor0.permute(Order3);
    MPS1.set(pos, Tensor0);
// 5.b)2. - compute error:
    NValue = Tensor0.scalarProduct();
    error1 = -MathAuxiliary::convertToDouble(NValue);
    epsAchieved = abs((error1-error0)/error1);
    if (epsAchieved < eps)
     return;
    error0 = error1;
// 5.b)2. - normalize tensor:
    if (DMPS1 == 1)
     tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
    else
    {
     Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
     Tensor0.LQDecompose(Indices1, Indices2, TensorL);
    }
    MPS1.set(pos, Tensor0);
// 5.b)2. - set bTensors:
    Tensor0.getShape(Shape3);
    Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape5[3]; Shape15[3] = Shape6[3]; Shape15[4] = Shape5[3];
    Tensor0.reshape(Shape15);
    Tensor0.complexConjugate();
    Indices04[0] = 0; Indices04[1] = 3; Indices04[2] = 5; Indices04[3] = 7; Indices14[0] = 1; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
    bTensor.contract(Indices04, Tensor0, Indices14);
    bTensors[pos] = bTensor;
   }
// 5.b)3. position 0:
   MPS1.get(0, Tensor0);
   if (DMPS1 == 1)
    Tensor0.multiply(tensorFrobeniusNorm);
   else
   {
    Indices0[0] = 1; Indices1[0] = 0;
    Tensor0.contract(Indices0, TensorL, Indices1);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
   }
   MPS1.set(0, Tensor0);
   SweepDirection = "right";
  }
  numSweepsDone++;
 }
}

template<class T> void PEPS<T>::getBoundaryMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary,
                                               double eps, unsigned int maxNumSweeps,
                                               double& epsAchieved, unsigned int& numSweepsDone,
                                               MPS<T>& MPS0) const
{
#ifdef DEBUG
 if ((this->BC != "open") || (this->N[0] == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getBoundaryMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary, " <<
                         "double eps, unsigned int maxNumSweeps, " <<
                         "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "((this->BC != open) || (this->N[0] == 0))." << endl;
  exit(1);
 }
 string BCPEPO0; PEPO0.getBC(BCPEPO0);
 vector<unsigned int> NPEPO0(2); NPEPO0[0] = PEPO0.getNrows(); NPEPO0[1] = PEPO0.getNcols();
 Matrix<unsigned int> dPEPO0; PEPO0.getd(dPEPO0);
 if ((this->BC != BCPEPO0) || (this->N != NPEPO0) || (this->d != dPEPO0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getBoundaryMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary, " <<
                         "double eps, unsigned int maxNumSweeps, " <<
                         "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "((this->BC != PEPO0.BC) || (this->N != PEPO0.N) || (this->d != PEPO0.d))." << endl;
  exit(1);
 }
 if ((this->BC != PEPS0.BC) || (this->N != PEPS0.N) || (this->d != PEPS0.d))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getBoundaryMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary, " <<
                         "double eps, unsigned int maxNumSweeps, " <<
                         "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "((this->BC != PEPS0.BC) || (this->N != PEPS0.N) || (this->d != PEPS0.d))." << endl;
  exit(1);
 }
 if ((Boundary != "left") && (Boundary != "top") && (Boundary != "right") && (Boundary != "bottom"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getBoundaryMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary, " <<
                         "double eps, unsigned int maxNumSweeps, " <<
                         "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "((Boundary != left) && (Boundary != top) && (Boundary != right) && (Boundary != bottom))." << endl;
  exit(1);
 }
 string BCMPS0; MPS0.getBC(BCMPS0);
 if ((BCMPS0 != "open") ||
     (((Boundary == "left") || (Boundary == "right")) && (MPS0.getN() != this->N[0])) ||
     (((Boundary == "top") || (Boundary == "bottom")) && (MPS0.getN() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getBoundaryMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary, " <<
                         "double eps, unsigned int maxNumSweeps, " <<
                         "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "((MPS0.BC != open) || (((Boundary == left) || (Boundary == right)) && (MPS0.N != this->N[0])) || " <<
           "(((Boundary == top) || (Boundary == bottom)) && (MPS0.N != this->N[1])))." << endl;
  exit(1);
 }
#endif
 unsigned int DMPS0, L;
 double error0, error1;
 T bValue, NValue, tensorFrobeniusNorm, tensorFrobeniusNormProduct;
 string SweepDirection;
 vector<unsigned int> dMPS0, Indices0(1), Indices02(2), Indices03(3), Indices04(4), Indices1(1), Indices12(2), Indices13(3), Indices14(4), Indices2(2);
 vector<unsigned int> Order3(3), Order5(5), Order6(6), Order9(9), Shape3(3), Shape4(4), Shape04(4), Shape5(5), Shape6(6), Shape7(7), Shape8(8), Shape9(9);
 Tensor<T> bTensor, Tensor0, Tensor1, TensorL, TensorR;
 vector< Tensor<T> > BoundaryTensors, BoundaryTensors0, BoundaryTensorsPEPO, bTensors;
// 0. initialize:
// 0. - set boundary-length L:
 if ((Boundary == "left") || (Boundary == "right"))
  L = this->N[0];
 else if ((Boundary == "top") || (Boundary == "bottom"))
  L = this->N[1];
// 0. - get boundary-tensors clockwise, permute them as "top" boundary, and write the ones from this PEPS in BoundaryTensors,
//      the ones from PEPS0 in BoundaryTensors0, and the ones from PEPO0 in BoundaryTensorsPEPO:
 BoundaryTensors = vector< Tensor<T> >(L);
 BoundaryTensors0 = vector< Tensor<T> >(L);
 BoundaryTensorsPEPO = vector< Tensor<T> >(L);
 if (Boundary == "left")
 {
  Order5[0] = 3; Order5[1] = 0; Order5[2] = 1; Order5[3] = 2; Order5[4] = 4;
  Order6[0] = 3; Order6[1] = 0; Order6[2] = 1; Order6[3] = 2; Order6[4] = 4; Order6[5] = 5;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[this->N[0]-1-pos];
   Tensor0.permute(Order5);
   BoundaryTensors[pos] = Tensor0;
   Tensor0 = PEPS0.Tensors[this->N[0]-1-pos];
   Tensor0.permute(Order5);
   BoundaryTensors0[pos] = Tensor0;
   PEPO0.get(this->N[0]-1-pos, 0, Tensor0);
   Tensor0.permute(Order6);
   BoundaryTensorsPEPO[pos] = Tensor0;
  }
 }
 else if (Boundary == "top")
 {
  for (int pos = 0; pos < L; pos++)
  {
   BoundaryTensors[pos] = this->Tensors[this->N[0]*pos];
   BoundaryTensors0[pos] = PEPS0.Tensors[this->N[0]*pos];
   PEPO0.get(0, pos, BoundaryTensorsPEPO[pos]);
  }
 }
 else if (Boundary == "right")
 {
  Order5[0] = 1; Order5[1] = 2; Order5[2] = 3; Order5[3] = 0; Order5[4] = 4;
  Order6[0] = 1; Order6[1] = 2; Order6[2] = 3; Order6[3] = 0; Order6[4] = 4; Order6[5] = 5;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[pos+this->N[0]*(this->N[1]-1)];
   Tensor0.permute(Order5);
   BoundaryTensors[pos] = Tensor0;
   Tensor0 = PEPS0.Tensors[pos+this->N[0]*(this->N[1]-1)];
   Tensor0.permute(Order5);
   BoundaryTensors0[pos] = Tensor0;
   PEPO0.get(pos, this->N[1]-1, Tensor0);
   Tensor0.permute(Order6);
   BoundaryTensorsPEPO[pos] = Tensor0;
  }
 }
 else if (Boundary == "bottom")
 {
  Order5[0] = 2; Order5[1] = 3; Order5[2] = 0; Order5[3] = 1; Order5[4] = 4;
  Order6[0] = 2; Order6[1] = 3; Order6[2] = 0; Order6[3] = 1; Order6[4] = 4; Order6[5] = 5;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[this->N[0]-1+this->N[0]*(this->N[1]-1-pos)];
   Tensor0.permute(Order5);
   BoundaryTensors[pos] = Tensor0;
   Tensor0 = PEPS0.Tensors[this->N[0]-1+this->N[0]*(this->N[1]-1-pos)];
   Tensor0.permute(Order5);
   BoundaryTensors0[pos] = Tensor0;
   PEPO0.get(this->N[0]-1, this->N[1]-1-pos, Tensor0);
   Tensor0.permute(Order6);
   BoundaryTensorsPEPO[pos] = Tensor0;
  }
 }
// 0. - reshape boundary-tensors eliminating index 1:
 for (int pos = 0; pos < L; pos++)
 {
  BoundaryTensors[pos].getShape(Shape5);
  Shape4[0] = Shape5[0]; Shape4[1] = Shape5[2]; Shape4[2] = Shape5[3]; Shape4[3] = Shape5[4];
  BoundaryTensors[pos].reshape(Shape4);
  BoundaryTensors0[pos].getShape(Shape5);
  Shape4[0] = Shape5[0]; Shape4[1] = Shape5[2]; Shape4[2] = Shape5[3]; Shape4[3] = Shape5[4];
  BoundaryTensors0[pos].reshape(Shape4);
  BoundaryTensorsPEPO[pos].getShape(Shape6);
  Shape5[0] = Shape6[0]; Shape5[1] = Shape6[2]; Shape5[2] = Shape6[3]; Shape5[3] = Shape6[4]; Shape5[4] = Shape6[5];
  BoundaryTensorsPEPO[pos].reshape(Shape5);
 }
// 0. - check MPS0.d:
 MPS0.getd(dMPS0);
#ifdef DEBUG
 for (int pos = 0; pos < L; pos++)
 {
  BoundaryTensors0[pos].getShape(Shape04);
  BoundaryTensorsPEPO[pos].getShape(Shape5);
  BoundaryTensors[pos].getShape(Shape4);
  if (dMPS0[pos] != Shape04[2]*Shape5[2]*Shape4[2])
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void PEPS<T>::" <<
           "getBoundaryMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary, " <<
                          "double eps, unsigned int maxNumSweeps, " <<
                          "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
           "(MPS0.d[" << pos << "] != PEPS0.D[" << pos << "]*PEPO0.D[" << pos << "]*this->D[" << pos << "])." << endl;
   exit(1);
  }
 }
#endif
// 1. MPS0.D >= PEPS0.D*PEPO0.D*this->D:
 DMPS0 = MPS0.getD();
 if (DMPS0 >= PEPS0.D*PEPO0.getD()*this->D)
 {
  MPS0 = MPS<T>("open", L, dMPS0, PEPS0.D*PEPO0.getD()*this->D);
  Indices1[0] = 3;
  Order9[0] = 0; Order9[1] = 3; Order9[2] = 6; Order9[3] = 1; Order9[4] = 4; Order9[5] = 7; Order9[6] = 2; Order9[7] = 5; Order9[8] = 8;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = BoundaryTensors0[pos];
   Tensor1 = BoundaryTensorsPEPO[pos];
   Indices0[0] = 3;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   BoundaryTensors[pos].complexConjugate(Tensor1);
   Indices0[0] = 6;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order9);
   Tensor0.getShape(Shape9);
   Shape3[0] = Shape9[0]*Shape9[1]*Shape9[2]; Shape3[1] = Shape9[3]*Shape9[4]*Shape9[5]; Shape3[2] = Shape9[6]*Shape9[7]*Shape9[8];
   Tensor0.reshape(Shape3);
   MPS0.set(pos, Tensor0);
  }
  MPS0.bringIntoNormalShape();
  MPS0.setD(DMPS0);
  epsAchieved = 0.0;
  numSweepsDone = 0;
  return;
 }
// 2. bring MPS0 into normal form from right to left:
// 2.a) MPS0.D == 1:
 if (DMPS0 == 1)
 {
  tensorFrobeniusNormProduct = 1.0;
  for (int pos = L-1; pos > 0; pos--)
  {
   MPS0.get(pos, Tensor0);
   tensorFrobeniusNormProduct *= Tensor0.frobeniusNormalize();
   MPS0.set(pos, Tensor0);
  }
  MPS0.get(0, Tensor0);
  Tensor0.multiply(tensorFrobeniusNormProduct);
  MPS0.set(0, Tensor0);
 }
// 2.b) MPS0.D != 1:
 else
 {
  Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
  Indices0[0] = 1;
  Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
  MPS0.get(L-1, Tensor0);
  for (int pos = L-1; pos > 0; pos--)
  {
   Tensor0.LQDecompose(Indices1, Indices2, TensorL);
   MPS0.set(pos, Tensor0);
   MPS0.get(pos-1, Tensor0);
   Tensor0.contract(Indices0, TensorL, Indices1);
   Tensor0.permute(Order3);
  }
  MPS0.set(0, Tensor0);
 }
// 3. construct right half of b for each site:
 bTensors = vector< Tensor<T> >(L);
 Tensor0 = BoundaryTensors0[L-1];
 Tensor1 = BoundaryTensorsPEPO[L-1];
 Indices0[0] = 3; Indices1[0] = 3;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 BoundaryTensors[L-1].complexConjugate(Tensor1);
 Indices0[0] = 6;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor0.getShape(Shape9);
 MPS0.get(L-1, Tensor1);
 Tensor1.getShape(Shape3);
 Shape5[0] = Shape3[0]; Shape5[1] = 1; Shape5[2] = Shape9[2]; Shape5[3] = Shape9[5]; Shape5[4] = Shape9[8];
 Tensor1.reshape(Shape5);
 Tensor1.complexConjugate();
 Indices03[0] = 2; Indices03[1] = 5; Indices03[2] = 8; Indices13[0] = 2; Indices13[1] = 3; Indices13[2] = 4;
 Tensor0.contract(Indices03, Tensor1, Indices13);
 Tensor0.getShape(Shape8);
 Shape4[0] = Shape8[0]; Shape4[1] = Shape8[2]; Shape4[2] = Shape8[4]; Shape4[3] = Shape8[6];
 Tensor0.reshape(Shape4);
 bTensors[L-1] = Tensor0;
 Indices0[0] = 0; Indices1[0] = 1;
 Indices02[0] = 0; Indices12[0] = 1; Indices12[1] = 3;
 Indices04[0] = 0; Indices04[1] = 2; Indices04[2] = 4; Indices04[3] = 6; Indices14[0] = 1; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
 for (int pos = L-2; pos > 0; pos--)
 {
  Tensor1 = BoundaryTensors0[pos];
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor1 = BoundaryTensorsPEPO[pos];
  Indices02[1] = 5;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  BoundaryTensors[pos].complexConjugate(Tensor1);
  Indices02[1] = 6;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape7);
  MPS0.get(pos, Tensor1);
  Tensor1.getShape(Shape3);
  Shape5[0] = Shape3[0]; Shape5[1] = Shape3[1]; Shape5[2] = Shape7[2]; Shape5[3] = Shape7[4]; Shape5[4] = Shape7[6];
  Tensor1.reshape(Shape5);
  Tensor1.complexConjugate();
  Tensor0.contract(Indices04, Tensor1, Indices14);
  bTensors[pos] = Tensor0;
 }
// 4. compute initial error:
 Tensor1 = BoundaryTensors0[0];
 Indices0[0] = 0; Indices1[0] = 1;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = BoundaryTensorsPEPO[0];
 Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 1; Indices12[1] = 3;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 BoundaryTensors[0].complexConjugate(Tensor1);
 Indices02[1] = 6;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 Tensor0.getShape(Shape7);
 MPS0.get(0, Tensor1);
 Tensor1.getShape(Shape3);
 Shape5[0] = 1; Shape5[1] = Shape3[1]; Shape5[2] = Shape7[2]; Shape5[3] = Shape7[4]; Shape5[4] = Shape7[6];
 Tensor1.reshape(Shape5);
 Tensor1.complexConjugate();
 Indices04[0] = 0; Indices04[1] = 2; Indices04[2] = 4; Indices04[3] = 6; Indices14[0] = 1; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
 Tensor0.contract(Indices04, Tensor1, Indices14);
 bValue = Tensor0.get(0);
 MPS0.get(0, Tensor0);
 NValue = Tensor0.scalarProduct();
 error0 = MathAuxiliary::convertToDouble(NValue)-2.0*MathAuxiliary::convertToDouble(bValue);
// 5. sweep:
 numSweepsDone = 0;
 SweepDirection = "right";
 while (numSweepsDone < maxNumSweeps)
 {
// 5.a) sweep to the right:
  if (SweepDirection == "right")
  {
// 5.a)1. position 0:
// 5.a)1. - compute bTensor:
   Tensor0 = bTensors[1];
   Tensor1 = BoundaryTensors0[0];
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BoundaryTensorsPEPO[0];
   Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 1; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   BoundaryTensors[0].complexConjugate(Tensor1);
   Indices02[1] = 6;
   Tensor0.contract(Indices02, Tensor1, Indices12);
// 5.a)1. - update tensor:
   Tensor0.getShape(Shape7);
   Shape3[0] = 1; Shape3[1] = Shape7[0]; Shape3[2] = Shape7[2]*Shape7[4]*Shape7[6];
   Tensor0.reshape(Shape3);
   MPS0.set(0, Tensor0);
// 5.a)1. - compute error:
   NValue = Tensor0.scalarProduct();
   error1 = -MathAuxiliary::convertToDouble(NValue);
   epsAchieved = abs((error1-error0)/error1);
   if (epsAchieved < eps)
    return;
   error0 = error1;
// 5.a)1. - normalize tensor:
   if (DMPS0 == 1)
    tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
   else
   {
    Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
    Tensor0.QRDecompose(Indices2, Indices1, TensorR);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
   }
   MPS0.set(0, Tensor0);
// 5.a)1. - set bTensors:
   Tensor0 = BoundaryTensors0[0];
   Tensor1 = BoundaryTensorsPEPO[0];
   Indices0[0] = 3; Indices1[0] = 3;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   BoundaryTensors[0].complexConjugate(Tensor1);
   Indices0[0] = 6;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.getShape(Shape9);
   MPS0.get(0, Tensor1);
   Tensor1.getShape(Shape3);
   Shape5[0] = 1; Shape5[1] = Shape3[1]; Shape5[2] = Shape9[2]; Shape5[3] = Shape9[5]; Shape5[4] = Shape9[8];
   Tensor1.reshape(Shape5);
   Tensor1.complexConjugate();
   Indices03[0] = 2; Indices03[1] = 5; Indices03[2] = 8; Indices13[0] = 2; Indices13[1] = 3; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape8);
   Shape4[0] = Shape8[1]; Shape4[1] = Shape8[3]; Shape4[2] = Shape8[5]; Shape4[3] = Shape8[7];
   Tensor0.reshape(Shape4);
   bTensors[0] = Tensor0;
// 5.a)2. 0 < position < L-1:
   for (int pos = 1; pos < L-1; pos++)
   {
// 5.a)2. - compute bTensor:
    Tensor0 = bTensors[pos-1];
    Tensor1 = BoundaryTensors0[pos];
    Indices0[0] = 0; Indices1[0] = 0;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = BoundaryTensorsPEPO[pos];
    Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 3;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    BoundaryTensors[pos].complexConjugate(Tensor1);
    Indices02[1] = 6;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    bTensor = Tensor0;
// 5.a)2. - update tensor:
    Tensor1 = bTensors[pos+1];
    Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 5; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 2;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    Tensor0.getShape(Shape5);
    Shape3[0] = Shape5[0]; Shape3[1] = Shape5[1]*Shape5[2]*Shape5[3]; Shape3[2] = Shape5[4];
    Tensor0.reshape(Shape3);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
    MPS0.set(pos, Tensor0);
// 5.a)2. - compute error:
    NValue = Tensor0.scalarProduct();
    error1 = -MathAuxiliary::convertToDouble(NValue);
    epsAchieved = abs((error1-error0)/error1);
    if (epsAchieved < eps)
     return;
    error0 = error1;
// 5.a)2. - normalize tensor:
    if (DMPS0 == 1)
     tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
    else
    {
     Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
     Tensor0.QRDecompose(Indices2, Indices1, TensorR);
     Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
     Tensor0.permute(Order3);
    }
    MPS0.set(pos, Tensor0);
// 5.a)2. - set bTensors:
    bTensor.getShape(Shape7);
    Tensor0.getShape(Shape3);
    Shape5[0] = Shape3[0]; Shape5[1] = Shape3[1]; Shape5[2] = Shape7[2]; Shape5[3] = Shape7[4]; Shape5[4] = Shape7[6];
    Tensor0.reshape(Shape5);
    Tensor0.complexConjugate();
    Indices04[0] = 0; Indices04[1] = 2; Indices04[2] = 4; Indices04[3] = 6; Indices14[0] = 0; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
    bTensor.contract(Indices04, Tensor0, Indices14);
    bTensors[pos] = bTensor;
   }
// 5.a)3. position L-1:
   MPS0.get(L-1, Tensor0);
   if (DMPS0 == 1)
    Tensor0.multiply(tensorFrobeniusNorm);
   else
   {
    Indices0[0] = 1; Indices1[0] = 0;
    TensorR.contract(Indices0, Tensor0, Indices1);
    Tensor0 = TensorR;
   }
   MPS0.set(L-1, Tensor0);
   SweepDirection = "left";
  }
// 5.b) sweep to the left:
  else if (SweepDirection == "left")
  {
// 5.b)1. position L-1:
// 5.b)1. - compute bTensor:
   Tensor0 = bTensors[L-2];
   Tensor1 = BoundaryTensors0[L-1];
   Indices0[0] = 0; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BoundaryTensorsPEPO[L-1];
   Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   BoundaryTensors[L-1].complexConjugate(Tensor1);
   Indices02[1] = 6;
   Tensor0.contract(Indices02, Tensor1, Indices12);
// 5.b)1. - update tensor:
   Tensor0.getShape(Shape7);
   Shape3[0] = Shape7[0]; Shape3[1] = 1; Shape3[2] = Shape7[2]*Shape7[4]*Shape7[6];
   Tensor0.reshape(Shape3);
   MPS0.set(L-1, Tensor0);
// 5.b)1. - compute error:
   NValue = Tensor0.scalarProduct();
   error1 = -MathAuxiliary::convertToDouble(NValue);
   epsAchieved = abs((error1-error0)/error1);
   if (epsAchieved < eps)
    return;
   error0 = error1;
// 5.b)1. - normalize tensor:
   if (DMPS0 == 1)
    tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
   else
   {
    Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
    Tensor0.LQDecompose(Indices1, Indices2, TensorL);
   }
   MPS0.set(L-1, Tensor0);
// 5.b)1. - set bTensors:
   Tensor0 = BoundaryTensors0[L-1];
   Tensor1 = BoundaryTensorsPEPO[L-1];
   Indices0[0] = 3; Indices1[0] = 3;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   BoundaryTensors[L-1].complexConjugate(Tensor1);
   Indices0[0] = 6;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.getShape(Shape9);
   MPS0.get(L-1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape5[0] = Shape3[0]; Shape5[1] = 1; Shape5[2] = Shape9[2]; Shape5[3] = Shape9[5]; Shape5[4] = Shape9[8];
   Tensor1.reshape(Shape5);
   Tensor1.complexConjugate();
   Indices03[0] = 2; Indices03[1] = 5; Indices03[2] = 8; Indices13[0] = 2; Indices13[1] = 3; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape8);
   Shape4[0] = Shape8[0]; Shape4[1] = Shape8[2]; Shape4[2] = Shape8[4]; Shape4[3] = Shape8[6];
   Tensor0.reshape(Shape4);
   bTensors[L-1] = Tensor0;
// 5.b)2. 0 < position < L-1:
   for (int pos = L-2; pos > 0; pos--)
   {
// 5.b)2. - compute bTensor:
    Tensor0 = bTensors[pos+1];
    Tensor1 = BoundaryTensors0[pos];
    Indices0[0] = 0; Indices1[0] = 1;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = BoundaryTensorsPEPO[pos];
    Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 1; Indices12[1] = 3;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    BoundaryTensors[pos].complexConjugate(Tensor1);
    Indices02[1] = 6;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    bTensor = Tensor0;
// 5.b)2. - update tensor:
    Tensor1 = bTensors[pos-1];
    Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 5; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 2;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    Tensor0.getShape(Shape5);
    Shape3[0] = Shape5[0]; Shape3[1] = Shape5[1]*Shape5[2]*Shape5[3]; Shape3[2] = Shape5[4];
    Tensor0.reshape(Shape3);
    Order3[0] = 2; Order3[1] = 0; Order3[2] = 1;
    Tensor0.permute(Order3);
    MPS0.set(pos, Tensor0);
// 5.b)2. - compute error:
    NValue = Tensor0.scalarProduct();
    error1 = -MathAuxiliary::convertToDouble(NValue);
    epsAchieved = abs((error1-error0)/error1);
    if (epsAchieved < eps)
     return;
    error0 = error1;
// 5.b)2. - normalize tensor:
    if (DMPS0 == 1)
     tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
    else
    {
     Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
     Tensor0.LQDecompose(Indices1, Indices2, TensorL);
    }
    MPS0.set(pos, Tensor0);
// 5.b)2. - set bTensors:
    bTensor.getShape(Shape7);
    Tensor0.getShape(Shape3);
    Shape5[0] = Shape3[0]; Shape5[1] = Shape3[1]; Shape5[2] = Shape7[2]; Shape5[3] = Shape7[4]; Shape5[4] = Shape7[6];
    Tensor0.reshape(Shape5);
    Tensor0.complexConjugate();
    Indices04[0] = 0; Indices04[1] = 2; Indices04[2] = 4; Indices04[3] = 6; Indices14[0] = 1; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
    bTensor.contract(Indices04, Tensor0, Indices14);
    bTensors[pos] = bTensor;
   }
// 5.b)3. position 0:
   MPS0.get(0, Tensor0);
   if (DMPS0 == 1)
    Tensor0.multiply(tensorFrobeniusNorm);
   else
   {
    Indices0[0] = 1; Indices1[0] = 0;
    Tensor0.contract(Indices0, TensorL, Indices1);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
   }
   MPS0.set(0, Tensor0);
   SweepDirection = "right";
  }
  numSweepsDone++;
 } 
}

template<class T> void PEPS<T>::multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0,
                                                   const string& Direction, unsigned int position, const MPS<T>& MPS0,
                                                   double eps, unsigned int maxNumSweeps,
                                                   double& epsAchieved, unsigned int& numSweepsDone,
                                                   MPS<T>& MPS1) const
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Direction, " <<
                             "unsigned int position, const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "(this->N[0] == 0)." << endl;
  exit(1);
 }
 string BCPEPO0; PEPO0.getBC(BCPEPO0);
 vector<unsigned int> NPEPO0(2); NPEPO0[0] = PEPO0.getNrows(); NPEPO0[1] = PEPO0.getNcols();
 Matrix<unsigned int> dPEPO0; PEPO0.getd(dPEPO0);
 if ((this->BC != BCPEPO0) || (this->N != NPEPO0) || (this->d != dPEPO0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Direction, " <<
                             "unsigned int position, const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((this->BC != PEPO0.BC) || (this->N != PEPO0.N) || (this->d != PEPO0.d))." << endl;
  exit(1);
 }
 if ((this->BC != PEPS0.BC) || (this->N != PEPS0.N) || (this->d != PEPS0.d))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Direction, " <<
                             "unsigned int position, const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((this->BC != PEPS0.BC) || (this->N != PEPS0.N) || (this->d != PEPS0.d))." << endl;
  exit(1);
 }
 if ((Direction != "right") && (Direction != "down") && (Direction != "left") && (Direction != "up"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Direction, " <<
                             "unsigned int position, const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((Direction != right) && (Direction != down) && (Direction != left) && (Direction != up))." << endl;
  exit(1);
 }
 if ((((Direction == "right") || (Direction == "left")) && (position > this->N[1]-1)) ||
     (((Direction == "down") || (Direction == "up")) && (position > this->N[0]-1)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Direction, " <<
                             "unsigned int position, const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((((Direction == right) || (Direction == left)) && (position > this->N[1]-1)) || " <<
           "(((Direction == down) || (Direction == up)) && (position > this->N[0]-1)))." << endl;
  exit(1);
 }
 string BCMPS;
 MPS0.getBC(BCMPS);
 if ((BCMPS != "open") ||
     (((Direction == "right") || (Direction == "left")) && (MPS0.getN() != this->N[0])) ||
     (((Direction == "down") || (Direction == "up")) && (MPS0.getN() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Direction, " <<
                             "unsigned int position, const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((MPS0.BC != open) || (((Direction == right) || (Direction == left)) && (MPS0.N != this->N[0])) || " <<
           "(((Direction == down) || (Direction == up)) && (MPS0.N != this->N[1])))." << endl;
  exit(1);
 }
 MPS1.getBC(BCMPS);
 if ((BCMPS != "open") ||
     (((Direction == "right") || (Direction == "left")) && (MPS1.getN() != this->N[0])) ||
     (((Direction == "down") || (Direction == "up")) && (MPS1.getN() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Direction, " <<
                             "unsigned int position, const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                             "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
          "((MPS1.BC != open) || (((Direction == right) || (Direction == left)) && (MPS1.N != this->N[0])) || " <<
           "(((Direction == down) || (Direction == up)) && (MPS1.N != this->N[1])))." << endl;
  exit(1);
 }
#endif
 unsigned int DMPS1, L;
 double error0, error1;
 T bValue, NValue, tensorFrobeniusNorm, tensorFrobeniusNormProduct;
 string SweepDirection;
 vector<unsigned int> dMPS, Indices0(1), Indices02(2), Indices03(3), Indices04(4), Indices1(1), Indices12(2), Indices13(3), Indices14(4), Indices2(2);
 vector<unsigned int> Order3(3), Order5(5), Order6(6), Order11(11);
 vector<unsigned int> Shape3(3), Shape5(5), Shape05(5), Shape15(5), Shape6(6), Shape8(8), Shape10(10), Shape11(11);
 Tensor<T> bTensor, Tensor0, Tensor1, TensorL, TensorR;
 vector< Tensor<T> > BulkTensors, BulkTensors0, BulkTensorsPEPO, bTensors;
// 0. initialize:
// 0. - set length L:
 if ((Direction == "right") || (Direction == "left"))
  L = this->N[0];
 else if ((Direction == "down") || (Direction == "up"))
  L = this->N[1];
// 0. - get bulk-tensors clockwise, permute them as Direction "down", and write the ones from this PEPS in BulkTensors, the
//      ones from PEPS0 in BulkTensors0, and the ones from PEPO0 in BulkTensorsPEPO:
 BulkTensors = vector< Tensor<T> >(L);
 BulkTensors0 = vector< Tensor<T> >(L);
 BulkTensorsPEPO = vector< Tensor<T> >(L);
 if (Direction == "right")
 {
  Order5[0] = 3; Order5[1] = 0; Order5[2] = 1; Order5[3] = 2; Order5[4] = 4;
  Order6[0] = 3; Order6[1] = 0; Order6[2] = 1; Order6[3] = 2; Order6[4] = 4; Order6[5] = 5;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[this->N[0]-1-pos+this->N[0]*position];
   Tensor0.permute(Order5);
   BulkTensors[pos] = Tensor0;
   Tensor0 = PEPS0.Tensors[this->N[0]-1-pos+this->N[0]*position];
   Tensor0.permute(Order5);
   BulkTensors0[pos] = Tensor0;
   PEPO0.get(this->N[0]-1-pos, position, Tensor0);
   Tensor0.permute(Order6);
   BulkTensorsPEPO[pos] = Tensor0;
  }
 }
 else if (Direction == "down")
 {
  for (int pos = 0; pos < L; pos++)
  {
   BulkTensors[pos] = this->Tensors[position+this->N[0]*pos];
   BulkTensors0[pos] = PEPS0.Tensors[position+this->N[0]*pos];
   PEPO0.get(position, pos, BulkTensorsPEPO[pos]);
  }
 }
 else if (Direction == "left")
 {
  Order5[0] = 1; Order5[1] = 2; Order5[2] = 3; Order5[3] = 0; Order5[4] = 4;
  Order6[0] = 1; Order6[1] = 2; Order6[2] = 3; Order6[3] = 0; Order6[4] = 4; Order6[5] = 5;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[pos+this->N[0]*position];
   Tensor0.permute(Order5);
   BulkTensors[pos] = Tensor0;
   Tensor0 = PEPS0.Tensors[pos+this->N[0]*position];
   Tensor0.permute(Order5);
   BulkTensors0[pos] = Tensor0;
   PEPO0.get(pos, position, Tensor0);
   Tensor0.permute(Order6);
   BulkTensorsPEPO[pos] = Tensor0;
  }
 }
 else if (Direction == "up")
 {
  Order5[0] = 2; Order5[1] = 3; Order5[2] = 0; Order5[3] = 1; Order5[4] = 4;
  Order6[0] = 2; Order6[1] = 3; Order6[2] = 0; Order6[3] = 1; Order6[4] = 4; Order6[5] = 5;
  for (int pos = 0; pos < L; pos++)
  {
   Tensor0 = this->Tensors[position+this->N[0]*(this->N[1]-1-pos)];
   Tensor0.permute(Order5);
   BulkTensors[pos] = Tensor0;
   Tensor0 = PEPS0.Tensors[position+this->N[0]*(this->N[1]-1-pos)];
   Tensor0.permute(Order5);
   BulkTensors0[pos] = Tensor0;
   PEPO0.get(position, this->N[1]-1-pos, Tensor0);
   Tensor0.permute(Order6);
   BulkTensorsPEPO[pos] = Tensor0;
  }
 }
// 0. - check MPS0.d:
#ifdef DEBUG
 MPS0.getd(dMPS);
 for (int pos = 0; pos < L; pos++)
 {
  BulkTensors0[pos].getShape(Shape05);
  BulkTensorsPEPO[pos].getShape(Shape6);
  BulkTensors[pos].getShape(Shape5);
  if (dMPS[pos] != Shape05[1]*Shape6[1]*Shape5[1])
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void PEPS<T>::" <<
           "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Direction, " <<
                              "unsigned int position, const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                              "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
           "(MPS0.d[" << pos << "] != PEPS0.D[" << pos << "]*PEPO0.D[" << pos << "]*this->D[" << pos << "])." << endl;
   exit(1);
  }
 }
#endif
// 0. - check MPS1.d:
 MPS1.getd(dMPS);
#ifdef DEBUG
 for (int pos = 0; pos < L; pos++)
 {
  BulkTensors0[pos].getShape(Shape05);
  BulkTensorsPEPO[pos].getShape(Shape6);
  BulkTensors[pos].getShape(Shape5);
  if (dMPS[pos] != Shape05[3]*Shape6[3]*Shape5[3])
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void PEPS<T>::" <<
           "multiplyBulkMPOMPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Direction, " <<
                              "unsigned int position, const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps, " <<
                              "double& epsAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) const: " <<
           "(MPS1.d[" << pos << "] != PEPS0.D[" << pos << "]*PEPO0.D[" << pos << "]*this->D[" << pos << "])." << endl;
   exit(1);
  }
 }
#endif
// 1. MPS1.D >= PEPS0.D*PEPO0.D*this->D*MPS0.D:
 DMPS1 = MPS1.getD();
 if (DMPS1 >= PEPS0.D*PEPO0.getD()*this->D*MPS0.getD())
 {
  MPS1 = MPS<T>("open", L, dMPS, PEPS0.D*PEPO0.getD()*this->D*MPS0.getD());
  Indices0[0] = 2; Indices1[0] = 1;
  Indices02[0] = 2; Indices12[0] = 1; Indices12[1] = 4;
  Order11[0] = 0; Order11[1] = 2; Order11[2] = 5; Order11[3] = 8; Order11[4] = 1; Order11[5] = 3;
  Order11[6] = 6; Order11[7] = 9; Order11[8] = 4; Order11[9] = 7; Order11[10] = 10;
  for (int pos = 0; pos < L; pos++)
  {
   MPS0.get(pos, Tensor0);
   Tensor0.getShape(Shape3);
   BulkTensors0[pos].getShape(Shape05);
   BulkTensorsPEPO[pos].getShape(Shape6);
   BulkTensors[pos].getShape(Shape5);
   Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
   Tensor0.reshape(Shape15);
   Tensor1 = BulkTensors0[pos];
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BulkTensorsPEPO[pos];
   Indices02[1] = 7;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   BulkTensors[pos].complexConjugate(Tensor1);
   Indices02[1] = 9;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.permute(Order11);
   Tensor0.getShape(Shape11);
   Shape3[0] = Shape11[0]*Shape11[1]*Shape11[2]*Shape11[3]; Shape3[1] = Shape11[4]*Shape11[5]*Shape11[6]*Shape11[7];
   Shape3[2] = Shape11[8]*Shape11[9]*Shape11[10];
   Tensor0.reshape(Shape3);
   MPS1.set(pos, Tensor0);
  }
  MPS1.bringIntoNormalShape();
  MPS1.setD(DMPS1);
  epsAchieved = 0.0;
  numSweepsDone = 0;
  return;
 }
// 2. bring MPS1 into normal form from right to left:
// 2.a) MPS1.D == 1:
 if (DMPS1 == 1)
 {
  tensorFrobeniusNormProduct = 1.0;
  for (int pos = L-1; pos > 0; pos--)
  {
   MPS1.get(pos, Tensor0);
   tensorFrobeniusNormProduct *= Tensor0.frobeniusNormalize();
   MPS1.set(pos, Tensor0);
  }
  MPS1.get(0, Tensor0);
  Tensor0.multiply(tensorFrobeniusNormProduct);
  MPS1.set(0, Tensor0);
 }
// 2.b) MPS1.D != 1:
 else
 {
  Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
  Indices0[0] = 1;
  Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
  MPS1.get(L-1, Tensor0);
  for (int pos = L-1; pos > 0; pos--)
  {
   Tensor0.LQDecompose(Indices1, Indices2, TensorL);
   MPS1.set(pos, Tensor0);
   MPS1.get(pos-1, Tensor0);
   Tensor0.contract(Indices0, TensorL, Indices1);
   Tensor0.permute(Order3);
  }
  MPS1.set(0, Tensor0);
 }
// 3. construct right half of b for each site:
 bTensors = vector< Tensor<T> >(L);
 MPS0.get(L-1, Tensor0);
 Tensor0.getShape(Shape3);
 BulkTensors0[L-1].getShape(Shape05);
 BulkTensorsPEPO[L-1].getShape(Shape6);
 BulkTensors[L-1].getShape(Shape5);
 Shape15[0] = Shape3[0]; Shape15[1] = 1; Shape15[2] = Shape05[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
 Tensor0.reshape(Shape15);
 Tensor1 = BulkTensors0[L-1];
 Indices0[0] = 2; Indices1[0] = 1;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = BulkTensorsPEPO[L-1];
 Indices02[0] = 2; Indices02[1] = 7; Indices12[0] = 1; Indices12[1] = 4;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 BulkTensors[L-1].complexConjugate(Tensor1);
 Indices02[1] = 9;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 MPS1.get(L-1, Tensor1);
 Tensor1.getShape(Shape3);
 Shape15[0] = Shape3[0]; Shape15[1] = 1; Shape15[2] = Shape05[3]; Shape15[3] = Shape6[3]; Shape15[4] = Shape5[3];
 Tensor1.reshape(Shape15);
 Tensor1.complexConjugate();
 Indices03[0] = 4; Indices03[1] = 7; Indices03[2] = 10; Indices13[0] = 2; Indices13[1] = 3; Indices13[2] = 4;
 Tensor0.contract(Indices03, Tensor1, Indices13);
 Tensor0.getShape(Shape10);
 Shape5[0] = Shape10[0]; Shape5[1] = Shape10[2]; Shape5[2] = Shape10[4]; Shape5[3] = Shape10[6]; Shape5[4] = Shape10[8];
 Tensor0.reshape(Shape5);
 bTensors[L-1] = Tensor0;
 Indices0[0] = 0; Indices1[0] = 1;
 Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 2; Indices12[1] = 1;
 Indices03[0] = 0; Indices03[2] = 8; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
 Indices04[0] = 0; Indices04[1] = 3; Indices04[2] = 5; Indices04[3] = 7; Indices14[0] = 1; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
 for (int pos = L-2; pos > 0; pos--)
 {
  MPS0.get(pos, Tensor1);
  Tensor1.getShape(Shape3);
  BulkTensors0[pos].getShape(Shape05);
  BulkTensorsPEPO[pos].getShape(Shape6);
  BulkTensors[pos].getShape(Shape5);
  Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
  Tensor1.reshape(Shape15);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor1 = BulkTensors0[pos];
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor1 = BulkTensorsPEPO[pos];
  Indices03[1] = 4;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  BulkTensors[pos].complexConjugate(Tensor1);
  Indices03[1] = 3;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  MPS1.get(pos, Tensor1);
  Tensor1.getShape(Shape3);
  Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[3]; Shape15[3] = Shape6[3]; Shape15[4] = Shape5[3];
  Tensor1.reshape(Shape15);
  Tensor1.complexConjugate();
  Tensor0.contract(Indices04, Tensor1, Indices14);
  bTensors[pos] = Tensor0;
 }
// 4. compute initial error:
 MPS0.get(0, Tensor1);
 Tensor1.getShape(Shape3);
 BulkTensors0[0].getShape(Shape05);
 BulkTensorsPEPO[0].getShape(Shape6);
 BulkTensors[0].getShape(Shape5);
 Shape15[0] = 1; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
 Tensor1.reshape(Shape15);
 Indices0[0] = 0; Indices1[0] = 1;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = BulkTensors0[0];
 Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 2; Indices12[1] = 1;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 Tensor1 = BulkTensorsPEPO[0];
 Indices03[0] = 0; Indices03[1] = 4; Indices03[2] = 8; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
 Tensor0.contract(Indices03, Tensor1, Indices13);
 BulkTensors[0].complexConjugate(Tensor1);
 Indices03[1] = 3;
 Tensor0.contract(Indices03, Tensor1, Indices13);
 MPS1.get(0, Tensor1);
 Tensor1.getShape(Shape3);
 Shape15[0] = 1; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[3]; Shape15[3] = Shape6[3]; Shape15[4] = Shape5[3];
 Tensor1.reshape(Shape15);
 Tensor1.complexConjugate();
 Indices04[0] = 0; Indices04[1] = 3; Indices04[2] = 5; Indices04[3] = 7; Indices14[0] = 1; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
 Tensor0.contract(Indices04, Tensor1, Indices14);
 bValue = Tensor0.get(0);
 MPS1.get(0, Tensor0);
 NValue = Tensor0.scalarProduct();
 error0 = MathAuxiliary::convertToDouble(NValue)-2.0*MathAuxiliary::convertToDouble(bValue);
// 5. sweep:
 numSweepsDone = 0;
 SweepDirection = "right";
 while (numSweepsDone < maxNumSweeps)
 {
// 5.a) sweep to the right:
  if (SweepDirection == "right")
  {
// 5.a)1. position 0:
// 5.a)1. - compute bTensor:
   Tensor0 = bTensors[1];
   MPS0.get(0, Tensor1);
   Tensor1.getShape(Shape3);
   BulkTensors0[0].getShape(Shape05);
   BulkTensorsPEPO[0].getShape(Shape6);
   BulkTensors[0].getShape(Shape5);
   Shape15[0] = 1; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
   Tensor1.reshape(Shape15);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BulkTensors0[0];
   Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 2; Indices12[1] = 1;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor1 = BulkTensorsPEPO[0];
   Indices03[0] = 0; Indices03[1] = 4; Indices03[2] = 8; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   BulkTensors[0].complexConjugate(Tensor1);
   Indices03[1] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
// 5.a)1. - update tensor:
   Tensor0.getShape(Shape8);
   Shape3[0] = 1; Shape3[1] = Shape8[0]; Shape3[2] = Shape8[3]*Shape8[5]*Shape8[7];
   Tensor0.reshape(Shape3);
   MPS1.set(0, Tensor0);
// 5.a)1. - compute error:
   NValue = Tensor0.scalarProduct();
   error1 = -MathAuxiliary::convertToDouble(NValue);
   epsAchieved = abs((error1-error0)/error1);
   if (epsAchieved < eps)
    return;
   error0 = error1;
// 5.a)1. - normalize tensor:
   if (DMPS1 == 1)
    tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
   else
   {
    Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
    Tensor0.QRDecompose(Indices2, Indices1, TensorR);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
   }
   MPS1.set(0, Tensor0);
// 5.a)1. - set bTensors:
   MPS0.get(0, Tensor0);
   Tensor0.getShape(Shape3);
   Shape15[0] = 1; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
   Tensor0.reshape(Shape15);
   Tensor1 = BulkTensors0[0];
   Indices0[0] = 2; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BulkTensorsPEPO[0];
   Indices02[0] = 2; Indices02[1] = 7; Indices12[0] = 1; Indices12[1] = 4;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   BulkTensors[0].complexConjugate(Tensor1);
   Indices02[1] = 9;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   MPS1.get(0, Tensor1);
   Tensor1.getShape(Shape3);
   Shape15[0] = 1; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[3]; Shape15[3] = Shape6[3]; Shape15[4] = Shape5[3];
   Tensor1.reshape(Shape15);
   Tensor1.complexConjugate();
   Indices03[0] = 4; Indices03[1] = 7; Indices03[2] = 10; Indices13[0] = 2; Indices13[1] = 3; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape10);
   Shape5[0] = Shape10[1]; Shape5[1] = Shape10[3]; Shape5[2] = Shape10[5]; Shape5[3] = Shape10[7]; Shape5[4] = Shape10[9];
   Tensor0.reshape(Shape5);
   bTensors[0] = Tensor0;
// 5.a)2. 0 < pos < L-1:
   for (int pos = 1; pos < L-1; pos++)
   {
// 5.a)2. - compute bTensor:
    Tensor0 = bTensors[pos-1];
    MPS0.get(pos, Tensor1);
    Tensor1.getShape(Shape3);
    BulkTensors0[pos].getShape(Shape05);
    BulkTensorsPEPO[pos].getShape(Shape6);
    BulkTensors[pos].getShape(Shape5);
    Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
    Tensor1.reshape(Shape15);
    Indices0[0] = 0; Indices1[0] = 0;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = BulkTensors0[pos];
    Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 1;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    Tensor1 = BulkTensorsPEPO[pos];
    Indices03[0] = 0; Indices03[1] = 4; Indices03[2] = 8; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 4;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    BulkTensors[pos].complexConjugate(Tensor1);
    Indices03[1] = 3;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    bTensor = Tensor0;
// 5.a)2. - update tensor:
    Tensor1 = bTensors[pos+1];
    Indices04[0] = 1; Indices04[1] = 2; Indices04[2] = 4; Indices04[3] = 6; Indices14[0] = 0; Indices14[1] = 1; Indices14[2] = 2; Indices14[3] = 3;
    Tensor0.contract(Indices04, Tensor1, Indices14);
    Tensor0.getShape(Shape15);
    Shape3[0] = Shape15[0]; Shape3[1] = Shape15[1]*Shape15[2]*Shape15[3]; Shape3[2] = Shape15[4];
    Tensor0.reshape(Shape3);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
    MPS1.set(pos, Tensor0);
// 5.a)2. - compute error:
    NValue = Tensor0.scalarProduct();
    error1 = -MathAuxiliary::convertToDouble(NValue);
    epsAchieved = abs((error1-error0)/error1);
    if (epsAchieved < eps)
     return;
    error0 = error1;
// 5.a)2. - normalize tensor:
    if (DMPS1 == 1)
     tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
    else
    {
     Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
     Tensor0.QRDecompose(Indices2, Indices1, TensorR);
     Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
     Tensor0.permute(Order3);
    }
    MPS1.set(pos, Tensor0);
// 5.a)2. - set bTensors:
    Tensor0.getShape(Shape3);
    Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[3]; Shape15[3] = Shape6[3]; Shape15[4] = Shape5[3];
    Tensor0.reshape(Shape15);
    Tensor0.complexConjugate();
    Indices04[0] = 0; Indices04[1] = 3; Indices04[2] = 5; Indices04[3] = 7; Indices14[0] = 0; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
    bTensor.contract(Indices04, Tensor0, Indices14);
    bTensors[pos] = bTensor;
   }
// 5.a)3. position L-1:
   MPS1.get(L-1, Tensor0);
   if (DMPS1 == 1)
    Tensor0.multiply(tensorFrobeniusNorm);
   else
   {
    Indices0[0] = 1; Indices1[0] = 0;
    TensorR.contract(Indices0, Tensor0, Indices1);
    Tensor0 = TensorR;
   }
   MPS1.set(L-1, Tensor0);
   SweepDirection = "left";
  }
// 5.b) sweep to the left:
  else if (SweepDirection == "left")
  {
// 5.b)1. position L-1:
// 5.b)1. - compute bTensor:
   Tensor0 = bTensors[L-2];
   MPS0.get(L-1, Tensor1);
   Tensor1.getShape(Shape3);
   BulkTensors0[L-1].getShape(Shape05);
   BulkTensorsPEPO[L-1].getShape(Shape6);
   BulkTensors[L-1].getShape(Shape5);
   Shape15[0] = Shape3[0]; Shape15[1] = 1; Shape15[2] = Shape05[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
   Tensor1.reshape(Shape15);
   Indices0[0] = 0; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BulkTensors0[L-1];
   Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 1;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor1 = BulkTensorsPEPO[L-1];
   Indices03[0] = 0; Indices03[1] = 4; Indices03[2] = 8; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   BulkTensors[L-1].complexConjugate(Tensor1);
   Indices03[1] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
// 5.b)1. - update tensor:
   Tensor0.getShape(Shape8);
   Shape3[0] = Shape8[0]; Shape3[1] = 1; Shape3[2] = Shape8[3]*Shape8[5]*Shape8[7];
   Tensor0.reshape(Shape3);
   MPS1.set(L-1, Tensor0);
// 5.b)1. - compute error:
   NValue = Tensor0.scalarProduct();
   error1 = -MathAuxiliary::convertToDouble(NValue);
   epsAchieved = abs((error1-error0)/error1);
   if (epsAchieved < eps)
    return;
   error0 = error1;
// 5.b)1. - normalize tensor:
   if (DMPS1 == 1)
    tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
   else
   {
    Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
    Tensor0.LQDecompose(Indices1, Indices2, TensorL);
   }
   MPS1.set(L-1, Tensor0);
// 5.b)1. - set bTensors:
   MPS0.get(L-1, Tensor0);
   Tensor0.getShape(Shape3);
   Shape15[0] = Shape3[0]; Shape15[1] = 1; Shape15[2] = Shape05[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
   Tensor0.reshape(Shape15);
   Tensor1 = BulkTensors0[L-1];
   Indices0[0] = 2; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = BulkTensorsPEPO[L-1];
   Indices02[0] = 2; Indices02[1] = 7; Indices12[0] = 1; Indices12[1] = 4;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   BulkTensors[L-1].complexConjugate(Tensor1);
   Indices02[1] = 9;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   MPS1.get(L-1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape15[0] = Shape3[0]; Shape15[1] = 1; Shape15[2] = Shape05[3]; Shape15[3] = Shape6[3]; Shape15[4] = Shape5[3];
   Tensor1.reshape(Shape15);
   Tensor1.complexConjugate();
   Indices03[0] = 4; Indices03[1] = 7; Indices03[2] = 10; Indices13[0] = 2; Indices13[1] = 3; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape10);
   Shape5[0] = Shape10[0]; Shape5[1] = Shape10[2]; Shape5[2] = Shape10[4]; Shape5[3] = Shape10[6]; Shape5[4] = Shape10[8];
   Tensor0.reshape(Shape5);
   bTensors[L-1] = Tensor0;
// 5.b)2. 0 < pos < L-1:
   for (int pos = L-2; pos > 0; pos--)
   {
// 5.b)2. - compute bTensor:
    Tensor0 = bTensors[pos+1];
    MPS0.get(pos, Tensor1);
    Tensor1.getShape(Shape3);
    BulkTensors0[pos].getShape(Shape05);
    BulkTensorsPEPO[pos].getShape(Shape6);
    BulkTensors[pos].getShape(Shape5);
    Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[1]; Shape15[3] = Shape6[1]; Shape15[4] = Shape5[1];
    Tensor1.reshape(Shape15);
    Indices0[0] = 0; Indices1[0] = 1;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Tensor1 = BulkTensors0[pos];
    Indices02[0] = 0; Indices02[1] = 5; Indices12[0] = 2; Indices12[1] = 1;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    Tensor1 = BulkTensorsPEPO[pos];
    Indices03[0] = 0; Indices03[1] = 4; Indices03[2] = 8; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    BulkTensors[pos].complexConjugate(Tensor1);
    Indices03[1] = 3;
    Tensor0.contract(Indices03, Tensor1, Indices13);
    bTensor = Tensor0;
// 5.b)2. - update tensor:
    Tensor1 = bTensors[pos-1];
    Indices04[0] = 1; Indices04[1] = 2; Indices04[2] = 4; Indices04[3] = 6; Indices14[0] = 0; Indices14[1] = 1; Indices14[2] = 2; Indices14[3] = 3;
    Tensor0.contract(Indices04, Tensor1, Indices14);
    Tensor0.getShape(Shape15);
    Shape3[0] = Shape15[0]; Shape3[1] = Shape15[1]*Shape15[2]*Shape15[3]; Shape3[2] = Shape15[4];
    Tensor0.reshape(Shape3);
    Order3[0] = 2; Order3[1] = 0; Order3[2] = 1;
    Tensor0.permute(Order3);
    MPS1.set(pos, Tensor0);
// 5.b)2. - compute error:
    NValue = Tensor0.scalarProduct();
    error1 = -MathAuxiliary::convertToDouble(NValue);
    epsAchieved = abs((error1-error0)/error1);
    if (epsAchieved < eps)
     return;
    error0 = error1;
// 5.b)2. - normalize tensor:
    if (DMPS1 == 1)
     tensorFrobeniusNorm = Tensor0.frobeniusNormalize();
    else
    {
     Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
     Tensor0.LQDecompose(Indices1, Indices2, TensorL);
    }
    MPS1.set(pos, Tensor0);
// 5.b)2. - set bTensors:
    Tensor0.getShape(Shape3);
    Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[3]; Shape15[3] = Shape6[3]; Shape15[4] = Shape5[3];
    Tensor0.reshape(Shape15);
    Tensor0.complexConjugate();
    Indices04[0] = 0; Indices04[1] = 3; Indices04[2] = 5; Indices04[3] = 7; Indices14[0] = 1; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
    bTensor.contract(Indices04, Tensor0, Indices14);
    bTensors[pos] = bTensor;
   }
// 5.b)3. position 0:
   MPS1.get(0, Tensor0);
   if (DMPS1 == 1)
    Tensor0.multiply(tensorFrobeniusNorm);
   else
   {
    Indices0[0] = 1; Indices1[0] = 0;
    Tensor0.contract(Indices0, TensorL, Indices1);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0.permute(Order3);
   }
   MPS1.set(0, Tensor0);
   SweepDirection = "right";
  }
  numSweepsDone++;
 }
}

template<class T> T PEPS<T>::scalarProduct(const string& Direction, const vector<unsigned int>& D2s,
                                           const vector<double>& Epss,
                                           const vector<unsigned int>& MaxNumsSweeps,
                                           vector<double>& ErrorsAchieved,
                                           vector<unsigned int>& NumsSweepsDone) const
{
#ifdef DEBUG
 if ((this->N[0] == 0) || ((Direction != "horizontal") && (Direction != "vertical")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
          "((this->N[0] == 0) || ((Direction != horizontal) && (Direction != vertical)))." << endl;
  exit(1);
 }
 if (((Direction == "horizontal") && (D2s.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (Epss.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (MaxNumsSweeps.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (ErrorsAchieved.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (NumsSweepsDone.size() != this->N[1]-3)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
          "(((Direction == horizontal) && (D2s.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (Epss.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (MaxNumsSweeps.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (ErrorsAchieved.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (NumsSweepsDone.size() != this->N[1]-3)))." << endl;
  exit(1);
 }
 if (((Direction == "vertical") && (D2s.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (Epss.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (MaxNumsSweeps.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (ErrorsAchieved.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (NumsSweepsDone.size() != this->N[0]-3)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
          "(((Direction == vertical) && (D2s.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (Epss.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (MaxNumsSweeps.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (ErrorsAchieved.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (NumsSweepsDone.size() != this->N[0]-3)))." << endl;
  exit(1);
 }
 for (int i = 0; i < D2s.size(); i++)
 {
  if (D2s[i] == 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> T PEPS<T>::" <<
           "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                         "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                         "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
           "(D2s[" << i << "] == 0)." << endl;
   exit(1);
  }
 }
#endif
 double eps, errorAchieved;
 unsigned int maxNumSweeps, numSweepsDone;
 MPS<T> MPS0, MPS1, MPS2; MPO<T> MPO0;
 unsigned int D1;
 if (this->BC == "open")
 {
  if (Direction == "horizontal")
  {
// compute MPS1 approximating the left half of the PEPS sandwich:
   this->getMPS("right", MPS1);
   for (int col = 1; col < this->N[1]/2; col++)
   {
    MPS0 = MPS1;
    this->getMPO("right", col, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[col-1]);
    MPS1.setD(D1);
    eps = Epss[col-1];
    maxNumSweeps = MaxNumsSweeps[col-1];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
    ErrorsAchieved[col-1] = errorAchieved;
    NumsSweepsDone[col-1] = numSweepsDone;
   }
// compute MPS2 approximating the right half of the PEPS sandwich:
   this->getMPS("left", MPS2);
   for (int col = this->N[1]-2; col > this->N[1]/2; col--)
   {
    MPS0 = MPS2;
    this->getMPO("left", col, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[col-2]);
    MPS2.setD(D1);
    eps = Epss[col-2];
    maxNumSweeps = MaxNumsSweeps[col-2];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS2);
    ErrorsAchieved[col-2] = errorAchieved;
    NumsSweepsDone[col-2] = numSweepsDone;
   }
// return the contraction value:
   if (this->N[1] == 2)
   {
    return MPS2.contractReverse(MPS1);
   }
   else
   {
    this->getMPO("right", this->N[1]/2, MPO0);
    return contractReverseMPOMPS(MPS2, MPO0, MPS1);
   }
  }
  else if (Direction == "vertical")
  {
// compute MPS1 approximating the upper half of the PEPS sandwich:
   this->getMPS("down", MPS1);
   for (int row = 1; row < this->N[0]/2; row++)
   {
    MPS0 = MPS1;
    this->getMPO("down", row, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[row-1]);
    MPS1.setD(D1);
    eps = Epss[row-1];
    maxNumSweeps = MaxNumsSweeps[row-1];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
    ErrorsAchieved[row-1] = errorAchieved;
    NumsSweepsDone[row-1] = numSweepsDone;
   }
// compute MPS2 approximating the lower half of the PEPS sandwich:
   this->getMPS("up", MPS2);
   for (int row = this->N[0]-2; row > this->N[0]/2; row--)
   {
    MPS0 = MPS2;
    this->getMPO("up", row, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[row-2]);
    MPS2.setD(D1);
    eps = Epss[row-2];
    maxNumSweeps = MaxNumsSweeps[row-2];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS2);
    ErrorsAchieved[row-2] = errorAchieved;
    NumsSweepsDone[row-2] = numSweepsDone;
   }
// return the contraction value:
   if (this->N[0] == 2)
   {
    return MPS2.contractReverse(MPS1);
   }
   else
   {
    this->getMPO("down", this->N[0]/2, MPO0);
    return contractReverseMPOMPS(MPS2, MPO0, MPS1);
   }
  }
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented for periodic BC: " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, " <<
                        "vector<unsigned int>& NumsSweepsDone) const." << endl;
  exit(1);
 }
}

template<class T> T PEPS<T>::scalarProduct(const string& Direction, const vector<unsigned int>& D2s,
                                           const vector<double>& Epss,
                                           const vector<unsigned int>& MaxNumsSweeps,
                                           vector<double>& ErrorsAchieved,
                                           vector<unsigned int>& NumsSweepsDone, vector<T>& Norms) const
{
#ifdef DEBUG
 if ((this->N[0] == 0) || ((Direction != "horizontal") && (Direction != "vertical")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "((this->N[0] == 0) || ((Direction != horizontal) && (Direction != vertical)))." << endl;
  exit(1);
 }
 if (((Direction == "horizontal") && (D2s.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (Epss.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (MaxNumsSweeps.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (ErrorsAchieved.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (NumsSweepsDone.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (Norms.size() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "(((Direction == horizontal) && (D2s.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (Epss.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (MaxNumsSweeps.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (ErrorsAchieved.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (NumsSweepsDone.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (Norms.size() != this->N[1])))." << endl;
  exit(1);
 }
 if (((Direction == "vertical") && (D2s.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (Epss.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (MaxNumsSweeps.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (ErrorsAchieved.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (NumsSweepsDone.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (Norms.size() != this->N[0])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "(((Direction == vertical) && (D2s.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (Epss.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (MaxNumsSweeps.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (ErrorsAchieved.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (NumsSweepsDone.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (Norms.size() != this->N[0])))." << endl;
  exit(1);
 }
 for (int i = 0; i < D2s.size(); i++)
 {
  if (D2s[i] == 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> T PEPS<T>::" <<
           "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                         "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                         "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                         "vector<T>& Norms) const: " <<
           "(D2s[" << i << "] == 0)." << endl;
   exit(1);
  }
 }
#endif
 double eps, errorAchieved;
 unsigned int maxNumSweeps, numSweepsDone;
 MPS<T> MPS0, MPS1, MPS2; MPO<T> MPO0;
 unsigned int D1;
 if (this->BC == "open")
 {
  if (Direction == "horizontal")
  {
// compute MPS1 approximating the left half of the PEPS sandwich:
   this->getMPS("right", MPS1);
   Norms[0] = MPS1.simplifiedNormalize();
   for (int col = 1; col < this->N[1]/2; col++)
   {
    MPS0 = MPS1;
    this->getMPO("right", col, MPO0);
    Norms[col] = MPO0.simplifiedNormalize();
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[col-1]);
    MPS1.setD(D1);
    eps = Epss[col-1];
    maxNumSweeps = MaxNumsSweeps[col-1];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
    ErrorsAchieved[col-1] = errorAchieved;
    NumsSweepsDone[col-1] = numSweepsDone;
   }
// compute MPS2 approximating the right half of the PEPS sandwich:
   this->getMPS("left", MPS2);
   Norms[this->N[1]-1] = MPS2.simplifiedNormalize();
   for (int col = this->N[1]-2; col > this->N[1]/2; col--)
   {
    MPS0 = MPS2;
    this->getMPO("left", col, MPO0);
    Norms[col] = MPO0.simplifiedNormalize();
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[col-2]);
    MPS2.setD(D1);
    eps = Epss[col-2];
    maxNumSweeps = MaxNumsSweeps[col-2];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS2);
    ErrorsAchieved[col-2] = errorAchieved;
    NumsSweepsDone[col-2] = numSweepsDone;
   }
// return the contraction value:
   if (this->N[1] == 2)
   {
    return MPS2.contractReverse(MPS1);
   }
   else
   {
    this->getMPO("right", this->N[1]/2, MPO0);
    Norms[this->N[1]/2] = MPO0.simplifiedNormalize();
    return contractReverseMPOMPS(MPS2, MPO0, MPS1);
   }
  }
  else if (Direction == "vertical")
  {
// compute MPS1 approximating the upper half of the PEPS sandwich:
   this->getMPS("down", MPS1);
   Norms[0] = MPS1.simplifiedNormalize();
   for (int row = 1; row < this->N[0]/2; row++)
   {
    MPS0 = MPS1;
    this->getMPO("down", row, MPO0);
    Norms[row] = MPO0.simplifiedNormalize();
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[row-1]);
    MPS1.setD(D1);
    eps = Epss[row-1];
    maxNumSweeps = MaxNumsSweeps[row-1];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
    ErrorsAchieved[row-1] = errorAchieved;
    NumsSweepsDone[row-1] = numSweepsDone;
   }
// compute MPS2 approximating the lower half of the PEPS sandwich:
   this->getMPS("up", MPS2);
   Norms[this->N[0]-1] = MPS2.simplifiedNormalize();
   for (int row = this->N[0]-2; row > this->N[0]/2; row--)
   {
    MPS0 = MPS2;
    this->getMPO("up", row, MPO0);
    Norms[row] = MPO0.simplifiedNormalize();
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[row-2]);
    MPS2.setD(D1);
    eps = Epss[row-2];
    maxNumSweeps = MaxNumsSweeps[row-2];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS2);
    ErrorsAchieved[row-2] = errorAchieved;
    NumsSweepsDone[row-2] = numSweepsDone;
   }
// return the contraction value:
   if (this->N[0] == 2)
   {
    return MPS2.contractReverse(MPS1);
   }
   else
   {
    this->getMPO("down", this->N[0]/2, MPO0);
    Norms[this->N[0]/2] = MPO0.simplifiedNormalize();
    return contractReverseMPOMPS(MPS2, MPO0, MPS1);
   }
  }
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented for periodic BC: " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, " <<
                        "vector<unsigned int>& NumsSweepsDone, vector<T>& Norms) const." << endl;
  exit(1);
 }
}

template<class T> T PEPS<T>::scalarProduct(const string& Direction,
                                           unsigned int D2, double eps, unsigned int maxNumSweeps,
                                           vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone,
                                           vector<T>& Norms) const
{
#ifdef DEBUG
 if ((this->BC != "open") || (this->N[0] == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, " <<
                        "unsigned int D2, double eps, unsigned int maxNumSweeps, " <<
                        "vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "((this->BC != open) || (this->N[0] == 0))." << endl;
 }
 if ((Direction != "horizontal") && (Direction != "vertical"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, " <<
                        "unsigned int D2, double eps, unsigned int maxNumSweeps, " <<
                        "vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "((Direction != horizontal) && (Direction != vertical))." << endl;
 }
 if ((D2 == 0) || (maxNumSweeps == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, " <<
                        "unsigned int D2, double eps, unsigned int maxNumSweeps, " <<
                        "vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "((D2 == 0) || (maxNumSweeps == 0))." << endl;
 }
 if (((Direction == "horizontal") && (EpssAchieved.size() != this->N[1]-1)) ||
     ((Direction == "horizontal") && (NumsSweepsDone.size() != this->N[1]-1)) ||
     ((Direction == "horizontal") && (Norms.size() != this->N[1]-1)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, " <<
                        "unsigned int D2, double eps, unsigned int maxNumSweeps, " <<
                        "vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "(((Direction == horizontal) && (EpssAchieved.size() != this->N[1]-1)) || " <<
           "((Direction == horizontal) && (NumsSweepsDone.size() != this->N[1]-1)) || " <<
           "((Direction == horizontal) && (Norms.size() != this->N[1]-1)))." << endl;
  exit(1);
 }
 if (((Direction == "vertical") && (EpssAchieved.size() != this->N[0]-1)) ||
     ((Direction == "vertical") && (NumsSweepsDone.size() != this->N[0]-1)) ||
     ((Direction == "vertical") && (Norms.size() != this->N[0]-1)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, " <<
                        "unsigned int D2, double eps, unsigned int maxNumSweeps, " <<
                        "vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "(((Direction == vertical) && (EpssAchieved.size() != this->N[0]-1)) || " <<
           "((Direction == vertical) && (NumsSweepsDone.size() != this->N[0]-1)) || " <<
           "((Direction == vertical) && (Norms.size() != this->N[0]-1)))." << endl;
  exit(1);
 }
#endif
 unsigned int L, Ncols = this->N[1], Nrows = this->N[0];
 T result;
 string BCMPS = "open", Boundary, SweepDirection;
 vector<unsigned int> dMPS, dMPS1, Seed, Shape(5);
 MPS<T> MPS0, MPS1, MPS2;
// a) horizontal contraction:
 if (Direction == "horizontal")
 {
  L = Nrows;
  dMPS = vector<unsigned int>(L);
  dMPS1 = vector<unsigned int>(L);
  Seed = vector<unsigned int>(L);
// a)1. compute MPS1 approximating the left half of the sandwich:
// a)1.1. boundary-MPS for left boundary:
  for (int i = 0; i < L; i++)
  {
   this->Tensors[L-1-i].getShape(Shape);
   dMPS[i] = Shape[2]*Shape[2];
  }
  MPS1 = MPS<T>(BCMPS, L, dMPS, min(D2, this->D*this->D));
  for (int i = 0; i < L; i++)
   Seed[i] = time(0)+13*i;
  MPS1.fillRandomly(Seed);
  Boundary = "left";
  this->getBoundaryMPS(Boundary, eps, maxNumSweeps, EpssAchieved[0], NumsSweepsDone[0], MPS1);
  Norms[0] = MPS1.normalize();
// a)1.2. boundary-MPSs for bulk:
  SweepDirection = "right";
  for (int pos = 1; pos < Ncols/2; pos++)
  {
   MPS0 = MPS1;
   for (int i = 0; i < L; i++)
   {
    this->Tensors[L-1-i+Nrows*pos].getShape(Shape);
    dMPS1[i] = Shape[2]*Shape[2];
   }
   if (dMPS1 == dMPS)
   {
    MPS1.setD(min(D2, this->D*this->D*MPS0.getD()));
   }
   else
   {
    dMPS = dMPS1;
    MPS1 = MPS<T>(BCMPS, L, dMPS, min(D2, this->D*this->D*MPS0.getD()));
    for (int i = 0; i < L; i++)
     Seed[i] = time(0)+13*(i+Nrows*pos);
    MPS1.fillRandomly(Seed);
   }
   this->multiplyBulkMPOMPS(SweepDirection, pos, MPS0, eps, maxNumSweeps, EpssAchieved[pos], NumsSweepsDone[pos], MPS1);
   Norms[pos] = MPS1.normalize();
  }
// a)2. compute MPS2 approximating the right half of the sandwich:
// a)2.1. boundary-MPS for right boundary:
  for (int i = 0; i < L; i++)
  {
   this->Tensors[i+Nrows*(Ncols-1)].getShape(Shape);
   dMPS[i] = Shape[0]*Shape[0];
  }
  MPS2 = MPS<T>(BCMPS, L, dMPS, min(D2, this->D*this->D));
  for (int i = 0; i < L; i++)
   Seed[i] = time(0)+13*(i+Nrows*Ncols/2);
  MPS2.fillRandomly(Seed);
  Boundary = "right";
  this->getBoundaryMPS(Boundary, eps, maxNumSweeps, EpssAchieved[Ncols-2], NumsSweepsDone[Ncols-2], MPS2);
  Norms[Ncols-2] = MPS2.normalize();
// a)2.2. boundary-MPSs for bulk:
  SweepDirection = "left";
  for (int pos = Ncols-2; pos > Ncols/2; pos--)
  {
   MPS0 = MPS2;
   for (int i = 0; i < L; i++)
   {
    this->Tensors[i+Nrows*pos].getShape(Shape);
    dMPS1[i] = Shape[0]*Shape[0];
   }
   if (dMPS1 == dMPS)
   {
    MPS2.setD(min(D2, this->D*this->D*MPS0.getD()));
   }
   else
   {
    dMPS = dMPS1;
    MPS2 = MPS<T>(BCMPS, L, dMPS, min(D2, this->D*this->D*MPS0.getD()));
    for (int i = 0; i < L; i++)
     Seed[i] = time(0)+13*(i+Nrows*(Ncols/2+Ncols-1-pos));
    MPS2.fillRandomly(Seed);
   }
   this->multiplyBulkMPOMPS(SweepDirection, pos, MPS0, eps, maxNumSweeps, EpssAchieved[pos-1], NumsSweepsDone[pos-1], MPS2);
   Norms[pos-1] = MPS2.normalize();
  }
// a)3. perform final contraction:
  result = contractReverseBulkMPOMPS(Direction, MPS1, Ncols/2, MPS2);
 }
// b) vertical contraction:
 else if (Direction == "vertical")
 {
  L = Ncols;
  dMPS = vector<unsigned int>(L);
  dMPS1 = vector<unsigned int>(L);
  Seed = vector<unsigned int>(L);
// b)1. compute MPS1 approximating the upper half of the sandwich:
// b)1.1. boundary-MPS for top boundary:
  for (int i = 0; i < L; i++)
  {
   this->Tensors[Nrows*i].getShape(Shape);
   dMPS[i] = Shape[3]*Shape[3];
  }
  MPS1 = MPS<T>(BCMPS, L, dMPS, min(D2, this->D*this->D));
  for (int i = 0; i < L; i++)
   Seed[i] = time(0)+13*i;
  MPS1.fillRandomly(Seed);
  Boundary = "top";
  this->getBoundaryMPS(Boundary, eps, maxNumSweeps, EpssAchieved[0], NumsSweepsDone[0], MPS1);
  Norms[0] = MPS1.normalize();
// b)1.2. boundary-MPSs for bulk:
  SweepDirection = "down";
  for (int pos = 1; pos < Nrows/2; pos++)
  {
   MPS0 = MPS1;
   for (int i = 0; i < L; i++)
   {
    this->Tensors[pos+Nrows*i].getShape(Shape);
    dMPS1[i] = Shape[3]*Shape[3];
   }
   if (dMPS1 == dMPS)
   {
    MPS1.setD(min(D2, this->D*this->D*MPS0.getD()));
   }
   else
   {
    dMPS = dMPS1;
    MPS1 = MPS<T>(BCMPS, L, dMPS, min(D2, this->D*this->D*MPS0.getD()));
    for (int i = 0; i < L; i++)
     Seed[i] = time(0)+13*(i+Ncols*pos);
    MPS1.fillRandomly(Seed);
   }
   this->multiplyBulkMPOMPS(SweepDirection, pos, MPS0, eps, maxNumSweeps, EpssAchieved[pos], NumsSweepsDone[pos], MPS1);
   Norms[pos] = MPS1.normalize();
  }
// b)2. compute MPS2 approximating the lower half of the sandwich:
// b)2.1. boundary-MPS for bottom boundary:
  for (int i = 0; i < L; i++)
  {
   this->Tensors[Nrows-1+Nrows*(L-1-i)].getShape(Shape);
   dMPS[i] = Shape[1]*Shape[1];
  }
  MPS2 = MPS<T>(BCMPS, L, dMPS, min(D2, this->D*this->D));
  for (int i = 0; i < L; i++)
   Seed[i] = time(0)+13*(i+Ncols*Nrows/2);
  MPS2.fillRandomly(Seed);
  Boundary = "bottom";
  this->getBoundaryMPS(Boundary, eps, maxNumSweeps, EpssAchieved[Nrows-2], NumsSweepsDone[Nrows-2], MPS2);
  Norms[Nrows-2] = MPS2.normalize();
// b)2.2. boundary-MPSs for bulk:
  SweepDirection = "up";
  for (int pos = Nrows-2; pos > Nrows/2; pos--)
  {
   MPS0 = MPS2;
   for (int i = 0; i < L; i++)
   {
    this->Tensors[pos+Nrows*(L-1-i)].getShape(Shape);
    dMPS1[i] = Shape[1]*Shape[1];
   }
   if (dMPS1 == dMPS)
   {
    MPS2.setD(min(D2, this->D*this->D*MPS0.getD()));
   }
   else
   {
    dMPS = dMPS1;
    MPS2 = MPS<T>(BCMPS, L, dMPS, min(D2, this->D*this->D*MPS0.getD()));
    for (int i = 0; i < L; i++)
     Seed[i] = time(0)+13*(i+Ncols*(Nrows/2+Nrows-1-pos));
    MPS2.fillRandomly(Seed);
   }
   this->multiplyBulkMPOMPS(SweepDirection, pos, MPS0, eps, maxNumSweeps, EpssAchieved[pos-1], NumsSweepsDone[pos-1], MPS2);
   Norms[pos-1] = MPS2.normalize();
  }
// b)3. perform final contraction:
  result = contractReverseBulkMPOMPS(Direction, MPS2, Nrows/2, MPS1);
 }
 return result;
}

template<class T> T PEPS<T>::scalarProduct(const string& Direction, const vector<unsigned int>& D2s,
                                           const vector<double>& Epss,
                                           const vector<unsigned int>& MaxNumsSweeps,
                                           vector<double>& ErrorsAchieved,
                                           vector<unsigned int>& NumsSweepsDone,
                                           const vector<unsigned int>& d2s, double eps,
                                           Matrix< vector<double> >& Eigenvalues) const
{
#ifdef DEBUG
 if ((this->N[0] == 0) || ((Direction != "horizontal") && (Direction != "vertical")) || (eps < 0.0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "const vector<unsigned int>& d2s, double eps, " <<
                        "Matrix< vector<double> >& Eigenvalues) const: " <<
          "((this->N[0] == 0) || ((Direction != horizontal) && (Direction != vertical)) || (eps < 0.0))." << endl;
  exit(1);
 }
 if (((Direction == "horizontal") && (D2s.size() != this->N[1]-2)) ||
     ((Direction == "horizontal") && (Epss.size() != this->N[1]-2)) ||
     ((Direction == "horizontal") && (MaxNumsSweeps.size() != this->N[1]-2)) ||
     ((Direction == "horizontal") && (ErrorsAchieved.size() != this->N[1]-2)) ||
     ((Direction == "horizontal") && (NumsSweepsDone.size() != this->N[1]-2)) ||
     ((Direction == "horizontal") && (d2s.size() != this->N[1]-2)) ||
     ((Direction == "horizontal") && (Eigenvalues.getDim0() != this->N[0])) ||
     ((Direction == "horizontal") && (Eigenvalues.getDim1() != this->N[1]-2)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "const vector<unsigned int>& d2s, double eps, " <<
                        "Matrix< vector<double> >& Eigenvalues) const: " <<
          "(((Direction == horizontal) && (D2s.size() != this->N[1]-2)) || " <<
           "((Direction == horizontal) && (Epss.size() != this->N[1]-2)) || " <<
           "((Direction == horizontal) && (MaxNumsSweeps.size() != this->N[1]-2)) || " <<
           "((Direction == horizontal) && (ErrorsAchieved.size() != this->N[1]-2)) || " <<
           "((Direction == horizontal) && (NumsSweepsDone.size() != this->N[1]-2)) || " <<
           "((Direction == horizontal) && (d2s.size() != this->N[1]-2)) || " <<
           "((Direction == horizontal) && (Eigenvalues.getDim0() != this->N[0])) || " <<
           "((Direction == horizontal) && (Eigenvalues.getDim1() != this->N[1]-2)))." << endl;
  exit(1);
 }
 if (((Direction == "vertical") && (D2s.size() != this->N[0]-2)) ||
     ((Direction == "vertical") && (Epss.size() != this->N[0]-2)) ||
     ((Direction == "vertical") && (MaxNumsSweeps.size() != this->N[0]-2)) ||
     ((Direction == "vertical") && (ErrorsAchieved.size() != this->N[0]-2)) ||
     ((Direction == "vertical") && (NumsSweepsDone.size() != this->N[0]-2)) ||
     ((Direction == "vertical") && (d2s.size() != this->N[0]-2)) ||
     ((Direction == "vertical") && (Eigenvalues.getDim0() != this->N[0]-2)) ||
     ((Direction == "vertical") && (Eigenvalues.getDim1() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "const vector<unsigned int>& d2s, double eps, " <<
                        "Matrix< vector<double> >& Eigenvalues) const: " <<
          "(((Direction == vertical) && (D2s.size() != this->N[0]-2)) || " <<
           "((Direction == vertical) && (Epss.size() != this->N[0]-2)) || " <<
           "((Direction == vertical) && (MaxNumsSweeps.size() != this->N[0]-2)) || " <<
           "((Direction == vertical) && (ErrorsAchieved.size() != this->N[0]-2)) || " <<
           "((Direction == vertical) && (NumsSweepsDone.size() != this->N[0]-2)) || " <<
           "((Direction == vertical) && (d2s.size() != this->N[0]-2)) || " <<
           "((Direction == vertical) && (Eigenvalues.getDim0() != this->N[0]-2)) || " <<
           "((Direction == vertical) && (Eigenvalues.getDim1() != this->N[1])))." << endl;
  exit(1);
 }
 for (int i = 0; i < D2s.size(); i++)
 {
  if ((D2s[i] == 0) || (d2s[i] == 0))
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> T PEPS<T>::" <<
           "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                         "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                         "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                         "const vector<unsigned int>& d2s, double eps, " <<
                         "Matrix< vector<double> >& Eigenvalues) const: " <<
           "((D2s[" << i << "] == 0) || (d2s[" << i << "] == 0))." << endl;
   exit(1);
  }
 }
#endif
 unsigned int row, col, position;
 unsigned int maxNumSweeps, numSweepsDone, D1, d0;
 double eps0, errorAchieved;
 vector<unsigned int> Shape4(4), Shape5(5), Order(5);
 vector< vector<double> > Eigenvalues0;
 Tensor<T> Tensor0;
 MPS<T> MPS0, MPS1; MPO<T> MPO0, MPO1, MPO2;
 Order[0] = 0; Order[1] = 1; Order[2] = 3; Order[3] = 2; Order[4] = 4;
 if (this->BC == "open")
 {
  if (Direction == "horizontal")
  {
// compute MPS0 approximating the left half of the PEPS sandwich:
   this->getMPOSL("right", MPO2);
   for (col = 1; col < (this->N[1]+1)/2; col++)
   {
    MPO0 = MPO2;
    this->getMPOSL("right", col, MPO1);
    D1 = min(MPO1.getD()*MPO0.getD(), D2s[col-1]);
    MPO2.setD(D1);
    for (position = 0; position < MPO2.getN(); position++)
    {
     MPO2.get(position, Tensor0);
     Tensor0.getShape(Shape4);
     d0 = this->d(this->N[0]-1-position, col);
     Shape4[3] = this->D*d0;
     Tensor0 = Tensor<T>(Shape4);
     Tensor0.fillRandomly();
     MPO2.set(position, Tensor0);
    }
    eps0 = Epss[col-1];
    maxNumSweeps = MaxNumsSweeps[col-1];
    MPO1.multiply(MPO0, eps0, maxNumSweeps, errorAchieved, numSweepsDone, MPO2);
    ErrorsAchieved[col-1] = errorAchieved;
    NumsSweepsDone[col-1] = numSweepsDone;
    for (position = 0; position < MPO2.getN(); position++)
    {
     MPO2.get(position, Tensor0);
     Tensor0.getShape(Shape4);
     d0 = this->d(this->N[0]-1-position, col);
     Shape5[0] = Shape4[0]; Shape5[1] = Shape4[1]; Shape5[2] = Shape4[2]; Shape5[3] = Shape4[3]/d0;
     Shape5[4] = d0;
     Tensor0.reshape(Shape5);
     Tensor0.permute(Order);
     Shape4[2] = Shape5[3]; Shape4[3] = Shape5[2]*d0;
     Tensor0.reshape(Shape4);
     MPO2.setd(position, Shape5[2]*d0);
     MPO2.set(position, Tensor0);
    }
    MPO2.setd2(d2s[col-1], eps, Eigenvalues0);
    MPO2.transpose();
    for (position = 0; position < MPO2.getN(); position++)
     MPO2.setd(position, this->D);
    for (position = 0; position < MPO2.getN(); position++)
     Eigenvalues(this->N[0]-1-position, col-1) = Eigenvalues0[position];
   }
   MPO2.adjoint(MPO1);
   MPO2.multiply(MPO1, MPO0);
   MPO0.getMPS(MPS0);
// compute MPS1 approximating the right half of the PEPS sandwich:
   this->getMPOSL("left", MPO2);
   for (col = this->N[1]-2; col >= (this->N[1]+1)/2; col--)
   {
    MPO0 = MPO2;
    this->getMPOSL("left", col, MPO1);
    D1 = min(MPO1.getD()*MPO0.getD(), D2s[col-1]);
    MPO2.setD(D1);
    for (position = 0; position < MPO2.getN(); position++)
    {
     MPO2.get(position, Tensor0);
     Tensor0.getShape(Shape4);
     d0 = this->d(position, col);
     Shape4[3] = this->D*d0;
     Tensor0 = Tensor<T>(Shape4);
     Tensor0.fillRandomly();
     MPO2.set(position, Tensor0);
    }
    eps0 = Epss[col-1];
    maxNumSweeps = MaxNumsSweeps[col-1];
    MPO1.multiply(MPO0, eps0, maxNumSweeps, errorAchieved, numSweepsDone, MPO2);
    ErrorsAchieved[col-1] = errorAchieved;
    NumsSweepsDone[col-1] = numSweepsDone;
    for (position = 0; position < MPO2.getN(); position++)
    {
     MPO2.get(position, Tensor0);
     Tensor0.getShape(Shape4);
     d0 = this->d(position, col);
     Shape5[0] = Shape4[0]; Shape5[1] = Shape4[1]; Shape5[2] = Shape4[2]; Shape5[3] = Shape4[3]/d0;
     Shape5[4] = d0;
     Tensor0.reshape(Shape5);
     Tensor0.permute(Order);
     Shape4[2] = Shape5[3]; Shape4[3] = Shape5[2]*d0;
     Tensor0.reshape(Shape4);
     MPO2.setd(position, Shape5[2]*d0);
     MPO2.set(position, Tensor0);
    }
    MPO2.setd2(d2s[col-1], eps, Eigenvalues0);
    MPO2.transpose();
    for (position = 0; position < MPO2.getN(); position++)
     MPO2.setd(position, this->D);
    for (position = 0; position < MPO2.getN(); position++)
     Eigenvalues(position, col-1) = Eigenvalues0[position];
   }
   MPO2.adjoint(MPO1);
   MPO2.multiply(MPO1, MPO0);
   MPO0.getMPS(MPS1);
  }
  else if (Direction == "vertical")
  {
// compute MPS0 approximating the upper half of the PEPS sandwich:
   this->getMPOSL("down", MPO2);
   for (row = 1; row < (this->N[0]+1)/2; row++)
   {
    MPO0 = MPO2;
    this->getMPOSL("down", row, MPO1);
    D1 = min(MPO1.getD()*MPO0.getD(), D2s[row-1]);
    MPO2.setD(D1);
    for (position = 0; position < MPO2.getN(); position++)
    {
     MPO2.get(position, Tensor0);
     Tensor0.getShape(Shape4);
     d0 = this->d(row, position);
     Shape4[3] = this->D*d0;
     Tensor0 = Tensor<T>(Shape4);
     Tensor0.fillRandomly();
     MPO2.set(position, Tensor0);
    }
    eps0 = Epss[row-1];
    maxNumSweeps = MaxNumsSweeps[row-1];
    MPO1.multiply(MPO0, eps0, maxNumSweeps, errorAchieved, numSweepsDone, MPO2);
    ErrorsAchieved[row-1] = errorAchieved;
    NumsSweepsDone[row-1] = numSweepsDone;
    for (position = 0; position < MPO2.getN(); position++)
    {
     MPO2.get(position, Tensor0);
     Tensor0.getShape(Shape4);
     d0 = this->d(row, position);
     Shape5[0] = Shape4[0]; Shape5[1] = Shape4[1]; Shape5[2] = Shape4[2]; Shape5[3] = Shape4[3]/d0;
     Shape5[4] = d0;
     Tensor0.reshape(Shape5);
     Tensor0.permute(Order);
     Shape4[2] = Shape5[3]; Shape4[3] = Shape5[2]*d0;
     Tensor0.reshape(Shape4);
     MPO2.setd(position, Shape5[2]*d0);
     MPO2.set(position, Tensor0);
    }
    MPO2.setd2(d2s[row-1], eps, Eigenvalues0);
    MPO2.transpose();
    for (position = 0; position < MPO2.getN(); position++)
     MPO2.setd(position, this->D);
    for (position = 0; position < MPO2.getN(); position++)
     Eigenvalues(row-1, position) = Eigenvalues0[position];
   }
   MPO2.adjoint(MPO1);
   MPO2.multiply(MPO1, MPO0);
   MPO0.getMPS(MPS0);
// compute MPS1 approximating the lower half of the PEPS sandwich:
   this->getMPOSL("up", MPO2);
   for (row = this->N[0]-2; row >= (this->N[0]+1)/2; row--)
   {
    MPO0 = MPO2;
    this->getMPOSL("up", row, MPO1);
    D1 = min(MPO1.getD()*MPO0.getD(), D2s[row-1]);
    MPO2.setD(D1);
    for (position = 0; position < MPO2.getN(); position++)
    {
     MPO2.get(position, Tensor0);
     Tensor0.getShape(Shape4);
     d0 = this->d(row, this->N[1]-1-position);
     Shape4[3] = this->D*d0;
     Tensor0 = Tensor<T>(Shape4);
     Tensor0.fillRandomly();
     MPO2.set(position, Tensor0);
    }
    eps0 = Epss[row-1];
    maxNumSweeps = MaxNumsSweeps[row-1];
    MPO1.multiply(MPO0, eps0, maxNumSweeps, errorAchieved, numSweepsDone, MPO2);
    ErrorsAchieved[row-1] = errorAchieved;
    NumsSweepsDone[row-1] = numSweepsDone;
    for (position = 0; position < MPO2.getN(); position++)
    {
     MPO2.get(position, Tensor0);
     Tensor0.getShape(Shape4);
     d0 = this->d(row, this->N[1]-1-position);
     Shape5[0] = Shape4[0]; Shape5[1] = Shape4[1]; Shape5[2] = Shape4[2]; Shape5[3] = Shape4[3]/d0;
     Shape5[4] = d0;
     Tensor0.reshape(Shape5);
     Tensor0.permute(Order);
     Shape4[2] = Shape5[3]; Shape4[3] = Shape5[2]*d0;
     Tensor0.reshape(Shape4);
     MPO2.setd(position, Shape5[2]*d0);
     MPO2.set(position, Tensor0);
    }
    MPO2.setd2(d2s[row-1], eps, Eigenvalues0);
    MPO2.transpose();
    for (position = 0; position < MPO2.getN(); position++)
     MPO2.setd(position, this->D);
    for (position = 0; position < MPO2.getN(); position++)
     Eigenvalues(row-1, this->N[1]-1-position) = Eigenvalues0[position];
   }
   MPO2.adjoint(MPO1);
   MPO2.multiply(MPO1, MPO0);
   MPO0.getMPS(MPS1);
  }
  return MPS0.contractReverse(MPS1);
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented for periodic BC: " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "const vector<unsigned int>& d2s, double eps, " <<
                        "Matrix< vector<double> >& Eigenvalues) const." << endl;
  exit(1);
 }
}

template<class T> T PEPS<T>::scalarProduct(const string& Direction, unsigned int D2, unsigned int d2,
                                           double eps, unsigned int maxNumSweeps, double cutoff,
                                           vector<double>& ErrorsAchieved1, vector<unsigned int>& NumsSweepsDone1,
                                           vector<double>& ErrorsAchieved2, vector<unsigned int>& NumsSweepsDone2)
                                          const
{
#ifdef DEBUG
 if ((this->N[0] == 0) || ((Direction != "horizontal") && (Direction != "vertical")) || (D2 == 0) || (d2 == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, unsigned int D2, unsigned int d2, " <<
                        "double eps, unsigned int maxNumSweeps, double cutoff, " <<
                        "vector<double>& ErrorsAchieved1, vector<unsigned int>& NumsSweepsDone1, " <<
                        "vector<double>& ErrorsAchieved2, vector<unsigned int>& NumsSweepsDone2) const: " <<
          "((this->N[0] == 0) || ((Direction != horizontal) && (Direction != vertical)) || " <<
           "(D2 == 0) || (d2 == 0))." << endl;
  exit(1);
 }
 if (((Direction == "horizontal") && (ErrorsAchieved1.size() != this->N[1]-2)) ||
     ((Direction == "horizontal") && (NumsSweepsDone1.size() != this->N[1]-2)) ||
     ((Direction == "horizontal") && (ErrorsAchieved2.size() != this->N[1]-2)) ||
     ((Direction == "horizontal") && (NumsSweepsDone2.size() != this->N[1]-2)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, unsigned int D2, unsigned int d2, " <<
                        "double eps, unsigned int maxNumSweeps, double cutoff, " <<
                        "vector<double>& ErrorsAchieved1, vector<unsigned int>& NumsSweepsDone1, " <<
                        "vector<double>& ErrorsAchieved2, vector<unsigned int>& NumsSweepsDone2) const: " <<
          "(((Direction == horizontal) && (ErrorsAchieved1.size() != this->N[1]-2)) || " <<
           "((Direction == horizontal) && (NumsSweepsDone1.size() != this->N[1]-2)) || " <<
           "((Direction == horizontal) && (ErrorsAchieved2.size() != this->N[1]-2)) || " <<
           "((Direction == horizontal) && (NumsSweepsDone2.size() != this->N[1]-2)))." << endl;
  exit(1);
 }
 else if (((Direction == "vertical") && (ErrorsAchieved1.size() != this->N[0]-2)) ||
          ((Direction == "vertical") && (NumsSweepsDone1.size() != this->N[0]-2)) ||
          ((Direction == "vertical") && (ErrorsAchieved2.size() != this->N[0]-2)) ||
          ((Direction == "vertical") && (NumsSweepsDone2.size() != this->N[0]-2)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, unsigned int D2, unsigned int d2, " <<
                        "double eps, unsigned int maxNumSweeps, double cutoff, " <<
                        "vector<double>& ErrorsAchieved1, vector<unsigned int>& NumsSweepsDone1, " <<
                        "vector<double>& ErrorsAchieved2, vector<unsigned int>& NumsSweepsDone2) const: " <<
          "(((Direction == vertical) && (ErrorsAchieved1.size() != this->N[0]-2)) || " <<
           "((Direction == vertical) && (NumsSweepsDone1.size() != this->N[0]-2)) || " <<
           "((Direction == vertical) && (ErrorsAchieved2.size() != this->N[0]-2)) || " <<
           "((Direction == vertical) && (NumsSweepsDone2.size() != this->N[0]-2)))." << endl;
  exit(1);
 }
#endif
 unsigned int row, col, position;
 unsigned int numSweepsDone, D1, d0;
 double errorAchieved;
 vector<unsigned int> Shape4(4), Shape5(5), Order(5);
 Tensor<T> Tensor0;
 MPS<T> MPS0, MPS1; MPO<T> MPO0, MPO1, MPO2;
 Order[0] = 0; Order[1] = 1; Order[2] = 3; Order[3] = 2; Order[4] = 4;
 if (this->BC == "open")
 {
  if (Direction == "horizontal")
  {
// compute MPS0 approximating the left half of the PEPS sandwich:
   this->getMPOSL("right", MPO2);
   for (col = 1; col < (this->N[1]+1)/2; col++)
   {
    MPO0 = MPO2;
    this->getMPOSL("right", col, MPO1);
    D1 = min(MPO1.getD()*MPO0.getD(), D2);
    MPO2.setD(D1);
    for (position = 0; position < MPO2.getN(); position++)
    {
     MPO2.get(position, Tensor0);
     Tensor0.getShape(Shape4);
     d0 = this->d(this->N[0]-1-position, col);
     Shape4[3] = this->D*d0;
     Tensor0 = Tensor<T>(Shape4);
     Tensor0.fillRandomly(time(0)+7*position);
     MPO2.set(position, Tensor0);
    }
    MPO1.multiply(MPO0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPO2);
    ErrorsAchieved1[col-1] = errorAchieved;
    NumsSweepsDone1[col-1] = numSweepsDone;
    for (position = 0; position < MPO2.getN(); position++)
    {
     MPO2.get(position, Tensor0);
     Tensor0.getShape(Shape4);
     d0 = this->d(this->N[0]-1-position, col);
     Shape5[0] = Shape4[0]; Shape5[1] = Shape4[1]; Shape5[2] = Shape4[2]; Shape5[3] = Shape4[3]/d0;
     Shape5[4] = d0;
     Tensor0.reshape(Shape5);
     Tensor0.permute(Order);
     Shape4[2] = Shape5[3]; Shape4[3] = Shape5[2]*d0;
     Tensor0.reshape(Shape4);
     MPO2.setd(position, Shape5[2]*d0);
     MPO2.set(position, Tensor0);
    }
    MPO2.setd2(d2, cutoff, eps, maxNumSweeps, errorAchieved, numSweepsDone);
    ErrorsAchieved2[col-1] = errorAchieved;
    NumsSweepsDone2[col-1] = numSweepsDone;
    MPO2.transpose();
    for (position = 0; position < MPO2.getN(); position++)
     MPO2.setd(position, this->D);
   }
   MPO2.adjoint(MPO1);
   MPO2.multiply(MPO1, MPO0);
   MPO0.getMPS(MPS0);
// compute MPS1 approximating the right half of the PEPS sandwich:
   this->getMPOSL("left", MPO2);
   for (col = this->N[1]-2; col >= (this->N[1]+1)/2; col--)
   {
    MPO0 = MPO2;
    this->getMPOSL("left", col, MPO1);
    D1 = min(MPO1.getD()*MPO0.getD(), D2);
    MPO2.setD(D1);
    for (position = 0; position < MPO2.getN(); position++)
    {
     MPO2.get(position, Tensor0);
     Tensor0.getShape(Shape4);
     d0 = this->d(position, col);
     Shape4[3] = this->D*d0;
     Tensor0 = Tensor<T>(Shape4);
     Tensor0.fillRandomly(time(0)+7*position);
     MPO2.set(position, Tensor0);
    }
    MPO1.multiply(MPO0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPO2);
    ErrorsAchieved1[col-1] = errorAchieved;
    NumsSweepsDone1[col-1] = numSweepsDone;
    for (position = 0; position < MPO2.getN(); position++)
    {
     MPO2.get(position, Tensor0);
     Tensor0.getShape(Shape4);
     d0 = this->d(position, col);
     Shape5[0] = Shape4[0]; Shape5[1] = Shape4[1]; Shape5[2] = Shape4[2]; Shape5[3] = Shape4[3]/d0;
     Shape5[4] = d0;
     Tensor0.reshape(Shape5);
     Tensor0.permute(Order);
     Shape4[2] = Shape5[3]; Shape4[3] = Shape5[2]*d0;
     Tensor0.reshape(Shape4);
     MPO2.setd(position, Shape5[2]*d0);
     MPO2.set(position, Tensor0);
    }
    MPO2.setd2(d2, cutoff, eps, maxNumSweeps, errorAchieved, numSweepsDone);
    ErrorsAchieved2[col-1] = errorAchieved;
    NumsSweepsDone2[col-1] = numSweepsDone;
    MPO2.transpose();
    for (position = 0; position < MPO2.getN(); position++)
     MPO2.setd(position, this->D);
   }
   MPO2.adjoint(MPO1);
   MPO2.multiply(MPO1, MPO0);
   MPO0.getMPS(MPS1);
  }
  else if (Direction == "vertical")
  {
// compute MPS0 approximating the upper half of the PEPS sandwich:
   this->getMPOSL("down", MPO2);
   for (row = 1; row < (this->N[0]+1)/2; row++)
   {
    MPO0 = MPO2;
    this->getMPOSL("down", row, MPO1);
    D1 = min(MPO1.getD()*MPO0.getD(), D2);
    MPO2.setD(D1);
    for (position = 0; position < MPO2.getN(); position++)
    {
     MPO2.get(position, Tensor0);
     Tensor0.getShape(Shape4);
     d0 = this->d(row, position);
     Shape4[3] = this->D*d0;
     Tensor0 = Tensor<T>(Shape4);
     Tensor0.fillRandomly(time(0)+7*position);
     MPO2.set(position, Tensor0);
    }
    MPO1.multiply(MPO0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPO2);
    ErrorsAchieved1[row-1] = errorAchieved;
    NumsSweepsDone1[row-1] = numSweepsDone;
    for (position = 0; position < MPO2.getN(); position++)
    {
     MPO2.get(position, Tensor0);
     Tensor0.getShape(Shape4);
     d0 = this->d(row, position);
     Shape5[0] = Shape4[0]; Shape5[1] = Shape4[1]; Shape5[2] = Shape4[2]; Shape5[3] = Shape4[3]/d0;
     Shape5[4] = d0;
     Tensor0.reshape(Shape5);
     Tensor0.permute(Order);
     Shape4[2] = Shape5[3]; Shape4[3] = Shape5[2]*d0;
     Tensor0.reshape(Shape4);
     MPO2.setd(position, Shape5[2]*d0);
     MPO2.set(position, Tensor0);
    }
    MPO2.setd2(d2, cutoff, eps, maxNumSweeps, errorAchieved, numSweepsDone);
    ErrorsAchieved2[row-1] = errorAchieved;
    NumsSweepsDone2[row-1] = numSweepsDone;
    MPO2.transpose();
    for (position = 0; position < MPO2.getN(); position++)
     MPO2.setd(position, this->D);
   }
   MPO2.adjoint(MPO1);
   MPO2.multiply(MPO1, MPO0);
   MPO0.getMPS(MPS0);
// compute MPS1 approximating the lower half of the PEPS sandwich:
   this->getMPOSL("up", MPO2);
   for (row = this->N[0]-2; row >= (this->N[0]+1)/2; row--)
   {
    MPO0 = MPO2;
    this->getMPOSL("up", row, MPO1);
    D1 = min(MPO1.getD()*MPO0.getD(), D2);
    MPO2.setD(D1);
    for (position = 0; position < MPO2.getN(); position++)
    {
     MPO2.get(position, Tensor0);
     Tensor0.getShape(Shape4);
     d0 = this->d(row, this->N[1]-1-position);
     Shape4[3] = this->D*d0;
     Tensor0 = Tensor<T>(Shape4);
     Tensor0.fillRandomly(time(0)+7*position);
     MPO2.set(position, Tensor0);
    }
    MPO1.multiply(MPO0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPO2);
    ErrorsAchieved1[row-1] = errorAchieved;
    NumsSweepsDone1[row-1] = numSweepsDone;
    for (position = 0; position < MPO2.getN(); position++)
    {
     MPO2.get(position, Tensor0);
     Tensor0.getShape(Shape4);
     d0 = this->d(row, this->N[1]-1-position);
     Shape5[0] = Shape4[0]; Shape5[1] = Shape4[1]; Shape5[2] = Shape4[2]; Shape5[3] = Shape4[3]/d0;
     Shape5[4] = d0;
     Tensor0.reshape(Shape5);
     Tensor0.permute(Order);
     Shape4[2] = Shape5[3]; Shape4[3] = Shape5[2]*d0;
     Tensor0.reshape(Shape4);
     MPO2.setd(position, Shape5[2]*d0);
     MPO2.set(position, Tensor0);
    }
    MPO2.setd2(d2, cutoff, eps, maxNumSweeps, errorAchieved, numSweepsDone);
    ErrorsAchieved2[row-1] = errorAchieved;
    NumsSweepsDone2[row-1] = numSweepsDone;
    MPO2.transpose();
    for (position = 0; position < MPO2.getN(); position++)
     MPO2.setd(position, this->D);
   }
   MPO2.adjoint(MPO1);
   MPO2.multiply(MPO1, MPO0);
   MPO0.getMPS(MPS1);
  }
  return MPS0.contractReverse(MPS1);
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented for periodic BC: " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, unsigned int D2, unsigned int d2, " <<
                        "double eps, unsigned int maxNumSweeps, double cutoff, " <<
                        "vector<double>& ErrorsAchieved1, vector<unsigned int>& NumsSweepsDone1, " <<
                        "vector<double>& ErrorsAchieved2, vector<unsigned int>& NumsSweepsDone2) const." << endl;
  exit(1);
 }
}

template<class T> T PEPS<T>::scalarProduct(const string& Direction, const vector<unsigned int>& D2s,
                                           const vector<unsigned int>& d2s,
                                           const vector<double>& Epss,
                                           const vector<unsigned int>& MaxNumsSweeps,
                                           vector<double>& ErrorsAchieved,
                                           vector<unsigned int>& NumsSweepsDone) const
{
#ifdef DEBUG
 if ((this->N[0] == 0) || ((Direction != "horizontal") && (Direction != "vertical")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<unsigned int>& d2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
          "((this->N[0] == 0) || ((Direction != horizontal) && (Direction != vertical)))." << endl;
  exit(1);
 }
 if (((Direction == "horizontal") && (D2s.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (d2s.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (Epss.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (MaxNumsSweeps.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (ErrorsAchieved.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (NumsSweepsDone.size() != this->N[1]-3)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<unsigned int>& d2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
          "(((Direction == horizontal) && (D2s.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (d2s.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (Epss.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (MaxNumsSweeps.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (ErrorsAchieved.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (NumsSweepsDone.size() != this->N[1]-3)))." << endl;
  exit(1);
 }
 if (((Direction == "vertical") && (D2s.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (d2s.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (Epss.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (MaxNumsSweeps.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (ErrorsAchieved.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (NumsSweepsDone.size() != this->N[0]-3)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<unsigned int>& d2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
          "(((Direction == vertical) && (D2s.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (d2s.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (Epss.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (MaxNumsSweeps.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (ErrorsAchieved.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (NumsSweepsDone.size() != this->N[0]-3)))." << endl;
  exit(1);
 }
 for (int i = 0; i < D2s.size(); i++)
 {
  if (D2s[i] == 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> T PEPS<T>::" <<
           "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                         "const vector<unsigned int>& d2s, " <<
                         "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                         "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
           "(D2s[" << i << "] == 0)." << endl;
   exit(1);
  }
 }
#endif
 unsigned int mode = 2; double cutoff = 1.0e-12;
 unsigned int d2, maxNumSweeps, numSweepsDone;
 double eps, errorAchieved;
 MPS<T> MPS0, MPS1, MPS2; MPO<T> MPO0;
 unsigned int D1;
 if (this->BC == "open")
 {
  if (Direction == "horizontal")
  {
// compute MPS1 approximating the left half of the PEPS sandwich:
   this->getMPS("right", MPS1);
   for (int col = 1; col < this->N[1]/2; col++)
   {
    MPS0 = MPS1;
    this->getMPO("right", col, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[col-1]);
    MPS1.setD(D1);
    d2 = d2s[col-1];
    eps = Epss[col-1];
    maxNumSweeps = MaxNumsSweeps[col-1];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1, cutoff, mode, d2);
    ErrorsAchieved[col-1] = errorAchieved;
    NumsSweepsDone[col-1] = numSweepsDone;
   }
// compute MPS2 approximating the right half of the PEPS sandwich:
   this->getMPS("left", MPS2);
   for (int col = this->N[1]-2; col > this->N[1]/2; col--)
   {
    MPS0 = MPS2;
    this->getMPO("left", col, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[col-2]);
    MPS2.setD(D1);
    d2 = d2s[col-2];
    eps = Epss[col-2];
    maxNumSweeps = MaxNumsSweeps[col-2];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS2, cutoff, mode, d2);
    ErrorsAchieved[col-2] = errorAchieved;
    NumsSweepsDone[col-2] = numSweepsDone;
   }
// return the contraction value:
   if (this->N[1] == 2)
   {
    return MPS2.contractReverse(MPS1);
   }
   else
   {
    this->getMPO("right", this->N[1]/2, MPO0);
    return contractReverseMPOMPS(MPS2, MPO0, MPS1);
   }
  }
  else if (Direction == "vertical")
  {
// compute MPS1 approximating the upper half of the PEPS sandwich:
   this->getMPS("down", MPS1);
   for (int row = 1; row < this->N[0]/2; row++)
   {
    MPS0 = MPS1;
    this->getMPO("down", row, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[row-1]);
    MPS1.setD(D1);
    d2 = d2s[row-1];
    eps = Epss[row-1];
    maxNumSweeps = MaxNumsSweeps[row-1];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1, cutoff, mode, d2);
    ErrorsAchieved[row-1] = errorAchieved;
    NumsSweepsDone[row-1] = numSweepsDone;
   }
// compute MPS2 approximating the lower half of the PEPS sandwich:
   this->getMPS("up", MPS2);
   for (int row = this->N[0]-2; row > this->N[0]/2; row--)
   {
    MPS0 = MPS2;
    this->getMPO("up", row, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[row-2]);
    MPS2.setD(D1);
    d2 = d2s[row-2];
    eps = Epss[row-2];
    maxNumSweeps = MaxNumsSweeps[row-2];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS2, cutoff, mode, d2);
    ErrorsAchieved[row-2] = errorAchieved;
    NumsSweepsDone[row-2] = numSweepsDone;
   }
// return the contraction value:
   if (this->N[0] == 2)
   {
    return MPS2.contractReverse(MPS1);
   }
   else
   {
    this->getMPO("down", this->N[0]/2, MPO0);
    return contractReverseMPOMPS(MPS2, MPO0, MPS1);
   }
  }
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented for periodic BC: " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<unsigned int>& d2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const." << endl;
  exit(1);
 }
}

template<class T> T PEPS<T>::normalize(const string& Direction,
                                       unsigned int D2, double eps, unsigned int maxNumSweeps)
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "normalize(const string& Direction, " <<
                    "unsigned int D2, double eps, unsigned int maxNumSweeps): " <<
          "(this->N[0] == 0)." << endl;
 }
 if ((Direction != "horizontal") && (Direction != "vertical"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "normalize(const string& Direction, " <<
                    "unsigned int D2, double eps, unsigned int maxNumSweeps): " <<
          "((Direction != horizontal) && (Direction != vertical))." << endl;
 }
 if ((D2 == 0) || (maxNumSweeps == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "normalize(const string& Direction, " <<
                    "unsigned int D2, double eps, unsigned int maxNumSweeps): " <<
          "((D2 == 0) || (maxNumSweeps == 0))." << endl;
 }
#endif
 unsigned int Ncols = this->N[1], Nrows = this->N[0];
 T result;
 vector<unsigned int> NumsSweepsDone;
 vector<double> EpssAchieved;
 vector<T> Norms;
// 1. divide each tensor by its largest element absolute:
 this->normalize();
// 2. compute norm and correct factor for each tensor to normalize this PEPS:
// 2.a) horizontal contraction:
 if (Direction == "horizontal")
 {
  NumsSweepsDone = vector<unsigned int>(Ncols-1);
  EpssAchieved = vector<double>(Ncols-1);
  Norms = vector<T>(Ncols-1);
  result = this->scalarProduct(Direction, D2, eps, maxNumSweeps, EpssAchieved, NumsSweepsDone, Norms);
  result = pow(abs(result), 1.0/(2.0*Nrows));
  for (int i = 0; i < Ncols-1; i++)
   result *= pow(abs(Norms[i]), 1.0/(2.0*Nrows));
  result = pow(result, 1.0/Ncols);
 }
// 2.b) vertical contraction:
 else if (Direction == "vertical")
 {
  NumsSweepsDone = vector<unsigned int>(Nrows-1);
  EpssAchieved = vector<double>(Nrows-1);
  Norms = vector<T>(Nrows-1);
  result = this->scalarProduct(Direction, D2, eps, maxNumSweeps, EpssAchieved, NumsSweepsDone, Norms);
  result = pow(abs(result), 1.0/(2.0*Ncols));
  for (int i = 0; i < Nrows-1; i++)
   result *= pow(abs(Norms[i]), 1.0/(2.0*Ncols));
  result = pow(result, 1.0/Nrows);
 }
 result = 1.0/result;
// 3. multiply each tensor with correct factor to normalize this PEPS:
 for (int col = 0; col < Ncols; col++){
  for (int row = 0; row < Nrows; row++){
   this->Tensors[row+Nrows*col].multiply(result);
  }
 }
 return result;
}

template<class T> T PEPS<T>::scalarProduct(const PEPS<T>& PEPS0, const string& Direction,
                                           const vector<unsigned int>& D2s,
                                           const vector<double>& Epss,
                                           const vector<unsigned int>& MaxNumsSweeps,
                                           vector<double>& ErrorsAchieved,
                                           vector<unsigned int>& NumsSweepsDone) const
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (this->BC != PEPS0.BC) || (this->N != PEPS0.N) || (this->d != PEPS0.d) ||
     ((Direction != "horizontal") && (Direction != "vertical")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const PEPS<T>& PEPS0, const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
          "((this->N[0] == 0) || (this->BC != PEPS0.BC) || (this->N != PEPS0.N) || (this->d != PEPS0.d) || " <<
           "((Direction != horizontal) && (Direction != vertical)))." << endl;
  exit(1);
 }
 if (((Direction == "horizontal") && (D2s.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (Epss.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (MaxNumsSweeps.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (ErrorsAchieved.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (NumsSweepsDone.size() != this->N[1]-3)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const PEPS<T>& PEPS0, const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
          "(((Direction == horizontal) && (D2s.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (Epss.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (MaxNumsSweeps.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (ErrorsAchieved.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (NumsSweepsDone.size() != this->N[1]-3)))." << endl;
  exit(1);
 }
 if (((Direction == "vertical") && (D2s.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (Epss.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (MaxNumsSweeps.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (ErrorsAchieved.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (NumsSweepsDone.size() != this->N[0]-3)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const PEPS<T>& PEPS0, const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
          "(((Direction == vertical) && (D2s.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (Epss.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (MaxNumsSweeps.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (ErrorsAchieved.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (NumsSweepsDone.size() != this->N[0]-3)))." << endl;
  exit(1);
 }
 for (int i = 0; i < D2s.size(); i++)
 {
  if (D2s[i] == 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> T PEPS<T>::" <<
           "scalarProduct(const PEPS<T>& PEPS0, const string& Direction, const vector<unsigned int>& D2s, " <<
                         "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                         "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
           "(D2s[" << i << "] == 0)." << endl;
   exit(1);
  }
 }
#endif
 double eps, errorAchieved;
 unsigned int maxNumSweeps, numSweepsDone;
 MPS<T> MPS0, MPS1, MPS2; MPO<T> MPO0;
 unsigned int D1;
 if (this->BC == "open")
 {
  if (Direction == "horizontal")
  {
// compute MPS1 approximating the left half of the PEPS sandwich:
   this->getMPS(PEPS0, "right", MPS1);
   for (int col = 1; col < this->N[1]/2; col++)
   {
    MPS0 = MPS1;
    this->getMPO(PEPS0, "right", col, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[col-1]);
    MPS1.setD(D1);
    eps = Epss[col-1];
    maxNumSweeps = MaxNumsSweeps[col-1];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
    ErrorsAchieved[col-1] = errorAchieved;
    NumsSweepsDone[col-1] = numSweepsDone;
   }
// compute MPS2 approximating the right half of the PEPS sandwich:
   this->getMPS(PEPS0, "left", MPS2);
   for (int col = this->N[1]-2; col > this->N[1]/2; col--)
   {
    MPS0 = MPS2;
    this->getMPO(PEPS0, "left", col, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[col-2]);
    MPS2.setD(D1);
    eps = Epss[col-2];
    maxNumSweeps = MaxNumsSweeps[col-2];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS2);
    ErrorsAchieved[col-2] = errorAchieved;
    NumsSweepsDone[col-2] = numSweepsDone;
   }
// return the contraction value:
   if (this->N[1] == 2)
   {
    return MPS2.contractReverse(MPS1);
   }
   else
   {
    this->getMPO(PEPS0, "right", this->N[1]/2, MPO0);
    return contractReverseMPOMPS(MPS2, MPO0, MPS1);
   }
  }
  else if (Direction == "vertical")
  {
// compute MPS1 approximating the upper half of the PEPS sandwich:
   this->getMPS(PEPS0, "down", MPS1);
   for (int row = 1; row < this->N[0]/2; row++)
   {
    MPS0 = MPS1;
    this->getMPO(PEPS0, "down", row, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[row-1]);
    MPS1.setD(D1);
    eps = Epss[row-1];
    maxNumSweeps = MaxNumsSweeps[row-1];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
    ErrorsAchieved[row-1] = errorAchieved;
    NumsSweepsDone[row-1] = numSweepsDone;
   }
// compute MPS2 approximating the lower half of the PEPS sandwich:
   this->getMPS(PEPS0, "up", MPS2);
   for (int row = this->N[0]-2; row > this->N[0]/2; row--)
   {
    MPS0 = MPS2;
    this->getMPO(PEPS0, "up", row, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[row-2]);
    MPS2.setD(D1);
    eps = Epss[row-2];
    maxNumSweeps = MaxNumsSweeps[row-2];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS2);
    ErrorsAchieved[row-2] = errorAchieved;
    NumsSweepsDone[row-2] = numSweepsDone;
   }
// return the contraction value:
   if (this->N[0] == 2)
   {
    return MPS2.contractReverse(MPS1);
   }
   else
   {
    this->getMPO(PEPS0, "down", this->N[0]/2, MPO0);
    return contractReverseMPOMPS(MPS2, MPO0, MPS1);
   }
  }
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented for periodic BC: " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const PEPS<T>& PEPS0, const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, " <<
                        "vector<unsigned int>& NumsSweepsDone) const." << endl;
  exit(1);
 }
}

template<class T> T PEPS<T>::scalarProduct(const PEPS<T>& PEPS0, const string& Direction,
                                           const vector<unsigned int>& D2s,
                                           const vector<double>& Epss,
                                           const vector<unsigned int>& MaxNumsSweeps,
                                           vector<double>& ErrorsAchieved,
                                           vector<unsigned int>& NumsSweepsDone,
                                           vector<T>& Norms) const
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (this->BC != PEPS0.BC) || (this->N != PEPS0.N) || (this->d != PEPS0.d) ||
     ((Direction != "horizontal") && (Direction != "vertical")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const PEPS<T>& PEPS0, const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "((this->N[0] == 0) || (this->BC != PEPS0.BC) || (this->N != PEPS0.N) || (this->d != PEPS0.d) || " <<
           "((Direction != horizontal) && (Direction != vertical)))." << endl;
  exit(1);
 }
 if (((Direction == "horizontal") && (D2s.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (Epss.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (MaxNumsSweeps.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (ErrorsAchieved.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (NumsSweepsDone.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (Norms.size() != this->N[1])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const PEPS<T>& PEPS0, const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "(((Direction == horizontal) && (D2s.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (Epss.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (MaxNumsSweeps.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (ErrorsAchieved.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (NumsSweepsDone.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (Norms.size() != this->N[1])))." << endl;
  exit(1);
 }
 if (((Direction == "vertical") && (D2s.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (Epss.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (MaxNumsSweeps.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (ErrorsAchieved.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (NumsSweepsDone.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (Norms.size() != this->N[0])))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const PEPS<T>& PEPS0, const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "(((Direction == vertical) && (D2s.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (Epss.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (MaxNumsSweeps.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (ErrorsAchieved.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (NumsSweepsDone.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (Norms.size() != this->N[0])))." << endl;
  exit(1);
 }
 for (int i = 0; i < D2s.size(); i++)
 {
  if (D2s[i] == 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> T PEPS<T>::" <<
           "scalarProduct(const PEPS<T>& PEPS0, const string& Direction, const vector<unsigned int>& D2s, " <<
                         "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                         "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                         "vector<T>& Norms) const: " <<
           "(D2s[" << i << "] == 0)." << endl;
   exit(1);
  }
 }
#endif
 double eps, errorAchieved;
 unsigned int maxNumSweeps, numSweepsDone;
 MPS<T> MPS0, MPS1, MPS2; MPO<T> MPO0;
 unsigned int D1;
 if (this->BC == "open")
 {
  if (Direction == "horizontal")
  {
// compute MPS1 approximating the left half of the PEPS sandwich:
   this->getMPS(PEPS0, "right", MPS1);
   Norms[0] = MPS1.simplifiedNormalize();
   for (int col = 1; col < this->N[1]/2; col++)
   {
    MPS0 = MPS1;
    this->getMPO(PEPS0, "right", col, MPO0);
    Norms[col] = MPO0.simplifiedNormalize();
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[col-1]);
    MPS1.setD(D1);
    eps = Epss[col-1];
    maxNumSweeps = MaxNumsSweeps[col-1];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
    ErrorsAchieved[col-1] = errorAchieved;
    NumsSweepsDone[col-1] = numSweepsDone;
   }
// compute MPS2 approximating the right half of the PEPS sandwich:
   this->getMPS(PEPS0, "left", MPS2);
   Norms[this->N[1]-1] = MPS2.simplifiedNormalize();
   for (int col = this->N[1]-2; col > this->N[1]/2; col--)
   {
    MPS0 = MPS2;
    this->getMPO(PEPS0, "left", col, MPO0);
    Norms[col] = MPO0.simplifiedNormalize();
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[col-2]);
    MPS2.setD(D1);
    eps = Epss[col-2];
    maxNumSweeps = MaxNumsSweeps[col-2];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS2);
    ErrorsAchieved[col-2] = errorAchieved;
    NumsSweepsDone[col-2] = numSweepsDone;
   }
// return the contraction value:
   if (this->N[1] == 2)
   {
    return MPS2.contractReverse(MPS1);
   }
   else
   {
    this->getMPO(PEPS0, "right", this->N[1]/2, MPO0);
    Norms[this->N[1]/2] = MPO0.simplifiedNormalize();
    return contractReverseMPOMPS(MPS2, MPO0, MPS1);
   }
  }
  else if (Direction == "vertical")
  {
// compute MPS1 approximating the upper half of the PEPS sandwich:
   this->getMPS(PEPS0, "down", MPS1);
   Norms[0] = MPS1.simplifiedNormalize();
   for (int row = 1; row < this->N[0]/2; row++)
   {
    MPS0 = MPS1;
    this->getMPO(PEPS0, "down", row, MPO0);
    Norms[row] = MPO0.simplifiedNormalize();
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[row-1]);
    MPS1.setD(D1);
    eps = Epss[row-1];
    maxNumSweeps = MaxNumsSweeps[row-1];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
    ErrorsAchieved[row-1] = errorAchieved;
    NumsSweepsDone[row-1] = numSweepsDone;
   }
// compute MPS2 approximating the lower half of the PEPS sandwich:
   this->getMPS(PEPS0, "up", MPS2);
   Norms[this->N[0]-1] = MPS2.simplifiedNormalize();
   for (int row = this->N[0]-2; row > this->N[0]/2; row--)
   {
    MPS0 = MPS2;
    this->getMPO(PEPS0, "up", row, MPO0);
    Norms[row] = MPO0.simplifiedNormalize();
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[row-2]);
    MPS2.setD(D1);
    eps = Epss[row-2];
    maxNumSweeps = MaxNumsSweeps[row-2];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS2);
    ErrorsAchieved[row-2] = errorAchieved;
    NumsSweepsDone[row-2] = numSweepsDone;
   }
// return the contraction value:
   if (this->N[0] == 2)
   {
    return MPS2.contractReverse(MPS1);
   }
   else
   {
    this->getMPO(PEPS0, "down", this->N[0]/2, MPO0);
    Norms[this->N[0]/2] = MPO0.simplifiedNormalize();
    return contractReverseMPOMPS(MPS2, MPO0, MPS1);
   }
  }
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented for periodic BC: " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const PEPS<T>& PEPS0, const string& Direction, const vector<unsigned int>& D2s, " <<
                        "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                        "vector<double>& ErrorsAchieved, " <<
                        "vector<unsigned int>& NumsSweepsDone, vector<T>& Norms) const." << endl;
  exit(1);
 }
}

template<class T> T PEPS<T>::scalarProduct(const PEPS<T>& PEPS0, const string& Direction,
                                           unsigned int D2, double eps, unsigned int maxNumSweeps,
                                           vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone,
                                           vector<T>& Norms) const
{
#ifdef DEBUG
 if ((this->BC != "open") || (this->N[0] == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const PEPS<T>& PEPS0, const string& Direction, " <<
                        "unsigned int D2, double eps, unsigned int maxNumSweeps, " <<
                        "vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "((this->BC != open) || (this->N[0] == 0))." << endl;
 }
 if ((this->BC != PEPS0.BC) || (this->N != PEPS0.N) || (this->d != PEPS0.d))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const PEPS<T>& PEPS0, const string& Direction, " <<
                        "unsigned int D2, double eps, unsigned int maxNumSweeps, " <<
                        "vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "((this->BC != PEPS0.BC) || (this->N != PEPS0.N) || (this->d != PEPS0.d))." << endl;
 }
 if ((Direction != "horizontal") && (Direction != "vertical"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const PEPS<T>& PEPS0, const string& Direction, " <<
                        "unsigned int D2, double eps, unsigned int maxNumSweeps, " <<
                        "vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "((Direction != horizontal) && (Direction != vertical))." << endl;
 }
 if ((D2 == 0) || (maxNumSweeps == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const PEPS<T>& PEPS0, const string& Direction, " <<
                        "unsigned int D2, double eps, unsigned int maxNumSweeps, " <<
                        "vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "((D2 == 0) || (maxNumSweeps == 0))." << endl;
 }
 if (((Direction == "horizontal") && (EpssAchieved.size() != this->N[1]-1)) ||
     ((Direction == "horizontal") && (NumsSweepsDone.size() != this->N[1]-1)) ||
     ((Direction == "horizontal") && (Norms.size() != this->N[1]-1)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const PEPS<T>& PEPS0, const string& Direction, " <<
                        "unsigned int D2, double eps, unsigned int maxNumSweeps, " <<
                        "vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "(((Direction == horizontal) && (EpssAchieved.size() != this->N[1]-1)) || " <<
           "((Direction == horizontal) && (NumsSweepsDone.size() != this->N[1]-1)) || " <<
           "((Direction == horizontal) && (Norms.size() != this->N[1]-1)))." << endl;
  exit(1);
 }
 if (((Direction == "vertical") && (EpssAchieved.size() != this->N[0]-1)) ||
     ((Direction == "vertical") && (NumsSweepsDone.size() != this->N[0]-1)) ||
     ((Direction == "vertical") && (Norms.size() != this->N[0]-1)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "scalarProduct(const PEPS<T>& PEPS0, const string& Direction, " <<
                        "unsigned int D2, double eps, unsigned int maxNumSweeps, " <<
                        "vector<double>& EpssAchieved, vector<unsigned int>& NumsSweepsDone, " <<
                        "vector<T>& Norms) const: " <<
          "(((Direction == vertical) && (EpssAchieved.size() != this->N[0]-1)) || " <<
           "((Direction == vertical) && (NumsSweepsDone.size() != this->N[0]-1)) || " <<
           "((Direction == vertical) && (Norms.size() != this->N[0]-1)))." << endl;
  exit(1);
 }
#endif
 unsigned int L, Ncols = this->N[1], Nrows = this->N[0];
 T result;
 string BCMPS = "open", Boundary, SweepDirection;
 vector<unsigned int> dMPS, dMPS1, Seed, Shape(5), Shape0(5);
 MPS<T> MPS0, MPS1, MPS2;
// a) horizontal contraction:
 if (Direction == "horizontal")
 {
  L = Nrows;
  dMPS = vector<unsigned int>(L);
  dMPS1 = vector<unsigned int>(L);
  Seed = vector<unsigned int>(L);
// a)1. compute MPS1 approximating the left half of the sandwich:
// a)1.1. boundary-MPS for left boundary:
  for (int i = 0; i < L; i++)
  {
   PEPS0.Tensors[L-1-i].getShape(Shape0);
   this->Tensors[L-1-i].getShape(Shape);
   dMPS[i] = Shape0[2]*Shape[2];
  }
  MPS1 = MPS<T>(BCMPS, L, dMPS, min(D2, this->D*this->D));
  for (int i = 0; i < L; i++)
   Seed[i] = time(0)+13*i;
  MPS1.fillRandomly(Seed);
  Boundary = "left";
  this->getBoundaryMPS(PEPS0, Boundary, eps, maxNumSweeps, EpssAchieved[0], NumsSweepsDone[0], MPS1);
  Norms[0] = MPS1.normalize();
// a)1.2. boundary-MPSs for bulk:
  SweepDirection = "right";
  for (int pos = 1; pos < Ncols/2; pos++)
  {
   MPS0 = MPS1;
   for (int i = 0; i < L; i++)
   {
    PEPS0.Tensors[L-1-i+Nrows*pos].getShape(Shape0);
    this->Tensors[L-1-i+Nrows*pos].getShape(Shape);
    dMPS1[i] = Shape0[2]*Shape[2];
   }
   if (dMPS1 == dMPS)
   {
    MPS1.setD(min(D2, this->D*this->D*MPS0.getD()));
   }
   else
   {
    dMPS = dMPS1;
    MPS1 = MPS<T>(BCMPS, L, dMPS, min(D2, this->D*this->D*MPS0.getD()));
    for (int i = 0; i < L; i++)
     Seed[i] = time(0)+13*(i+Nrows*pos);
    MPS1.fillRandomly(Seed);
   }
   this->multiplyBulkMPOMPS(PEPS0, SweepDirection, pos, MPS0, eps, maxNumSweeps, EpssAchieved[pos], NumsSweepsDone[pos], MPS1);
   Norms[pos] = MPS1.normalize();
  }
// a)2. compute MPS2 approximating the right half of the sandwich:
// a)2.1. boundary-MPS for right boundary:
  for (int i = 0; i < L; i++)
  {
   PEPS0.Tensors[i+Nrows*(Ncols-1)].getShape(Shape0);
   this->Tensors[i+Nrows*(Ncols-1)].getShape(Shape);
   dMPS[i] = Shape0[0]*Shape[0];
  }
  MPS2 = MPS<T>(BCMPS, L, dMPS, min(D2, this->D*this->D));
  for (int i = 0; i < L; i++)
   Seed[i] = time(0)+13*(i+Nrows*Ncols/2);
  MPS2.fillRandomly(Seed);
  Boundary = "right";
  this->getBoundaryMPS(PEPS0, Boundary, eps, maxNumSweeps, EpssAchieved[Ncols-2], NumsSweepsDone[Ncols-2], MPS2);
  Norms[Ncols-2] = MPS2.normalize();
// a)2.2. boundary-MPSs for bulk:
  SweepDirection = "left";
  for (int pos = Ncols-2; pos > Ncols/2; pos--)
  {
   MPS0 = MPS2;
   for (int i = 0; i < L; i++)
   {
    PEPS0.Tensors[i+Nrows*pos].getShape(Shape0);
    this->Tensors[i+Nrows*pos].getShape(Shape);
    dMPS1[i] = Shape0[0]*Shape[0];
   }
   if (dMPS1 == dMPS)
   {
    MPS2.setD(min(D2, this->D*this->D*MPS0.getD()));
   }
   else
   {
    dMPS = dMPS1;
    MPS2 = MPS<T>(BCMPS, L, dMPS, min(D2, this->D*this->D*MPS0.getD()));
    for (int i = 0; i < L; i++)
     Seed[i] = time(0)+13*(i+Nrows*(Ncols/2+Ncols-1-pos));
    MPS2.fillRandomly(Seed);
   }
   this->multiplyBulkMPOMPS(PEPS0, SweepDirection, pos, MPS0, eps, maxNumSweeps, EpssAchieved[pos-1], NumsSweepsDone[pos-1], MPS2);
   Norms[pos-1] = MPS2.normalize();
  }
// a)3. perform final contraction:
  result = contractReverseBulkMPOMPS(PEPS0, Direction, MPS1, Ncols/2, MPS2);
 }
// b) vertical contraction:
 else if (Direction == "vertical")
 {
  L = Ncols;
  dMPS = vector<unsigned int>(L);
  dMPS1 = vector<unsigned int>(L);
  Seed = vector<unsigned int>(L);
// b)1. compute MPS1 approximating the upper half of the sandwich:
// b)1.1. boundary-MPS for top boundary:
  for (int i = 0; i < L; i++)
  {
   PEPS0.Tensors[Nrows*i].getShape(Shape0);
   this->Tensors[Nrows*i].getShape(Shape);
   dMPS[i] = Shape0[3]*Shape[3];
  }
  MPS1 = MPS<T>(BCMPS, L, dMPS, min(D2, this->D*this->D));
  for (int i = 0; i < L; i++)
   Seed[i] = time(0)+13*i;
  MPS1.fillRandomly(Seed);
  Boundary = "top";
  this->getBoundaryMPS(PEPS0, Boundary, eps, maxNumSweeps, EpssAchieved[0], NumsSweepsDone[0], MPS1);
  Norms[0] = MPS1.normalize();
// b)1.2. boundary-MPSs for bulk:
  SweepDirection = "down";
  for (int pos = 1; pos < Nrows/2; pos++)
  {
   MPS0 = MPS1;
   for (int i = 0; i < L; i++)
   {
    PEPS0.Tensors[pos+Nrows*i].getShape(Shape0);
    this->Tensors[pos+Nrows*i].getShape(Shape);
    dMPS1[i] = Shape0[3]*Shape[3];
   }
   if (dMPS1 == dMPS)
   {
    MPS1.setD(min(D2, this->D*this->D*MPS0.getD()));
   }
   else
   {
    dMPS = dMPS1;
    MPS1 = MPS<T>(BCMPS, L, dMPS, min(D2, this->D*this->D*MPS0.getD()));
    for (int i = 0; i < L; i++)
     Seed[i] = time(0)+13*(i+Ncols*pos);
    MPS1.fillRandomly(Seed);
   }
   this->multiplyBulkMPOMPS(PEPS0, SweepDirection, pos, MPS0, eps, maxNumSweeps, EpssAchieved[pos], NumsSweepsDone[pos], MPS1);
   Norms[pos] = MPS1.normalize();
  }
// b)2. compute MPS2 approximating the lower half of the sandwich:
// b)2.1. boundary-MPS for bottom boundary:
  for (int i = 0; i < L; i++)
  {
   PEPS0.Tensors[Nrows-1+Nrows*(L-1-i)].getShape(Shape0);
   this->Tensors[Nrows-1+Nrows*(L-1-i)].getShape(Shape);
   dMPS[i] = Shape0[1]*Shape[1];
  }
  MPS2 = MPS<T>(BCMPS, L, dMPS, min(D2, this->D*this->D));
  for (int i = 0; i < L; i++)
   Seed[i] = time(0)+13*(i+Ncols*Nrows/2);
  MPS2.fillRandomly(Seed);
  Boundary = "bottom";
  this->getBoundaryMPS(PEPS0, Boundary, eps, maxNumSweeps, EpssAchieved[Nrows-2], NumsSweepsDone[Nrows-2], MPS2);
  Norms[Nrows-2] = MPS2.normalize();
// b)2.2. boundary-MPSs for bulk:
  SweepDirection = "up";
  for (int pos = Nrows-2; pos > Nrows/2; pos--)
  {
   MPS0 = MPS2;
   for (int i = 0; i < L; i++)
   {
    PEPS0.Tensors[pos+Nrows*(L-1-i)].getShape(Shape0);
    this->Tensors[pos+Nrows*(L-1-i)].getShape(Shape);
    dMPS1[i] = Shape0[1]*Shape[1];
   }
   if (dMPS1 == dMPS)
   {
    MPS2.setD(min(D2, this->D*this->D*MPS0.getD()));
   }
   else
   {
    dMPS = dMPS1;
    MPS2 = MPS<T>(BCMPS, L, dMPS, min(D2, this->D*this->D*MPS0.getD()));
    for (int i = 0; i < L; i++)
     Seed[i] = time(0)+13*(i+Ncols*(Nrows/2+Nrows-1-pos));
    MPS2.fillRandomly(Seed);
   }
   this->multiplyBulkMPOMPS(PEPS0, SweepDirection, pos, MPS0, eps, maxNumSweeps, EpssAchieved[pos-1], NumsSweepsDone[pos-1], MPS2);
   Norms[pos-1] = MPS2.normalize();
  }
// b)3. perform final contraction:
  result = contractReverseBulkMPOMPS(PEPS0, Direction, MPS2, Nrows/2, MPS1);
 }
 return result;
}

template<class T> T PEPS<T>::expectationValue(const PEPO<T>& PEPO0, const string& Direction,
                                              const vector<unsigned int>& D2s,
                                              const vector<double>& Epss,
                                              const vector<unsigned int>& MaxNumsSweeps,
                                              vector<double>& ErrorsAchieved,
                                              vector<unsigned int>& NumsSweepsDone) const
{
#ifdef DEBUG
 string BCPEPO0; PEPO0.getBC(BCPEPO0);
 vector<unsigned int> NPEPO0(2); NPEPO0[0] = PEPO0.getNrows(); NPEPO0[1] = PEPO0.getNcols();
 Matrix<unsigned int> dPEPO0; PEPO0.getd(dPEPO0);
 if ((this->N[0] == 0) || (this->BC != BCPEPO0) || (this->N != NPEPO0) || (this->d != dPEPO0) ||
     ((Direction != "horizontal") && (Direction != "vertical")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "expectationValue(const PEPO<T>& PEPO0, const string& Direction, const vector<unsigned int>& D2s, " <<
                           "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                           "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
          "((this->N[0] == 0) || (this->BC != BCPEPO0) || (this->N != NPEPO0) || (this->d != dPEPO0) || " <<
           "((Direction != horizontal) && (Direction != vertical)))." << endl;
  exit(1);
 }
 if (((Direction == "horizontal") && (D2s.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (Epss.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (MaxNumsSweeps.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (ErrorsAchieved.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (NumsSweepsDone.size() != this->N[1]-3)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "expectationValue(const PEPO<T>& PEPO0, const string& Direction, const vector<unsigned int>& D2s, " <<
                           "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                           "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
          "(((Direction == horizontal) && (D2s.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (Epss.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (MaxNumsSweeps.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (ErrorsAchieved.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (NumsSweepsDone.size() != this->N[1]-3)))." << endl;
  exit(1);
 }
 if (((Direction == "vertical") && (D2s.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (Epss.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (MaxNumsSweeps.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (ErrorsAchieved.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (NumsSweepsDone.size() != this->N[0]-3)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "expectationValue(const PEPO<T>& PEPO0, const string& Direction, const vector<unsigned int>& D2s, " <<
                           "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                           "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
          "(((Direction == vertical) && (D2s.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (Epss.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (MaxNumsSweeps.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (ErrorsAchieved.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (NumsSweepsDone.size() != this->N[0]-3)))." << endl;
  exit(1);
 }
 for (int i = 0; i < D2s.size(); i++)
 {
  if (D2s[i] == 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> T PEPS<T>::" <<
           "expectationValue(const PEPO<T>& PEPO0, const string& Direction, const vector<unsigned int>& D2s, " <<
                            "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                            "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
           "(D2s[" << i << "] == 0)." << endl;
   exit(1);
  }
 }
#endif
// compute |PEPS0>=PEPO0*|this>:
 unsigned D0 = this->D*PEPO0.getD();
 PEPS<T> PEPS0(this->BC, this->N[0], this->N[1], this->d, D0);
 Tensor<T> Tensor0, Tensor1;
 vector<unsigned int> Indices0(1), Indices1(1);
 Indices0[0] = 4; Indices1[0] = 4;
 vector<unsigned int> Order(9), Shape(5);
 Order[0] = 0; Order[1] = 4; Order[2] = 1; Order[3] = 5; Order[4] = 2; Order[5] = 6; Order[6] = 3; Order[7] = 7; Order[8] = 8;
 for (int col = 0; col < this->N[1]; col++)
 {
  for (int row = 0; row < this->N[0]; row++)
  {
   this->get(row, col, Tensor0); PEPO0.get(row, col, Tensor1);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order);
   PEPS0.getOpenBCShape(row, col, Shape);
   Tensor0.reshape(Shape);
   PEPS0.set(row, col, Tensor0);
  }
 }
// perform the contraction:
 double eps, errorAchieved;
 unsigned int maxNumSweeps, numSweepsDone;
 MPS<T> MPS0, MPS1, MPS2; MPO<T> MPO0;
 unsigned int D1;
 if (this->BC == "open")
 {
  if (Direction == "horizontal")
  {
// compute MPS1 approximating the left half of the PEPS sandwich:
   this->getMPS(PEPS0, "right", MPS1);
   for (int col = 1; col < this->N[1]/2; col++)
   {
    MPS0 = MPS1;
    this->getMPO(PEPS0, "right", col, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[col-1]);
    MPS1.setD(D1);
    eps = Epss[col-1];
    maxNumSweeps = MaxNumsSweeps[col-1];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
    ErrorsAchieved[col-1] = errorAchieved;
    NumsSweepsDone[col-1] = numSweepsDone;
   }
// compute MPS2 approximating the right half of the PEPS sandwich:
   this->getMPS(PEPS0, "left", MPS2);
   for (int col = this->N[1]-2; col > this->N[1]/2; col--)
   {
    MPS0 = MPS2;
    this->getMPO(PEPS0, "left", col, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[col-2]);
    MPS2.setD(D1);
    eps = Epss[col-2];
    maxNumSweeps = MaxNumsSweeps[col-2];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS2);
    ErrorsAchieved[col-2] = errorAchieved;
    NumsSweepsDone[col-2] = numSweepsDone;
   }
// return the contraction value:
   if (this->N[1] == 2)
   {
    return MPS2.contractReverse(MPS1);
   }
   else
   {
    this->getMPO(PEPS0, "right", this->N[1]/2, MPO0);
    return contractReverseMPOMPS(MPS2, MPO0, MPS1);
   }
  }
  else if (Direction == "vertical")
  {
// compute MPS1 approximating the upper half of the PEPS sandwich:
   this->getMPS(PEPS0, "down", MPS1);
   for (int row = 1; row < this->N[0]/2; row++)
   {
    MPS0 = MPS1;
    this->getMPO(PEPS0, "down", row, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[row-1]);
    MPS1.setD(D1);
    eps = Epss[row-1];
    maxNumSweeps = MaxNumsSweeps[row-1];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
    ErrorsAchieved[row-1] = errorAchieved;
    NumsSweepsDone[row-1] = numSweepsDone;
   }
// compute MPS2 approximating the lower half of the PEPS sandwich:
   this->getMPS(PEPS0, "up", MPS2);
   for (int row = this->N[0]-2; row > this->N[0]/2; row--)
   {
    MPS0 = MPS2;
    this->getMPO(PEPS0, "up", row, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2s[row-2]);
    MPS2.setD(D1);
    eps = Epss[row-2];
    maxNumSweeps = MaxNumsSweeps[row-2];
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS2);
    ErrorsAchieved[row-2] = errorAchieved;
    NumsSweepsDone[row-2] = numSweepsDone;
   }
// return the contraction value:
   if (this->N[0] == 2)
   {
    return MPS2.contractReverse(MPS1);
   }
   else
   {
    this->getMPO(PEPS0, "down", this->N[0]/2, MPO0);
    return contractReverseMPOMPS(MPS2, MPO0, MPS1);
   }
  }
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented for periodic BC: " <<
          "template<class T> T PEPS<T>::" <<
          "expectationValue(const PEPO<T>& PEPO0, const string& Direction, const vector<unsigned int>& D2s, " <<
                           "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                           "vector<double>& ErrorsAchieved, " <<
                           "vector<unsigned int>& NumsSweepsDone) const." << endl;
  exit(1);
 }
}

template<class T> T PEPS<T>::expectationValue(const vector< vector<unsigned int> >& PositionsRow,
                                              const vector< vector<unsigned int> >& PositionsCol,
                                              const vector< vector< Matrix<T> > >& Interactions,
                                              const string& Direction, const vector<unsigned int>& D2s,
                                              const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps,
                                              vector<double>& ErrorsAchieved,
                                              vector<unsigned int>& NumsSweepsDone) const
{
 unsigned int numInteractions = Interactions.size();
#ifdef DEBUG
 if ((this->N[0] == 0) || (numInteractions == 0) || (PositionsRow.size() != numInteractions) ||
     (PositionsCol.size() != numInteractions) || ((Direction != "horizontal") && (Direction != "vertical")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "expectationValue(const vector< vector<unsigned int> >& PositionsRow, " <<
                           "const vector< vector<unsigned int> >& PositionsCol, " <<
                           "const vector< vector< Matrix<T> > >& Interactions, " <<
                           "const string& Direction, const vector<unsigned int>& D2s, " <<
                           "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                           "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
          "((this->N[0] == 0) || (Interactions.size() == 0) || (PositionsRow.size() != Interactions.size()) || " <<
           "(PositionsCol.size() != numInteractions) || " <<
           "((Direction != horizontal) && (Direction != vertical)))." << endl;
  exit(1);
 }
 if (((Direction == "horizontal") && (D2s.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (Epss.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (MaxNumsSweeps.size() != this->N[1]-3)) ||
     ((Direction == "horizontal") && (ErrorsAchieved.size() != (this->N[1]-3)*numInteractions)) ||
     ((Direction == "horizontal") && (NumsSweepsDone.size() != (this->N[1]-3)*numInteractions)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "expectationValue(const vector< vector<unsigned int> >& PositionsRow, " <<
                           "const vector< vector<unsigned int> >& PositionsCol, " <<
                           "const vector< vector< Matrix<T> > >& Interactions, " <<
                           "const string& Direction, const vector<unsigned int>& D2s, " <<
                           "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                           "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
          "(((Direction == horizontal) && (D2s.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (Epss.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (MaxNumsSweeps.size() != this->N[1]-3)) || " <<
           "((Direction == horizontal) && (ErrorsAchieved.size() != (this->N[1]-3)*Interactions.size())) || " <<
           "((Direction == horizontal) && (NumsSweepsDone.size() != (this->N[1]-3)*Interactions.size())))." << endl;
  exit(1);
 }
 if (((Direction == "vertical") && (D2s.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (Epss.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (MaxNumsSweeps.size() != this->N[0]-3)) ||
     ((Direction == "vertical") && (ErrorsAchieved.size() != (this->N[0]-3)*numInteractions)) ||
     ((Direction == "vertical") && (NumsSweepsDone.size() != (this->N[0]-3)*numInteractions)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "expectationValue(const vector< vector<unsigned int> >& PositionsRow, " <<
                           "const vector< vector<unsigned int> >& PositionsCol, " <<
                           "const vector< vector< Matrix<T> > >& Interactions, " <<
                           "const string& Direction, const vector<unsigned int>& D2s, " <<
                           "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                           "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
          "(((Direction == vertical) && (D2s.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (Epss.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (MaxNumsSweeps.size() != this->N[0]-3)) || " <<
           "((Direction == vertical) && (ErrorsAchieved.size() != (this->N[0]-3)*Interactions.size())) || " <<
           "((Direction == vertical) && (NumsSweepsDone.size() != (this->N[0]-3)*Interactions.size())))." << endl;
  exit(1);
 }
 for (int i = 0; i < D2s.size(); i++)
 {
  if (D2s[i] == 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> T PEPS<T>::" <<
           "expectationValue(const vector< vector<unsigned int> >& PositionsRow, " <<
                            "const vector< vector<unsigned int> >& PositionsCol, " <<
                            "const vector< vector< Matrix<T> > >& Interactions, " <<
                            "const string& Direction, const vector<unsigned int>& D2s, " <<
                            "const vector<double>& Epss, const vector<unsigned int>& MaxNumsSweeps, " <<
                            "vector<double>& ErrorsAchieved, vector<unsigned int>& NumsSweepsDone) const: " <<
           "(D2s[" << i << "] == 0)." << endl;
   exit(1);
  }
 }
#endif
 vector<double> ErrorsAchieved0; vector<unsigned int> NumsSweepsDone0;
 vector<T> NormsNorm, NormsExp;
 if (Direction == "horizontal")
 {
  ErrorsAchieved0 = vector<double>(this->N[1]-3);
  NumsSweepsDone0 = vector<unsigned int>(this->N[1]-3);
  NormsNorm = vector<T>(this->N[1]);
  NormsExp = vector<T>(this->N[1]);
 }
 else if (Direction == "vertical")
 {
  ErrorsAchieved0 = vector<double>(this->N[0]-3);
  NumsSweepsDone0 = vector<unsigned int>(this->N[0]-3);
  NormsNorm = vector<T>(this->N[0]);
  NormsExp = vector<T>(this->N[0]);
 }
// compute norm of this PEPS:
 T norm = this->scalarProduct(Direction, D2s, Epss, MaxNumsSweeps, ErrorsAchieved0, NumsSweepsDone0, NormsNorm);
// compute expectation values of Interactions:
 PEPS<T> PEPS0;
 unsigned int positionRow, positionCol;
 Matrix<T> Matrix0; Tensor<T> Tensor0;
 vector<unsigned int> Indices0(1), Indices1(1);
 Indices0[0] = 4; Indices1[0] = 1;
 T exp, result0, result = 0.0;
 for (int i = 0; i < numInteractions; i++)
 {
  PEPS0 = *this;
  for (int j = 0; j < Interactions[i].size(); j++)
  {
   positionRow = PositionsRow[i][j]; positionCol = PositionsCol[i][j];
   Matrix0 = Interactions[i][j];
   PEPS0.get(positionRow, positionCol, Tensor0);
   Tensor0.contract(Indices0, Matrix0, Indices1);
   PEPS0.set(positionRow, positionCol, Tensor0);
  }
  exp = this->scalarProduct(PEPS0, Direction, D2s, Epss, MaxNumsSweeps, ErrorsAchieved0, NumsSweepsDone0, NormsExp);
  result0 = 1.0;
  for (int j = 0; j < NormsExp.size(); j++)
   result0 *= NormsExp[j]/NormsNorm[j];
  result0 *= exp/norm;
  result += result0;
  for (int j = 0; j < ErrorsAchieved0.size(); j++)
  {
   ErrorsAchieved[i*ErrorsAchieved0.size()+j] = ErrorsAchieved0[j];
   NumsSweepsDone[i*ErrorsAchieved0.size()+j] = NumsSweepsDone0[j];
  }
 }
 return result;
}

template<class T> T PEPS<T>::expectationValue(const vector< vector<unsigned int> >& PositionsRow,
                                              const vector< vector<unsigned int> >& PositionsCol,
                                              const vector< vector< Matrix<T> > >& Interactions,
                                              const string& Direction,
                                              unsigned int D2, double eps, unsigned int maxNumSweeps,
                                              unsigned int M) const
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "expectationValue(const vector< vector<unsigned int> >& PositionsRow, " <<
                           "const vector< vector<unsigned int> >& PositionsCol, " <<
                           "const vector< vector< Matrix<T> > >& Interactions, " <<
                           "const string& Direction, unsigned int D2, " <<
                           "double eps, unsigned int maxNumSweeps, unsigned int M) const: " <<
          "(this->N[0] == 0)." << endl;
 }
 if ((Interactions.size() == 0) || (PositionsRow.size() != Interactions.size()) ||
     (PositionsCol.size() != Interactions.size()))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "expectationValue(const vector< vector<unsigned int> >& PositionsRow, " <<
                           "const vector< vector<unsigned int> >& PositionsCol, " <<
                           "const vector< vector< Matrix<T> > >& Interactions, " <<
                           "const string& Direction, unsigned int D2, " <<
                           "double eps, unsigned int maxNumSweeps, unsigned int M) const: " <<
          "((Interactions.size() == 0) || (PositionsRow.size() != Interactions.size()) || " <<
           "(PositionsCol.size() != Interactions.size()))." << endl;
  exit(1);
 }
 if ((Direction != "horizontal") && (Direction != "vertical"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "expectationValue(const vector< vector<unsigned int> >& PositionsRow, " <<
                           "const vector< vector<unsigned int> >& PositionsCol, " <<
                           "const vector< vector< Matrix<T> > >& Interactions, " <<
                           "const string& Direction, unsigned int D2, " <<
                           "double eps, unsigned int maxNumSweeps, unsigned int M) const: " <<
          "((Direction != horizontal) && (Direction != vertical))." << endl;
 }
 if ((D2 == 0) || (maxNumSweeps == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "expectationValue(const vector< vector<unsigned int> >& PositionsRow, " <<
                           "const vector< vector<unsigned int> >& PositionsCol, " <<
                           "const vector< vector< Matrix<T> > >& Interactions, " <<
                           "const string& Direction, unsigned int D2, " <<
                           "double eps, unsigned int maxNumSweeps, unsigned int M) const: " <<
          "((D2 == 0) || (maxNumSweeps == 0))." << endl;
 }
 if (M == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T PEPS<T>::" <<
          "expectationValue(const vector< vector<unsigned int> >& PositionsRow, " <<
                           "const vector< vector<unsigned int> >& PositionsCol, " <<
                           "const vector< vector< Matrix<T> > >& Interactions, " <<
                           "const string& Direction, unsigned int D2, " <<
                           "double eps, unsigned int maxNumSweeps, unsigned int M) const: " <<
          "(M == 0)." << endl;
 }
#endif
 unsigned int Ncols = this->N[1], Nrows = this->N[0], positionRow, positionCol;
 T exp, norm, result, result0;
 vector<unsigned int> Indices0(1), Indices1(1), NumsSweepsDone;
 vector<double> EpssAchieved;
 vector<T> NormsExp, NormsNorm;
 Tensor<T> Tensor0;
 Matrix<T> Matrix0;
 PEPS<T> PEPS0;
// 0. initialize:
 if (Direction == "horizontal")
 {
  NumsSweepsDone = vector<unsigned int>(Ncols-1);
  EpssAchieved = vector<double>(Ncols-1);
  NormsExp = vector<T>(Ncols-1);
  NormsNorm = vector<T>(Ncols-1);
 }
 else if (Direction == "vertical")
 {
  NumsSweepsDone = vector<unsigned int>(Nrows-1);
  EpssAchieved = vector<double>(Nrows-1);
  NormsExp = vector<T>(Nrows-1);
  NormsNorm = vector<T>(Nrows-1);
 }
// 1. compute norm-sandwich <thisPEPS|thisPEPS>:
 norm = this->scalarProduct(Direction, D2, eps, maxNumSweeps, EpssAchieved, NumsSweepsDone, NormsNorm);
// 2. compute expectation values of Interactions <thisPEPS|Interaction|thisPEPS>=:<thisPEPS|PEPS0>:
 result = 0.0;
 Indices0[0] = 4; Indices1[0] = 1;
 for (int i = 0; i < Interactions.size(); i++)
 {
// 2.1. construct |PEPS0>:=Interactions[i]*|thisPEPS>:
  PEPS0 = *this;
  for (int j = 0; j < Interactions[i].size(); j++)
  {
   positionRow = PositionsRow[i][j]*M;
   positionCol = PositionsCol[i][j]*M;
   Matrix0 = Interactions[i][j];
   PEPS0.Tensors[positionRow+Nrows*positionCol].contract(Indices0, Matrix0, Indices1);
  }
// 2.2. compute expectation value <thisPEPS|Interactions[i]|thisPEPS>=:<thisPEPS|PEPS0>:
  exp = this->scalarProduct(PEPS0, Direction, D2, eps, maxNumSweeps, EpssAchieved, NumsSweepsDone, NormsExp);
// 2.3. compute ratio <thisPEPS|Interactions[i]|thisPEPS>/<thisPEPS|thisPEPS>=:<thisPEPS|PEPS0>/<thisPEPS|thisPEPS>:
  result0 = exp/norm;
  for (int j = 0; j < NormsExp.size(); j++)
   result0 *= NormsExp[j]/NormsNorm[j];
// 2.4. add <thisPEPS|Interactions[i]|thisPEPS>/<thisPEPS|thisPEPS> to result:
  result += result0;
 }
// 3. return result:
 return result;
}

template<class T> double distancePEPOPEPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const PEPS<T>& PEPS1,
                                          const string& Direction, const vector<unsigned int>& D2s,
                                          const vector<double>& Epss,
                                          const vector<unsigned int>& MaxNumsSweeps, vector<double>& ErrorsAchieved,
                                          vector<unsigned int>& NumsSweepsDone)
{
#ifdef DEBUG
 if ((PEPO0.BC != PEPS0.BC) || (PEPO0.BC != PEPS1.BC) || (PEPO0.N != PEPS0.N) || (PEPO0.N != PEPS1.N) ||
     (PEPO0.N[0] == 0) || (PEPO0.d != PEPS0.d) || (PEPO0.d != PEPS1.d))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> double " <<
          "distancePEPOPEPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const PEPS<T>& PEPS1, " <<
                           "const string& Direction, const vector<unsigned int>& D2s, const vector<double>& Epss, " <<
                           "const vector<unsigned int>& MaxNumsSweeps, vector<double>& ErrorsAchieved, " <<
                           "vector<unsigned int>& NumsSweepsDone): " <<
          "PEPO0, PEPS0 and PEPS1 are not all of the same form or (PEPO0.N[0] == 0)." << endl;
  exit(1);
 }
 if ((Direction != "horizontal") && (Direction != "vertical"))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> double " <<
          "distancePEPOPEPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const PEPS<T>& PEPS1, " <<
                           "const string& Direction, const vector<unsigned int>& D2s, const vector<double>& Epss, " <<
                           "const vector<unsigned int>& MaxNumsSweeps, vector<double>& ErrorsAchieved, " <<
                           "vector<unsigned int>& NumsSweepsDone): " <<
          "((Direction != horizontal) && (Direction != vertical))." << endl;
  exit(1);
 }
 if (((Direction == "horizontal") && (D2s.size() != 3*(PEPO0.N[1]-3))) ||
     ((Direction == "horizontal") && (Epss.size() != 3*(PEPO0.N[1]-3))) ||
     ((Direction == "horizontal") && (MaxNumsSweeps.size() != 3*(PEPO0.N[1]-3))) ||
     ((Direction == "horizontal") && (ErrorsAchieved.size() != 3*(PEPO0.N[1]-3))) ||
     ((Direction == "horizontal") && (NumsSweepsDone.size() != 3*(PEPO0.N[1]-3))))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> double " <<
          "distancePEPOPEPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const PEPS<T>& PEPS1, " <<
                           "const string& Direction, const vector<unsigned int>& D2s, const vector<double>& Epss, " <<
                           "const vector<unsigned int>& MaxNumsSweeps, vector<double>& ErrorsAchieved, " <<
                           "vector<unsigned int>& NumsSweepsDone): " <<
          "(((Direction == horizontal) && (D2s.size() != 3*(PEPO0.N[1]-3))) || " <<
           "((Direction == horizontal) && (Epss.size() != 3*(PEPO0.N[1]-3))) || " <<
           "((Direction == horizontal) && (MaxNumsSweeps.size() != 3*(PEPO0.N[1]-3))) || " <<
           "((Direction == horizontal) && (ErrorsAchieved.size() != 3*(PEPO0.N[1]-3))) || " <<
           "((Direction == horizontal) && (NumsSweepsDone.size() != 3*(PEPO0.N[1]-3))))." << endl;
  exit(1);
 }
 if (((Direction == "vertical") && (D2s.size() != 3*(PEPO0.N[0]-3))) ||
     ((Direction == "vertical") && (Epss.size() != 3*(PEPO0.N[0]-3))) ||
     ((Direction == "vertical") && (MaxNumsSweeps.size() != 3*(PEPO0.N[0]-3))) ||
     ((Direction == "vertical") && (ErrorsAchieved.size() != 3*(PEPO0.N[0]-3))) ||
     ((Direction == "vertical") && (NumsSweepsDone.size() != 3*(PEPO0.N[0]-3))))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> double " <<
          "distancePEPOPEPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const PEPS<T>& PEPS1, " <<
                           "const string& Direction, const vector<unsigned int>& D2s, const vector<double>& Epss, " <<
                           "const vector<unsigned int>& MaxNumsSweeps, vector<double>& ErrorsAchieved, " <<
                           "vector<unsigned int>& NumsSweepsDone): " <<
          "(((Direction == vertical) && (D2s.size() != 3*(PEPO0.N[0]-3))) || " <<
           "((Direction == vertical) && (Epss.size() != 3*(PEPO0.N[0]-3))) || " <<
           "((Direction == vertical) && (MaxNumsSweeps.size() != 3*(PEPO0.N[0]-3))) || " <<
           "((Direction == vertical) && (ErrorsAchieved.size() != 3*(PEPO0.N[0]-3))) || " <<
           "((Direction == vertical) && (NumsSweepsDone.size() != 3*(PEPO0.N[0]-3))))." << endl;
  exit(1);
 }
 for (int i = 0; i < D2s.size(); i++)
 {
  if (D2s[i] == 0)
  {
   cerr << "Program terminated because of error in friend function " <<
           "template<class T> double " <<
           "distancePEPOPEPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const PEPS<T>& PEPS1, " <<
                            "const string& Direction, const vector<unsigned int>& D2s, const vector<double>& Epss, " <<
                            "const vector<unsigned int>& MaxNumsSweeps, vector<double>& ErrorsAchieved, " <<
                            "vector<unsigned int>& NumsSweepsDone): " <<
           "(D2s[" << i << "] == 0)." << endl;
   exit(1);
  }
 }
#endif
 unsigned int numMPOMPSMult;
 if (Direction == "horizontal")
  numMPOMPSMult = PEPO0.N[1]-3;
 else if (Direction == "vertical")
  numMPOMPSMult = PEPO0.N[0]-3;
 vector<unsigned int> D2s0(numMPOMPSMult), MaxNumsSweeps0(numMPOMPSMult), NumsSweepsDone0(numMPOMPSMult);
 vector<double> Epss0(numMPOMPSMult), ErrorsAchieved0(numMPOMPSMult);
// |PEPS2> = PEPO0*|PEPS0>:
 PEPS<T> PEPS2(PEPO0.BC, PEPO0.N[0], PEPO0.N[1], PEPO0.d, PEPO0.D*PEPS0.D);
 multiplyPEPOPEPSExact(PEPO0, PEPS0, PEPS2);
// normSquaredPEPO0PEPS0 = <PEPS0|PEPO0**{+}PEPO0|PEPS0>:
 for (int i = 0; i < numMPOMPSMult; i++)
 {
  D2s0[i] = D2s[i];
  Epss0[i] = Epss[i];
  MaxNumsSweeps0[i] = MaxNumsSweeps[i];
 }
 double normSquaredPEPO0PEPS0 = MathAuxiliary::convertToDouble(PEPS2.scalarProduct(Direction, D2s0, Epss0, MaxNumsSweeps0,
                                                                                   ErrorsAchieved0, NumsSweepsDone0));
 for (int i = 0; i < numMPOMPSMult; i++)
 {
  ErrorsAchieved[i] = ErrorsAchieved0[i];
  NumsSweepsDone[i] = NumsSweepsDone0[i];
 }
// realScalarProductPEPS1PEPO0PEPS0 = real(<PEPS1|PEPO0|PEPS0>):
 for (int i = 0; i < numMPOMPSMult; i++)
 {
  D2s0[i] = D2s[numMPOMPSMult+i];
  Epss0[i] = Epss[numMPOMPSMult+i];
  MaxNumsSweeps0[i] = MaxNumsSweeps[numMPOMPSMult+i];
 }
 double realScalarProductPEPS1PEPO0PEPS0 = MathAuxiliary::convertToDouble(PEPS1.scalarProduct(PEPS2, Direction, D2s0, Epss0, MaxNumsSweeps0,
                                                                                              ErrorsAchieved0, NumsSweepsDone0));
 for (int i = 0; i < numMPOMPSMult; i++)
 {
  ErrorsAchieved[numMPOMPSMult+i] = ErrorsAchieved0[i];
  NumsSweepsDone[numMPOMPSMult+i] = NumsSweepsDone0[i];
 }
// normSquaredPEPS1=<PEPS1|PEPS1>:
 for (int i = 0; i < numMPOMPSMult; i++)
 {
  D2s0[i] = D2s[2*numMPOMPSMult+i];
  Epss0[i] = Epss[2*numMPOMPSMult+i];
  MaxNumsSweeps0[i] = MaxNumsSweeps[2*numMPOMPSMult+i];
 }
 double normSquaredPEPS1 = MathAuxiliary::convertToDouble(PEPS1.scalarProduct(Direction, D2s0, Epss0, MaxNumsSweeps0, ErrorsAchieved0,
                                                                              NumsSweepsDone0));
 for (int i = 0; i < numMPOMPSMult; i++)
 {
  ErrorsAchieved[2*numMPOMPSMult+i] = ErrorsAchieved0[i];
  NumsSweepsDone[2*numMPOMPSMult+i] = NumsSweepsDone0[i];
 }
// return ||PEPO0*|PEPS0>-|PEPS1>||:
 return sqrt(abs(normSquaredPEPO0PEPS0 - 2.0*realScalarProductPEPS1PEPO0PEPS0 + normSquaredPEPS1));
}

template<class T> void multiplyPEPOPEPSExact(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, PEPS<T>& PEPS1)
{
#ifdef DEBUG
 if ((PEPO0.BC != PEPS0.BC) || (PEPO0.BC != PEPS1.BC) || (PEPO0.N != PEPS0.N) || (PEPO0.N != PEPS1.N) ||
     (PEPO0.N[0] == 0) || (PEPO0.d != PEPS0.d) || (PEPO0.d != PEPS1.d))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> void " <<
          "multiplyPEPOPEPSExact(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, PEPS<T>& PEPS1): " <<
          "((PEPO0.BC != PEPS0.BC) || (PEPO0.BC != PEPS1.BC) || (PEPO0.N != PEPS0.N) || (PEPO0.N != PEPS1.N) || " <<
           "(PEPO0.N[0] == 0) || (PEPO0.d != PEPS0.d) || (PEPO0.d != PEPS1.d))." << endl;
  exit(1);
 }
 if (PEPS1.D < PEPO0.D*PEPS0.D)
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> void " <<
          "multiplyPEPOPEPSExact(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, PEPS<T>& PEPS1): " <<
          "(PEPS1.D < PEPO0.D*PEPS0.D)." << endl;
  exit(1);
 }
#endif
 unsigned int DPEPS1 = PEPS1.D, Ncols = PEPO0.N[1], Nrows = PEPO0.N[0], position;
 vector<unsigned int> Indices0(1), Indices1(1), Order(9), Shape(5), Shape0(5), Shape1(6);
 Tensor<T> Tensor0, Tensor1;
 PEPS1.D = PEPS0.D*PEPO0.D;
 Indices0[0] = 4; Indices1[0] = 4;
 Order[0] = 0; Order[1] = 4; Order[2] = 1; Order[3] = 5; Order[4] = 2; Order[5] = 6; Order[6] = 3; Order[7] = 7; Order[8] = 8;
 for (int col = 0; col < Ncols; col++)
 {
  for (int row = 0; row < Nrows; row++)
  {
   position = row+Nrows*col;
   Tensor0 = PEPS0.Tensors[position]; Tensor1 = PEPO0.Tensors[position];
   Tensor0.getShape(Shape0); Tensor1.getShape(Shape1);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor0.permute(Order);
   Shape[0] = Shape0[0]*Shape1[0];
   Shape[1] = Shape0[1]*Shape1[1];
   Shape[2] = Shape0[2]*Shape1[2];
   Shape[3] = Shape0[3]*Shape1[3];
   Shape[4] = Shape1[5];
   Tensor0.reshape(Shape);
   PEPS1.Tensors[position] = Tensor0;
  }
 }
 PEPS1.setD(DPEPS1);
}

template<class T> void multiplyPEPOPEPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0,
                                        double eps, unsigned int maxNumSweeps,
                                        unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv,
                                        const string& UpdateTensor, const string& UpdateMode, double cutoff,
                                        double& epsAchieved, unsigned int& numSweepsDone,
                                        PEPS<T>& PEPS1)
{
#ifdef DEBUG
 if ((PEPO0.BC != PEPS0.BC) || (PEPO0.BC != PEPS1.BC) || (PEPO0.N != PEPS0.N) || (PEPO0.N != PEPS1.N) ||
     (PEPO0.N[0] == 0) || (PEPO0.d != PEPS0.d) || (PEPO0.d != PEPS1.d))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> void " <<
          "multiplyPEPOPEPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, double eps, unsigned int maxNumSweeps, " <<
                           "unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv, " <<
                           "const string& UpdateTensor, const string& UpdateMode, double cutoff, " <<
                           "double& epsAchieved, unsigned int& numSweepsDone, PEPS<T>& PEPS1): " <<
          "((PEPO0.BC != PEPS0.BC) || (PEPO0.BC != PEPS1.BC) || (PEPO0.N != PEPS0.N) || (PEPO0.N != PEPS1.N) || " <<
           "(PEPO0.N[0] == 0) || (PEPO0.d != PEPS0.d) || (PEPO0.d != PEPS1.d))." << endl;
  exit(1);
 }
 if ((UpdateTensor != "full") && (UpdateTensor != "reduced"))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> void " <<
          "multiplyPEPOPEPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, double eps, unsigned int maxNumSweeps, " <<
                           "unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv, " <<
                           "const string& UpdateTensor, const string& UpdateMode, double cutoff, " <<
                           "double& epsAchieved, unsigned int& numSweepsDone, PEPS<T>& PEPS1): " <<
          "((UpdateTensor != full) && (UpdateTensor != reduced))." << endl;
  exit(1);
 }
 if ((UpdateMode != "Pseudoinverse") && (UpdateMode != "PositiveNormMatrix"))
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> void " <<
          "multiplyPEPOPEPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, double eps, unsigned int maxNumSweeps, " <<
                           "unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv, " <<
                           "const string& UpdateTensor, const string& UpdateMode, double cutoff, " <<
                           "double& epsAchieved, unsigned int& numSweepsDone, PEPS<T>& PEPS1): " <<
          "((UpdateMode != Pseudoinverse) && (UpdateMode != PositiveNormMatrix))." << endl;
  exit(1);
 }
#endif
 unsigned int col, Ncols = PEPS1.N[1];
 string Boundary, SweepDirection;
 vector< MPS<T> > ColumnbMPSs, ColumnNormMPSs;
// a) PEPS1.D >= PEPO0.D*PEPS0.D:
 if (PEPS1.D >= PEPO0.D*PEPS0.D)
 {
  multiplyPEPOPEPSExact(PEPO0, PEPS0, PEPS1);
  epsAchieved = 0.0;
  numSweepsDone = 0;
  return;
 }
// b) PEPS1.D < PEPO0.D*PEPS0.D:
// b)0. initialize sweeping:
 PEPS1.getInitialColumnNormMPSs(D2Env, epsEnv, maxNumSweepsEnv, ColumnNormMPSs);
 PEPS1.getInitialColumnbMPSs(PEPO0, PEPS0, D2Env, epsEnv, maxNumSweepsEnv, ColumnbMPSs);
// b)1. sweep:
 numSweepsDone = 0;
 SweepDirection = "right";
 while (numSweepsDone < maxNumSweeps)
 {
// b)1.a) sweep right:
  if (SweepDirection == "right")
  {
// b)1.a)1. update column 0:
   Boundary = "left";
   PEPS1.updateBoundaryColumn(PEPO0, PEPS0, Boundary, eps, ColumnNormMPSs, ColumnbMPSs, UpdateTensor, UpdateMode, cutoff, epsAchieved);
   if (epsAchieved < eps)
    return;
   PEPS1.updateColumnNormMPSsBoundary(Boundary, D2Env, epsEnv, maxNumSweepsEnv, ColumnNormMPSs);
   PEPS1.updateColumnbMPSsBoundary(PEPO0, PEPS0, Boundary, D2Env, epsEnv, maxNumSweepsEnv, ColumnbMPSs);
// b)1.a)2. update columns 1 to Ncols-2:
   for (col = 1; col < Ncols-1; col++)
   {
    PEPS1.updateBulkColumn(PEPO0, PEPS0, col, eps, ColumnNormMPSs, ColumnbMPSs, UpdateTensor, UpdateMode, cutoff, epsAchieved);
    if (epsAchieved < eps)
     return;
    PEPS1.updateColumnNormMPSsBulk(col, SweepDirection, D2Env, epsEnv, maxNumSweepsEnv, ColumnNormMPSs);
    PEPS1.updateColumnbMPSsBulk(PEPO0, PEPS0, col, SweepDirection, D2Env, epsEnv, maxNumSweepsEnv, ColumnbMPSs);
   }
   SweepDirection = "left";
  }
// b)1.b) sweep left:
  else if (SweepDirection == "left")
  {
// b)1.b)1. update column Ncols-1:
   Boundary = "right";
   PEPS1.updateBoundaryColumn(PEPO0, PEPS0, Boundary, eps, ColumnNormMPSs, ColumnbMPSs, UpdateTensor, UpdateMode, cutoff, epsAchieved);
   if (epsAchieved < eps)
    return;
   PEPS1.updateColumnNormMPSsBoundary(Boundary, D2Env, epsEnv, maxNumSweepsEnv, ColumnNormMPSs);
   PEPS1.updateColumnbMPSsBoundary(PEPO0, PEPS0, Boundary, D2Env, epsEnv, maxNumSweepsEnv, ColumnbMPSs);
// b)1.b)2. update columns Ncols-2 to 1:
   for (col = Ncols-2; col > 0; col--)
   {
    PEPS1.updateBulkColumn(PEPO0, PEPS0, col, eps, ColumnNormMPSs, ColumnbMPSs, UpdateTensor, UpdateMode, cutoff, epsAchieved);
    if (epsAchieved < eps)
     return;
    PEPS1.updateColumnNormMPSsBulk(col, SweepDirection, D2Env, epsEnv, maxNumSweepsEnv, ColumnNormMPSs);
    PEPS1.updateColumnbMPSsBulk(PEPO0, PEPS0, col, SweepDirection, D2Env, epsEnv, maxNumSweepsEnv, ColumnbMPSs);
   }
   SweepDirection = "right";
  }
  numSweepsDone++;
 }
}

template<class T> void PEPS<T>::getMPS(const string& Direction, const vector<unsigned int>& dIndices, MPS<T>& MPS0) const
{
 MPS0 = MPS<T>(this->BC, this->N[0], this->D, this->D);
 Tensor<T> TensorPEPS, TensorMPS;
 vector<unsigned int> Indices(3), Index(2);
 if (Direction == "right")
 {
  Indices[0] = 3; Indices[1] = 1; Indices[2] = 2;
  Index[0] = 0;
  for (int position = 0; position < this->N[0]; position++)
  {
   TensorPEPS = this->Tensors[this->N[0]-1-position];
   Index[1] = dIndices[position];
   TensorPEPS.getSubtensor(Indices, Index, TensorMPS);
   MPS0.set(position, TensorMPS);
  }
 }
 else if (Direction == "left")
 {
  Indices[0] = 1; Indices[1] = 3; Indices[2] = 0;
  Index[0] = 0;
  for (int position = 0; position < this->N[0]; position++)
  {
   TensorPEPS = this->Tensors[position+this->N[0]*(this->N[1]-1)];
   Index[1] = dIndices[position];
   TensorPEPS.getSubtensor(Indices, Index, TensorMPS);
   MPS0.set(position, TensorMPS);
  }
 }
}

template<class T> void PEPS<T>::getMPO(const string& Direction, unsigned int positionCol, const vector<unsigned int>& dIndices, MPO<T>& MPO0) const
{
 MPO0 = MPO<T>(this->BC, this->N[0], this->D, this->D);
 Tensor<T> TensorPEPS, TensorMPO;
 vector<unsigned int> Indices(4), Index(1);
 if (Direction == "right")
 {
  Indices[0] = 3; Indices[1] = 1; Indices[2] = 0; Indices[3] = 2;
  for (int position = 0; position < this->N[0]; position++)
  {
   TensorPEPS = this->Tensors[this->N[0]-1-position+this->N[0]*positionCol];
   Index[0] = dIndices[position];
   TensorPEPS.getSubtensor(Indices, Index, TensorMPO);
   MPO0.set(position, TensorMPO);
  }
 }
 else if (Direction == "left")
 {
  Indices[0] = 1; Indices[1] = 3; Indices[2] = 2; Indices[3] = 0;
  for (int position = 0; position < this->N[0]; position++)
  {
   TensorPEPS = this->Tensors[position+this->N[0]*positionCol];
   Index[0] = dIndices[position];
   TensorPEPS.getSubtensor(Indices, Index, TensorMPO);
   MPO0.set(position, TensorMPO);
  }
 }
}

template<class T> void PEPS<T>::sampleExpectationValue(unsigned int positionRow, unsigned int positionCol,
                                                       const Matrix<T>& Observable,
                                                       unsigned int D1, const string& Filling, double eps, unsigned int maxNumSweeps,
                                                       unsigned int MEqu, unsigned int M, unsigned int MCSteps,
                                                       const string& FileName) const
{
// input parameters:
 double errorAchieved; unsigned int numSweepsDone;
 fstream File;
 File.open(FileName.c_str(), fstream::out);
 File.precision(15);
 srand(time(0));
 unsigned int row, col, Nrows = this->N[0], Ncols = this->N[1];
 Matrix<unsigned int> Sample(Nrows, Ncols);
// initialize sample:
 for (col = 0; col < Ncols; col++)
 {
  for (row = 0; row < Nrows; row++)
  {
   Sample(row, col) = rand()%2;
  }
 }
// initialize boundary MPSs:
 vector< MPS<T> > MPSs(Ncols);
 vector<unsigned int> dIndices(Nrows);
 MPS<T> MPS0, MPS1; MPO<T> MPO0;
 string Direction = "left";
 T coeff0, coeff1; double prob0, prob1;
 double randomNum;
 col = Ncols-1;
 for (row = 0; row < Nrows; row++)
  dIndices[row] = Sample(row, col);
 this->getMPS(Direction, dIndices, MPS1);
 MPSs[col] = MPS1;
 for (col = Ncols-2; col > 0; col--)
 {
  for (row = 0; row < Nrows; row++)
   dIndices[row] = Sample(row, col);
  this->getMPO(Direction, col, dIndices, MPO0);
  MPS0 = MPSs[col+1];
  MPS1 = MPS0;
  MPS1.setD(D1, Filling);
  multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
  MPSs[col] = MPS1;
 }
// equilibrate the Metropolis sampling:
 for (int Conf = 0; Conf < MEqu; Conf++)
 {
// Direction right:
  Direction = "right";
// col 0:
  col = 0;
  MPS1 = MPSs[col+1];
  for (row = 0; row < Nrows; row++)
  {
   for (int i = 0; i < Nrows; i++)
    dIndices[i] = Sample(Nrows-1-i, col);
   this->getMPS(Direction, dIndices, MPS0);
   coeff0 = MPS0.contractReverse(MPS1);
   prob0 = coeff0*MathAuxiliary::complexConjugate(coeff0);
   dIndices[row] = abs(int(dIndices[row])-1);
   this->getMPS(Direction, dIndices, MPS0);
   coeff1 = MPS0.contractReverse(MPS1);
   prob1 = coeff1*MathAuxiliary::complexConjugate(coeff1);
   randomNum = double(rand())/double(RAND_MAX+1.0);
   if (randomNum < prob1/prob0)
    Sample(Nrows-1-row, col) = dIndices[row];
  }
  for (row = 0; row < Nrows; row++)
   dIndices[row] = Sample(Nrows-1-row, col);
  this->getMPS(Direction, dIndices, MPS0);
  MPSs[col] = MPS0;
// 0 < col < Ncols-1:
  for (col = 1; col < Ncols-1; col++)
  {
   MPS0 = MPSs[col-1]; MPS1 = MPSs[col+1];
   for (row = 0; row < Nrows; row++)
   {
    for (int i = 0; i < Nrows; i++)
     dIndices[i] = Sample(Nrows-1-i, col);
    this->getMPO(Direction, col, dIndices, MPO0);
    coeff0 = contractReverseMPOMPS(MPS1, MPO0, MPS0);
    prob0 = coeff0*MathAuxiliary::complexConjugate(coeff0);
    dIndices[row] = abs(int(dIndices[row])-1);
    this->getMPO(Direction, col, dIndices, MPO0);
    coeff1 = contractReverseMPOMPS(MPS1, MPO0, MPS0);
    prob1 = coeff1*MathAuxiliary::complexConjugate(coeff1);
    randomNum = double(rand())/double(RAND_MAX+1.0);
    if (randomNum < prob1/prob0)
     Sample(Nrows-1-row, col) = dIndices[row];
   }
   for (row = 0; row < Nrows; row++)
    dIndices[row] = Sample(Nrows-1-row, col);
   this->getMPO(Direction, col, dIndices, MPO0);
   MPS1 = MPS0;
   MPS1.setD(D1, Filling);
   multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
   MPSs[col] = MPS1;
  }
// Direction left:
  Direction = "left";
// col Ncols-1:
  col = Ncols-1;
  MPS0 = MPSs[col-1];
  for (row = 0; row < Nrows; row++)
  {
   for (int i = 0; i < Nrows; i++)
    dIndices[i] = Sample(i, col);
   this->getMPS(Direction, dIndices, MPS1);
   coeff0 = MPS0.contractReverse(MPS1);
   prob0 = coeff0*MathAuxiliary::complexConjugate(coeff0);
   dIndices[row] = abs(int(dIndices[row])-1);
   this->getMPS(Direction, dIndices, MPS1);
   coeff1 = MPS0.contractReverse(MPS1);
   prob1 = coeff1*MathAuxiliary::complexConjugate(coeff1);
   randomNum = double(rand())/double(RAND_MAX+1.0);
   if (randomNum < prob1/prob0)
    Sample(row, col) = dIndices[row];
  }
  for (row = 0; row < Nrows; row++)
   dIndices[row] = Sample(row, col);
  this->getMPS(Direction, dIndices, MPS1);
  MPSs[col] = MPS1;
// 0 < col < Ncols-1:
  for (col = Ncols-2; col > 0; col--)
  {
   MPS0 = MPSs[col-1]; MPS1 = MPSs[col+1];
   for (row = 0; row < Nrows; row++)
   {
    for (int i = 0; i < Nrows; i++)
     dIndices[i] = Sample(i, col);
    this->getMPO(Direction, col, dIndices, MPO0);
    coeff0 = contractReverseMPOMPS(MPS0, MPO0, MPS1);
    prob0 = coeff0*MathAuxiliary::complexConjugate(coeff0);
    dIndices[row] = abs(int(dIndices[row])-1);
    this->getMPO(Direction, col, dIndices, MPO0);
    coeff1 = contractReverseMPOMPS(MPS0, MPO0, MPS1);
    prob1 = coeff1*MathAuxiliary::complexConjugate(coeff1);
    randomNum = double(rand())/double(RAND_MAX+1.0);
    if (randomNum < prob1/prob0)
     Sample(row, col) = dIndices[row];
   }
   for (row = 0; row < Nrows; row++)
    dIndices[row] = Sample(row, col);
   this->getMPO(Direction, col, dIndices, MPO0);
   MPS0 = MPS1;
   MPS0.setD(D1, Filling);
   multiplyMPOMPS(MPO0, MPS1, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS0);
   MPSs[col] = MPS0;
  }
 }
// compute observable:
 T result;
 unsigned int MCStep = 0;
 while (MCStep < MCSteps)
 {
  result = 0.0;
  for (int Conf = 0; Conf < M; Conf++)
  {
// Direction right:
   Direction = "right";
// col 0:
   col = 0;
   MPS1 = MPSs[col+1];
   for (row = 0; row < Nrows; row++)
   {
    for (int i = 0; i < Nrows; i++)
     dIndices[i] = Sample(Nrows-1-i, col);
    this->getMPS(Direction, dIndices, MPS0);
    coeff0 = MPS0.contractReverse(MPS1);
    prob0 = coeff0*MathAuxiliary::complexConjugate(coeff0);
    dIndices[row] = abs(int(dIndices[row])-1);
    this->getMPS(Direction, dIndices, MPS0);
    coeff1 = MPS0.contractReverse(MPS1);
    prob1 = coeff1*MathAuxiliary::complexConjugate(coeff1);
    randomNum = double(rand())/double(RAND_MAX+1.0);
    if (randomNum < prob1/prob0)
     Sample(Nrows-1-row, col) = dIndices[row];
   }
   for (row = 0; row < Nrows; row++)
    dIndices[row] = Sample(Nrows-1-row, col);
   this->getMPS(Direction, dIndices, MPS0);
   MPSs[col] = MPS0;
// 0 < col < Ncols-1:
   for (col = 1; col < Ncols-1; col++)
   {
    MPS0 = MPSs[col-1]; MPS1 = MPSs[col+1];
    for (row = 0; row < Nrows; row++)
    {
     for (int i = 0; i < Nrows; i++)
      dIndices[i] = Sample(Nrows-1-i, col);
     this->getMPO(Direction, col, dIndices, MPO0);
     coeff0 = contractReverseMPOMPS(MPS1, MPO0, MPS0);
     prob0 = coeff0*MathAuxiliary::complexConjugate(coeff0);
     dIndices[row] = abs(int(dIndices[row])-1);
     this->getMPO(Direction, col, dIndices, MPO0);
     coeff1 = contractReverseMPOMPS(MPS1, MPO0, MPS0);
     prob1 = coeff1*MathAuxiliary::complexConjugate(coeff1);
     randomNum = double(rand())/double(RAND_MAX+1.0);
     if (randomNum < prob1/prob0)
      Sample(Nrows-1-row, col) = dIndices[row];
    }
    for (row = 0; row < Nrows; row++)
     dIndices[row] = Sample(Nrows-1-row, col);
    this->getMPO(Direction, col, dIndices, MPO0);
    MPS1 = MPS0;
    MPS1.setD(D1, Filling);
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
    MPSs[col] = MPS1;
   }
// Direction left:
   Direction = "left";
// col Ncols-1:
   col = Ncols-1;
   MPS0 = MPSs[col-1];
   for (row = 0; row < Nrows; row++)
   {
    for (int i = 0; i < Nrows; i++)
     dIndices[i] = Sample(i, col);
    this->getMPS(Direction, dIndices, MPS1);
    coeff0 = MPS0.contractReverse(MPS1);
    prob0 = coeff0*MathAuxiliary::complexConjugate(coeff0);
    dIndices[row] = abs(int(dIndices[row])-1);
    this->getMPS(Direction, dIndices, MPS1);
    coeff1 = MPS0.contractReverse(MPS1);
    prob1 = coeff1*MathAuxiliary::complexConjugate(coeff1);
    randomNum = double(rand())/double(RAND_MAX+1.0);
    if (randomNum < prob1/prob0)
     Sample(row, col) = dIndices[row];
   }
   for (row = 0; row < Nrows; row++)
    dIndices[row] = Sample(row, col);
   this->getMPS(Direction, dIndices, MPS1);
   MPSs[col] = MPS1;
// 0 < col < Ncols-1:
   for (col = Ncols-2; col > 0; col--)
   {
    MPS0 = MPSs[col-1]; MPS1 = MPSs[col+1];
    for (row = 0; row < Nrows; row++)
    {
     for (int i = 0; i < Nrows; i++)
      dIndices[i] = Sample(i, col);
     this->getMPO(Direction, col, dIndices, MPO0);
     coeff0 = contractReverseMPOMPS(MPS0, MPO0, MPS1);
     prob0 = coeff0*MathAuxiliary::complexConjugate(coeff0);
     dIndices[row] = abs(int(dIndices[row])-1);
     this->getMPO(Direction, col, dIndices, MPO0);
     coeff1 = contractReverseMPOMPS(MPS0, MPO0, MPS1);
     prob1 = coeff1*MathAuxiliary::complexConjugate(coeff1);
     randomNum = double(rand())/double(RAND_MAX+1.0);
     if (randomNum < prob1/prob0)
      Sample(row, col) = dIndices[row];
    }
    for (row = 0; row < Nrows; row++)
     dIndices[row] = Sample(row, col);
    this->getMPO(Direction, col, dIndices, MPO0);
    MPS0 = MPS1;
    MPS0.setD(D1, Filling);
    multiplyMPOMPS(MPO0, MPS1, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS0);
    MPSs[col] = MPS0;
   }
// compute observable:
   result += Observable(Sample(positionRow, positionCol), Sample(positionRow, positionCol));
  }
  File << MCStep << "   " << result/double(M) << endl;
  MCStep++;
 }
 File.close();
}

template<class T> void PEPS<T>::sampleExpectationValue(const Matrix<T>& Observable,
                                                       unsigned int D1, const string& Filling, double eps, unsigned int maxNumSweeps,
                                                       unsigned int MEqu, unsigned int M, unsigned int MCSteps,
                                                       const string& FileName) const
{
// input parameters:
 double errorAchieved; unsigned int numSweepsDone;
 fstream File;
 File.open(FileName.c_str(), fstream::out);
 File.precision(15);
 srand(time(0));
 unsigned int row, col, Nrows = this->N[0], Ncols = this->N[1];
 Matrix<unsigned int> Sample(Nrows, Ncols);
// initialize sample:
 for (col = 0; col < Ncols; col++)
 {
  for (row = 0; row < Nrows; row++)
  {
   Sample(row, col) = rand()%2;
  }
 }
// initialize boundary MPSs:
 vector< MPS<T> > MPSs(Ncols);
 vector<unsigned int> dIndices(Nrows);
 MPS<T> MPS0, MPS1; MPO<T> MPO0;
 string Direction = "left";
 T coeff0, coeff1; double prob0, prob1;
 double randomNum;
 col = Ncols-1;
 for (row = 0; row < Nrows; row++)
  dIndices[row] = Sample(row, col);
 this->getMPS(Direction, dIndices, MPS1);
 MPSs[col] = MPS1;
 for (col = Ncols-2; col > 0; col--)
 {
  for (row = 0; row < Nrows; row++)
   dIndices[row] = Sample(row, col);
  this->getMPO(Direction, col, dIndices, MPO0);
  MPS0 = MPSs[col+1];
  MPS1 = MPS0;
  MPS1.setD(D1, Filling);
  multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
  MPSs[col] = MPS1;
 }
// equilibrate the Metropolis sampling:
 for (int Conf = 0; Conf < MEqu; Conf++)
 {
// Direction right:
  Direction = "right";
// col 0:
  col = 0;
  MPS1 = MPSs[col+1];
  for (row = 0; row < Nrows; row++)
  {
   for (int i = 0; i < Nrows; i++)
    dIndices[i] = Sample(Nrows-1-i, col);
   this->getMPS(Direction, dIndices, MPS0);
   coeff0 = MPS0.contractReverse(MPS1);
   prob0 = coeff0*MathAuxiliary::complexConjugate(coeff0);
   dIndices[row] = abs(int(dIndices[row])-1);
   this->getMPS(Direction, dIndices, MPS0);
   coeff1 = MPS0.contractReverse(MPS1);
   prob1 = coeff1*MathAuxiliary::complexConjugate(coeff1);
   randomNum = double(rand())/double(RAND_MAX+1.0);
   if (randomNum < prob1/prob0)
    Sample(Nrows-1-row, col) = dIndices[row];
  }
  for (row = 0; row < Nrows; row++)
   dIndices[row] = Sample(Nrows-1-row, col);
  this->getMPS(Direction, dIndices, MPS0);
  MPSs[col] = MPS0;
// 0 < col < Ncols-1:
  for (col = 1; col < Ncols-1; col++)
  {
   MPS0 = MPSs[col-1]; MPS1 = MPSs[col+1];
   for (row = 0; row < Nrows; row++)
   {
    for (int i = 0; i < Nrows; i++)
     dIndices[i] = Sample(Nrows-1-i, col);
    this->getMPO(Direction, col, dIndices, MPO0);
    coeff0 = contractReverseMPOMPS(MPS1, MPO0, MPS0);
    prob0 = coeff0*MathAuxiliary::complexConjugate(coeff0);
    dIndices[row] = abs(int(dIndices[row])-1);
    this->getMPO(Direction, col, dIndices, MPO0);
    coeff1 = contractReverseMPOMPS(MPS1, MPO0, MPS0);
    prob1 = coeff1*MathAuxiliary::complexConjugate(coeff1);
    randomNum = double(rand())/double(RAND_MAX+1.0);
    if (randomNum < prob1/prob0)
     Sample(Nrows-1-row, col) = dIndices[row];
   }
   for (row = 0; row < Nrows; row++)
    dIndices[row] = Sample(Nrows-1-row, col);
   this->getMPO(Direction, col, dIndices, MPO0);
   MPS1 = MPS0;
   MPS1.setD(D1, Filling);
   multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
   MPSs[col] = MPS1;
  }
// Direction left:
  Direction = "left";
// col Ncols-1:
  col = Ncols-1;
  MPS0 = MPSs[col-1];
  for (row = 0; row < Nrows; row++)
  {
   for (int i = 0; i < Nrows; i++)
    dIndices[i] = Sample(i, col);
   this->getMPS(Direction, dIndices, MPS1);
   coeff0 = MPS0.contractReverse(MPS1);
   prob0 = coeff0*MathAuxiliary::complexConjugate(coeff0);
   dIndices[row] = abs(int(dIndices[row])-1);
   this->getMPS(Direction, dIndices, MPS1);
   coeff1 = MPS0.contractReverse(MPS1);
   prob1 = coeff1*MathAuxiliary::complexConjugate(coeff1);
   randomNum = double(rand())/double(RAND_MAX+1.0);
   if (randomNum < prob1/prob0)
    Sample(row, col) = dIndices[row];
  }
  for (row = 0; row < Nrows; row++)
   dIndices[row] = Sample(row, col);
  this->getMPS(Direction, dIndices, MPS1);
  MPSs[col] = MPS1;
// 0 < col < Ncols-1:
  for (col = Ncols-2; col > 0; col--)
  {
   MPS0 = MPSs[col-1]; MPS1 = MPSs[col+1];
   for (row = 0; row < Nrows; row++)
   {
    for (int i = 0; i < Nrows; i++)
     dIndices[i] = Sample(i, col);
    this->getMPO(Direction, col, dIndices, MPO0);
    coeff0 = contractReverseMPOMPS(MPS0, MPO0, MPS1);
    prob0 = coeff0*MathAuxiliary::complexConjugate(coeff0);
    dIndices[row] = abs(int(dIndices[row])-1);
    this->getMPO(Direction, col, dIndices, MPO0);
    coeff1 = contractReverseMPOMPS(MPS0, MPO0, MPS1);
    prob1 = coeff1*MathAuxiliary::complexConjugate(coeff1);
    randomNum = double(rand())/double(RAND_MAX+1.0);
    if (randomNum < prob1/prob0)
     Sample(row, col) = dIndices[row];
   }
   for (row = 0; row < Nrows; row++)
    dIndices[row] = Sample(row, col);
   this->getMPO(Direction, col, dIndices, MPO0);
   MPS0 = MPS1;
   MPS0.setD(D1, Filling);
   multiplyMPOMPS(MPO0, MPS1, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS0);
   MPSs[col] = MPS0;
  }
 }
// compute observable:
 T result;
 unsigned int MCStep = 0;
 while (MCStep < MCSteps)
 {
  result = 0.0;
  for (int Conf = 0; Conf < M; Conf++)
  {
// Direction right:
   Direction = "right";
// col 0:
   col = 0;
   MPS1 = MPSs[col+1];
   for (row = 0; row < Nrows; row++)
   {
    for (int i = 0; i < Nrows; i++)
     dIndices[i] = Sample(Nrows-1-i, col);
    this->getMPS(Direction, dIndices, MPS0);
    coeff0 = MPS0.contractReverse(MPS1);
    prob0 = coeff0*MathAuxiliary::complexConjugate(coeff0);
    dIndices[row] = abs(int(dIndices[row])-1);
    this->getMPS(Direction, dIndices, MPS0);
    coeff1 = MPS0.contractReverse(MPS1);
    prob1 = coeff1*MathAuxiliary::complexConjugate(coeff1);
    randomNum = double(rand())/double(RAND_MAX+1.0);
    if (randomNum < prob1/prob0)
     Sample(Nrows-1-row, col) = dIndices[row];
   }
   for (row = 0; row < Nrows; row++)
    dIndices[row] = Sample(Nrows-1-row, col);
   this->getMPS(Direction, dIndices, MPS0);
   MPSs[col] = MPS0;
// 0 < col < Ncols-1:
   for (col = 1; col < Ncols-1; col++)
   {
    MPS0 = MPSs[col-1]; MPS1 = MPSs[col+1];
    for (row = 0; row < Nrows; row++)
    {
     for (int i = 0; i < Nrows; i++)
      dIndices[i] = Sample(Nrows-1-i, col);
     this->getMPO(Direction, col, dIndices, MPO0);
     coeff0 = contractReverseMPOMPS(MPS1, MPO0, MPS0);
     prob0 = coeff0*MathAuxiliary::complexConjugate(coeff0);
     dIndices[row] = abs(int(dIndices[row])-1);
     this->getMPO(Direction, col, dIndices, MPO0);
     coeff1 = contractReverseMPOMPS(MPS1, MPO0, MPS0);
     prob1 = coeff1*MathAuxiliary::complexConjugate(coeff1);
     randomNum = double(rand())/double(RAND_MAX+1.0);
     if (randomNum < prob1/prob0)
      Sample(Nrows-1-row, col) = dIndices[row];
    }
    for (row = 0; row < Nrows; row++)
     dIndices[row] = Sample(Nrows-1-row, col);
    this->getMPO(Direction, col, dIndices, MPO0);
    MPS1 = MPS0;
    MPS1.setD(D1, Filling);
    multiplyMPOMPS(MPO0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS1);
    MPSs[col] = MPS1;
   }
// Direction left:
   Direction = "left";
// col Ncols-1:
   col = Ncols-1;
   MPS0 = MPSs[col-1];
   for (row = 0; row < Nrows; row++)
   {
    for (int i = 0; i < Nrows; i++)
     dIndices[i] = Sample(i, col);
    this->getMPS(Direction, dIndices, MPS1);
    coeff0 = MPS0.contractReverse(MPS1);
    prob0 = coeff0*MathAuxiliary::complexConjugate(coeff0);
    dIndices[row] = abs(int(dIndices[row])-1);
    this->getMPS(Direction, dIndices, MPS1);
    coeff1 = MPS0.contractReverse(MPS1);
    prob1 = coeff1*MathAuxiliary::complexConjugate(coeff1);
    randomNum = double(rand())/double(RAND_MAX+1.0);
    if (randomNum < prob1/prob0)
     Sample(row, col) = dIndices[row];
   }
   for (row = 0; row < Nrows; row++)
    dIndices[row] = Sample(row, col);
   this->getMPS(Direction, dIndices, MPS1);
   MPSs[col] = MPS1;
// 0 < col < Ncols-1:
   for (col = Ncols-2; col > 0; col--)
   {
    MPS0 = MPSs[col-1]; MPS1 = MPSs[col+1];
    for (row = 0; row < Nrows; row++)
    {
     for (int i = 0; i < Nrows; i++)
      dIndices[i] = Sample(i, col);
     this->getMPO(Direction, col, dIndices, MPO0);
     coeff0 = contractReverseMPOMPS(MPS0, MPO0, MPS1);
     prob0 = coeff0*MathAuxiliary::complexConjugate(coeff0);
     dIndices[row] = abs(int(dIndices[row])-1);
     this->getMPO(Direction, col, dIndices, MPO0);
     coeff1 = contractReverseMPOMPS(MPS0, MPO0, MPS1);
     prob1 = coeff1*MathAuxiliary::complexConjugate(coeff1);
     randomNum = double(rand())/double(RAND_MAX+1.0);
     if (randomNum < prob1/prob0)
      Sample(row, col) = dIndices[row];
    }
    for (row = 0; row < Nrows; row++)
     dIndices[row] = Sample(row, col);
    this->getMPO(Direction, col, dIndices, MPO0);
    MPS0 = MPS1;
    MPS0.setD(D1, Filling);
    multiplyMPOMPS(MPO0, MPS1, eps, maxNumSweeps, errorAchieved, numSweepsDone, MPS0);
    MPSs[col] = MPS0;
   }
// compute observable:
   for (int j = 0; j < Ncols; j++)
   {
    for (int i = 0; i < Nrows; i++)
    {
     result += Observable(Sample(i, j), Sample(i, j))/(double(Nrows*Ncols));
    }
   }
  }
  File << MCStep << "   " << result/double(M) << endl;
  MCStep++;
 }
 File.close();
}

template<class T> void PEPS<T>::getOpenBCShape(unsigned int positionRow, unsigned int positionCol,
                                               vector<unsigned int>& Shape)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (this->BC != "open") || (positionRow > this->N[0]-1) ||
     (positionCol > this->N[1]-1) || (Shape.size() != 5))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "getOpenBCShape(unsigned int positionRow, unsigned int positionCol, " <<
                         "vector<unsigned int>& Shape) const: " <<
          "((this->N[0] == 0) || (this->BC != open) || (positionRow > this->N[0]-1) || " <<
           "(positionCol > this->N[1]-1) || (Shape.size() != 5))." << endl;
  exit(1);
 }
#endif
 Shape[0] = this->D; Shape[1] = this->D; Shape[2] = this->D; Shape[3] = this->D;
 Shape[4] = this->d(positionRow, positionCol);
 if (positionRow == 0)
  Shape[1] = 1;
 if (positionRow == this->N[0]-1)
  Shape[3] = 1;
 if (positionCol == 0)
  Shape[0] = 1;
 if (positionCol == this->N[1]-1)
  Shape[2] = 1;
}

template<class T> void PEPS<T>::getInitialColumnNormMPSs(unsigned int D2, double eps, unsigned int maxNumSweeps,
                                                         vector< MPS<T> >& ColumnNormMPSs) const
{
 unsigned int Ncols = this->N[1], Nrows = this->N[0], numSweepsDone;
 double epsAchieved;
 string BC = "open", Boundary = "right", Direction = "left";
 vector<unsigned int> dMPS(Nrows), dMPS1(Nrows), Seed(Nrows), Shape(5);
 ColumnNormMPSs = vector< MPS<T> >(Ncols);
// 1. boundary-MPS for right boundary:
 for (int i = 0; i < Nrows; i++)
 {
  this->Tensors[i+Nrows*(Ncols-1)].getShape(Shape);
  dMPS[i] = Shape[0]*Shape[0];
 }
 ColumnNormMPSs[Ncols-1] = MPS<T>(BC, Nrows, dMPS, min(D2, this->D*this->D));
 for (int i = 0; i < Nrows; i++)
  Seed[i] = time(0)+13*i;
 ColumnNormMPSs[Ncols-1].fillRandomly(Seed);
 this->getBoundaryMPS(Boundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, ColumnNormMPSs[Ncols-1]);
// 2. boundary-MPSs for bulk:
 for (int pos = Ncols-2; pos > 0; pos--)
 {
  for (int i = 0; i < Nrows; i++)
  {
   this->Tensors[i+Nrows*pos].getShape(Shape);
   dMPS1[i] = Shape[0]*Shape[0];
  }
  if (dMPS1 == dMPS)
  {
   ColumnNormMPSs[pos] = ColumnNormMPSs[pos+1];
   ColumnNormMPSs[pos].setD(min(D2, this->D*this->D*ColumnNormMPSs[pos+1].getD()));
  }
  else
  {
   dMPS = dMPS1;
   ColumnNormMPSs[pos] = MPS<T>(BC, Nrows, dMPS, min(D2, this->D*this->D*ColumnNormMPSs[pos+1].getD()));
   for (int i = 0; i < Nrows; i++)
    Seed[i] = time(0)+13*((Ncols-1-pos)*Nrows+i);
   ColumnNormMPSs[pos].fillRandomly(Seed);
  }
  this->multiplyBulkMPOMPS(Direction, pos, ColumnNormMPSs[pos+1], eps, maxNumSweeps, epsAchieved, numSweepsDone, ColumnNormMPSs[pos]);
 }
}

template<class T> void PEPS<T>::getInitialColumnbMPSs(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0,
                                                      unsigned int D2, double eps, unsigned int maxNumSweeps,
                                                      vector< MPS<T> >& ColumnbMPSs) const
{
 unsigned int Ncols = this->N[1], Nrows = this->N[0], numSweepsDone;
 double epsAchieved;
 string BC = "open", Boundary = "right", Direction = "left";
 vector<unsigned int> dMPS(Nrows), dMPS1(Nrows), Seed(Nrows), Shape(5), Shape0(5), ShapePEPO(6);
 Tensor<T> Tensor0;
 ColumnbMPSs = vector< MPS<T> >(Ncols);
// 1. boundary-MPS for right boundary:
 for (int i = 0; i < Nrows; i++)
 {
  PEPS0.Tensors[i+Nrows*(Ncols-1)].getShape(Shape0);
  PEPO0.get(i, Ncols-1, Tensor0); Tensor0.getShape(ShapePEPO);
  this->Tensors[i+Nrows*(Ncols-1)].getShape(Shape);
  dMPS[i] = Shape0[0]*ShapePEPO[0]*Shape[0];
 }
 ColumnbMPSs[Ncols-1] = MPS<T>(BC, Nrows, dMPS, min(D2, PEPS0.D*PEPO0.getD()*this->D));
 for (int i = 0; i < Nrows; i++)
  Seed[i] = time(0)+13*i;
 ColumnbMPSs[Ncols-1].fillRandomly(Seed);
 this->getBoundaryMPS(PEPO0, PEPS0, Boundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, ColumnbMPSs[Ncols-1]);
// 2. boundary-MPSs for bulk:
 for (int pos = Ncols-2; pos > 0; pos--)
 {
  for (int i = 0; i < Nrows; i++)
  {
   PEPS0.Tensors[i+Nrows*pos].getShape(Shape0);
   PEPO0.get(i, pos, Tensor0); Tensor0.getShape(ShapePEPO);
   this->Tensors[i+Nrows*pos].getShape(Shape);
   dMPS1[i] = Shape0[0]*ShapePEPO[0]*Shape[0];
  }
  if (dMPS1 == dMPS)
  {
   ColumnbMPSs[pos] = ColumnbMPSs[pos+1];
   ColumnbMPSs[pos].setD(min(D2, PEPS0.D*PEPO0.getD()*this->D*ColumnbMPSs[pos+1].getD()));
  }
  else
  {
   dMPS = dMPS1;
   ColumnbMPSs[pos] = MPS<T>(BC, Nrows, dMPS, min(D2, PEPS0.D*PEPO0.getD()*this->D*ColumnbMPSs[pos+1].getD()));
   for (int i = 0; i < Nrows; i++)
    Seed[i] = time(0)+13*((Ncols-1-pos)*Nrows+i);
   ColumnbMPSs[pos].fillRandomly(Seed);
  }
  this->multiplyBulkMPOMPS(PEPO0, PEPS0, Direction, pos, ColumnbMPSs[pos+1], eps, maxNumSweeps, epsAchieved, numSweepsDone, ColumnbMPSs[pos]);
 }
}

template<class T> void PEPS<T>::updateColumnNormMPSsBoundary(const string& Boundary,
                                                             unsigned int D2, double eps, unsigned int maxNumSweeps,
                                                             vector< MPS<T> >& ColumnNormMPSs) const
{
 unsigned int Ncols = this->N[1], Nrows = this->N[0], numSweepsDone;
 double epsAchieved;
 string BC = "open";
 vector<unsigned int> dMPS(Nrows), Seed(Nrows), Shape(5);
 for (int i = 0; i < Nrows; i++)
  Seed[i] = time(0)+13*i;
// a) left boundary:
 if (Boundary == "left")
 {
  for (int i = 0; i < Nrows; i++)
  {
   this->Tensors[Nrows-1-i].getShape(Shape);
   dMPS[i] = Shape[2]*Shape[2];
  }
  ColumnNormMPSs[0] = MPS<T>(BC, Nrows, dMPS, min(D2, this->D*this->D));
  ColumnNormMPSs[0].fillRandomly(Seed);
  this->getBoundaryMPS(Boundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, ColumnNormMPSs[0]);
 }
// b) right boundary:
 else if (Boundary == "right")
 {
  for (int i = 0; i < Nrows; i++)
  {
   this->Tensors[i+Nrows*(Ncols-1)].getShape(Shape);
   dMPS[i] = Shape[0]*Shape[0];
  }
  ColumnNormMPSs[Ncols-1] = MPS<T>(BC, Nrows, dMPS, min(D2, this->D*this->D));
  ColumnNormMPSs[Ncols-1].fillRandomly(Seed);
  this->getBoundaryMPS(Boundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, ColumnNormMPSs[Ncols-1]);
 }
}

template<class T> void PEPS<T>::updateColumnbMPSsBoundary(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary,
                                                          unsigned int D2, double eps, unsigned int maxNumSweeps,
                                                          vector< MPS<T> >& ColumnbMPSs) const
{
 unsigned int Ncols = this->N[1], Nrows = this->N[0], numSweepsDone;
 double epsAchieved;
 string BC = "open";
 vector<unsigned int> dMPS(Nrows), Seed(Nrows), Shape(5), Shape0(5), ShapePEPO(6);
 Tensor<T> Tensor0;
 for (int i = 0; i < Nrows; i++)
  Seed[i] = time(0)+13*i;
// left boundary:
 if (Boundary == "left")
 {
  for (int i = 0; i < Nrows; i++)
  {
   PEPS0.Tensors[Nrows-1-i].getShape(Shape0);
   PEPO0.get(Nrows-1-i, 0, Tensor0); Tensor0.getShape(ShapePEPO);
   this->Tensors[Nrows-1-i].getShape(Shape);
   dMPS[i] = Shape0[2]*ShapePEPO[2]*Shape[2];
  }
  ColumnbMPSs[0] = MPS<T>(BC, Nrows, dMPS, min(D2, PEPS0.D*PEPO0.getD()*this->D));
  ColumnbMPSs[0].fillRandomly(Seed);
  this->getBoundaryMPS(PEPO0, PEPS0, Boundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, ColumnbMPSs[0]);
 }
// right boundary:
 else if (Boundary == "right")
 {
  for (int i = 0; i < Nrows; i++)
  {
   PEPS0.Tensors[i+Nrows*(Ncols-1)].getShape(Shape0);
   PEPO0.get(i, Ncols-1, Tensor0); Tensor0.getShape(ShapePEPO);
   this->Tensors[i+Nrows*(Ncols-1)].getShape(Shape);
   dMPS[i] = Shape0[0]*ShapePEPO[0]*Shape[0];
  }
  ColumnbMPSs[Ncols-1] = MPS<T>(BC, Nrows, dMPS, min(D2, PEPS0.D*PEPO0.getD()*this->D));
  ColumnbMPSs[Ncols-1].fillRandomly(Seed);
  this->getBoundaryMPS(PEPO0, PEPS0, Boundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, ColumnbMPSs[Ncols-1]);
 }
}

template<class T> void PEPS<T>::updateColumnNormMPSsBulk(unsigned int position, const string& Direction,
                                                         unsigned int D2, double eps, unsigned int maxNumSweeps,
                                                         vector< MPS<T> >& ColumnNormMPSs) const
{
 unsigned int Ncols = this->N[1], Nrows = this->N[0], numSweepsDone;
 double epsAchieved;
 string BC = "open";
 vector<unsigned int> dMPS(Nrows), dMPS1(Nrows), Seed, Shape(5);
// a) direction right:
 if (Direction == "right")
 {
  ColumnNormMPSs[position-1].getd(dMPS);
  for (int i = 0; i < Nrows; i++)
  {
   this->Tensors[Nrows-1-i+Nrows*position].getShape(Shape);
   dMPS1[i] = Shape[2]*Shape[2];
  }
  if (dMPS1 == dMPS)
  {
   ColumnNormMPSs[position] = ColumnNormMPSs[position-1];
   ColumnNormMPSs[position].setD(min(D2, this->D*this->D*ColumnNormMPSs[position-1].getD()));
  }
  else
  {
   ColumnNormMPSs[position] = MPS<T>(BC, Nrows, dMPS1, min(D2, this->D*this->D*ColumnNormMPSs[position-1].getD()));
   Seed = vector<unsigned int>(Nrows);
   for (int i = 0; i < Nrows; i++)
    Seed[i] = time(0)+13*(position*Nrows+i);
   ColumnNormMPSs[position].fillRandomly(Seed);
  }
  this->multiplyBulkMPOMPS(Direction, position, ColumnNormMPSs[position-1], eps, maxNumSweeps, epsAchieved, numSweepsDone, ColumnNormMPSs[position]);
 }
// b) direction left:
 else if (Direction == "left")
 {
  ColumnNormMPSs[position+1].getd(dMPS);
  for (int i = 0; i < Nrows; i++)
  {
   this->Tensors[i+Nrows*position].getShape(Shape);
   dMPS1[i] = Shape[0]*Shape[0];
  }
  if (dMPS1 == dMPS)
  {
   ColumnNormMPSs[position] = ColumnNormMPSs[position+1];
   ColumnNormMPSs[position].setD(min(D2, this->D*this->D*ColumnNormMPSs[position+1].getD()));
  }
  else
  {
   ColumnNormMPSs[position] = MPS<T>(BC, Nrows, dMPS1, min(D2, this->D*this->D*ColumnNormMPSs[position+1].getD()));
   Seed = vector<unsigned int>(Nrows);
   for (int i = 0; i < Nrows; i++)
    Seed[i] = time(0)+13*(position*Nrows+i);
   ColumnNormMPSs[position].fillRandomly(Seed);
  }
  this->multiplyBulkMPOMPS(Direction, position, ColumnNormMPSs[position+1], eps, maxNumSweeps, epsAchieved, numSweepsDone, ColumnNormMPSs[position]);
 }
}

template<class T> void PEPS<T>::updateColumnbMPSsBulk(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, unsigned int position, const string& Direction,
                                                      unsigned int D2, double eps, unsigned int maxNumSweeps,
                                                      vector< MPS<T> >& ColumnbMPSs) const
{
 unsigned int Ncols = this->N[1], Nrows = this->N[0], numSweepsDone;
 double epsAchieved;
 string BC = "open";
 vector<unsigned int> dMPS(Nrows), dMPS1(Nrows), Seed, Shape(5), Shape0(5), ShapePEPO(6);
 Tensor<T> Tensor0;
// a) direction right:
 if (Direction == "right")
 {
  ColumnbMPSs[position-1].getd(dMPS);
  for (int i = 0; i < Nrows; i++)
  {
   PEPS0.Tensors[Nrows-1-i+Nrows*position].getShape(Shape0);
   PEPO0.get(Nrows-1-i, position, Tensor0); Tensor0.getShape(ShapePEPO);
   this->Tensors[Nrows-1-i+Nrows*position].getShape(Shape);
   dMPS1[i] = Shape0[2]*ShapePEPO[2]*Shape[2];
  }
  if (dMPS1 == dMPS)
  {
   ColumnbMPSs[position] = ColumnbMPSs[position-1];
   ColumnbMPSs[position].setD(min(D2, PEPS0.D*PEPO0.getD()*this->D*ColumnbMPSs[position-1].getD()));
  }
  else
  {
   ColumnbMPSs[position] = MPS<T>(BC, Nrows, dMPS1, min(D2, PEPS0.D*PEPO0.getD()*this->D*ColumnbMPSs[position-1].getD()));
   Seed = vector<unsigned int>(Nrows);
   for (int i = 0; i < Nrows; i++)
    Seed[i] = time(0)+13*(position*Nrows+i);
   ColumnbMPSs[position].fillRandomly(Seed);
  }
  this->multiplyBulkMPOMPS(PEPO0, PEPS0, Direction, position, ColumnbMPSs[position-1], eps, maxNumSweeps, epsAchieved, numSweepsDone, ColumnbMPSs[position]);
 }
// b) direction left:
 else if (Direction == "left")
 {
  ColumnbMPSs[position+1].getd(dMPS);
  for (int i = 0; i < Nrows; i++)
  {
   PEPS0.Tensors[i+Nrows*position].getShape(Shape0);
   PEPO0.get(i, position, Tensor0); Tensor0.getShape(ShapePEPO);
   this->Tensors[i+Nrows*position].getShape(Shape);
   dMPS1[i] = Shape0[0]*ShapePEPO[0]*Shape[0];
  }
  if (dMPS1 == dMPS)
  {
   ColumnbMPSs[position] = ColumnbMPSs[position+1];
   ColumnbMPSs[position].setD(min(D2, PEPS0.D*PEPO0.getD()*this->D*ColumnbMPSs[position+1].getD()));
  }
  else
  {
   ColumnbMPSs[position] = MPS<T>(BC, Nrows, dMPS1, min(D2, PEPS0.D*PEPO0.getD()*this->D*ColumnbMPSs[position+1].getD()));
   Seed = vector<unsigned int>(Nrows);
   for (int i = 0; i < Nrows; i++)
    Seed[i] = time(0)+13*(position*Nrows+i);
   ColumnbMPSs[position].fillRandomly(Seed);
  }
  this->multiplyBulkMPOMPS(PEPO0, PEPS0, Direction, position, ColumnbMPSs[position+1], eps, maxNumSweeps, epsAchieved, numSweepsDone, ColumnbMPSs[position]);
 }
}

template<class T> void PEPS<T>::getInitialColumnNormTensorsBoundary(const string& Boundary, const MPS<T>& ColumnNormMPS,
                                                                    vector< Tensor<T> >& ColumnNormTensors) const
{
 unsigned int L = this->N[0], Ncols = this->N[1], Nrows = this->N[0];
 vector<unsigned int> Indices0(1), Indices02(2), Indices03(3), Indices1(1), Indices12(2), Indices13(3);
 vector<unsigned int> Shape3(3), Shape4(4), Shape5(5), Shape6(6), Shape8(8);
 Tensor<T> Tensor0, Tensor1;
 ColumnNormTensors = vector< Tensor<T> >(L);
// a) left boundary:
 if (Boundary == "left")
 {
// a)1. position L-1:
  ColumnNormMPS.get(L-1, Tensor0);
  Tensor0.getShape(Shape3);
  this->Tensors[L-1].getShape(Shape5);
  Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[2];
  Tensor0.reshape(Shape4);
  Tensor1 = this->Tensors[L-1];
  Indices0[0] = 2; Indices1[0] = 2;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  this->Tensors[L-1].complexConjugate(Tensor1);
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 2; Indices12[1] = 4;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape8);
  Shape4[0] = Shape8[0]; Shape4[1] = 1; Shape4[2] = Shape8[3]; Shape4[3] = Shape8[3];
  Tensor0.reshape(Shape4);
  ColumnNormTensors[L-1] = Tensor0;
// a)2. L-1 > position > 0:
  Indices0[0] = 0; Indices1[0] = 1;
  Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 3; Indices12[1] = 2;
  Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 3; Indices13[1] = 2; Indices13[2] = 4;
  for (int pos = L-2; pos > 0; pos--)
  {
   ColumnNormMPS.get(pos, Tensor1);
   Tensor1.getShape(Shape3);
   this->Tensors[pos].getShape(Shape5);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[2];
   Tensor1.reshape(Shape4);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = this->Tensors[pos];
   Tensor0.contract(Indices02, Tensor1, Indices12);
   this->Tensors[pos].complexConjugate(Tensor1);
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = Shape6[1]; Shape4[1] = 1; Shape4[2] = Shape6[3]; Shape4[3] = Shape6[3];
   Tensor0.reshape(Shape4);
   ColumnNormTensors[pos] = Tensor0;
  }
 }
// b) right boundary:
 else if (Boundary == "right")
 {
// b)1. position L-1:
  ColumnNormMPS.get(0, Tensor0);
  Tensor0.getShape(Shape3);
  this->Tensors[L-1+Nrows*(Ncols-1)].getShape(Shape5);
  Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[0]; Shape4[3] = Shape5[0];
  Tensor0.reshape(Shape4);
  Tensor1 = this->Tensors[L-1+Nrows*(Ncols-1)];
  Indices0[0] = 2; Indices1[0] = 0;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  this->Tensors[L-1+Nrows*(Ncols-1)].complexConjugate(Tensor1);
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 4;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape8);
  Shape4[0] = 1; Shape4[1] = Shape8[1]; Shape4[2] = Shape8[2]; Shape4[3] = Shape8[2];
  Tensor0.reshape(Shape4);
  ColumnNormTensors[L-1] = Tensor0;
// b)2. L-1 > position > 0:
  Indices0[0] = 1; Indices1[0] = 0;
  Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 3; Indices12[1] = 0;
  Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 3; Indices13[1] = 0; Indices13[2] = 4;
  for (int pos = L-2; pos > 0; pos--)
  {
   ColumnNormMPS.get(L-1-pos, Tensor1);
   Tensor1.getShape(Shape3);
   this->Tensors[pos+Nrows*(Ncols-1)].getShape(Shape5);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[0]; Shape4[3] = Shape5[0];
   Tensor1.reshape(Shape4);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = this->Tensors[pos+Nrows*(Ncols-1)];
   Tensor0.contract(Indices02, Tensor1, Indices12);
   this->Tensors[pos+Nrows*(Ncols-1)].complexConjugate(Tensor1);
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = 1; Shape4[1] = Shape6[1]; Shape4[2] = Shape6[2]; Shape4[3] = Shape6[2];
   Tensor0.reshape(Shape4);
   ColumnNormTensors[pos] = Tensor0;
  }
 }
}

template<class T> void PEPS<T>::getInitialColumnbTensorsBoundary(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary,
                                                                 const MPS<T>& ColumnbMPS,
                                                                 vector< Tensor<T> >& ColumnbTensors) const
{
 unsigned int L = this->N[0], Ncols = this->N[1], Nrows = this->N[0];
 vector<unsigned int> Indices0(1), Indices02(2), Indices03(3), Indices1(1), Indices12(2), Indices13(3);
 vector<unsigned int> Shape3(3), Shape5(5), Shape05(5), Shape15(5), Shape6(6), Shape8(8), Shape11(11);
 Tensor<T> Tensor0, Tensor1, TensorPEPO;
 ColumnbTensors = vector< Tensor<T> >(L);
// a) left boundary:
 if (Boundary == "left")
 {
// a)1. position L-1:
  ColumnbMPS.get(L-1, Tensor0);
  Tensor0.getShape(Shape3);
  PEPS0.Tensors[L-1].getShape(Shape05);
  PEPO0.get(L-1, 0, TensorPEPO); TensorPEPO.getShape(Shape6);
  this->Tensors[L-1].getShape(Shape5);
  Shape15[0] = Shape3[0]; Shape15[1] = 1; Shape15[2] = Shape05[2]; Shape15[3] = Shape6[2]; Shape15[4] = Shape5[2];
  Tensor0.reshape(Shape15);
  Tensor1 = PEPS0.Tensors[L-1];
  Indices0[0] = 2; Indices1[0] = 2;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Indices02[0] = 2; Indices02[1] = 7; Indices12[0] = 2; Indices12[1] = 4;
  Tensor0.contract(Indices02, TensorPEPO, Indices12);
  this->Tensors[L-1].complexConjugate(Tensor1);
  Indices02[1] = 9;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape11);
  Shape5[0] = Shape11[0]; Shape5[1] = 1; Shape5[2] = Shape11[3]; Shape5[3] = Shape11[6]; Shape5[4] = Shape11[9];
  Tensor0.reshape(Shape5);
  ColumnbTensors[L-1] = Tensor0;
// a)2. L-1 > position > 0:
  Indices0[0] = 0; Indices1[0] = 1;
  Indices02[0] = 1; Indices02[1] = 5; Indices12[0] = 3; Indices12[1] = 2;
  Indices03[0] = 1; Indices03[2] = 8; Indices13[0] = 3; Indices13[1] = 2; Indices13[2] = 4;
  for (int pos = L-2; pos > 0; pos--)
  {
   ColumnbMPS.get(pos, Tensor1);
   Tensor1.getShape(Shape3);
   PEPS0.Tensors[pos].getShape(Shape05);
   PEPO0.get(pos, 0, TensorPEPO); TensorPEPO.getShape(Shape6);
   this->Tensors[pos].getShape(Shape5);
   Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[2]; Shape15[3] = Shape6[2]; Shape15[4] = Shape5[2];
   Tensor1.reshape(Shape15);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = PEPS0.Tensors[pos];
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Indices03[1] = 4;
   Tensor0.contract(Indices03, TensorPEPO, Indices13);
   this->Tensors[pos].complexConjugate(Tensor1);
   Indices03[1] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape8);
   Shape5[0] = Shape8[1]; Shape5[1] = 1; Shape5[2] = Shape8[3]; Shape5[3] = Shape8[5]; Shape5[4] = Shape8[7];
   Tensor0.reshape(Shape5);
   ColumnbTensors[pos] = Tensor0;
  }
 }
// b) right boundary:
 else if (Boundary == "right")
 {
// b)1. position L-1:
  ColumnbMPS.get(0, Tensor0);
  Tensor0.getShape(Shape3);
  PEPS0.Tensors[L-1+Nrows*(Ncols-1)].getShape(Shape05);
  PEPO0.get(L-1, Ncols-1, TensorPEPO); TensorPEPO.getShape(Shape6);
  this->Tensors[L-1+Nrows*(Ncols-1)].getShape(Shape5);
  Shape15[0] = 1; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[0]; Shape15[3] = Shape6[0]; Shape15[4] = Shape5[0];
  Tensor0.reshape(Shape15);
  Tensor1 = PEPS0.Tensors[L-1+Nrows*(Ncols-1)];
  Indices0[0] = 2; Indices1[0] = 0;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Indices02[0] = 2; Indices02[1] = 7; Indices12[0] = 0; Indices12[1] = 4;
  Tensor0.contract(Indices02, TensorPEPO, Indices12);
  this->Tensors[L-1+Nrows*(Ncols-1)].complexConjugate(Tensor1);
  Indices02[1] = 9;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape11);
  Shape5[0] = 1; Shape5[1] = Shape11[1]; Shape5[2] = Shape11[2]; Shape5[3] = Shape11[5]; Shape5[4] = Shape11[8];
  Tensor0.reshape(Shape5);
  ColumnbTensors[L-1] = Tensor0;
// b)2. L-1 > position > 0:
  Indices0[0] = 1; Indices1[0] = 0;
  Indices02[0] = 1; Indices02[1] = 5; Indices12[0] = 3; Indices12[1] = 0;
  Indices03[0] = 1; Indices03[2] = 8; Indices13[0] = 3; Indices13[1] = 0; Indices13[2] = 4;
  for (int pos = L-2; pos > 0; pos--)
  {
   ColumnbMPS.get(L-1-pos, Tensor1);
   Tensor1.getShape(Shape3);
   PEPS0.Tensors[pos+Nrows*(Ncols-1)].getShape(Shape05);
   PEPO0.get(pos, Ncols-1, TensorPEPO); TensorPEPO.getShape(Shape6);
   this->Tensors[pos+Nrows*(Ncols-1)].getShape(Shape5);
   Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[0]; Shape15[3] = Shape6[0]; Shape15[4] = Shape5[0];
   Tensor1.reshape(Shape15);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = PEPS0.Tensors[pos+Nrows*(Ncols-1)];
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Indices03[1] = 4;
   Tensor0.contract(Indices03, TensorPEPO, Indices13);
   this->Tensors[pos+Nrows*(Ncols-1)].complexConjugate(Tensor1);
   Indices03[1] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape8);
   Shape5[0] = 1; Shape5[1] = Shape8[1]; Shape5[2] = Shape8[2]; Shape5[3] = Shape8[4]; Shape5[4] = Shape8[6];
   Tensor0.reshape(Shape5);
   ColumnbTensors[pos] = Tensor0;
  }
 }
}

template<class T> void PEPS<T>::getInitialColumnNormTensorsBulk(unsigned int position, const MPS<T>& ColumnNormMPSL, const MPS<T>& ColumnNormMPSR,
                                                                vector< Tensor<T> >& ColumnNormTensors) const
{
 unsigned int L = this->N[0], Nrows = this->N[0];
 vector<unsigned int> Indices0(1), Indices02(2), Indices03(3), Indices1(1), Indices12(2), Indices13(3);
 vector<unsigned int> Order4(4), Shape3(3), Shape4(4), Shape5(5), Shape8(8);
 Tensor<T> Tensor0, Tensor1;
 ColumnNormTensors = vector< Tensor<T> >(L);
// 1. position L-1:
 ColumnNormMPSL.get(0, Tensor0);
 Tensor0.getShape(Shape3);
 this->Tensors[L-1+Nrows*position].getShape(Shape5);
 Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[0]; Shape4[3] = Shape5[0];
 Tensor0.reshape(Shape4);
 Tensor1 = this->Tensors[L-1+Nrows*position];
 Indices0[0] = 2; Indices1[0] = 0;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 this->Tensors[L-1+Nrows*position].complexConjugate(Tensor1);
 Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 4;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 ColumnNormMPSR.get(L-1, Tensor1);
 Tensor1.getShape(Shape3);
 Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[2];
 Tensor1.reshape(Shape4);
 Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 2; Indices12[1] = 3;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 Tensor0.getShape(Shape8);
 Shape4[0] = Shape8[1]; Shape4[1] = Shape8[2]; Shape4[2] = Shape8[4]; Shape4[3] = Shape8[6];
 Tensor0.reshape(Shape4);
 Order4[0] = 3; Order4[1] = 0; Order4[2] = 1; Order4[3] = 2;
 Tensor0.permute(Order4);
 ColumnNormTensors[L-1] = Tensor0;
// b)2. L-1 > position > 0:
 Indices0[0] = 1; Indices1[0] = 0;
 Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 3; Indices12[1] = 0;
 Indices03[1] = 3;
 for (int pos = L-2; pos > 0; pos--)
 {
  ColumnNormMPSL.get(L-1-pos, Tensor1);
  Tensor1.getShape(Shape3);
  this->Tensors[pos+Nrows*position].getShape(Shape5);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[0]; Shape4[3] = Shape5[0];
  Tensor1.reshape(Shape4);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor1 = this->Tensors[pos+Nrows*position];
  Tensor0.contract(Indices02, Tensor1, Indices12);
  this->Tensors[pos+Nrows*position].complexConjugate(Tensor1);
  Indices03[0] = 1; Indices03[2] = 6; Indices13[0] = 3; Indices13[1] = 0; Indices13[2] = 4;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  ColumnNormMPSR.get(pos, Tensor1);
  Tensor1.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[2];
  Tensor1.reshape(Shape4);
  Indices03[0] = 0; Indices03[2] = 5; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  Tensor0.permute(Order4);
  ColumnNormTensors[pos] = Tensor0;
 }
}

template<class T> void PEPS<T>::getInitialColumnbTensorsBulk(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, unsigned int position,
                                                             const MPS<T>& ColumnbMPSL, const MPS<T>& ColumnbMPSR,
                                                             vector< Tensor<T> >& ColumnbTensors) const
{
 unsigned int L = this->N[0], Nrows = this->N[0];
 vector<unsigned int> Indices0(1), Indices02(2), Indices03(3), Indices04(4), Indices1(1), Indices12(2), Indices13(3), Indices14(4);
 vector<unsigned int> Order5(5), Shape3(3), Shape5(5), Shape05(5), Shape15(5), Shape6(6), Shape10(10);
 Tensor<T> Tensor0, Tensor1, TensorPEPO;
 ColumnbTensors = vector< Tensor<T> >(L);
// 1. position L-1:
 ColumnbMPSL.get(0, Tensor0);
 Tensor0.getShape(Shape3);
 PEPS0.Tensors[L-1+Nrows*position].getShape(Shape05);
 PEPO0.get(L-1, position, TensorPEPO); TensorPEPO.getShape(Shape6);
 this->Tensors[L-1+Nrows*position].getShape(Shape5);
 Shape15[0] = 1; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[0]; Shape15[3] = Shape6[0]; Shape15[4] = Shape5[0];
 Tensor0.reshape(Shape15);
 Tensor1 = PEPS0.Tensors[L-1+Nrows*position];
 Indices0[0] = 2; Indices1[0] = 0;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Indices02[0] = 2; Indices02[1] = 7; Indices12[0] = 0; Indices12[1] = 4;
 Tensor0.contract(Indices02, TensorPEPO, Indices12);
 this->Tensors[L-1+Nrows*position].complexConjugate(Tensor1);
 Indices02[1] = 9;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 ColumnbMPSR.get(L-1, Tensor1);
 Tensor1.getShape(Shape3);
 Shape15[0] = Shape3[0]; Shape15[1] = 1; Shape15[2] = Shape05[2]; Shape15[3] = Shape6[2]; Shape15[4] = Shape5[2];
 Tensor1.reshape(Shape15);
 Indices03[0] = 3; Indices03[1] = 6; Indices03[2] = 9; Indices13[0] = 2; Indices13[1] = 3; Indices13[2] = 4;
 Tensor0.contract(Indices03, Tensor1, Indices13);
 Tensor0.getShape(Shape10);
 Shape5[0] = Shape10[1]; Shape5[1] = Shape10[2]; Shape5[2] = Shape10[4]; Shape5[3] = Shape10[6]; Shape5[4] = Shape10[8];
 Tensor0.reshape(Shape5);
 Order5[0] = 4; Order5[1] = 0; Order5[2] = 1; Order5[3] = 2; Order5[4] = 3;
 Tensor0.permute(Order5);
 ColumnbTensors[L-1] = Tensor0;
// b)2. L-1 > position > 0:
 Indices0[0] = 1; Indices1[0] = 0;
 Indices02[0] = 1; Indices02[1] = 5; Indices12[0] = 3; Indices12[1] = 0;
 Indices03[0] = 1; Indices03[2] = 8; Indices13[0] = 3; Indices13[1] = 0; Indices13[2] = 4;
 Indices04[0] = 0; Indices04[1] = 3; Indices04[2] = 5; Indices04[3] = 7; Indices14[0] = 1; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
 for (int pos = L-2; pos > 0; pos--)
 {
  ColumnbMPSL.get(L-1-pos, Tensor1);
  Tensor1.getShape(Shape3);
  PEPS0.Tensors[pos+Nrows*position].getShape(Shape05);
  PEPO0.get(pos, position, TensorPEPO); TensorPEPO.getShape(Shape6);
  this->Tensors[pos+Nrows*position].getShape(Shape5);
  Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[0]; Shape15[3] = Shape6[0]; Shape15[4] = Shape5[0];
  Tensor1.reshape(Shape15);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor1 = PEPS0.Tensors[pos+Nrows*position];
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Indices03[1] = 4;
  Tensor0.contract(Indices03, TensorPEPO, Indices13);
  this->Tensors[pos+Nrows*position].complexConjugate(Tensor1);
  Indices03[1] = 3;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  ColumnbMPSR.get(pos, Tensor1);
  Tensor1.getShape(Shape3);
  Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[2]; Shape15[3] = Shape6[2]; Shape15[4] = Shape5[2];
  Tensor1.reshape(Shape15);
  Tensor0.contract(Indices04, Tensor1, Indices14);
  Tensor0.permute(Order5);
  ColumnbTensors[pos] = Tensor0;
 }
}

template<class T> void PEPS<T>::updateColumnNormTensorsBoundary(const string& Boundary, unsigned int tensorPosition, const MPS<T>& ColumnNormMPS,
                                                                vector< Tensor<T> >& ColumnNormTensors) const
{
 unsigned int L = this->N[0], Ncols = this->N[1], Nrows = this->N[0];
 vector<unsigned int> Indices0(1), Indices02(2), Indices03(3), Indices1(1), Indices12(2), Indices13(3);
 vector<unsigned int> Shape3(3), Shape4(4), Shape5(5), Shape6(6), Shape8(8);
 Tensor<T> Tensor0, Tensor1;
// a) left boundary:
 if (Boundary == "left")
 {
// a)a) tensor-position 0:
  if (tensorPosition == 0)
  {
   ColumnNormMPS.get(0, Tensor0);
   Tensor0.getShape(Shape3);
   this->Tensors[0].getShape(Shape5);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[2];
   Tensor0.reshape(Shape4);
   Tensor1 = this->Tensors[0];
   Indices0[0] = 2; Indices1[0] = 2;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   this->Tensors[0].complexConjugate(Tensor1);
   Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 2; Indices12[1] = 4;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape8);
   Shape4[0] = 1; Shape4[1] = Shape8[1]; Shape4[2] = Shape8[4]; Shape4[3] = Shape8[4];
   Tensor0.reshape(Shape4);
   ColumnNormTensors[0] = Tensor0;
  }
// a)b) 0 < tensor-position < L-1:
  else if ((tensorPosition > 0) && (tensorPosition < L-1))
  {
   Tensor0 = ColumnNormTensors[tensorPosition-1];
   ColumnNormMPS.get(tensorPosition, Tensor1);
   Tensor1.getShape(Shape3);
   this->Tensors[tensorPosition].getShape(Shape5);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[2];
   Tensor1.reshape(Shape4);
   Indices0[0] = 1; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = this->Tensors[tensorPosition];
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   this->Tensors[tensorPosition].complexConjugate(Tensor1);
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = 1; Shape4[1] = Shape6[1]; Shape4[2] = Shape6[3]; Shape4[3] = Shape6[3];
   Tensor0.reshape(Shape4);
   ColumnNormTensors[tensorPosition] = Tensor0;
  }
 }
// b) right boundary:
 else if (Boundary == "right")
 {
// b)a) tensor-position 0:
  if (tensorPosition == 0)
  {
   ColumnNormMPS.get(L-1, Tensor0);
   Tensor0.getShape(Shape3);
   this->Tensors[Nrows*(Ncols-1)].getShape(Shape5);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape5[0]; Shape4[3] = Shape5[0];
   Tensor0.reshape(Shape4);
   Tensor1 = this->Tensors[Nrows*(Ncols-1)];
   Indices0[0] = 2; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   this->Tensors[Nrows*(Ncols-1)].complexConjugate(Tensor1);
   Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 4;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape8);
   Shape4[0] = Shape8[0]; Shape4[1] = 1; Shape4[2] = Shape8[4]; Shape4[3] = Shape8[4];
   Tensor0.reshape(Shape4);
   ColumnNormTensors[0] = Tensor0;
  }
// b)b) 0 < tensor-position < L-1:
  else if ((tensorPosition > 0) && (tensorPosition < L-1))
  {
   Tensor0 = ColumnNormTensors[tensorPosition-1];
   ColumnNormMPS.get(L-1-tensorPosition, Tensor1);
   Tensor1.getShape(Shape3);
   this->Tensors[tensorPosition+Nrows*(Ncols-1)].getShape(Shape5);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[0]; Shape4[3] = Shape5[0];
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = this->Tensors[tensorPosition+Nrows*(Ncols-1)];
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 0;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   this->Tensors[tensorPosition+Nrows*(Ncols-1)].complexConjugate(Tensor1);
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 1; Indices13[1] = 0; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = Shape6[1]; Shape4[1] = 1; Shape4[2] = Shape6[3]; Shape4[3] = Shape6[3];
   Tensor0.reshape(Shape4);
   ColumnNormTensors[tensorPosition] = Tensor0;
  }
 }
}

template<class T> void PEPS<T>::updateColumnbTensorsBoundary(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary,
                                                             unsigned int tensorPosition, const MPS<T>& ColumnbMPS,
                                                             vector< Tensor<T> >& ColumnbTensors) const
{
 unsigned int L = this->N[0], Ncols = this->N[1], Nrows = this->N[0];
 vector<unsigned int> Indices0(1), Indices02(2), Indices03(3), Indices1(1), Indices12(2), Indices13(3);
 vector<unsigned int> Shape3(3), Shape5(5), Shape05(5), Shape15(5), Shape6(6), Shape8(8), Shape11(11);
 Tensor<T> Tensor0, Tensor1, TensorPEPO;
// a) left boundary:
 if (Boundary == "left")
 {
// a)a) tensor-position 0:
  if (tensorPosition == 0)
  {
   ColumnbMPS.get(0, Tensor0);
   Tensor0.getShape(Shape3);
   PEPS0.Tensors[0].getShape(Shape05);
   PEPO0.get(0, 0, TensorPEPO); TensorPEPO.getShape(Shape6);
   this->Tensors[0].getShape(Shape5);
   Shape15[0] = 1; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[2]; Shape15[3] = Shape6[2]; Shape15[4] = Shape5[2];
   Tensor0.reshape(Shape15);
   Tensor1 = PEPS0.Tensors[0];
   Indices0[0] = 2; Indices1[0] = 2;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Indices02[0] = 2; Indices02[1] = 7; Indices12[0] = 2; Indices12[1] = 4;
   Tensor0.contract(Indices02, TensorPEPO, Indices12);
   this->Tensors[0].complexConjugate(Tensor1);
   Indices02[1] = 9;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape11);
   Shape5[0] = 1; Shape5[1] = Shape11[1]; Shape5[2] = Shape11[4]; Shape5[3] = Shape11[7]; Shape5[4] = Shape11[10];
   Tensor0.reshape(Shape5);
   ColumnbTensors[0] = Tensor0;
  }
// a)b) 0 < tensor-position < L-1:
  else if ((tensorPosition > 0) && (tensorPosition < L-1))
  {
   Tensor0 = ColumnbTensors[tensorPosition-1];
   ColumnbMPS.get(tensorPosition, Tensor1);
   Tensor1.getShape(Shape3);
   PEPS0.Tensors[tensorPosition].getShape(Shape05);
   PEPO0.get(tensorPosition, 0, TensorPEPO); TensorPEPO.getShape(Shape6);
   this->Tensors[tensorPosition].getShape(Shape5);
   Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[2]; Shape15[3] = Shape6[2]; Shape15[4] = Shape5[2];
   Tensor1.reshape(Shape15);
   Indices0[0] = 1; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = PEPS0.Tensors[tensorPosition];
   Indices02[0] = 1; Indices02[1] = 5; Indices12[0] = 1; Indices12[1] = 2;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Indices03[0] = 1; Indices03[1] = 4; Indices03[2] = 8; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 4;
   Tensor0.contract(Indices03, TensorPEPO, Indices13);
   this->Tensors[tensorPosition].complexConjugate(Tensor1);
   Indices03[1] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape8);
   Shape5[0] = 1; Shape5[1] = Shape8[1]; Shape5[2] = Shape8[3]; Shape5[3] = Shape8[5]; Shape5[4] = Shape8[7];
   Tensor0.reshape(Shape5);
   ColumnbTensors[tensorPosition] = Tensor0;
  }
 }
// b) right boundary:
 else if (Boundary == "right")
 {
// b)a) tensor-position 0:
  if (tensorPosition == 0)
  {
   ColumnbMPS.get(L-1, Tensor0);
   Tensor0.getShape(Shape3);
   PEPS0.Tensors[Nrows*(Ncols-1)].getShape(Shape05);
   PEPO0.get(0, Ncols-1, TensorPEPO); TensorPEPO.getShape(Shape6);
   this->Tensors[Nrows*(Ncols-1)].getShape(Shape5);
   Shape15[0] = Shape3[0]; Shape15[1] = 1; Shape15[2] = Shape05[0]; Shape15[3] = Shape6[0]; Shape15[4] = Shape5[0];
   Tensor0.reshape(Shape15);
   Tensor1 = PEPS0.Tensors[Nrows*(Ncols-1)];
   Indices0[0] = 2; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Indices02[0] = 2; Indices02[1] = 7; Indices12[0] = 0; Indices12[1] = 4;
   Tensor0.contract(Indices02, TensorPEPO, Indices12);
   this->Tensors[Nrows*(Ncols-1)].complexConjugate(Tensor1);
   Indices02[1] = 9;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape11);
   Shape5[0] = Shape11[0]; Shape5[1] = 1; Shape5[2] = Shape11[4]; Shape5[3] = Shape11[7]; Shape5[4] = Shape11[10];
   Tensor0.reshape(Shape5);
   ColumnbTensors[0] = Tensor0;
  }
// b)b) 0 < tensor-position < L-1:
  else if ((tensorPosition > 0) && (tensorPosition < L-1))
  {
   Tensor0 = ColumnbTensors[tensorPosition-1];
   ColumnbMPS.get(L-1-tensorPosition, Tensor1);
   Tensor1.getShape(Shape3);
   PEPS0.Tensors[tensorPosition+Nrows*(Ncols-1)].getShape(Shape05);
   PEPO0.get(tensorPosition, Ncols-1, TensorPEPO); TensorPEPO.getShape(Shape6);
   this->Tensors[tensorPosition+Nrows*(Ncols-1)].getShape(Shape5);
   Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[0]; Shape15[3] = Shape6[0]; Shape15[4] = Shape5[0];
   Tensor1.reshape(Shape15);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Tensor1 = PEPS0.Tensors[tensorPosition+Nrows*(Ncols-1)];
   Indices02[0] = 1; Indices02[1] = 5; Indices12[0] = 1; Indices12[1] = 0;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Indices03[0] = 1; Indices03[1] = 4; Indices03[2] = 8; Indices13[0] = 1; Indices13[1] = 0; Indices13[2] = 4;
   Tensor0.contract(Indices03, TensorPEPO, Indices13);
   this->Tensors[tensorPosition+Nrows*(Ncols-1)].complexConjugate(Tensor1);
   Indices03[1] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape8);
   Shape5[0] = Shape8[1]; Shape5[1] = 1; Shape5[2] = Shape8[3]; Shape5[3] = Shape8[5]; Shape5[4] = Shape8[7];
   Tensor0.reshape(Shape5);
   ColumnbTensors[tensorPosition] = Tensor0;
  }
 }
}

template<class T> void PEPS<T>::updateColumnNormTensorsBulk(unsigned int position, unsigned int tensorPosition,
                                                            const MPS<T>& ColumnNormMPSL, const MPS<T>& ColumnNormMPSR,
                                                            vector< Tensor<T> >& ColumnNormTensors) const
{
 unsigned int L = this->N[0], Nrows = this->N[0];
 vector<unsigned int> Indices0(1), Indices02(2), Indices03(3), Indices1(1), Indices12(2), Indices13(3);
 vector<unsigned int> Order4(4), Shape3(3), Shape4(4), Shape5(5), Shape8(8);
 Tensor<T> Tensor0, Tensor1;
// a) tensor-position 0:
 if (tensorPosition == 0)
 {
  ColumnNormMPSL.get(L-1, Tensor0);
  Tensor0.getShape(Shape3);
  this->Tensors[Nrows*position].getShape(Shape5);
  Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape5[0]; Shape4[3] = Shape5[0];
  Tensor0.reshape(Shape4);
  Tensor1 = this->Tensors[Nrows*position];
  Indices0[0] = 2; Indices1[0] = 0;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  this->Tensors[Nrows*position].complexConjugate(Tensor1);
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 4;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  ColumnNormMPSR.get(0, Tensor1);
  Tensor1.getShape(Shape3);
  Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[2];
  Tensor1.reshape(Shape4);
  Indices02[0] = 3; Indices12[0] = 2; Indices12[1] = 3;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape8);
  Shape4[0] = Shape8[0]; Shape4[1] = Shape8[3]; Shape4[2] = Shape8[5]; Shape4[3] = Shape8[7];
  Tensor0.reshape(Shape4);
  Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
  Tensor0.permute(Order4);
  ColumnNormTensors[0] = Tensor0;
 }
// b) 0 < tensor-position < L-1:
 else if ((tensorPosition > 0) && (tensorPosition < L-1))
 {
  Tensor0 = ColumnNormTensors[tensorPosition-1];
  ColumnNormMPSL.get(L-1-tensorPosition, Tensor1);
  Tensor1.getShape(Shape3);
  this->Tensors[tensorPosition+Nrows*position].getShape(Shape5);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[0]; Shape4[3] = Shape5[0];
  Tensor1.reshape(Shape4);
  Indices0[0] = 0; Indices1[0] = 1;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor1 = this->Tensors[tensorPosition+Nrows*position];
  Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 0;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  this->Tensors[tensorPosition+Nrows*position].complexConjugate(Tensor1);
  Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 1; Indices13[1] = 0; Indices13[2] = 4;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  ColumnNormMPSR.get(tensorPosition, Tensor1);
  Tensor1.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[2];
  Tensor1.reshape(Shape4);
  Indices03[0] = 0; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
  Tensor0.permute(Order4);
  ColumnNormTensors[tensorPosition] = Tensor0;
 }
}

template<class T> void PEPS<T>::updateColumnbTensorsBulk(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, unsigned int position,
                                                         unsigned int tensorPosition, const MPS<T>& ColumnbMPSL, const MPS<T>& ColumnbMPSR,
                                                         vector< Tensor<T> >& ColumnbTensors) const
{
 unsigned int L = this->N[0], Nrows = this->N[0];
 vector<unsigned int> Indices0(1), Indices02(2), Indices03(3), Indices04(4), Indices1(1), Indices12(2), Indices13(3), Indices14(4);
 vector<unsigned int> Order5(5), Shape3(3), Shape5(5), Shape05(5), Shape15(5), Shape6(6), Shape10(10);
 Tensor<T> Tensor0, Tensor1, TensorPEPO;
// a) tensor-position 0:
 if (tensorPosition == 0)
 {
  ColumnbMPSL.get(L-1, Tensor0);
  Tensor0.getShape(Shape3);
  PEPS0.Tensors[Nrows*position].getShape(Shape05);
  PEPO0.get(0, position, TensorPEPO); TensorPEPO.getShape(Shape6);
  this->Tensors[Nrows*position].getShape(Shape5);
  Shape15[0] = Shape3[0]; Shape15[1] = 1; Shape15[2] = Shape05[0]; Shape15[3] = Shape6[0]; Shape15[4] = Shape5[0];
  Tensor0.reshape(Shape15);
  Tensor1 = PEPS0.Tensors[Nrows*position];
  Indices0[0] = 2; Indices1[0] = 0;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Indices02[0] = 2; Indices02[1] = 7; Indices12[0] = 0; Indices12[1] = 4;
  Tensor0.contract(Indices02, TensorPEPO, Indices12);
  this->Tensors[Nrows*position].complexConjugate(Tensor1);
  Indices02[1] = 9;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  ColumnbMPSR.get(0, Tensor1);
  Tensor1.getShape(Shape3);
  Shape15[0] = 1; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[2]; Shape15[3] = Shape6[2]; Shape15[4] = Shape5[2];
  Tensor1.reshape(Shape15);
  Indices03[0] = 3; Indices03[1] = 6; Indices03[2] = 9; Indices13[0] = 2; Indices13[1] = 3; Indices13[2] = 4;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  Tensor0.getShape(Shape10);
  Shape5[0] = Shape10[0]; Shape5[1] = Shape10[3]; Shape5[2] = Shape10[5]; Shape5[3] = Shape10[7]; Shape5[4] = Shape10[9];
  Tensor0.reshape(Shape5);
  Order5[0] = 0; Order5[1] = 4; Order5[2] = 1; Order5[3] = 2; Order5[4] = 3;
  Tensor0.permute(Order5);
  ColumnbTensors[0] = Tensor0;
 }
// b) 0 < tensor-position < L-1:
 else if ((tensorPosition > 0) && (tensorPosition < L-1))
 {
  Tensor0 = ColumnbTensors[tensorPosition-1];
  ColumnbMPSL.get(L-1-tensorPosition, Tensor1);
  Tensor1.getShape(Shape3);
  PEPS0.Tensors[tensorPosition+Nrows*position].getShape(Shape05);
  PEPO0.get(tensorPosition, position, TensorPEPO); TensorPEPO.getShape(Shape6);
  this->Tensors[tensorPosition+Nrows*position].getShape(Shape5);
  Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[0]; Shape15[3] = Shape6[0]; Shape15[4] = Shape5[0];
  Tensor1.reshape(Shape15);
  Indices0[0] = 0; Indices1[0] = 1;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor1 = PEPS0.Tensors[tensorPosition+Nrows*position];
  Indices02[0] = 1; Indices02[1] = 5; Indices12[0] = 1; Indices12[1] = 0;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Indices03[0] = 1; Indices03[1] = 4; Indices03[2] = 8; Indices13[0] = 1; Indices13[1] = 0; Indices13[2] = 4;
  Tensor0.contract(Indices03, TensorPEPO, Indices13);
  this->Tensors[tensorPosition+Nrows*position].complexConjugate(Tensor1);
  Indices03[1] = 3;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  ColumnbMPSR.get(tensorPosition, Tensor1);
  Tensor1.getShape(Shape3);
  Shape15[0] = Shape3[0]; Shape15[1] = Shape3[1]; Shape15[2] = Shape05[2]; Shape15[3] = Shape6[2]; Shape15[4] = Shape5[2];
  Tensor1.reshape(Shape15);
  Indices04[0] = 0; Indices04[1] = 2; Indices04[2] = 4; Indices04[3] = 6; Indices14[0] = 0; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
  Tensor0.contract(Indices04, Tensor1, Indices14);
  Order5[0] = 0; Order5[1] = 4; Order5[2] = 1; Order5[3] = 2; Order5[4] = 3;
  Tensor0.permute(Order5);  
  ColumnbTensors[tensorPosition] = Tensor0;
 }
}

template<class T> void PEPS<T>::getColumnNormMPOBoundary(const string& Boundary, unsigned int tensorPosition,
                                                         const MPS<T>& ColumnNormMPS, const vector< Tensor<T> >& ColumnNormTensors,
                                                         MPO<T>& ColumnNormMPO) const
{
 unsigned int L = this->N[0], Ncols = this->N[1], Nrows = this->N[0];
 vector<unsigned int> Shape3(3), Shape4(4), Shape5(5);
 Tensor<T> OneTensor, Tensor0;
// OneTensor:
 Shape4[0] = 1; Shape4[1] = 1; Shape4[2] = 1; Shape4[3] = 1;
 OneTensor = Tensor<T>(Shape4);
 OneTensor.set(0, 1.0);
// ColumnNormMPO:
 ColumnNormMPO = MPO<T>("periodic", 4, this->D, ColumnNormMPS.getD());
// a) left boundary:
 if (Boundary == "left")
 {
// a)a) tensor-position 0:
  if (tensorPosition == 0)
  {
   ColumnNormMPO.set(0, OneTensor);
   ColumnNormMPO.set(1, OneTensor);
   ColumnNormMPS.get(0, Tensor0);
   Tensor0.getShape(Shape3);
   this->Tensors[0].getShape(Shape5);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[2];
   Tensor0.reshape(Shape4);
   ColumnNormMPO.set(2, Tensor0);
   ColumnNormMPO.set(3, ColumnNormTensors[1]);
  }
// a)b) 0 < tensor-position < L-1:
  else if ((tensorPosition > 0) && (tensorPosition < L-1))
  {
   ColumnNormMPO.set(0, OneTensor);
   ColumnNormMPO.set(1, ColumnNormTensors[tensorPosition-1]);
   ColumnNormMPS.get(tensorPosition, Tensor0);
   Tensor0.getShape(Shape3);
   this->Tensors[tensorPosition].getShape(Shape5);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[2];
   Tensor0.reshape(Shape4);
   ColumnNormMPO.set(2, Tensor0);
   ColumnNormMPO.set(3, ColumnNormTensors[tensorPosition+1]);
  }
// a)c) tensor-position L-1:
  else if (tensorPosition == L-1)
  {
   ColumnNormMPO.set(0, OneTensor);
   ColumnNormMPO.set(1, ColumnNormTensors[L-2]);
   ColumnNormMPS.get(L-1, Tensor0);
   Tensor0.getShape(Shape3);
   this->Tensors[L-1].getShape(Shape5);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[2];
   Tensor0.reshape(Shape4);
   ColumnNormMPO.set(2, Tensor0);
   ColumnNormMPO.set(3, OneTensor);
  }
 }
// b) right boundary:
 else if (Boundary == "right")
 {
// b)a) tensor-position 0:
  if (tensorPosition == 0)
  {
   ColumnNormMPS.get(L-1, Tensor0);
   Tensor0.getShape(Shape3);
   this->Tensors[Nrows*(Ncols-1)].getShape(Shape5);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape5[0]; Shape4[3] = Shape5[0];
   Tensor0.reshape(Shape4);
   ColumnNormMPO.set(0, Tensor0);
   ColumnNormMPO.set(1, OneTensor);
   ColumnNormMPO.set(2, OneTensor);
   ColumnNormMPO.set(3, ColumnNormTensors[1]);
  }
// b)b) 0 < tensor-position < L-1:
  else if ((tensorPosition > 0) && (tensorPosition < L-1))
  {
   ColumnNormMPS.get(L-1-tensorPosition, Tensor0);
   Tensor0.getShape(Shape3);
   this->Tensors[tensorPosition+Nrows*(Ncols-1)].getShape(Shape5);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[0]; Shape4[3] = Shape5[0];
   Tensor0.reshape(Shape4);
   ColumnNormMPO.set(0, Tensor0);
   ColumnNormMPO.set(1, ColumnNormTensors[tensorPosition-1]);
   ColumnNormMPO.set(2, OneTensor);
   ColumnNormMPO.set(3, ColumnNormTensors[tensorPosition+1]);
  }
// b)c) tensor-position L-1:
  else if (tensorPosition == L-1)
  {
   ColumnNormMPS.get(0, Tensor0);
   Tensor0.getShape(Shape3);
   this->Tensors[L-1+Nrows*(Ncols-1)].getShape(Shape5);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[0]; Shape4[3] = Shape5[0];
   Tensor0.reshape(Shape4);
   ColumnNormMPO.set(0, Tensor0);
   ColumnNormMPO.set(1, ColumnNormTensors[L-2]);
   ColumnNormMPO.set(2, OneTensor);
   ColumnNormMPO.set(3, OneTensor);
  }
 }
}

template<class T> void PEPS<T>::getColumnbMPOBoundary(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary, unsigned int tensorPosition,
                                                      const MPS<T>& ColumnbMPS, const vector< Tensor<T> >& ColumnbTensors,
                                                      MPO<T>& ColumnbMPO) const
{
 unsigned int L = this->N[0], Ncols = this->N[1], Nrows = this->N[0];
 vector<unsigned int> Shape3(3), Shape4(4), Shape5(5), Shape05(5), Shape15(5), Shape6(6);
 Tensor<T> OneTensor, Tensor0, Tensor1;
// OneTensor:
 Shape4[0] = 1; Shape4[1] = 1; Shape4[2] = 1; Shape4[3] = 1;
 OneTensor = Tensor<T>(Shape4);
 OneTensor.set(0, 1.0);
// ColumnbMPO:
 ColumnbMPO = MPO<T>("periodic", 4, max(PEPS0.D*PEPO0.getD(), this->D), ColumnbMPS.getD());
// a) left boundary:
 if (Boundary == "left")
 {
// a)a) tensor-position 0:
  if (tensorPosition == 0)
  {
   ColumnbMPO.set(0, OneTensor);
   ColumnbMPO.set(1, OneTensor);
   ColumnbMPS.get(0, Tensor0);
   Tensor0.getShape(Shape3);
   PEPS0.Tensors[0].getShape(Shape05);
   PEPO0.get(0, 0, Tensor1); Tensor1.getShape(Shape6);
   this->Tensors[0].getShape(Shape5);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[2]*Shape6[2]; Shape4[3] = Shape5[2];
   Tensor0.reshape(Shape4);
   ColumnbMPO.set(2, Tensor0);
   Tensor0 = ColumnbTensors[1];
   Tensor0.getShape(Shape15);
   Shape4[0] = Shape15[0]; Shape4[1] = 1; Shape4[2] = Shape05[3]*Shape6[3]; Shape4[3] = Shape5[3];
   Tensor0.reshape(Shape4);
   ColumnbMPO.set(3, Tensor0);
  }
// a)b) 0 < tensor-position < L-1:
  else if ((tensorPosition > 0) && (tensorPosition < L-1))
  {
   ColumnbMPO.set(0, OneTensor);
   Tensor0 = ColumnbTensors[tensorPosition-1];
   Tensor0.getShape(Shape15);
   PEPS0.Tensors[tensorPosition].getShape(Shape05);
   PEPO0.get(tensorPosition, 0, Tensor1); Tensor1.getShape(Shape6);
   this->Tensors[tensorPosition].getShape(Shape5);
   Shape4[0] = 1; Shape4[1] = Shape15[1]; Shape4[2] = Shape05[1]*Shape6[1]; Shape4[3] = Shape5[1];
   Tensor0.reshape(Shape4);
   ColumnbMPO.set(1, Tensor0);
   ColumnbMPS.get(tensorPosition, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[2]*Shape6[2]; Shape4[3] = Shape5[2];
   Tensor0.reshape(Shape4);
   ColumnbMPO.set(2, Tensor0);
   Tensor0 = ColumnbTensors[tensorPosition+1];
   Tensor0.getShape(Shape15);
   Shape4[0] = Shape15[0]; Shape4[1] = 1; Shape4[2] = Shape05[3]*Shape6[3]; Shape4[3] = Shape5[3];
   Tensor0.reshape(Shape4);
   ColumnbMPO.set(3, Tensor0);
  }
// a)c) tensor-position L-1:
  else if (tensorPosition == L-1)
  {
   ColumnbMPO.set(0, OneTensor);
   Tensor0 = ColumnbTensors[L-2];
   Tensor0.getShape(Shape15);
   PEPS0.Tensors[L-1].getShape(Shape05);
   PEPO0.get(L-1, 0, Tensor1); Tensor1.getShape(Shape6);
   this->Tensors[L-1].getShape(Shape5);
   Shape4[0] = 1; Shape4[1] = Shape15[1]; Shape4[2] = Shape05[1]*Shape6[1]; Shape4[3] = Shape5[1];
   Tensor0.reshape(Shape4);
   ColumnbMPO.set(1, Tensor0);
   ColumnbMPS.get(L-1, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape05[2]*Shape6[2]; Shape4[3] = Shape5[2];
   Tensor0.reshape(Shape4);
   ColumnbMPO.set(2, Tensor0);
   ColumnbMPO.set(3, OneTensor);
  }
 }
// b) right boundary:
 else if (Boundary == "right")
 {
// b)a) tensor-position 0:
  if (tensorPosition == 0)
  {
   ColumnbMPS.get(L-1, Tensor0);
   Tensor0.getShape(Shape3);
   PEPS0.Tensors[Nrows*(Ncols-1)].getShape(Shape05);
   PEPO0.get(0, Ncols-1, Tensor1); Tensor1.getShape(Shape6);
   this->Tensors[Nrows*(Ncols-1)].getShape(Shape5);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape05[0]*Shape6[0]; Shape4[3] = Shape5[0];
   Tensor0.reshape(Shape4);
   ColumnbMPO.set(0, Tensor0);
   ColumnbMPO.set(1, OneTensor);
   ColumnbMPO.set(2, OneTensor);
   Tensor0 = ColumnbTensors[1];
   Tensor0.getShape(Shape15);
   Shape4[0] = 1; Shape4[1] = Shape15[1]; Shape4[2] = Shape05[3]*Shape6[3]; Shape4[3] = Shape5[3];
   Tensor0.reshape(Shape4);
   ColumnbMPO.set(3, Tensor0);
  }
// b)b) 0 < tensor-position < L-1:
  else if ((tensorPosition > 0) && (tensorPosition < L-1))
  {
   ColumnbMPS.get(L-1-tensorPosition, Tensor0);
   Tensor0.getShape(Shape3);
   PEPS0.Tensors[tensorPosition+Nrows*(Ncols-1)].getShape(Shape05);
   PEPO0.get(tensorPosition, Ncols-1, Tensor1); Tensor1.getShape(Shape6);
   this->Tensors[tensorPosition+Nrows*(Ncols-1)].getShape(Shape5);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[0]*Shape6[0]; Shape4[3] = Shape5[0];
   Tensor0.reshape(Shape4);
   ColumnbMPO.set(0, Tensor0);
   Tensor0 = ColumnbTensors[tensorPosition-1];
   Tensor0.getShape(Shape15);
   Shape4[0] = Shape15[0]; Shape4[1] = 1; Shape4[2] = Shape05[1]*Shape6[1]; Shape4[3] = Shape5[1];
   Tensor0.reshape(Shape4);
   ColumnbMPO.set(1, Tensor0);
   ColumnbMPO.set(2, OneTensor);
   Tensor0 = ColumnbTensors[tensorPosition+1];
   Tensor0.getShape(Shape15);
   Shape4[0] = 1; Shape4[1] = Shape15[1]; Shape4[2] = Shape05[3]*Shape6[3]; Shape4[3] = Shape5[3];
   Tensor0.reshape(Shape4);
   ColumnbMPO.set(3, Tensor0);
  }
// b)c) tensor-position L-1:
  else if (tensorPosition == L-1)
  {
   ColumnbMPS.get(0, Tensor0);
   Tensor0.getShape(Shape3);
   PEPS0.Tensors[L-1+Nrows*(Ncols-1)].getShape(Shape05);
   PEPO0.get(L-1, Ncols-1, Tensor1); Tensor1.getShape(Shape6);
   this->Tensors[L-1+Nrows*(Ncols-1)].getShape(Shape5);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[0]*Shape6[0]; Shape4[3] = Shape5[0];
   Tensor0.reshape(Shape4);
   ColumnbMPO.set(0, Tensor0);
   Tensor0 = ColumnbTensors[L-2];
   Tensor0.getShape(Shape15);
   Shape4[0] = Shape15[0]; Shape4[1] = 1; Shape4[2] = Shape05[1]*Shape6[1]; Shape4[3] = Shape5[1];
   Tensor0.reshape(Shape4);
   ColumnbMPO.set(1, Tensor0);
   ColumnbMPO.set(2, OneTensor);
   ColumnbMPO.set(3, OneTensor);
  }
 }
}

template<class T> void PEPS<T>::getColumnNormMPOBulk(unsigned int position, unsigned int tensorPosition,
                                                     const MPS<T>& ColumnNormMPSL, const MPS<T>& ColumnNormMPSR, const vector< Tensor<T> >& ColumnNormTensors,
                                                     MPO<T>& ColumnNormMPO) const
{
 unsigned int L = this->N[0], Nrows = this->N[0];
 vector<unsigned int> Shape3(3), Shape4(4), Shape5(5);
 Tensor<T> OneTensor, Tensor0;
// OneTensor:
 Shape4[0] = 1; Shape4[1] = 1; Shape4[2] = 1; Shape4[3] = 1;
 OneTensor = Tensor<T>(Shape4);
 OneTensor.set(0, 1.0);
// ColumnNormMPO:
 ColumnNormMPO = MPO<T>("periodic", 4, this->D, max(ColumnNormMPSL.getD(), ColumnNormMPSR.getD()));
// a) tensor-position 0:
 if (tensorPosition == 0)
 {
  ColumnNormMPSL.get(L-1, Tensor0);
  Tensor0.getShape(Shape3);
  this->Tensors[Nrows*position].getShape(Shape5);
  Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape5[0]; Shape4[3] = Shape5[0];
  Tensor0.reshape(Shape4);
  ColumnNormMPO.set(0, Tensor0);
  ColumnNormMPO.set(1, OneTensor);
  ColumnNormMPSR.get(0, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[2];
  Tensor0.reshape(Shape4);
  ColumnNormMPO.set(2, Tensor0);
  ColumnNormMPO.set(3, ColumnNormTensors[1]);
 }
// b) 0 < tensor-position < L-1:
 else if ((tensorPosition > 0) && (tensorPosition < L-1))
 {
  ColumnNormMPSL.get(L-1-tensorPosition, Tensor0);
  Tensor0.getShape(Shape3);
  this->Tensors[tensorPosition+Nrows*position].getShape(Shape5);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[0]; Shape4[3] = Shape5[0];
  Tensor0.reshape(Shape4);
  ColumnNormMPO.set(0, Tensor0);
  ColumnNormMPO.set(1, ColumnNormTensors[tensorPosition-1]);
  ColumnNormMPSR.get(tensorPosition, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[2];
  Tensor0.reshape(Shape4);
  ColumnNormMPO.set(2, Tensor0);
  ColumnNormMPO.set(3, ColumnNormTensors[tensorPosition+1]);
 }
// c) tensor-position L-1:
 else if (tensorPosition == L-1)
 {
  ColumnNormMPSL.get(0, Tensor0);
  Tensor0.getShape(Shape3);
  this->Tensors[L-1+Nrows*position].getShape(Shape5);
  Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape5[0]; Shape4[3] = Shape5[0];
  Tensor0.reshape(Shape4);
  ColumnNormMPO.set(0, Tensor0);
  ColumnNormMPO.set(1, ColumnNormTensors[L-2]);
  ColumnNormMPSR.get(L-1, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape5[2]; Shape4[3] = Shape5[2];
  Tensor0.reshape(Shape4);
  ColumnNormMPO.set(2, Tensor0);
  ColumnNormMPO.set(3, OneTensor);
 }
}

template<class T> void PEPS<T>::getColumnbMPOBulk(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, unsigned int position, unsigned int tensorPosition,
                                                  const MPS<T>& ColumnbMPSL, const MPS<T>& ColumnbMPSR, const vector< Tensor<T> >& ColumnbTensors,
                                                  MPO<T>& ColumnbMPO) const
{
 unsigned int L = this->N[0], Nrows = this->N[0];
 vector<unsigned int> Shape3(3), Shape4(4), Shape5(5), Shape05(5), Shape15(5), Shape6(6);
 Tensor<T> OneTensor, Tensor0, Tensor1;
// OneTensor:
 Shape4[0] = 1; Shape4[1] = 1; Shape4[2] = 1; Shape4[3] = 1;
 OneTensor = Tensor<T>(Shape4);
 OneTensor.set(0, 1.0);
// ColumnbMPO:
 ColumnbMPO = MPO<T>("periodic", 4, max(PEPS0.D*PEPO0.getD(), this->D), max(ColumnbMPSL.getD(), ColumnbMPSR.getD()));
// a) tensor-position 0:
 if (tensorPosition == 0)
 {
  ColumnbMPSL.get(L-1, Tensor0);
  Tensor0.getShape(Shape3);
  PEPS0.Tensors[Nrows*position].getShape(Shape05);
  PEPO0.get(0, position, Tensor1); Tensor1.getShape(Shape6);
  this->Tensors[Nrows*position].getShape(Shape5);
  Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape05[0]*Shape6[0]; Shape4[3] = Shape5[0];
  Tensor0.reshape(Shape4);
  ColumnbMPO.set(0, Tensor0);
  ColumnbMPO.set(1, OneTensor);
  ColumnbMPSR.get(0, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[2]*Shape6[2]; Shape4[3] = Shape5[2];
  Tensor0.reshape(Shape4);
  ColumnbMPO.set(2, Tensor0);
  Tensor0 = ColumnbTensors[1];
  Tensor0.getShape(Shape15);
  Shape4[0] = Shape15[0]; Shape4[1] = Shape15[1]; Shape4[2] = Shape05[3]*Shape6[3]; Shape4[3] = Shape5[3];
  Tensor0.reshape(Shape4);
  ColumnbMPO.set(3, Tensor0);
 }
// b) 0 < tensor-position < L-1:
 else if ((tensorPosition > 0) && (tensorPosition < L-1))
 {
  ColumnbMPSL.get(L-1-tensorPosition, Tensor0);
  Tensor0.getShape(Shape3);
  PEPS0.Tensors[tensorPosition+Nrows*position].getShape(Shape05);
  PEPO0.get(tensorPosition, position, Tensor1); Tensor1.getShape(Shape6);
  this->Tensors[tensorPosition+Nrows*position].getShape(Shape5);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[0]*Shape6[0]; Shape4[3] = Shape5[0];
  Tensor0.reshape(Shape4);
  ColumnbMPO.set(0, Tensor0);
  Tensor0 = ColumnbTensors[tensorPosition-1];
  Tensor0.getShape(Shape15);
  Shape4[0] = Shape15[0]; Shape4[1] = Shape15[1]; Shape4[2] = Shape05[1]*Shape6[1]; Shape4[3] = Shape5[1];
  Tensor0.reshape(Shape4);  
  ColumnbMPO.set(1, Tensor0);
  ColumnbMPSR.get(tensorPosition, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[2]*Shape6[2]; Shape4[3] = Shape5[2];
  Tensor0.reshape(Shape4);
  ColumnbMPO.set(2, Tensor0);
  Tensor0 = ColumnbTensors[tensorPosition+1];
  Tensor0.getShape(Shape15);
  Shape4[0] = Shape15[0]; Shape4[1] = Shape15[1]; Shape4[2] = Shape05[3]*Shape6[3]; Shape4[3] = Shape5[3];
  Tensor0.reshape(Shape4);
  ColumnbMPO.set(3, Tensor0);
 }
// c) tensor-position L-1:
 else if (tensorPosition == L-1)
 {
  ColumnbMPSL.get(0, Tensor0);
  Tensor0.getShape(Shape3);
  PEPS0.Tensors[L-1+Nrows*position].getShape(Shape05);
  PEPO0.get(L-1, position, Tensor1); Tensor1.getShape(Shape6);
  this->Tensors[L-1+Nrows*position].getShape(Shape5);
  Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = Shape05[0]*Shape6[0]; Shape4[3] = Shape5[0];
  Tensor0.reshape(Shape4);
  ColumnbMPO.set(0, Tensor0);
  Tensor0 = ColumnbTensors[L-2];
  Tensor0.getShape(Shape15);
  Shape4[0] = Shape15[0]; Shape4[1] = Shape15[1]; Shape4[2] = Shape05[1]*Shape6[1]; Shape4[3] = Shape5[1];
  Tensor0.reshape(Shape4);  
  ColumnbMPO.set(1, Tensor0);
  ColumnbMPSR.get(L-1, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = Shape05[2]*Shape6[2]; Shape4[3] = Shape5[2];
  Tensor0.reshape(Shape4);
  ColumnbMPO.set(2, Tensor0);
  ColumnbMPO.set(3, OneTensor);
 }
}

template<class T> void PEPS<T>::updateTensor(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, unsigned int positionRow, unsigned int positionCol,
                                             const MPO<T>& NormMPO, const MPO<T>& bMPO,
                                             const string& UpdateTensor, const string& UpdateMode, double cutoff,
                                             double& epsAchieved)
{
 if (UpdateTensor == "reduced")
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPS<T>::" <<
          "updateTensor(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, unsigned int positionRow, unsigned int positionCol, " <<
                       "const MPO<T>& NormMPO, const MPO<T>& bMPO, const string& UpdateTensor, const string& UpdateMode, double cutoff, double& epsAchieved): " <<
          "(UpdateTensor == reduced) is not implemented yet." << endl;
  exit(1);
 }
 unsigned int dim, dimPos, index0, index1, Nrows = this->N[0];
 double error0, error1;
 T bValue, NValue;
 vector<unsigned int> Index5(5), Index8(8), Indices0(1), Indices02(2), Indices03(3), Indices06(6), Indices1(1), Indices12(2), Indices13(3), Indices16(6);
 vector<unsigned int> Shape(5), Shape0(5), Shape4(4), Shape5(5), ShapePEPO(6);
 vector<double> W;
 vector<T> bVector, Vector0;
 Tensor<T> Tensor0, Tensor1, Tensor2;
 Matrix<T> NMatrix, NMatrixAdj, NMatrixH, NMatrixInv, Vr;
// 1. dimension:
 this->Tensors[positionRow+Nrows*positionCol].getShape(Shape);
 dim = 1;
 for (int i = 0; i < 5; i++)
  dim *= Shape[i];
// 2. NMatrix:
 NormMPO.get(0, Tensor0);
 NormMPO.get(1, Tensor1);
 Indices0[0] = 1; Indices1[0] = 0;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 NormMPO.get(2, Tensor1);
 NormMPO.get(3, Tensor2);
 Tensor1.contract(Indices0, Tensor2, Indices1);
 Indices02[0] = 0; Indices02[1] = 3; Indices12[0] = 3; Indices12[1] = 0;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 NMatrix = Matrix<T>(dim, dim);
 NMatrix.fillZeroes();
 for (int k = 0; k < Shape[4]; k++){
  for (int j3 = 0; j3 < Shape[3]; j3++){
   Index8[6] = j3;
  for (int j2 = 0; j2 < Shape[2]; j2++){
   Index8[4] = j2;
  for (int j1 = 0; j1 < Shape[1]; j1++){
   Index8[2] = j1;
  for (int j0 = 0; j0 < Shape[0]; j0++){
   Index8[0] = j0;
   index1 = j0 + j1*Shape[0] + j2*Shape[0]*Shape[1] + j3*Shape[0]*Shape[1]*Shape[2] + k*Shape[0]*Shape[1]*Shape[2]*Shape[3];
   for (int i3 = 0; i3 < Shape[3]; i3++){
    Index8[7] = i3;
   for (int i2 = 0; i2 < Shape[2]; i2++){
    Index8[5] = i2;
   for (int i1 = 0; i1 < Shape[1]; i1++){
    Index8[3] = i1;
   for (int i0 = 0; i0 < Shape[0]; i0++){
    Index8[1] = i0;
    index0 = i0 + i1*Shape[0] + i2*Shape[0]*Shape[1] + i3*Shape[0]*Shape[1]*Shape[2] + k*Shape[0]*Shape[1]*Shape[2]*Shape[3];
    NMatrix(index0, index1) = Tensor0.get(Index8);
   }}}}
  }}}}
 }
// 3. bVector:
 bMPO.get(0, Tensor0);
 Tensor0.getShape(Shape4);
 PEPS0.Tensors[positionRow+Nrows*positionCol].getShape(Shape0);
 PEPO0.get(positionRow, positionCol, Tensor1); Tensor1.getShape(ShapePEPO);
 Shape5[0] = Shape4[0]; Shape5[1] = Shape4[1]; Shape5[2] = Shape0[0]; Shape5[3] = ShapePEPO[0]; Shape5[4] = Shape[0];
 Tensor0.reshape(Shape5);
 bMPO.get(1, Tensor1);
 Tensor1.getShape(Shape4);
 Shape5[0] = Shape4[0]; Shape5[1] = Shape4[1]; Shape5[2] = Shape0[1]; Shape5[3] = ShapePEPO[1]; Shape5[4] = Shape[1];
 Tensor1.reshape(Shape5);
 Indices0[0] = 1; Indices1[0] = 0;
 Tensor0.contract(Indices0, Tensor1, Indices1);
 Tensor1 = PEPS0.Tensors[positionRow+Nrows*positionCol];
 Indices02[0] = 1; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 1;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 PEPO0.get(positionRow, positionCol, Tensor1);
 Indices03[0] = 1; Indices03[1] = 4; Indices03[2] = 8; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 4;
 Tensor0.contract(Indices03, Tensor1, Indices13);
 bMPO.get(2, Tensor1);
 Tensor1.getShape(Shape4);
 Shape5[0] = Shape4[0]; Shape5[1] = Shape4[1]; Shape5[2] = Shape0[2]; Shape5[3] = ShapePEPO[2]; Shape5[4] = Shape[2];
 Tensor1.reshape(Shape5);
 bMPO.get(3, Tensor2);
 Tensor2.getShape(Shape4);
 Shape5[0] = Shape4[0]; Shape5[1] = Shape4[1]; Shape5[2] = Shape0[3]; Shape5[3] = ShapePEPO[3]; Shape5[4] = Shape[3];
 Tensor2.reshape(Shape5);
 Tensor1.contract(Indices0, Tensor2, Indices1);
 Indices06[0] = 0; Indices06[1] = 2; Indices06[2] = 4; Indices06[3] = 5; Indices06[4] = 6; Indices06[5] = 7;
 Indices16[0] = 4; Indices16[1] = 0; Indices16[2] = 1; Indices16[3] = 5; Indices16[4] = 2; Indices16[5] = 6;
 Tensor0.contract(Indices06, Tensor1, Indices16);
 bVector = vector<T>(dim, 0.0);
 for (int k = 0; k < Shape[4]; k++){
  Index5[2] = k;
  for (int i3 = 0; i3 < Shape[3]; i3++){
   Index5[4] = i3;
  for (int i2 = 0; i2 < Shape[2]; i2++){
   Index5[3] = i2;
  for (int i1 = 0; i1 < Shape[1]; i1++){
   Index5[1] = i1;
  for (int i0 = 0; i0 < Shape[0]; i0++){
   Index5[0] = i0;
   index0 = i0 + i1*Shape[0] + i2*Shape[0]*Shape[1] + i3*Shape[0]*Shape[1]*Shape[2] + k*Shape[0]*Shape[1]*Shape[2]*Shape[3];
   bVector[index0] = Tensor0.get(Index5);
  }}}}
 }
// 4. initial error0:
 Vector0 = vector<T>(dim);
 for (int i = 0; i < dim; i++)
  Vector0[i] = this->Tensors[positionRow+Nrows*positionCol].get(i);
 NValue = expectationValueMatrix(Vector0, NMatrix);
 bValue = scalarProductVector(Vector0, bVector);
 error0 = MathAuxiliary::convertToDouble(NValue)-2.0*MathAuxiliary::convertToDouble(bValue);
// 5. NMatrixInv:
// 5.a) via pseudoinverse:
 if (UpdateMode == "Pseudoinverse")
 {
  NMatrix.pseudoinvert(cutoff, NMatrixInv);
 }
// 5.b) via closest positive norm-matrix:
 else if (UpdateMode == "PositiveNormMatrix")
 {
  NMatrix.adjoint(NMatrixAdj);
  NMatrixH = NMatrix;
  NMatrixH.add(NMatrixAdj);
  NMatrixH.multiply(0.5);
  NMatrixH.setType("hermitian");
  W = vector<double>(dim); Vr = Matrix<T>(dim, dim);
  NMatrixH.eigenDecompose(W, Vr);
  NMatrixInv = Matrix<T>(dim, dim);
  NMatrixInv.fillZeroes();
  dimPos = 0;
  for (int k = 0; k < dim; k++){
   if (W[k]/abs(W[dim-1]) > cutoff){
    dimPos++;
    for (int j = 0; j < dim; j++){
     for (int i = 0; i < dim; i++){
      NMatrixInv(i, j) += Vr(i, k)*MathAuxiliary::complexConjugate(Vr(j, k))/W[k];
     }
    }
   }
  }
  if (dimPos == 0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void PEPS<T>::" <<
           "updateTensor(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, unsigned int positionRow, unsigned int positionCol, " <<
                        "const MPO<T>& NormMPO, const MPO<T>& bMPO, const string& UpdateTensor, const string& UpdateMode, double cutoff, double& epsAchieved): " <<
           "All eigenvalues of the norm-matrix were negative. Increase the bond dimension of the environment!" << endl;
   exit(1);
  }
 }
 NMatrixInv.multiply(bVector, Vector0);
// 6. write solution into this PEPS:
 for (int i = 0; i < dim; i++)
  this->Tensors[positionRow+Nrows*positionCol].set(i, Vector0[i]);
// 7. compute final error:
 NValue = expectationValueMatrix(Vector0, NMatrix);
 bValue = scalarProductVector(Vector0, bVector);
 error1 = MathAuxiliary::convertToDouble(NValue)-2.0*MathAuxiliary::convertToDouble(bValue);
 epsAchieved = abs((error1-error0)/error1);
}

template<class T> void PEPS<T>::updateBoundaryColumn(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, const string& Boundary,
                                                     double eps,
                                                     const vector< MPS<T> >& ColumnNormMPSs, const vector< MPS<T> >& ColumnbMPSs,
                                                     const string& UpdateTensor, const string& UpdateMode, double cutoff,
                                                     double& epsAchieved)
{
 unsigned int col, Ncols = this->N[1], Nrows = this->N[0], row;
 MPO<T> bMPO, NormMPO;
 vector< Tensor<T> > bTensors, NormTensors;
// a) left boundary-column:
 if (Boundary == "left")
 {
  col = 0;
// a)1. get initial norm-tensors and b-tensors:
  this->getInitialColumnNormTensorsBoundary(Boundary, ColumnNormMPSs[1], NormTensors);
  this->getInitialColumnbTensorsBoundary(PEPO0, PEPS0, Boundary, ColumnbMPSs[1], bTensors);
// a)2. update boundary-column:
  for (row = 0; row < Nrows; row++)
  {
// a)2.1. get norm-MPO and b-MPO:
   this->getColumnNormMPOBoundary(Boundary, row, ColumnNormMPSs[1], NormTensors, NormMPO);
   this->getColumnbMPOBoundary(PEPO0, PEPS0, Boundary, row, ColumnbMPSs[1], bTensors, bMPO);
// a)2.2. update tensor:
   this->updateTensor(PEPO0, PEPS0, row, col, NormMPO, bMPO, UpdateTensor, UpdateMode, cutoff, epsAchieved);
   if (epsAchieved < eps)
    return;
// a)2.3. update norm-tensors and b-tensors:
   this->updateColumnNormTensorsBoundary(Boundary, row, ColumnNormMPSs[1], NormTensors);
   this->updateColumnbTensorsBoundary(PEPO0, PEPS0, Boundary, row, ColumnbMPSs[1], bTensors);
  }
 }
// b) right boundary-column:
 else if (Boundary == "right")
 {
  col = Ncols-1;
// b)1. get initial norm-tensors and b-tensors:
  this->getInitialColumnNormTensorsBoundary(Boundary, ColumnNormMPSs[Ncols-2], NormTensors);
  this->getInitialColumnbTensorsBoundary(PEPO0, PEPS0, Boundary, ColumnbMPSs[Ncols-2], bTensors);
// b)2. update boundary-column:
  for (row = 0; row < Nrows; row++)
  {
// b)2.1. get norm-MPO and b-MPO:
   this->getColumnNormMPOBoundary(Boundary, row, ColumnNormMPSs[Ncols-2], NormTensors, NormMPO);
   this->getColumnbMPOBoundary(PEPO0, PEPS0, Boundary, row, ColumnbMPSs[Ncols-2], bTensors, bMPO);
// b)2.2. update tensor:
   this->updateTensor(PEPO0, PEPS0, row, col, NormMPO, bMPO, UpdateTensor, UpdateMode, cutoff, epsAchieved);
   if (epsAchieved < eps)
    return;
// b)2.3. update norm-tensors and b-tensors:
   this->updateColumnNormTensorsBoundary(Boundary, row, ColumnNormMPSs[Ncols-2], NormTensors);
   this->updateColumnbTensorsBoundary(PEPO0, PEPS0, Boundary, row, ColumnbMPSs[Ncols-2], bTensors);
  }
 }
}

template<class T> void PEPS<T>::updateBulkColumn(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0, unsigned int position,
                                                 double eps,
                                                 const vector< MPS<T> >& ColumnNormMPSs, const vector< MPS<T> >& ColumnbMPSs,
                                                 const string& UpdateTensor, const string& UpdateMode, double cutoff,
                                                 double& epsAchieved)
{
 unsigned int col = position, Nrows = this->N[0], row;
 MPO<T> bMPO, NormMPO;
 vector< Tensor<T> > bTensors, NormTensors;
// 1. get initial norm-tensors and b-tensors:
 this->getInitialColumnNormTensorsBulk(col, ColumnNormMPSs[col-1], ColumnNormMPSs[col+1], NormTensors);
 this->getInitialColumnbTensorsBulk(PEPO0, PEPS0, col, ColumnbMPSs[col-1], ColumnbMPSs[col+1], bTensors);
// 2. update bulk-column:
 for (row = 0; row < Nrows; row++)
 {
// 2.1. get norm-MPO and b-MPO:
  this->getColumnNormMPOBulk(col, row, ColumnNormMPSs[col-1], ColumnNormMPSs[col+1], NormTensors, NormMPO);
  this->getColumnbMPOBulk(PEPO0, PEPS0, col, row, ColumnbMPSs[col-1], ColumnbMPSs[col+1], bTensors, bMPO);
// 2.2. update tensor:
  this->updateTensor(PEPO0, PEPS0, row, col, NormMPO, bMPO, UpdateTensor, UpdateMode, cutoff, epsAchieved);
  if (epsAchieved < eps)
   return;
// 2.3. update norm-tensors and b-tensors:
  this->updateColumnNormTensorsBulk(col, row, ColumnNormMPSs[col-1], ColumnNormMPSs[col+1], NormTensors);
  this->updateColumnbTensorsBulk(PEPO0, PEPS0, col, row, ColumnbMPSs[col-1], ColumnbMPSs[col+1], bTensors);
 }
}
