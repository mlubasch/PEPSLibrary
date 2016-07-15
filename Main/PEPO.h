/// Template class PEPO implements Projected Entangled-Pair Operators.
/** The template class PEPO implements Projected Entangled-Pair Operators in 2D.
    PEPOs consist of Tensors in a matrix.
    \param BC string, the boundary conditions of the PEPO, is "open" or "periodic"
    \param N vector<unsigned int>(2), the number of Tensors of the PEPO, fulfills N.size()==2
    \param d Matrix<unsigned int>(N[0], N[1]), the physical dimensions of the PEPO, fulfills
             d.getDim0()==N[0] and d.getDim1()==N[1]
    \param D unsigned int, the maximal virtual bond dimension of the PEPO
    \param Tensors Tensor<T>*, the tensors of the PEPO, stored in Fortran's column-major order
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

template<class T> class PEPO
{
 public:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0.
    \sa PEPO(const string& BC0, unsigned int Nrows, unsigned int Ncols, const Matrix<unsigned int>& d0,
             unsigned int D0)
    \sa PEPO(const string& BC0, unsigned int Nrows, unsigned int Ncols, unsigned int d0,
             unsigned int D0) */
  PEPO();

/// Constructor for PEPO with specific BC, N, d and D.
/** This constructor initializes a PEPO of a specific form with Nrows rows and Ncols columns.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param Nrows input: unsigned int, the number of rows, must be > 1
    \param Ncols input: unsigned int, the number of columns, must be > 1
    \param d0 input: const Matrix<unsigned int>&, the physical dimensions, must fulfill
                     d0.getDim0()==Nrows and d0.getDim1()==Ncols, and all entries must be > 1
    \param D0 input: unsigned int, the maximal virtual bond dimension, must be > 0
    \sa PEPO()
    \sa PEPO(const string& BC0, unsigned int Nrows, unsigned int Ncols, unsigned int d0,
             unsigned int D0) */
  PEPO(const string& BC0, unsigned int Nrows, unsigned int Ncols, const Matrix<unsigned int>& d0,
       unsigned int D0);

/// Constructor for PEPO with specific BC, N, d and D.
/** This constructor initializes a PEPO of a specific form with Nrows rows, Ncols columns and all
    physical dimensions equal to d0.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param Nrows input: unsigned int, the number of rows, must be > 1
    \param Ncols input: unsigned int, the number of columns, must be > 1
    \param d0 input: unsigned int, the physical dimension, must be > 1
    \param D0 input: unsigned int, the maximal virtual bond dimension, must be > 0
    \sa PEPO()
    \sa PEPO(const string& BC0, unsigned int Nrows, unsigned int Ncols, const Matrix<unsigned int>& d0,
             unsigned int D0) */
  PEPO(const string& BC0, unsigned int Nrows, unsigned int Ncols, unsigned int d0, unsigned int D0);

/// Standard copy constructor.
/** The standard copy constructor copies the input PEPO into this.
    \param PEPO0 input: const PEPO<T>&, to be copied into this
    \sa PEPO<T>& operator=(const PEPO<T>& PEPO0) */
  PEPO(const PEPO<T>& PEPO0);

/// Standard destructor.
/** The standard destructor deletes the elements of the PEPO. */
  ~PEPO();

/// Assigns PEPO to this.
/** The operator= allows to assign a PEPO to this. Hereby this is destroyed and newly constructed to be
    a copy of the right-hand side PEPO.
    \param PEPO0 input: const PEPO<T>&, to be copied into this
    \return PEPO<T>&, a reference to the new this
    \sa PEPO(const PEPO<T>& PEPO0) */
  PEPO<T>& operator=(const PEPO<T>& PEPO0);

/// Returns boundary conditions.
/** The returned boundary conditions are either "open" or "periodic".
    \param BC0 output: string&, the boundary conditions of this PEPO */
  void getBC(string& BC0) const { BC0 = this->BC; }

/// Returns number of rows Nrows.
/** This function returns the number of rows of this PEPO.
    \return unsigned int, the number of rows of this PEPO */
  unsigned int getNrows() const { return this->N[0]; }

/// Returns number of columns Ncols.
/** This function returns the number of columns of this PEPO.
    \return unsigned int, the number of columns of this PEPO */
  unsigned int getNcols() const { return this->N[1]; }

/// Returns physical dimensions d.
/** This function returns the physical dimensions of this PEPO.
    \param d0 output: Matrix<unsigned int>&, the physical dimensions of this PEPO */
  void getd(Matrix<unsigned int>& d0) const { d0 = this->d; }

/// Returns physical dimension d.
/** This function returns the physical dimension of this PEPO at site (0, 0). It is useful, if this
    PEPO has all physical dimensions equal to each other.
    \return unsigned int, the physical dimension of this PEPO at site (0, 0) */
  unsigned int getd() const { return this->d(0, 0); }

/// Returns maximal virtual bond dimension D.
/** This function returns the maximal virtual bond dimension of this PEPO.
    \return unsigned int, the maximal virtual bond dimension of this PEPO */
  unsigned int getD() const { return this->D; }

/// Sets Tensor at position.
/** This function sets Tensor0 at a given row number positionRow and column number positionCol in this
    PEPO.
    \param positionRow input: unsigned int, the row position
    \param positionCol input: unsigned int, the column position
    \param Tensor0 input: const Tensor<T>&, to be written at position, must have the correct Shape */
  void set(unsigned int positionRow, unsigned int positionCol, const Tensor<T>& Tensor0);

/// Returns Tensor at position.
/** This function returns as Tensor0 a copy of the tensor at a given row number positionRow and column
    number positionCol in this PEPO.
    \param positionRow input: unsigned int, the row position
    \param positionCol input: unsigned int, the column position
    \param Tensor0 output: Tensor<T>&, a copy of the tensor at position */
  void get(unsigned int positionRow, unsigned int positionCol, Tensor<T>& Tensor0) const;

/// Writes this PEPO to binary file.
/** Given a file name FileName, a new binary file is constructed into which this PEPO is written.
    A PEPO is represented in a binary file by:
    {BCSize, BC[0], ..., BC[BCSize-1], N[0], N[1], d[0], ..., d[N[0]*N[1]-1], D, Tensors[0], ..., Tensors[N[0]*N[1]-1]}
    where each Tensor is represented by:
    {rank, Shape[0], ..., Shape[rank-1], size, Elements[0], ..., Elements[size-1]}.
    Note that, as always, d is stored in Fortran's column-major order.
    \param FileName input: const string&, the name for the new binary file to which this PEPO is written
    \sa void read(const string& FileName) */
  void write(const string& FileName) const;

/// Reads PEPO from binary file.
/** Given a binary file called FileName, this PEPO is replaced by the PEPO in FileName.
    A PEPO is represented in a binary file by:
    {BCSize, BC[0], ..., BC[BCSize-1], N[0], N[1], d[0], ..., d[N[0]*N[1]-1], D, Tensors[0], ..., Tensors[N[0]*N[1]-1]}
    where each Tensor is represented by:
    {rank, Shape[0], ..., Shape[rank-1], size, Elements[0], ..., Elements[size-1]}.
    Note that, as always, d is stored in Fortran's column-major order.
    \param FileName input: const string&, the binary file from which this PEPO is read
    \sa void write(const string& FileName) const */
  void read(const string& FileName);

/// Returns this PEPO as matrix.
/** This PEPO is returned as a matrix in the standard basis.
    \param Matrix0 output: Matrix<T>&, this PEPO as a matrix,
                           must fulfill ((Matrix0.getDim0()==d^(N[0]*N[1])) && (Matrix0.getDim1()==d^(N[0]*N[1]))) */
  void getMatrix(Matrix<T>& Matrix0) const;

/// This PEPO is filled with random entries.
/** This PEPO is filled with uniformly distributed random entries. LAPACK's XLARNV is used.
    \param Seed input: Matrix<unsigned int>&, is not changed, contains the seed values for each tensor,
                       must fulfill ((Seed.getDim0()==this->N[0]) && (Seed.getDim1()==this->N[1]))
    \sa void PEPO<T>::fillZeroes()
    \sa void PEPO<T>::setOne() */
  void fillRandomly(Matrix<unsigned int>& Seed);

/// Fills this PEPO with zeroes.
/** This PEPO is filled with zeroes.
    \sa void PEPO<T>::fillRandomly(const Matrix<unsigned int>& Seed)
    \sa void PEPO<T>::setOne() */
  void fillZeroes();

/// Sets this PEPO as unity.
/** This PEPO is set as unity 1.
    \sa void PEPO<T>::fillRandomly(const Matrix<unsigned int>& Seed)
    \sa void PEPO<T>::fillZeroes() */
  void setOne();

/// This PEPO is initialized as Concatenated PEPO.
/** This PEPO is initialized as a Concatenated PEPO. Each tensor is replaced by a concatenation of tensors,
    consisting of the initial tensor and M^{2}-1 auxiliary tensors having physical dimension 1. Each tensor
    is replaced by a tensor network resembling a PEPO.
    M denotes the concatenation level, such that M==1 denotes this PEPO with no auxiliary tensors, M==2
    denotes this PEPO with 3 auxiliary tensors per physical tensor, M==3 denotes this PEPO with 8 auxiliary
    tensors per physical tensor, and so on. The auxiliary tensors are chosen as Deltas to guarantee that
    the resulting Concatenated PEPO is equivalent to this PEPO for all concatenation levels.
    \param M input: unsigned int, the concatenation level, must be > 0 */
  void setConcatenated(unsigned int M);

/// Multiplies this PEPO with element.
/** This PEPO is multiplied with element, by multiplying the first tensor with element.
    \param element input: T, the scalar with which this PEPO is multiplied */
  void multiply(T element);

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
    Either the pseudoinverse or the closest positive N-matrix is used in the local update algorithm.
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
                             "PositiveNMatrix" for the closest positive N-matrix
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

 protected:

/// Boundary conditions BC.
/** The boundary conditions are either "open" or "periodic". */
  string BC;

/// Number of tensors N.
/** N[0]=Nrows is the number of rows and N[1]=Ncols is the number of columns of this PEPO. */
  vector<unsigned int> N;

/// Physical dimensions d.
/** The Matrix d has dimensions d.getDim0()==N[0] and d.getDim1()==N[1]. */
  Matrix<unsigned int> d;

/// Maximal virtual bond dimension D.
  unsigned int D;

/// Tensors.
/** The Tensors of this PEPO are stored in Fortran's column-major order, such that the Tensor at matrix
    position (row, col) is located in Tensors at linear position (row + Nrows*col). Each PEPO tensor has the
    Shape {Dleft, Dup, Dright, Ddown, dup, ddown}. */
  Tensor<T>* Tensors;

/// Returns Shape to position for open boundary conditions.
/** Given a row position positionRow and a column position positionCol in this PEPO with open boundary
    conditions, the correct Shape of the corresponding Tensor is returned.
    This PEPO is not changed.
    \param positionRow input: unsigned int, the row position
    \param positionCol input: unsigned int, the column position
    \param Shape output: vector<unsigned int>&, the resulting Shape, on input it must fulfill
                          Shape.size()==6, on output it will contain {Dleft, Dup, Dright, Ddown, dup, ddown} */
  void getOpenBCShape(unsigned int positionRow, unsigned int positionCol,
                      vector<unsigned int>& Shape);

};

template<class T> PEPO<T>::PEPO()
{
 this->N = vector<unsigned int>(2);
 this->N[0] = 0; this->N[1] = 0;
 this->D = 0;
 this->Tensors = 0;
}

template<class T> PEPO<T>::PEPO(const string& BC0, unsigned int Nrows, unsigned int Ncols,
                                const Matrix<unsigned int>& d0, unsigned int D0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (Nrows < 2) || (Ncols < 2) ||
     (d0.getDim0() != Nrows) || (d0.getDim1() != Ncols) || (D0 == 0))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> PEPO<T>::" <<
          "PEPO(const string& BC0, unsigned int Nrows, unsigned int Ncols, " <<
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
            "template<class T> PEPO<T>::" <<
            "PEPO(const string& BC0, unsigned int Nrows, unsigned int Ncols, " <<
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
 vector<unsigned int> Shape(6);
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
    Shape[4] = this->d(i, j); Shape[5] = this->d(i, j);
    position = i + j*this->N[0];
    this->Tensors[position] = Tensor<T>(Shape);
   }
  }
 }
}

template<class T> PEPO<T>::PEPO(const string& BC0, unsigned int Nrows, unsigned int Ncols,
                                unsigned int d0, unsigned int D0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (Nrows < 2) || (Ncols < 2) || (d0 < 2) || (D0 == 0))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> PEPO<T>::" <<
          "PEPO(const string& BC0, unsigned int Nrows, unsigned int Ncols, unsigned int d0, " <<
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
 vector<unsigned int> Shape(6);
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
  Shape[4] = this->d(0, 0); Shape[5] = this->d(0, 0);
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

template<class T> PEPO<T>::PEPO(const PEPO<T>& PEPO0)
{
 this->BC = PEPO0.BC;
 this->N = PEPO0.N;
 this->d = PEPO0.d;
 this->D = PEPO0.D;
 unsigned int numTensors = this->N[0] * this->N[1];
 this->Tensors = new Tensor<T>[numTensors];
 for (int i = 0; i < numTensors; i++)
 {
  this->Tensors[i] = PEPO0.Tensors[i];
 }
}

template<class T> PEPO<T>::~PEPO()
{
 delete[] this->Tensors;
}

template<class T> PEPO<T>& PEPO<T>::operator=(const PEPO<T>& PEPO0)
{
 if (this != &PEPO0)
 {
  this->BC = PEPO0.BC;
  this->N = PEPO0.N;
  this->d = PEPO0.d;
  this->D = PEPO0.D;
  delete[] this->Tensors;
  unsigned int numTensors = this->N[0] * this->N[1];
  this->Tensors = new Tensor<T>[numTensors];
  for (int i = 0; i < numTensors; i++)
  {
   this->Tensors[i] = PEPO0.Tensors[i];
  }
 }
 return *this;
}

template<class T> inline void PEPO<T>::set(unsigned int positionRow, unsigned int positionCol,
                                           const Tensor<T>& Tensor0)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (positionRow >= this->N[0]) || (positionCol >= this->N[1]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> inline void PEPO<T>::" <<
          "set(unsigned int positionRow, unsigned int positionCol, const Tensor<T>& Tensor0): " <<
          "((this->N[0] == 0) || (positionRow >= this->N[0]) || (positionCol >= this->N[1]))." << endl;
  exit(1);
 }
 vector<unsigned int> Shape0;
 Tensor0.getShape(Shape0);
 if (this->BC == "open")
 {
  vector<unsigned int> Shape1(6);
  this->getOpenBCShape(positionRow, positionCol, Shape1);
  if ((Shape0[0] > Shape1[0]) || (Shape0[1] > Shape1[1]) || (Shape0[2] > Shape1[2]) ||
      (Shape0[3] > Shape1[3]) || (Shape0[4] != this->d(positionRow, positionCol)) || (Shape0[5] != this->d(positionRow, positionCol)))
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> inline void PEPO<T>::" <<
           "set(unsigned int positionRow, unsigned int positionCol, const Tensor<T>& Tensor0): " <<
           "Tensor0 has an incorrect shape." << endl;
   exit(1);
  }
 }
 else if (this->BC == "periodic")
 {
  if ((Shape0[0] > this->D) || (Shape0[1] > this->D) || (Shape0[2] > this->D) ||
      (Shape0[3] > this->D) || (Shape0[4] != this->d(positionRow, positionCol)) || (Shape0[5] != this->d(positionRow, positionCol)))
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> inline void PEPO<T>::" <<
           "set(unsigned int positionRow, unsigned int positionCol, const Tensor<T>& Tensor0): " <<
           "Tensor0 has an incorrect shape." << endl;
   exit(1);
  }
 }
#endif
 this->Tensors[positionRow + this->N[0]*positionCol] = Tensor0;
}

template<class T> inline void PEPO<T>::get(unsigned int positionRow, unsigned int positionCol,
                                           Tensor<T>& Tensor0) const
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (positionRow >= this->N[0]) || (positionCol >= this->N[1]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> inline void PEPO<T>::" <<
          "get(unsigned int positionRow, unsigned int positionCol, Tensor<T>& Tensor0) const: " <<
          "((this->N[0] == 0) || (positionRow >= this->N[0]) || (positionCol >= this->N[1]))." << endl;
  exit(1);
 }
#endif
 Tensor0 = this->Tensors[positionRow + this->N[0]*positionCol];
}

template<class T> void PEPO<T>::write(const string& FileName) const
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
          "template<class T> void PEPO<T>::" <<
          "write(const string& FileName) const: " <<
          "Binary output file " << FileName << " could not be opened." << endl;
  exit(1);
 }
}

template<class T> void PEPO<T>::read(const string& FileName)
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
          "template<class T> void PEPO<T>::" <<
          "read(const string& FileName): " <<
          "Binary input file " << FileName << " could not be opened." << endl;
  exit(1);
 }
}

template<class T> void PEPO<T>::getMatrix(Matrix<T>& Matrix0) const
{}

template<class T> void PEPO<T>::fillRandomly(Matrix<unsigned int>& Seed)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (Seed.getDim0() != this->N[0]) || (Seed.getDim1() != this->N[1]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPO<T>::" <<
          "fillRandomly(const Matrix<unsigned int>& Seed): " <<
          "((this->N[0] == 0) || (Seed.getDim0() != this->N[0]) || " <<
           "(Seed.getDim1() != this->N[1]))." << endl;
  exit(1);
 }
#endif
 for (int col = 0; col < this->N[1]; col++)
 {
  for (int row = 0; row < this->N[0]; row++)
  {
   this->Tensors[row+this->N[0]*col].fillRandomly(Seed(row, col));
  }
 }
}

template<class T> void PEPO<T>::fillZeroes()
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPO<T>::" <<
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

template<class T> void PEPO<T>::setOne()
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPO<T>::" <<
          "setOne(): " <<
          "(this->N[0] == 0)." << endl;
  exit(1);
 }
#endif
 this->fillZeroes();
 vector<unsigned int> Index(6);
 Index[0] = 0; Index[1] = 0; Index[2] = 0; Index[3] = 0;
 T element = 1.0;
 for (int col = 0; col < this->N[1]; col++)
 {
  for (int row = 0; row < this->N[0]; row++)
  {
   for (int j = 0; j < this->d(row, col); j++)
   {
    Index[4] = j; Index[5] = j;
    this->Tensors[row+this->N[0]*col].set(Index, element);
   }
  }
 }
}

template<class T> void PEPO<T>::setConcatenated(unsigned int M)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (M == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPO<T>::" <<
          "setConcatenated(unsigned int M): " <<
          "((this->N[0] == 0) || (M == 0))." << endl;
  exit(1);
 }
#endif
 if (M == 1)
  return;
 unsigned int positionCol, positionRow;
 vector<unsigned int> Index(6), N0, Shape(6);
 Tensor<T> Tensor0;
 Tensor<T>* TensorsConc;
 Matrix<unsigned int> d0;
// 1. get N and d of Concatenated PEPO:
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
// 2. get tensors of Concatenated PEPO:
 TensorsConc = new Tensor<T>[this->N[0]*this->N[1]];
 Shape[0] = this->D; Shape[1] = this->D; Shape[2] = this->D; Shape[3] = this->D;
 for (int col0 = 0; col0 < N0[1]; col0++){
  for (int col1 = 0; col1 < M; col1++){
   positionCol = col0*M+col1;
   for (int row0 = 0; row0 < N0[0]; row0++){
    for (int row1 = 0; row1 < M; row1++){
     positionRow = row0*M+row1;
     if (this->BC == "open")
     {
      this->getOpenBCShape(positionRow, positionCol, Shape);
     }
     else if (this->BC == "periodic")
     {
      Shape[4] = this->d(positionRow, positionCol);
      Shape[5] = this->d(positionRow, positionCol);
     }
// 2.a) ((row1 == 0) && (col1 == 0)): physical tensors taken from this PEPO
     if ((row1 == 0) && (col1 == 0))
     {
      TensorsConc[positionRow+this->N[0]*positionCol] = this->Tensors[row0+N0[0]*col0];
      TensorsConc[positionRow+this->N[0]*positionCol].setShape(Shape);
     }
// 2.b) !((row1 == 0) && (col1 == 0)): auxiliary tensors set as Deltas
     else
     {
      Tensor0 = Tensor<T>(Shape);
      Tensor0.fillZeroes();
      Index[4] = 0; Index[5] = 0;
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
// 3. write Concatenated PEPO into this PEPO:
 delete[] this->Tensors;
 this->Tensors = TensorsConc;
}

template<class T> void PEPO<T>::multiply(T element)
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPO<T>::" <<
          "multiply(T element): " <<
          "(this->N[0] == 0)." << endl;
  exit(1);
 }
#endif
 this->Tensors[0].multiply(element);
}

template<class T> void PEPO<T>::getOpenBCShape(unsigned int positionRow, unsigned int positionCol,
                                               vector<unsigned int>& Shape)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (this->BC != "open") || (positionRow >= this->N[0]) ||
     (positionCol >= this->N[1]) || (Shape.size() != 6))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void PEPO<T>::" <<
          "getOpenBCShape(unsigned int positionRow, unsigned int positionCol, " <<
                         "vector<unsigned int>& Shape) const: " <<
          "((this->N[0] == 0) || (this->BC != open) || (positionRow >= this->N[0]) || " <<
           "(positionCol >= this->N[1]) || (Shape.size() != 6))." << endl;
  exit(1);
 }
#endif
 Shape[0] = this->D; Shape[1] = this->D; Shape[2] = this->D; Shape[3] = this->D;
 Shape[4] = this->d(positionRow, positionCol); Shape[5] = this->d(positionRow, positionCol);
 if (positionRow == 0)
  Shape[1] = 1;
 if (positionRow == this->N[0]-1)
  Shape[3] = 1;
 if (positionCol == 0)
  Shape[0] = 1;
 if (positionCol == this->N[1]-1)
  Shape[2] = 1;
}
