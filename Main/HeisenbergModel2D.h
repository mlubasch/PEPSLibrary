/// Template class HeisenbergModel2D implements quantum Heisenberg model in 2D.
/** The template class HeisenbergModel2D implements the quantum Heisenberg model
         H = J \sum_{<l,m>} \vec{S}_{l} \vec{S}_{m}
           = J \sum_{<l,m>} 0.5 * (S_{l}^{+} S_{m}^{-} + S_{l}^{-} S_{m}^{+}) + S_{l}^{z} S_{m}^{z}
    in 2D, where S^{+}, S^{-}, and S^{z} are spin operators.
    It has the parameter J.
    This template class is child class to the abstract template classes Hamiltonian2D and Operator2D.
    \param Representation string, the representation, is "Interactions", "PEPO" or "Matrix"
    \param BC string, the boundary conditions, is "open" or "periodic"
    \param N vector<unsigned int>(2), the number of sites, fulfills N.size()==2
    \param d Matrix<unsigned int>(N[0], N[1]), the physical dimensions, fulfills
             d.getDim0()==N[0] and d.getDim1()==N[1] and is fixed to 2 everywhere
    \param Parameters vector<T>, the parameters {J}
    \param timeDependent bool, the time dependence
    \param TimeFunctions vector<PointerToFunction>, the functions defining the time dependence
    \param time double, the time
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

template<class T> class HeisenbergModel2D : public Hamiltonian2D<T>
{
 typedef T (*PointerToFunction)(double time);

 public:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  HeisenbergModel2D();

/// Constructor for time-independent HeisenbergModel2D.
/** This constructor initializes a time-independent HeisenbergModel2D with specific BC, N, and Parameters.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param Nrows input: unsigned int, the number of rows, must be > 1
    \param Ncols input: unsigned int, the number of columns, must be > 1
    \param Parameters0 input: const vector<T>&, the parameters {J},
                              must fulfill Parameters0.size()==1 */
  HeisenbergModel2D(const string& BC0, unsigned int Nrows, unsigned int Ncols,
                    const vector<T>& Parameters0);

/// Constructor for time-dependent HeisenbergModel2D.
/** This constructor initializes a time-dependent HeisenbergModel2D with specific BC, N, TimeFunctions and time.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param Nrows input: unsigned int, the number of rows, must be > 1
    \param Ncols input: unsigned int, the number of columns, must be > 1
    \param TimeFunctions0 input: const vector<PointerToFunction>&, the time functions,
                                 must fulfill TimeFunctions0.size()==1
    \param time0 input: double, the time */
  HeisenbergModel2D(const string& BC0, unsigned int Nrows, unsigned int Ncols,
                    const vector<PointerToFunction>& TimeFunctions0, double time0);

/// Standard copy constructor.
/** The standard copy constructor copies the input HeisenbergModel2D into this.
    \param HeisenbergModel2D0 input: const HeisenbergModel2D<T>&, to be copied into this
    \sa HeisenbergModel2D<T>& operator=(const HeisenbergModel2D<T>& HeisenbergModel2D0) */
  HeisenbergModel2D(const HeisenbergModel2D<T>& HeisenbergModel2D0);

/// Standard destructor.
/** The standard destructor deletes the elements of HeisenbergModel2D. */
  ~HeisenbergModel2D();

/// Assigns HeisenbergModel2D to this.
/** The operator= allows to assign HeisenbergModel2D0 to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side HeisenbergModel2D0.
    \param HeisenbergModel2D0 input: const HeisenbergModel2D<T>&, to be copied into this
    \return HeisenbergModel2D<T>&, a reference to the new this
    \sa HeisenbergModel2D(const HeisenbergModel2D<T>& HeisenbergModel2D0) */
  HeisenbergModel2D<T>& operator=(const HeisenbergModel2D<T>& HeisenbergModel2D0);

/// Returns Interactions representation.
/** This function returns the Interactions of this HeisenbergModel2D.
    The Interactions are implemented as a vector of row positions, a vector of column positions, and a vector of
    matrices, where each vector entry corresponds to one Interaction, i.e. one term of the sum making up this
    HeisenbergModel2D.
    \param PositionsRow output: vector< vector<unsigned int> >&, the row positions of the Interactions
    \param PositionsCol output: vector< vector<unsigned int> >&, the column positions of the Interactions
    \param Interactions output: vector< vector< Matrix<T> > >&, the Interactions */
  void getInteractions(vector< vector<unsigned int> >& PositionsRow,
                       vector< vector<unsigned int> >& PositionsCol,
                       vector< vector< Matrix<T> > >& Interactions) const;

/// Returns PEPO representation.
/** This function returns this HeisenbergModel2D as a PEPO.
    \param PEPO0 output: PEPO<T>&, the PEPO representing this HeisenbergModel2D */
  void getPEPO(PEPO<T>& PEPO0) const;

/// Returns Matrix representation.
/** This function returns this HeisenbergModel2D as a matrix.
    \param Matrix0 output: Matrix<T>&, the matrix representing this HeisenbergModel2D,
                           must fulfill ((Matrix0.getDim0()==2^{Nrows*Ncols}) && (Matrix0.getDim1()==2^{Nrows*Ncols})),
                           will be of type "hermitian" */
  void getMatrix(Matrix<T>& Matrix0) const;

/// Returns two-site Hamiltonian MPOs.
/** This function returns the two-site Hamiltonian MPOs corresponding to the two-site Trotter gates.
    \param HMPOs output: vector< MPO<T> >&, the two-site Hamiltonian MPOs */
  void getHMPOs(vector< MPO<T> >& HMPOs) const;

/// Returns two-site Trotter gates.
/** This function returns the two-site Trotter gates.
    \param RealImaginaryTE input: const string&, must be "real" or "imaginary"
    \param timeStep input: double, the time step
    \param TEMPOs output: vector< MPO<T> >&, the two-site Trotter gates */
  void getTEMPOs(const string& RealImaginaryTE, double timeStep, vector< MPO<T> >& TEMPOs);

/// Returns PEPOs for time evolution.
/** This function returns the PEPOs resulting from a Trotter decomposition of the evolution operator.
    \param RealImaginaryTE input: const string&, must be "real" or "imaginary"
    \param timeStep input: double, the time step
    \param TEPEPOs output: vector< PEPO<T> >&, the time evolution PEPOs */
  void getTEPEPOs(const string& RealImaginaryTE, double timeStep, vector< PEPO<T> >& TEPEPOs);

};

template<class T> HeisenbergModel2D<T>::HeisenbergModel2D()
{
 this->N = vector<unsigned int>(2);
 this->N[0] = 0; this->N[1] = 0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> HeisenbergModel2D<T>::HeisenbergModel2D(const string& BC0, unsigned int Nrows, unsigned int Ncols,
                                                          const vector<T>& Parameters0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (Nrows < 2) || (Ncols < 2) ||
     (Parameters0.size() != 1))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> HeisenbergModel2D<T>::" <<
          "HeisenbergModel2D(const string& BC0, unsigned int Nrows, unsigned int Ncols, " <<
                            "const vector<T>& Parameters0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (Nrows < 2) || (Ncols < 2) || " <<
           "(Parameters0.size() != 1))." << endl;
  exit(1);
 }
#endif
 this->Representation = "PEPO";
 this->BC = BC0;
 this->N = vector<unsigned int>(2);
 this->N[0] = Nrows; this->N[1] = Ncols;
 this->d = Matrix<unsigned int>(this->N[0], this->N[1]);
 for (int j = 0; j < this->N[1]; j++)
 {
  for (int i = 0; i < this->N[0]; i++)
  {
   this->d(i, j) = 2;
  }
 }
 this->Parameters = Parameters0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> HeisenbergModel2D<T>::HeisenbergModel2D(const string& BC0, unsigned int Nrows, unsigned int Ncols,
                                                          const vector<PointerToFunction>& TimeFunctions0,
                                                          double time0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (Nrows < 2) || (Ncols < 2) ||
     (TimeFunctions0.size() != 1))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> HeisenbergModel2D<T>::" <<
          "HeisenbergModel2D(const string& BC0, unsigned int Nrows, unsigned int Ncols, " <<
                            "const vector<PointerToFunction>& TimeFunctions0, double time0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (Nrows < 2) || (Ncols < 2) || " <<
           "(TimeFunctions0.size() != 1))." << endl;
  exit(1);
 }
#endif
 this->Representation = "PEPO";
 this->BC = BC0;
 this->N = vector<unsigned int>(2);
 this->N[0] = Nrows; this->N[1] = Ncols;
 this->d = Matrix<unsigned int>(this->N[0], this->N[1]);
 for (int j = 0; j < this->N[1]; j++)
 {
  for (int i = 0; i < this->N[0]; i++)
  {
   this->d(i, j) = 2;
  }
 }
 this->Parameters = vector<T>(1);
 this->timeDependent = true;
 this->TimeFunctions = TimeFunctions0;
 this->time = time0;
 for (int i = 0; i < this->Parameters.size(); i++)
  this->Parameters[i] = this->TimeFunctions[i](this->time);
}

template<class T> HeisenbergModel2D<T>::HeisenbergModel2D(const HeisenbergModel2D<T>& HeisenbergModel2D0)
{
 this->Representation = HeisenbergModel2D0.Representation;
 this->BC = HeisenbergModel2D0.BC;
 this->N = HeisenbergModel2D0.N;
 this->d = HeisenbergModel2D0.d;
 this->Parameters = HeisenbergModel2D0.Parameters;
 this->timeDependent = HeisenbergModel2D0.timeDependent;
 this->TimeFunctions = HeisenbergModel2D0.TimeFunctions;
 this->time = HeisenbergModel2D0.time;
}

template<class T> HeisenbergModel2D<T>::~HeisenbergModel2D() {}

template<class T> HeisenbergModel2D<T>& HeisenbergModel2D<T>::operator=(const HeisenbergModel2D<T>& HeisenbergModel2D0)
{
 if (this != &HeisenbergModel2D0)
 {
  this->Representation = HeisenbergModel2D0.Representation;
  this->BC = HeisenbergModel2D0.BC;
  this->N = HeisenbergModel2D0.N;
  this->d = HeisenbergModel2D0.d;
  this->Parameters = HeisenbergModel2D0.Parameters;
  this->timeDependent = HeisenbergModel2D0.timeDependent;
  this->TimeFunctions = HeisenbergModel2D0.TimeFunctions;
  this->time = HeisenbergModel2D0.time;
 }
 return *this;
}

template<class T> void HeisenbergModel2D<T>::getInteractions(vector< vector<unsigned int> >& PositionsRow,
                                                             vector< vector<unsigned int> >& PositionsCol,
                                                             vector< vector< Matrix<T> > >& Interactions) const
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void HeisenbergModel2D<T>::" <<
          "getInteractions(vector< vector<unsigned int> >& PositionsRow, " <<
                          "vector< vector<unsigned int> >& PositionsCol, " <<
                          "vector< vector< Matrix<T> > >& Interactions) const: " <<
          "(this->N[0] == 0)." << endl;
  exit(1);
 }
#endif
 unsigned int numInteractions = 3*(2*this->N[0]*this->N[1]-this->N[0]-this->N[1]);
 PositionsRow = vector< vector<unsigned int> >(numInteractions);
 PositionsCol = vector< vector<unsigned int> >(numInteractions);
 Interactions = vector< vector< Matrix<T> > >(numInteractions);
// parameters:
 T J = this->Parameters[0];
// spin operators:
// - S0 := 0.5*J*S^{+}, S1 := S^{-}:
 Matrix<T> S0(2, 2), S1(2, 2);
 S0.fillZeroes(); S1.fillZeroes();
 S0(0, 1) = 0.5*J; S1(1, 0) = 1.0;
// - S2 := 0.5*J*S^{-}, S3 := S^{+}:
 Matrix<T> S2(2, 2), S3(2, 2);
 S2.fillZeroes(); S3.fillZeroes();
 S2(1, 0) = 0.5*J; S3(0, 1) = 1.0;
// - S4 := J*S^{z}, S5 := S^{z}:
 Matrix<T> S4(2, 2), S5(2, 2);
 S4.fillZeroes(); S5.fillZeroes();
 S4(0, 0) = 0.5*J; S4(1, 1) = -0.5*J; S5(0, 0) = 0.5; S5(1, 1) = -0.5;
// Interaction vectors:
 vector< Matrix<T> > VectorS0S1(2), VectorS2S3(2), VectorS4S5(2);
 VectorS0S1[0] = S0; VectorS0S1[1] = S1;
 VectorS2S3[0] = S2; VectorS2S3[1] = S3;
 VectorS4S5[0] = S4; VectorS4S5[1] = S5;
// vertical Interactions:
 vector<unsigned int> VectorRow(2), VectorCol(2);
 unsigned int count = 0;
 for (int col = 0; col < this->N[1]; col++)
 {
  VectorCol[0] = col; VectorCol[1] = col;
  for (int row = 0; row < this->N[0]-1; row++)
  {
   VectorRow[0] = row; VectorRow[1] = row+1;
   PositionsRow[count] = VectorRow;
   PositionsCol[count] = VectorCol;
   Interactions[count] = VectorS0S1;
   PositionsRow[count+1] = VectorRow;
   PositionsCol[count+1] = VectorCol;
   Interactions[count+1] = VectorS2S3;
   PositionsRow[count+2] = VectorRow;
   PositionsCol[count+2] = VectorCol;
   Interactions[count+2] = VectorS4S5;
   count += 3;
  }
 }
// horizontal Interactions:
 for (int col = 0; col < this->N[1]-1; col++)
 {
  VectorCol[0] = col; VectorCol[1] = col+1;
  for (int row = 0; row < this->N[0]; row++)
  {
   VectorRow[0] = row; VectorRow[1] = row;
   PositionsRow[count] = VectorRow;
   PositionsCol[count] = VectorCol;
   Interactions[count] = VectorS0S1;
   PositionsRow[count+1] = VectorRow;
   PositionsCol[count+1] = VectorCol;
   Interactions[count+1] = VectorS2S3;
   PositionsRow[count+2] = VectorRow;
   PositionsCol[count+2] = VectorCol;
   Interactions[count+2] = VectorS4S5;
   count += 3;
  }
 }
}

template<class T> void HeisenbergModel2D<T>::getPEPO(PEPO<T>& PEPO0) const
{
 cerr << "The following function is not implemented: " <<
         "template<class T> void HeisenbergModel2D<T>::" <<
         "getPEPO(PEPO<T>& PEPO0) const." << endl;
 exit(1);
}

template<class T> void HeisenbergModel2D<T>::getMatrix(Matrix<T>& Matrix0) const
{
#ifdef DEBUG
 unsigned int dim = pow(2.0, int(this->N[0]*this->N[1]));
 if ((this->N[0] == 0) || (Matrix0.getDim0() != dim) || (Matrix0.getDim1() != dim))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void HeisenbergModel2D<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const: " <<
          "((this->N[0] == 0) || (Matrix0.getDim0() != 2^{Nrows*Ncols}) || " <<
           "(Matrix0.getDim1() != 2^{Nrows*Ncols}))." << endl;
  exit(1);
 }
#endif
 Matrix0.fillZeroes();
// parameters:
 T J = this->Parameters[0];
// Interactions:
// - Identity:
 Matrix<T> Identity(2, 2);
 Identity.fillZeroes();
 for (int i = 0; i < 2; i++)
  Identity(i, i) = 1.0;
// - S0 := S^{+}, S1 := S^{-}, S2 := S^{z}:
 Matrix<T> S0(2, 2), S1(2, 2), S2(2, 2);
 S0.fillZeroes(); S1.fillZeroes(); S2.fillZeroes();
 S0(0, 1) = 1.0; S1(1, 0) = 1.0; S2(0, 0) = 0.5; S2(1, 1) = -0.5;
// - S0S1 := S^{+}S^{-}, S1S0 := S^{-}S^{+}, S2S2 := S^{z}S^{z}:
 Matrix<T> S0S1(4, 4), S1S0(4, 4), S2S2(4, 4);
 S0.multiplyDirectProduct(S1, S0S1);
 S1.multiplyDirectProduct(S0, S1S0);
 S2.multiplyDirectProduct(S2, S2S2);
// S := J*(0.5*(S^{+}S^{-}+S^{-}S^{+})+S^{z}S^{z}):
 Matrix<T> S(4, 4);
 S = S0S1;
 S.add(S1S0);
 S.multiply(0.5);
 S.add(S2S2);
 S.multiply(J);
 Matrix<T> X, Y;
 unsigned int position;
// open boundary conditions:
 if (this->BC == "open")
 {
// Hamiltonian sum:
// 1. \sum_{<l,m>} S_{l,m}, vertical terms:
  for (int col = 0; col < this->N[1]; col++)
  {
   for (int row = 0; row < this->N[0]-1; row++)
   {
    position = row + this->N[0]*col;
    X = S;
    for (int j = position+2; j < this->N[0]*this->N[1]; j++)
    {
     X.multiplyDirectProduct(Identity, Y);
     X = Y;
    }
    for (int j = position-1; j >= 0; j--)
    {
     Identity.multiplyDirectProduct(X, Y);
     X = Y;
    }
    Matrix0.add(X);
   }
  }
// 2. \sum_{<l,m>} S_{l,m}, horizontal terms:
  for (int col = 0; col < this->N[1]-1; col++)
  {
   for (int row = 0; row < this->N[0]; row++)
   {
    position = row + this->N[0]*col;
// 2.a) J*0.5*S^{+}S^{-}:
    X = S0;
    X.multiply(J*0.5);
    for (int j = 1; j < this->N[0]; j++)
    {
     X.multiplyDirectProduct(Identity, Y);
     X = Y;
    }
    X.multiplyDirectProduct(S1, Y);
    X = Y;
    for (int j = position+this->N[0]+1; j < this->N[0]*this->N[1]; j++)
    {
     X.multiplyDirectProduct(Identity, Y);
     X = Y;
    }
    for (int j = position-1; j >= 0; j--)
    {
     Identity.multiplyDirectProduct(X, Y);
     X = Y;
    }
    Matrix0.add(X);
// 2.b) J*0.5*S^{-}S^{+}:
    X = S1;
    X.multiply(J*0.5);
    for (int j = 1; j < this->N[0]; j++)
    {
     X.multiplyDirectProduct(Identity, Y);
     X = Y;
    }
    X.multiplyDirectProduct(S0, Y);
    X = Y;
    for (int j = position+this->N[0]+1; j < this->N[0]*this->N[1]; j++)
    {
     X.multiplyDirectProduct(Identity, Y);
     X = Y;
    }
    for (int j = position-1; j >= 0; j--)
    {
     Identity.multiplyDirectProduct(X, Y);
     X = Y;
    }
    Matrix0.add(X);
// 2.c) J*S^{z}S^{z}:
    X = S2;
    X.multiply(J);
    for (int j = 1; j < this->N[0]; j++)
    {
     X.multiplyDirectProduct(Identity, Y);
     X = Y;
    }
    X.multiplyDirectProduct(S2, Y);
    X = Y;
    for (int j = position+this->N[0]+1; j < this->N[0]*this->N[1]; j++)
    {
     X.multiplyDirectProduct(Identity, Y);
     X = Y;
    }
    for (int j = position-1; j >= 0; j--)
    {
     Identity.multiplyDirectProduct(X, Y);
     X = Y;
    }
    Matrix0.add(X);
   }
  }
 }
// periodic boundary conditions:
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented for periodic BC: " <<
          "template<class T> void HeisenbergModel2D<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const." << endl;
  exit(1);
 }
 Matrix0.setType("hermitian");
}

template<class T> void HeisenbergModel2D<T>::getHMPOs(vector< MPO<T> >& HMPOs) const
{
 T element, J = this->Parameters[0];
 vector<unsigned int> Index4(4), Shape4(4);
 Tensor<T> TensorA;
 MPO<T> HMPO;
// 1. HMPO := J*(0.5(S^{+}S^{-}+S^{-}S^{+})+S^{z}S^{z}):
 HMPO = MPO<T>("open", 2, 2, 3);
// 1.1. position 0 in HMPO:
 Shape4[0] = 1; Shape4[1] = 3; Shape4[2] = 2; Shape4[3] = 2;
 TensorA = Tensor<T>(Shape4);
 TensorA.fillZeroes();
 Index4[0] = 0;
 element = J*0.5;
 Index4[1] = 0; Index4[2] = 1; Index4[3] = 0;
 TensorA.set(Index4, element);
 Index4[1] = 1; Index4[2] = 0; Index4[3] = 1;
 TensorA.set(Index4, element);
 Index4[1] = 2; Index4[2] = 0; Index4[3] = 0;
 TensorA.set(Index4, element);
 element = -J*0.5;
 Index4[2] = 1; Index4[3] = 1;
 TensorA.set(Index4, element);
 HMPO.set(0, TensorA);
// 1.2. position 1 in HMPO:
 Shape4[0] = 3; Shape4[1] = 1; Shape4[2] = 2; Shape4[3] = 2;
 TensorA = Tensor<T>(Shape4);
 TensorA.fillZeroes();
 Index4[1] = 0;
 element = 1.0;
 Index4[0] = 0; Index4[2] = 0; Index4[3] = 1;
 TensorA.set(Index4, element);
 Index4[0] = 1; Index4[2] = 1; Index4[3] = 0;
 TensorA.set(Index4, element);
 element = 0.5;
 Index4[0] = 2; Index4[2] = 0; Index4[3] = 0;
 TensorA.set(Index4, element);
 element = -0.5;
 Index4[2] = 1; Index4[3] = 1;
 TensorA.set(Index4, element);
 HMPO.set(1, TensorA);
// 2. HMPOs:
 HMPOs = vector< MPO<T> >(4);
 HMPOs[0] = HMPO; HMPOs[1] = HMPO; HMPOs[2] = HMPO; HMPOs[3] = HMPO;
}

template<class T> void HeisenbergModel2D<T>::getTEMPOs(const string& RealImaginaryTE, double timeStep,
                                                       vector< MPO<T> >& TEMPOs)
{
#ifdef DEBUG
 if ((RealImaginaryTE != "real") && (RealImaginaryTE != "imaginary"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void HeisenbergModel2D<T>::" <<
          "getTEMPOs(const string& RealImaginaryTE, double timeStep, " <<
                    "vector< MPO<T> >& TEMPOs): " <<
          "((RealImaginaryTE != real) && (RealImaginaryTE != imaginary))." << endl;
  exit(1);
 }
#endif
 TEMPOs = vector< MPO<T> >(4);
 T delta;
 if (RealImaginaryTE == "real")
  MathAuxiliary::convertComplex(complex<double>(0.0, 1.0)*T(timeStep), delta);
 else if (RealImaginaryTE == "imaginary")
  delta = T(timeStep);
 if ((RealImaginaryTE == "real") && (this->timeDependent))
  this->setTime(this->time+timeStep);
 T J = this->Parameters[0];
// S0 := S^{+}, S1 := S^{-}, S2 := S^{z}:
 Matrix<T> S0(2, 2), S1(2, 2), S2(2, 2);
 S0.fillZeroes(); S1.fillZeroes(); S2.fillZeroes();
 S0(0, 1) = 1.0; S1(1, 0) = 1.0; S2(0, 0) = 0.5; S2(1, 1) = -0.5;
// S0S1 := S^{+}S^{-}, S1S0 := S^{-}S^{+}, S2S2 := S^{z}S^{z}:
 Matrix<T> S0S1(4, 4), S1S0(4, 4), S2S2(4, 4);
 S0.multiplyDirectProduct(S1, S0S1);
 S1.multiplyDirectProduct(S0, S1S0);
 S2.multiplyDirectProduct(S2, S2S2);
// H := 0.5(S^{+}S^{-}+S^{-}S^{+})+S^{z}S^{z}:
 Matrix<T> H(4, 4);
 H = S0S1;
 H.add(S1S0);
 H.multiply(0.5);
 H.add(S2S2);
 H.multiply(J);
// get TensorLeft and TensorRight for exp(-delta*H):
 vector<unsigned int> ShapeLeft(4), ShapeRight(4);
 ShapeLeft[0] = 1; ShapeLeft[1] = 4; ShapeLeft[2] = 2; ShapeLeft[3] = 2;
 ShapeRight[0] = 4; ShapeRight[1] = 1; ShapeRight[2] = 2; ShapeRight[3] = 2;
 Tensor<T> TensorLeft(ShapeLeft), TensorRight(ShapeRight);
 H.setType("hermitian");
 H.twoBodyExponentSVD(delta, 2, 2, TensorLeft, TensorRight);
// TEMPO = exp(-delta*H):
 MPO<T> TEMPO("open", 2, 2, 4);
 TEMPO.set(0, TensorLeft);
 TEMPO.set(1, TensorRight);
 TEMPOs[0] = TEMPO; TEMPOs[1] = TEMPO; TEMPOs[2] = TEMPO; TEMPOs[3] = TEMPO;
 if ((RealImaginaryTE == "real") && (this->timeDependent))
  this->setTime(this->time-timeStep);
}

template<class T> void HeisenbergModel2D<T>::getTEPEPOs(const string& RealImaginaryTE, double timeStep,
                                                        vector< PEPO<T> >& TEPEPOs)
{
#ifdef DEBUG
 if ((RealImaginaryTE != "real") && (RealImaginaryTE != "imaginary"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void HeisenbergModel2D<T>::" <<
          "getTEPEPOs(const string& RealImaginaryTE, double timeStep, " <<
                     "vector< PEPO<T> >& TEPEPOs): " <<
          "((RealImaginaryTE != real) && (RealImaginaryTE != imaginary))." << endl;
  exit(1);
 }
#endif
 unsigned int col, Ncols = this->N[1], Nrows = this->N[0], row;
 vector<unsigned int> Index6(6), Shape3(3), Shape6(6);
 Tensor<T> OneTensor, TETensor0, TETensor1, TETensorB, TETensorL, TETensorR, TETensorT;
 vector< MPO<T> > TEMPOs;
// 0. initialize:
 TEPEPOs = vector< PEPO<T> >(4);
 for (int i = 0; i < 4; i++)
  TEPEPOs[i] = PEPO<T>(this->BC, Nrows, Ncols, this->d, 4);
// 1. get TEMPOs:
 this->getTEMPOs(RealImaginaryTE, timeStep, TEMPOs);
// 2. get TETensors:
 TEMPOs[0].get(0, TETensor0); TEMPOs[0].get(1, TETensor1);
 Shape3[0] = 4; Shape3[1] = 2; Shape3[2] = 2;
 TETensor0.reshape(Shape3);
 TETensor1.reshape(Shape3);
// 2. set TEPEPOs:
// 2.1 set vertical-even Trotter gates in TEPEPOs[0]:
 TETensorT = TETensor0;
 Shape6[0] = 1; Shape6[1] = 1; Shape6[2] = 1; Shape6[3] = 4; Shape6[4] = 2; Shape6[5] = 2;
 TETensorT.reshape(Shape6);
 TETensorB = TETensor1;
 Shape6[0] = 1; Shape6[1] = 4; Shape6[2] = 1; Shape6[3] = 1; Shape6[4] = 2; Shape6[5] = 2;
 TETensorB.reshape(Shape6);
 for (col = 0; col < Ncols; col++)
 {
  for (row = 0; row < Nrows-1; row+=2)
  {
   TEPEPOs[0].set(row, col, TETensorT);
   TEPEPOs[0].set(row+1, col, TETensorB);
  }
 }
 if (Nrows%2 == 1)
 {
  row = Nrows-1;
  Shape6[0] = 1; Shape6[1] = 1; Shape6[2] = 1; Shape6[3] = 1;
  for (col = 0; col < Ncols; col++)
  {
   Shape6[4] = this->d(row, col); Shape6[5] = this->d(row, col);
   OneTensor = Tensor<T>(Shape6);
   OneTensor.fillZeroes();
   Index6[0] = 0; Index6[1] = 0; Index6[2] = 0; Index6[3] = 0;
   for (int i = 0; i < this->d(row, col); i++)
   {
    Index6[4] = i; Index6[5] = i;
    OneTensor.set(Index6, 1.0);
   }
   TEPEPOs[0].set(row, col, OneTensor);
  }
 }
// 2.2 set vertical-odd Trotter gates in TEPEPOs[1]:
 row = 0;
 Shape6[0] = 1; Shape6[1] = 1; Shape6[2] = 1; Shape6[3] = 1;
 for (col = 0; col < Ncols; col++)
 {
  Shape6[4] = this->d(row, col); Shape6[5] = this->d(row, col);
  OneTensor = Tensor<T>(Shape6);
  OneTensor.fillZeroes();
  Index6[0] = 0; Index6[1] = 0; Index6[2] = 0; Index6[3] = 0;
  for (int i = 0; i < this->d(row, col); i++)
  {
   Index6[4] = i; Index6[5] = i;
   OneTensor.set(Index6, 1.0);
  }
  TEPEPOs[1].set(row, col, OneTensor);
 }
 for (col = 0; col < Ncols; col++)
 {
  for (row = 1; row < Nrows-1; row+=2)
  {
   TEPEPOs[1].set(row, col, TETensorT);
   TEPEPOs[1].set(row+1, col, TETensorB);
  }
 }
 if (Nrows%2 == 0)
 {
  row = Nrows-1;
  Shape6[0] = 1; Shape6[1] = 1; Shape6[2] = 1; Shape6[3] = 1;
  for (col = 0; col < Ncols; col++)
  {
   Shape6[4] = this->d(row, col); Shape6[5] = this->d(row, col);
   OneTensor = Tensor<T>(Shape6);
   OneTensor.fillZeroes();
   Index6[0] = 0; Index6[1] = 0; Index6[2] = 0; Index6[3] = 0;
   for (int i = 0; i < this->d(row, col); i++)
   {
    Index6[4] = i; Index6[5] = i;
    OneTensor.set(Index6, 1.0);
   }
   TEPEPOs[1].set(row, col, OneTensor);
  }
 }
// 2.3 set horizontal-even Trotter gates in TEPEPOs[2]:
 TETensorL = TETensor0;
 Shape6[0] = 1; Shape6[1] = 1; Shape6[2] = 4; Shape6[3] = 1; Shape6[4] = 2; Shape6[5] = 2;
 TETensorL.reshape(Shape6);
 TETensorR = TETensor1;
 Shape6[0] = 4; Shape6[1] = 1; Shape6[2] = 1; Shape6[3] = 1; Shape6[4] = 2; Shape6[5] = 2;
 TETensorR.reshape(Shape6);
 for (col = 0; col < Ncols-1; col+=2)
 {
  for (row = 0; row < Nrows; row++)
  {
   TEPEPOs[2].set(row, col, TETensorL);
   TEPEPOs[2].set(row, col+1, TETensorR);
  }
 }
 if (Ncols%2 == 1)
 {
  col = Ncols-1;
  Shape6[0] = 1; Shape6[1] = 1; Shape6[2] = 1; Shape6[3] = 1;
  for (row = 0; row < Nrows; row++)
  {
   Shape6[4] = this->d(row, col); Shape6[5] = this->d(row, col);
   OneTensor = Tensor<T>(Shape6);
   OneTensor.fillZeroes();
   Index6[0] = 0; Index6[1] = 0; Index6[2] = 0; Index6[3] = 0;
   for (int i = 0; i < this->d(row, col); i++)
   {
    Index6[4] = i; Index6[5] = i;
    OneTensor.set(Index6, 1.0);
   }
   TEPEPOs[2].set(row, col, OneTensor);
  }
 }
// 2.4 set horizontal-odd Trotter gates in TEPEPOs[3]:
 col = 0;
 Shape6[0] = 1; Shape6[1] = 1; Shape6[2] = 1; Shape6[3] = 1;
 for (row = 0; row < Nrows; row++)
 {
  Shape6[4] = this->d(row, col); Shape6[5] = this->d(row, col);
  OneTensor = Tensor<T>(Shape6);
  OneTensor.fillZeroes();
  Index6[0] = 0; Index6[1] = 0; Index6[2] = 0; Index6[3] = 0;
  for (int i = 0; i < this->d(row, col); i++)
  {
   Index6[4] = i; Index6[5] = i;
   OneTensor.set(Index6, 1.0);
  }
  TEPEPOs[3].set(row, col, OneTensor);
 }
 for (col = 1; col < Ncols-1; col+=2)
 {
  for (row = 0; row < Nrows; row++)
  {
   TEPEPOs[3].set(row, col, TETensorL);
   TEPEPOs[3].set(row, col+1, TETensorR);
  }
 }
 if (Ncols%2 == 0)
 {
  col = Ncols-1;
  Shape6[0] = 1; Shape6[1] = 1; Shape6[2] = 1; Shape6[3] = 1;
  for (row = 0; row < Nrows; row++)
  {
   Shape6[4] = this->d(row, col); Shape6[5] = this->d(row, col);
   OneTensor = Tensor<T>(Shape6);
   OneTensor.fillZeroes();
   Index6[0] = 0; Index6[1] = 0; Index6[2] = 0; Index6[3] = 0;
   for (int i = 0; i < this->d(row, col); i++)
   {
    Index6[4] = i; Index6[5] = i;
    OneTensor.set(Index6, 1.0);
   }
   TEPEPOs[3].set(row, col, OneTensor);
  }
 }
}
