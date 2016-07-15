/// Template class IsingModel2D implements quantum Ising model in 2D.
/** The template class IsingModel2D implements the quantum Ising model
         H = B \sum_{l} sx_{l} + J \sum_{<l,m>} sz_{l} sz_{m}
    in 2D, where sx and sz are Pauli matrices.
    It has the parameters magnetic field B and spin interaction J.
    This template class is child class to the abstract template classes Hamiltonian2D and Operator2D.
    \param Representation string, the representation, is "Interactions", "PEPO" or "Matrix"
    \param BC string, the boundary conditions, is "open" or "periodic"
    \param N vector<unsigned int>(2), the number of sites, fulfills N.size()==2
    \param d Matrix<unsigned int>(N[0], N[1]), the physical dimensions, fulfills
             d.getDim0()==N[0] and d.getDim1()==N[1] and is fixed to 2 everywhere
    \param Parameters vector<T>, the parameters {B, J}
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

template<class T> class IsingModel2D : public Hamiltonian2D<T>
{
 typedef T (*PointerToFunction)(double time);

 public:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  IsingModel2D();

/// Constructor for time-independent IsingModel2D.
/** This constructor initializes a time-independent IsingModel2D with specific BC, N, and Parameters.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param Nrows input: unsigned int, the number of rows, must be > 1
    \param Ncols input: unsigned int, the number of columns, must be > 1
    \param Parameters0 input: const vector<T>&, the parameters {B, J},
                              must fulfill Parameters0.size()==2 */
  IsingModel2D(const string& BC0, unsigned int Nrows, unsigned int Ncols,
               const vector<T>& Parameters0);

/// Constructor for time-dependent IsingModel2D.
/** This constructor initializes a time-dependent IsingModel2D with specific BC, N, TimeFunctions and time.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param Nrows input: unsigned int, the number of rows, must be > 1
    \param Ncols input: unsigned int, the number of columns, must be > 1
    \param TimeFunctions0 input: const vector<PointerToFunction>&, the time functions,
                                 must fulfill TimeFunctions0.size()==2
    \param time0 input: double, the time */
  IsingModel2D(const string& BC0, unsigned int Nrows, unsigned int Ncols,
               const vector<PointerToFunction>& TimeFunctions0, double time0);

/// Standard copy constructor.
/** The standard copy constructor copies the input IsingModel2D into this.
    \param IsingModel2D0 input: const IsingModel2D<T>&, to be copied into this
    \sa IsingModel2D<T>& operator=(const IsingModel2D<T>& IsingModel2D0) */
  IsingModel2D(const IsingModel2D<T>& IsingModel2D0);

/// Standard destructor.
/** The standard destructor deletes the elements of IsingModel2D. */
  ~IsingModel2D();

/// Assigns IsingModel2D to this.
/** The operator= allows to assign IsingModel2D0 to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side IsingModel2D0.
    \param IsingModel2D0 input: const IsingModel2D<T>&, to be copied into this
    \return IsingModel2D<T>&, a reference to the new this
    \sa IsingModel2D(const IsingModel2D<T>& IsingModel2D0) */
  IsingModel2D<T>& operator=(const IsingModel2D<T>& IsingModel2D0);

/// Returns Interactions representation.
/** This function returns the Interactions of this IsingModel2D.
    The Interactions are implemented as a vector of row positions, a vector of column positions, and a vector of
    matrices, where each vector entry corresponds to one Interaction, i.e. one term of the sum making up this
    IsingModel2D.
    \param PositionsRow output: vector< vector<unsigned int> >&, the row positions of the Interactions
    \param PositionsCol output: vector< vector<unsigned int> >&, the column positions of the Interactions
    \param Interactions output: vector< vector< Matrix<T> > >&, the Interactions */
  void getInteractions(vector< vector<unsigned int> >& PositionsRow,
                       vector< vector<unsigned int> >& PositionsCol,
                       vector< vector< Matrix<T> > >& Interactions) const;

/// Returns PEPO representation.
/** This function returns this IsingModel2D as a PEPO.
    \param PEPO0 output: PEPO<T>&, the PEPO representing this IsingModel2D */
  void getPEPO(PEPO<T>& PEPO0) const;

/// Returns Matrix representation.
/** This function returns this IsingModel2D as a matrix.
    \param Matrix0 output: Matrix<T>&, the matrix representing this IsingModel2D,
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

template<class T> IsingModel2D<T>::IsingModel2D()
{
 this->N = vector<unsigned int>(2);
 this->N[0] = 0; this->N[1] = 0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> IsingModel2D<T>::IsingModel2D(const string& BC0, unsigned int Nrows, unsigned int Ncols,
                                                const vector<T>& Parameters0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (Nrows < 2) || (Ncols < 2) ||
     (Parameters0.size() != 2))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> IsingModel2D<T>::" <<
          "IsingModel2D(const string& BC0, unsigned int Nrows, unsigned int Ncols, " <<
                       "const vector<T>& Parameters0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (Nrows < 2) || (Ncols < 2) || " <<
           "(Parameters0.size() != 2))." << endl;
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

template<class T> IsingModel2D<T>::IsingModel2D(const string& BC0, unsigned int Nrows, unsigned int Ncols,
                                                const vector<PointerToFunction>& TimeFunctions0,
                                                double time0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (Nrows < 2) || (Ncols < 2) ||
     (TimeFunctions0.size() != 2))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> IsingModel2D<T>::" <<
          "IsingModel2D(const string& BC0, unsigned int Nrows, unsigned int Ncols, " <<
                       "const vector<PointerToFunction>& TimeFunctions0, double time0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (Nrows < 2) || (Ncols < 2) || " <<
           "(TimeFunctions0.size() != 2))." << endl;
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
 this->Parameters = vector<T>(2);
 this->timeDependent = true;
 this->TimeFunctions = TimeFunctions0;
 this->time = time0;
 for (int i = 0; i < this->Parameters.size(); i++)
  this->Parameters[i] = this->TimeFunctions[i](this->time);
}

template<class T> IsingModel2D<T>::IsingModel2D(const IsingModel2D<T>& IsingModel2D0)
{
 this->Representation = IsingModel2D0.Representation;
 this->BC = IsingModel2D0.BC;
 this->N = IsingModel2D0.N;
 this->d = IsingModel2D0.d;
 this->Parameters = IsingModel2D0.Parameters;
 this->timeDependent = IsingModel2D0.timeDependent;
 this->TimeFunctions = IsingModel2D0.TimeFunctions;
 this->time = IsingModel2D0.time;
}

template<class T> IsingModel2D<T>::~IsingModel2D() {}

template<class T> IsingModel2D<T>& IsingModel2D<T>::operator=(const IsingModel2D<T>& IsingModel2D0)
{
 if (this != &IsingModel2D0)
 {
  this->Representation = IsingModel2D0.Representation;
  this->BC = IsingModel2D0.BC;
  this->N = IsingModel2D0.N;
  this->d = IsingModel2D0.d;
  this->Parameters = IsingModel2D0.Parameters;
  this->timeDependent = IsingModel2D0.timeDependent;
  this->TimeFunctions = IsingModel2D0.TimeFunctions;
  this->time = IsingModel2D0.time;
 }
 return *this;
}

template<class T> void IsingModel2D<T>::getInteractions(vector< vector<unsigned int> >& PositionsRow,
                                                        vector< vector<unsigned int> >& PositionsCol,
                                                        vector< vector< Matrix<T> > >& Interactions) const
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void IsingModel2D<T>::" <<
          "getInteractions(vector< vector<unsigned int> >& PositionsRow, " <<
                          "vector< vector<unsigned int> >& PositionsCol, " <<
                          "vector< vector< Matrix<T> > >& Interactions) const: " <<
          "(this->N[0] == 0)." << endl;
  exit(1);
 }
#endif
 unsigned int numInteractions = 3*this->N[0]*this->N[1]-this->N[0]-this->N[1];
 PositionsRow = vector< vector<unsigned int> >(numInteractions);
 PositionsCol = vector< vector<unsigned int> >(numInteractions);
 Interactions = vector< vector< Matrix<T> > >(numInteractions);
// parameters:
 T B = this->Parameters[0], J = this->Parameters[1];
// Pauli matrices:
 Matrix<T> sx(2, 2), sz0(2, 2), sz1(2, 2);
 sx.fillZeroes();
 sx(0, 1) = B; sx(1, 0) = B;
 sz0.fillZeroes();
 sz0(0, 0) = J; sz0(1, 1) = -J;
 sz1.fillZeroes();
 sz1(0, 0) = 1.0; sz1(1, 1) = -1.0;
// Interaction vectors:
 vector< Matrix<T> > VectorSx(1), VectorSzSz(2);
 VectorSx[0] = sx;
 VectorSzSz[0] = sz0; VectorSzSz[1] = sz1;
// sx Interactions:
 vector<unsigned int> VectorRow(1), VectorCol(1);
 unsigned int count = 0;
 for (int col = 0; col < this->N[1]; col++)
 {
  VectorCol[0] = col;
  for (int row = 0; row < this->N[0]; row++)
  {
   VectorRow[0] = row;
   PositionsRow[count] = VectorRow;
   PositionsCol[count] = VectorCol;
   Interactions[count] = VectorSx;
   count++;
  }
 }
// vertical sz-sz Interactions:
 VectorRow = vector<unsigned int>(2);
 VectorCol = vector<unsigned int>(2);
 for (int col = 0; col < this->N[1]; col++)
 {
  VectorCol[0] = col; VectorCol[1] = col;
  for (int row = 0; row < this->N[0]-1; row++)
  {
   VectorRow[0] = row; VectorRow[1] = row+1;
   PositionsRow[count] = VectorRow;
   PositionsCol[count] = VectorCol;
   Interactions[count] = VectorSzSz;
   count++;
  }
 }
// horizontal sz-sz Interactions:
 for (int col = 0; col < this->N[1]-1; col++)
 {
  VectorCol[0] = col; VectorCol[1] = col+1;
  for (int row = 0; row < this->N[0]; row++)
  {
   VectorRow[0] = row; VectorRow[1] = row;
   PositionsRow[count] = VectorRow;
   PositionsCol[count] = VectorCol;
   Interactions[count] = VectorSzSz;
   count++;
  }
 }
}

template<class T> void IsingModel2D<T>::getPEPO(PEPO<T>& PEPO0) const
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void IsingModel2D<T>::" <<
          "getPEPO(PEPO<T>& PEPO0) const: " <<
          "(this->N[0] == 0)." << endl;
  exit(1);
 }
#endif
// PEPO operators O^{i}:
 vector<unsigned int> ShapeO(3);
 ShapeO[0] = 3; ShapeO[1] = this->d(0, 0); ShapeO[2] = this->d(0, 0);
 Tensor<T> O(ShapeO);
 O.fillZeroes();
 vector<unsigned int> IndexO(3);
// O^{0} = Identity:
 T element = 1.0;
 IndexO[0] = 0;
 IndexO[1] = 0; IndexO[2] = 0;
 O.set(IndexO, element);
 IndexO[1] = 1; IndexO[2] = 1;
 O.set(IndexO, element);
// O^{1} = sx:
 IndexO[0] = 1;
 IndexO[1] = 1; IndexO[2] = 0;
 O.set(IndexO, element);
 IndexO[1] = 0; IndexO[2] = 1;
 O.set(IndexO, element);
// O^{2} = sz:
 IndexO[0] = 2;
 IndexO[1] = 0; IndexO[2] = 0;
 O.set(IndexO, element);
 element = -1.0;
 IndexO[1] = 1; IndexO[2] = 1;
 O.set(IndexO, element);
 unsigned int D = 4;
 PEPO0 = PEPO<T>(this->BC, this->N[0], this->N[1], this->d, D);
// parameters:
 T B = this->Parameters[0], J = this->Parameters[1];
// open boundary conditions:
 if (this->BC == "open")
 {
// general A:
  vector<unsigned int> ShapeA(5), IndexA(5);
  ShapeA[0] = D; ShapeA[1] = D; ShapeA[2] = D; ShapeA[3] = D; ShapeA[4] = 3;
  Tensor<T> A(ShapeA);
  A.fillZeroes();
// for O^{0}=Identity:
  IndexA[4] = 0;
  element = 1.0;
  IndexA[0] = 0; IndexA[1] = 0; IndexA[2] = 0; IndexA[3] = 0;
  A.set(IndexA, element);
  IndexA[0] = 0; IndexA[1] = 2; IndexA[2] = 3; IndexA[3] = 2;
  A.set(IndexA, element);
  IndexA[0] = 2; IndexA[1] = 0; IndexA[2] = 2; IndexA[3] = 3;
  A.set(IndexA, element);
  IndexA[0] = 3; IndexA[1] = 3; IndexA[2] = 3; IndexA[3] = 3;
  A.set(IndexA, element);
// for O^{1}=sx:
  IndexA[4] = 1;
  element = B;
  IndexA[0] = 0; IndexA[1] = 0; IndexA[2] = 2; IndexA[3] = 2;
  A.set(IndexA, element);
// for O^{2}=sz:
  IndexA[4] = 2;
  element = J;
  IndexA[0] = 0; IndexA[1] = 0; IndexA[2] = 2; IndexA[3] = 1;
  A.set(IndexA, element);
  IndexA[0] = 0; IndexA[1] = 0; IndexA[2] = 1; IndexA[3] = 2;
  A.set(IndexA, element);
  element = 1.0;
  IndexA[0] = 0; IndexA[1] = 1; IndexA[2] = 3; IndexA[3] = 2;
  A.set(IndexA, element);
  IndexA[0] = 1; IndexA[1] = 0; IndexA[2] = 2; IndexA[3] = 3;
  A.set(IndexA, element);
  vector<unsigned int> IndexA0(1), IndexO0(1);
  IndexA0[0] = 4; IndexO0[0] = 0;
  A.contract(IndexA0, O, IndexO0);
  Tensor<T> AC(A);
// on the boundary of PEPO0 project A according to its position with projectors P0 := <0|,
// P1 := <0|+<2|+<3|, and P2 := <2|+<3|:
  vector<unsigned int> ShapeP(2), IndexP(2);
  ShapeP[0] = 1; ShapeP[1] = D;
  Tensor<T> P0(ShapeP), P1(ShapeP), P2(ShapeP);
  P0.fillZeroes(); P1.fillZeroes(); P2.fillZeroes();
  element = 1.0;
  IndexP[0] = 0; IndexP[1] = 0;
  P0.set(IndexP, element);
  P1.set(IndexP, element);
  IndexP[1] = 2;
  P1.set(IndexP, element);
  P2.set(IndexP, element);
  IndexP[1] = 3;
  P1.set(IndexP, element);
  P2.set(IndexP, element);
  Tensor<T> P0C(P0), P1C(P1), P2C(P2);
// put A into PEPO0:
  vector<unsigned int> IndicesA(1), IndicesP(1), OrderA(6);
  IndicesP[0] = 1;
  for (int col = 0; col < this->N[1]; col++)
  {
   for (int row = 0; row < this->N[0]; row++)
   {
    A = AC;
    if (col == 0)
    {
     IndicesA[0] = 0;
     A.contract(IndicesA, P0, IndicesP);
     P0 = P0C;
     OrderA[0] = 5; OrderA[1] = 0; OrderA[2] = 1; OrderA[3] = 2; OrderA[4] = 3; OrderA[5] = 4;
     A.permute(OrderA);
    }
    if (row == 0)
    {
     IndicesA[0] = 1;
     A.contract(IndicesA, P0, IndicesP);
     P0 = P0C;
     OrderA[0] = 0; OrderA[1] = 5; OrderA[2] = 1; OrderA[3] = 2; OrderA[4] = 3; OrderA[5] = 4;
     A.permute(OrderA);
    }
    if ((col == this->N[1]-1) && (row != this->N[0]-1))
    {
     IndicesA[0] = 2;
     A.contract(IndicesA, P1, IndicesP);
     P1 = P1C;
     OrderA[0] = 0; OrderA[1] = 1; OrderA[2] = 5; OrderA[3] = 2; OrderA[4] = 3; OrderA[5] = 4;
     A.permute(OrderA);
    }
    if ((row == this->N[0]-1) && (col != this->N[1]-1))
    {
     IndicesA[0] = 3;
     A.contract(IndicesA, P1, IndicesP);
     P1 = P1C;
     OrderA[0] = 0; OrderA[1] = 1; OrderA[2] = 2; OrderA[3] = 5; OrderA[4] = 3; OrderA[5] = 4;
     A.permute(OrderA);
    }
    if ((col == this->N[1]-1) && (row == this->N[0]-1))
    {
     IndicesA[0] = 2;
     A.contract(IndicesA, P2, IndicesP);
     P2 = P2C;
     A.contract(IndicesA, P2, IndicesP);
     OrderA[0] = 0; OrderA[1] = 1; OrderA[2] = 4; OrderA[3] = 5; OrderA[4] = 2; OrderA[5] = 3;
     A.permute(OrderA);
    }
    PEPO0.set(row, col, A);
   }
  }
 }
// periodic boundary conditions:
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented for periodic BC: " <<
          "template<class T> void IsingModel2D<T>::" <<
          "getPEPO(PEPO<T>& PEPO0) const." << endl;
  exit(1);
 }
}

template<class T> void IsingModel2D<T>::getMatrix(Matrix<T>& Matrix0) const
{
#ifdef DEBUG
 unsigned int dim = pow(2.0, int(this->N[0]*this->N[1]));
 if ((this->N[0] == 0) || (Matrix0.getDim0() != dim) || (Matrix0.getDim1() != dim))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void IsingModel2D<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const: " <<
          "((this->N[0] == 0) || (Matrix0.getDim0() != 2^{Nrows*Ncols}) || " <<
           "(Matrix0.getDim1() != 2^{Nrows*Ncols}))." << endl;
  exit(1);
 }
#endif
 Matrix0.fillZeroes();
// Interactions:
// - Identity:
 Matrix<T> Identity(2, 2);
 Identity.fillZeroes();
 for (int i = 0; i < 2; i++)
  Identity(i, i) = 1.0;
// - Pauli matrix sx:
 Matrix<T> sx(2, 2);
 sx.fillZeroes();
 sx(1, 0) = 1.0;
 sx(0, 1) = 1.0;
// - Pauli matrix sz:
 Matrix<T> sz(2, 2);
 sz.fillZeroes();
 sz(0, 0) = 1.0;
 sz(1, 1) = -1.0;
// - S := sz \otimes sz:
 Matrix<T> S;
 sz.multiplyDirectProduct(sz, S);
// parameters:
 T B = this->Parameters[0], J = this->Parameters[1];
 Matrix<T> X, Y;
 unsigned int position;
// open boundary conditions:
 if (this->BC == "open")
 {
// Hamiltonian sum:
// - B \sum_{l} sx_{l}:
  for (int i = 0; i < this->N[0]*this->N[1]; i++)
  {
   X = sx;
   X.multiply(B);
   for (int j = i+1; j < this->N[0]*this->N[1]; j++)
   {
    X.multiplyDirectProduct(Identity, Y);
    X = Y;
   }
   for (int j = i-1; j >= 0; j--)
   {
    Identity.multiplyDirectProduct(X, Y);
    X = Y;
   }
   Matrix0.add(X);
  }
// - J \sum_{<l,m>} sz_{l} sz_{m}, vertical terms:
  for (int col = 0; col < this->N[1]; col++)
  {
   for (int row = 0; row < this->N[0]-1; row++)
   {
    position = row + this->N[0]*col;
    X = S;
    X.multiply(J);
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
// - J \sum_{<l,m>} sz_{l} sz_{m}, horizontal terms:
  for (int col = 0; col < this->N[1]-1; col++)
  {
   for (int row = 0; row < this->N[0]; row++)
   {
    position = row + this->N[0]*col;
    X = sz;
    X.multiply(J);
    for (int j = 1; j < this->N[0]; j++)
    {
     X.multiplyDirectProduct(Identity, Y);
     X = Y;
    }
    X.multiplyDirectProduct(sz, Y);
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
          "template<class T> void IsingModel2D<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const." << endl;
  exit(1);
 }
 Matrix0.setType("hermitian");
}

template<class T> void IsingModel2D<T>::getHMPOs(vector< MPO<T> >& HMPOs) const
{
 HMPOs = vector< MPO<T> >(4);
 T element, B = this->Parameters[0], J = this->Parameters[1];
 vector<unsigned int> Index4(4), Shape4(4);
 Tensor<T> TensorA;
// 1. HMPOs[0] := J*s^{z}s^{z}:
 HMPOs[0] = MPO<T>("open", 2, 2, 1);
 Shape4[0] = 1; Shape4[1] = 1; Shape4[2] = 2; Shape4[3] = 2;
 TensorA = Tensor<T>(Shape4);
 TensorA.fillZeroes();
 Index4[0] = 0; Index4[1] = 0;
 element = J;
 Index4[2] = 0; Index4[3] = 0;
 TensorA.set(Index4, element);
 element = -J;
 Index4[2] = 1; Index4[3] = 1;
 TensorA.set(Index4, element);
 HMPOs[0].set(0, TensorA);
 TensorA.fillZeroes();
 element = 1.0;
 Index4[2] = 0; Index4[3] = 0;
 TensorA.set(Index4, element);
 element = -1.0;
 Index4[2] = 1; Index4[3] = 1;
 TensorA.set(Index4, element);
 HMPOs[0].set(1, TensorA);
// 2. HMPOs[1] := J*s^{z}s^{z} + B*s^{x}*one:
 HMPOs[1] = MPO<T>("open", 2, 2, 2);
 Shape4[0] = 1; Shape4[1] = 2; Shape4[2] = 2; Shape4[3] = 2;
 TensorA = Tensor<T>(Shape4);
 TensorA.fillZeroes();
 Index4[0] = 0;
 element = J;
 Index4[1] = 0; Index4[2] = 0; Index4[3] = 0;
 TensorA.set(Index4, element);
 element = -J;
 Index4[2] = 1; Index4[3] = 1;
 TensorA.set(Index4, element);
 element = B;
 Index4[1] = 1; Index4[2] = 1; Index4[3] = 0;
 TensorA.set(Index4, element);
 Index4[2] = 0; Index4[3] = 1;
 TensorA.set(Index4, element);
 HMPOs[1].set(0, TensorA);
 Shape4[0] = 2; Shape4[1] = 1; Shape4[2] = 2; Shape4[3] = 2;
 TensorA = Tensor<T>(Shape4);
 TensorA.fillZeroes();
 Index4[1] = 0;
 element = 1.0;
 Index4[0] = 0; Index4[2] = 0; Index4[3] = 0;
 TensorA.set(Index4, element);
 element = -1.0;
 Index4[2] = 1; Index4[3] = 1;
 TensorA.set(Index4, element);
 element = 1.0;
 Index4[0] = 1; Index4[2] = 0; Index4[3] = 0;
 TensorA.set(Index4, element);
 Index4[2] = 1; Index4[3] = 1;
 TensorA.set(Index4, element);
 HMPOs[1].set(1, TensorA);
// 3. HMPOs[2] := J*s^{z}s^{z} + B*one*s^{x}:
 HMPOs[2] = MPO<T>("open", 2, 2, 2);
 Shape4[0] = 1; Shape4[1] = 2; Shape4[2] = 2; Shape4[3] = 2;
 TensorA = Tensor<T>(Shape4);
 TensorA.fillZeroes();
 Index4[0] = 0;
 element = J;
 Index4[1] = 0; Index4[2] = 0; Index4[3] = 0;
 TensorA.set(Index4, element);
 element = -J;
 Index4[2] = 1; Index4[3] = 1;
 TensorA.set(Index4, element);
 element = B;
 Index4[1] = 1; Index4[2] = 0; Index4[3] = 0;
 TensorA.set(Index4, element);
 Index4[2] = 1; Index4[3] = 1;
 TensorA.set(Index4, element);
 HMPOs[2].set(0, TensorA);
 Shape4[0] = 2; Shape4[1] = 1; Shape4[2] = 2; Shape4[3] = 2;
 TensorA = Tensor<T>(Shape4);
 TensorA.fillZeroes();
 Index4[1] = 0;
 element = 1.0;
 Index4[0] = 0; Index4[2] = 0; Index4[3] = 0;
 TensorA.set(Index4, element);
 element = -1.0;
 Index4[2] = 1; Index4[3] = 1;
 TensorA.set(Index4, element);
 element = 1.0;
 Index4[0] = 1; Index4[2] = 1; Index4[3] = 0;
 TensorA.set(Index4, element);
 Index4[2] = 0; Index4[3] = 1;
 TensorA.set(Index4, element);
 HMPOs[2].set(1, TensorA);
// 4. HMPOs[3] := J*s^{z}s^{z} + B*s^{x}*one + B*one*s^{x}:
 HMPOs[3] = MPO<T>("open", 2, 2, 3);
 Shape4[0] = 1; Shape4[1] = 3; Shape4[2] = 2; Shape4[3] = 2;
 TensorA = Tensor<T>(Shape4);
 TensorA.fillZeroes();
 Index4[0] = 0;
 element = J;
 Index4[1] = 0; Index4[2] = 0; Index4[3] = 0;
 TensorA.set(Index4, element);
 element = -J;
 Index4[2] = 1; Index4[3] = 1;
 TensorA.set(Index4, element);
 element = B;
 Index4[1] = 1; Index4[2] = 1; Index4[3] = 0;
 TensorA.set(Index4, element);
 Index4[2] = 0; Index4[3] = 1;
 TensorA.set(Index4, element);
 Index4[1] = 2; Index4[2] = 0; Index4[3] = 0;
 TensorA.set(Index4, element);
 Index4[2] = 1; Index4[3] = 1;
 TensorA.set(Index4, element);
 HMPOs[3].set(0, TensorA);
 Shape4[0] = 3; Shape4[1] = 1; Shape4[2] = 2; Shape4[3] = 2;
 TensorA = Tensor<T>(Shape4);
 TensorA.fillZeroes();
 Index4[1] = 0;
 element = 1.0;
 Index4[0] = 0; Index4[2] = 0; Index4[3] = 0;
 TensorA.set(Index4, element);
 element = -1.0;
 Index4[2] = 1; Index4[3] = 1;
 TensorA.set(Index4, element);
 element = 1.0;
 Index4[0] = 1; Index4[2] = 0; Index4[3] = 0;
 TensorA.set(Index4, element);
 Index4[2] = 1; Index4[3] = 1;
 TensorA.set(Index4, element);
 Index4[0] = 2; Index4[2] = 1; Index4[3] = 0;
 TensorA.set(Index4, element);
 Index4[2] = 0; Index4[3] = 1;
 TensorA.set(Index4, element);
 HMPOs[3].set(1, TensorA);
}

template<class T> void IsingModel2D<T>::getTEMPOs(const string& RealImaginaryTE, double timeStep,
                                                  vector< MPO<T> >& TEMPOs)
{
#ifdef DEBUG
 if ((RealImaginaryTE != "real") && (RealImaginaryTE != "imaginary"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void IsingModel2D<T>::" <<
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
 T B = this->Parameters[0], J = this->Parameters[1];
// TensorX = ExpX = exp(-delta*HX) = exp(-delta*B*SigmaX):
 Matrix<T> HX(2, 2), ExpX(2, 2);
 HX.fillZeroes(); HX(1, 0) = B; HX(0, 1) = B;
 HX.setType("hermitian");
 HX.exponentialize(-delta, ExpX);
 vector<unsigned int> ShapeX(2), IndexX(2);
 ShapeX[0] = 2; ShapeX[1] = 2;
 Tensor<T> TensorX(ShapeX);
 T element;
 for (int j = 0; j < 2; j++)
 {
  IndexX[0] = j;
  for (int i = 0; i < 2; i++)
  {
   IndexX[1] = i;
   element = ExpX(i, j);
   TensorX.set(IndexX, element);
  }
 }
// TensorZLeft*TensorZRight = ExpZ = exp(-delta*HZ) = exp(-delta*J*SigmaZ*SigmaZ):
 Matrix<T> HZ0(2, 2), HZ1(2, 2), HZ(4, 4), ExpZ(4, 4);
 HZ0.fillZeroes(); HZ0(0, 0) = J; HZ0(1, 1) = -J;
 HZ1.fillZeroes(); HZ1(0, 0) = 1.0; HZ1(1, 1) = -1.0;
 HZ0.multiplyDirectProduct(HZ1, HZ);
 vector<unsigned int> ShapeLeft(4), ShapeRight(4);
 ShapeLeft[0] = 1; ShapeLeft[1] = 4; ShapeLeft[2] = 2; ShapeLeft[3] = 2;
 ShapeRight[0] = 4; ShapeRight[1] = 1; ShapeRight[2] = 2; ShapeRight[3] = 2;
 Tensor<T> TensorLeft(ShapeLeft), TensorRight(ShapeRight);
 HZ.setType("hermitian");
 HZ.twoBodyExponentSVD(delta, 2, 2, TensorLeft, TensorRight);
// - cut the virtual bond range from 4 to 2:
 vector<unsigned int> ShapeZLeft(4), ShapeZRight(4), Index4(4);
 ShapeZLeft[0] = 1; ShapeZLeft[1] = 2; ShapeZLeft[2] = 2; ShapeZLeft[3] = 2;
 ShapeZRight[0] = 2; ShapeZRight[1] = 1; ShapeZRight[2] = 2; ShapeZRight[3] = 2;
 Tensor<T> TensorZLeft(ShapeZLeft), TensorZRight(ShapeZRight);
 for (int i3 = 0; i3 < ShapeZLeft[3]; i3++)
 {
  Index4[3] = i3;
  for (int i2 = 0; i2 < ShapeZLeft[2]; i2++)
  {
   Index4[2] = i2;
   for (int i1 = 0; i1 < ShapeZLeft[1]; i1++)
   {
    Index4[1] = i1;
    for (int i0 = 0; i0 < ShapeZLeft[0]; i0++)
    {
     Index4[0] = i0;
     element = TensorLeft.get(Index4);
     TensorZLeft.set(Index4, element);
    }
   }
  }
 }
 for (int i3 = 0; i3 < ShapeZRight[3]; i3++)
 {
  Index4[3] = i3;
  for (int i2 = 0; i2 < ShapeZRight[2]; i2++)
  {
   Index4[2] = i2;
   for (int i1 = 0; i1 < ShapeZRight[1]; i1++)
   {
    Index4[1] = i1;
    for (int i0 = 0; i0 < ShapeZRight[0]; i0++)
    {
     Index4[0] = i0;
     element = TensorRight.get(Index4);
     TensorZRight.set(Index4, element);
    }
   }
  }
 }
// TEMPOs[0] = TEMPO0 = exp(-delta*H_{Z,Z}):
 MPO<T> TEMPO0("open", 2, 2, 2);
 TEMPO0.set(0, TensorZLeft);
 TEMPO0.set(1, TensorZRight);
 TEMPOs[0] = TEMPO0;
// TEMPOs[1] = TEMPO1 = exp(-delta*H_{XZ,Z}):
 Tensor<T> TensorXZLeft(TensorZLeft), Tensor1(TensorX);
 vector<unsigned int> Indices0(1), Indices1(1);
 Indices0[0] = 3; Indices1[0] = 0;
 TensorXZLeft.contract(Indices0, Tensor1, Indices1);
 MPO<T> TEMPO1("open", 2, 2, 2);
 TEMPO1.set(0, TensorXZLeft);
 TEMPO1.set(1, TensorZRight);
 TEMPOs[1] = TEMPO1;
// TEMPOs[2] = TEMPO2 = exp(-delta*H_{Z,XZ}):
 Tensor<T> TensorXZRight(TensorZRight);
 TensorXZRight.contract(Indices0, TensorX, Indices1);
 MPO<T> TEMPO2("open", 2, 2, 2);
 TEMPO2.set(0, TensorZLeft);
 TEMPO2.set(1, TensorXZRight);
 TEMPOs[2] = TEMPO2;
// TEMPOs[3] = TEMPO3 = exp(-delta*H_{XZ,XZ}):
 MPO<T> TEMPO3("open", 2, 2, 2);
 TEMPO3.set(0, TensorXZLeft);
 TEMPO3.set(1, TensorXZRight);
 TEMPOs[3] = TEMPO3;
 if ((RealImaginaryTE == "real") && (this->timeDependent))
  this->setTime(this->time-timeStep);
}

template<class T> void IsingModel2D<T>::getTEPEPOs(const string& RealImaginaryTE, double timeStep,
                                                   vector< PEPO<T> >& TEPEPOs)
{
#ifdef DEBUG
 if ((RealImaginaryTE != "real") && (RealImaginaryTE != "imaginary"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void IsingModel2D<T>::" <<
          "getTEPEPOs(const string& RealImaginaryTE, double timeStep, " <<
                     "vector< PEPO<T> >& TEPEPOs): " <<
          "((RealImaginaryTE != real) && (RealImaginaryTE != imaginary))." << endl;
  exit(1);
 }
#endif
 unsigned int col, Ncols = this->N[1], Nrows = this->N[0], row;
 vector<unsigned int> Index6(6), Shape6(6);
 Tensor<T> OneTensor, TETensor0, TETensor1;
 vector< MPO<T> > TEMPOs;
// 0. initialize:
 TEPEPOs = vector< PEPO<T> >(4);
 for (int i = 0; i < 4; i++)
  TEPEPOs[i] = PEPO<T>(this->BC, Nrows, Ncols, this->d, 2);
// 1. get TEMPOs:
 this->getTEMPOs(RealImaginaryTE, timeStep, TEMPOs);
// 2. set TEPEPOs:
// 2.1 set vertical-even Trotter gates in TEPEPOs[0]:
 TEMPOs[0].get(0, TETensor0);
 Shape6[0] = 1; Shape6[1] = 1; Shape6[2] = 1; Shape6[3] = 2; Shape6[4] = 2; Shape6[5] = 2;
 TETensor0.reshape(Shape6);
 TEMPOs[0].get(1, TETensor1);
 Shape6[0] = 1; Shape6[1] = 2; Shape6[2] = 1; Shape6[3] = 1; Shape6[4] = 2; Shape6[5] = 2;
 TETensor1.reshape(Shape6);
 for (col = 0; col < Ncols; col++)
 {
  for (row = 0; row < Nrows-1; row+=2)
  {
   TEPEPOs[0].set(row, col, TETensor0);
   TEPEPOs[0].set(row+1, col, TETensor1);
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
   TEPEPOs[1].set(row, col, TETensor0);
   TEPEPOs[1].set(row+1, col, TETensor1);
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
// 2.3 - column 0:
 col = 0;
 TEMPOs[1].get(0, TETensor0);
 Shape6[0] = 1; Shape6[1] = 1; Shape6[2] = 2; Shape6[3] = 1; Shape6[4] = 2; Shape6[5] = 2;
 TETensor0.reshape(Shape6);
 TEMPOs[1].get(1, TETensor1);
 Shape6[0] = 2; Shape6[1] = 1; Shape6[2] = 1; Shape6[3] = 1; Shape6[4] = 2; Shape6[5] = 2;
 TETensor1.reshape(Shape6);
 for (row = 0; row < Nrows; row++)
 {
  TEPEPOs[2].set(row, col, TETensor0);
  TEPEPOs[2].set(row, col+1, TETensor1);
 }
// 2.3 - column > 0:
 TEMPOs[0].get(0, TETensor0);
 Shape6[0] = 1; Shape6[1] = 1; Shape6[2] = 2; Shape6[3] = 1; Shape6[4] = 2; Shape6[5] = 2;
 TETensor0.reshape(Shape6);
 TEMPOs[0].get(1, TETensor1);
 Shape6[0] = 2; Shape6[1] = 1; Shape6[2] = 1; Shape6[3] = 1; Shape6[4] = 2; Shape6[5] = 2;
 TETensor1.reshape(Shape6);
 for (col = 2; col < Ncols-2; col+=2)
 {
  for (row = 0; row < Nrows; row++)
  {
   TEPEPOs[2].set(row, col, TETensor0);
   TEPEPOs[2].set(row, col+1, TETensor1);
  }
 }
 if (Ncols%2 == 0)
 {
// 2.3 - column Ncols-2:
  col = Ncols-2;
  TEMPOs[2].get(0, TETensor0);
  Shape6[0] = 1; Shape6[1] = 1; Shape6[2] = 2; Shape6[3] = 1; Shape6[4] = 2; Shape6[5] = 2;
  TETensor0.reshape(Shape6);
  TEMPOs[2].get(1, TETensor1);
  Shape6[0] = 2; Shape6[1] = 1; Shape6[2] = 1; Shape6[3] = 1; Shape6[4] = 2; Shape6[5] = 2;
  TETensor1.reshape(Shape6);
  for (row = 0; row < Nrows; row++)
  {
   TEPEPOs[2].set(row, col, TETensor0);
   TEPEPOs[2].set(row, col+1, TETensor1);
  }
 }
 else if (Ncols%2 == 1)
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
 TEMPOs[3].get(0, TETensor0);
 Shape6[0] = 1; Shape6[1] = 1; Shape6[2] = 2; Shape6[3] = 1; Shape6[4] = 2; Shape6[5] = 2;
 TETensor0.reshape(Shape6);
 TEMPOs[3].get(1, TETensor1);
 Shape6[0] = 2; Shape6[1] = 1; Shape6[2] = 1; Shape6[3] = 1; Shape6[4] = 2; Shape6[5] = 2;
 TETensor1.reshape(Shape6);
 for (col = 1; col < Ncols-1; col+=2)
 {
  for (row = 0; row < Nrows; row++)
  {
   TEPEPOs[3].set(row, col, TETensor0);
   TEPEPOs[3].set(row, col+1, TETensor1);
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
