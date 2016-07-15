/// Template class BoseHubbardModel implements Bose-Hubbard model.
/** The template class BoseHubbardModel implements the Bose-Hubbard model:
    H = - \sum_{l} t_{l} (a_{l+1}^{+} a_{l} + a_{l}^{+} a_{l+1})
        + \sum_{l} 0.5*U_{l} n_{l} (n_{l} - 1)
        + \sum_{l} V_{l} n_{l}
        - mu \sum_{l} n_{l}   .
    We have site-dependent parameters tunneling t_{l}, on-site interaction U_{l} and offset V_{l}.
    \param Representation string, the representation, must be "interactions", "MPO" or "matrix"
    \param BC string, the boundary conditions, must be "open" or "periodic"
    \param N unsigned int, the number of spins, i.e. lattice sites
    \param d vector<unsigned int>, the physical dimensions of the lattice sites,
             i.e. d[l]-1 is the maximal number of bosons allowed on lattice site l
    \param Parameters vector<T>, the parameters {t_{0}, ..., t_{N-2}, U_{0}, ..., U_{N-1}, V_{0}, ..., V_{N-1}, mu}
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

template<class T> class BoseHubbardModel: public Hamiltonian<T>
{
 typedef T (*PointerToFunction)(double time);

 public:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  BoseHubbardModel();

/// Constructor for time-independent BoseHubbardModel with specific BC, N, d and Parameters.
/** This constructor initializes a time-independent BoseHubbardModel.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N0 input: unsigned int, the number of lattice sites
    \param d0 input: unsigned int, d0-1 is the maximal number of bosons per lattice site
    \param Parameters0 input: const vector<T>&, the parameters, must have Parameters0.size()==3*N0 */
  BoseHubbardModel(const string& BC0, unsigned int N0, unsigned int d0, const vector<T>& Parameters0);

/// Constructor for time-dependent BoseHubbardModel with specific BC, N, d, TimeFunctions and time.
/** This constructor initializes a time-dependent BoseHubbardModel.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N0 input: unsigned int, the number of lattice sites
    \param d0 input: unsigned int, d0-1 is the maximal number of bosons per lattice site
    \param TimeFunctions0 input: const vector<PointerToFunction>&, the time functions, must have
                                 TimeFunctions.size()==3*N0
    \param time0 input: double, the time */
  BoseHubbardModel(const string& BC0, unsigned int N0, unsigned int d0,
                   const vector<PointerToFunction>& TimeFunctions0, double time0);

/// Standard copy constructor.
/** The standard copy constructor copies the input BoseHubbardModel into this.
    \param BoseHubbardModel0 input: const BoseHubbardModel<T>&, to be copied into this
    \sa BoseHubbardModel<T>& operator=(const BoseHubbardModel<T>& BoseHubbardModel0) */
  BoseHubbardModel(const BoseHubbardModel<T>& BoseHubbardModel0);

/// Standard destructor.
/** The standard destructor deletes the elements of BoseHubbardModel. */
  ~BoseHubbardModel();

/// Assigns BoseHubbardModel to this.
/** The operator= allows to assign BoseHubbardModel0 to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side BoseHubbardModel0.
    \param BoseHubbardModel0 input: const BoseHubbardModel<T>&, to be copied into this
    \return BoseHubbardModel<T>&, a reference to the new this
    \sa BoseHubbardModel(const BoseHubbardModel<T>& BoseHubbardModel0) */
  BoseHubbardModel<T>& operator=(const BoseHubbardModel<T>& BoseHubbardModel0);

/// Returns interactions representation of BoseHubbardModel.
/** This function returns the interactions of this BoseHubbardModel. The interactions
    are implemented as a vector of positions and a vector of matrices where each vector entry
    corresponds to one term of the sum making up the Hamiltonian.
    \param Positions output: vector< vector<unsigned int> >&, the positions of the matrices
    \param Matrices output: vector< vector< Matrix<T> > >&, the matrices */
  void getInteractions(vector< vector<unsigned int> >& Positions,
                       vector< vector< Matrix<T> > >& Matrices) const;

/// Returns MPO representation of BoseHubbardModel.
/** This function returns this BoseHubbardModel as a MPO.
    \param MPO0 output: MPO<T>&, the MPO representing this BoseHubbardModel */
  void getMPO(MPO<T>& MPO0) const;

/// Returns matrix representation of BoseHubbardModel.
/** This function returns this BoseHubbardModel as a matrix.
    \param Matrix0 output: Matrix<T>&, the matrix representing this BoseHubbardModel,
                           must fulfill ((Matrix0.getDim0()==d^{N}) && (Matrix0.getDim1()==d^{N})),
                           will be of type "hermitian" */
  void getMatrix(Matrix<T>& Matrix0) const;

/// Returns local particle number operator as MPO.
/** This function returns the local particle number operator n_{x} at position x as a MPO.
    \param x input: unsigned int, the position
    \param nMPO output: MPO<T>&, the local particle number operator as a MPO
    \sa void getNMPO(MPO<T>& NMPO) const */
  void getnMPO(unsigned int x, MPO<T>& nMPO) const;

/// Returns total particle number operator as MPO.
/** This function returns the total particle number operator N = \sum_{l} n_{l} as a MPO.
    \param NMPO output: MPO<T>&, the total particle number operator as a MPO
    \sa void getnMPO(unsigned int x, MPO<T>& nMPO) const */
  void getNMPO(MPO<T>& NMPO) const;

 private:

/// Returns two body Hamiltonian for Trotter terms.
/** This function computes the hermitian part H of the exponent for the Trotter terms exp(-delta*H) for
    the even-odd Trotter decomposition and returns it as a Matrix TwoBodyHamiltonian. This function has
    to be implemented for Hamiltonian<T>::getTEMPOs to work.
    \param position input: unsigned int, the left position of the two body Hamiltonian in the
                           Hamiltonian sum, must be out of {0, 1, ..., N-2}
    \param TwoBodyHamiltonian output: Matrix<T>&, the two body Hamiltonian, must have the correct
                                      form */
  void getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian) const;
};

template<class T> BoseHubbardModel<T>::BoseHubbardModel()
{
 this->N = 0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> BoseHubbardModel<T>::BoseHubbardModel(const string& BC0, unsigned int N0, unsigned int d0,
                                                        const vector<T>& Parameters0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (N0 == 0) || (Parameters0.size() != 3*N0))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> BoseHubbardModel<T>::" <<
          "BoseHubbardModel(const string& BC0, unsigned int N0, unsigned int d0, " <<
                           "const vector<T>& Parameters0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (N0 == 0) || (Parameters0.size() != 3*N0))." << endl;
  exit(1);
 }
#endif
 this->Representation = "MPO";
 this->BC = BC0;
 this->N = N0;
 this->d = vector<unsigned int>(this->N);
 for (int i = 0; i < this->N; i++)
 {
  this->d[i] = d0;
 }
 this->Parameters = Parameters0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> BoseHubbardModel<T>::BoseHubbardModel(const string& BC0, unsigned int N0, unsigned int d0,
                                              const vector<PointerToFunction>& TimeFunctions0, double time0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (N0 == 0) || (TimeFunctions0.size() != 3*N0))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> BoseHubbardModel<T>::" <<
          "BoseHubbardModel(const string& BC0, unsigned int N0, unsigned int d0, " <<
                           "const vector<PointerToFunction>& TimeFunctions0, double time0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (N0 == 0) || (TimeFunctions0.size() != 3*N0))." << endl;
  exit(1);
 }
#endif
 this->Representation = "MPO";
 this->BC = BC0;
 this->N = N0;
 this->d = vector<unsigned int>(this->N);
 for (int i = 0; i < this->N; i++)
 {
  this->d[i] = d0;
 }
 this->Parameters = vector<T>(3*this->N);
 this->timeDependent = true;
 this->TimeFunctions = TimeFunctions0;
 this->time = time0;
 for (int i = 0; i < this->Parameters.size(); i++)
 {
  this->Parameters[i] = this->TimeFunctions[i](this->time);
 }
}

template<class T> BoseHubbardModel<T>::BoseHubbardModel(const BoseHubbardModel<T>& BoseHubbardModel0)
{
 this->Representation = BoseHubbardModel0.Representation;
 this->BC = BoseHubbardModel0.BC;
 this->N = BoseHubbardModel0.N;
 this->d = BoseHubbardModel0.d;
 this->Parameters = BoseHubbardModel0.Parameters;
 this->timeDependent = BoseHubbardModel0.timeDependent;
 this->TimeFunctions = BoseHubbardModel0.TimeFunctions;
 this->time = BoseHubbardModel0.time;
}

template<class T> BoseHubbardModel<T>::~BoseHubbardModel() {}

template<class T> BoseHubbardModel<T>& BoseHubbardModel<T>::operator=(const BoseHubbardModel<T>& BoseHubbardModel0)
{
 if (this != &BoseHubbardModel0)
 {
  this->Representation = BoseHubbardModel0.Representation;
  this->BC = BoseHubbardModel0.BC;
  this->N = BoseHubbardModel0.N;
  this->d = BoseHubbardModel0.d;
  this->Parameters = BoseHubbardModel0.Parameters;
  this->timeDependent = BoseHubbardModel0.timeDependent;
  this->TimeFunctions = BoseHubbardModel0.TimeFunctions;
  this->time = BoseHubbardModel0.time;
 }
 return *this;
}

template<class T> void BoseHubbardModel<T>::getInteractions(vector< vector<unsigned int> >& Positions,
                                                            vector< vector< Matrix<T> > >& Matrices) const
{
 cerr << "The following function is not implemented yet: " <<
         "template<class T> void BoseHubbardModel<T>::" <<
         "getInteractions(vector< vector<unsigned int> >& Positions, " <<
                         "vector< vector< Matrix<T> > >& Matrices) const." << endl;
 exit(1);
}

template<class T> void BoseHubbardModel<T>::getMPO(MPO<T>& MPO0) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void BoseHubbardModel<T>::" <<
          "getMPO(MPO<T>& MPO0) const: " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
 if (this->BC == "open")
 {
// the MPO operators O^{i}:
  vector<unsigned int> Shape(3);
  Shape[0] = 5; Shape[1] = this->d[0]; Shape[2] = this->d[0];
  Tensor<T> O(Shape);
  O.fillZeroes();
  vector<unsigned int> Index(3);
// O^{0} = identity:
  Index[0] = 0;
  T element = 1.0;
  for (int i = 0; i < this->d[0]; i++)
  {
   Index[1] = i; Index[2] = i;
   O.set(Index, element);
  }
// O^{1} = N:
  Index[0] = 1;
  for (int i = 0; i < this->d[0]; i++)
  {
   element = T(i);
   Index[1] = i; Index[2] = i;
   O.set(Index, element);
  }
// O^{2} = 0.5*N(N-1):
  Index[0] = 2;
  for (int i = 0; i < this->d[0]; i++)
  {
   element = 0.5*T(i*(i-1));
   Index[1] = i; Index[2] = i;
   O.set(Index, element);
  }
// O^{3} = A^{+}:
  Index[0] = 3;
  for (int i = 1; i < this->d[0]; i++)
  {
   element = sqrt(T(i));
   Index[1] = i-1; Index[2] = i;
   O.set(Index, element);
  }
// O^{4} = A^{-}:
  Index[0] = 4;
  for (int i = 1; i < this->d[0]; i++)
  {
   element = sqrt(T(i));
   Index[1] = i; Index[2] = i-1;
   O.set(Index, element);
  }
  unsigned int D = 4;
  MPO0 = MPO<T>(this->BC, this->N, this->d, D);
// the parameters:
  vector<T> t(this->N-1);
  for (int i = 0; i < this->N-1; i++)
  {
   t[i] = this->Parameters[i];
  }
  vector<T> U(this->N);
  for (int i = 0; i < this->N; i++)
  {
   U[i] = this->Parameters[this->N-1+i];
  }
  vector<T> V(this->N);
  for (int i = 0; i < this->N; i++)
  {
   V[i] = this->Parameters[2*this->N-1+i];
  }
  T mu = this->Parameters[3*this->N-1];
// A for position 0:
  Shape[0] = 1; Shape[1] = D; Shape[2] = 5;
  Tensor<T> A(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  element = V[0] - mu;
  Index[0] = 0; Index[1] = 3; Index[2] = 1;
  A.set(Index, element);
  element = U[0];
  Index[0] = 0; Index[1] = 3; Index[2] = 2;
  A.set(Index, element);
  element = -t[0];
  Index[0] = 0; Index[1] = 1; Index[2] = 3;
  A.set(Index, element);
  Index[0] = 0; Index[1] = 2; Index[2] = 4;
  A.set(Index, element);
  Tensor<T> OC(O);
  vector<unsigned int> IndexA(1), IndexOC(1);
  IndexA[0] = 2; IndexOC[0] = 0;
  A.contract(IndexA, OC, IndexOC);
  unsigned int position = 0;
  MPO0.set(position, A);
// A for position 1 <= l <= N-2:
  Shape[0] = D; Shape[1] = D; Shape[2] = 5;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 3; Index[1] = 3; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 2; Index[1] = 3; Index[2] = 3;
  A.set(Index, element);
  Index[0] = 1; Index[1] = 3; Index[2] = 4;
  A.set(Index, element);
  Tensor<T> AC(A);
  for (position = 1; position < this->N-1; position++)
  {
   A = AC;
   element = V[position] - mu;
   Index[0] = 0; Index[1] = 3; Index[2] = 1;
   A.set(Index, element);
   element = U[position];
   Index[0] = 0; Index[1] = 3; Index[2] = 2;
   A.set(Index, element);
   element = -t[position];
   Index[0] = 0; Index[1] = 1; Index[2] = 3;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 2; Index[2] = 4;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   MPO0.set(position, A);
  }
// A for position N-1:
  Shape[0] = D; Shape[1] = 1; Shape[2] = 5;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 3; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  element = V[this->N-1] - mu;
  Index[0] = 0; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  element = U[this->N-1];
  Index[0] = 0; Index[1] = 0; Index[2] = 2;
  A.set(Index, element);
  element = 1.0;
  Index[0] = 2; Index[1] = 0; Index[2] = 3;
  A.set(Index, element);
  Index[0] = 1; Index[1] = 0; Index[2] = 4;
  A.set(Index, element);
  OC = O;
  A.contract(IndexA, OC, IndexOC);
  position = this->N-1;
  MPO0.set(position, A);
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void BoseHubbardModel<T>::" <<
          "getMPO(MPO<T>& MPO0) const." << endl;
  exit(1);
 }
}

template<class T> void BoseHubbardModel<T>::getMatrix(Matrix<T>& Matrix0) const
{
#ifdef DEBUG
 unsigned int dim = 1;
 for (int i = 0; i < this->N; i++)
 {
  dim *= this->d[0];
 }
 if ((Matrix0.getDim0() != dim) || (Matrix0.getDim1() != dim))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void BoseHubbardModel<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const: " <<
          "((Matrix0.getDim0() != this->d[0]^{N}) || (Matrix0.getDim1() != this->d[0]^{N}))." << endl;
  exit(1);
 }
#endif
 Matrix0.fillZeroes();
// the interactions:
// - the identity:
 Matrix<T> Identity(this->d[0], this->d[0]);
 Identity.fillZeroes();
 T element = 1.0;
 for (int i = 0; i < this->d[0]; i++)
 {
  Identity(i, i) = element;
 }
// - particle number N =: N0:
 Matrix<T> N0(this->d[0], this->d[0]);
 N0.fillZeroes();
 for (int i = 0; i < this->d[0]; i++)
 {
  element = T(i);
  N0(i, i) = element;
 }
// - interaction particle number 0.5*N(N-1) =: N1:
 Matrix<T> N1(this->d[0], this->d[0]);
 N1.fillZeroes();
 for (int i = 0; i < this->d[0]; i++)
 {
  element = 0.5*T(i*(i-1));
  N1(i, i) = element;
 }
// - A^{+} =: A0:
 Matrix<T> A0(this->d[0], this->d[0]);
 A0.fillZeroes();
 for (int i = 1; i < this->d[0]; i++)
 {
  element = sqrt(T(i));
  A0(i, i-1) = element;
 }
// - A^{-} =: A1:
 Matrix<T> A1(this->d[0], this->d[0]);
 A1.fillZeroes();
 for (int i = 1; i < this->d[0]; i++)
 {
  element = sqrt(T(i));
  A1(i-1, i) = element;
 }
// - A := A^{+} \otimes A^{-} + A^{-} \otimes A^{+}:
 Matrix<T> A(this->d[0]*this->d[0], this->d[0]*this->d[0]), X, Y;
 A.fillZeroes();
 A0.multiplyDirectProduct(A1, X);
 A.add(X);
 A1.multiplyDirectProduct(A0, X);
 A.add(X);
// the parameters:
 vector<T> t(this->N-1);
 for (int i = 0; i < this->N-1; i++)
 {
  t[i] = this->Parameters[i];
 }
 vector<T> U(this->N);
 for (int i = 0; i < this->N; i++)
 {
  U[i] = this->Parameters[this->N-1+i];
 }
 vector<T> V(this->N);
 for (int i = 0; i < this->N; i++)
 {
  V[i] = this->Parameters[2*this->N-1+i];
 }
 T mu = this->Parameters[3*this->N-1];
 if (this->BC == "open")
 {
// the Hamiltonian sum:
// - \sum_{l} -t_{l} (a_{l+1}^{+} a_{l} + a_{l}^{+} a_{l+1}):
  for (int i = 0; i < this->N-1; i++)
  {
   X = A;
   for (int j = i+2; j < this->N; j++)
   {
    X.multiplyDirectProduct(Identity, Y);
    X = Y;
   }
   for (int j = i-1; j >= 0; j--)
   {
    Identity.multiplyDirectProduct(X, Y);
    X = Y;
   }
   element = -t[i];
   X.multiply(element);
   Matrix0.add(X);
  }
// - \sum_{l} 0.5*U_{l} n_{l} (n_{l} - 1):
  for (int i = 0; i < this->N; i++)
  {
   X = N1;
   for (int j = i+1; j < this->N; j++)
   {
    X.multiplyDirectProduct(Identity, Y);
    X = Y;
   }
   for (int j = i-1; j >= 0; j--)
   {
    Identity.multiplyDirectProduct(X, Y);
    X = Y;
   }
   element = U[i];
   X.multiply(element);
   Matrix0.add(X);
  }
// - \sum_{l} (V_{l}-mu) n_{l}:
  for (int i = 0; i < this->N; i++)
  {
   X = N0;
   for (int j = i+1; j < this->N; j++)
   {
    X.multiplyDirectProduct(Identity, Y);
    X = Y;
   }
   for (int j = i-1; j >= 0; j--)
   {
    Identity.multiplyDirectProduct(X, Y);
    X = Y;
   }
   element = V[i]-mu;
   X.multiply(element);
   Matrix0.add(X);
  }
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void BoseHubbardModel<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const." << endl;
  exit(1);
 }
 Matrix0.setType("hermitian");
}

template<class T> void BoseHubbardModel<T>::getnMPO(unsigned int x, MPO<T>& nMPO) const
{
#ifdef DEBUG
 if ((this->N == 0) || (x > this->N-1))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void BoseHubbardModel<T>::" <<
          "getnMPO(unsigned int x, MPO<T>& nMPO) const: " <<
          "((this->N == 0) || (x > this->N-1))." << endl;
  exit(1);
 }
#endif
// the particle number operator O = n_{l}:
 Matrix<T> O(this->d[0], this->d[0]);
 O.fillZeroes();
 for (int i = 0; i < this->d[0]; i++)
  O(i, i) = T(i);
 getMPOFromLocalOperator(this->BC, this->N, x, O, nMPO);
}

template<class T> void BoseHubbardModel<T>::getNMPO(MPO<T>& NMPO) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void BoseHubbardModel<T>::" <<
          "getNMPO(MPO<T>& NMPO) const: " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
// the particle number operator O = n_{l}:
 Matrix<T> O(this->d[0], this->d[0]);
 O.fillZeroes();
 for (int i = 0; i < this->d[0]; i++)
  O(i, i) = T(i);
 getMPOFromSumLocalOperator(this->BC, this->N, O, NMPO);
}

template<class T> void BoseHubbardModel<T>::getTwoBodyHamiltonian(unsigned int position,
                                                                  Matrix<T>& TwoBodyHamiltonian) const
{
#ifdef DEBUG
 if ((position > this->N-2) || (TwoBodyHamiltonian.getDim0() != this->d[0]*this->d[0]) ||
     (TwoBodyHamiltonian.getDim1() != this->d[0]*this->d[0]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void BoseHubbardModel<T>::" <<
          "getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian) const: " <<
          "((position > this->N-2) || (TwoBodyHamiltonian.getDim0() != this->d[0]*this->d[0]) || " <<
           "(TwoBodyHamiltonian.getDim1() != this->d[0]*this->d[0]))." << endl;
  exit(1);
 }
#endif
 TwoBodyHamiltonian.fillZeroes();
// the interactions:
// - the identity:
 Matrix<T> Identity(this->d[0], this->d[0]);
 Identity.fillZeroes();
 T element = 1.0;
 for (int i = 0; i < this->d[0]; i++)
 {
  Identity(i, i) = element;
 }
// - particle number N =: N0:
 Matrix<T> N0(this->d[0], this->d[0]);
 N0.fillZeroes();
 for (int i = 0; i < this->d[0]; i++)
 {
  element = T(i);
  N0(i, i) = element;
 }
// - interaction particle number 0.5*N(N-1) =: N1:
 Matrix<T> N1(this->d[0], this->d[0]);
 N1.fillZeroes();
 for (int i = 0; i < this->d[0]; i++)
 {
  element = 0.5*T(i*(i-1));
  N1(i, i) = element;
 }
// - A^{+} =: A0:
 Matrix<T> A0(this->d[0], this->d[0]);
 A0.fillZeroes();
 for (int i = 1; i < this->d[0]; i++)
 {
  element = sqrt(T(i));
  A0(i, i-1) = element;
 }
// - A^{-} =: A1:
 Matrix<T> A1(this->d[0], this->d[0]);
 A1.fillZeroes();
 for (int i = 1; i < this->d[0]; i++)
 {
  element = sqrt(T(i));
  A1(i-1, i) = element;
 }
// - A := A^{+} \otimes A^{-} + A^{-} \otimes A^{+}:
 Matrix<T> A(this->d[0]*this->d[0], this->d[0]*this->d[0]), X, Y;
 A.fillZeroes();
 A0.multiplyDirectProduct(A1, X);
 A.add(X);
 A1.multiplyDirectProduct(A0, X);
 A.add(X);
// the parameters:
 vector<T> t(this->N-1);
 for (int i = 0; i < this->N-1; i++)
 {
  t[i] = this->Parameters[i];
 }
 vector<T> U(this->N);
 for (int i = 0; i < this->N; i++)
 {
  U[i] = this->Parameters[this->N-1+i];
 }
 vector<T> V(this->N);
 for (int i = 0; i < this->N; i++)
 {
  V[i] = this->Parameters[2*this->N-1+i];
 }
 T mu = this->Parameters[3*this->N-1];
 if (this->BC == "open")
 {
  if (position == 0)
  {
   N0.multiplyDirectProduct(Identity, X);
   element = V[position]-mu;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N0, X);
   element = 0.5*(V[position+1]-mu);
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   N1.multiplyDirectProduct(Identity, X);
   element = U[position];
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N1, X);
   element = 0.5*U[position+1];
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   X = A;
   element = -t[position];
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
  }
  else if ((position > 0) && (position < this->N-2))
  {
   N0.multiplyDirectProduct(Identity, X);
   element = 0.5*(V[position]-mu);
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N0, X);
   element = 0.5*(V[position+1]-mu);
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   N1.multiplyDirectProduct(Identity, X);
   element = 0.5*U[position];
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N1, X);
   element = 0.5*U[position+1];
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   X = A;
   element = -t[position];
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
  }
  else if (position == this->N-2)
  {
   N0.multiplyDirectProduct(Identity, X);
   element = 0.5*(V[position]-mu);
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N0, X);
   element = V[position+1]-mu;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   N1.multiplyDirectProduct(Identity, X);
   element = 0.5*U[position];
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N1, X);
   element = U[position+1];
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   X = A;
   element = -t[position];
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
  }
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void BoseHubbardModel<T>::" <<
          "getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian) const." << endl;
  exit(1);
 }
}
