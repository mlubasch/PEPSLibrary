/// Template class SpinlessCoulombModel implements spinless fermions with soft-Coulomb interaction.
/** The template class SpinlessCoulombModel implements spinless fermions with a soft-Coulomb interaction:
         H = - \sum_{l} t_{l} (c_{l}^{+} c_{l+1} + c_{l+1}^{+} c_{l})
             + \sum_{m>l} U_{l,m} n_{l} n_{m}
             + \sum_{l} (V_{l}-\mu) n_{l}   .
    The soft-Coulomb interaction
         c_{0}/sqrt((m\Delta-l\Delta)^2+c_{1})
    is approximated by a sum of r exponentials:
         U_{l,m} = \sum_{i=1}^{r} p_{i} q_{i}^{m-l-1}   .
    We have the following parameters:
    - site-dependent tunneling t_{0}, t_{1}, ..., t_{N-2}
    - soft-Coulomb interaction U_{l,m} defined via c_{0}, c_{1}, \Delta, r
    - site-dependent external potential V_{0}, V_{1}, ..., V_{N-1}
    - chemical potential \mu   .
    \param Representation string, the representation, must be "interactions", "MPO" or "matrix"
    \param BC string, the boundary conditions, must be "open" or "periodic"
    \param N unsigned int, the number of spins, i.e. lattice sites
    \param d vector<unsigned int>, the physical dimensions of the lattice sites, is fixed to 2 everywhere
    \param Parameters vector<T>, the parameters {t_{0}, ..., t_{N-2}, c_{0}, c_{1}, \Delta, r,
                                                 V_{0}, ..., V_{N-1}, mu}
    \param timeDependent bool, the time-dependence
    \param TimeFunctions vector<PointerToFunction>, the functions defining the time-dependence
    \param time double, the time
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

template<class T> class SpinlessCoulombModel: public Hamiltonian<T>
{
 typedef T (*PointerToFunction)(double time);

 public:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  SpinlessCoulombModel();

/// Constructor for time-independent SpinlessCoulombModel with specific BC, N, and Parameters.
/** This constructor initializes a time-independent SpinlessCoulombModel.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N0 input: unsigned int, the number of lattice sites
    \param Parameters0 input: const vector<T>&, the parameters {t_{0}, ..., t_{N-2},
                                                                c_{0}, c_{1}, \Delta, r,
                                                                V_{0}, ..., V_{N-1}, mu},
                              must fulfill Parameters0.size()==2*N0+4 */
  SpinlessCoulombModel(const string& BC0, unsigned int N0, const vector<T>& Parameters0);

/// Constructor for time-dependent SpinlessCoulombModel with specific BC, N, TimeFunctions and time.
/** This constructor initializes a time-dependent SpinlessCoulombModel.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N0 input: unsigned int, the number of lattice sites
    \param TimeFunctions0 input: const vector<PointerToFunction>&, the time functions for
                                 {t_{0}, ..., t_{N-2}, c_{0}, c_{1}, \Delta, r, V_{0}, ..., V_{N-1}, mu},
                                 must fulfill TimeFunctions.size()==2*N0+4
    \param time0 input: double, the time */
  SpinlessCoulombModel(const string& BC0, unsigned int N0,
                       const vector<PointerToFunction>& TimeFunctions0, double time0);

/// Standard copy constructor.
/** The standard copy constructor copies the input SpinlessCoulombModel into this.
    \param SpinlessCoulombModel0 input: const SpinlessCoulombModel<T>&, to be copied into this
    \sa SpinlessCoulombModel<T>& operator=(const SpinlessCoulombModel<T>& SpinlessCoulombModel0) */
  SpinlessCoulombModel(const SpinlessCoulombModel<T>& SpinlessCoulombModel0);

/// Standard destructor.
/** The standard destructor deletes the elements of SpinlessCoulombModel. */
  ~SpinlessCoulombModel();

/// Assigns SpinlessCoulombModel to this.
/** The operator= allows to assign SpinlessCoulombModel0 to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side SpinlessCoulombModel0.
    \param SpinlessCoulombModel0 input: const SpinlessCoulombModel<T>&, to be copied into this
    \return SpinlessCoulombModel<T>&, a reference to the new this
    \sa SpinlessCoulombModel(const SpinlessCoulombModel<T>& SpinlessCoulombModel0) */
  SpinlessCoulombModel<T>& operator=(const SpinlessCoulombModel<T>& SpinlessCoulombModel0);

/// Approximates arbitrary function y by sum of exponentials.
/** An arbitrary function y with m values is approximated by a sum of r exponentials,
         y[x] \approx \sum_{i=1}^{r} p[i] q[i]^{x}   ,
    where the p[i] are the weights and the q[i] are the exponentials.
    The number r is specified via the length of p and q, and it must fulfill (m > 2*r).
    The approximation error is returned as
         error = sqrt(\sum_{x=0}^{m-1}(y[x] - \sum_{i=1}^{r} p[i] q[i]^{x})^{2})   .
    \param y input: const vector<T>&, the function that is to be approximated
    \param p input/output: vector<T>&, on input its length r specifies the number of exponentials,
                                       must fulfill ((p.size() == q.size()) && (y.size() > 2*p.size())),
                                       on output it contains the weights
    \param q input/output: vector<T>&, on input its length r specifies the number of exponentials,
                                       must fulfill ((p.size() == q.size()) && (y.size() > 2*q.size())),
                                       on output it contains the exponentials
    \param error output: double&, the approximation distance */
  static void getSumExp(const vector<T>& y, vector<T>& p, vector<T>& q, double& error);

/// Returns interactions representation of SpinlessCoulombModel.
/** This function returns the interactions of this SpinlessCoulombModel. The interactions
    are implemented as a vector of positions and a vector of matrices where each vector entry
    corresponds to one term of the sum making up the Hamiltonian.
    \param Positions output: vector< vector<unsigned int> >&, the positions of the matrices
    \param Matrices output: vector< vector< Matrix<T> > >&, the matrices */
  void getInteractions(vector< vector<unsigned int> >& Positions,
                       vector< vector< Matrix<T> > >& Matrices) const;

/// Returns MPO representation of SpinlessCoulombModel.
/** This function returns this SpinlessCoulombModel as a MPO.
    \param MPO0 output: MPO<T>&, the MPO representing this SpinlessCoulombModel */
  void getMPO(MPO<T>& MPO0) const;

/// Returns matrix representation of SpinlessCoulombModel.
/** This function returns this SpinlessCoulombModel as a matrix.
    \param Matrix0 output: Matrix<T>&, the matrix representing this SpinlessCoulombModel,
                           must fulfill (Matrix0.getDim0()==2^{N} && Matrix0.getDim1()==2^{N}),
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
    \param TwoBodyHamiltonian output: Matrix<T>&, the two body Hamiltonian, must have the correct form */
  void getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian) const;
};

template<class T> SpinlessCoulombModel<T>::SpinlessCoulombModel()
{
 this->N = 0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> SpinlessCoulombModel<T>::SpinlessCoulombModel(const string& BC0, unsigned int N0,
                                                                const vector<T>& Parameters0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (N0 == 0) || (Parameters0.size() != 2*N0+4))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> SpinlessCoulombModel<T>::" <<
          "SpinlessCoulombModel(const string& BC0, unsigned int N0, const vector<T>& Parameters0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (N0 == 0) || (Parameters0.size() != 2*N0+4))." << endl;
  exit(1);
 }
#endif
 this->Representation = "MPO";
 this->BC = BC0;
 this->N = N0;
 this->d = vector<unsigned int>(this->N);
 for (int i = 0; i < this->N; i++)
  this->d[i] = 2;
 this->Parameters = Parameters0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> SpinlessCoulombModel<T>::SpinlessCoulombModel(const string& BC0, unsigned int N0,
                                              const vector<PointerToFunction>& TimeFunctions0, double time0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (N0 == 0) || (TimeFunctions0.size() != 2*N0+4))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> SpinlessCoulombModel<T>::" <<
          "SpinlessCoulombModel(const string& BC0, unsigned int N0, " <<
                               "const vector<PointerToFunction>& TimeFunctions0, double time0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (N0 == 0) || " <<
           "(TimeFunctions0.size() != 2*N0+4))." << endl;
  exit(1);
 }
#endif
 this->Representation = "MPO";
 this->BC = BC0;
 this->N = N0;
 this->d = vector<unsigned int>(this->N);
 for (int i = 0; i < this->N; i++)
  this->d[i] = 2;
 this->Parameters = vector<T>(2*this->N+4);
 this->timeDependent = true;
 this->TimeFunctions = TimeFunctions0;
 this->time = time0;
 for (int i = 0; i < this->Parameters.size(); i++)
  this->Parameters[i] = this->TimeFunctions[i](this->time);
}

template<class T> SpinlessCoulombModel<T>::SpinlessCoulombModel(const SpinlessCoulombModel<T>& 
                                                                SpinlessCoulombModel0)
{
 this->Representation = SpinlessCoulombModel0.Representation;
 this->BC = SpinlessCoulombModel0.BC;
 this->N = SpinlessCoulombModel0.N;
 this->d = SpinlessCoulombModel0.d;
 this->Parameters = SpinlessCoulombModel0.Parameters;
 this->timeDependent = SpinlessCoulombModel0.timeDependent;
 this->TimeFunctions = SpinlessCoulombModel0.TimeFunctions;
 this->time = SpinlessCoulombModel0.time;
}

template<class T> SpinlessCoulombModel<T>::~SpinlessCoulombModel() {}

template<class T> SpinlessCoulombModel<T>& SpinlessCoulombModel<T>::operator=(const
                                                             SpinlessCoulombModel<T>& SpinlessCoulombModel0)
{
 if (this != &SpinlessCoulombModel0)
 {
  this->Representation = SpinlessCoulombModel0.Representation;
  this->BC = SpinlessCoulombModel0.BC;
  this->N = SpinlessCoulombModel0.N;
  this->d = SpinlessCoulombModel0.d;
  this->Parameters = SpinlessCoulombModel0.Parameters;
  this->timeDependent = SpinlessCoulombModel0.timeDependent;
  this->TimeFunctions = SpinlessCoulombModel0.TimeFunctions;
  this->time = SpinlessCoulombModel0.time;
 }
 return *this;
}

template<class T> void SpinlessCoulombModel<T>::getSumExp(const vector<T>& y, vector<T>& p, vector<T>& q,
                                                          double& error)
{
#ifdef DEBUG
 if ((y.size() == 0) || (p.size() == 0) || (y.size() <= 2*p.size()) || (p.size() != q.size()))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void SpinlessCoulombModel<T>::" <<
          "getSumExp(const vector<T>& y, vector<T>& p, vector<T>& q, double& error): " <<
          "((y.size() == 0) || (p.size() == 0) || (y.size() <= 2*p.size()) || " <<
           "(p.size() != q.size()))." << endl;
  exit(1);
 }
#endif
 unsigned int m = y.size(), r = p.size();
 double rcond = 1.0e-12; unsigned int rank; vector<double> S(r);
// 1. Find optimal s values:
 Matrix<T> A(m-r, r);
 for (int i = 0; i < m-r; i++)
 {
  for (int j = 0; j < r; j++)
  {
   A(i, j) = y[i+j];
  }
 }
 Matrix<T> b(m-r, 1);
 for (int i = 0; i < m-r; i++)
  b(i, 0) = -y[r+i];
 A.linearLeastSquares(rcond, b, rank, S);
// 2. Find optimal exponentials q[i]:
 Matrix<T> C(r, r);
 C.fillZeroes();
 for (int j = 0; j < r; j++)
  C(0, j) = -b(r-1-j, 0);
 for (int i = 1; i < r; i++)
  C(i, i-1) = 1.0;
 vector< complex<double> > W(r);
 Matrix<T> Vr(r, r), Vl(r, r);
 C.eigenDecompose(W, Vr, Vl);
 if ((typeid(T) == typeid(float)) || (typeid(T) == typeid(double)))
 {
  for (int i = 0; i < r; i++)
   q[i] = abs(W[i]);
 }
 else
 {
  for (int i = 0; i < r; i++)
   MathAuxiliary::convertComplex(W[i], q[i]);
 }
// 3. Find optimal weights p[i]:
 Matrix<T> B(m, r);
 for (int i = 0; i < m; i++)
 {
  for (int j = 0; j < r; j++)
  {
   B(i, j) = pow(q[j], i);
  }
 }
 Matrix<T> Y(m, 1);
 for (int i = 0; i < m; i++)
  Y(i, 0) = y[i];
 B.linearLeastSquares(rcond, Y, rank, S);
 for (int i = 0; i < r; i++)
  p[i] = Y(i, 0);
// compute error:
 T sumExp;
 error = 0.0;
 for (int i = 0; i < m; i++)
 {
  sumExp = 0.0;
  for (int j = 0; j < r; j++)
   sumExp += p[j]*pow(q[j], i);
  error += pow(abs(y[i]-sumExp), 2);
 }
 error = sqrt(error);
}

template<class T> void SpinlessCoulombModel<T>::getInteractions(vector< vector<unsigned int> >& Positions,
                                                              vector< vector< Matrix<T> > >& Matrices) const
{
 cerr << "The following function is not implemented: " <<
         "template<class T> void SpinlessCoulombModel<T>::" <<
         "getInteractions(vector< vector<unsigned int> >& Positions, " <<
                         "vector< vector< Matrix<T> > >& Matrices) const." << endl;
 exit(1);
}

template<class T> void SpinlessCoulombModel<T>::getMPO(MPO<T>& MPO0) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void SpinlessCoulombModel<T>::" <<
          "getMPO(MPO<T>& MPO0) const: " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
// get approximation of soft-Coulomb interaction by sum of r exponentials:
 T c0 = this->Parameters[this->N-1], c1 = this->Parameters[this->N], Delta = this->Parameters[this->N+1];
 unsigned int r = this->Parameters[this->N+2];
 unsigned int m = max(this->N-1, 2*r+1);
 vector<T> y(m);
 for (int i = 0; i < m; i++)
  y[i] = c0 / sqrt(pow((i+1)*Delta, 2)+c1);
 vector<T> p(r), q(r);
 double error;
 SpinlessCoulombModel<T>::getSumExp(y, p, q, error);
// open boundary conditions:
 if (this->BC == "open")
 {
// MPO operators O^{i}:
  vector<unsigned int> Shape(3);
  Shape[0] = 3*r+4; Shape[1] = this->d[0]; Shape[2] = this->d[0];
  Tensor<T> O(Shape);
  O.fillZeroes();
  vector<unsigned int> Index(3);
// O^{0} = Identity:
  T element = 1.0;
  Index[0] = 0;
  Index[1] = 0; Index[2] = 0;
  O.set(Index, element);
  Index[1] = 1; Index[2] = 1;
  O.set(Index, element);
// O^{1} = N:
  Index[0] = 1;
  Index[1] = 1; Index[2] = 1;
  O.set(Index, element);
// O^{2} = C^{+}:
  Index[0] = 2;
  Index[1] = 0; Index[2] = 1;
  O.set(Index, element);
// O^{3} = C:
  Index[0] = 3;
  Index[1] = 1; Index[2] = 0;
  O.set(Index, element);
// O^{4}, O^{7}, ..., O^{3r+1} = p_{i}*N:
  for (int i = 0; i < r; i++)
  {
   element = p[i];
   Index[0] = 4+3*i;
   Index[1] = 1; Index[2] = 1;
   O.set(Index, element);
  }
// O^{5}, O^{8}, ..., O^{3r+2} = q_{i}*Identity:
  for (int i = 0; i < r; i++)
  {
   element = q[i];
   Index[0] = 5+3*i;
   Index[1] = 0; Index[2] = 0;
   O.set(Index, element);
   Index[1] = 1; Index[2] = 1;
   O.set(Index, element);
  }
// O^{6}, O^{9}, ..., O^{3r+3} = N:
  element = 1.0;
  for (int i = 0; i < r; i++)
  {
   Index[0] = 6+3*i;
   Index[1] = 1; Index[2] = 1;
   O.set(Index, element);
  }
  unsigned int D = r+4;
  MPO0 = MPO<T>(this->BC, this->N, this->d, D);
// parameters:
  vector<T> t(this->N-1);
  for (int i = 0; i < this->N-1; i++)
   t[i] = this->Parameters[i];
  vector<T> V(this->N);
  for (int i = 0; i < this->N; i++)
   V[i] = this->Parameters[this->N+3+i];
  T mu = this->Parameters[2*this->N+3];
// A for position 0:
  unsigned int position = 0;
  Shape[0] = 1; Shape[1] = D; Shape[2] = 3*r+4;
  Tensor<T> A(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  element = V[position] - mu;
  Index[0] = 0; Index[1] = r+3; Index[2] = 1;
  A.set(Index, element);
  element = -t[position];
  Index[0] = 0; Index[1] = 1; Index[2] = 2;
  A.set(Index, element);
  Index[0] = 0; Index[1] = 2; Index[2] = 3;
  A.set(Index, element);
  element = 1.0;
  for (int i = 0; i < r; i++)
  {
   Index[0] = 0; Index[1] = 3+i; Index[2] = 4+3*i;
   A.set(Index, element);
  }
  Tensor<T> OC(O);
  vector<unsigned int> IndexA(1), IndexOC(1);
  IndexA[0] = 2; IndexOC[0] = 0;
  A.contract(IndexA, OC, IndexOC);
  MPO0.set(position, A);
// A for position 0 < l < N-1:
  Shape[0] = D; Shape[1] = D; Shape[2] = 3*r+4;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  Index[0] = r+3; Index[1] = r+3; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 2; Index[1] = r+3; Index[2] = 2;
  A.set(Index, element);
  Index[0] = 1; Index[1] = r+3; Index[2] = 3;
  A.set(Index, element);
  for (int i = 0; i < r; i++)
  {
   Index[0] = 0; Index[1] = 3+i; Index[2] = 4+3*i;
   A.set(Index, element);
   Index[0] = 3+i; Index[1] = 3+i; Index[2] = 5+3*i;
   A.set(Index, element);
   Index[0] = 3+i; Index[1] = r+3; Index[2] = 6+3*i;
   A.set(Index, element);
  }
  Tensor<T> AC(A);
  for (position = 1; position < this->N-1; position++)
  {
   A = AC;
   element = V[position] - mu;
   Index[0] = 0; Index[1] = r+3; Index[2] = 1;
   A.set(Index, element);
   element = -t[position];
   Index[0] = 0; Index[1] = 1; Index[2] = 2;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 2; Index[2] = 3;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   MPO0.set(position, A);
  }
// A for position N-1:
  position = this->N-1;
  Shape[0] = D; Shape[1] = 1; Shape[2] = 3*r+4;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = r+3; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  element = V[position] - mu;
  Index[0] = 0; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  element = 1.0;
  Index[0] = 2; Index[1] = 0; Index[2] = 2;
  A.set(Index, element);
  Index[0] = 1; Index[1] = 0; Index[2] = 3;
  A.set(Index, element);
  for (int i = 0; i < r; i++)
  {
   Index[0] = 3+i; Index[1] = 0; Index[2] = 6+3*i;
   A.set(Index, element);
  }
  OC = O;
  A.contract(IndexA, OC, IndexOC);
  MPO0.set(position, A);
 }
// periodic boundary conditions:
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented for periodic BC: " <<
          "template<class T> void SpinlessCoulombModel<T>::" <<
          "getMPO(MPO<T>& MPO0) const." << endl;
  exit(1);
 }
}

template<class T> void SpinlessCoulombModel<T>::getMatrix(Matrix<T>& Matrix0) const
{
#ifdef DEBUG
 unsigned int dim = 1;
 for (int i = 0; i < this->N; i++)
  dim *= this->d[0];
 if ((this->N == 0) || (Matrix0.getDim0() != dim) || (Matrix0.getDim1() != dim))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void SpinlessCoulombModel<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const: " <<
          "((this->N == 0) || (Matrix0.getDim0() != dim) || (Matrix0.getDim1() != dim))." << endl;
  exit(1);
 }
#endif
 Matrix0.fillZeroes();
// get approximation of soft-Coulomb interaction by sum of r exponentials:
 T c0 = this->Parameters[this->N-1], c1 = this->Parameters[this->N], Delta = this->Parameters[this->N+1];
 unsigned int r = this->Parameters[this->N+2];
 unsigned int m = max(this->N-1, 2*r+1);
 vector<T> y(m);
 for (int i = 0; i < m; i++)
  y[i] = c0 / sqrt(pow((i+1)*Delta, 2)+c1);
 vector<T> p(r), q(r);
 double error;
 SpinlessCoulombModel<T>::getSumExp(y, p, q, error);
// Interactions:
// - Identity:
 Matrix<T> Identity(this->d[0], this->d[0]);
 Identity.fillZeroes();
 T element = 1.0;
 Identity(0, 0) = element;
 Identity(1, 1) = element;
// - particle number N =: N0:
 Matrix<T> N0(this->d[0], this->d[0]);
 N0.fillZeroes();
 N0(1, 1) = element;
// - C^{+} =: C0:
 Matrix<T> C0(this->d[0], this->d[0]);
 C0.fillZeroes();
 C0(1, 0) = element;
// - C =: C1:
 Matrix<T> C1(this->d[0], this->d[0]);
 C1.fillZeroes();
 C1(0, 1) = element;
// - C := C^{+} \otimes C + C \otimes C^{+}:
 Matrix<T> C(this->d[0]*this->d[0], this->d[0]*this->d[0]), X, Y;
 C.fillZeroes();
 C0.multiplyDirectProduct(C1, X);
 C.add(X);
 C1.multiplyDirectProduct(C0, X);
 C.add(X);
// parameters:
 vector<T> t(this->N-1);
 for (int i = 0; i < this->N-1; i++)
  t[i] = this->Parameters[i];
 vector<T> V(this->N);
 for (int i = 0; i < this->N; i++)
  V[i] = this->Parameters[this->N+3+i];
 T mu = this->Parameters[2*this->N+3];
// open boundary conditions:
 if (this->BC == "open")
 {
// Hamiltonian sum:
// - \sum_{l} -t_{l} (c_{l}^{+} c_{l+1} + c_{l+1}^{+} c_{l}):
  for (int i = 0; i < this->N-1; i++)
  {
   X = C;
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
// - \sum_{m>l} U_{l,m} n_{l} n_{m} with U_{l,m} = \sum_{i=1}^{r} p_{i} q_{i}^{m-l-1}:
  for (int i = 0; i < r; i++)
  {
   for (int l = 0; l < this->N-1; l++)
   {
    for (int m = l+1; m < this->N; m++)
    {
     X = N0;
     for (int j = l+1; j < m; j++)
     {
      X.multiplyDirectProduct(Identity, Y);
      X = Y;
     }
     X.multiplyDirectProduct(N0, Y);
     X = Y;
     for (int j = m+1; j < this->N; j++)
     {
      X.multiplyDirectProduct(Identity, Y);
      X = Y;
     }
     for (int j = l-1; j >= 0; j--)
     {
      Identity.multiplyDirectProduct(X, Y);
      X = Y;
     }
     element = p[i]*pow(q[i], m-l-1);
     X.multiply(element);
     Matrix0.add(X);
    }
   }
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
// periodic boundary conditions:
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented for periodic BC: " <<
          "template<class T> void SpinlessCoulombModel<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const." << endl;
  exit(1);
 }
 Matrix0.setType("hermitian");
}

template<class T> void SpinlessCoulombModel<T>::getnMPO(unsigned int x, MPO<T>& nMPO) const
{
#ifdef DEBUG
 if ((this->N == 0) || (x > this->N-1))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void SpinlessCoulombModel<T>::" <<
          "getnMPO(unsigned int x, MPO<T>& nMPO) const: " <<
          "((this->N == 0) || (x > this->N-1))." << endl;
  exit(1);
 }
#endif
// the particle number operator O = n_{l}:
 Matrix<T> O(this->d[0], this->d[0]);
 O.fillZeroes();
 O(1, 1) = 1.0;
 getMPOFromLocalOperator(this->BC, this->N, x, O, nMPO);
}

template<class T> void SpinlessCoulombModel<T>::getNMPO(MPO<T>& NMPO) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void SpinlessCoulombModel<T>::" <<
          "getNMPO(MPO<T>& NMPO) const: " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
// the particle number operator O = n_{l}:
 Matrix<T> O(this->d[0], this->d[0]);
 O.fillZeroes();
 O(1, 1) = 1.0;
 getMPOFromSumLocalOperator(this->BC, this->N, O, NMPO);
}

template<class T> void SpinlessCoulombModel<T>::getTwoBodyHamiltonian(unsigned int position,
                                                                  Matrix<T>& TwoBodyHamiltonian) const
{
 cerr << "The following function is not implemented: " <<
         "template<class T> void SpinlessCoulombModel<T>::" <<
         "getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian) const." << endl;
 exit(1);
}
