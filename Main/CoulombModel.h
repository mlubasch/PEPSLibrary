/// Template class CoulombModel implements spin-1/2 fermions with soft-Coulomb interaction.
/** The template class CoulombModel implements spin-1/2 fermions with a soft-Coulomb interaction:
         H = - \sum_{l, \sigma} t_{l} (c_{l, \sigma}^{+} c_{l+1, \sigma} +
                                       c_{l+1, \sigma}^{+} c_{l, \sigma})
             + \sum_{l} U_{l} n_{l, \spindown} n_{l, \spinup}
             + \sum_{m>l} U_{l,m} n_{l} n_{m}
             + \sum_{l} (V_{l}-\mu) n_{l}   .
    The soft-Coulomb interaction
         c_{0}/sqrt((m\Delta-l\Delta)^2+c_{1})
    is approximated by a sum of r exponentials:
         U_{l,m} = \sum_{i=1}^{r} p_{i} q_{i}^{m-l-1}   .
    We have the following parameters:
    - site-dependent tunneling t_{0}, t_{1}, ..., t_{N-2}
    - site-dependent interaction U_{l}=c_{0}/sqrt(c_{1})
    - soft-Coulomb interaction U_{l,m} defined via c_{0}, c_{1}, \Delta, r
    - site-dependent external potential V_{0}, V_{1}, ..., V_{N-1}
    - chemical potential \mu   .
    \param Representation string, the representation, must be "interactions", "MPO" or "matrix"
    \param BC string, the boundary conditions, must be "open" or "periodic"
    \param N unsigned int, the number of spins, i.e. lattice sites
    \param d vector<unsigned int>, the physical dimensions of the lattice sites, is fixed to 4 everywhere
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

template<class T> class CoulombModel: public Hamiltonian<T>
{
 typedef T (*PointerToFunction)(double time);

 public:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  CoulombModel();

/// Constructor for time-independent CoulombModel with specific BC, N, and Parameters.
/** This constructor initializes a time-independent CoulombModel.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N0 input: unsigned int, the number of lattice sites
    \param Parameters0 input: const vector<T>&, the parameters {t_{0}, ..., t_{N-2},
                                                                c_{0}, c_{1}, \Delta, r,
                                                                V_{0}, ..., V_{N-1}, mu},
                              must fulfill Parameters0.size()==2*N0+4 */
  CoulombModel(const string& BC0, unsigned int N0, const vector<T>& Parameters0);

/// Constructor for time-dependent CoulombModel with specific BC, N, TimeFunctions and time.
/** This constructor initializes a time-dependent CoulombModel.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N0 input: unsigned int, the number of lattice sites
    \param TimeFunctions0 input: const vector<PointerToFunction>&, the time functions for
                                 {t_{0}, ..., t_{N-2}, c_{0}, c_{1}, \Delta, r, V_{0}, ..., V_{N-1}, mu},
                                 must fulfill TimeFunctions.size()==2*N0+4
    \param time0 input: double, the time */
  CoulombModel(const string& BC0, unsigned int N0,
               const vector<PointerToFunction>& TimeFunctions0, double time0);

/// Standard copy constructor.
/** The standard copy constructor copies the input CoulombModel into this.
    \param CoulombModel0 input: const CoulombModel<T>&, to be copied into this
    \sa CoulombModel<T>& operator=(const CoulombModel<T>& CoulombModel0) */
  CoulombModel(const CoulombModel<T>& CoulombModel0);

/// Standard destructor.
/** The standard destructor deletes the elements of CoulombModel. */
  ~CoulombModel();

/// Assigns CoulombModel to this.
/** The operator= allows to assign CoulombModel0 to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side CoulombModel0.
    \param CoulombModel0 input: const CoulombModel<T>&, to be copied into this
    \return CoulombModel<T>&, a reference to the new this
    \sa CoulombModel(const CoulombModel<T>& CoulombModel0) */
  CoulombModel<T>& operator=(const CoulombModel<T>& CoulombModel0);

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

/// Returns interactions representation of CoulombModel.
/** This function returns the interactions of this CoulombModel. The interactions
    are implemented as a vector of positions and a vector of matrices where each vector entry
    corresponds to one term of the sum making up the Hamiltonian.
    \param Positions output: vector< vector<unsigned int> >&, the positions of the matrices
    \param Matrices output: vector< vector< Matrix<T> > >&, the matrices */
  void getInteractions(vector< vector<unsigned int> >& Positions,
                       vector< vector< Matrix<T> > >& Matrices) const;

/// Returns MPO representation of CoulombModel.
/** This function returns this CoulombModel as a MPO.
    \param MPO0 output: MPO<T>&, the MPO representing this CoulombModel */
  void getMPO(MPO<T>& MPO0) const;

/// Returns matrix representation of CoulombModel.
/** This function returns this CoulombModel as a matrix.
    \param Matrix0 output: Matrix<T>&, the matrix representing this CoulombModel,
                           must fulfill (Matrix0.getDim0()==4^{N} && Matrix0.getDim1()==4^{N}),
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

/// Computes ground state of this Hamiltonian without interaction.
/** This function computes the ground state of this Hamiltonian without interaction.
    The matrix of the non-interacting Hamiltonian has dimension 2*this->N and we use its NP lowest-lying
    eigenstates to construct the ground state for NP fermions.
    \param NP input: unsigned int, the total number of particles
    \param energy output: double&, the ground state energy
    \param nP output: vector<double>&, the ground state density */
  void computeNonInteractingGroundState(unsigned int NP, double& energy, vector<double>& nP) const;

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

template<class T> CoulombModel<T>::CoulombModel()
{
 this->N = 0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> CoulombModel<T>::CoulombModel(const string& BC0, unsigned int N0,
                                                const vector<T>& Parameters0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (N0 == 0) || (Parameters0.size() != 2*N0+4))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> CoulombModel<T>::" <<
          "CoulombModel(const string& BC0, unsigned int N0, const vector<T>& Parameters0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (N0 == 0) || (Parameters0.size() != 2*N0+4))." << endl;
  exit(1);
 }
#endif
 this->Representation = "MPO";
 this->BC = BC0;
 this->N = N0;
 this->d = vector<unsigned int>(this->N);
 for (int i = 0; i < this->N; i++)
  this->d[i] = 4;
 this->Parameters = Parameters0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> CoulombModel<T>::CoulombModel(const string& BC0, unsigned int N0,
                                                const vector<PointerToFunction>& TimeFunctions0,
                                                double time0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (N0 == 0) || (TimeFunctions0.size() != 2*N0+4))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> CoulombModel<T>::" <<
          "CoulombModel(const string& BC0, unsigned int N0, " <<
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
  this->d[i] = 4;
 this->Parameters = vector<T>(2*this->N+4);
 this->timeDependent = true;
 this->TimeFunctions = TimeFunctions0;
 this->time = time0;
 for (int i = 0; i < this->Parameters.size(); i++)
  this->Parameters[i] = this->TimeFunctions[i](this->time);
}

template<class T> CoulombModel<T>::CoulombModel(const CoulombModel<T>& CoulombModel0)
{
 this->Representation = CoulombModel0.Representation;
 this->BC = CoulombModel0.BC;
 this->N = CoulombModel0.N;
 this->d = CoulombModel0.d;
 this->Parameters = CoulombModel0.Parameters;
 this->timeDependent = CoulombModel0.timeDependent;
 this->TimeFunctions = CoulombModel0.TimeFunctions;
 this->time = CoulombModel0.time;
}

template<class T> CoulombModel<T>::~CoulombModel() {}

template<class T> CoulombModel<T>& CoulombModel<T>::operator=(const CoulombModel<T>& CoulombModel0)
{
 if (this != &CoulombModel0)
 {
  this->Representation = CoulombModel0.Representation;
  this->BC = CoulombModel0.BC;
  this->N = CoulombModel0.N;
  this->d = CoulombModel0.d;
  this->Parameters = CoulombModel0.Parameters;
  this->timeDependent = CoulombModel0.timeDependent;
  this->TimeFunctions = CoulombModel0.TimeFunctions;
  this->time = CoulombModel0.time;
 }
 return *this;
}

template<class T> void CoulombModel<T>::getSumExp(const vector<T>& y, vector<T>& p, vector<T>& q,
                                                  double& error)
{
#ifdef DEBUG
 if ((y.size() == 0) || (p.size() == 0) || (y.size() <= 2*p.size()) || (p.size() != q.size()))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void CoulombModel<T>::" <<
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

template<class T> void CoulombModel<T>::getInteractions(vector< vector<unsigned int> >& Positions,
                                                        vector< vector< Matrix<T> > >& Matrices) const
{
 cerr << "The following function is not implemented: " <<
         "template<class T> void CoulombModel<T>::" <<
         "getInteractions(vector< vector<unsigned int> >& Positions, " <<
                         "vector< vector< Matrix<T> > >& Matrices) const." << endl;
 exit(1);
}

template<class T> void CoulombModel<T>::getMPO(MPO<T>& MPO0) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void CoulombModel<T>::" <<
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
 CoulombModel<T>::getSumExp(y, p, q, error);
// open boundary conditions:
 if (this->BC == "open")
 {
// MPO operators O^{i}:
  vector<unsigned int> Shape(3);
  Shape[0] = 11; Shape[1] = this->d[0]; Shape[2] = this->d[0];
  Tensor<T> O(Shape);
  O.fillZeroes();
  vector<unsigned int> Index(3);
// - O^{0} = Identity:
  T element = 1.0;
  Index[0] = 0;
  for (int i = 0; i < this->d[0]; i++)
  {
   Index[1] = i; Index[2] = i;
   O.set(Index, element);
  }
// - O^{1} = N_{\spindown}+N_{\spinup}:
  Index[0] = 1;
  Index[1] = 1; Index[2] = 1;
  O.set(Index, element);
  Index[1] = 2; Index[2] = 2;
  O.set(Index, element);
  element = 2.0;
  Index[1] = 3; Index[2] = 3;
  O.set(Index, element);
// - O^{2} = N_{\spindown, \spinup}:
  Index[0] = 2;
  element = 1.0;
  Index[1] = 3; Index[2] = 3;
  O.set(Index, element);
// - O^{3} = C_{\spindown}^{+}:
  Index[0] = 3;
  Index[1] = 0; Index[2] = 1;
  O.set(Index, element);
  Index[1] = 2; Index[2] = 3;
  O.set(Index, element);
// - O^{4} = \tilde{C}_{\spindown}^{+}:
  Index[0] = 4;
  Index[1] = 0; Index[2] = 1;
  O.set(Index, element);
  element = -1.0;
  Index[1] = 2; Index[2] = 3;
  O.set(Index, element);
// - O^{5} = C_{\spindown}:
  Index[0] = 5;
  element = 1.0;
  Index[1] = 1; Index[2] = 0;
  O.set(Index, element);
  Index[1] = 3; Index[2] = 2;
  O.set(Index, element);
// - O^{6} = \tilde{C}_{\spindown}:
  Index[0] = 6;
  Index[1] = 1; Index[2] = 0;
  O.set(Index, element);
  element = -1.0;
  Index[1] = 3; Index[2] = 2;
  O.set(Index, element);
// - O^{7} = C_{\spinup}^{+}:
  Index[0] = 7;
  element = 1.0;
  Index[1] = 0; Index[2] = 2;
  O.set(Index, element);
  Index[1] = 1; Index[2] = 3;
  O.set(Index, element);
// - O^{8} = \tilde{C}_{\spinup}^{+}:
  Index[0] = 8;
  Index[1] = 0; Index[2] = 2;
  O.set(Index, element);
  element = -1.0;
  Index[1] = 1; Index[2] = 3;
  O.set(Index, element);
// - O^{9} = C_{\spinup}:
  Index[0] = 9;
  element = 1.0;
  Index[1] = 2; Index[2] = 0;
  O.set(Index, element);
  Index[1] = 3; Index[2] = 1;
  O.set(Index, element);
// - O^{10} = \tilde{C}_{\spinup}:
  Index[0] = 10;
  Index[1] = 2; Index[2] = 0;
  O.set(Index, element);
  element = -1.0;
  Index[1] = 3; Index[2] = 1;
  O.set(Index, element);
  Tensor<T> OC(O);
// define MPO0:
  unsigned int D = r+6;
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
  Shape[0] = 1; Shape[1] = D; Shape[2] = 11;
  Tensor<T> A(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[2] = 0;
  Index[0] = 0; Index[1] = 0;
  A.set(Index, element);
  Index[2] = 1;
  element = V[position]-mu;
  Index[0] = 0; Index[1] = r+5;
  A.set(Index, element);
  for (int i = 0; i < r; i++)
  {
   element = p[i];
   Index[0] = 0; Index[1] = 5+i;
   A.set(Index, element);
  }
  Index[2] = 2;
  element = c0/sqrt(c1);
  Index[0] = 0; Index[1] = r+5;
  A.set(Index, element);
  element = -t[position];
  Index[2] = 4;
  Index[0] = 0; Index[1] = 1;
  A.set(Index, element);
  Index[2] = 6;
  Index[0] = 0; Index[1] = 2;
  A.set(Index, element);
  Index[2] = 7;
  Index[0] = 0; Index[1] = 3;
  A.set(Index, element);
  Index[2] = 9;
  Index[0] = 0; Index[1] = 4;
  A.set(Index, element);
  vector<unsigned int> IndexA(1), IndexO(1);
  IndexA[0] = 2; IndexO[0] = 0;
  A.contract(IndexA, O, IndexO);
  MPO0.set(position, A);
// A for position 0 < l < N-1:
  Shape[0] = D; Shape[1] = D; Shape[2] = 11;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[2] = 0;
  Index[0] = 0; Index[1] = 0;
  A.set(Index, element);
  for (int i = 0; i < r; i++)
  {
   element = q[i];
   Index[0] = 5+i; Index[1] = 5+i;
   A.set(Index, element);
  }
  element = 1.0;
  Index[0] = r+5; Index[1] = r+5;
  A.set(Index, element);
  Index[2] = 1;
  for (int i = 0; i < r; i++)
  {
   element = p[i];
   Index[0] = 0; Index[1] = 5+i;
   A.set(Index, element);
  }
  element = 1.0;
  for (int i = 0; i < r; i++)
  {
   Index[0] = 5+i; Index[1] = r+5;
   A.set(Index, element);
  }
  Index[2] = 2;
  element = c0/sqrt(c1);
  Index[0] = 0; Index[1] = r+5;
  A.set(Index, element);
  element = 1.0;
  Index[2] = 3;
  Index[0] = 2; Index[1] = r+5;
  A.set(Index, element);
  Index[2] = 5;
  Index[0] = 1; Index[1] = r+5;
  A.set(Index, element);
  Index[2] = 8;
  Index[0] = 4; Index[1] = r+5;
  A.set(Index, element);
  Index[2] = 10;
  Index[0] = 3; Index[1] = r+5;
  A.set(Index, element);
  Tensor<T> AC(A);
  for (position = 1; position < this->N-1; position++)
  {
   A = AC;
   Index[2] = 1;
   element = V[position]-mu;
   Index[0] = 0; Index[1] = r+5;
   A.set(Index, element);
   element = -t[position];
   Index[2] = 4;
   Index[0] = 0; Index[1] = 1;
   A.set(Index, element);
   Index[2] = 6;
   Index[0] = 0; Index[1] = 2;
   A.set(Index, element);
   Index[2] = 7;
   Index[0] = 0; Index[1] = 3;
   A.set(Index, element);
   Index[2] = 9;
   Index[0] = 0; Index[1] = 4;
   A.set(Index, element);
   O = OC;
   A.contract(IndexA, O, IndexO);
   MPO0.set(position, A);
  }
// A for position N-1:
  position = this->N-1;
  Shape[0] = D; Shape[1] = 1; Shape[2] = 11;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[2] = 0;
  Index[0] = r+5; Index[1] = 0;
  A.set(Index, element);
  Index[2] = 1;
  element = V[position]-mu;
  Index[0] = 0; Index[1] = 0;
  A.set(Index, element);
  element = 1.0;
  for (int i = 0; i < r; i++)
  {
   Index[0] = 5+i; Index[1] = 0;
   A.set(Index, element);
  }
  Index[2] = 2;
  element = c0/sqrt(c1);
  Index[0] = 0; Index[1] = 0;
  A.set(Index, element);
  element = 1.0;
  Index[2] = 3;
  Index[0] = 2; Index[1] = 0;
  A.set(Index, element);
  Index[2] = 5;
  Index[0] = 1; Index[1] = 0;
  A.set(Index, element);
  Index[2] = 8;
  Index[0] = 4; Index[1] = 0;
  A.set(Index, element);
  Index[2] = 10;
  Index[0] = 3; Index[1] = 0;
  A.set(Index, element);
  O = OC;
  A.contract(IndexA, O, IndexO);
  MPO0.set(position, A);
 }
// periodic boundary conditions:
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented for periodic BC: " <<
          "template<class T> void CoulombModel<T>::" <<
          "getMPO(MPO<T>& MPO0) const." << endl;
  exit(1);
 }
}

template<class T> void CoulombModel<T>::getMatrix(Matrix<T>& Matrix0) const
{
#ifdef DEBUG
 unsigned int dim = 1;
 for (int i = 0; i < this->N; i++)
  dim *= this->d[0];
 if ((this->N == 0) || (Matrix0.getDim0() != dim) || (Matrix0.getDim1() != dim))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void CoulombModel<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const: " <<
          "((this->N == 0) || (Matrix0.getDim0() != 4^{N}) || (Matrix0.getDim1() != 4^{N}))." << endl;
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
 CoulombModel<T>::getSumExp(y, p, q, error);
// Interactions:
// - Identity:
 Matrix<T> Identity(4, 4);
 Identity.fillZeroes();
 for (int i = 0; i < 4; i++)
  Identity(i, i) = 1.0;
// - particle number N =: N1:
 Matrix<T> N1(4, 4);
 N1.fillZeroes();
 N1(1, 1) = 1.0;
 N1(2, 2) = 1.0;
 N1(3, 3) = 2.0;
// - N_{\spindown \spinup} =: N2:
 Matrix<T> N2(4, 4);
 N2.fillZeroes();
 N2(3, 3) = 1.0;
// - C_{\spindown}^{+} =: C0:
 Matrix<T> C0(4, 4);
 C0.fillZeroes();
 C0(1, 0) = 1.0;
 C0(3, 2) = 1.0;
// - \tilde{C}_{\spindown}^{+} =: C1:
 Matrix<T> C1(4, 4);
 C1.fillZeroes();
 C1(1, 0) = 1.0;
 C1(3, 2) = -1.0;
// - C_{\spindown} =: C2:
 Matrix<T> C2(4, 4);
 C2.fillZeroes();
 C2(0, 1) = 1.0;
 C2(2, 3) = 1.0;
// - \tilde{C}_{\spindown} =: C3:
 Matrix<T> C3(4, 4);
 C3.fillZeroes();
 C3(0, 1) = 1.0;
 C3(2, 3) = -1.0;
// - C_{\spinup}^{+} =: C4:
 Matrix<T> C4(4, 4);
 C4.fillZeroes();
 C4(2, 0) = 1.0;
 C4(3, 1) = 1.0;
// - \tilde{C}_{\spinup}^{+} =: C5:
 Matrix<T> C5(4, 4);
 C5.fillZeroes();
 C5(2, 0) = 1.0;
 C5(3, 1) = -1.0;
// - C_{\spinup} =: C6:
 Matrix<T> C6(4, 4);
 C6.fillZeroes();
 C6(0, 2) = 1.0;
 C6(1, 3) = 1.0;
// - \tilde{C}_{\spinup} =: C7:
 Matrix<T> C7(4, 4);
 C7.fillZeroes();
 C7(0, 2) = 1.0;
 C7(1, 3) = -1.0;
// - C := \tilde{C}_{\spindown}^{+} \otimes C_{\spindown} +
//        \tilde{C}_{\spindown} \otimes C_{\spindown}^{+} +
//        C_{\spinup}^{+} \otimes \tilde{C}_{\spinup} +
//        C_{\spinup} \otimes \tilde{C}_{\spinup}^{+}:
 Matrix<T> C(16, 16), X, Y;
 C.fillZeroes();
 C1.multiplyDirectProduct(C2, X);
 C.add(X);
 C3.multiplyDirectProduct(C0, X);
 C.add(X);
 C4.multiplyDirectProduct(C7, X);
 C.add(X);
 C6.multiplyDirectProduct(C5, X);
 C.add(X);
// parameters:
 vector<T> t(this->N-1);
 for (int i = 0; i < this->N-1; i++)
  t[i] = this->Parameters[i];
 vector<T> V(this->N);
 for (int i = 0; i < this->N; i++)
  V[i] = this->Parameters[this->N+3+i];
 T mu = this->Parameters[2*this->N+3];
 T element;
// open boundary conditions:
 if (this->BC == "open")
 {
// Hamiltonian sum:
// - \sum_{l} -t_{l} (c_{l}^{+} c_{l+1} + c_{l+1}^{+} c_{l}):
  for (int i = 0; i < this->N-1; i++)
  {
   X = C;
   element = -t[i];
   X.multiply(element);
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
   Matrix0.add(X);
  }
// - \sum_{l} U_{l} n_{l, \spindown} n_{l, \spinup}:
  for (int i = 0; i < this->N; i++)
  {
   X = N2;
   element = c0/sqrt(c1);
   X.multiply(element);
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
   Matrix0.add(X);
  }
// - \sum_{m>l} U_{l,m} n_{l} n_{m} with U_{l,m} = \sum_{i=1}^{r} p_{i} q_{i}^{m-l-1}:
  for (int i = 0; i < r; i++)
  {
   for (int l = 0; l < this->N-1; l++)
   {
    for (int m = l+1; m < this->N; m++)
    {
     X = N1;
     element = p[i]*pow(q[i], m-l-1);
     X.multiply(element);
     for (int j = l+1; j < m; j++)
     {
      X.multiplyDirectProduct(Identity, Y);
      X = Y;
     }
     X.multiplyDirectProduct(N1, Y);
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
     Matrix0.add(X);
    }
   }
  }
// - \sum_{l} (V_{l}-mu) n_{l}:
  for (int i = 0; i < this->N; i++)
  {
   X = N1;
   element = V[i]-mu;
   X.multiply(element);
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
   Matrix0.add(X);
  }
 }
// periodic boundary conditions:
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented for periodic BC: " <<
          "template<class T> void CoulombModel<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const." << endl;
  exit(1);
 }
 Matrix0.setType("hermitian");
}

template<class T> void CoulombModel<T>::getnMPO(unsigned int x, MPO<T>& nMPO) const
{
#ifdef DEBUG
 if ((this->N == 0) || (x > this->N-1))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void CoulombModel<T>::" <<
          "getnMPO(unsigned int x, MPO<T>& nMPO) const: " <<
          "((this->N == 0) || (x > this->N-1))." << endl;
  exit(1);
 }
#endif
// the particle number operator O = n_{l}:
 Matrix<T> O(this->d[0], this->d[0]);
 O.fillZeroes();
 O(1, 1) = 1.0;
 O(2, 2) = 1.0;
 O(3, 3) = 2.0;
 getMPOFromLocalOperator(this->BC, this->N, x, O, nMPO);
}

template<class T> void CoulombModel<T>::getNMPO(MPO<T>& NMPO) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void CoulombModel<T>::" <<
          "getNMPO(MPO<T>& NMPO) const: " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
// the particle number operator O = n_{l}:
 Matrix<T> O(this->d[0], this->d[0]);
 O.fillZeroes();
 O(1, 1) = 1.0;
 O(2, 2) = 1.0;
 O(3, 3) = 2.0;
 getMPOFromSumLocalOperator(this->BC, this->N, O, NMPO);
}

template<class T> void CoulombModel<T>::computeNonInteractingGroundState(unsigned int NP, double& energy, vector<double>& nP) const
{
#ifdef DEBUG
 if ((this->BC != "open") || (this->N == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void CoulombModel<T>::" <<
          "computeNonInteractingGroundState(unsigned int NP, double& energy, vector<double>& nP) const: " <<
          "((this->BC != open) || (this->N == 0))." << endl;
  exit(1);
 }
#endif
// parameters:
 vector<T> t(this->N-1);
 for (int i = 0; i < this->N-1; i++)
  t[i] = this->Parameters[i];
 vector<T> V(this->N);
 for (int i = 0; i < this->N; i++)
  V[i] = this->Parameters[this->N+3+i];
 T mu = this->Parameters[2*this->N+3];
// non-interacting Hamiltonian matrix NonIntH:
 unsigned int dim = 2*this->N;
 Matrix<T> NonIntH(dim, dim);
 NonIntH.fillZeroes();
// - diagonal terms for external potential:
 for (int l = 0; l < this->N; l++)
 {
  NonIntH(2*l, 2*l) = V[l]-mu;
  NonIntH(2*l+1, 2*l+1) = V[l]-mu;
 }
// - offdiagonal terms for tunneling:
 for (int l = 0; l < this->N-1; l++)
 {
  NonIntH(2*l+2, 2*l) = -t[l];
  NonIntH(2*l+3, 2*l+1) = -t[l];
  NonIntH(2*l, 2*l+2) = -t[l];
  NonIntH(2*l+1, 2*l+3) = -t[l];
 }
// diagonalize NonIntH:
 NonIntH.setType("hermitian");
 vector<double> W(dim); Matrix<T> Vr(dim, dim);
 NonIntH.eigenDecompose(W, Vr);
// solution:
 energy = 0.0;
 for (int x = 0; x < NP; x++)
  energy += W[x];
 nP = vector<double>(this->N);
 for (int l = 0; l < this->N; l++)
 {
  nP[l] = 0.0;
  for (int x = 0; x < NP; x++)
  {
   nP[l] += Vr(2*l, x)*MathAuxiliary::complexConjugate(Vr(2*l, x));
   nP[l] += Vr(2*l+1, x)*MathAuxiliary::complexConjugate(Vr(2*l+1, x));
  }
 }
}

template<class T> void CoulombModel<T>::getTwoBodyHamiltonian(unsigned int position,
                                                              Matrix<T>& TwoBodyHamiltonian) const
{
 cerr << "The following function is not implemented: " <<
         "template<class T> void CoulombModel<T>::" <<
         "getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian) const." << endl;
 exit(1);
}
