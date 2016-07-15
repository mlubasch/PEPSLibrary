/// Template class HubbardModel implements fermionic Hubbard model.
/** The template class HubbardModel implements the fermionic Hubbard model:
         H = - \sum_{l, \sigma} t_{l} (c_{l, \sigma}^{+} c_{l+1, \sigma} +
                                       c_{l+1, \sigma}^{+} c_{l, \sigma})
             + \sum_{l} U_{l} n_{l, \spindown} n_{l, \spinup}
             + \sum_{l} (V_{l}-\mu) n_{l}   .
    We have the following parameters:
    - site-dependent tunneling t_{0}, t_{1}, ..., t_{N-2}
    - site-dependent on-site interaction U_{0}, U_{1}, ..., U_{N-1}
    - site-dependent external potential V_{0}, V_{1}, ..., V_{N-1}
    - chemical potential \mu   .
    \param Representation string, the representation, must be "interactions", "MPO" or "matrix"
    \param BC string, the boundary conditions, must be "open" or "periodic"
    \param N unsigned int, the number of spins, i.e. lattice sites
    \param d vector<unsigned int>, the physical dimensions of the lattice sites, is fixed to 4 everywhere
    \param Parameters vector<T>, the parameters {t_{0}, ..., t_{N-2}, U_{0}, ..., U_{N-1},
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

template<class T> class HubbardModel: public Hamiltonian<T>
{
 typedef T (*PointerToFunction)(double time);

 public:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  HubbardModel();

/// Constructor for time-independent HubbardModel with specific BC, N, and Parameters.
/** This constructor initializes a time-independent HubbardModel.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N0 input: unsigned int, the number of lattice sites
    \param Parameters0 input: const vector<T>&, the parameters {t_{0}, ..., t_{N-2},
                                                                U_{0}, ..., U_{N-1},
                                                                V_{0}, ..., V_{N-1}, mu},
                              must fulfill Parameters0.size()==3*N0 */
  HubbardModel(const string& BC0, unsigned int N0, const vector<T>& Parameters0);

/// Constructor for time-dependent HubbardModel with specific BC, N, TimeFunctions and time.
/** This constructor initializes a time-dependent HubbardModel.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N0 input: unsigned int, the number of lattice sites
    \param TimeFunctions0 input: const vector<PointerToFunction>&, the time functions for
                                 {t_{0}, ..., t_{N-2}, U_{0}, ..., U_{N-1}, V_{0}, ..., V_{N-1}, mu},
                                 must fulfill TimeFunctions.size()==3*N0
    \param time0 input: double, the time */
  HubbardModel(const string& BC0, unsigned int N0,
               const vector<PointerToFunction>& TimeFunctions0, double time0);

/// Standard copy constructor.
/** The standard copy constructor copies the input HubbardModel into this.
    \param HubbardModel0 input: const HubbardModel<T>&, to be copied into this
    \sa HubbardModel<T>& operator=(const HubbardModel<T>& HubbardModel0) */
  HubbardModel(const HubbardModel<T>& HubbardModel0);

/// Standard destructor.
/** The standard destructor deletes the elements of HubbardModel. */
  ~HubbardModel();

/// Assigns HubbardModel to this.
/** The operator= allows to assign HubbardModel0 to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side HubbardModel0.
    \param HubbardModel0 input: const HubbardModel<T>&, to be copied into this
    \return HubbardModel<T>&, a reference to the new this
    \sa HubbardModel(const HubbardModel<T>& HubbardModel0) */
  HubbardModel<T>& operator=(const HubbardModel<T>& HubbardModel0);

/// Returns interactions representation of HubbardModel.
/** This function returns the interactions of this HubbardModel. The interactions
    are implemented as a vector of positions and a vector of matrices where each vector entry
    corresponds to one term of the sum making up the Hamiltonian.
    \param Positions output: vector< vector<unsigned int> >&, the positions of the matrices
    \param Matrices output: vector< vector< Matrix<T> > >&, the matrices */
  void getInteractions(vector< vector<unsigned int> >& Positions,
                       vector< vector< Matrix<T> > >& Matrices) const;

/// Returns MPO representation of HubbardModel.
/** This function returns this HubbardModel as a MPO.
    \param MPO0 output: MPO<T>&, the MPO representing this HubbardModel */
  void getMPO(MPO<T>& MPO0) const;

/// Returns matrix representation of HubbardModel.
/** This function returns this HubbardModel as a matrix.
    \param Matrix0 output: Matrix<T>&, the matrix representing this HubbardModel,
                           must fulfill (Matrix0.getDim0()==4^{N} && Matrix0.getDim1()==4^{N}),
                           will be of type "hermitian" */
  void getMatrix(Matrix<T>& Matrix0) const;

/// Returns sparse matrix representation of Hubbard model.
/** This function returns the Hubbard model as a sparse matrix.
    Given homogeneous tunneling t, on-site interaction U, and chemical potential mu, the non-zero elements
    of the corresponding Hubbard Hamiltonian are returned.
    The non-zero offdiagonal elements of the upper half of the kinetic energy Ht will be indexed by Indexti
    and Indextj such that the kth element in Ht, i.e. Ht(k), is located at row Indexti(k) and column
    Indextj(k).
    The non-zero diagonal elements of the on-site interaction HU will be indexed by IndexU such that the
    kth element in HU, i.e. HU(k), is located at row and column IndexU(k). The same convention applies to
    Hmu and Indexmu.
    \param BC input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N input: unsigned int, the number of lattice sites
    \param t input: T, the tunneling
    \param U input: T, the on-site interaction
    \param mu input: T, the chemical potential
    \param Ht output: vector<T>&, the non-zero offdiagonal elements of the upper half of the kinetic energy,
                      indexed by Indexti and Indextj
    \param Indexti output: vector<int>&, the row index to Ht
    \param Indextj output: vector<int>&, the column index to Ht
    \param HU output: vector<T>&, the non-zero diagonal elements of the on-site interaction,
                      indexed by IndexU
    \param IndexU output: vector<int>&, the diagonal index to HU
    \param Hmu output: vector<T>&, the non-zero diagonal elements of the chemical potential,
                       indexed by Indexmu
    \param Indexmu output: vector<int>&, the diagonal index to Hmu */
  static void getSparseHubbardMatrix(const string& BC, unsigned int N, T t, T U, T mu,
                                     vector<T>& Ht, vector<int>& Indexti, vector<int>& Indextj,
                                     vector<T>& HU, vector<int>& IndexU,
                                     vector<T>& Hmu, vector<int>& Indexmu);

/// Returns FHK for Hubbard model.
/** This function computes FHK for the Hubbard model by exact diagonalization of the sparse Hubbard matrix
    from
       static void HubbardModel<T>::getSparseHubbardMatrix(const string& BC, unsigned int N, T t, T U, T mu,
                                                           vector<T>& Ht, vector<int>& Indexti,
                                                           vector<int>& Indextj,
                                                           vector<T>& HU, vector<int>& IndexU,
                                                           vector<T>& Hmu, vector<int>& Indexmu)   .
    For every desired site occupation n_{l}^{desired}, the interacting inversion is implemented as the
    following iteration:
       v_{l}(run+1) = v_{l}(run)+alpha(run)(x)*(n_{l}(run)-n_{l}^{desired})   .
    In each run we perform two diagonalizations: one with alpha(run)(x)=x*alpha(run) and one with
    alpha(run)(x)=alpha(run)/x, and we take the alpha with the smaller occupation error
       error := \sqrt{\sum_{l}(n_{l}-n_{l}^{desired})^{2}}
    for the next run.
    The function evaluates the occupation error after each run, and stops when it is smaller than eps.
    It always stops when the number of runs exceeds maxNumRuns.
    The template parameter T must be real.
    \param BC input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N input: unsigned int, the number of lattice sites
    \param t input: T, the tunneling
    \param U input: T, the on-site interaction
    \param mu input: T, the chemical potential
    \param NP input: unsigned int, the total number of particles
    \param nPStart input: T, the first grid point of the local particle number nP
    \param nPStop input: T, the last grid point of the local particle number nP
    \param R input: unsigned int, R+1 is the number of grid points in the local particle number nP
    \param alphaStart input: T, the first alpha used in the interacting inversion
    \param x input: T, the factor multiplying alpha
    \param eps input: T, the occupation precision
    \param maxNumRuns input: unsigned int, the maximal number of runs allowed per desired site occupation
                             n_{l}^{desired}
    \param ErrorAchieved output: Tensor<T>&, the achieved occupation errors, must have the correct shape
    \param NumRunsDone output: Tensor<unsigned int>&, the number of runs done, must have the correct shape
    \param FHK output: Tensor<T>&, the Hohenberg-Kohn functional FHK, must have the correct shape */
  static void getFHKHubbard(const string& BC, unsigned int N, T t, T U, T mu, unsigned int NP,
                            T nPStart, T nPStop, unsigned int R,
                            T alphaStart, T x, T eps, unsigned int maxNumRuns,
                            Tensor<T>& ErrorAchieved, Tensor<unsigned int>& NumRunsDone, Tensor<T>& FHK);

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

template<class T> HubbardModel<T>::HubbardModel()
{
 this->N = 0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> HubbardModel<T>::HubbardModel(const string& BC0, unsigned int N0,
                                                const vector<T>& Parameters0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (N0 == 0) || (Parameters0.size() != 3*N0))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> HubbardModel<T>::" <<
          "HubbardModel(const string& BC0, unsigned int N0, const vector<T>& Parameters0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (N0 == 0) || (Parameters0.size() != 3*N0))." << endl;
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

template<class T> HubbardModel<T>::HubbardModel(const string& BC0, unsigned int N0,
                                                const vector<PointerToFunction>& TimeFunctions0,
                                                double time0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (N0 == 0) || (TimeFunctions0.size() != 3*N0))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> HubbardModel<T>::" <<
          "HubbardModel(const string& BC0, unsigned int N0, " <<
                       "const vector<PointerToFunction>& TimeFunctions0, double time0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (N0 == 0) || " <<
           "(TimeFunctions0.size() != 3*N0))." << endl;
  exit(1);
 }
#endif
 this->Representation = "MPO";
 this->BC = BC0;
 this->N = N0;
 this->d = vector<unsigned int>(this->N);
 for (int i = 0; i < this->N; i++)
  this->d[i] = 4;
 this->Parameters = vector<T>(3*this->N);
 this->timeDependent = true;
 this->TimeFunctions = TimeFunctions0;
 this->time = time0;
 for (int i = 0; i < this->Parameters.size(); i++)
  this->Parameters[i] = this->TimeFunctions[i](this->time);
}

template<class T> HubbardModel<T>::HubbardModel(const HubbardModel<T>& HubbardModel0)
{
 this->Representation = HubbardModel0.Representation;
 this->BC = HubbardModel0.BC;
 this->N = HubbardModel0.N;
 this->d = HubbardModel0.d;
 this->Parameters = HubbardModel0.Parameters;
 this->timeDependent = HubbardModel0.timeDependent;
 this->TimeFunctions = HubbardModel0.TimeFunctions;
 this->time = HubbardModel0.time;
}

template<class T> HubbardModel<T>::~HubbardModel() {}

template<class T> HubbardModel<T>& HubbardModel<T>::operator=(const HubbardModel<T>& HubbardModel0)
{
 if (this != &HubbardModel0)
 {
  this->Representation = HubbardModel0.Representation;
  this->BC = HubbardModel0.BC;
  this->N = HubbardModel0.N;
  this->d = HubbardModel0.d;
  this->Parameters = HubbardModel0.Parameters;
  this->timeDependent = HubbardModel0.timeDependent;
  this->TimeFunctions = HubbardModel0.TimeFunctions;
  this->time = HubbardModel0.time;
 }
 return *this;
}

template<class T> void HubbardModel<T>::getInteractions(vector< vector<unsigned int> >& Positions,
                                                        vector< vector< Matrix<T> > >& Matrices) const
{
 cerr << "The following function is not implemented: " <<
         "template<class T> void HubbardModel<T>::" <<
         "getInteractions(vector< vector<unsigned int> >& Positions, " <<
                         "vector< vector< Matrix<T> > >& Matrices) const." << endl;
 exit(1);
}

template<class T> void HubbardModel<T>::getMPO(MPO<T>& MPO0) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void HubbardModel<T>::" <<
          "getMPO(MPO<T>& MPO0) const: " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
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
  unsigned int D = 6;
  MPO0 = MPO<T>(this->BC, this->N, this->d, D);
// parameters:
  vector<T> t(this->N-1), U(this->N), V(this->N);
  for (int i = 0; i < this->N-1; i++)
   t[i] = this->Parameters[i];
  for (int i = 0; i < this->N; i++)
  {
   U[i] = this->Parameters[this->N-1+i];
   V[i] = this->Parameters[2*this->N-1+i];
  }
  T mu = this->Parameters[3*this->N-1];
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
  Index[0] = 0; Index[1] = 5;
  A.set(Index, element);
  Index[2] = 2;
  element = U[position];
  Index[0] = 0; Index[1] = 5;
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
  Index[0] = 5; Index[1] = 5;
  A.set(Index, element);
  Index[2] = 3;
  Index[0] = 2; Index[1] = 5;
  A.set(Index, element);
  Index[2] = 5;
  Index[0] = 1; Index[1] = 5;
  A.set(Index, element);
  Index[2] = 8;
  Index[0] = 4; Index[1] = 5;
  A.set(Index, element);
  Index[2] = 10;
  Index[0] = 3; Index[1] = 5;
  A.set(Index, element);
  Tensor<T> AC(A);
  for (position = 1; position < this->N-1; position++)
  {
   A = AC;
   Index[2] = 1;
   element = V[position]-mu;
   Index[0] = 0; Index[1] = 5;
   A.set(Index, element);
   Index[2] = 2;
   element = U[position];
   Index[0] = 0; Index[1] = 5;
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
  Index[0] = 5; Index[1] = 0;
  A.set(Index, element);
  Index[2] = 1;
  element = V[position]-mu;
  Index[0] = 0; Index[1] = 0;
  A.set(Index, element);
  Index[2] = 2;
  element = U[position];
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
          "template<class T> void HubbardModel<T>::" <<
          "getMPO(MPO<T>& MPO0) const." << endl;
  exit(1);
 }
}

template<class T> void HubbardModel<T>::getMatrix(Matrix<T>& Matrix0) const
{
#ifdef DEBUG
 unsigned int dim = 1;
 for (int i = 0; i < this->N; i++)
  dim *= this->d[0];
 if ((this->N == 0) || (Matrix0.getDim0() != dim) || (Matrix0.getDim1() != dim))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void HubbardModel<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const: " <<
          "((this->N == 0) || (Matrix0.getDim0() != 4^{N}) || (Matrix0.getDim1() != 4^{N}))." << endl;
  exit(1);
 }
#endif
 Matrix0.fillZeroes();
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
 vector<T> t(this->N-1), U(this->N), V(this->N);
 for (int i = 0; i < this->N-1; i++)
  t[i] = this->Parameters[i];
 for (int i = 0; i < this->N; i++)
 {
  U[i] = this->Parameters[this->N-1+i];
  V[i] = this->Parameters[2*this->N-1+i];
 }
 T mu = this->Parameters[3*this->N-1];
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
   element = U[i];
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
          "template<class T> void HubbardModel<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const." << endl;
  exit(1);
 }
 Matrix0.setType("hermitian");
}

template<class T> void HubbardModel<T>::getSparseHubbardMatrix(const string& BC, unsigned int N, T t, T U,
                                                               T mu,
                                                               vector<T>& Ht, vector<int>& Indexti,
                                                               vector<int>& Indextj,
                                                               vector<T>& HU, vector<int>& IndexU,
                                                               vector<T>& Hmu, vector<int>& Indexmu)
{
#ifdef DEBUG
 if (((BC != "open") && (BC != "periodic")) || (N == 0))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void HubbardModel<T>::" <<
          "getSparseHubbardMatrix(const string& BC, unsigned int N, T t, T U, T mu, " <<
                                 "vector<T>& Ht, vector<int>& Indexti, vector<int>& Indextj, " <<
                                 "vector<T>& HU, vector<int>& IndexU, vector<T>& Hmu, " <<
                                 "vector<int>& Indexmu): " <<
          "(((BC != open) && (BC != periodic)) || (N == 0))." << endl;
  exit(1);
 }
#endif
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
// - C_{\spindown}^{-} =: C2:
 Matrix<T> C2(4, 4);
 C2.fillZeroes();
 C2(0, 1) = 1.0;
 C2(2, 3) = 1.0;
// - \tilde{C}_{\spindown}^{-} =: C3:
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
// - C_{\spinup}^{-} =: C6:
 Matrix<T> C6(4, 4);
 C6.fillZeroes();
 C6(0, 2) = 1.0;
 C6(1, 3) = 1.0;
// - \tilde{C}_{\spinup}^{-} =: C7:
 Matrix<T> C7(4, 4);
 C7.fillZeroes();
 C7(0, 2) = 1.0;
 C7(1, 3) = -1.0;
// - C := \tilde{C}_{\spindown}^{+} \otimes C_{\spindown}^{-} +
//        \tilde{C}_{\spindown}^{-} \otimes C_{\spindown}^{+} +
//        C_{\spinup}^{+} \otimes \tilde{C}_{\spinup}^{-} +
//        C_{\spinup}^{-} \otimes \tilde{C}_{\spinup}^{+}:
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
// multiply Interactions with parameters:
 N1.multiply(-mu);
 N2.multiply(U);
 C.multiply(-t);
// obtain non-zero elements:
 unsigned int dim = pow(double(4), double(N));
 Matrix<T> Matrix0(dim, dim);
 unsigned int count;
// open boundary conditions:
 if (BC == "open")
 {
// Hamiltonian sum:
// - -t \sum_{l} (c_{l}^{+} c_{l+1} + c_{l+1}^{+} c_{l}):
  Matrix0.fillZeroes();
  for (int i = 0; i < N-1; i++)
  {
   X = C;
   for (int j = i+2; j < N; j++)
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
  count = 0;
  for (int j = 0; j < dim; j++)
  {
   for (int i = 0; i < j; i++)
   {
    if (Matrix0(i, j) != 0.0)
     count++;
   }
  }
  Ht = vector<T>(count);
  Indexti = vector<int>(count);
  Indextj = vector<int>(count);
  count = 0;
  for (int j = 0; j < dim; j++)
  {
   for (int i = 0; i < j; i++)
   {
    if (Matrix0(i, j) != 0.0)
    {
     Ht[count] = Matrix0(i, j);
     Indexti[count] = i;
     Indextj[count] = j;
     count++;
    }
   }
  }
// - U \sum_{l} n_{l, \spindown} n_{l, \spinup}:
  Matrix0.fillZeroes();
  for (int i = 0; i < N; i++)
  {
   X = N2;
   for (int j = i+1; j < N; j++)
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
  count = 0;
  for (int i = 0; i < dim; i++)
  {
   if (Matrix0(i, i) != 0.0)
    count++;
  }
  HU = vector<T>(count);
  IndexU = vector<int>(count);
  count = 0;
  for (int i = 0; i < dim; i++)
  {
   if (Matrix0(i, i) != 0.0)
   {
    HU[count] = Matrix0(i, i);
    IndexU[count] = i;
    count++;
   }
  }
// - -mu \sum_{l} n_{l}:
  Matrix0.fillZeroes();
  for (int i = 0; i < N; i++)
  {
   X = N1;
   for (int j = i+1; j < N; j++)
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
  count = 0;
  for (int i = 0; i < dim; i++)
  {
   if (Matrix0(i, i) != 0.0)
    count++;
  }
  Hmu = vector<T>(count);
  Indexmu = vector<int>(count);
  count = 0;
  for (int i = 0; i < dim; i++)
  {
   if (Matrix0(i, i) != 0.0)
   {
    Hmu[count] = Matrix0(i, i);
    Indexmu[count] = i;
    count++;
   }
  }
 }
// periodic boundary conditions:
 else if (BC == "periodic")
 {
  cerr << "The following static function is not implemented for periodic BC: " <<
          "template<class T> void HubbardModel<T>::" <<
          "getSparseHubbardMatrix(const string& BC, unsigned int N, T t, T U, T mu, " <<
                                 "vector<T>& Ht, vector<int>& Indexti, vector<int>& Indextj, " <<
                                 "vector<T>& HU, vector<int>& IndexU, " <<
                                 "vector<T>& Hmu, vector<int>& Indexmu)." << endl;
  exit(1);
 }
}

template<class T> void HubbardModel<T>::getFHKHubbard(const string& BC, unsigned int N, T t, T U, T mu,
                                                      unsigned int NP,
                                                      T nPStart, T nPStop, unsigned int R,
                                                      T alphaStart, T x, T eps, unsigned int maxNumRuns,
                                                      Tensor<T>& ErrorAchieved,
                                                      Tensor<unsigned int>& NumRunsDone, Tensor<T>& FHK)
{
#ifdef DEBUG
 if (((BC != "open") && (BC != "periodic")) || (N == 0))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void HubbardModel<T>::" <<
          "getFHKHubbard(const string& BC, unsigned int N, T t, T U, T mu, unsigned int NP, " <<
                        "T nPStart, T nPStop, unsigned int R, " <<
                        "T alphaStart, T x, T eps, unsigned int maxNumRuns, Tensor<T>& ErrorAchieved, " <<
                        "Tensor<unsigned int>& NumRunsDone, Tensor<T>& FHK): " <<
          "(((BC != open) && (BC != periodic)) || (N == 0))." << endl;
  exit(1);
 }
 if ((nPStart < 0.0) || (nPStart >= 2.0) || (nPStop <= nPStart) || (nPStop > 2.0))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void HubbardModel<T>::" <<
          "getFHKHubbard(const string& BC, unsigned int N, T t, T U, T mu, unsigned int NP, " <<
                        "T nPStart, T nPStop, unsigned int R, " <<
                        "T alphaStart, T x, T eps, unsigned int maxNumRuns, Tensor<T>& ErrorAchieved, " <<
                        "Tensor<unsigned int>& NumRunsDone, Tensor<T>& FHK): " <<
          "((nPStart < 0.0) || (nPStart >= 2.0) || (nPStop <= nPStart) || (nPStop > 2.0))." << endl;
  exit(1);
 }
 if ((typeid(T) != typeid(float)) && (typeid(T) != typeid(double)))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void HubbardModel<T>::" <<
          "getFHKHubbard(const string& BC, unsigned int N, T t, T U, T mu, unsigned int NP, " <<
                        "T nPStart, T nPStop, unsigned int R, " <<
                        "T alphaStart, T x, T eps, unsigned int maxNumRuns, Tensor<T>& ErrorAchieved, " <<
                        "Tensor<unsigned int>& NumRunsDone, Tensor<T>& FHK): " <<
          "This function is not implemented for a complex template parameter." << endl;
  exit(1);
 }
 vector<unsigned int> Shape;
 ErrorAchieved.getShape(Shape);
 if (Shape.size() != N)
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void HubbardModel<T>::" <<
          "getFHKHubbard(const string& BC, unsigned int N, T t, T U, T mu, unsigned int NP, " <<
                        "T nPStart, T nPStop, unsigned int R, " <<
                        "T alphaStart, T x, T eps, unsigned int maxNumRuns, Tensor<T>& ErrorAchieved, " <<
                        "Tensor<unsigned int>& NumRunsDone, Tensor<T>& FHK): " <<
          "(ErrorAchieved.Shape.size() != N)." << endl;
  exit(1);
 }
 for (int i = 0; i < Shape.size(); i++)
 {
  if (Shape[i] != R+1)
  {
   cerr << "Program terminated because of error in static function " <<
           "template<class T> void HubbardModel<T>::" <<
           "getFHKHubbard(const string& BC, unsigned int N, T t, T U, T mu, unsigned int NP, " <<
                         "T nPStart, T nPStop, unsigned int R, " <<
                         "T alphaStart, T x, T eps, unsigned int maxNumRuns, Tensor<T>& ErrorAchieved, " <<
                         "Tensor<unsigned int>& NumRunsDone, Tensor<T>& FHK): " <<
           "(ErrorAchieved.Shape[" << i << "] != R+1)." << endl;
   exit(1);
  }
 }
 NumRunsDone.getShape(Shape);
 if (Shape.size() != N)
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void HubbardModel<T>::" <<
          "getFHKHubbard(const string& BC, unsigned int N, T t, T U, T mu, unsigned int NP, " <<
                        "T nPStart, T nPStop, unsigned int R, " <<
                        "T alphaStart, T x, T eps, unsigned int maxNumRuns, Tensor<T>& ErrorAchieved, " <<
                        "Tensor<unsigned int>& NumRunsDone, Tensor<T>& FHK): " <<
          "(NumRunsDone.Shape.size() != N)." << endl;
  exit(1);
 }
 for (int i = 0; i < Shape.size(); i++)
 {
  if (Shape[i] != R+1)
  {
   cerr << "Program terminated because of error in static function " <<
           "template<class T> void HubbardModel<T>::" <<
           "getFHKHubbard(const string& BC, unsigned int N, T t, T U, T mu, unsigned int NP, " <<
                         "T nPStart, T nPStop, unsigned int R, " <<
                         "T alphaStart, T x, T eps, unsigned int maxNumRuns, Tensor<T>& ErrorAchieved, " <<
                         "Tensor<unsigned int>& NumRunsDone, Tensor<T>& FHK): " <<
           "(NumRunsDone.Shape[" << i << "] != R+1)." << endl;
   exit(1);
  }
 }
 FHK.getShape(Shape);
 if (Shape.size() != N)
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void HubbardModel<T>::" <<
          "getFHKHubbard(const string& BC, unsigned int N, T t, T U, T mu, unsigned int NP, " <<
                        "T nPStart, T nPStop, unsigned int R, " <<
                        "T alphaStart, T x, T eps, unsigned int maxNumRuns, Tensor<T>& ErrorAchieved, " <<
                        "Tensor<unsigned int>& NumRunsDone, Tensor<T>& FHK): " <<
          "(FHK.Shape.size() != N)." << endl;
  exit(1);
 }
 for (int i = 0; i < Shape.size(); i++)
 {
  if (Shape[i] != R+1)
  {
   cerr << "Program terminated because of error in static function " <<
           "template<class T> void HubbardModel<T>::" <<
           "getFHKHubbard(const string& BC, unsigned int N, T t, T U, T mu, unsigned int NP, " <<
                         "T nPStart, T nPStop, unsigned int R, " <<
                         "T alphaStart, T x, T eps, unsigned int maxNumRuns, Tensor<T>& ErrorAchieved, " <<
                         "Tensor<unsigned int>& NumRunsDone, Tensor<T>& FHK): " <<
           "(FHK.Shape[" << i << "] != R+1)." << endl;
   exit(1);
  }
 }
#endif
 ErrorAchieved.fillZeroes();
 NumRunsDone.fillZeroes();
 FHK.fillZeroes();
// initialize HubbardModel:
 unsigned int d = 4;
 int n = pow(double(d), double(N)), nev = 1, ncv = min(30, n);
 unsigned int dim = pow(double(d), double(N));
 vector<T> Evals(nev); Matrix<T> Evecs(n, nev);
 vector<T> Parameters(3*N);
 for (int i = 0; i < 3*N; i++)
  Parameters[i] = 1.0;
 HubbardModel<T> HM(BC, N, Parameters);
 MPO<T> nPMPO; vector< Matrix<T> > nPMatrices(N);
 for (int l = 0; l < N; l++)
 {
  HM.getnMPO(l, nPMPO);
  nPMatrices[l] = Matrix<T>(dim, dim);
  nPMPO.getMatrix(nPMatrices[l]);
 }
// initialize sparse Hubbard matrix:
 vector<T> Ht, HU, Hmu, HD, HOD;
 vector<int> Indexti, Indextj, IndexU, Indexmu, IndexD, IndexODi, IndexODj;
 HubbardModel<T>::getSparseHubbardMatrix(BC, N, t, U, mu, Ht, Indexti, Indextj, HU, IndexU, Hmu, Indexmu);
// offdiagonal elements:
 HOD = Ht; IndexODi = Indexti; IndexODj = Indextj;
// diagonal elements:
 vector<T> VectorD(dim), VectorDC(dim);
 fillZeroes(VectorD);
 for (int k = 0; k < IndexU.size(); k++)
  VectorD[IndexU[k]] = HU[k];
 VectorDC = VectorD;
 vector<T> VectorI(4);
 for (int k = 0; k < 4; k++)
  VectorI[k] = 1.0;
 vector<T> VectorN(4);
 VectorN[0] = 0.0; VectorN[1] = 1.0; VectorN[2] = 1.0; VectorN[3] = 2.0;
// potentials:
 vector< vector<T> > VectorsV(N), VectorsVC;
 vector<T> VectorX, VectorY;
 for (int l = 0; l < N; l++)
 {
  VectorX = VectorN;
  for (int m = l; m > 0; m--)
  {
   multiplyDirectProduct(VectorI, VectorX, VectorY);
   VectorX = VectorY;
  }
  for (int m = l; m < N-1; m++)
  {
   multiplyDirectProduct(VectorX, VectorI, VectorY);
   VectorX = VectorY;
  }
  VectorsV[l] = VectorX;
 }
 VectorsVC = VectorsV;
// iterative interacting inversion:
 unsigned int count, run;
 T NPPos, deltanP, alpha, energy, energy0, energy1, distance, distance0, distance1, fhk;
 deltanP = (nPStop-nPStart)/double(R);
 vector<unsigned int> Index(N);
 vector<T> V(N), nP(N), nP0(N), nP1(N), nPDes(N), State(dim);
 for (int l = 0; l < N; l++)
 {
  V[l] = 0.0;
  nP[l] = 0.0;
 }
 for (int pos = 0; pos < FHK.getSize(); pos++)
 {
  FHK.getIndex(pos, Index);
  NPPos = 0.0;
  for (int l = 0; l < N; l++)
   NPPos += nPStart+Index[l]*deltanP;
  if ((NPPos > T(NP)-0.5*deltanP) && (NPPos < T(NP)+0.5*deltanP))
  {
   for (int l = 0; l < N; l++)
    nPDes[l] = nPStart+Index[l]*deltanP;
   alpha = alphaStart;
   distance = 1.0e300; run = 0;
   while ((distance > eps) && (run < maxNumRuns))
   {
// alpha_{0} = x*alpha:
    VectorsV = VectorsVC;
    for (int l = 0; l < N; l++)
     multiply(V[l]+x*alpha*(nP[l]-nPDes[l])-mu, VectorsV[l]);
    VectorD = VectorDC;
    for (int l = 0; l < N; l++)
     add(VectorD, VectorsV[l]);
    count = 0;
    for (int k = 0; k < dim; k++)
    {
     if (VectorD[k] != 0.0)
      count++;
    }
    HD = vector<T>(count);
    IndexD = vector<int>(count);
    count = 0;
    for (int k = 0; k < dim; k++)
    {
     if (VectorD[k] != 0.0)
     {
      HD[count] = VectorD[k];
      IndexD[count] = k;
      count++;
     }
    }
    Matrix<T>::computeLowestEigenstates(HD, IndexD, HOD, IndexODi, IndexODj, n, nev, ncv, Evals, Evecs);
    for (int i = 0; i < dim; i++)
     State[i] = Evecs(i, 0);
    energy0 = Evals[0];
    distance0 = 0.0;
    for (int l = 0; l < N; l++)
    {
     nP0[l] = expectationValueMatrix(State, nPMatrices[l]);
     distance0 += pow(nP0[l]-nPDes[l], 2);
    }
    distance0 = sqrt(abs(distance0));
// alpha_{1} = alpha/x:
    VectorsV = VectorsVC;
    for (int l = 0; l < N; l++)
     multiply(V[l]+alpha*(nP[l]-nPDes[l])/x-mu, VectorsV[l]);
    VectorD = VectorDC;
    for (int l = 0; l < N; l++)
     add(VectorD, VectorsV[l]);
    count = 0;
    for (int k = 0; k < dim; k++)
    {
     if (VectorD[k] != 0.0)
      count++;
    }
    HD = vector<T>(count);
    IndexD = vector<int>(count);
    count = 0;
    for (int k = 0; k < dim; k++)
    {
     if (VectorD[k] != 0.0)
     {
      HD[count] = VectorD[k];
      IndexD[count] = k;
      count++;
     }
    }
    Matrix<T>::computeLowestEigenstates(HD, IndexD, HOD, IndexODi, IndexODj, n, nev, ncv, Evals, Evecs);
    for (int i = 0; i < dim; i++)
     State[i] = Evecs(i, 0);
    energy1 = Evals[0];
    distance1 = 0.0;
    for (int l = 0; l < N; l++)
    {
     nP1[l] = expectationValueMatrix(State, nPMatrices[l]);
     distance1 += pow(nP1[l]-nPDes[l], 2);
    }
    distance1 = sqrt(abs(distance1));
// choose better alpha:
    if (distance0 <= distance1)
    {
     alpha *= x;
     for (int l = 0; l < N; l++)
      V[l] += alpha*(nP[l]-nPDes[l]);
     energy = energy0;
     for (int l = 0; l < N; l++)
      nP[l] = nP0[l];
     distance = distance0;
    }
    else
    {
     alpha /= x;
     for (int l = 0; l < N; l++)
      V[l] += alpha*(nP[l]-nPDes[l]);
     energy = energy1;
     for (int l = 0; l < N; l++)
      nP[l] = nP1[l];
     distance = distance1;
    }
    run++;
   }
   ErrorAchieved.set(pos, distance);
   NumRunsDone.set(pos, run);
   fhk = energy;
   for (int l = 0; l < N; l++)
    fhk -= (V[l]-mu)*nP[l];
   FHK.set(pos, fhk);
  }
 }
}

template<class T> void HubbardModel<T>::getnMPO(unsigned int x, MPO<T>& nMPO) const
{
#ifdef DEBUG
 if ((this->N == 0) || (x > this->N-1))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void HubbardModel<T>::" <<
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

template<class T> void HubbardModel<T>::getNMPO(MPO<T>& NMPO) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void HubbardModel<T>::" <<
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

template<class T> void HubbardModel<T>::getTwoBodyHamiltonian(unsigned int position,
                                                              Matrix<T>& TwoBodyHamiltonian) const
{
 cerr << "The following function is not implemented: " <<
         "template<class T> void HubbardModel<T>::" <<
         "getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian) const." << endl;
 exit(1);
}
