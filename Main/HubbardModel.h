/// Template class HubbardModel implements Fermi-Hubbard model.
/** The template class HubbardModel implements the Fermi-Hubbard model on a double well lattice:
    H = - t0 \sum_{\sigma,l} (c_{\sigma,2l+1}^{+} c_{\sigma,2l} + c_{\sigma,2l}^{+} c_{\sigma,2l+1})
        - t1 \sum_{\sigma,l} (c_{\sigma,2l+2}^{+} c_{\sigma,2l+1} + c_{\sigma,2l+1}^{+} c_{\sigma,2l+2})
        + U0 \sum_{l} n_{\spindown,l} n_{\spinup,l}
        + U1 \sum_{l} n_{2l} n_{2l+1}
        - \mu \sum_{l} n_{l}
        + V_{0} \sum_{l} (l-l_{0})^{2} n_{l}   .
    \param Representation string, the representation, must be "interactions", "MPO" or "matrix"
    \param BC string, the boundary conditions, must be "open" or "periodic"
    \param N unsigned int, the number of spins
    \param d vector<unsigned int>, the physical dimensions of the lattice sites, is fixed to 4
             everywhere
    \param Parameters vector<T>, the parameters {t0, t1, U0, U1, \mu, V_{0}, l_{0}}
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

template<class T> class HubbardModel: public Hamiltonian<T>
{
 typedef T (*PointerToFunction)(double time);

 public:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  HubbardModel();

/// Constructor for time-independent HubbardModel with specific BC, N and Parameters.
/** This constructor initializes a time-independent HubbardModel.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N0 input: unsigned int, the number of lattice sites
    \param Parameters0 input: const vector<T>&, the parameters, must have Parameters0.size()==7 */
  HubbardModel(const string& BC0, unsigned int N0, const vector<T>& Parameters0);

/// Constructor for time-dependent HubbardModel with specific BC, N, TimeFunctions and time.
/** This constructor initializes a time-dependent HubbardModel.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N0 input: unsigned int, the number of lattice sites
    \param TimeFunctions0 input: const vector<PointerToFunction>&, the time functions, must have
                                 TimeFunctions.size()==7
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
                           must fulfill ((Matrix0.getDim0()==4^{N}) && (Matrix0.getDim1()==4^{N})),
                           will be of type "hermitian" */
  void getMatrix(Matrix<T>& Matrix0) const;

/// Returns local particle number operator as MPO.
/** This function returns the local particle number operator n_{x} = n_{\spindown, x} + n_{\spinup, x}
    at position x as a MPO.
    \param x input: unsigned int, the position
    \param nMPO output: MPO<T>&, the local particle number operator as a MPO
    \sa void getNMPO(MPO<T>& NMPO) const */
  void getnMPO(unsigned int x, MPO<T>& nMPO) const;

/// Returns total particle number operator as MPO.
/** This function returns the total particle number operator N = \sum_{l} n_{l} where
    n_{l} = n_{\spindown, l} + n_{\spinup, l} as a MPO.
    \param NMPO output: MPO<T>&, the total particle number operator as a MPO
    \sa void getnMPO(unsigned int x, MPO<T>& nMPO) const */
  void getNMPO(MPO<T>& NMPO) const;

/// Returns local spin operator as MPO.
/** This function returns the local spin operator s_{x} = n_{\spinup, x} - n_{\spindown, x}
    at position x as a MPO.
    \param x input: unsigned int, the position
    \param sMPO output: MPO<T>&, the local spin operator as a MPO
    \sa void getSMPO(MPO<T>& SMPO) const */
  void getsMPO(unsigned int x, MPO<T>& sMPO) const;

/// Returns total spin operator as MPO.
/** This function returns the total spin operator S = \sum_{l} s_{l} where
    s_{l} = n_{\spinup, l} - n_{\spindown, l} as a MPO.
    \param SMPO output: MPO<T>&, the total spin operator as a MPO
    \sa void getsMPO(unsigned int x, MPO<T>& sMPO) const */
  void getSMPO(MPO<T>& SMPO) const;

/// Returns double well state as MPS.
/** This function returns the double well state
    |DW> = \phi^{\otimes N/2}   where   \phi = (|\spinup \spindown> - |\spindown \spinup>) / sqrt(2)
    is a singlet as a MPS. This->N must be even.
    \param D input: unsigned int, the bond dimension of DWMPS
    \param DWMPS output: MPS<T>&, the double well MPS
    \sa void getAFMPS(unsigned int D, MPS<T>& AFMPS) const */
  void getDWMPS(unsigned int D, MPS<T>& DWMPS) const;

/// Returns antiferromagnetic state as MPS.
/** This function returns the antiferromagnetic state
    |AF> = (|\spinup \spindown \spinup ...> + |\spindown \spinup \spindown ...>) / sqrt(2)
    as a MPS.
    \param D input: unsigned int, the bond dimension of AFMPS
    \param AFMPS output: MPS<T>&, the antiferromagnetic MPS
    \sa void getDWMPS(unsigned int D, MPS<T>& DWMPS) const */
  void getAFMPS(unsigned int D, MPS<T>& AFMPS) const;

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
 if (((BC0 != "open") && (BC0 != "periodic")) || (N0 == 0) || (Parameters0.size() != 7))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> HubbardModel<T>::" <<
          "HubbardModel(const string& BC0, unsigned int N0, const vector<T>& Parameters0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (N0 == 0) || (Parameters0.size() != 7))." << endl;
  exit(1);
 }
#endif
 this->Representation = "MPO";
 this->BC = BC0;
 this->N = N0;
 this->d = vector<unsigned int>(this->N);
// in the case of HubbardModel d=4 everywhere:
 for (int i = 0; i < this->N; i++)
 {
  this->d[i] = 4;
 }
 this->Parameters = Parameters0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> HubbardModel<T>::HubbardModel(const string& BC0, unsigned int N0,
                                                const vector<PointerToFunction>& TimeFunctions0,
                                                double time0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (N0 == 0) || (TimeFunctions0.size() != 7))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> HubbardModel<T>::" <<
          "HubbardModel(const string& BC0, unsigned int N0, " <<
          "const vector<PointerToFunction>& TimeFunctions0, double time0): " <<
          "(((BC0 != open) && (BC0 != periodic)) ||  (N0 == 0) || " <<
          "(TimeFunctions0.size() != 7))." << endl;
  exit(1);
 }
#endif
 this->Representation = "MPO";
 this->BC = BC0;
 this->N = N0;
 this->d = vector<unsigned int>(this->N);
// in the case of HubbardModel d=4 everywhere:
 for (int i = 0; i < this->N; i++)
 {
  this->d[i] = 4;
 }
 this->timeDependent = true;
 this->TimeFunctions = TimeFunctions0;
 this->time = time0;
 for (int i = 0; i < this->Parameters.size(); i++)
 {
  this->Parameters[i] = this->TimeFunctions[i](this->time);
 }
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
 cerr << "The following function is not implemented yet: " <<
         "template<class T> void HubbardModel<T>::" <<
         "getInteractions(vector< vector<unsigned int> >& Positions, " <<
                         "vector< vector< Matrix<T> > >& Matrices)." << endl;
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
 if (this->BC == "open")
 {
// the MPO operators O^{i}:
  vector<unsigned int> Shape(3);
  Shape[0] = 11; Shape[1] = 4; Shape[2] = 4;
  Tensor<T> O(Shape);
  O.fillZeroes();
  vector<unsigned int> Index(3);
// O^{0} = identity:
  T element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  O.set(Index, element);
  Index[0] = 0; Index[1] = 1; Index[2] = 1;
  O.set(Index, element);
  Index[0] = 0; Index[1] = 2; Index[2] = 2;
  O.set(Index, element);
  Index[0] = 0; Index[1] = 3; Index[2] = 3;
  O.set(Index, element);
// O^{1} = N:
  element = 1.0;
  Index[0] = 1; Index[1] = 1; Index[2] = 1;
  O.set(Index, element);
  Index[0] = 1; Index[1] = 2; Index[2] = 2;
  O.set(Index, element);
  element = 2.0;
  Index[0] = 1; Index[1] = 3; Index[2] = 3;
  O.set(Index, element);
// O^{2} = N_{\spindown \spinup}:
  element = 1.0;
  Index[0] = 2; Index[1] = 3; Index[2] = 3;
  O.set(Index, element);
// O^{3} = S_{\spinup}^{+}:
  element = 1.0;
  Index[0] = 3; Index[1] = 0; Index[2] = 1;
  O.set(Index, element);
  Index[0] = 3; Index[1] = 2; Index[2] = 3;
  O.set(Index, element);
// O^{4} = \tilde{S}_{\spinup}^{+}:
  element = 1.0;
  Index[0] = 4; Index[1] = 0; Index[2] = 1;
  O.set(Index, element);
  element = -1.0;
  Index[0] = 4; Index[1] = 2; Index[2] = 3;
  O.set(Index, element);
// O^{5} = S_{\spinup}^{-}:
  element = 1.0;
  Index[0] = 5; Index[1] = 1; Index[2] = 0;
  O.set(Index, element);
  Index[0] = 5; Index[1] = 3; Index[2] = 2;
  O.set(Index, element);
// O^{6} = \tilde{S}_{\spinup}^{-}:
  element = 1.0;
  Index[0] = 6; Index[1] = 1; Index[2] = 0;
  O.set(Index, element);
  element = -1.0;
  Index[0] = 6; Index[1] = 3; Index[2] = 2;
  O.set(Index, element);
// O^{7} = S_{\spindown}^{+}:
  element = 1.0;
  Index[0] = 7; Index[1] = 0; Index[2] = 2;
  O.set(Index, element);
  Index[0] = 7; Index[1] = 1; Index[2] = 3;
  O.set(Index, element);
// O^{8} = \tilde{S}_{\spindown}^{+}:
  element = 1.0;
  Index[0] = 8; Index[1] = 0; Index[2] = 2;
  O.set(Index, element);
  element = -1.0;
  Index[0] = 8; Index[1] = 1; Index[2] = 3;
  O.set(Index, element);
// O^{9} = S_{\spindown}^{-}:
  element = 1.0;
  Index[0] = 9; Index[1] = 2; Index[2] = 0;
  O.set(Index, element);
  Index[0] = 9; Index[1] = 3; Index[2] = 1;
  O.set(Index, element);
// O^{10} = \tilde{S}_{\spindown}^{-}:
  element = 1.0;
  Index[0] = 10; Index[1] = 2; Index[2] = 0;
  O.set(Index, element);
  element = -1.0;
  Index[0] = 10; Index[1] = 3; Index[2] = 1;
  O.set(Index, element);
  unsigned int D = 7;
  MPO0 = MPO<T>(this->BC, this->N, this->d, D);
// the parameters:
  T t0 = this->Parameters[0];
  T t1 = this->Parameters[1];
  T U0 = this->Parameters[2];
  T U1 = this->Parameters[3];
  T mu = this->Parameters[4];
  T V0 = this->Parameters[5];
  T l0 = this->Parameters[6];
// A for position 0:
  Shape[0] = 1; Shape[1] = 7; Shape[2] = 11;
  Tensor<T> A(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  element = -mu + V0 * pow(l0, 2);
  Index[0] = 0; Index[1] = 6; Index[2] = 1;
  A.set(Index, element);
  element = U1;
  Index[0] = 0; Index[1] = 1; Index[2] = 1;
  A.set(Index, element);
  element = U0;
  Index[0] = 0; Index[1] = 6; Index[2] = 2;
  A.set(Index, element);
  element = -t0;
  Index[0] = 0; Index[1] = 2; Index[2] = 3;
  A.set(Index, element);
  element = -t0;
  Index[0] = 0; Index[1] = 3; Index[2] = 5;
  A.set(Index, element);
  element = -t0;
  Index[0] = 0; Index[1] = 4; Index[2] = 8;
  A.set(Index, element);
  element = -t0;
  Index[0] = 0; Index[1] = 5; Index[2] = 10;
  A.set(Index, element);
  Tensor<T> OC(O);
  vector<unsigned int> IndexA(1), IndexOC(1);
  IndexA[0] = 2; IndexOC[0] = 0;
  A.contract(IndexA, OC, IndexOC);
  unsigned int position = 0;
  MPO0.set(position, A);
// A for position 1 <= l <= N-2:
  Shape[0] = 7; Shape[1] = 7; Shape[2] = 11;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 6; Index[1] = 6; Index[2] = 0;
  A.set(Index, element);
  element = U0;
  Index[0] = 0; Index[1] = 6; Index[2] = 2;
  A.set(Index, element);
  element = 1.0;
  Index[0] = 3; Index[1] = 6; Index[2] = 4;
  A.set(Index, element);
  Index[0] = 2; Index[1] = 6; Index[2] = 6;
  A.set(Index, element);
  Index[0] = 5; Index[1] = 6; Index[2] = 7;
  A.set(Index, element);
  Index[0] = 4; Index[1] = 6; Index[2] = 9;
  A.set(Index, element);
  Tensor<T> AC(A);
  for (position = 1; position < this->N-1; position++)
  {
   A = AC;
   element = -mu + V0 * pow(T(position) - l0, 2);
   Index[0] = 0; Index[1] = 6; Index[2] = 1;
   A.set(Index, element);
// at even sites use t0 and include U1:
   if ((position % 2) == 0)
   {
    element = U1;
    Index[0] = 0; Index[1] = 1; Index[2] = 1;
    A.set(Index, element);
    element = -t0;
    Index[0] = 0; Index[1] = 2; Index[2] = 3;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 3; Index[2] = 5;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 4; Index[2] = 8;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 5; Index[2] = 10;
    A.set(Index, element);
   }
// at odd sites use t1:
   else if ((position % 2) == 1)
   {
    element = 1.0;
    Index[0] = 1; Index[1] = 6; Index[2] = 1;
    A.set(Index, element);
    element = -t1;
    Index[0] = 0; Index[1] = 2; Index[2] = 3;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 3; Index[2] = 5;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 4; Index[2] = 8;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 5; Index[2] = 10;
    A.set(Index, element);
   }
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   MPO0.set(position, A);
  }
// A for position N-1:
  Shape[0] = 7; Shape[1] = 1; Shape[2] = 11;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 6; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  element = -mu + V0 * pow(T(this->N-1) - l0, 2);
  Index[0] = 0; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  element = 1.0;
  Index[0] = 1; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  element = U0;
  Index[0] = 0; Index[1] = 0; Index[2] = 2;
  A.set(Index, element);
  element = 1.0;
  Index[0] = 3; Index[1] = 0; Index[2] = 4;
  A.set(Index, element);
  Index[0] = 2; Index[1] = 0; Index[2] = 6;
  A.set(Index, element);
  Index[0] = 5; Index[1] = 0; Index[2] = 7;
  A.set(Index, element);
  Index[0] = 4; Index[1] = 0; Index[2] = 9;
  A.set(Index, element);
  OC = O;
  A.contract(IndexA, OC, IndexOC);
  position = this->N-1;
  MPO0.set(position, A);
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
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
 {
  dim *= 4;
 }
 if ((Matrix0.getDim0() != dim) || (Matrix0.getDim1() != dim))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void HubbardModel<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const: " <<
          "((Matrix0.getDim0() != 4^{N}) || (Matrix0.getDim1() != 4^{N}))." << endl;
  exit(1);
 }
#endif
 Matrix0.fillZeroes();
// the interactions:
// - the identity:
 Matrix<T> Identity(4, 4);
 Identity.fillZeroes();
 for (int i = 0; i < 4; i++)
 {
  Identity(i, i) = 1.0;
 }
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
// - S_{\spinup}^{+} =: S0:
 Matrix<T> S0(4, 4);
 S0.fillZeroes();
 S0(1, 0) = 1.0;
 S0(3, 2) = 1.0;
// - \tilde{S}_{\spinup}^{+} =: S1:
 Matrix<T> S1(4, 4);
 S1.fillZeroes();
 S1(1, 0) = 1.0;
 S1(3, 2) = -1.0;
// - S_{\spinup}^{-} =: S2:
 Matrix<T> S2(4, 4);
 S2.fillZeroes();
 S2(0, 1) = 1.0;
 S2(2, 3) = 1.0;
// - \tilde{S}_{\spinup}^{-} =: S3:
 Matrix<T> S3(4, 4);
 S3.fillZeroes();
 S3(0, 1) = 1.0;
 S3(2, 3) = -1.0;
// - S_{\spindown}^{+} =: S4:
 Matrix<T> S4(4, 4);
 S4.fillZeroes();
 S4(2, 0) = 1.0;
 S4(3, 1) = 1.0;
// - \tilde{S}_{\spindown}^{+} =: S5:
 Matrix<T> S5(4, 4);
 S5.fillZeroes();
 S5(2, 0) = 1.0;
 S5(3, 1) = -1.0;
// - S_{\spindown}^{-} =: S6:
 Matrix<T> S6(4, 4);
 S6.fillZeroes();
 S6(0, 2) = 1.0;
 S6(1, 3) = 1.0;
// - \tilde{S}_{\spindown}^{-} =: S7:
 Matrix<T> S7(4, 4);
 S7.fillZeroes();
 S7(0, 2) = 1.0;
 S7(1, 3) = -1.0;
// the parameters:
 T t0 = this->Parameters[0];
 T t1 = this->Parameters[1];
 T U0 = this->Parameters[2];
 T U1 = this->Parameters[3];
 T mu = this->Parameters[4];
 T V0 = this->Parameters[5];
 T l0 = this->Parameters[6];
 T element;
 Matrix<T> X, Y;
// the Hamiltonian sum:
// - -\mu \sum_{l} N_{l}:
 for (int i = 0; i < this->N; i++)
 {
  X = N1;
  for (int j = i; j < this->N-1; j++)
  {
   X.multiplyDirectProduct(Identity, Y);
   X = Y;
  }
  for (int j = i; j > 0; j--)
  {
   Identity.multiplyDirectProduct(X, Y);
   X = Y;
  }
  element = -mu;
  X.multiply(element);
  Matrix0.add(X);
 }
// - V0 \sum_{l} (l-l0)^{2} N_{l}:
 for (int i = 0; i < this->N; i++)
 {
  X = N1;
  for (int j = i; j < this->N-1; j++)
  {
   X.multiplyDirectProduct(Identity, Y);
   X = Y;
  }
  for (int j = i; j > 0; j--)
  {
   Identity.multiplyDirectProduct(X, Y);
   X = Y;
  }
  element = V0 * pow(T(i)-l0, 2);
  X.multiply(element);
  Matrix0.add(X);
 }
// - U0 \sum_{l} N_{\spindown \spinup}:
 for (int i = 0; i < this->N; i++)
 {
  X = N2;
  for (int j = i; j < this->N-1; j++)
  {
   X.multiplyDirectProduct(Identity, Y);
   X = Y;
  }
  for (int j = i; j > 0; j--)
  {
   Identity.multiplyDirectProduct(X, Y);
   X = Y;
  }
  element = U0;
  X.multiply(element);
  Matrix0.add(X);
 }
// - U1 \sum_{l} N_{2l} \otimes N_{2l+1}:
// N3 := N_{2l} \otimes N_{2l+1}:
 Matrix<T> N3;
 N1.multiplyDirectProduct(N1, N3);
 for (int i = 0; i < this->N-1; i+=2)
 {
  X = N3;
  for (int j = i; j < this->N-2; j++)
  {
   X.multiplyDirectProduct(Identity, Y);
   X = Y;
  }
  for (int j = i; j > 0; j--)
  {
   Identity.multiplyDirectProduct(X, Y);
   X = Y;
  }
  element = U1;
  X.multiply(element);
  Matrix0.add(X);
 }
// - -t0 \sum_{l} S_{2l}
// S := \tilde{S}_{\spindown}^{-} \otimes S_{\spindown}^{+} +
//      \tilde{S}_{\spindown}^{+} \otimes S_{\spindown}^{-} +
//      S_{\spinup}^{-} \otimes \tilde{S}_{\spinup}^{+} +
//      S_{\spinup}^{+} \otimes \tilde{S}_{\spinup}^{-}:
 Matrix<T> S(16, 16);
 S.fillZeroes();
 S7.multiplyDirectProduct(S4, X);
 S.add(X);
 S5.multiplyDirectProduct(S6, X);
 S.add(X);
 S2.multiplyDirectProduct(S1, X);
 S.add(X);
 S0.multiplyDirectProduct(S3, X);
 S.add(X);
 for (int i = 0; i < this->N-1; i+=2)
 {
  X = S;
  for (int j = i; j < this->N-2; j++)
  {
   X.multiplyDirectProduct(Identity, Y);
   X = Y;
  }
  for (int j = i; j > 0; j--)
  {
   Identity.multiplyDirectProduct(X, Y);
   X = Y;
  }
  element = -t0;
  X.multiply(element);
  Matrix0.add(X);
 }
// - -t1 \sum_{l} S_{2l+1}
 for (int i = 1; i < this->N-1; i+=2)
 {
  X = S;
  for (int j = i; j < this->N-2; j++)
  {
   X.multiplyDirectProduct(Identity, Y);
   X = Y;
  }
  for (int j = i; j > 0; j--)
  {
   Identity.multiplyDirectProduct(X, Y);
   X = Y;
  }
  element = -t1;
  X.multiply(element);
  Matrix0.add(X);
 }
 if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void HubbardModel<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const." << endl;
  exit(1);
 }
 Matrix0.setType("hermitian");
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
// the particle number operator O = n_{l} = n_{\spindown, l} + n_{\spinup, l}:
 Matrix<T> O(4, 4);
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
// the particle number operator O = n_{l} = n_{\spindown, l} + n_{\spinup, l}:
 Matrix<T> O(4, 4);
 O.fillZeroes();
 O(1, 1) = 1.0;
 O(2, 2) = 1.0;
 O(3, 3) = 2.0;
 getMPOFromSumLocalOperator(this->BC, this->N, O, NMPO);
}

template<class T> void HubbardModel<T>::getsMPO(unsigned int x, MPO<T>& sMPO) const
{
#ifdef DEBUG
 if ((this->N == 0) || (x > this->N-1))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void HubbardModel<T>::" <<
          "getsMPO(unsigned int x, MPO<T>& sMPO) const: " <<
          "((this->N == 0) || (x > this->N-1))." << endl;
  exit(1);
 }
#endif
// the spin operator O = s_{l} = n_{\spinup, l} - n_{\spindown, l}:
 Matrix<T> O(4, 4);
 O.fillZeroes();
 O(1, 1) = 1.0;
 O(2, 2) = -1.0;
 getMPOFromLocalOperator(this->BC, this->N, x, O, sMPO);
}

template<class T> void HubbardModel<T>::getSMPO(MPO<T>& SMPO) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void HubbardModel<T>::" <<
          "getSMPO(MPO<T>& SMPO) const: " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
// the spin operator O = s_{l} = n_{\spinup, l} - n_{\spindown, l}:
 Matrix<T> O(4, 4);
 O.fillZeroes();
 O(1, 1) = 1.0;
 O(2, 2) = -1.0;
 getMPOFromSumLocalOperator(this->BC, this->N, O, SMPO);
}

template<class T> void HubbardModel<T>::getDWMPS(unsigned int D, MPS<T>& DWMPS) const
{
#ifdef DEBUG
 if ((this->N == 0) || ((this->N % 2) != 0) || (D < 2))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void tJModel<T>::" <<
          "getDWMPS(unsigned int D, MPS<T>& DWMPS) const: " <<
          "((this->N == 0) || ((this->N % 2) != 0) || (D < 2))." << endl;
  exit(1);
 }
#endif
 DWMPS = MPS<T>(this->BC, this->N, this->d, D);
 vector<unsigned int> Shape(3), Index(3);
 T element;
 unsigned int position;
 if (this->BC == "open")
 {
// A for position 0:
  position = 0;
  DWMPS.getOpenBCShape(position, Shape);
  Tensor<T> A(Shape);
  A.fillZeroes();
  element = 1.0 / sqrt(2.0);
  Index[0] = 0; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  element = -1.0 / sqrt(2.0);
  Index[0] = 0; Index[1] = 1; Index[2] = 2;
  A.set(Index, element);
  DWMPS.set(position, A);
// A for odd position 1 <= l <= N-2:
  element = 1.0;
  for (position = 1; position < this->N-1; position += 2)
  {
   DWMPS.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   Index[0] = 1; Index[1] = 0; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 0; Index[2] = 2;
   A.set(Index, element);
   DWMPS.set(position, A);
  }
// A for even position 1 <= l <= N-2:
  for (position = 2; position < this->N-1; position += 2)
  {
   DWMPS.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0 / sqrt(2.0);
   Index[0] = 0; Index[1] = 0; Index[2] = 1;
   A.set(Index, element);
   element = -1.0 / sqrt(2.0);
   Index[0] = 0; Index[1] = 1; Index[2] = 2;
   A.set(Index, element);
   DWMPS.set(position, A);
  }
// A for position N-1:
  position = this->N-1;
  DWMPS.getOpenBCShape(position, Shape);
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 1; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  Index[0] = 0; Index[1] = 0; Index[2] = 2;
  A.set(Index, element);
  DWMPS.set(position, A);
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void tJModel<T>::" <<
          "getDWMPS(unsigned int D, MPS<T>& DWMPS) const." << endl;
  exit(1);
 }
}

template<class T> void HubbardModel<T>::getAFMPS(unsigned int D, MPS<T>& AFMPS) const
{
#ifdef DEBUG
 if ((this->N == 0) || (D < 2))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void HubbardModel<T>::" <<
          "getAFMPS(unsigned int D, MPS<T>& AFMPS) const: " <<
          "((this->N == 0) || (D < 2))." << endl;
  exit(1);
 }
#endif
 AFMPS = MPS<T>(this->BC, this->N, this->d, D);
 vector<unsigned int> Shape(3), Index(3);
 T element;
 unsigned int position;
 if (this->BC == "open")
 {
// A for position 0:
  position = 0;
  AFMPS.getOpenBCShape(position, Shape);
  Tensor<T> A(Shape);
  A.fillZeroes();
  element = 1.0 / sqrt(2.0);
  Index[0] = 0; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  Index[0] = 0; Index[1] = 1; Index[2] = 2;
  A.set(Index, element);
  AFMPS.set(position, A);
// A for position 1 <= l <= N-2:
  element = 1.0;
  for (position = 1; position < this->N-1; position++)
  {
   AFMPS.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   Index[0] = 1; Index[1] = 0; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 1; Index[2] = 2;
   A.set(Index, element);
   AFMPS.set(position, A);
  }
// A for position N-1:
  position = this->N-1;
  AFMPS.getOpenBCShape(position, Shape);
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 1; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  Index[0] = 0; Index[1] = 0; Index[2] = 2;
  A.set(Index, element);
  AFMPS.set(position, A);
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void HubbardModel<T>::" <<
          "getAFMPS(unsigned int D, MPS<T>& AFMPS) const." << endl;
  exit(1);
 }
}

template<class T> void HubbardModel<T>::getTwoBodyHamiltonian(unsigned int position,
                                                              Matrix<T>& TwoBodyHamiltonian) const
{
#ifdef DEBUG
 if ((position >= this->N-1) || (TwoBodyHamiltonian.getDim0() != 16) ||
     (TwoBodyHamiltonian.getDim1() != 16))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void HubbardModel<T>::" <<
          "getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian) const: " <<
          "((position >= this->N-1) || (TwoBodyHamiltonian.getDim0() != 16) || " <<
          "(TwoBodyHamiltonian.getDim1() != 16))." << endl;
  exit(1);
 }
#endif
 TwoBodyHamiltonian.fillZeroes();
// the interactions:
// - the identity:
 Matrix<T> Identity(4, 4);
 Identity.fillZeroes();
 for (int i = 0; i < 4; i++)
 {
  Identity(i, i) = 1.0;
 }
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
// - S_{\spinup}^{+} =: S0:
 Matrix<T> S0(4, 4);
 S0.fillZeroes();
 S0(1, 0) = 1.0;
 S0(3, 2) = 1.0;
// - \tilde{S}_{\spinup}^{+} =: S1:
 Matrix<T> S1(4, 4);
 S1.fillZeroes();
 S1(1, 0) = 1.0;
 S1(3, 2) = -1.0;
// - S_{\spinup}^{-} =: S2:
 Matrix<T> S2(4, 4);
 S2.fillZeroes();
 S2(0, 1) = 1.0;
 S2(2, 3) = 1.0;
// - \tilde{S}_{\spinup}^{-} =: S3:
 Matrix<T> S3(4, 4);
 S3.fillZeroes();
 S3(0, 1) = 1.0;
 S3(2, 3) = -1.0;
// - S_{\spindown}^{+} =: S4:
 Matrix<T> S4(4, 4);
 S4.fillZeroes();
 S4(2, 0) = 1.0;
 S4(3, 1) = 1.0;
// - \tilde{S}_{\spindown}^{+} =: S5:
 Matrix<T> S5(4, 4);
 S5.fillZeroes();
 S5(2, 0) = 1.0;
 S5(3, 1) = -1.0;
// - S_{\spindown}^{-} =: S6:
 Matrix<T> S6(4, 4);
 S6.fillZeroes();
 S6(0, 2) = 1.0;
 S6(1, 3) = 1.0;
// - \tilde{S}_{\spindown}^{-} =: S7:
 Matrix<T> S7(4, 4);
 S7.fillZeroes();
 S7(0, 2) = 1.0;
 S7(1, 3) = -1.0;
// S := \tilde{S}_{\spindown}^{-} \otimes S_{\spindown}^{+} +
//      \tilde{S}_{\spindown}^{+} \otimes S_{\spindown}^{-} +
//      S_{\spinup}^{-} \otimes \tilde{S}_{\spinup}^{+} +
//      S_{\spinup}^{+} \otimes \tilde{S}_{\spinup}^{-}:
 Matrix<T> S(16, 16), X;
 S.fillZeroes();
 S7.multiplyDirectProduct(S4, X);
 S.add(X);
 S5.multiplyDirectProduct(S6, X);
 S.add(X);
 S2.multiplyDirectProduct(S1, X);
 S.add(X);
 S0.multiplyDirectProduct(S3, X);
 S.add(X);
// the parameters:
 T t0 = this->Parameters[0];
 T t1 = this->Parameters[1];
 T U0 = this->Parameters[2];
 T U1 = this->Parameters[3];
 T mu = this->Parameters[4];
 T V0 = this->Parameters[5];
 T l0 = this->Parameters[6];
 T element;
 if (this->BC == "open")
 {
  if (position == 0)
  {
   N1.multiplyDirectProduct(Identity, X);
   element = -mu;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N1, X);
   element = -0.5 * mu;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   N1.multiplyDirectProduct(Identity, X);
   element = V0 * pow(l0, 2);
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N1, X);
   element = 0.5 * V0 * pow(l0-1.0, 2);
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   N1.multiplyDirectProduct(N1, X);
   element = U1;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   N2.multiplyDirectProduct(Identity, X);
   element = U0;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N2, X);
   element = 0.5 * U0;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   X = S;
   element = -t0;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
  }
  else if ((position >= 1) && (position <= this->N-3))
  {
   N1.multiplyDirectProduct(Identity, X);
   element = -0.5 * mu;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N1, X);
   element = -0.5 * mu;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   N1.multiplyDirectProduct(Identity, X);
   element = 0.5 * V0 * pow(T(position)-l0, 2);
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N1, X);
   element = 0.5 * V0 * pow(T(position+1)-l0, 2);
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   N2.multiplyDirectProduct(Identity, X);
   element = 0.5 * U0;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N2, X);
   element = 0.5 * U0;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   if ((position % 2) == 0)
   {
    X = S;
    element = -t0;
    X.multiply(element);
    TwoBodyHamiltonian.add(X);
    N1.multiplyDirectProduct(N1, X);
    element = U1;
    X.multiply(element);
    TwoBodyHamiltonian.add(X);
   }
   else if ((position % 2) == 1)
   {
    X = S;
    element = -t1;
    X.multiply(element);
    TwoBodyHamiltonian.add(X);
   }
  }
  else if (position == this->N-2)
  {
   N1.multiplyDirectProduct(Identity, X);
   element = -0.5 * mu;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N1, X);
   element = -mu;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   N1.multiplyDirectProduct(Identity, X);
   element = 0.5 * V0 * pow(T(this->N-2)-l0, 2);
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N1, X);
   element = V0 * pow(T(this->N-1)-l0, 2);
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   N2.multiplyDirectProduct(Identity, X);
   element = 0.5 * U0;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N2, X);
   element = U0;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   if (((this->N-2) % 2) == 0)
   {
    X = S;
    element = -t0;
    X.multiply(element);
    TwoBodyHamiltonian.add(X);
    N1.multiplyDirectProduct(N1, X);
    element = U1;
    X.multiply(element);
    TwoBodyHamiltonian.add(X);
   }
   else if (((this->N-2) % 2) == 1)
   {
    X = S;
    element = -t1;
    X.multiply(element);
    TwoBodyHamiltonian.add(X);
   }
  }
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void HubbardModel<T>::" <<
          "getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian) const." << endl;
  exit(1);
 }
}
