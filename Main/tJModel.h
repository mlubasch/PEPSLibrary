/// Template class tJModel implements t-J model.
/** The template class tJModel implements the t-J model on a double well lattice:
    H = - t0 \sum_{\sigma,l} (c_{\sigma,2l+1}^{+} c_{\sigma,2l} + c_{\sigma,2l}^{+} c_{\sigma,2l+1})
        + J0 \sum_{l} (0.5 * (S_{2l}^{+} S_{2l+1}^{-} + S_{2l}^{-} S_{2l+1}^{+}) +
                       + S_{2l}^{z} S_{2l+1}^{z} - 0.25 n_{2l} n_{2l+1})
        - t1 \sum_{\sigma,l} (c_{\sigma,2l+2}^{+} c_{\sigma,2l+1} + c_{\sigma,2l+1}^{+} c_{\sigma,2l+2})
        + J1 \sum_{l} (0.5 * (S_{2l+1}^{+} S_{2l+2}^{-} + S_{2l+1}^{-} S_{2l+2}^{+}) +
                       + S_{2l+1}^{z} S_{2l+2}^{z} - 0.25 n_{2l+1} n_{2l+2})
        + U \sum_{l} n_{2l} n_{2l+1}
        - \mu \sum_{l} n_{l}
        + V_{0} \sum_{l} (l-l_{0})^{2} n_{l}   .
    We note that this Hamiltonian is projected onto the subspace of only one fermion per lattice site.
    \param Representation string, the representation, must be "interactions", "MPO" or "matrix"
    \param BC string, the boundary conditions, must be "open" or "periodic"
    \param N unsigned int, the number of spins
    \param d vector<unsigned int>, the physical dimensions of the lattice sites, is fixed to 3
             everywhere
    \param Parameters vector<T>, the parameters {t0, J0, t1, J1, U, \mu, V_{0}, l_{0}}
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

template<class T> class tJModel: public Hamiltonian<T>
{
 typedef T (*PointerToFunction)(double time);

 public:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  tJModel();

/// Constructor for time-independent tJModel with specific BC, N and Parameters.
/** This constructor initializes a time-independent tJModel.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N0 input: unsigned int, the number of lattice sites
    \param Parameters0 input: const vector<T>&, the parameters, must have Parameters0.size()==8 */
  tJModel(const string& BC0, unsigned int N0, const vector<T>& Parameters0);

/// Constructor for time-dependent tJModel with specific BC, N, TimeFunctions and time.
/** This constructor initializes a time-dependent tJModel.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N0 input: unsigned int, the number of lattice sites
    \param TimeFunctions0 input: const vector<PointerToFunction>&, the time functions, must have
                                 TimeFunctions.size()==8
    \param time0 input: double, the time */
  tJModel(const string& BC0, unsigned int N0, const vector<PointerToFunction>& TimeFunctions0,
          double time0);

/// Standard copy constructor.
/** The standard copy constructor copies the input tJModel into this.
    \param tJModel0 input: const tJModel<T>&, to be copied into this
    \sa tJModel<T>& operator=(const tJModel<T>& tJModel0) */
  tJModel(const tJModel<T>& tJModel0);

/// Standard destructor.
/** The standard destructor deletes the elements of tJModel. */
  ~tJModel();

/// Assigns tJModel to this.
/** The operator= allows to assign tJModel0 to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side tJModel0.
    \param tJModel0 input: const tJModel<T>&, to be copied into this
    \return tJModel<T>&, a reference to the new this
    \sa tJModel(const tJModel<T>& tJModel0) */
  tJModel<T>& operator=(const tJModel<T>& tJModel0);

/// Returns interactions representation of tJModel.
/** This function returns the interactions of this tJModel. The interactions
    are implemented as a vector of positions and a vector of matrices where each vector entry
    corresponds to one term of the sum making up the Hamiltonian.
    \param Positions output: vector< vector<unsigned int> >&, the positions of the matrices
    \param Matrices output: vector< vector< Matrix<T> > >&, the matrices */
  void getInteractions(vector< vector<unsigned int> >& Positions,
                       vector< vector< Matrix<T> > >& Matrices) const;

/// Returns MPO representation of tJModel.
/** This function returns this tJModel as a MPO.
    \param MPO0 output: MPO<T>&, the MPO representing this tJModel */
  void getMPO(MPO<T>& MPO0) const;

/// Returns matrix representation of tJModel.
/** This function returns this tJModel as a matrix.
    \param Matrix0 output: Matrix<T>&, the matrix representing this tJModel,
                           must fulfill ((Matrix0.getDim0()==3^{N}) && (Matrix0.getDim1()==3^{N})),
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

/// Returns squared staggered magnetization operator as MPO.
/** This function returns the squared staggered magnetization operator
    M^{2} = (\sum_{l} (-1)^{l} \vec{S}_{l})^{2} ,
    where the sum over l runs from leftEnd to rightEnd, as a MPO.
    \param leftEnd input: unsigned int, the left end position of this operator
    \param rightEnd input: unsigned int, the right end position of this operator
    \param M2MPO output: MPO<T>&, the squared staggered magnetization operator as a MPO */
  void getM2MPO(unsigned int leftEnd, unsigned int rightEnd, MPO<T>& M2MPO) const;

/// Returns double well Heisenberg model as MPO.
/** This function returns the Heisenberg model on a double well lattice
    H_{DWHM} = J0 \sum_{l} (0.5 * (S_{2l}^{+} S_{2l+1}^{-} + S_{2l}^{-} S_{2l+1}^{+}) +
                            + S_{2l}^{z} S_{2l+1}^{z}) +
             + J1 \sum_{l} (0.5 * (S_{2l+1}^{+} S_{2l+2}^{-} + S_{2l+1}^{-} S_{2l+2}^{+}) +
                            + S_{2l+1}^{z} S_{2l+2}^{z}),
    where the sum over l runs from leftEnd to rightEnd, as a MPO.
    \param leftEnd input: unsigned int, the left end position of this operator
    \param rightEnd input: unsigned int, the right end position of this operator
    \param DWHMMPO output: MPO<T>&, the double well Heisenberg model as a MPO */
  void getDWHMMPO(unsigned int leftEnd, unsigned int rightEnd, MPO<T>& DWHMMPO) const;

/// Returns local hole number operator as MPO.
/** This function returns the local hole number operator n_{x} = n_{\hole, x} at position x as a MPO.
    \param x input: unsigned int, the position
    \param holeMPO output: MPO<T>&, the local hole number operator as a MPO
    \sa void getHoleMPO(MPO<T>& HoleMPO) const */
  void getholeMPO(unsigned int x, MPO<T>& holeMPO) const;

/// Returns total hole number operator as MPO.
/** This function returns the total hole number operator N = \sum_{l} n_{l} where
    n_{l} = n_{\hole, l} as a MPO.
    \param HoleMPO output: MPO<T>&, the total hole number operator as a MPO
    \sa void getholeMPO(unsigned int x, MPO<T>& holeMPO) const */
  void getHoleMPO(MPO<T>& HoleMPO) const;

/// Sets hole at position.
/** This function sets a hole at position in MPS0 by applying
    c_{\spinup, position}+c_{\spindown, position}. MPS0 is not normalized afterwards.
    \param position input: unsigned int, the position of the hole
    \param MPS0 input/output: MPS<T>&, the MPS */
  static void setHole(unsigned int position, MPS<T>& MPS0);

/// Returns single particle state.
/** This static function returns the MPS MPS0 describing a single particle in the state
    |\phi> = Coefficients[0]|0> + Coefficients[1]|1> + Coefficients[2]|2>
    located at position.
    Boundary conditions BC, length N and virtual dimension D are given by MPS0.
    MPS0 must have physical dimension d==3 everywhere.
    \param position input: unsigned int, the position
    \param Coefficients input: const vector<T>&, the coefficients of the single particle state,
                               must fulfill Coefficients.size()==3
    \param MPS0 input/output: MPS<T>&, on input its shape specifies the shape of the single particle state and it
                              must fulfill d==3 everywhere, on output it will be the single particle state */
  static void getSingleParticle(unsigned int position, const vector<T>& Coefficients, MPS<T>& MPS0);

/// Returns double well single particle eigenstate.
/** This static function returns the two site MPS MPS0 describing the single particle eigenstates of a double well:
    |\phi_{GS}> = (|x 0> + |0 x>)/sqrt(2) with energy -t
    and
    |\phi_{ES}> = (|x 0> - |0 x>)/sqrt(2) with energy t
    where
    |x> = Cup |\spinup> + Cdown |\spindown> .
    \param WhichEigenstate input: const string&, must be GS for ground state and ES for excited state
    \param Cup input: T, the coefficient of the spin-up component
    \param Cdown input: T, the coefficient of the spin-down component
    \param MPS0 output: MPS<T>&, the resulting double well single particle eigenstate */
  static void getDWSingleParticle(const string& WhichEigenstate, T Cup, T Cdown, MPS<T>& MPS0);

/// Returns dimer state.
/** This static function returns the dimer state |Dimer> = |\phi>^{\otimes N/2}, where
    |\phi> = (|\spinup \spindown> - |\spindown \spinup>) / sqrt(2) is a singlet, in MPS0. Holes are set
    according to HolePositions and a single particle is put into the state
    |\chi> = (|\spinup> + |\spindown>) / sqrt(2). Boundary conditions BC, length N and virtual
    dimension D of the dimer are given by MPS0. MPS0 must have physical dimension d==3 everywhere.
    \param HolePositions input: const vector<unsigned int>&, the hole positions
    \param MPS0 input/output: MPS<T>&, on input its shape specifies the shape of the dimer state and it
                              must fulfill d==3 everywhere, on output it will be the dimer state with
                              holes */
  static void getDimer(const vector<unsigned int>& HolePositions, MPS<T>& MPS0);

/// Returns Néel state.
/** This static function returns the Néel state
    |Néel> = (|\spinup \spindown \spinup ...> + |\spindown \spinup \spindown ...>) / sqrt(2)
    in MPS0. Boundary conditions BC, length N and virtual dimension D are given by MPS0. MPS0 must have
    physical dimension d==3 everywhere.
    \param MPS0 input/output: MPS<T>&, on input its shape specifies the shape of the Néel state and it
                              must fulfill d==3 everywhere, on output it will be the Néel state */
  static void getNeel(MPS<T>& MPS0);

/// Sets holes at left and right boundary.
/** This static function sets numHolesLeft holes left of MPS0 and numHolesRight right of MPS0 and
    stores the resulting state in MPS1.
    \param MPS0 input: const MPS<T>&, the initial MPS
    \param numHolesLeft input: unsigned int, the number of holes to be set left of MPS0
    \param numHolesRight input: unsigned int, the number of holes to be set right of MPS0
    \param MPS1 output: MPS<T>&, the resulting MPS */
  static void setHoles(const MPS<T>& MPS0, unsigned int numHolesLeft, unsigned int numHolesRight,
                       MPS<T>& MPS1);

/// Returns spin-spin correlation function.
/** This static function returns the spin-spin correlation function
    \vec{S}_{position0} \cdot \vec{S}_{position1}
    in MPO0. Boundary conditions BC and length N are given by MPO0.
    \param position0 input: unsigned int, the position of the left spin operator
    \param position1 input: unsigned int, the position of the right spin operator
    \param MPO0 input/output: MPO<T>&, on input its shape specifies the shape of the spin-spin correlator,
                              on output it will be the spin-spin correlator */
  static void getSpinSpinCorrelator(unsigned int position0, unsigned int position1, MPO<T>& MPO0);

/// Computes lowest lying eigenstates of double well Heisenberg model.
/** This static function computes the nev lowest lying eigenstates of the double well Heisenberg model
    H_{DWHM} = J0 \sum_{l} (0.5 * (S_{2l}^{+} S_{2l+1}^{-} + S_{2l}^{-} S_{2l+1}^{+}) +
                            + S_{2l}^{z} S_{2l+1}^{z}) +
             + J1 \sum_{l} (0.5 * (S_{2l+1}^{+} S_{2l+2}^{-} + S_{2l+1}^{-} S_{2l+2}^{+}) +
                            + S_{2l+1}^{z} S_{2l+2}^{z}),
    for spin-1/2 particles. T must be real.
    \param BC input: const string&, the boundary conditions, must be either "open" or "periodic"
    \param N input: unsigned int, the length of the chain, must be even
    \param J0 input: T, parameter
    \param J1 input: T, parameter
    \param nev input: unsigned int, the number of eigenvalues and eigenstates requested
    \param Evals output: vector<T>&, the resulting eigenvalues
    \param Evecs output: Matrix<T>&, the resulting eigenvectors in the columns */
  static void computeDWHMEigenstates(const string& BC, unsigned int N, T J0, T J1, unsigned int nev,
                                     vector<T>& Evals, Matrix<T>& Evecs);

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

template<class T> tJModel<T>::tJModel()
{
 this->N = 0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> tJModel<T>::tJModel(const string& BC0, unsigned int N0, const vector<T>& Parameters0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (N0 == 0) || (Parameters0.size() != 8))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> tJModel<T>::" <<
          "tJModel(const string& BC0, unsigned int N0, const vector<T>& Parameters0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (N0 == 0) || (Parameters0.size() != 8))." << endl;
  exit(1);
 }
#endif
 this->Representation = "MPO";
 this->BC = BC0;
 this->N = N0;
 this->d = vector<unsigned int>(this->N);
// in the case of tJModel d=3 everywhere:
 for (int i = 0; i < this->N; i++)
 {
  this->d[i] = 3;
 }
 this->Parameters = Parameters0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> tJModel<T>::tJModel(const string& BC0, unsigned int N0,
                                      const vector<PointerToFunction>& TimeFunctions0, double time0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (N0 == 0) || (TimeFunctions0.size() != 8))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> tJModel<T>::" <<
          "tJModel(const string& BC0, unsigned int N0, " <<
          "const vector<PointerToFunction>& TimeFunctions0, double time0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (N0 == 0) || " <<
          "(TimeFunctions0.size() != 8))." << endl;
  exit(1);
 }
#endif
 this->Representation = "MPO";
 this->BC = BC0;
 this->N = N0;
 this->d = vector<unsigned int>(this->N);
// in the case of tJModel d=3 everywhere:
 for (int i = 0; i < this->N; i++)
 {
  this->d[i] = 3;
 }
 this->Parameters = vector<T>(8);
 this->timeDependent = true;
 this->TimeFunctions = TimeFunctions0;
 this->time = time0;
 for (int i = 0; i < this->Parameters.size(); i++)
 {
  this->Parameters[i] = this->TimeFunctions[i](this->time);
 }
}

template<class T> tJModel<T>::tJModel(const tJModel<T>& tJModel0)
{
 this->Representation = tJModel0.Representation;
 this->BC = tJModel0.BC;
 this->N = tJModel0.N;
 this->d = tJModel0.d;
 this->Parameters = tJModel0.Parameters;
 this->timeDependent = tJModel0.timeDependent;
 this->TimeFunctions = tJModel0.TimeFunctions;
 this->time = tJModel0.time;
}

template<class T> tJModel<T>::~tJModel() {}

template<class T> tJModel<T>& tJModel<T>::operator=(const tJModel<T>& tJModel0)
{
 if (this != &tJModel0)
 {
  this->Representation = tJModel0.Representation;
  this->BC = tJModel0.BC;
  this->N = tJModel0.N;
  this->d = tJModel0.d;
  this->Parameters = tJModel0.Parameters;
  this->timeDependent = tJModel0.timeDependent;
  this->TimeFunctions = tJModel0.TimeFunctions;
  this->time = tJModel0.time;
 }
 return *this;
}

template<class T> void tJModel<T>::getInteractions(vector< vector<unsigned int> >& Positions,
                                                   vector< vector< Matrix<T> > >& Matrices) const
{
 cerr << "The following function is not implemented yet: " <<
         "template<class T> void tJModel<T>::" <<
         "getInteractions(vector< vector<unsigned int> >& Positions, " <<
                         "vector< vector< Matrix<T> > >& Matrices)." << endl;
 exit(1);
}

template<class T> void tJModel<T>::getMPO(MPO<T>& MPO0) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void tJModel<T>::" <<
          "getMPO(MPO<T>& MPO0) const: " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
 if (this->BC == "open")
 {
// the MPO operators O^{i}:
  vector<unsigned int> Shape(3);
  Shape[0] = 9; Shape[1] = 3; Shape[2] = 3;
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
// O^{1} = N:
  element = 1.0;
  Index[0] = 1; Index[1] = 1; Index[2] = 1;
  O.set(Index, element);
  Index[0] = 1; Index[1] = 2; Index[2] = 2;
  O.set(Index, element);
// O^{2} = S^{+}:
  element = 1.0;
  Index[0] = 2; Index[1] = 2; Index[2] = 1;
  O.set(Index, element);
// O^{3} = S^{-}:
  element = 1.0;
  Index[0] = 3; Index[1] = 1; Index[2] = 2;
  O.set(Index, element);
// O^{4} = S^{z}:
  element = 0.5;
  Index[0] = 4; Index[1] = 1; Index[2] = 1;
  O.set(Index, element);
  element = -0.5;
  Index[0] = 4; Index[1] = 2; Index[2] = 2;
  O.set(Index, element);
// O^{5} = C_{\spinup}^{+}:
  element = 1.0;
  Index[0] = 5; Index[1] = 0; Index[2] = 1;
  O.set(Index, element);
// O^{6} = C_{\spinup}^{-}:
  element = 1.0;
  Index[0] = 6; Index[1] = 1; Index[2] = 0;
  O.set(Index, element);
// O^{7} = C_{\spindown}^{+}:
  element = 1.0;
  Index[0] = 7; Index[1] = 0; Index[2] = 2;
  O.set(Index, element);
// O^{8} = C_{\spindown}^{-}:
  element = 1.0;
  Index[0] = 8; Index[1] = 2; Index[2] = 0;
  O.set(Index, element);
  unsigned int D = 10;
  MPO0 = MPO<T>(this->BC, this->N, this->d, D);
// the parameters:
  T t0 = this->Parameters[0];
  T J0 = this->Parameters[1];
  T t1 = this->Parameters[2];
  T J1 = this->Parameters[3];
  T U = this->Parameters[4];
  T mu = this->Parameters[5];
  T V0 = this->Parameters[6];
  T l0 = this->Parameters[7];
// A for position 0:
  Shape[0] = 1; Shape[1] = D; Shape[2] = 9;
  Tensor<T> A(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  element = V0 * pow(l0, 2) - mu;
  Index[0] = 0; Index[1] = 9; Index[2] = 1;
  A.set(Index, element);
  element = U - T(0.25) * J0;
  Index[0] = 0; Index[1] = 1; Index[2] = 1;
  A.set(Index, element);
  element = T(0.5) * J0;
  Index[0] = 0; Index[1] = 2; Index[2] = 2;
  A.set(Index, element);
  element = T(0.5) * J0;
  Index[0] = 0; Index[1] = 3; Index[2] = 3;
  A.set(Index, element);
  element = J0;
  Index[0] = 0; Index[1] = 4; Index[2] = 4;
  A.set(Index, element);
  element = -t0;
  Index[0] = 0; Index[1] = 5; Index[2] = 5;
  A.set(Index, element);
  Index[0] = 0; Index[1] = 6; Index[2] = 6;
  A.set(Index, element);
  Index[0] = 0; Index[1] = 7; Index[2] = 7;
  A.set(Index, element);
  Index[0] = 0; Index[1] = 8; Index[2] = 8;
  A.set(Index, element);
  Tensor<T> OC(O);
  vector<unsigned int> IndexA(1), IndexOC(1);
  IndexA[0] = 2; IndexOC[0] = 0;
  A.contract(IndexA, OC, IndexOC);
  unsigned int position = 0;
  MPO0.set(position, A);
// A for position 1 <= l <= N-2:
  Shape[0] = D; Shape[1] = D; Shape[2] = 9;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 9; Index[1] = 9; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 1; Index[1] = 9; Index[2] = 1;
  A.set(Index, element);
  Index[0] = 3; Index[1] = 9; Index[2] = 2;
  A.set(Index, element);
  Index[0] = 2; Index[1] = 9; Index[2] = 3;
  A.set(Index, element);
  Index[0] = 4; Index[1] = 9; Index[2] = 4;
  A.set(Index, element);
  Index[0] = 6; Index[1] = 9; Index[2] = 5;
  A.set(Index, element);
  Index[0] = 5; Index[1] = 9; Index[2] = 6;
  A.set(Index, element);
  Index[0] = 8; Index[1] = 9; Index[2] = 7;
  A.set(Index, element);
  Index[0] = 7; Index[1] = 9; Index[2] = 8;
  A.set(Index, element);
  Tensor<T> AC(A);
  for (position = 1; position < this->N-1; position++)
  {
   A = AC;
   element = V0 * pow(T(position) - l0, 2) - mu;
   Index[0] = 0; Index[1] = 9; Index[2] = 1;
   A.set(Index, element);
// at even sites use t0, J0 and include U:
   if ((position % 2) == 0)
   {
    element = U - T(0.25) * J0;
    Index[0] = 0; Index[1] = 1; Index[2] = 1;
    A.set(Index, element);
    element = T(0.5) * J0;
    Index[0] = 0; Index[1] = 2; Index[2] = 2;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 3; Index[2] = 3;
    A.set(Index, element);
    element = J0;
    Index[0] = 0; Index[1] = 4; Index[2] = 4;
    A.set(Index, element);
    element = -t0;
    Index[0] = 0; Index[1] = 5; Index[2] = 5;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 6; Index[2] = 6;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 7; Index[2] = 7;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 8; Index[2] = 8;
    A.set(Index, element);
   }
// at odd sites use t1 and J1:
   else if ((position % 2) == 1)
   {
    element = -T(0.25) * J1;
    Index[0] = 0; Index[1] = 1; Index[2] = 1;
    A.set(Index, element);
    element = T(0.5) * J1;
    Index[0] = 0; Index[1] = 2; Index[2] = 2;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 3; Index[2] = 3;
    A.set(Index, element);
    element = J1;
    Index[0] = 0; Index[1] = 4; Index[2] = 4;
    A.set(Index, element);
    element = -t1;
    Index[0] = 0; Index[1] = 5; Index[2] = 5;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 6; Index[2] = 6;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 7; Index[2] = 7;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 8; Index[2] = 8;
    A.set(Index, element);
   }
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   MPO0.set(position, A);
  }
// A for position N-1:
  Shape[0] = D; Shape[1] = 1; Shape[2] = 9;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 9; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  element = V0 * pow(T(this->N-1) - l0, 2) - mu;
  Index[0] = 0; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  element = 1.0;
  Index[0] = 1; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  Index[0] = 3; Index[1] = 0; Index[2] = 2;
  A.set(Index, element);
  Index[0] = 2; Index[1] = 0; Index[2] = 3;
  A.set(Index, element);
  Index[0] = 4; Index[1] = 0; Index[2] = 4;
  A.set(Index, element);
  Index[0] = 6; Index[1] = 0; Index[2] = 5;
  A.set(Index, element);
  Index[0] = 5; Index[1] = 0; Index[2] = 6;
  A.set(Index, element);
  Index[0] = 8; Index[1] = 0; Index[2] = 7;
  A.set(Index, element);
  Index[0] = 7; Index[1] = 0; Index[2] = 8;
  A.set(Index, element);
  OC = O;
  A.contract(IndexA, OC, IndexOC);
  position = this->N-1;
  MPO0.set(position, A);
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void tJModel<T>::" <<
          "getMPO(MPO<T>& MPO0) const." << endl;
  exit(1);
 }
}

template<class T> void tJModel<T>::getMatrix(Matrix<T>& Matrix0) const
{
#ifdef DEBUG
 unsigned int dim = 1;
 for (int i = 0; i < this->N; i++)
 {
  dim *= 3;
 }
 if ((Matrix0.getDim0() != dim) || (Matrix0.getDim1() != dim))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void tJModel<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const: " <<
          "((Matrix0.getDim0() != 3^{N}) || (Matrix0.getDim1() != 3^{N}))." << endl;
  exit(1);
 }
#endif
 Matrix0.fillZeroes();
// the interactions:
// - the identity:
 Matrix<T> Identity(3, 3);
 Identity.fillZeroes();
 for (int i = 0; i < 3; i++)
 {
  Identity(i, i) = 1.0;
 }
// - particle number N =: N0:
 Matrix<T> N0(3, 3);
 N0.fillZeroes();
 N0(1, 1) = 1.0;
 N0(2, 2) = 1.0;
// - S^{+} =: S0:
 Matrix<T> S0(3, 3);
 S0.fillZeroes();
 S0(1, 2) = 1.0;
// - S^{-} =: S1:
 Matrix<T> S1(3, 3);
 S1.fillZeroes();
 S1(2, 1) = 1.0;
// - S^{z} =: S2:
 Matrix<T> S2(3, 3);
 S2.fillZeroes();
 S2(1, 1) = 0.5;
 S2(2, 2) = -0.5;
// - C_{\spinup}^{+} =: C0:
 Matrix<T> C0(3, 3);
 C0.fillZeroes();
 C0(1, 0) = 1.0;
// - C_{\spinup}^{-} =: C1:
 Matrix<T> C1(3, 3);
 C1.fillZeroes();
 C1(0, 1) = 1.0;
// - C_{\spindown}^{+} =: C2:
 Matrix<T> C2(3, 3);
 C2.fillZeroes();
 C2(2, 0) = 1.0;
// - C_{\spindown}^{-} =: C3:
 Matrix<T> C3(3, 3);
 C3.fillZeroes();
 C3(0, 2) = 1.0;
// the parameters:
 T t0 = this->Parameters[0];
 T J0 = this->Parameters[1];
 T t1 = this->Parameters[2];
 T J1 = this->Parameters[3];
 T U = this->Parameters[4];
 T mu = this->Parameters[5];
 T V0 = this->Parameters[6];
 T l0 = this->Parameters[7];
 T element;
 Matrix<T> X, Y;
// the Hamiltonian sum:
// - V0 \sum_{l} (l-l0)^{2} N_{l}:
 for (int i = 0; i < this->N; i++)
 {
  X = N0;
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
// - -\mu \sum_{l} N_{l}:
 for (int i = 0; i < this->N; i++)
 {
  X = N0;
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
// - U \sum_{l} N_{2l} \otimes N_{2l+1}:
// N1 := N_{2l} \otimes N_{2l+1}:
 Matrix<T> N1;
 N0.multiplyDirectProduct(N0, N1);
 for (int i = 0; i < this->N-1; i += 2)
 {
  X = N1;
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
  element = U;
  X.multiply(element);
  Matrix0.add(X);
 }
// - -t0 \sum_{l} C_{2l}:
// C := C_{\spinup}^{+} \otimes C_{\spinup}^{-} +
//      C_{\spinup}^{-} \otimes C_{\spinup}^{+} +
//      C_{\spindown}^{+} \otimes C_{\spindown}^{-} +
//      C_{\spindown}^{-} \otimes C_{\spindown}^{+}:
 Matrix<T> C(9, 9);
 C.fillZeroes();
 C0.multiplyDirectProduct(C1, X);
 C.add(X);
 C1.multiplyDirectProduct(C0, X);
 C.add(X);
 C2.multiplyDirectProduct(C3, X);
 C.add(X);
 C3.multiplyDirectProduct(C2, X);
 C.add(X);
 for (int i = 0; i < this->N-1; i += 2)
 {
  X = C;
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
// - J0 \sum_{l} S_{2l}:
// S := 0.5 * (S^{+} \otimes S^{-} + S^{-} \otimes S^{+}) +
//      S^{z} \otimes S^{z} -
//      0.25 * N \otimes N:
 Matrix<T> S(9, 9);
 S.fillZeroes();
 S0.multiplyDirectProduct(S1, X);
 element = 0.5;
 X.multiply(element);
 S.add(X);
 S1.multiplyDirectProduct(S0, X);
 X.multiply(element);
 S.add(X);
 S2.multiplyDirectProduct(S2, X);
 S.add(X);
 N0.multiplyDirectProduct(N0, X);
 element = -0.25;
 X.multiply(element);
 S.add(X);
 for (int i = 0; i < this->N-1; i += 2)
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
  element = J0;
  X.multiply(element);
  Matrix0.add(X);
 }
// - -t1 \sum_{l} C_{2l+1}
 for (int i = 1; i < this->N-1; i += 2)
 {
  X = C;
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
// - J1 \sum_{l} S_{2l+1}
 for (int i = 1; i < this->N-1; i += 2)
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
  element = J1;
  X.multiply(element);
  Matrix0.add(X);
 }
 if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void tJModel<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const." << endl;
  exit(1);
 }
 Matrix0.setType("hermitian");
}

template<class T> void tJModel<T>::getnMPO(unsigned int x, MPO<T>& nMPO) const
{
#ifdef DEBUG
 if ((this->N == 0) || (x > this->N-1))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void tJModel<T>::" <<
          "getnMPO(unsigned int x, MPO<T>& nMPO) const: " <<
          "((this->N == 0) || (x > this->N-1))." << endl;
  exit(1);
 }
#endif
// the particle number operator O = n_{l} = n_{\spindown, l} + n_{\spinup, l}:
 Matrix<T> O(3, 3);
 O.fillZeroes();
 O(1, 1) = 1.0;
 O(2, 2) = 1.0;
 getMPOFromLocalOperator(this->BC, this->N, x, O, nMPO);
}

template<class T> void tJModel<T>::getNMPO(MPO<T>& NMPO) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void tJModel<T>::" <<
          "getNMPO(MPO<T>& NMPO) const: " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
// the particle number operator O = n_{l} = n_{\spindown, l} + n_{\spinup, l}:
 Matrix<T> O(3, 3);
 O.fillZeroes();
 O(1, 1) = 1.0;
 O(2, 2) = 1.0;
 getMPOFromSumLocalOperator(this->BC, this->N, O, NMPO);
}

template<class T> void tJModel<T>::getsMPO(unsigned int x, MPO<T>& sMPO) const
{
#ifdef DEBUG
 if ((this->N == 0) || (x > this->N-1))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void tJModel<T>::" <<
          "getsMPO(unsigned int x, MPO<T>& sMPO) const: " <<
          "((this->N == 0) || (x > this->N-1))." << endl;
  exit(1);
 }
#endif
// the spin operator O = s_{l} = n_{\spinup, l} - n_{\spindown, l}:
 Matrix<T> O(3, 3);
 O.fillZeroes();
 O(1, 1) = 1.0;
 O(2, 2) = -1.0;
 getMPOFromLocalOperator(this->BC, this->N, x, O, sMPO);
}

template<class T> void tJModel<T>::getSMPO(MPO<T>& SMPO) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void tJModel<T>::" <<
          "getSMPO(MPO<T>& SMPO) const: " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
// the spin operator O = s_{l} = n_{\spinup, l} - n_{\spindown, l}:
 Matrix<T> O(3, 3);
 O.fillZeroes();
 O(1, 1) = 1.0;
 O(2, 2) = -1.0;
 getMPOFromSumLocalOperator(this->BC, this->N, O, SMPO);
}

template<class T> void tJModel<T>::getM2MPO(unsigned int leftEnd, unsigned int rightEnd,
                                            MPO<T>& M2MPO) const
{
#ifdef DEBUG
 if ((this->N == 0) || (leftEnd > this->N-2) || (rightEnd > this->N-1) ||
     (leftEnd >= rightEnd))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void tJModel<T>::" <<
          "getM2MPO(unsigned int leftEnd, unsigned int rightEnd, " <<
                   "MPO<T>& M2MPO) const: " <<
          "((this->N == 0) || (leftEnd > this->N-2) || (rightEnd > this->N-1) || " <<
           "(leftEnd >= rightEnd))." << endl;
  exit(1);
 }
#endif
 if (this->BC == "open")
 {
// the MPO operators O^{i}:
  vector<unsigned int> Shape(3);
  Shape[0] = 5; Shape[1] = 3; Shape[2] = 3;
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
// O^{1} = S^{+}:
  element = 1.0;
  Index[0] = 1; Index[1] = 2; Index[2] = 1;
  O.set(Index, element);
// O^{2} = S^{-}:
  element = 1.0;
  Index[0] = 2; Index[1] = 1; Index[2] = 2;
  O.set(Index, element);
// O^{3} = S^{z}:
  element = 0.5;
  Index[0] = 3; Index[1] = 1; Index[2] = 1;
  O.set(Index, element);
  element = -0.5;
  Index[0] = 3; Index[1] = 2; Index[2] = 2;
  O.set(Index, element);
// O^{4} = 3/4SpinIdentity:
  element = 0.75;
  Index[0] = 4; Index[1] = 1; Index[2] = 1;
  O.set(Index, element);
  Index[0] = 4; Index[1] = 2; Index[2] = 2;
  O.set(Index, element);
  unsigned int D = 5;
  M2MPO = MPO<T>(this->BC, this->N, this->d, D);
// A for position 0:
  Shape[0] = 1; Shape[1] = D; Shape[2] = 5;
  Tensor<T> A(Shape);
  A.fillZeroes();
  Tensor<T> OC(O);
  vector<unsigned int> IndexA(1), IndexOC(1);
  IndexA[0] = 2; IndexOC[0] = 0;
  int position = 0;
  if (position < leftEnd)
  {
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
  }
  else if (position == leftEnd)
  {
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 1; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 2; Index[2] = 2;
   A.set(Index, element);
   element = 2.0;
   Index[0] = 0; Index[1] = 3; Index[2] = 3;
   A.set(Index, element);
   element = 1.0;
   Index[0] = 0; Index[1] = 4; Index[2] = 4;
   A.set(Index, element);
  }
  A.contract(IndexA, OC, IndexOC);
  M2MPO.set(position, A);
// A for position 0 < l < leftEnd:
  if (1 < leftEnd)
  {
   Shape[0] = D; Shape[1] = D; Shape[2] = 5;
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   for (int position = 1; position < leftEnd; position++)
   {
    M2MPO.set(position, A);
   }
  }
// A for position leftEnd > 0:
  if (0 < leftEnd)
  {
   position = leftEnd;
   Shape[0] = D; Shape[1] = D; Shape[2] = 5;
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 1; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 2; Index[2] = 2;
   A.set(Index, element);
   element = 2.0;
   Index[0] = 0; Index[1] = 3; Index[2] = 3;
   A.set(Index, element);
   element = 1.0;
   Index[0] = 0; Index[1] = 4; Index[2] = 4;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   M2MPO.set(position, A);
  }
// A for position leftEnd < l < rightEnd:
  Shape[0] = D; Shape[1] = D; Shape[2] = 5;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 1; Index[1] = 1; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 2; Index[1] = 2; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 3; Index[1] = 3; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 4; Index[1] = 4; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 0; Index[1] = 4; Index[2] = 4;
  A.set(Index, element);
  Tensor<T> AC(A);
  for (position = leftEnd+1; position < rightEnd; position++)
  {
   A = AC;
   element = pow(-1.0, int(position-leftEnd));
   Index[0] = 0; Index[1] = 1; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 2; Index[2] = 2;
   A.set(Index, element);
   Index[0] = 2; Index[1] = 4; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 1; Index[1] = 4; Index[2] = 2;
   A.set(Index, element);
   element = 2.0 * pow(-1.0, int(position-leftEnd));
   Index[0] = 0; Index[1] = 3; Index[2] = 3;
   A.set(Index, element);
   element = pow(-1.0, int(position-leftEnd));
   Index[0] = 3; Index[1] = 4; Index[2] = 3;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   M2MPO.set(position, A);
  }
// A for position rightEnd < this->N-1:
  if (rightEnd < this->N-1)
  {
   position = rightEnd;
   Shape[0] = D; Shape[1] = D; Shape[2] = 5;
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0;
   Index[0] = 4; Index[1] = 4; Index[2] = 0;
   A.set(Index, element);
   element = pow(-1.0, int(rightEnd-leftEnd));
   Index[0] = 2; Index[1] = 4; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 1; Index[1] = 4; Index[2] = 2;
   A.set(Index, element);
   Index[0] = 3; Index[1] = 4; Index[2] = 3;
   A.set(Index, element);
   element = 1.0;
   Index[0] = 0; Index[1] = 4; Index[2] = 4;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   M2MPO.set(position, A);
  }
// A for position rightEnd < l < this->N-1:
  if (rightEnd < this->N-2)
  {
   Shape[0] = D; Shape[1] = D; Shape[2] = 5;
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0;
   Index[0] = 4; Index[1] = 4; Index[2] = 0;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   for (position = rightEnd+1; position < this->N-1; position++)
   {
    M2MPO.set(position, A);
   }
  }
// A for position this->N-1:
  Shape[0] = D; Shape[1] = 1; Shape[2] = 5;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  position = this->N-1;
  if (rightEnd < position)
  {
   element = 1.0;
   Index[0] = 4; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
  }
  else if (rightEnd == position)
  {
   element = 1.0;
   Index[0] = 4; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
   element = pow(-1.0, int(rightEnd-leftEnd));
   Index[0] = 2; Index[1] = 0; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 1; Index[1] = 0; Index[2] = 2;
   A.set(Index, element);
   Index[0] = 3; Index[1] = 0; Index[2] = 3;
   A.set(Index, element);
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 4;
   A.set(Index, element);
  }
  OC = O;
  A.contract(IndexA, OC, IndexOC);
  M2MPO.set(position, A);
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void tJModel<T>::" <<
          "getM2MPO(unsigned int leftEnd, unsigned int rightEnd, MPO<T>& M2MPO) const." << endl;
  exit(1);
 }
}

template<class T> void tJModel<T>::getDWHMMPO(unsigned int leftEnd, unsigned int rightEnd,
                                              MPO<T>& DWHMMPO) const
{
#ifdef DEBUG
 if ((this->N == 0) || (leftEnd > this->N-2) || (rightEnd > this->N-1) ||
     (leftEnd >= rightEnd))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void tJModel<T>::" <<
          "getDWHMMPO(unsigned int leftEnd, unsigned int rightEnd, " <<
                     "MPO<T>& DWHMMPO) const: " <<
          "((this->N == 0) || (leftEnd > this->N-2) || (rightEnd > this->N-1) || " <<
           "(leftEnd >= rightEnd))." << endl;
  exit(1);
 }
#endif
 if (this->BC == "open")
 {
// the MPO operators O^{i}:
  vector<unsigned int> Shape(3);
  Shape[0] = 4; Shape[1] = 3; Shape[2] = 3;
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
// O^{1} = S^{+}:
  element = 1.0;
  Index[0] = 1; Index[1] = 2; Index[2] = 1;
  O.set(Index, element);
// O^{2} = S^{-}:
  element = 1.0;
  Index[0] = 2; Index[1] = 1; Index[2] = 2;
  O.set(Index, element);
// O^{3} = S^{z}:
  element = 0.5;
  Index[0] = 3; Index[1] = 1; Index[2] = 1;
  O.set(Index, element);
  element = -0.5;
  Index[0] = 3; Index[1] = 2; Index[2] = 2;
  O.set(Index, element);
// the parameters:
  T J0 = this->Parameters[1];
  T J1 = this->Parameters[3];
  T J;
  unsigned int D = 5;
  DWHMMPO = MPO<T>(this->BC, this->N, this->d, D);
// A for position 0:
  Shape[0] = 1; Shape[1] = D; Shape[2] = 4;
  Tensor<T> A(Shape);
  A.fillZeroes();
  Tensor<T> OC(O);
  vector<unsigned int> IndexA(1), IndexOC(1);
  IndexA[0] = 2; IndexOC[0] = 0;
  int position = 0;
  if (position < leftEnd)
  {
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
  }
  else if (position == leftEnd)
  {
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
   element = 0.5 * J0;
   Index[0] = 0; Index[1] = 1; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 2; Index[2] = 2;
   A.set(Index, element);
   element = J0;
   Index[0] = 0; Index[1] = 3; Index[2] = 3;
   A.set(Index, element);
  }
  A.contract(IndexA, OC, IndexOC);
  DWHMMPO.set(position, A);
// A for position 0 < l < leftEnd:
  if (1 < leftEnd)
  {
   Shape[0] = D; Shape[1] = D; Shape[2] = 4;
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   for (int position = 1; position < leftEnd; position++)
   {
    DWHMMPO.set(position, A);
   }
  }
// A for position leftEnd > 0:
  if (0 < leftEnd)
  {
   position = leftEnd;
   Shape[0] = D; Shape[1] = D; Shape[2] = 4;
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
   if (leftEnd % 2 == 0)
    J = J0;
   else if (leftEnd % 2 == 1)
    J = J1;
   element = 0.5 * J;
   Index[0] = 0; Index[1] = 1; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 2; Index[2] = 2;
   A.set(Index, element);
   element = J;
   Index[0] = 0; Index[1] = 3; Index[2] = 3;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   DWHMMPO.set(position, A);
  }
// A for position leftEnd < l < rightEnd:
  Shape[0] = D; Shape[1] = D; Shape[2] = 4;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 4; Index[1] = 4; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 2; Index[1] = 4; Index[2] = 1;
  A.set(Index, element);
  Index[0] = 1; Index[1] = 4; Index[2] = 2;
  A.set(Index, element);
  Index[0] = 3; Index[1] = 4; Index[2] = 3;
  A.set(Index, element);
  Tensor<T> AC(A);
  for (position = leftEnd+1; position < rightEnd; position++)
  {
   A = AC;
   if (position % 2 == 0)
    J = J0;
   else if (position % 2 == 1)
    J = J1;
   element = 0.5 * J;
   Index[0] = 0; Index[1] = 1; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 2; Index[2] = 2;
   A.set(Index, element);
   element = J;
   Index[0] = 0; Index[1] = 3; Index[2] = 3;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   DWHMMPO.set(position, A);
  }
// A for position rightEnd < this->N-1:
  if (rightEnd < this->N-1)
  {
   position = rightEnd;
   Shape[0] = D; Shape[1] = D; Shape[2] = 4;
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0;
   Index[0] = 4; Index[1] = 4; Index[2] = 0;
   A.set(Index, element);
   Index[0] = 2; Index[1] = 4; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 1; Index[1] = 4; Index[2] = 2;
   A.set(Index, element);
   Index[0] = 3; Index[1] = 4; Index[2] = 3;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   DWHMMPO.set(position, A);
  }
// A for position rightEnd < l < this->N-1:
  if (rightEnd < this->N-2)
  {
   Shape[0] = D; Shape[1] = D; Shape[2] = 4;
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0;
   Index[0] = 4; Index[1] = 4; Index[2] = 0;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   for (position = rightEnd+1; position < this->N-1; position++)
   {
    DWHMMPO.set(position, A);
   }
  }
// A for position this->N-1:
  Shape[0] = D; Shape[1] = 1; Shape[2] = 4;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  position = this->N-1;
  if (rightEnd < position)
  {
   element = 1.0;
   Index[0] = 4; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
  }
  else if (rightEnd == position)
  {
   element = 1.0;
   Index[0] = 4; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
   Index[0] = 2; Index[1] = 0; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 1; Index[1] = 0; Index[2] = 2;
   A.set(Index, element);
   Index[0] = 3; Index[1] = 0; Index[2] = 3;
   A.set(Index, element);
  }
  OC = O;
  A.contract(IndexA, OC, IndexOC);
  DWHMMPO.set(position, A);
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void tJModel<T>::" <<
          "getDWHMMPO(unsigned int leftEnd, unsigned int rightEnd, MPO<T>& DWHMMPO) const." << endl;
  exit(1);
 }
}

template<class T> void tJModel<T>::getholeMPO(unsigned int x, MPO<T>& holeMPO) const
{
#ifdef DEBUG
 if ((this->N == 0) || (x > this->N-1))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void tJModel<T>::" <<
          "getholeMPO(unsigned int x, MPO<T>& holeMPO) const: " <<
          "((this->N == 0) || (x > this->N-1))." << endl;
  exit(1);
 }
#endif
// the hole number operator O = n_{l} = n_{\hole, l}:
 Matrix<T> O(3, 3);
 O.fillZeroes();
 O(0, 0) = 1.0;
 getMPOFromLocalOperator(this->BC, this->N, x, O, holeMPO);
}

template<class T> void tJModel<T>::getHoleMPO(MPO<T>& HoleMPO) const
{
#ifdef DEBUG
 if (this->N == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void tJModel<T>::" <<
          "getHoleMPO(MPO<T>& HoleMPO) const: " <<
          "(this->N == 0)." << endl;
  exit(1);
 }
#endif
// the hole number operator O = n_{l} = n_{\hole, l}:
 Matrix<T> O(3, 3);
 O.fillZeroes();
 O(0, 0) = 1.0;
 getMPOFromSumLocalOperator(this->BC, this->N, O, HoleMPO);
}

template<class T> void tJModel<T>::setHole(unsigned int position, MPS<T>& MPS0)
{
#ifdef DEBUG
 if (position > MPS0.getN()-1)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void tJModel<T>::" <<
          "setHole(unsigned int position, MPS<T>& MPS0): " <<
          "(position > MPS0.getN()-1)." << endl;
  exit(1);
 }
#endif
 Tensor<T> Tensor0;
 vector<unsigned int> Shape0, Index0(3);
 T element, element0;
 for (int l = 0; l < position; l++)
 {
  MPS0.get(l, Tensor0);
  Tensor0.getShape(Shape0);
  for (int k = 1; k < Shape0[2]; k++)
  {
   Index0[2] = k;
   for (int j = 0; j < Shape0[1]; j++)
   {
    Index0[1] = j;
    for (int i = 0; i < Shape0[0]; i++)
    {
     Index0[0] = i;
     element = Tensor0.get(Index0);
     element = -element;
     Tensor0.set(Index0, element);
    }
   }
  }
  MPS0.set(l, Tensor0);
 }
 MPS0.get(position, Tensor0);
 Tensor0.getShape(Shape0);
 for (int j = 0; j < Shape0[1]; j++)
 {
  Index0[1] = j;
  for (int i = 0; i < Shape0[0]; i++)
  {
   Index0[0] = i;
   Index0[2] = 1;
   element = Tensor0.get(Index0);
   Tensor0.set(Index0, 0.0);
   Index0[2] = 2;
   element0 = Tensor0.get(Index0);
   Tensor0.set(Index0, 0.0);
   element = element+element0;
   Index0[2] = 0;
   Tensor0.set(Index0, element);
  }
 }
 MPS0.set(position, Tensor0);
}

template<class T> void tJModel<T>::getSingleParticle(unsigned int position, const vector<T>& Coefficients, MPS<T>& MPS0)
{
 string BC0; MPS0.getBC(BC0);
 unsigned int N0 = MPS0.getN();
#ifdef DEBUG
 vector<unsigned int> d0(N0); MPS0.getd(d0);
 if ((N0 == 0) || (position > N0-1) || (Coefficients.size() != 3))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void tJModel<T>::" <<
          "getSingleParticle(unsigned int position, const vector<T>& Coefficients, MPS<T>& MPS0): " <<
          "((MPS0.getN() == 0) || (position > MPS0.getN()-1) || (Coefficients.size() != 3))." << endl;
  exit(1);
 }
 for (int i = 0; i < N0; i++)
 {
  if (d0[i] != 3)
  {
   cerr << "Program terminated because of error in static function " <<
           "template<class T> void tJModel<T>::" <<
           "getSingleParticle(unsigned int position, const vector<T>& Coefficients, MPS<T>& MPS0): " <<
           "(MPS0.d != 3)." << endl;
   exit(1);
  }
 }
#endif
 vector<unsigned int> Shape(3), Index(3);
 T element;
 Tensor<T> A;
 if (BC0 == "open")
 {
  for (int x = 0; x < N0; x++)
  {
   MPS0.getOpenBCShape(x, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   if (x != position)
   {
    element = 1.0;
    Index[0] = 0; Index[1] = 0; Index[2] = 0;
    A.set(Index, element);
    MPS0.set(x, A);
   }
   else if (x == position)
   {
    element = Coefficients[0];
    Index[0] = 0; Index[1] = 0; Index[2] = 0;
    A.set(Index, element);
    element = Coefficients[1];
    Index[0] = 0; Index[1] = 0; Index[2] = 1;
    A.set(Index, element);
    element = Coefficients[2];
    Index[0] = 0; Index[1] = 0; Index[2] = 2;
    A.set(Index, element);
    MPS0.set(x, A);
   }
  }
 }
 else if (BC0 == "periodic")
 {
  cerr << "The following static function is not implemented yet for periodic BC: " <<
          "template<class T> void tJModel<T>::" <<
          "getSingleParticle(unsigned int position, const vector<T>& Coefficients, MPS<T>& MPS0)." << endl;
  exit(1);
 }
}

template<class T> void tJModel<T>::getDWSingleParticle(const string& WhichEigenstate, T Cup, T Cdown, MPS<T>& MPS0)
{
#ifdef DEBUG
 if ((WhichEigenstate != "GS") && (WhichEigenstate != "ES"))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void tJModel<T>::" <<
          "getDWSingleParticle(const string& WhichEigenstate, T Cup, T Cdown, MPS<T>& MPS0): " <<
          "((WhichEigenstate != GS) && (WhichEigenstate != ES))." << endl;
  exit(1);
 }
#endif
 string BC0 = "open";
 unsigned int N0 = 2;
 unsigned int d0 = 3;
 unsigned int D0 = 3;
 MPS0 = MPS<T>(BC0, N0, d0, D0);
 unsigned int position;
 vector<unsigned int> Shape(3), Index(3);
 Tensor<T> A;
 T element;
 if (WhichEigenstate == "GS")
 {
  position = 0;
  MPS0.getOpenBCShape(position, Shape);
  A = Tensor<T>(Shape);
  A.fillZeroes();
  Index[0] = 0; Index[1] = 1; Index[2] = 0;
  element = Cup / sqrt(2.0);
  A.set(Index, element);
  Index[0] = 0; Index[1] = 2; Index[2] = 0;
  element = Cdown / sqrt(2.0);
  A.set(Index, element);
  Index[0] = 0; Index[1] = 0; Index[2] = 1;
  element = Cup / sqrt(2.0);
  A.set(Index, element);
  Index[0] = 0; Index[1] = 0; Index[2] = 2;
  element = Cdown / sqrt(2.0);
  A.set(Index, element);
  MPS0.set(position, A);
  position = 1;
  MPS0.getOpenBCShape(position, Shape);
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 1; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  Index[0] = 2; Index[1] = 0; Index[2] = 2;
  A.set(Index, element);
  MPS0.set(position, A);
 }
 else if (WhichEigenstate == "ES")
 {
  position = 0;
  MPS0.getOpenBCShape(position, Shape);
  A = Tensor<T>(Shape);
  A.fillZeroes();
  Index[0] = 0; Index[1] = 1; Index[2] = 0;
  element = -Cup / sqrt(2.0);
  A.set(Index, element);
  Index[0] = 0; Index[1] = 2; Index[2] = 0;
  element = -Cdown / sqrt(2.0);
  A.set(Index, element);
  Index[0] = 0; Index[1] = 0; Index[2] = 1;
  element = Cup / sqrt(2.0);
  A.set(Index, element);
  Index[0] = 0; Index[1] = 0; Index[2] = 2;
  element = Cdown / sqrt(2.0);
  A.set(Index, element);
  MPS0.set(position, A);
  position = 1;
  MPS0.getOpenBCShape(position, Shape);
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 1; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  Index[0] = 2; Index[1] = 0; Index[2] = 2;
  A.set(Index, element);
  MPS0.set(position, A);
 }
}

template<class T> void tJModel<T>::getDimer(const vector<unsigned int>& HolePositions, MPS<T>& MPS0)
{
 string BC0; MPS0.getBC(BC0);
 unsigned int N0 = MPS0.getN();
#ifdef DEBUG
 vector<unsigned int> d0(N0); MPS0.getd(d0);
 unsigned int D0 = MPS0.getD();
 if ((N0 == 0) || (D0 < 2))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void tJModel<T>::" <<
          "getDimer(const vector<unsigned int>& HolePositions, MPS<T>& MPS0): " <<
          "((MPS0.getN() == 0) || (MPS0.getD() < 2))." << endl;
  exit(1);
 }
 for (int i = 0; i < N0; i++)
 {
  if (d0[i] != 3)
  {
   cerr << "Program terminated because of error in static function " <<
           "template<class T> void tJModel<T>::" <<
           "getDimer(const vector<unsigned int>& HolePositions, MPS<T>& MPS0): " <<
           "(MPS0.d != 3)." << endl;
   exit(1);
  }
 }
 for (int i = 0; i < HolePositions.size(); i++)
 {
  if (HolePositions[i] > N0-1)
  {
   cerr << "Program terminated because of error in static function " <<
           "template<class T> void tJModel<T>::" <<
           "getDimer(const vector<unsigned int>& HolePositions, MPS<T>& MPS0): " <<
           "(HolePositions > MPS0.getN()-1)." << endl;
   exit(1);
  }
 }
#endif
// NumHoles will have 0 on sites with no hole and 1 on sites with hole:
 vector<unsigned int> NumHoles(N0);
 for (int i = 0; i < N0; i++)
 {
  NumHoles[i] = 0;
 }
 for (int i = 0; i < HolePositions.size(); i++)
 {
  NumHoles[HolePositions[i]] = 1;
 }
 vector<unsigned int> Shape(3), Index(3);
 T element;
 Tensor<T> A;
 unsigned int position;
 if (BC0 == "open")
 {
// A for position 0 and position 1:
  position = 0;
  MPS0.getOpenBCShape(position, Shape);
  A = Tensor<T>(Shape);
  A.fillZeroes();
  if ((NumHoles[position] == 0) && (NumHoles[position+1] == 0))
  {
   element = 1.0 / sqrt(2.0);
   Index[0] = 0; Index[1] = 0; Index[2] = 1;
   A.set(Index, element);
   element = -1.0 / sqrt(2.0);
   Index[0] = 0; Index[1] = 1; Index[2] = 2;
   A.set(Index, element);
   MPS0.set(position, A);
   position = 1;
   MPS0.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0;
   Index[0] = 1; Index[1] = 0; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 0; Index[2] = 2;
   A.set(Index, element);
   MPS0.set(position, A);
  }
  else if ((NumHoles[position] == 0) && (NumHoles[position+1] == 1))
  {
   element = 1.0 / sqrt(2.0);
   Index[0] = 0; Index[1] = 0; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 0; Index[2] = 2;
   A.set(Index, element);
   MPS0.set(position, A);
   position = 1;
   MPS0.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
   MPS0.set(position, A);
  }
  else if ((NumHoles[position] == 1) && (NumHoles[position+1] == 0))
  {
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
   MPS0.set(position, A);
   position = 1;
   MPS0.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0 / sqrt(2.0);
   Index[0] = 0; Index[1] = 0; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 0; Index[2] = 2;
   A.set(Index, element);
   MPS0.set(position, A);
  }
  else if ((NumHoles[position] == 1) && (NumHoles[position+1] == 1))
  {
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
   MPS0.set(position, A);
   position = 1;
   MPS0.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
   MPS0.set(position, A);
  }
// A for positions 1 < position < N0-2:
  for (position = 2; position < N0-2; position = position+2)
  {
   MPS0.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   if ((NumHoles[position] == 0) && (NumHoles[position+1] == 0))
   {
    element = 1.0 / sqrt(2.0);
    Index[0] = 0; Index[1] = 0; Index[2] = 1;
    A.set(Index, element);
    element = -1.0 / sqrt(2.0);
    Index[0] = 0; Index[1] = 1; Index[2] = 2;
    A.set(Index, element);
    MPS0.set(position, A);
    MPS0.getOpenBCShape(position+1, Shape);
    A = Tensor<T>(Shape);
    A.fillZeroes();
    element = 1.0;
    Index[0] = 1; Index[1] = 0; Index[2] = 1;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 0; Index[2] = 2;
    A.set(Index, element);
    MPS0.set(position+1, A);
   }
   else if ((NumHoles[position] == 0) && (NumHoles[position+1] == 1))
   {
    element = 1.0 / sqrt(2.0);
    Index[0] = 0; Index[1] = 0; Index[2] = 1;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 0; Index[2] = 2;
    A.set(Index, element);
    MPS0.set(position, A);
    MPS0.getOpenBCShape(position+1, Shape);
    A = Tensor<T>(Shape);
    A.fillZeroes();
    element = 1.0;
    Index[0] = 0; Index[1] = 0; Index[2] = 0;
    A.set(Index, element);
    MPS0.set(position+1, A);
   }
   else if ((NumHoles[position] == 1) && (NumHoles[position+1] == 0))
   {
    element = 1.0;
    Index[0] = 0; Index[1] = 0; Index[2] = 0;
    A.set(Index, element);
    MPS0.set(position, A);
    MPS0.getOpenBCShape(position+1, Shape);
    A = Tensor<T>(Shape);
    A.fillZeroes();
    element = 1.0 / sqrt(2.0);
    Index[0] = 0; Index[1] = 0; Index[2] = 1;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 0; Index[2] = 2;
    A.set(Index, element);
    MPS0.set(position+1, A);
   }
   else if ((NumHoles[position] == 1) && (NumHoles[position+1] == 1))
   {
    element = 1.0;
    Index[0] = 0; Index[1] = 0; Index[2] = 0;
    A.set(Index, element);
    MPS0.set(position, A);
    MPS0.getOpenBCShape(position+1, Shape);
    A = Tensor<T>(Shape);
    A.fillZeroes();
    element = 1.0;
    Index[0] = 0; Index[1] = 0; Index[2] = 0;
    A.set(Index, element);
    MPS0.set(position+1, A);
   }
  }
// If N0 is even, A for remaining positions N0-2 and N0-1:
  if (N0 % 2 == 0)
  {
   position = N0-2;
   MPS0.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   if ((NumHoles[position] == 0) && (NumHoles[position+1] == 0))
   {
    element = 1.0 / sqrt(2.0);
    Index[0] = 0; Index[1] = 0; Index[2] = 1;
    A.set(Index, element);
    element = -1.0 / sqrt(2.0);
    Index[0] = 0; Index[1] = 1; Index[2] = 2;
    A.set(Index, element);
    MPS0.set(position, A);
    position = N0-1;
    MPS0.getOpenBCShape(position, Shape);
    A = Tensor<T>(Shape);
    A.fillZeroes();
    element = 1.0;
    Index[0] = 1; Index[1] = 0; Index[2] = 1;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 0; Index[2] = 2;
    A.set(Index, element);
    MPS0.set(position, A);
   }
   else if ((NumHoles[position] == 0) && (NumHoles[position+1] == 1))
   {
    element = 1.0 / sqrt(2.0);
    Index[0] = 0; Index[1] = 0; Index[2] = 1;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 0; Index[2] = 2;
    A.set(Index, element);
    MPS0.set(position, A);
    position = N0-1;
    MPS0.getOpenBCShape(position, Shape);
    A = Tensor<T>(Shape);
    A.fillZeroes();
    element = 1.0;
    Index[0] = 0; Index[1] = 0; Index[2] = 0;
    A.set(Index, element);
    MPS0.set(position, A);
   }
   else if ((NumHoles[position] == 1) && (NumHoles[position+1] == 0))
   {
    element = 1.0;
    Index[0] = 0; Index[1] = 0; Index[2] = 0;
    A.set(Index, element);
    MPS0.set(position, A);
    position = N0-1;
    MPS0.getOpenBCShape(position, Shape);
    A = Tensor<T>(Shape);
    A.fillZeroes();
    element = 1.0 / sqrt(2.0);
    Index[0] = 0; Index[1] = 0; Index[2] = 1;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 0; Index[2] = 2;
    A.set(Index, element);
    MPS0.set(position, A);
   }
   else if ((NumHoles[position] == 1) && (NumHoles[position+1] == 1))
   {
    element = 1.0;
    Index[0] = 0; Index[1] = 0; Index[2] = 0;
    A.set(Index, element);
    MPS0.set(position, A);
    position = N0-1;
    MPS0.getOpenBCShape(position, Shape);
    A = Tensor<T>(Shape);
    A.fillZeroes();
    element = 1.0;
    Index[0] = 0; Index[1] = 0; Index[2] = 0;
    A.set(Index, element);
    MPS0.set(position, A);
   }
  }
// else if N0 is odd, A for remaining position N0-1:
  else if (N0 % 2 == 1)
  {
   position = N0-1;
   MPS0.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   if (NumHoles[position] == 0)
   {
    element = 1.0 / sqrt(2.0);
    Index[0] = 0; Index[1] = 0; Index[2] = 1;
    A.set(Index, element);
    Index[0] = 0; Index[1] = 0; Index[2] = 2;
    A.set(Index, element);
    MPS0.set(position, A);
   }
   else if (NumHoles[position] == 1)
   {
    element = 1.0;
    Index[0] = 0; Index[1] = 0; Index[2] = 0;
    A.set(Index, element);
    MPS0.set(position, A);
   }
  }
 }
 else if (BC0 == "periodic")
 {
  cerr << "The following static function is not implemented yet for periodic BC: " <<
          "template<class T> void tJModel<T>::" <<
          "getDimer(const vector<unsigned int>& HolePositions, MPS<T>& MPS0)." << endl;
  exit(1);
 }
}

template<class T> void tJModel<T>::getNeel(MPS<T>& MPS0)
{
 unsigned int N0 = MPS0.getN();
#ifdef DEBUG
 if ((N0 == 0) || (MPS0.getD() < 2) || (MPS0.getd() != 3))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void tJModel<T>::" <<
          "getNeel(MPS<T>& MPS0): " <<
          "((MPS0.getN() == 0) || (MPS0.getD() < 2) || (MPS0.getd() != 3))." << endl;
  exit(1);
 }
#endif
 string BC0; MPS0.getBC(BC0);
 vector<unsigned int> Shape(3), Index(3);
 T element;
 unsigned int position;
 if (BC0 == "open")
 {
// A for position 0:
  position = 0;
  MPS0.getOpenBCShape(position, Shape);
  Tensor<T> A(Shape);
  A.fillZeroes();
  element = 1.0 / sqrt(2.0);
  Index[0] = 0; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  Index[0] = 0; Index[1] = 1; Index[2] = 2;
  A.set(Index, element);
  MPS0.set(position, A);
// A for position 1 <= l <= N0-2:
  element = 1.0;
  for (position = 1; position < N0-1; position++)
  {
   MPS0.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   Index[0] = 1; Index[1] = 0; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 1; Index[2] = 2;
   A.set(Index, element);
   MPS0.set(position, A);
  }
// A for position N0-1:
  position = N0-1;
  MPS0.getOpenBCShape(position, Shape);
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 1; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  Index[0] = 0; Index[1] = 0; Index[2] = 2;
  A.set(Index, element);
  MPS0.set(position, A);
 }
 else if (BC0 == "periodic")
 {
  cerr << "The following static function is not implemented yet for periodic BC: " <<
          "template<class T> void tJModel<T>::" <<
          "getNeel(MPS<T>& MPS0)." << endl;
  exit(1);
 }
}

template<class T> void tJModel<T>::setHoles(const MPS<T>& MPS0, unsigned int numHolesLeft,
                                            unsigned int numHolesRight, MPS<T>& MPS1)
{
 unsigned int N0 = MPS0.getN();
#ifdef DEBUG
 if (N0 == 0)
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void tJModel<T>::" <<
          "setHoles(const MPS<T>& MPS0, unsigned int numHolesLeft, " <<
                   "unsigned int numHolesRight, MPS<T>& MPS1): " <<
          "(MPS0.getN() == 0)." << endl;
  exit(1);
 }
#endif
 unsigned int N1 = numHolesLeft + N0 + numHolesRight;
 unsigned int posLeft = numHolesLeft, posRight = numHolesLeft + N0 - 1;
 string BC0;
 MPS0.getBC(BC0);
 unsigned int D0 = MPS0.getD();
 vector<unsigned int> d0(N0), d1(N1);
 MPS0.getd(d0);
 for (int i = 0; i < posLeft; i++)
 {
  d1[i] = 3;
 }
 for (int i = posLeft; i <= posRight; i++)
 {
  d1[i] = d0[i-posLeft];
 }
 for (int i = posRight+1; i < N1; i++)
 {
  d1[i] = 3;
 }
 MPS1 = MPS<T>(BC0, N1, d1, D0);
 vector<unsigned int> Shape(3), Index(3), Shape0(3);
 Tensor<T> A, Tensor0;
 T element;
 unsigned int position;
 if (BC0 == "open")
 {
// A for position 0:
  position = 0;
  MPS1.getOpenBCShape(position, Shape);
  A = Tensor<T>(Shape);
// if position 0 < posLeft:
  if (position < posLeft)
  {
   A.fillZeroes();
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
  }
// else if position == posLeft:
  else if (position == posLeft)
  {
   MPS0.get(position, Tensor0);
   A = Tensor0;
  }
  MPS1.set(position, A);
// A for position 0 < l < posLeft:
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  for (position = 1; position < posLeft; position++)
  {
   MPS1.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   A.set(Index, element);
   MPS1.set(position, A);
  }
// A for position posLeft > 0:
  if (posLeft > 0)
  {
   position = posLeft;
   MPS1.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   MPS0.get(0, Tensor0);
   Tensor0.getShape(Shape0);
   Index[0] = 0;
   for (int j = 0; j < Shape0[2]; j++)
   {
    Index[2] = j;
    for (int i = 0; i < Shape0[1]; i++)
    {
     Index[1] = i;
     element = Tensor0.get(Index);
     A.set(Index, element);
    }
   }
   MPS1.set(position, A);
  }
// A for position posLeft < l < posRight:
  for (position = posLeft+1; position < posRight; position++)
  {
   MPS1.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   MPS0.get(position-posLeft, Tensor0);
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
// A for position posRight < N1-1:
  if (posRight < N1-1)
  {
   position = posRight;
   MPS1.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   MPS0.get(N0-1, Tensor0);
   Tensor0.getShape(Shape0);
   Index[1] = 0;
   for (int j = 0; j < Shape0[2]; j++)
   {
    Index[2] = j;
    for (int i = 0; i < Shape0[0]; i++)
    {
     Index[0] = i;
     element = Tensor0.get(Index);
     A.set(Index, element);
    }
   }
   MPS1.set(position, A);
  }
// A for position posRight < l < N1-1:
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  for (position = posRight+1; position < N1-1; position++)
  {
   MPS1.getOpenBCShape(position, Shape);
   A = Tensor<T>(Shape);
   A.fillZeroes();
   A.set(Index, element);
   MPS1.set(position, A);
  }
// A for position N1-1:
  position = N1-1;
  MPS1.getOpenBCShape(position, Shape);
  A = Tensor<T>(Shape);
// if position N1-1 > posRight:
  if (position > posRight)
  {
   A.fillZeroes();
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
  }
// else if position == posRight:
  else if (position == posRight)
  {
   MPS0.get(N0-1, Tensor0);
   A = Tensor0;
  }
  MPS1.set(position, A);
 }
 else if (BC0 == "periodic")
 {
  cerr << "The following static function is not implemented yet for periodic BC: " <<
          "template<class T> void tJModel<T>::" <<
          "setHoles(const MPS<T>& MPS0, unsigned int numHolesLeft, " <<
                   "unsigned int numHolesRight, MPS<T>& MPS1)." << endl;
  exit(1);
 }
}

template<class T> void tJModel<T>::getSpinSpinCorrelator(unsigned int position0, unsigned int position1, MPO<T>& MPO0)
{
 unsigned int N0 = MPO0.getN();
 string BC0; MPO0.getBC(BC0);
#ifdef DEBUG
 if ((N0 < 2) || (position0 > N0-2) || (position1 > N0-1) || (position0 >= position1))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void tJModel<T>::" <<
          "getSpinSpinCorrelator(unsigned int position0, unsigned int position1, MPO<T>& MPO0): " <<
          "((MPO0.getN() < 2) || (position0 > MPO0.getN()-2) || (position1 > MPO0.getN()-1) || (position0 >= position1))." << endl;
  exit(1);
 }
#endif
 if (BC0 == "open")
 {
// the MPO operators O^{i}:
  vector<unsigned int> Shape(3);
  Shape[0] = 4; Shape[1] = 3; Shape[2] = 3;
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
// O^{1} = S^{+}:
  element = 1.0;
  Index[0] = 1; Index[1] = 2; Index[2] = 1;
  O.set(Index, element);
// O^{2} = S^{-}:
  element = 1.0;
  Index[0] = 2; Index[1] = 1; Index[2] = 2;
  O.set(Index, element);
// O^{3} = S^{z}:
  element = 0.5;
  Index[0] = 3; Index[1] = 1; Index[2] = 1;
  O.set(Index, element);
  element = -0.5;
  Index[0] = 3; Index[1] = 2; Index[2] = 2;
  O.set(Index, element);
  unsigned int d0 = 3;
  unsigned int D = 3;
  MPO0 = MPO<T>(BC0, N0, d0, D);
// A for position 0:
  Shape[0] = 1; Shape[1] = D; Shape[2] = 4;
  Tensor<T> A(Shape);
  A.fillZeroes();
  Tensor<T> OC(O);
  vector<unsigned int> IndexA(1), IndexOC(1);
  IndexA[0] = 2; IndexOC[0] = 0;
  int position = 0;
  if (position < position0)
  {
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
  }
  else if (position == position0)
  {
   element = 0.5;
   Index[0] = 0; Index[1] = 0; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 1; Index[2] = 2;
   A.set(Index, element);
   element = 1.0;
   Index[0] = 0; Index[1] = 2; Index[2] = 3;
   A.set(Index, element);
  }
  A.contract(IndexA, OC, IndexOC);
  MPO0.set(position, A);
// A for position 0 < l < position0:
  if (1 < position0)
  {
   Shape[0] = D; Shape[1] = D; Shape[2] = 4;
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   for (position = 1; position < position0; position++)
   {
    MPO0.set(position, A);
   }
  }
// A for position position0 > 0:
  if (0 < position0)
  {
   position = position0;
   Shape[0] = D; Shape[1] = D; Shape[2] = 4;
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 0.5;
   Index[0] = 0; Index[1] = 0; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 1; Index[2] = 2;
   A.set(Index, element);
   element = 1.0;
   Index[0] = 0; Index[1] = 2; Index[2] = 3;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   MPO0.set(position, A);
  }
// A for position position0 < l < position1:
  if (position0+1 < position1)
  {
   Shape[0] = D; Shape[1] = D; Shape[2] = 4;
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
   Index[0] = 1; Index[1] = 1; Index[2] = 0;
   A.set(Index, element);
   Index[0] = 2; Index[1] = 2; Index[2] = 0;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   for (position = position0+1; position < position1; position++)
   {
    MPO0.set(position, A);
   }
  }
// A for position position1 < MPO0.getN()-1:
  if (position1 < N0-1)
  {
   position = position1;
   Shape[0] = D; Shape[1] = D; Shape[2] = 4;
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0;
   Index[0] = 1; Index[1] = 0; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 0; Index[2] = 2;
   A.set(Index, element);
   Index[0] = 2; Index[1] = 0; Index[2] = 3;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   MPO0.set(position, A);
  }
// A for position position1 < l < MPO0.getN()-1:
  if (position1 < N0-2)
  {
   Shape[0] = D; Shape[1] = D; Shape[2] = 4;
   A = Tensor<T>(Shape);
   A.fillZeroes();
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
   OC = O;
   A.contract(IndexA, OC, IndexOC);
   for (position = position1+1; position < N0-1; position++)
   {
    MPO0.set(position, A);
   }
  }
// A for position MPO0.getN()-1:
  Shape[0] = D; Shape[1] = 1; Shape[2] = 4;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  position = N0-1;
  if (position1 < position)
  {
   element = 1.0;
   Index[0] = 0; Index[1] = 0; Index[2] = 0;
   A.set(Index, element);
  }
  else if (position1 == position)
  {
   element = 1.0;
   Index[0] = 1; Index[1] = 0; Index[2] = 1;
   A.set(Index, element);
   Index[0] = 0; Index[1] = 0; Index[2] = 2;
   A.set(Index, element);
   Index[0] = 2; Index[1] = 0; Index[2] = 3;
   A.set(Index, element);
  }
  OC = O;
  A.contract(IndexA, OC, IndexOC);
  MPO0.set(position, A);
 }
 else if (BC0 == "periodic")
 {
  cerr << "The following static function is not implemented yet for periodic BC: " <<
          "template<class T> void tJModel<T>::" <<
          "getSpinSpinCorrelator(unsigned int position0, unsigned int position1, MPO<T>& MPO0)." << endl;
  exit(1);
 }
}

template<class T> void tJModel<T>::computeDWHMEigenstates(const string& BC, unsigned int N, T J0, T J1,
                                                          unsigned int nev, vector<T>& Evals,
                                                          Matrix<T>& Evecs)
{
#ifdef DEBUG
 if (((BC != "open") && (BC != "periodic")) || (N%2 != 0))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void tJModel<T>::" <<
          "computeDWHMEigenstates(const string& BC, unsigned int N, T J0, T J1, " <<
                                 "unsigned int nev, vector<T>& Evals, Matrix<T>& Evecs): " <<
          "(((BC != open) && (BC != periodic)) || (N%2 != 0))." << endl;
  exit(1);
 }
 else if ((typeid(T) != typeid(float)) && (typeid(T) != typeid(double)))
 {
  cerr << "Program terminated because of error in static function " <<
          "template<class T> void tJModel<T>::" <<
          "computeDWHMEigenstates(const string& BC, unsigned int N, T J0, T J1, " <<
                                 "unsigned int nev, vector<T>& Evals, Matrix<T>& Evecs): " <<
          "This function is not implemented yet for complex matrices." << endl;
  exit(1);
 }
#endif
// initialization:
 Matrix<int> Identity(2, 2), S0(2, 2), S1(2, 2), S2(2, 2);
 vector<T> AD, AOD; vector<int> IndexD, IndexODi, IndexODj;
 vector<int> Indexi(N), Indexj(N);
 int i0, j0, resultJ0, resultJ1, result;
 int dim = 1;
 for (int i = 0; i < N; i++)
  dim *= 2;
 int n0, nev0, ncv;
 n0 = dim; nev0 = nev;
 if (nev <= 10)
  ncv = min(60, dim);
 else
  ncv = min(nev0+50, dim);
 Evals = vector<T>(nev);
 Evecs = Matrix<T>(dim, nev);
// the single site matrices:
// - the identity:
 Identity.fillZeroes();
 Identity(0, 0) = 1;
 Identity(1, 1) = 1;
// - S^{+} =: S0:
 S0.fillZeroes();
 S0(0, 1) = 1;
// - S^{-} =: S1:
 S1.fillZeroes();
 S1(1, 0) = 1;
// - S^{z} =: S2:
 S2.fillZeroes();
 S2(0, 0) = 1;
 S2(1, 1) = -1;
// computing the diagonal part: Jx S^{z} \otimes S^{z}
 for (int i = 0; i < dim; i++)
 {
  i0 = i;
  for (int k = N-1; k >= 0; k--)
  {
   Indexi[k] = i0 % 2;
   i0 -= i0 % 2;
   i0 /= 2;
  }
// coupling J0 for even sites:
  resultJ0 = 0;
  for (int j = 0; j < N-1; j += 2)
  {
   result = 1;
   for (int k = 0; k < j; k++)
    result *= Identity(Indexi[k], Indexi[k]);
   result *= S2(Indexi[j], Indexi[j]);
   result *= S2(Indexi[j+1], Indexi[j+1]);
   for (int k = j+2; k < N; k++)
    result *= Identity(Indexi[k], Indexi[k]);
   resultJ0 += result;
  }
// coupling J1 for odd sites:
  resultJ1 = 0;
  for (int j = 1; j < N-1; j += 2)
  {
   result = 1;
   for (int k = 0; k < j; k++)
    result *= Identity(Indexi[k], Indexi[k]);
   result *= S2(Indexi[j], Indexi[j]);
   result *= S2(Indexi[j+1], Indexi[j+1]);
   for (int k = j+2; k < N; k++)
    result *= Identity(Indexi[k], Indexi[k]);
   resultJ1 += result;
  }
  if (BC == "periodic")
  {
   result = 1;
   result *= S2(Indexi[0], Indexi[0]);
   for (int k = 1; k < N-1; k++)
    result *= Identity(Indexi[k], Indexi[k]);
   result *= S2(Indexi[N-1], Indexi[N-1]);
   resultJ1 += result;
  }
  if ((resultJ0 != 0) || (resultJ1 != 0))
  {
   AD.push_back(0.25*J0*T(resultJ0)+0.25*J1*T(resultJ1));
   IndexD.push_back(i+1);
  }
 }
// computing the offdiagonal part: 0.5*Jx S^{+/-} \otimes S^{-/+}
 for (int i = 0; i < dim; i++)
 {
  i0 = i;
  for (int k = N-1; k >= 0; k--)
  {
   Indexi[k] = i0 % 2;
   i0 -= i0 % 2;
   i0 /= 2;
  }
  for (int j = i+1; j < dim; j++)
  {
   j0 = j;
   for (int k = N-1; k >= 0; k--)
   {
    Indexj[k] = j0 % 2;
    j0 -= j0 % 2;
    j0 /= 2;
   }
// coupling J0 for even sites:
   resultJ0 = 0;
   for (int k = 0; k < N-1; k += 2)
   {
    result = 1;
    for (int l = 0; l < k; l++)
     result *= Identity(Indexi[l], Indexj[l]);
    result *= S0(Indexi[k], Indexj[k]);
    result *= S1(Indexi[k+1], Indexj[k+1]);
    for (int l = k+2; l < N; l++)
     result *= Identity(Indexi[l], Indexj[l]);
    resultJ0 += result;
    result = 1;
    for (int l = 0; l < k; l++)
     result *= Identity(Indexi[l], Indexj[l]);
    result *= S1(Indexi[k], Indexj[k]);
    result *= S0(Indexi[k+1], Indexj[k+1]);
    for (int l = k+2; l < N; l++)
     result *= Identity(Indexi[l], Indexj[l]);
    resultJ0 += result;
   }
// coupling J1 for odd sites:
   resultJ1 = 0;
   for (int k = 1; k < N-1; k += 2)
   {
    result = 1;
    for (int l = 0; l < k; l++)
     result *= Identity(Indexi[l], Indexj[l]);
    result *= S0(Indexi[k], Indexj[k]);
    result *= S1(Indexi[k+1], Indexj[k+1]);
    for (int l = k+2; l < N; l++)
     result *= Identity(Indexi[l], Indexj[l]);
    resultJ1 += result;
    result = 1;
    for (int l = 0; l < k; l++)
     result *= Identity(Indexi[l], Indexj[l]);
    result *= S1(Indexi[k], Indexj[k]);
    result *= S0(Indexi[k+1], Indexj[k+1]);
    for (int l = k+2; l < N; l++)
     result *= Identity(Indexi[l], Indexj[l]);
    resultJ1 += result;
   }
   if (BC == "periodic")
   {
    result = 1;
    result *= S0(Indexi[0], Indexj[0]);
    for (int k = 1; k < N-1; k++)
     result *= Identity(Indexi[k], Indexj[k]);
    result *= S1(Indexi[N-1], Indexj[N-1]);
    resultJ1 += result;
    result = 1;
    result *= S1(Indexi[0], Indexj[0]);
    for (int k = 1; k < N-1; k++)
     result *= Identity(Indexi[k], Indexj[k]);
    result *= S0(Indexi[N-1], Indexj[N-1]);
    resultJ1 += result;
   }
   if ((resultJ0 != 0) || (resultJ1 != 0))
   {
    AOD.push_back(0.5*J0*T(resultJ0)+0.5*J1*T(resultJ1));
    IndexODi.push_back(i+1);
    IndexODj.push_back(j+1);
   }
  }
 }
// computing the nev lowest lying eigenstates:
 Matrix<T>::computeLowestEigenstates(AD, IndexD, AOD, IndexODi, IndexODj, n0, nev0, ncv, Evals, Evecs);
}

template<class T> void tJModel<T>::getTwoBodyHamiltonian(unsigned int position,
                                                         Matrix<T>& TwoBodyHamiltonian) const
{
#ifdef DEBUG
 if ((position > this->N-2) || (TwoBodyHamiltonian.getDim0() != 9) ||
     (TwoBodyHamiltonian.getDim1() != 9))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void tJModel<T>::" <<
          "getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian) const: " <<
          "((position > this->N-2) || (TwoBodyHamiltonian.getDim0() != 9) || " <<
          "(TwoBodyHamiltonian.getDim1() != 9))." << endl;
  exit(1);
 }
#endif
 TwoBodyHamiltonian.fillZeroes();
// the interactions:
// - the identity:
 Matrix<T> Identity(3, 3);
 Identity.fillZeroes();
 for (int i = 0; i < 3; i++)
 {
  Identity(i, i) = 1.0;
 }
// - particle number N =: N0:
 Matrix<T> N0(3, 3);
 N0.fillZeroes();
 N0(1, 1) = 1.0;
 N0(2, 2) = 1.0;
// - S^{+} =: S0:
 Matrix<T> S0(3, 3);
 S0.fillZeroes();
 S0(1, 2) = 1.0;
// - S^{-} =: S1:
 Matrix<T> S1(3, 3);
 S1.fillZeroes();
 S1(2, 1) = 1.0;
// - S^{z} =: S2:
 Matrix<T> S2(3, 3);
 S2.fillZeroes();
 S2(1, 1) = 0.5;
 S2(2, 2) = -0.5;
// - C_{\spinup}^{+} =: C0:
 Matrix<T> C0(3, 3);
 C0.fillZeroes();
 C0(1, 0) = 1.0;
// - C_{\spinup}^{-} =: C1:
 Matrix<T> C1(3, 3);
 C1.fillZeroes();
 C1(0, 1) = 1.0;
// - C_{\spindown}^{+} =: C2:
 Matrix<T> C2(3, 3);
 C2.fillZeroes();
 C2(2, 0) = 1.0;
// - C_{\spindown}^{-} =: C3:
 Matrix<T> C3(3, 3);
 C3.fillZeroes();
 C3(0, 2) = 1.0;
// - C := C_{\spinup}^{+} \otimes C_{\spinup}^{-} +
//        C_{\spinup}^{-} \otimes C_{\spinup}^{+} +
//        C_{\spindown}^{+} \otimes C_{\spindown}^{-} +
//        C_{\spindown}^{-} \otimes C_{\spindown}^{+}:
 Matrix<T> C(9, 9), X;
 C.fillZeroes();
 C0.multiplyDirectProduct(C1, X);
 C.add(X);
 C1.multiplyDirectProduct(C0, X);
 C.add(X);
 C2.multiplyDirectProduct(C3, X);
 C.add(X);
 C3.multiplyDirectProduct(C2, X);
 C.add(X);
// - S := 0.5 * (S^{+} \otimes S^{-} + S^{-} \otimes S^{+}) +
//        S^{z} \otimes S^{z} -
//        0.25 * N \otimes N:
 Matrix<T> S(9, 9);
 S.fillZeroes();
 S0.multiplyDirectProduct(S1, X);
 T element = 0.5;
 X.multiply(element);
 S.add(X);
 S1.multiplyDirectProduct(S0, X);
 X.multiply(element);
 S.add(X);
 S2.multiplyDirectProduct(S2, X);
 S.add(X);
 N0.multiplyDirectProduct(N0, X);
 element = -0.25;
 X.multiply(element);
 S.add(X);
// the parameters:
 T t0 = this->Parameters[0];
 T J0 = this->Parameters[1];
 T t1 = this->Parameters[2];
 T J1 = this->Parameters[3];
 T U = this->Parameters[4];
 T mu = this->Parameters[5];
 T V0 = this->Parameters[6];
 T l0 = this->Parameters[7];
 if (this->BC == "open")
 {
  if (position == 0)
  {
   N0.multiplyDirectProduct(Identity, X);
   element = V0 * pow(l0, 2);
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N0, X);
   element = T(0.5) * V0 * pow(l0-T(1.0), 2);
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   N0.multiplyDirectProduct(Identity, X);
   element = -mu;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N0, X);
   element = -T(0.5) * mu;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   N0.multiplyDirectProduct(N0, X);
   element = U;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   X = S;
   element = J0;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   X = C;
   element = -t0;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
  }
  else if ((position > 0) && (position < this->N-2))
  {
   N0.multiplyDirectProduct(Identity, X);
   element = T(0.5) * V0 * pow(T(position)-l0, 2);
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N0, X);
   element = T(0.5) * V0 * pow(T(position+1)-l0, 2);
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   N0.multiplyDirectProduct(Identity, X);
   element = -T(0.5) * mu;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N0, X);
   element = -T(0.5) * mu;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   if ((position % 2) == 0)
   {
    N0.multiplyDirectProduct(N0, X);
    element = U;
    X.multiply(element);
    TwoBodyHamiltonian.add(X);
    X = S;
    element = J0;
    X.multiply(element);
    TwoBodyHamiltonian.add(X);
    X = C;
    element = -t0;
    X.multiply(element);
    TwoBodyHamiltonian.add(X);
   }
   else if ((position % 2) == 1)
   {
    X = S;
    element = J1;
    X.multiply(element);
    TwoBodyHamiltonian.add(X);
    X = C;
    element = -t1;
    X.multiply(element);
    TwoBodyHamiltonian.add(X);
   }
  }
  else if (position == this->N-2)
  {
   N0.multiplyDirectProduct(Identity, X);
   element = T(0.5) * V0 * pow(T(this->N-2)-l0, 2);
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N0, X);
   element = V0 * pow(T(this->N-1)-l0, 2);
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   N0.multiplyDirectProduct(Identity, X);
   element = -T(0.5) * mu;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   Identity.multiplyDirectProduct(N0, X);
   element = -mu;
   X.multiply(element);
   TwoBodyHamiltonian.add(X);
   if (((this->N-2) % 2) == 0)
   {
    N0.multiplyDirectProduct(N0, X);
    element = U;
    X.multiply(element);
    TwoBodyHamiltonian.add(X);
    X = S;
    element = J0;
    X.multiply(element);
    TwoBodyHamiltonian.add(X);
    X = C;
    element = -t0;
    X.multiply(element);
    TwoBodyHamiltonian.add(X);
   }
   else if (((this->N-2) % 2) == 1)
   {
    X = S;
    element = J1;
    X.multiply(element);
    TwoBodyHamiltonian.add(X);
    X = C;
    element = -t1;
    X.multiply(element);
    TwoBodyHamiltonian.add(X);
   }
  }
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void tJModel<T>::" <<
          "getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian) const." << endl;
  exit(1);
 }
}
