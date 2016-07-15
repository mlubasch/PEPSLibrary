/// Template class SpinlessFermions implements spinless fermions on a lattice.
/** The template class SpinlessFermions implements spinless fermions on a lattice:
    H = -t \sum_{l} (c_{l+1}^{+} c_{l} + c_{l}^{+} c_{l+1}) + U \sum_{l} n_{l} n_{l+1}
        -\mu \sum_{l} n_{l} + V_{0} \sum_{l} (l-l_{0})^{2} n_{l}
    \param Representation string, the representation, must be "interactions", "MPO" or "matrix"
    \param BC string, the boundary conditions, must be "open" or "periodic"
    \param N unsigned int, the number of lattice sites
    \param d vector<unsigned int>, the physical dimensions of the lattice sites, is fixed to 2
             everywhere
    \param Parameters vector<T>, the parameters {t, U, \mu, V_{0}, l_{0}}
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

template<class T> class SpinlessFermions: public Hamiltonian<T>
{
 typedef T (*PointerToFunction)(double time);

 public:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  SpinlessFermions();

/// Constructor for time-independent SpinlessFermions with specific BC, N and Parameters.
/** This constructor initializes a time-independent Hamiltonian of spinless fermions on a lattice.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N0 input: unsigned int, the number of lattice sites
    \param Parameters0 input: const vector<T>&, the parameters, must have Parameters0.size()==5 */
  SpinlessFermions(const string& BC0, unsigned int N0, const vector<T>& Parameters0);

/// Constructor for time-dependent SpinlessFermions with specific BC, N, TimeFunctions and time.
/** This constructor initializes a time-dependent Hamiltonian of spinless fermions on a lattice.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N0 input: unsigned int, the number of lattice sites
    \param TimeFunctions0 input: const vector<PointerToFunction>&, the time functions
    \param time0 input: double, the time */
  SpinlessFermions(const string& BC0, unsigned int N0,
                   const vector<PointerToFunction>& TimeFunctions0, double time0);

/// Standard copy constructor.
/** The standard copy constructor copies the input SpinlessFermions into this.
    \param SpinlessFermions0 input: const SpinlessFermions<T>&, to be copied into this
    \sa SpinlessFermions<T>& operator=(const SpinlessFermions<T>& SpinlessFermions0) */
  SpinlessFermions(const SpinlessFermions<T>& SpinlessFermions0);

/// Standard destructor.
/** The standard destructor deletes the elements of SpinlessFermions. */
  ~SpinlessFermions();

/// Assigns SpinlessFermions to this.
/** The operator= allows to assign SpinlessFermions0 to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side SpinlessFermions0.
    \param SpinlessFermions0 input: const SpinlessFermions<T>&, to be copied into this
    \return SpinlessFermions<T>&, a reference to the new this
    \sa SpinlessFermions(const SpinlessFermions<T>& SpinlessFermions0) */
  SpinlessFermions<T>& operator=(const SpinlessFermions<T>& SpinlessFermions0);

/// Returns interactions representation of SpinlessFermions.
/** This function returns the interactions of this SpinlessFermions. The interactions
    are implemented as a vector of positions and a vector of matrices where each vector entry
    corresponds to one term of the sum making up the Hamiltonian.
    \param Positions output: vector< vector<unsigned int> >&, the positions of the matrices
    \param Matrices output: vector< vector< Matrix<T> > >&, the matrices */
  void getInteractions(vector< vector<unsigned int> >& Positions,
                       vector< vector< Matrix<T> > >& Matrices) const;

/// Returns MPO representation of SpinlessFermions.
/** This function returns this SpinlessFermions as a MPO.
    \param MPO0 output: MPO<T>&, the MPO representing this SpinlessFermions */
  void getMPO(MPO<T>& MPO0) const;

/// Returns matrix representation of SpinlessFermions.
/** This function returns this SpinlessFermions as a matrix.
    \param Matrix0 output: Matrix<T>&, the matrix representing this SpinlessFermions,
                           must fulfill ((Matrix0.getDim0()==2^{N}) && (Matrix0.getDim1()==2^{N})),
                           will be of type "hermitian" */
  void getMatrix(Matrix<T>& Matrix0) const;

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

template<class T> SpinlessFermions<T>::SpinlessFermions()
{
 this->N = 0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> SpinlessFermions<T>::SpinlessFermions(const string& BC0, unsigned int N0,
                                                        const vector<T>& Parameters0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (Parameters0.size() != 5))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> SpinlessFermions<T>::" <<
          "SpinlessFermions(const string& BC0, unsigned int N0, const vector<T>& Parameters0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (Parameters0.size() != 5))." << endl;
  exit(1);
 }
#endif
 this->Representation = "MPO";
 this->BC = BC0;
 this->N = N0;
 this->d = vector<unsigned int>(this->N);
// in the case of spinless fermions d=2 everywhere:
 for (int i = 0; i < this->N; i++)
 {
  this->d[i] = 2;
 }
 this->Parameters = Parameters0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> SpinlessFermions<T>::SpinlessFermions(const string& BC0, unsigned int N0,
                                                        const vector<PointerToFunction>&
                                                        TimeFunctions0, double time0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (TimeFunctions0.size() != 5))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> SpinlessFermions<T>::" <<
          "SpinlessFermions(const string& BC0, unsigned int N0, " <<
          "const vector<PointerToFunction>& TimeFunctions0, double time0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (TimeFunctions0.size() != 5))." << endl;
  exit(1);
 }
#endif
 this->Representation = "MPO";
 this->BC = BC0;
 this->N = N0;
 this->d = vector<unsigned int>(this->N);
// in the case of spinless fermions d=2 everywhere:
 for (int i = 0; i < this->N; i++)
 {
  this->d[i] = 2;
 }
 this->timeDependent = true;
 this->TimeFunctions = TimeFunctions0;
 this->time = time0;
 for (int i = 0; i < this->Parameters.size(); i++)
 {
  this->Parameters[i] = this->TimeFunctions[i](this->time);
 }
}

template<class T> SpinlessFermions<T>::SpinlessFermions(const SpinlessFermions<T>&
                                                        SpinlessFermions0)
{
 this->Representation = SpinlessFermions0.Representation;
 this->BC = SpinlessFermions0.BC;
 this->N = SpinlessFermions0.N;
 this->d = SpinlessFermions0.d;
 this->Parameters = SpinlessFermions0.Parameters;
 this->timeDependent = SpinlessFermions0.timeDependent;
 this->TimeFunctions = SpinlessFermions0.TimeFunctions;
 this->time = SpinlessFermions0.time;
}

template<class T> SpinlessFermions<T>::~SpinlessFermions() {}

template<class T> SpinlessFermions<T>& SpinlessFermions<T>::operator=(const SpinlessFermions<T>& 
                                                                      SpinlessFermions0)
{
 if (this != &SpinlessFermions0)
 {
  this->Representation = SpinlessFermions0.Representation;
  this->BC = SpinlessFermions0.BC;
  this->N = SpinlessFermions0.N;
  this->d = SpinlessFermions0.d;
  this->Parameters = SpinlessFermions0.Parameters;
  this->timeDependent = SpinlessFermions0.timeDependent;
  this->TimeFunctions = SpinlessFermions0.TimeFunctions;
  this->time = SpinlessFermions0.time;
 }
 return *this;
}

template<class T> void SpinlessFermions<T>::getInteractions(vector< vector<unsigned int> >& 
                                                            Positions,
                                                            vector< vector< Matrix<T> > >& Matrices)
                                            const
{
 cerr << "The following function is not implemented yet: " <<
         "template<class T> void SpinlessFermions<T>::" <<
         "getInteractions(vector< vector<unsigned int> >& Positions, " <<
         "vector< vector< Matrix<T> > >& Matrices)." << endl;
 exit(1);
}

template<class T> void SpinlessFermions<T>::getMPO(MPO<T>& MPO0) const
{
 if (this->BC == "open")
 {
  unsigned int D = 5;
  MPO0 = MPO<T>(this->BC, this->N, this->d, D);
// the sigmas:
  vector<unsigned int> Shape(3);
  Shape[0] = 4; Shape[1] = 2; Shape[2] = 2;
  Tensor<T> Sigmas(Shape);
  Sigmas.fillZeroes();
  vector<unsigned int> Index(3);
  T element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  Sigmas.set(Index, element);
  Index[0] = 0; Index[1] = 1; Index[2] = 1;
  Sigmas.set(Index, element);
  Index[0] = 1; Index[1] = 1; Index[2] = 1;
  Sigmas.set(Index, element);
  Index[0] = 2; Index[1] = 0; Index[2] = 1;
  Sigmas.set(Index, element);
  Index[0] = 3; Index[1] = 1; Index[2] = 0;
  Sigmas.set(Index, element);
  Tensor<T> Sigmas0(Sigmas);
// A for the left end:
  Shape[0] = 1; Shape[1] = D; Shape[2] = 4;
  Tensor<T> A(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  element = this->Parameters[1];
  Index[0] = 0; Index[1] = 1; Index[2] = 1;
  A.set(Index, element);
  element = -this->Parameters[2] + this->Parameters[3]*pow(this->Parameters[4], 2);
  Index[0] = 0; Index[1] = 4; Index[2] = 1;
  A.set(Index, element);
  element = -this->Parameters[0];
  Index[0] = 0; Index[1] = 2; Index[2] = 2;
  A.set(Index, element);
  element = -this->Parameters[0];
  Index[0] = 0; Index[1] = 3; Index[2] = 3;
  A.set(Index, element);
  vector<unsigned int> IndexA(1), Index0(1);
  IndexA[0] = 2; Index0[0] = 0;
  A.contract(IndexA, Sigmas0, Index0);
  unsigned int position = 0;
  MPO0.set(position, A);
// A for the intermediate lattice sites:
  Shape[0] = D; Shape[1] = D; Shape[2] = 4;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 4; Index[1] = 4; Index[2] = 0;
  A.set(Index, element);
  element = this->Parameters[1];
  Index[0] = 0; Index[1] = 1; Index[2] = 1;
  A.set(Index, element);
  element = 1.0;
  Index[0] = 1; Index[1] = 4; Index[2] = 1;
  A.set(Index, element);
  element = -this->Parameters[0];
  Index[0] = 0; Index[1] = 2; Index[2] = 2;
  A.set(Index, element);
  element = 1.0;
  Index[0] = 3; Index[1] = 4; Index[2] = 2;
  A.set(Index, element);
  element = -this->Parameters[0];
  Index[0] = 0; Index[1] = 3; Index[2] = 3;
  A.set(Index, element);
  element = 1.0;
  Index[0] = 2; Index[1] = 4; Index[2] = 3;
  A.set(Index, element);
  Tensor<T> AC(A);
  Index[0] = 0; Index[1] = 4; Index[2] = 1;
  for (position = 1; position < this->N-1; position++)
  {
   element = this->Parameters[3]*pow(this->Parameters[4] - T(position), 2) - this->Parameters[2];
   A.set(Index, element);
   Sigmas0 = Sigmas;
   A.contract(IndexA, Sigmas0, Index0);
   MPO0.set(position, A);
   A = AC;
  }
// A for the right end:
  Shape[0] = D; Shape[1] = 1; Shape[2] = 4;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 4; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  element = -this->Parameters[2] + this->Parameters[3]*pow(this->Parameters[4]-double(this->N)+1.0, 2);
  Index[0] = 0; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  element = 1.0;
  Index[0] = 1; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  Index[0] = 3; Index[1] = 0; Index[2] = 2;
  A.set(Index, element);
  Index[0] = 2; Index[1] = 0; Index[2] = 3;
  A.set(Index, element);
  Sigmas0 = Sigmas;
  A.contract(IndexA, Sigmas0, Index0);
  position = this->N-1;
  MPO0.set(position, A);
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void SpinlessFermions<T>::" <<
          "getMPO(MPO<T>& MPO0) const." << endl;
  exit(1);
 }
}

template<class T> void SpinlessFermions<T>::getMatrix(Matrix<T>& Matrix0) const
{
#ifdef DEBUG
 unsigned int dim = 1;
 for (int i = 0; i < this->N; i++)
 {
  dim *= 2;
 }
 if ((Matrix0.getDim0() != dim) || (Matrix0.getDim1() != dim))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void SpinlessFermions<T>::" <<
          "getMatrix(Matrix<T>& Matrix0) const: " <<
          "((Matrix0.getDim0() != 2^{N}) || (Matrix0.getDim1() != 2^{N}))." << endl;
  exit(1);
 }
#endif
 MPO<T> MPO0;
 this->getMPO(MPO0);
 MPO0.getMatrix(Matrix0);
 Matrix0.setType("hermitian");
}

template<class T> void SpinlessFermions<T>::getTwoBodyHamiltonian(unsigned int position,
                                                                  Matrix<T>& TwoBodyHamiltonian) const
{
#ifdef DEBUG
 if ((position >= this->N-1) || (TwoBodyHamiltonian.getDim0() != 4) ||
     (TwoBodyHamiltonian.getDim1() != 4))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void SpinlessFermions<T>::" <<
          "getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian) const: " <<
          "((position >= this->N-1) || (TwoBodyHamiltonian.getDim0() != 4) || " <<
          "(TwoBodyHamiltonian.getDim1() != 4))." << endl;
  exit(1);
 }
#endif
 TwoBodyHamiltonian.fillZeroes();
 T t = this->Parameters[0];
 T U = this->Parameters[1];
 T mu = this->Parameters[2];
 T V0 = this->Parameters[3];
 T l0 = this->Parameters[4];
 T l = position;
 TwoBodyHamiltonian(1, 2) = -t;
 TwoBodyHamiltonian(2, 1) = -t;
 if (this->BC == "open")
 {
  if (position == 0)
  {
   TwoBodyHamiltonian(1, 1) = 0.5*V0*pow(l0-1.0, 2) - 0.5*mu;
   TwoBodyHamiltonian(2, 2) = V0*pow(l0, 2) - mu;
   TwoBodyHamiltonian(3, 3) = U - 1.5*mu + 1.5*V0*pow(l0, 2) - V0*l0 + 0.5*V0;
  }
  else if ((position >= 1) && (position <= this->N-3))
  {
   TwoBodyHamiltonian(1, 1) = 0.5*V0*pow(l+1.0-l0, 2) - 0.5*mu;
   TwoBodyHamiltonian(2, 2) = 0.5*V0*pow(l-l0, 2) - 0.5*mu;
   TwoBodyHamiltonian(3, 3) = U - mu + V0*pow(l-l0, 2) + V0*(l-l0) + 0.5*V0;
  }
  else if (position == this->N-2)
  {
   TwoBodyHamiltonian(1, 1) = V0*pow(this->N-1.0-l0, 2) - mu;
   TwoBodyHamiltonian(2, 2) = 0.5*V0*pow(this->N-2.0-l0, 2) - 0.5*mu;
   TwoBodyHamiltonian(3, 3) = U - 1.5*mu + 0.5*V0*pow(this->N-2.0-l0, 2) +
                              V0*pow(this->N-1.0-l0, 2);
  }
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void SpinlessFermions<T>::" <<
          "getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian) const." << endl;
  exit(1);
 }
}
