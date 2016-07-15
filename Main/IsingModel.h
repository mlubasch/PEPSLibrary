/// Template class IsingModel implements quantum Ising model.
/** The template class IsingModel implements the quantum Ising model
    H = - B \sum_{i} sx_{i} -J \sum_{<i,j>} sz_{i} sz_{j}.
    \param Representation string, the representation, must be "interactions", "MPO" or "matrix"
    \param BC string, the boundary conditions, must be "open" or "periodic"
    \param N unsigned int, the number of spins
    \param d vector<unsigned int>, the physical dimensions of the lattice sites, is fixed to 2
             everywhere
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

template<class T> class IsingModel: public Hamiltonian<T>
{
 typedef T (*PointerToFunction)(double time);

 public:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  IsingModel();

/// Constructor for time-independent IsingModel with specific BC, N and Parameters.
/** This constructor initializes a time-independent IsingModel.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N0 input: unsigned int, the number of lattice sites
    \param Parameters0 input: const vector<T>&, the parameters, must have Parameters0.size()==2 */
  IsingModel(const string& BC0, unsigned int N0, const vector<T>& Parameters0);

/// Constructor for time-dependent IsingModel with specific BC, N, TimeFunctions and time.
/** This constructor initializes a time-dependent IsingModel.
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic"
    \param N0 input: unsigned int, the number of lattice sites
    \param TimeFunctions0 input: const vector<PointerToFunction>&, the time functions, must have
                                 TimeFunctions.size()==2
    \param time0 input: double, the time */
  IsingModel(const string& BC0, unsigned int N0,
                   const vector<PointerToFunction>& TimeFunctions0, double time0);

/// Standard copy constructor.
/** The standard copy constructor copies the input IsingModel into this.
    \param IsingModel0 input: const IsingModel<T>&, to be copied into this
    \sa IsingModel<T>& operator=(const IsingModel<T>& IsingModel0) */
  IsingModel(const IsingModel<T>& IsingModel0);

/// Standard destructor.
/** The standard destructor deletes the elements of IsingModel. */
  ~IsingModel();

/// Assigns IsingModel to this.
/** The operator= allows to assign IsingModel0 to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side IsingModel0.
    \param IsingModel0 input: const IsingModel<T>&, to be copied into this
    \return IsingModel<T>&, a reference to the new this
    \sa IsingModel(const IsingModel<T>& IsingModel0) */
  IsingModel<T>& operator=(const IsingModel<T>& IsingModel0);

/// Returns interactions representation of IsingModel.
/** This function returns the interactions of this IsingModel. The interactions
    are implemented as a vector of positions and a vector of matrices where each vector entry
    corresponds to one term of the sum making up the Hamiltonian.
    \param Positions output: vector< vector<unsigned int> >&, the positions of the matrices
    \param Matrices output: vector< vector< Matrix<T> > >&, the matrices */
  void getInteractions(vector< vector<unsigned int> >& Positions,
                       vector< vector< Matrix<T> > >& Matrices) const;

/// Returns MPO representation of IsingModel.
/** This function returns this IsingModel as a MPO.
    \param MPO0 output: MPO<T>&, the MPO representing this IsingModel */
  void getMPO(MPO<T>& MPO0) const;

/// Returns matrix representation of IsingModel.
/** This function returns this IsingModel as a matrix.
    \param Matrix0 output: Matrix<T>&, the matrix representing this IsingModel,
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

template<class T> IsingModel<T>::IsingModel()
{
 this->N = 0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> IsingModel<T>::IsingModel(const string& BC0, unsigned int N0,
                                            const vector<T>& Parameters0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (Parameters0.size() != 2))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> IsingModel<T>::" <<
          "IsingModel(const string& BC0, unsigned int N0, const vector<T>& Parameters0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (Parameters0.size() != 2))." << endl;
  exit(1);
 }
#endif
 this->Representation = "MPO";
 this->BC = BC0;
 this->N = N0;
 this->d = vector<unsigned int>(this->N);
// in the case of IsingModel d=2 everywhere:
 for (int i = 0; i < this->N; i++)
 {
  this->d[i] = 2;
 }
 this->Parameters = Parameters0;
 this->timeDependent = false;
 this->time = 0.0;
}

template<class T> IsingModel<T>::IsingModel(const string& BC0, unsigned int N0,
                                            const vector<PointerToFunction>& TimeFunctions0,
                                            double time0)
{
#ifdef DEBUG
 if (((BC0 != "open") && (BC0 != "periodic")) || (TimeFunctions0.size() != 2))
 {
  cerr << "Program terminated because of error in constructor " <<
          "template<class T> IsingModel<T>::" <<
          "IsingModel(const string& BC0, unsigned int N0, " <<
          "const vector<PointerToFunction>& TimeFunctions0, double time0): " <<
          "(((BC0 != open) && (BC0 != periodic)) || (TimeFunctions0.size() != 2))." << endl;
  exit(1);
 }
#endif
 this->Representation = "MPO";
 this->BC = BC0;
 this->N = N0;
 this->d = vector<unsigned int>(this->N);
// in the case of IsingModel d=2 everywhere:
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

template<class T> IsingModel<T>::IsingModel(const IsingModel<T>& IsingModel0)
{
 this->Representation = IsingModel0.Representation;
 this->BC = IsingModel0.BC;
 this->N = IsingModel0.N;
 this->d = IsingModel0.d;
 this->Parameters = IsingModel0.Parameters;
 this->timeDependent = IsingModel0.timeDependent;
 this->TimeFunctions = IsingModel0.TimeFunctions;
 this->time = IsingModel0.time;
}

template<class T> IsingModel<T>::~IsingModel() {}

template<class T> IsingModel<T>& IsingModel<T>::operator=(const IsingModel<T>& IsingModel0)
{
 if (this != &IsingModel0)
 {
  this->Representation = IsingModel0.Representation;
  this->BC = IsingModel0.BC;
  this->N = IsingModel0.N;
  this->d = IsingModel0.d;
  this->Parameters = IsingModel0.Parameters;
  this->timeDependent = IsingModel0.timeDependent;
  this->TimeFunctions = IsingModel0.TimeFunctions;
  this->time = IsingModel0.time;
 }
 return *this;
}

template<class T> void IsingModel<T>::getInteractions(vector< vector<unsigned int> >& Positions,
                                                      vector< vector< Matrix<T> > >& Matrices) const
{
 cerr << "The following function is not implemented yet: " <<
         "template<class T> void IsingModel<T>::" <<
         "getInteractions(vector< vector<unsigned int> >& Positions, " <<
         "vector< vector< Matrix<T> > >& Matrices)." << endl;
 exit(1);
}

template<class T> void IsingModel<T>::getMPO(MPO<T>& MPO0) const
{
 if (this->BC == "open")
 {
// the sigmas:
  vector<unsigned int> Shape(3);
  Shape[0] = 3; Shape[1] = 2; Shape[2] = 2;
  Tensor<T> Sigmas(Shape);
  Sigmas.fillZeroes();
  vector<unsigned int> Index(3);
  T element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  Sigmas.set(Index, element);
  Index[0] = 0; Index[1] = 1; Index[2] = 1;
  Sigmas.set(Index, element);
  Index[0] = 1; Index[1] = 1; Index[2] = 0;
  Sigmas.set(Index, element);
  Index[0] = 1; Index[1] = 0; Index[2] = 1;
  Sigmas.set(Index, element);
  Index[0] = 2; Index[1] = 0; Index[2] = 0;
  Sigmas.set(Index, element);
  element = -1.0;
  Index[0] = 2; Index[1] = 1; Index[2] = 1;
  Sigmas.set(Index, element);
  unsigned int D = 3;
  MPO0 = MPO<T>(this->BC, this->N, this->d, D);
// A for the left end:
  Shape[0] = 1; Shape[1] = 3; Shape[2] = 3;
  Tensor<T> A(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  element = -this->Parameters[0];
  Index[0] = 0; Index[1] = 2; Index[2] = 1;
  A.set(Index, element);
  element = -this->Parameters[1];
  Index[0] = 0; Index[1] = 1; Index[2] = 2;
  A.set(Index, element);
  Tensor<T> Sigmas0(Sigmas);
  vector<unsigned int> IndexA(1), Index0(1);
  IndexA[0] = 2; Index0[0] = 0;
  A.contract(IndexA, Sigmas0, Index0);
  unsigned int position = 0;
  MPO0.set(position, A);
// A for the intermediate lattice sites:
  Shape[0] = 3; Shape[1] = 3; Shape[2] = 3;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 0; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  Index[0] = 2; Index[1] = 2; Index[2] = 0;
  A.set(Index, element);
  element = -this->Parameters[0];
  Index[0] = 0; Index[1] = 2; Index[2] = 1;
  A.set(Index, element);
  element = -this->Parameters[1];
  Index[0] = 0; Index[1] = 1; Index[2] = 2;
  A.set(Index, element);
  element = 1.0;
  Index[0] = 1; Index[1] = 2; Index[2] = 2;
  A.set(Index, element);
  Sigmas0 = Sigmas;
  A.contract(IndexA, Sigmas0, Index0);
  for (position = 1; position < this->N-1; position++)
  {
   MPO0.set(position, A);
  }
// A for the right end:
  Shape[0] = 3; Shape[1] = 1; Shape[2] = 3;
  A = Tensor<T>(Shape);
  A.fillZeroes();
  element = 1.0;
  Index[0] = 2; Index[1] = 0; Index[2] = 0;
  A.set(Index, element);
  element = -this->Parameters[0];
  Index[0] = 0; Index[1] = 0; Index[2] = 1;
  A.set(Index, element);
  element = 1.0;
  Index[0] = 1; Index[1] = 0; Index[2] = 2;
  A.set(Index, element);
  Sigmas0 = Sigmas;
  A.contract(IndexA, Sigmas0, Index0);
  position = this->N-1;
  MPO0.set(position, A);
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void IsingModel<T>::" <<
          "getMPO(MPO<T>& MPO0) const." << endl;
  exit(1);
 }
}

template<class T> void IsingModel<T>::getMatrix(Matrix<T>& Matrix0) const
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
          "template<class T> void IsingModel<T>::" <<
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

template<class T> void IsingModel<T>::getTwoBodyHamiltonian(unsigned int position,
                                                            Matrix<T>& TwoBodyHamiltonian) const
{
#ifdef DEBUG
 if ((position >= this->N-1) || (TwoBodyHamiltonian.getDim0() != 4) ||
     (TwoBodyHamiltonian.getDim1() != 4))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void IsingModel<T>::" <<
          "getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian) const: " <<
          "((position >= this->N-1) || (TwoBodyHamiltonian.getDim0() != 4) || " <<
          "(TwoBodyHamiltonian.getDim1() != 4))." << endl;
  exit(1);
 }
#endif
 TwoBodyHamiltonian.fillZeroes();
 T B = this->Parameters[0];
 T J = this->Parameters[1];
 TwoBodyHamiltonian(0, 0) = -J;
 TwoBodyHamiltonian(1, 1) = J;
 TwoBodyHamiltonian(2, 2) = J;
 TwoBodyHamiltonian(3, 3) = -J;
 if (this->BC == "open")
 {
  if (position == 0)
  {
   TwoBodyHamiltonian(0, 1) = -0.5*B;
   TwoBodyHamiltonian(0, 2) = -B;
   TwoBodyHamiltonian(1, 0) = -0.5*B;
   TwoBodyHamiltonian(1, 3) = -B;
   TwoBodyHamiltonian(2, 0) = -B;
   TwoBodyHamiltonian(2, 3) = -0.5*B;
   TwoBodyHamiltonian(3, 1) = -B;
   TwoBodyHamiltonian(3, 2) = -0.5*B;
  }
  else if ((position >= 1) && (position <= this->N-3))
  {
   TwoBodyHamiltonian(0, 1) = -0.5*B;
   TwoBodyHamiltonian(0, 2) = -0.5*B;
   TwoBodyHamiltonian(1, 0) = -0.5*B;
   TwoBodyHamiltonian(1, 3) = -0.5*B;
   TwoBodyHamiltonian(2, 0) = -0.5*B;
   TwoBodyHamiltonian(2, 3) = -0.5*B;
   TwoBodyHamiltonian(3, 1) = -0.5*B;
   TwoBodyHamiltonian(3, 2) = -0.5*B;
  }
  else if (position == this->N-2)
  {
   TwoBodyHamiltonian(0, 1) = -B;
   TwoBodyHamiltonian(0, 2) = -0.5*B;
   TwoBodyHamiltonian(1, 0) = -B;
   TwoBodyHamiltonian(1, 3) = -0.5*B;
   TwoBodyHamiltonian(2, 0) = -0.5*B;
   TwoBodyHamiltonian(2, 3) = -B;
   TwoBodyHamiltonian(3, 1) = -0.5*B;
   TwoBodyHamiltonian(3, 2) = -B;
  }
 }
 else if (this->BC == "periodic")
 {
  cerr << "The following function is not implemented yet for periodic BC: " <<
          "template<class T> void IsingModel<T>::" <<
          "getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian) const." << endl;
  exit(1);
 }
}
