/// Abstract template class Hamiltonian implements interface for Hamiltonians.
/** The abstract template class Hamiltonian implements the interface for Hamiltonians and
    related operators in 1D.
    \param BC string, the boundary conditions
    \param N unsigned int, the number of spins
    \param d vector<unsigned int>, the physical dimensions of the spins
    \param Parameters vector<T>, the parameters
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

template<class T> class Hamiltonian : public Operator<T>
{
 typedef T (*PointerToFunction)(double time);

 public:

/// Sets boundary conditions of Hamiltonian.
/** This function sets the boundary conditions of this Hamiltonian. They must be
    either "open" or "periodic".
    \param BC0 input: const string&, the boundary conditions of this Hamiltonian,
                      must be "open" or "periodic" */
  inline void setBC(const string& BC0);

/// Returns boundary conditions of Hamiltonian.
/** The returned boundary conditions are either "open" or "periodic".
    \param BC0 output: string&, the boundary conditions of this Hamiltonian */
  inline void getBC(string& BC0) const;

/// Sets number of spins of Hamiltonian.
/** This function sets the number of spins of this Hamiltonian.
    \param N0 input: unsigned int, the number of spins of this Hamiltonian */
  inline void setN(unsigned int N0);

/// Returns number of spins of Hamiltonian.
/** This function returns the number of spins of this Hamiltonian.
    \return unsigned int, the number of spins of this Hamiltonian */
  inline unsigned int getN() const;

/// Sets physical dimensions of Hamiltonian.
/** This function sets the physical dimensions of this Hamiltonian.
    \param d0 input: const vector<unsigned int>&, the physical dimensions of this Hamiltonian */
  inline void setd(const vector<unsigned int>& d0);

/// Sets physical dimensions of Hamiltonian.
/** This function sets all the physical dimensions of this Hamiltonian to the same value.
    \param d0 input: unsigned int, the physical dimensions of this Hamiltonian */
  inline void setd(unsigned int d0);

/// Returns physical dimensions of Hamiltonian.
/** This function returns the physical dimensions of this Hamiltonian.
    \param d0 output: vector<unsigned int>&, the physical dimensions of this Hamiltonian */
  inline void getd(vector<unsigned int>& d0) const;

/// Sets parameters of Hamiltonian.
/** This function sets the parameters of this Hamiltonian.
    \param Parameters0 input: const vector<T>&, the parameters of this Hamiltonian */
  inline void setParameters(const vector<T>& Parameters0);

/// Returns parameters of Hamiltonian.
/** This function returns the parameters of this Hamiltonian.
    \param Parameters0 output: vector<T>&, the parameters of this Hamiltonian */
  inline void getParameters(vector<T>& Parameters0) const;

/// Sets Hamiltonian time-dependent.
/** This function sets this Hamiltonian time-dependent and passes the functions defining the
    time dependence.
    \param TimeFunctions0 input: const vector<PointerToFunction>&, the functions defining
                                 the time dependence of Parameters */
  inline void setTimeDependent(const vector<PointerToFunction>& TimeFunctions0);

/// Sets Hamiltonian time-independent.
/** This function sets this Hamiltonian time-independent. */
  inline void unsetTimeDependent();

/// Returns time dependence of Hamiltonian.
/** This function returns true if this Hamiltonian is time-dependent, it returns false if
    this Hamiltonian is time-independent. */
  inline bool isTimeDependent() const;

/// Sets functions defining time dependence of Hamiltonian.
/** This function sets the functions that define the time dependence of Parameters.
    \param TimeFunctions0 input: const vector<PointerToFunction>&, the time functions of
                                 this Hamiltonian */
  inline void setTimeFunctions(const vector<PointerToFunction>& TimeFunctions0);

/// Returns functions defining time dependence of Hamiltonian.
/** This function returns the functions that define the time dependence of Parameters.
    \param TimeFunctions0 output: vector<PointerToFunction>&, the time functions of
                                  this Hamiltonian */
  inline void getTimeFunctions(vector<PointerToFunction>& TimeFunctions0) const;

/// Sets time of Hamiltonian.
/** This function sets the time of this Hamiltonian and updates Parameters according to their
    time dependence defined in TimeFunctions.
    \param time0 input: double, the time of this Hamiltonian */
  inline void setTime(double time0);

/// Returns time of Hamiltonian.
/** This function returns the time of this Hamiltonian.
    \return double, the time of this Hamiltonian */
  inline double getTime() const;

/// Returns identity tensor at position.
/** This function returns the identity tensor at position in a specific Hamiltonian.
    \param position input: unsigned int, the position
    \param IdentityTensor output: Tensor<T>&, the IdentityTensor */
  void getIdentityTensor(unsigned int position, Tensor<T>& IdentityTensor) const;

/// Returns two body Hamiltonian for Trotter terms.
/** This function computes the hermitian part H of the exponent for the Trotter terms exp(-delta*H) for
    the even-odd Trotter decomposition and returns it as a Matrix TwoBodyHamiltonian. This function has
    to be implemented for Hamiltonian<T>::getTEMPOs to work.
    \param position input: unsigned int, the left position of the two body Hamiltonian in the
                           Hamiltonian sum, must be out of {0, 1, ..., N-2}
    \param TwoBodyHamiltonian output: Matrix<T>&, the two body Hamiltonian, must have the correct form
    \sa void getTEMPOs(const string& RealImaginaryTE, const string& TrotterOrder,
                       const string& TrotterDecomposition, double timeStep,
                       vector< MPO<T> >& TEMPOs, MPO<T>& LastFirstTEMPOsMerged) const */
  virtual void getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian) const = 0;

/// Returns MPOs for time evolution.
/** This function returns the MPOs resulting from a Trotter decomposition of the evolution operator. It
    needs getTwoBodyHamiltonian to work.
    \param RealImaginaryTE input: const string&, specifying if it is real or imaginary
                                  time evolution, must be "real" or "imaginary"
    \param TrotterOrder input: const string&, the Trotter order, must be "2nd" or "4th"
    \param TrotterDecomposition input: const string&, the Trotter decomposition,
                                       must be "even-odd" for the
                                       even-odd decomposition or "oneMPO" for both
                                       noncommuting parts put together into one MPO
    \param timeStep input: double, the time step
    \param TEMPOs output: vector< MPO<T> >&, the time evolution MPOs
    \param LastFirstTEMPOsMerged output: MPO<T>&, the MPO resulting from merging the
                                         last and the first TEMPO of
                                         the Trotter decomposition
    \sa virtual void getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian)
                                                                                         const = 0 */
  void getTEMPOs(const string& RealImaginaryTE, const string& TrotterOrder,
                 const string& TrotterDecomposition, double timeStep,
                 vector< MPO<T> >& TEMPOs, MPO<T>& LastFirstTEMPOsMerged);

/// Returns Matrices for time evolution.
/** This function returns the Matrices resulting from a Trotter decomposition of the evolution
    operator. It needs getTwoBodyHamiltonian to work. This function is implemented only for debugging
    purposes of void Hamiltonian<T>::getTEMPOs.
    \param RealImaginaryTE input: const string&, specifying if it is real or imaginary
                                  time evolution, must be "real" or "imaginary"
    \param TrotterOrder input: const string&, the Trotter order, must be "2nd" or "4th"
    \param TrotterDecomposition input: const string&, the Trotter decomposition,
                                       must be "even-odd" for the
                                       even-odd decomposition or "oneMPO" for both
                                       noncommuting parts put together into one MPO
    \param timeStep input: double, the time step
    \param TEMatrices output: vector< Matrix<T> >&, the time evolution Matrices
    \param LastFirstTEMatricesMerged output: Matrix<T>&, the TEMatrix resulting from merging the
                                             last and the first TEMatrix of
                                             the Trotter decomposition
    \sa virtual void getTwoBodyHamiltonian(unsigned int position, Matrix<T>& TwoBodyHamiltonian)
                                                                                         const = 0 */
  void getTEMatrices(const string& RealImaginaryTE, const string& TrotterOrder,
                     const string& TrotterDecomposition, double timeStep,
                     vector< Matrix<T> >& TEMatrices, Matrix<T>& LastFirstTEMatricesMerged);

 protected:

/// Boundary conditions.
 string BC;

/// Number of spins.
 unsigned int N;

/// Physical dimensions of spins.
 vector<unsigned int> d;

/// Parameters.
 vector<T> Parameters;

/// Time dependence.
 bool timeDependent;

/// Time dependence of Parameters.
 vector<PointerToFunction> TimeFunctions;

/// Time.
 double time;
};

template<class T> inline void Hamiltonian<T>::setBC(const string& BC0)
{
#ifdef DEBUG
 if ((BC0 != "open") && (BC0 != "periodic"))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline void Hamiltonian<T>::" <<
          "setBC(const string& BC0)" << endl;
  exit(1);
 }
#endif
 this->BC = BC0;
}

template<class T> inline void Hamiltonian<T>::getBC(string& BC0) const
{
 BC0 = this->BC;
}

template<class T> inline void Hamiltonian<T>::setN(unsigned int N0)
{
 this->N = N0;
}

template<class T> inline unsigned int Hamiltonian<T>::getN() const
{
 return this->N;
}

template<class T> inline void Hamiltonian<T>::setd(const vector<unsigned int>& d0)
{
#ifdef DEBUG
 if (d0.size() != this->N)
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline void Hamiltonian<T>::" <<
          "setd(const vector<unsigned int>& d0)" << endl;
  exit(1);
 }
#endif
 this->d = d0;
}

template<class T> inline void Hamiltonian<T>::setd(unsigned int d0)
{
 this->d = vector<unsigned int>(this->N);
 for (int i = 0; i < this->N; i++)
 {
  this->d[i] = d0;
 }
}

template<class T> inline void Hamiltonian<T>::getd(vector<unsigned int>& d0) const
{
 d0 = this->d;
}

template<class T> inline void Hamiltonian<T>::setParameters(const vector<T>& Parameters0)
{
#ifdef DEBUG
 if (Parameters0.size() != this->Parameters.size())
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline void Hamiltonian<T>::" <<
          "setParameters(const vector<T>& Parameters0)" << endl;
  exit(1);
 }
#endif
 this->Parameters = Parameters0;
}

template<class T> inline void Hamiltonian<T>::getParameters(vector<T>& Parameters0) const
{
 Parameters0 = this->Parameters;
}

template<class T> inline void Hamiltonian<T>::setTimeDependent(const vector<PointerToFunction>&
                                                               TimeFunctions0)
{
#ifdef DEBUG
 if (TimeFunctions0.size() != this->Parameters.size())
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline void Hamiltonian<T>::" <<
          "setTimeDependent(const vector<PointerToFunction>& TimeFunctions0)" << endl;
  exit(1);
 }
#endif
 this->timeDependent = true;
 this->TimeFunctions = TimeFunctions0;
}

template<class T> inline void Hamiltonian<T>::unsetTimeDependent()
{
 this->timeDependent = false;
 for (int i = 0; i < this->TimeFunctions.size(); i++)
 {
  this->TimeFunctions[i] = 0;
 }
}

template<class T> inline bool Hamiltonian<T>::isTimeDependent() const
{
 return this->timeDependent;
}

template<class T> inline void Hamiltonian<T>::setTimeFunctions(const vector<PointerToFunction>&
                                                               TimeFunctions0)
{
#ifdef DEBUG
 if (TimeFunctions0.size() != this->Parameters.size())
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline void Hamiltonian<T>::" <<
          "setTimeFunctions(const vector<PointerToFunction>& TimeFunctions0)" << endl;
  exit(1);
 }
#endif
 this->timeDependent = true;
 this->TimeFunctions = TimeFunctions0;
}

template<class T> inline void Hamiltonian<T>::getTimeFunctions(vector<PointerToFunction>&
                                                               TimeFunctions0) const
{
 TimeFunctions0 = this->TimeFunctions;
}

template<class T> inline void Hamiltonian<T>::setTime(double time0)
{
#ifdef DEBUG
 if (this->timeDependent == false)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> inline void Hamiltonian<T>::" <<
          "setTime(double time0): " <<
          "(this->timeDependent == false)" << endl;
  exit(1);
 }
#endif
 this->time = time0;
 for (int i = 0; i < this->Parameters.size(); i++)
 {
  this->Parameters[i] = this->TimeFunctions[i](time0);
 }
}

template<class T> inline double Hamiltonian<T>::getTime() const
{
 return this->time;
}

template<class T> void Hamiltonian<T>::getIdentityTensor(unsigned int position,
                                                         Tensor<T>& IdentityTensor) const
{
#ifdef DEBUG
 if (position >= this->N)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Hamiltonian<T>::" <<
          "getIdentityTensor(unsigned int position, Tensor<T>& IdentityTensor) const: " <<
          "(position >= this->N)." << endl;
  exit(1);
 }
#endif
 vector<unsigned int> Shape(4), Index(4);
 Shape[0] = 1; Shape[1] = 1; Shape[2] = this->d[position]; Shape[3] = this->d[position];
 IdentityTensor = Tensor<T>(Shape);
 IdentityTensor.fillZeroes();
 Index[0] = 0; Index[1] = 0;
 T element = 1.0;
 for (int i = 0; i < this->d[position]; i++)
 {
  Index[2] = i; Index[3] = i;
  IdentityTensor.set(Index, element);
 }
}

template<class T> void Hamiltonian<T>::getTEMPOs(const string& RealImaginaryTE,
                                                 const string& TrotterOrder,
                                                 const string& TrotterDecomposition,
                                                 double timeStep, vector< MPO<T> >& TEMPOs,
                                                 MPO<T>& LastFirstTEMPOsMerged)
{
#ifdef DEBUG
 if ((this->N == 0) || ((RealImaginaryTE != "real") && (RealImaginaryTE != "imaginary")) ||
     ((TrotterOrder != "2nd") && (TrotterOrder != "4th")) ||
     ((TrotterDecomposition != "even-odd") && (TrotterDecomposition != "oneMPO")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Hamiltonian<T>::" <<
          "getTEMPOs(const string& RealImaginaryTE, const string& TrotterOrder, " <<
                    "const string& TrotterDecomposition, double timeStep, vector< MPO<T> >& TEMPOs, " <<
                    "MPO<T>& LastFirstTEMPOsMerged) const: " <<
          "The input arguments are incorrect." << endl;
  exit(1);
 }
#endif
 Matrix<T> TwoBodyHamiltonian;
 unsigned int d0, d1, numSV;
// we compute DEven of TEMPOEven and DOdd of TEMPOOdd:
 unsigned int DEven = 0;
 for (int position = 0; position <= this->N-2; position += 2)
 {
  d0 = this->d[position]; d1 = this->d[position+1]; numSV = min(d0*d0, d1*d1);
  if (numSV > DEven)
   DEven = numSV;
 }
 unsigned int DOdd = 0;
 for (int position = 1; position <= this->N-2; position += 2)
 {
  d0 = this->d[position]; d1 = this->d[position+1]; numSV = min(d0*d0, d1*d1);
  if (numSV > DOdd)
   DOdd = numSV;
 }
 MPO<T> TEMPOEven(this->BC, this->N, this->d, DEven), TEMPOOdd(this->BC, this->N, this->d, DOdd);
 Tensor<T> TensorLeft, TensorRight, IdentityTensor;
 vector<unsigned int> Shape(4), Index(4);
 T deltaTime, delta;
 if (RealImaginaryTE == "real")
 {
  deltaTime = T(0.0, 1.0) * T(timeStep);
 }
 else if (RealImaginaryTE == "imaginary")
 {
  deltaTime = T(timeStep);
 }
 if ((TrotterOrder == "2nd") && (TrotterDecomposition == "even-odd"))
 {
  TEMPOs = vector< MPO<T> >(3);
// in the case of real time evolution with a time-dependent Hamiltonian we put
// this->time = this->time + 0.5*timeStep:
  if ((RealImaginaryTE == "real") && (this->timeDependent))
  {
   this->setTime(this->time + 0.5*timeStep);
  }
// 0. MPO for position 0 and 2 in TEMPOs: TEMPOEven=exp(-0.5*deltaTime*H_{EVEN})
  delta = T(0.5) * deltaTime;
  for (int position = 0; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1]; numSV = min(d0*d0, d1*d1);
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   TwoBodyHamiltonian.setType("hermitian");
   Shape[0] = 1; Shape[1] = numSV; Shape[2] = d0; Shape[3] = d0;
   TensorLeft = Tensor<T>(Shape);
   Shape[0] = numSV; Shape[1] = 1; Shape[2] = d1; Shape[3] = d1;
   TensorRight = Tensor<T>(Shape);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   TwoBodyHamiltonian.twoBodyExponentSVD(delta, d0, d1, TensorLeft, TensorRight);
   TEMPOEven.set(position, TensorLeft);
   TEMPOEven.set(position+1, TensorRight);
  }
  if (this->N % 2 == 1)
  {
   this->getIdentityTensor(this->N-1, IdentityTensor);
   TEMPOEven.set(this->N-1, IdentityTensor);
  }
  TEMPOs[0] = TEMPOEven;
  TEMPOs[2] = TEMPOEven;
// 1. MPO for position 1 in TEMPOs: TEMPOOdd=exp(-deltaTime*H_{ODD})
  delta = deltaTime;
  this->getIdentityTensor(0, IdentityTensor);
  TEMPOOdd.set(0, IdentityTensor);
  for (int position = 1; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1]; numSV = min(d0*d0, d1*d1);
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   TwoBodyHamiltonian.setType("hermitian");
   Shape[0] = 1; Shape[1] = numSV; Shape[2] = d0; Shape[3] = d0;
   TensorLeft = Tensor<T>(Shape);
   Shape[0] = numSV; Shape[1] = 1; Shape[2] = d1; Shape[3] = d1;
   TensorRight = Tensor<T>(Shape);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   TwoBodyHamiltonian.twoBodyExponentSVD(delta, d0, d1, TensorLeft, TensorRight);
   TEMPOOdd.set(position, TensorLeft);
   TEMPOOdd.set(position+1, TensorRight);
  }
  if (this->N % 2 == 0)
  {
   this->getIdentityTensor(this->N-1, IdentityTensor);
   TEMPOOdd.set(this->N-1, IdentityTensor);
  }
  TEMPOs[1] = TEMPOOdd;
// 2. LastFirstTEMPOsMerged: TEMPOEven=exp(-deltaTime*H_{EVEN})
  delta = deltaTime;
  for (int position = 0; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1]; numSV = min(d0*d0, d1*d1);
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   TwoBodyHamiltonian.setType("hermitian");
   Shape[0] = 1; Shape[1] = numSV; Shape[2] = d0; Shape[3] = d0;
   TensorLeft = Tensor<T>(Shape);
   Shape[0] = numSV; Shape[1] = 1; Shape[2] = d1; Shape[3] = d1;
   TensorRight = Tensor<T>(Shape);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   TwoBodyHamiltonian.twoBodyExponentSVD(delta, d0, d1, TensorLeft, TensorRight);
   TEMPOEven.set(position, TensorLeft);
   TEMPOEven.set(position+1, TensorRight);
  }
  if (this->N % 2 == 1)
  {
   this->getIdentityTensor(this->N-1, IdentityTensor);
   TEMPOEven.set(this->N-1, IdentityTensor);
  }
  LastFirstTEMPOsMerged = TEMPOEven;
// in the case of real time evolution with a time-dependent Hamiltonian we undo the time change, so
// this->time = this->time - 0.5*timeStep:
  if ((RealImaginaryTE == "real") && (this->timeDependent))
  {
   this->setTime(this->time - 0.5*timeStep);
  }
 }
 else if ((TrotterOrder == "4th") && (TrotterDecomposition == "even-odd") && (!this->timeDependent))
 {
  TEMPOs = vector< MPO<T> >(5);
  T p;
// 0. MPO for position 0 in TEMPOs: TEMPOEven=exp(-p1*deltaTime*H_{EVEN})
//    (T. Prosen & I. Pizorn, J. Phys. A: Math. Gen. 39 (2006), 5957-5964):
  p = T(0.25, sqrt(3.0)/12.0);
  delta = p*deltaTime;
  for (int position = 0; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1]; numSV = min(d0*d0, d1*d1);
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   TwoBodyHamiltonian.setType("hermitian");
   Shape[0] = 1; Shape[1] = numSV; Shape[2] = d0; Shape[3] = d0;
   TensorLeft = Tensor<T>(Shape);
   Shape[0] = numSV; Shape[1] = 1; Shape[2] = d1; Shape[3] = d1;
   TensorRight = Tensor<T>(Shape);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   TwoBodyHamiltonian.twoBodyExponentSVD(delta, d0, d1, TensorLeft, TensorRight);
   TEMPOEven.set(position, TensorLeft);
   TEMPOEven.set(position+1, TensorRight);
  }
  if (this->N % 2 == 1)
  {
   this->getIdentityTensor(this->N-1, IdentityTensor);
   TEMPOEven.set(this->N-1, IdentityTensor);
  }
  TEMPOs[0] = TEMPOEven;
// 1. MPO for position 1 in TEMPOs: TEMPOOdd=exp(-p2*deltaTime*H_{ODD})
//    (T. Prosen & I. Pizorn, J. Phys. A: Math. Gen. 39 (2006), 5957-5964):
  p = T(0.5, sqrt(3.0)/6.0);
  delta = p*deltaTime;
  this->getIdentityTensor(0, IdentityTensor);
  TEMPOOdd.set(0, IdentityTensor);
  for (int position = 1; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1]; numSV = min(d0*d0, d1*d1);
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   TwoBodyHamiltonian.setType("hermitian");
   Shape[0] = 1; Shape[1] = numSV; Shape[2] = d0; Shape[3] = d0;
   TensorLeft = Tensor<T>(Shape);
   Shape[0] = numSV; Shape[1] = 1; Shape[2] = d1; Shape[3] = d1;
   TensorRight = Tensor<T>(Shape);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   TwoBodyHamiltonian.twoBodyExponentSVD(delta, d0, d1, TensorLeft, TensorRight);
   TEMPOOdd.set(position, TensorLeft);
   TEMPOOdd.set(position+1, TensorRight);
  }
  if (this->N % 2 == 0)
  {
   this->getIdentityTensor(this->N-1, IdentityTensor);
   TEMPOOdd.set(this->N-1, IdentityTensor);
  }
  TEMPOs[1] = TEMPOOdd;
// 2. MPO for position 2 in TEMPOs: TEMPOEven=exp(-p3*deltaTime*H_{EVEN})
//    (T. Prosen & I. Pizorn, J. Phys. A: Math. Gen. 39 (2006), 5957-5964):
  p = 0.5;
  delta = p*deltaTime;
  for (int position = 0; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1]; numSV = min(d0*d0, d1*d1);
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   TwoBodyHamiltonian.setType("hermitian");
   Shape[0] = 1; Shape[1] = numSV; Shape[2] = d0; Shape[3] = d0;
   TensorLeft = Tensor<T>(Shape);
   Shape[0] = numSV; Shape[1] = 1; Shape[2] = d1; Shape[3] = d1;
   TensorRight = Tensor<T>(Shape);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   TwoBodyHamiltonian.twoBodyExponentSVD(delta, d0, d1, TensorLeft, TensorRight);
   TEMPOEven.set(position, TensorLeft);
   TEMPOEven.set(position+1, TensorRight);
  }
  if (this->N % 2 == 1)
  {
   this->getIdentityTensor(this->N-1, IdentityTensor);
   TEMPOEven.set(this->N-1, IdentityTensor);
  }
  TEMPOs[2] = TEMPOEven;
// 3. MPO for position 3 in TEMPOs: TEMPOOdd=exp(-p4*deltaTime*H_{ODD})
//    (T. Prosen & I. Pizorn, J. Phys. A: Math. Gen. 39 (2006), 5957-5964):
  p = T(0.5, -sqrt(3.0)/6.0);
  delta = p*deltaTime;
  this->getIdentityTensor(0, IdentityTensor);
  TEMPOOdd.set(0, IdentityTensor);
  for (int position = 1; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1]; numSV = min(d0*d0, d1*d1);
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   TwoBodyHamiltonian.setType("hermitian");
   Shape[0] = 1; Shape[1] = numSV; Shape[2] = d0; Shape[3] = d0;
   TensorLeft = Tensor<T>(Shape);
   Shape[0] = numSV; Shape[1] = 1; Shape[2] = d1; Shape[3] = d1;
   TensorRight = Tensor<T>(Shape);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   TwoBodyHamiltonian.twoBodyExponentSVD(delta, d0, d1, TensorLeft, TensorRight);
   TEMPOOdd.set(position, TensorLeft);
   TEMPOOdd.set(position+1, TensorRight);
  }
  if (this->N % 2 == 0)
  {
   this->getIdentityTensor(this->N-1, IdentityTensor);
   TEMPOOdd.set(this->N-1, IdentityTensor);
  }
  TEMPOs[3] = TEMPOOdd;
// 4. MPO for position 4 in TEMPOs: TEMPOEven=exp(-p5*deltaTime*H_{EVEN})
//    (T. Prosen & I. Pizorn, J. Phys. A: Math. Gen. 39 (2006), 5957-5964):
  p = T(0.25, -sqrt(3.0)/12.0);
  delta = p*deltaTime;
  for (int position = 0; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1]; numSV = min(d0*d0, d1*d1);
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   TwoBodyHamiltonian.setType("hermitian");
   Shape[0] = 1; Shape[1] = numSV; Shape[2] = d0; Shape[3] = d0;
   TensorLeft = Tensor<T>(Shape);
   Shape[0] = numSV; Shape[1] = 1; Shape[2] = d1; Shape[3] = d1;
   TensorRight = Tensor<T>(Shape);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   TwoBodyHamiltonian.twoBodyExponentSVD(delta, d0, d1, TensorLeft, TensorRight);
   TEMPOEven.set(position, TensorLeft);
   TEMPOEven.set(position+1, TensorRight);
  }
  if (this->N % 2 == 1)
  {
   this->getIdentityTensor(this->N-1, IdentityTensor);
   TEMPOEven.set(this->N-1, IdentityTensor);
  }
  TEMPOs[4] = TEMPOEven;
// 5. LastFirstTEMPOsMerged: TEMPOEven=exp(-(p5+p1)*deltaTime*H_{EVEN})
//    (T. Prosen & I. Pizorn, J. Phys. A: Math. Gen. 39 (2006), 5957-5964):
  p = 0.5;
  delta = p*deltaTime;
  for (int position = 0; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1]; numSV = min(d0*d0, d1*d1);
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   TwoBodyHamiltonian.setType("hermitian");
   Shape[0] = 1; Shape[1] = numSV; Shape[2] = d0; Shape[3] = d0;
   TensorLeft = Tensor<T>(Shape);
   Shape[0] = numSV; Shape[1] = 1; Shape[2] = d1; Shape[3] = d1;
   TensorRight = Tensor<T>(Shape);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   TwoBodyHamiltonian.twoBodyExponentSVD(delta, d0, d1, TensorLeft, TensorRight);
   TEMPOEven.set(position, TensorLeft);
   TEMPOEven.set(position+1, TensorRight);
  }
  if (this->N % 2 == 1)
  {
   this->getIdentityTensor(this->N-1, IdentityTensor);
   TEMPOEven.set(this->N-1, IdentityTensor);
  }
  LastFirstTEMPOsMerged = TEMPOEven;
 }
}

template<class T> void Hamiltonian<T>::getTEMatrices(const string& RealImaginaryTE,
                                                     const string& TrotterOrder,
                                                     const string& TrotterDecomposition,
                                                     double timeStep, vector< Matrix<T> >& TEMatrices,
                                                     Matrix<T>& LastFirstTEMatricesMerged)
{
#ifdef DEBUG
 if ((this->N == 0) || ((RealImaginaryTE != "real") && (RealImaginaryTE != "imaginary")) ||
     ((TrotterOrder != "2nd") && (TrotterOrder != "4th")) ||
     ((TrotterDecomposition != "even-odd") && (TrotterDecomposition != "oneMPO")))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Hamiltonian<T>::" <<
          "getTEMatrices(const string& RealImaginaryTE, const string& TrotterOrder, " <<
                        "const string& TrotterDecomposition, double timeStep, " <<
                        "vector< Matrix<T> >& TEMatrices, " <<
                        "Matrix<T>& LastFirstTEMatricesMerged) const: " <<
          "The input arguments are incorrect." << endl;
  exit(1);
 }
#endif
 unsigned int dim;
 vector< Matrix<T> > Identities(this->N); Matrix<T> Identity;
 for (int i = 0; i < this->N; i++)
 {
  dim = this->d[i];
  Identity = Matrix<T>(dim, dim);
  Identity.fillZeroes();
  for (int j = 0; j < dim; j++)
   Identity(j, j) = T(1.0);
  Identities[i] = Identity;
 }
 dim = 1;
 for (int i = 0; i < this->N; i++)
  dim *= this->d[i];
 Matrix<T> TwoBodyHamiltonian, X, Y, Z(dim, dim);
 Z.setType("hermitian");
 unsigned int d0, d1;
 T deltaTime, delta;
 if (RealImaginaryTE == "real")
 {
  deltaTime = T(0.0, 1.0) * T(timeStep);
 }
 else if (RealImaginaryTE == "imaginary")
 {
  deltaTime = T(timeStep);
 }
 if ((TrotterOrder == "2nd") && (TrotterDecomposition == "even-odd"))
 {
  TEMatrices = vector< Matrix<T> >(3);
// in the case of real time evolution with a time-dependent Hamiltonian we put
// this->time = this->time + 0.5*timeStep:
  if ((RealImaginaryTE == "real") && (this->timeDependent))
  {
   this->setTime(this->time + 0.5*timeStep);
  }
// 0. Matrix for position 0 and 2 in TEMatrices: X=exp(-0.5*deltaTime*H_{EVEN})
  delta = T(0.5) * deltaTime;
  Z.fillZeroes();
  for (int position = 0; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1];
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   X = TwoBodyHamiltonian;
   for (int i = position; i < this->N-2; i++)
   {
    Identity = Identities[i+2];
    X.multiplyDirectProduct(Identity, Y);
    X = Y;
   }
   for (int i = position; i > 0; i--)
   {
    Identity = Identities[i-1];
    Identity.multiplyDirectProduct(X, Y);
    X = Y;
   }
   Z.add(X);
  }
  Z.exponentialize(-delta, X);
  TEMatrices[0] = X;
  TEMatrices[2] = X;
// 1. Matrix for position 1 in TEMatrices: X=exp(-deltaTime*H_{ODD})
  delta = deltaTime;
  Z.fillZeroes();
  for (int position = 1; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1];
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   X = TwoBodyHamiltonian;
   for (int i = position; i < this->N-2; i++)
   {
    Identity = Identities[i+2];
    X.multiplyDirectProduct(Identity, Y);
    X = Y;
   }
   for (int i = position; i > 0; i--)
   {
    Identity = Identities[i-1];
    Identity.multiplyDirectProduct(X, Y);
    X = Y;
   }
   Z.add(X);
  }
  Z.exponentialize(-delta, X);
  TEMatrices[1] = X;
// 2. LastFirstTEMatricesMerged: X=exp(-deltaTime*H_{EVEN})
  delta = deltaTime;
  Z.fillZeroes();
  for (int position = 0; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1];
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   X = TwoBodyHamiltonian;
   for (int i = position; i < this->N-2; i++)
   {
    Identity = Identities[i+2];
    X.multiplyDirectProduct(Identity, Y);
    X = Y;
   }
   for (int i = position; i > 0; i--)
   {
    Identity = Identities[i-1];
    Identity.multiplyDirectProduct(X, Y);
    X = Y;
   }
   Z.add(X);
  }
  Z.exponentialize(-delta, X);
  LastFirstTEMatricesMerged = X;
// in the case of real time evolution with a time-dependent Hamiltonian we undo the time change, so
// this->time = this->time - 0.5*timeStep:
  if ((RealImaginaryTE == "real") && (this->timeDependent))
  {
   this->setTime(this->time - 0.5*timeStep);
  }
 }
 else if ((TrotterOrder == "4th") && (TrotterDecomposition == "even-odd") && (!this->timeDependent))
 {
  TEMatrices = vector< Matrix<T> >(5);
  T p;
// 0. Matrix for position 0 in TEMatrices: X=exp(-p1*deltaTime*H_{EVEN})
//    (T. Prosen & I. Pizorn, J. Phys. A: Math. Gen. 39 (2006), 5957-5964):
  p = T(0.25, sqrt(3.0)/12.0);
  delta = p*deltaTime;
  Z.fillZeroes();
  for (int position = 0; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1];
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   X = TwoBodyHamiltonian;
   for (int i = position; i < this->N-2; i++)
   {
    Identity = Identities[i+2];
    X.multiplyDirectProduct(Identity, Y);
    X = Y;
   }
   for (int i = position; i > 0; i--)
   {
    Identity = Identities[i-1];
    Identity.multiplyDirectProduct(X, Y);
    X = Y;
   }
   Z.add(X);
  }
  Z.exponentialize(-delta, X);
  TEMatrices[0] = X;
// 1. Matrix for position 1 in TEMatrices: X=exp(-p2*deltaTime*H_{ODD})
//    (T. Prosen & I. Pizorn, J. Phys. A: Math. Gen. 39 (2006), 5957-5964):
  p = T(0.5, sqrt(3.0)/6.0);
  delta = p*deltaTime;
  Z.fillZeroes();
  for (int position = 1; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1];
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   X = TwoBodyHamiltonian;
   for (int i = position; i < this->N-2; i++)
   {
    Identity = Identities[i+2];
    X.multiplyDirectProduct(Identity, Y);
    X = Y;
   }
   for (int i = position; i > 0; i--)
   {
    Identity = Identities[i-1];
    Identity.multiplyDirectProduct(X, Y);
    X = Y;
   }
   Z.add(X);
  }
  Z.exponentialize(-delta, X);
  TEMatrices[1] = X;
// 2. Matrix for position 2 in TEMatrices: X=exp(-p3*deltaTime*H_{EVEN})
//    (T. Prosen & I. Pizorn, J. Phys. A: Math. Gen. 39 (2006), 5957-5964):
  p = 0.5;
  delta = p*deltaTime;
  Z.fillZeroes();
  for (int position = 0; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1];
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   X = TwoBodyHamiltonian;
   for (int i = position; i < this->N-2; i++)
   {
    Identity = Identities[i+2];
    X.multiplyDirectProduct(Identity, Y);
    X = Y;
   }
   for (int i = position; i > 0; i--)
   {
    Identity = Identities[i-1];
    Identity.multiplyDirectProduct(X, Y);
    X = Y;
   }
   Z.add(X);
  }
  Z.exponentialize(-delta, X);
  TEMatrices[2] = X;
// 3. Matrix for position 3 in TEMatrices: X=exp(-p4*deltaTime*H_{ODD})
//    (T. Prosen & I. Pizorn, J. Phys. A: Math. Gen. 39 (2006), 5957-5964):
  p = T(0.5, -sqrt(3.0)/6.0);
  delta = p*deltaTime;
  Z.fillZeroes();
  for (int position = 1; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1];
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   X = TwoBodyHamiltonian;
   for (int i = position; i < this->N-2; i++)
   {
    Identity = Identities[i+2];
    X.multiplyDirectProduct(Identity, Y);
    X = Y;
   }
   for (int i = position; i > 0; i--)
   {
    Identity = Identities[i-1];
    Identity.multiplyDirectProduct(X, Y);
    X = Y;
   }
   Z.add(X);
  }
  Z.exponentialize(-delta, X);
  TEMatrices[3] = X;
// 4. Matrix for position 4 in TEMatrices: X=exp(-p5*deltaTime*H_{EVEN})
//    (T. Prosen & I. Pizorn, J. Phys. A: Math. Gen. 39 (2006), 5957-5964):
  p = T(0.25, -sqrt(3.0)/12.0);
  delta = p*deltaTime;
  Z.fillZeroes();
  for (int position = 0; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1];
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   X = TwoBodyHamiltonian;
   for (int i = position; i < this->N-2; i++)
   {
    Identity = Identities[i+2];
    X.multiplyDirectProduct(Identity, Y);
    X = Y;
   }
   for (int i = position; i > 0; i--)
   {
    Identity = Identities[i-1];
    Identity.multiplyDirectProduct(X, Y);
    X = Y;
   }
   Z.add(X);
  }
  Z.exponentialize(-delta, X);
  TEMatrices[4] = X;
// 5. LastFirstTEMatricesMerged: X=exp(-(p5+p1)*deltaTime*H_{EVEN})
//    (T. Prosen & I. Pizorn, J. Phys. A: Math. Gen. 39 (2006), 5957-5964):
  p = 0.5;
  delta = p*deltaTime;
  Z.fillZeroes();
  for (int position = 0; position <= this->N-2; position += 2)
  {
   d0 = this->d[position]; d1 = this->d[position+1];
   TwoBodyHamiltonian = Matrix<T>(d0*d1, d0*d1);
   this->getTwoBodyHamiltonian(position, TwoBodyHamiltonian);
   X = TwoBodyHamiltonian;
   for (int i = position; i < this->N-2; i++)
   {
    Identity = Identities[i+2];
    X.multiplyDirectProduct(Identity, Y);
    X = Y;
   }
   for (int i = position; i > 0; i--)
   {
    Identity = Identities[i-1];
    Identity.multiplyDirectProduct(X, Y);
    X = Y;
   }
   Z.add(X);
  }
  Z.exponentialize(-delta, X);
  LastFirstTEMatricesMerged = X;
 }
}
