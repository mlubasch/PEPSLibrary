/// Abstract template class Hamiltonian2D implements interface for Hamiltonians in 2D.
/** The abstract template class Hamiltonian2D implements the interface for Hamiltonians in 2D.
    It is child class to the abstract template class Operator2D.
    \param BC string, the boundary conditions, is "open" or "periodic"
    \param N vector<unsigned int>(2), the number of sites, fulfills N.size()==2
    \param d Matrix<unsigned int>(N[0], N[1]), the physical dimensions, fulfills
             d.getDim0()==N[0] and d.getDim1()==N[1]
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

template<class T> class Hamiltonian2D : public Operator2D<T>
{

 typedef T (*PointerToFunction)(double time);

 public:

/// Sets boundary conditions.
/** This function sets the boundary conditions of this Hamiltonian2D.
    They must be either "open" or "periodic".
    \param BC0 input: const string&, the boundary conditions, must be "open" or "periodic" */
  void setBC(const string& BC0);

/// Returns boundary conditions.
/** The returned boundary conditions are either "open" or "periodic".
    \param BC0 output: string&, the boundary conditions of this Hamiltonian2D */
  void getBC(string& BC0) const { BC0 = this->BC; }

/// Sets number of rows and columns.
/** This function sets the number of rows and columns of this Hamiltonian2D, and both must be > 1.
    \param Nrows input: unsigned int, the number of rows, must be > 1
    \param Ncols input: unsigned int, the number of columns, must be > 1 */
  void setN(unsigned int Nrows, unsigned int Ncols);

/// Returns number of rows.
/** This function returns the number of rows of this Hamiltonian2D.
    \return unsigned int, the number of rows of this Hamiltonian2D */
  unsigned int getNrows() const { return this->N[0]; }

/// Returns number of columns.
/** This function returns the number of columns of this Hamiltonian2D.
    \return unsigned int, the number of columns of this Hamiltonian2D */
  unsigned int getNcols() const { return this->N[1]; }

/// Sets physical dimensions d.
/** This function sets the physical dimensions of this Hamiltonian2D.
    \param d0 input: const Matrix<unsigned int>&, the physical dimensions, must fulfill
                     d0.getDim0()==Nrows and d0.getDim1()==Ncols, and all entries must be > 1 */
  void setd(const Matrix<unsigned int>& d0);

/// Returns physical dimensions d.
/** This function returns the physical dimensions d of this Hamiltonian2D.
    \param d0 output: Matrix<unsigned int>&, the physical dimensions of this Hamiltonian2D */
  void getd(Matrix<unsigned int>& d0) const { d0 = this->d; }

/// Sets physical dimension d.
/** This function sets all the physical dimensions of this Hamiltonian2D to the same value d.
    \param d0 input: unsigned int, the physical dimension, must be > 1 */
  void setd(unsigned int d0);

/// Returns physical dimension d.
/** This function returns the physical dimension of this Hamiltonian2D at site (0, 0). It is useful
    if this Hamiltonian2D has all physical dimensions equal to each other.
    \return unsigned int, the physical dimension of this Hamiltonian2D at site (0, 0) */
  unsigned int getd() const { return this->d(0, 0); }

/// Sets parameters.
/** This function sets the parameters.
    \param Parameters0 input: const vector<T>&, the parameters, must fulfill
                              Parameters0.size() == this->Parameters.size() */
  void setParameters(const vector<T>& Parameters0);

/// Returns parameters.
/** This function returns the parameters of this Hamiltonian2D.
    \param Parameters0 output: vector<T>&, the parameters of this Hamiltonian2D */
  void getParameters(vector<T>& Parameters0) const { Parameters0 = this->Parameters; }

/// Sets Hamiltonian2D time-dependent.
/** This function sets this Hamiltonian2D time-dependent and passes the functions defining the
    time dependence.
    \param TimeFunctions0 input: const vector<PointerToFunction>&, the functions defining
                                 the time dependence of Parameters, must fulfill
                                 TimeFunctions0.size() == this->Parameters.size() */
  void setTimeDependent(const vector<PointerToFunction>& TimeFunctions0);

/// Sets Hamiltonian time-independent.
/** This function sets this Hamiltonian time-independent. */
  void unsetTimeDependent();

/// Returns time dependence.
/** This function returns true if this Hamiltonian2D is time-dependent, it returns false if
    this Hamiltonian2D is time-independent. */
  bool isTimeDependent() const { return this->timeDependent; }

/// Sets functions defining time dependence.
/** This function sets this Hamiltonian2D time-dependent and it sets the functions that define the
    time dependence of Parameters.
    \param TimeFunctions0 input: const vector<PointerToFunction>&, the time functions */
  void setTimeFunctions(const vector<PointerToFunction>& TimeFunctions0);

/// Returns functions defining time dependence.
/** This function returns the functions that define the time dependence of Parameters.
    \param TimeFunctions0 output: vector<PointerToFunction>&, the time functions of
                                  this Hamiltonian2D */
  void getTimeFunctions(vector<PointerToFunction>& TimeFunctions0) const
  { TimeFunctions0 = this->TimeFunctions; }

/// Sets time.
/** This function sets the time and updates Parameters according to their time dependence defined
    in TimeFunctions.
    \param time0 input: double, the time */
  void setTime(double time0);

/// Returns time.
/** This function returns the time of this Hamiltonian2D.
    \return double, the time of this Hamiltonian2D */
  double getTime() const { return this->time; }

/// Returns two-site Hamiltonian MPOs.
/** This function returns the two-site Hamiltonian MPOs corresponding to the two-site Trotter gates.
    \param HMPOs output: vector< MPO<T> >&, the two-site Hamiltonian MPOs */
  virtual void getHMPOs(vector< MPO<T> >& HMPOs) const = 0;

/// Returns two-site Trotter gates.
/** This function returns the two-site Trotter gates.
    \param RealImaginaryTE input: const string&, must be "real" or "imaginary"
    \param timeStep input: double, the time step
    \param TEMPOs output: vector< MPO<T> >&, the two-site Trotter gates */
  virtual void getTEMPOs(const string& RealImaginaryTE, double timeStep, vector< MPO<T> >& TEMPOs) = 0;

/// Returns PEPOs for time evolution.
/** This function returns the PEPOs resulting from a Trotter decomposition of the evolution operator.
    \param RealImaginaryTE input: const string&, must be "real" or "imaginary"
    \param timeStep input: double, the time step
    \param TEPEPOs output: vector< PEPO<T> >&, the time evolution PEPOs */
  virtual void getTEPEPOs(const string& RealImaginaryTE, double timeStep, vector< PEPO<T> >& TEPEPOs) = 0;

 protected:

/// Boundary conditions BC.
/** The boundary conditions are either "open" or "periodic". */
  string BC;

/// Number of sites N.
/** N[0]=Nrows is the number of rows and N[1]=Ncols is the number of columns of this Hamiltonian2D. */
  vector<unsigned int> N;

/// Physical dimensions d.
/** The Matrix d has dimensions d.getDim0()==N[0] and d.getDim1()==N[1]. */
  Matrix<unsigned int> d;

/// Parameters.
/** The Parameters of this Hamiltonian2D. */
  vector<T> Parameters;

/// Time dependence.
/** timeDependent==false for a time-independent Hamiltonian2D, and timeDependent==true for a
    time-dependent Hamiltonian2D. */
  bool timeDependent;

/// Time dependence of Parameters.
/** The functions TimeFunctions define the time dependence of Parameters when this Hamiltonian2D is
    time-dependent. */
  vector<PointerToFunction> TimeFunctions;

/// Time.
/** The time of this Hamiltonian2D. */
  double time;

};

template<class T> void Hamiltonian2D<T>::setBC(const string& BC0)
{
#ifdef DEBUG
 if ((BC0 != "open") && (BC0 != "periodic"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Hamiltonian2D<T>::" <<
          "setBC(const string& BC0): " <<
          "((BC0 != open) && (BC0 != periodic))." << endl;
  exit(1);
 }
#endif
 this->BC = BC0;
}

template<class T> void Hamiltonian2D<T>::setN(unsigned int Nrows, unsigned int Ncols)
{
#ifdef DEBUG
 if ((Nrows < 2) || (Ncols < 2))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Hamiltonian2D<T>::" <<
          "setN(unsigned int Nrows, unsigned int Ncols): " <<
          "((Nrows < 2) || (Ncols < 2))." << endl;
  exit(1);
 }
#endif
 this->N[0] = Nrows; this->N[1] = Ncols;
}

template<class T> void Hamiltonian2D<T>::setd(const Matrix<unsigned int>& d0)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (d0.getDim0() != this->N[0]) || (d0.getDim1() != this->N[1]))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Hamiltonian2D<T>::" <<
          "setd(const Matrix<unsigned int>& d0): " <<
      "((this->N[0] == 0) || (d0.getDim0() != this->N[0]) || (d0.getDim1() != this->N[1]))." << endl;
  exit(1);
 }
 for (int j = 0; j < d0.getDim1(); j++)
 {
  for (int i = 0; i < d0.getDim0(); i++)
  {
   if (d0(i, j) < 2)
   {
    cerr << "Program terminated because of error in function " <<
            "template<class T> void Hamiltonian2D<T>::" <<
            "setd(const Matrix<unsigned int>& d0): " <<
            "(d0(" << i << ", " << j << ") < 2)." << endl;
    exit(1);
   }
  }
 }
#endif
 this->d = d0;
}

template<class T> void Hamiltonian2D<T>::setd(unsigned int d0)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (d0 < 2))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Hamiltonian2D<T>::" <<
          "setd(unsigned int d0): " <<
          "((this->N[0] == 0) || (d0 < 2))." << endl;
  exit(1);
 }
#endif
 this->d = Matrix<unsigned int>(this->N[0], this->N[1]);
 for (int j = 0; j < this->N[1]; j++)
 {
  for (int i = 0; i < this->N[0]; i++)
  {
   this->d(i, j) = d0;
  }
 }
}

template<class T> void Hamiltonian2D<T>::setParameters(const vector<T>& Parameters0)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (Parameters0.size() != this->Parameters.size()))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Hamiltonian2D<T>::" <<
          "setParameters(const vector<T>& Parameters0): " <<
          "((this->N[0] == 0) || (Parameters0.size() != this->Parameters.size()))." << endl;
  exit(1);
 }
#endif
 this->Parameters = Parameters0;
}

template<class T> void Hamiltonian2D<T>::setTimeDependent(const vector<PointerToFunction>&
                                                          TimeFunctions0)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (TimeFunctions0.size() != this->Parameters.size()))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Hamiltonian2D<T>::" <<
          "setTimeDependent(const vector<PointerToFunction>& TimeFunctions0): " <<
          "((this->N[0] == 0) || (TimeFunctions0.size() != this->Parameters.size()))." << endl;
  exit(1);
 }
#endif
 this->timeDependent = true;
 this->TimeFunctions = TimeFunctions0;
}

template<class T> void Hamiltonian2D<T>::unsetTimeDependent()
{
#ifdef DEBUG
 if (this->N[0] == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Hamiltonian2D<T>::" <<
          "unsetTimeDependent(): " <<
          "(this->N[0] == 0)." << endl;
  exit(1);
 }
#endif
 this->timeDependent = false;
 for (int i = 0; i < this->TimeFunctions.size(); i++)
  this->TimeFunctions[i] = 0;
}

template<class T> void Hamiltonian2D<T>::setTimeFunctions(const vector<PointerToFunction>&
                                                          TimeFunctions0)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (TimeFunctions0.size() != this->Parameters.size()))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Hamiltonian2D<T>::" <<
          "setTimeFunctions(const vector<PointerToFunction>& TimeFunctions0): " <<
          "((this->N[0] == 0) || (TimeFunctions0.size() != this->Parameters.size()))." << endl;
  exit(1);
 }
#endif
 this->timeDependent = true;
 this->TimeFunctions = TimeFunctions0;
}

template<class T> void Hamiltonian2D<T>::setTime(double time0)
{
#ifdef DEBUG
 if ((this->N[0] == 0) || (this->timeDependent == false))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Hamiltonian2D<T>::" <<
          "setTime(double time0): " <<
          "((this->N[0] == 0) || (this->timeDependent == false))." << endl;
  exit(1);
 }
#endif
 this->time = time0;
 for (int i = 0; i < this->Parameters.size(); i++)
  this->Parameters[i] = this->TimeFunctions[i](time0);
}
