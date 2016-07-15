/// Singleton template class TimeEvolution implements time evolution for MPS.
/** The singleton template class TimeEvolution implements all functions for time evolution of MPS.
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

template<class T> class TimeEvolution
{
 public:

/// Constructs Singleton TimeEvolution<T> and returns reference.
  static TimeEvolution<T>& reference()
  {
   static TimeEvolution<T> TE;
   return TE;
  }

/// Evolves MPS in time.
/** This function evolves an initial MPS with a Hamiltonian in real or imaginary time. It uses a
    Trotter decomposition of the evolution operator of second or fourth order. The approximation has an
    error per time step and the program stops at timeStopAchieved if it can not fulfill this error.
    The evolved state is not normalized during evolution.
    \param MPS0 input: const MPS<T>&, the initial MPS
    \param H0 input: Hamiltonian<T>&, the Hamiltonian
    \param RealImaginaryTE input: const string&, must be "real" or "imaginary"
    \param TrotterOrder input: const string&, must be "2nd" or "4th"
    \param TrotterDecomposition input: const string&, must be "even-odd" or "oneMPO"
    \param timeStart input: double, the start time
    \param timeStep input: double, the time step
    \param timeStop input: double, the stop time
    \param eps input: double, the desired change in error in the MPS-MPO approximation
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed in the MPS-MPO approximation
    \param error input: double, the maximally allowed error for approximating a MPS-MPO product by a MPS
    \param errorAchieved output: double&, the largest errorAchieved in the MPS-MPO approximation
    \param numSweepsDone output: unsigned int&, the largest numSweepsDone in the MPS-MPO approximation
    \param timeStopAchieved output: double&, the final time possible with error
    \param MPS1 output: MPS<T>&, the final MPS */
  void timeEvolve(const MPS<T>& MPS0, Hamiltonian<T>& H0, const string& RealImaginaryTE,
                  const string& TrotterOrder, const string& TrotterDecomposition, double timeStart,
                  double timeStep, double timeStop, double eps, unsigned int maxNumSweeps, double error,
                  double& errorAchieved, unsigned int& numSweepsDone, double& timeStopAchieved, MPS<T>& MPS1) const;

/// Evolves MPS in time with TEBD.
/** This function evolves an initial MPS with a Hamiltonian in real or imaginary time with TEBD. It
    assumes a even-odd Trotter decomposition of the evolution operator. The approximation has an
    error per time step and the program stops at timeStopAchieved if it can not fulfill this error.
    \param MPS0 input: const MPS<T>&, the initial MPS
    \param H0 input: Hamiltonian<T>&, the Hamiltonian
    \param RealImaginaryTE input: const string&, must be "real" or "imaginary"
    \param TrotterOrder input: const string&, must be "2nd" or "4th"
    \param timeStart input: double, the start time
    \param timeStep input: double, the time step
    \param timeStop input: double, the stop time
    \param error input: double, the allowed error for approximating a MPS-MPO product by a MPS
    \param timeStopAchieved output: double&, the final time possible with error
    \param MPS1 output: MPS<T>&, the final MPS */
  void timeEvolveTEBD(const MPS<T>& MPS0, Hamiltonian<T>& H0, const string& RealImaginaryTE,
                      const string& TrotterOrder, double timeStart, double timeStep,
                      double timeStop, double error, double& timeStopAchieved, MPS<T>& MPS1) const;

/// Evolves vector in time exactly.
/** This function evolves an initial state vector Vector0 with a Hamiltonian in real or imaginary time.
    Initially, Vector0 will be normalized and the resulting vector Vector1 will also be normalized. For
    a time-independent Hamiltonian, the evolution operator is obtained from diagonalization. For a
    time-dependent Hamiltonian, a Runge-Kutta scheme of fourth order is used.
    \param Vector0 input: const vector<T>&, the initial vector, must be normalized
    \param H0 input: Hamiltonian<T>&, the Hamiltonian
    \param RealImaginaryTE input: const string&, must be "real" or "imaginary"
    \param timeStart input: double, the start time
    \param timeStep input: double, the time step
    \param timeStop input: double, the stop time
    \param Vector1 output: vector<T>&, the final vector, must have the same size as Vector0, will be
                           normalized */
  void timeEvolveExactly(const vector<T>& Vector0, Hamiltonian<T>& H0,
                         const string& RealImaginaryTE, double timeStart, double timeStep,
                         double timeStop, vector<T>& Vector1) const;

 private:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  TimeEvolution();

/// Standard copy constructor. Not implemented for this Singleton.
/** The standard copy constructor copies the input TimeEvolution into this. This constructor must not
    be implemented for this Singleton class.
    \param TimeEvolution0 input: const TimeEvolution<T>&, to be copied into this
    \sa TimeEvolution<T>& operator=(const TimeEvolution<T>& TimeEvolution0) */
  TimeEvolution(const TimeEvolution<T>& TimeEvolution0);

/// Assigns TimeEvolution to this. Not implemented for this Singleton.
/** The operator= allows to assign a TimeEvolution to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side TimeEvolution. This function must not be
    implemented for this Singleton class.
    \param TimeEvolution0 input: const TimeEvolution<T>&, to be copied into this
    \return TimeEvolution<T>&, a reference to the new this
    \sa TimeEvolution(const TimeEvolution<T>& TimeEvolution0) */
  TimeEvolution<T>& operator=(const TimeEvolution<T>& TimeEvolution0);
};

template<class T> void TimeEvolution<T>::timeEvolve(const MPS<T>& MPS0, Hamiltonian<T>& H0,
                                                    const string& RealImaginaryTE,
                                                    const string& TrotterOrder,
                                                    const string& TrotterDecomposition,
                                                    double timeStart, double timeStep, double timeStop,
                                                    double eps, unsigned int maxNumSweeps, double error,
                                                    double& errorAchieved, unsigned int& numSweepsDone, double& timeStopAchieved,
                                                    MPS<T>& MPS1) const
{
 vector< MPO<T> > TEMPOs;
 MPO<T> LastFirstTEMPOsMerged;
 double time = timeStart;
 double errorAchieved0, errorAchievedMax = 0.0; unsigned int numSweepsDone0, numSweepsDoneMax = 0;
 MPS1 = MPS0;
 MPS<T> MPS2(MPS1), MPS3(MPS1);
 errorAchieved = 1.0e300; numSweepsDone = 0; timeStopAchieved = time;
 if (H0.isTimeDependent())
 {
  while (time < timeStop)
  {
   H0.setTime(time);
   H0.getTEMPOs(RealImaginaryTE, TrotterOrder, TrotterDecomposition, timeStep, TEMPOs,
                LastFirstTEMPOsMerged);
   for (int j = 0; j < TEMPOs.size(); j++)
   {
    multiplyMPOMPS(TEMPOs[j], MPS1, eps, maxNumSweeps, errorAchieved0, numSweepsDone0, MPS2);
    if (errorAchieved0 > error)
    {
     MPS1 = MPS3;
     return;
    }
    if (errorAchieved0 > errorAchievedMax)
     errorAchievedMax = errorAchieved0;
    if (numSweepsDone0 > numSweepsDoneMax)
     numSweepsDoneMax = numSweepsDone0;
// in the case of real time evolution the evolution operator is unitary and the evolved state remains
// normalized within the Trotter error:
//    MPS2.normalize();
    MPS1 = MPS2;
   }
   errorAchieved = errorAchievedMax;
   numSweepsDone = numSweepsDoneMax;
   timeStopAchieved += timeStep;
   MPS3 = MPS1;
   time += timeStep;
  }
  H0.setTime(timeStart);
 }
 else if (!(H0.isTimeDependent()))
 {
  H0.getTEMPOs(RealImaginaryTE, TrotterOrder, TrotterDecomposition, timeStep, TEMPOs,
               LastFirstTEMPOsMerged);
  while (time < timeStop)
  {
   for (int j = 0; j < TEMPOs.size(); j++)
   {
    multiplyMPOMPS(TEMPOs[j], MPS1, eps, maxNumSweeps, errorAchieved0, numSweepsDone0, MPS2);
    if (errorAchieved0 > error)
    {
     MPS1 = MPS3;
     return;
    }
    if (errorAchieved0 > errorAchievedMax)
     errorAchievedMax = errorAchieved0;
    if (numSweepsDone0 > numSweepsDoneMax)
     numSweepsDoneMax = numSweepsDone0;
// in the case of real time evolution the evolution operator is unitary and the evolved state remains
// normalized within the Trotter error:
//    MPS2.normalize();
    MPS1 = MPS2;
   }
   errorAchieved = errorAchievedMax;
   numSweepsDone = numSweepsDoneMax;
   timeStopAchieved += timeStep;
   MPS3 = MPS1;
   time += timeStep;
  }
 }
}

template<class T> void TimeEvolution<T>::timeEvolveTEBD(const MPS<T>& MPS0, Hamiltonian<T>& H0,
                                                        const string& RealImaginaryTE,
                                                        const string& TrotterOrder,
                                                        double timeStart, double timeStep,
                                                        double timeStop, double error,
                                                        double& timeStopAchieved, MPS<T>& MPS1) const
{
 vector< MPO<T> > TEMPOs;
 MPO<T> LastFirstTEMPOsMerged;
 string TrotterDecomposition = "even-odd";
 H0.getTEMPOs(RealImaginaryTE, TrotterOrder, TrotterDecomposition, timeStep, TEMPOs,
              LastFirstTEMPOsMerged);
 double numTimeSteps = (timeStop - timeStart) / timeStep;
 double errorAchieved;
 MPS1 = MPS0;
 MPS1.normalize();
 MPS<T> MPS2(MPS1), MPS3(MPS1);
 timeStopAchieved = 0.0;
 for (int i = 1; i <= numTimeSteps; i++)
 {
  for (int j = 0; j < TEMPOs.size(); j++)
  {
   multiplyEvenOddMPOMPS(TEMPOs[j], MPS1, MPS2);
   errorAchieved = distanceMPOMPS(TEMPOs[j], MPS1, MPS2);
   if (errorAchieved > error)
   {
    MPS1 = MPS3;
    return;
   }
   MPS2.normalize();
   MPS1 = MPS2;
  }
  timeStopAchieved += timeStep;
  MPS3 = MPS1;
 }
}

template<class T> void TimeEvolution<T>::timeEvolveExactly(const vector<T>& Vector0, Hamiltonian<T>& H0,
                                                           const string& RealImaginaryTE,
                                                           double timeStart, double timeStep,
                                                           double timeStop,
                                                           vector<T>& Vector1) const
{
 unsigned int N = H0.getN();
 vector<unsigned int> d(N); H0.getd(d);
 unsigned int dim = 1;
 for (int i = 0; i < N; i++)
 {
  dim *= d[i];
 }
#ifdef DEBUG
 if ((Vector0.size() != dim) || (N == 0) ||
     ((RealImaginaryTE != "real") && (RealImaginaryTE != "imaginary")) || (Vector1.size() != dim))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void TimeEvolution<T>::" <<
          "timeEvolveExactly(const vector<T>& Vector0, const Hamiltonian<T>& H0, " <<
                            "const string& RealImaginaryTE, double timeStart, double timeStop, " <<
                            "vector<T>& Vector1) const: " <<
          "((Vector0.size() != dim) || (N == 0) || " <<
           "((RealImaginaryTE != real) && (RealImaginaryTE != imaginary)) || " <<
           "(Vector1.size() != dim))." << endl;
  exit(1);
 }
#endif
 double time = timeStart;
 if (!(H0.isTimeDependent()))
// if H0 is time-independent, we determine the exact time evolution operator by diagonalizing H0:
 {
  Matrix<T> MatrixH0(dim, dim);
  H0.getMatrix(MatrixH0);
  vector<double> W(dim); Matrix<T> Vr(dim, dim);
  MatrixH0.eigenDecompose(W, Vr);
  if (RealImaginaryTE == "real")
  {
   complex<double> prefactor = -complex<double>(0.0,1.0) * (timeStop-timeStart);
   for (int i = 0; i < dim; i++)
   {
    Vector1[i] = 0.0;
    for (int j = 0; j < dim; j++)
    {
     for (int k = 0; k < dim; k++)
     {
      Vector1[i] += exp(prefactor*W[k]) * Vector0[j] * Vr(i,k) *
                   MathAuxiliary::complexConjugate(Vr(j,k));
     }
    }
   }
  }
  else if (RealImaginaryTE == "imaginary")
  {
   double prefactor = -timeStep;
   vector<T> Vector2(Vector0);
   while (time < timeStop)
   {
    for (int i = 0; i < dim; i++)
    {
     Vector1[i] = 0.0;
     for (int j = 0; j < dim; j++)
     {
      for (int k = 0; k < dim; k++)
      {
       Vector1[i] += exp(prefactor*W[k]) * Vector2[j] * Vr(i,k) *
                     MathAuxiliary::complexConjugate(Vr(j,k));
      }
     }
    }
    normalize(Vector1);
    Vector2 = Vector1;
    time += timeStep;
   }
  }
 }
 else if (H0.isTimeDependent())
// if H0 is time-dependent, we do time evolution with a Runge-Kutta algorithm of fourth order:
 {
  double numTimeSteps = (timeStop - timeStart) / timeStep;
  double time = timeStart;
  Matrix< complex<double> > MatrixH00(dim, dim), MatrixH01(dim, dim), MatrixH02(dim, dim);
  H0.setTime(timeStart);
  H0.getMatrix(MatrixH02);
  complex<double>* Y = new complex<double>[dim];
  complex<double>* H00 = new complex<double>[dim*dim];
  complex<double>* H01 = new complex<double>[dim*dim];
  complex<double>* H02 = new complex<double>[dim*dim];
  int n = dim; double h = timeStep;
  complex<double>* Yout = new complex<double>[dim];
  for (int j = 0; j < dim; j++)
  {
   Yout[j] = Vector0[j];
  }
  for (int k = 0; k < dim; k++)
  {
   for (int j = 0; j < dim; j++)
   {
    H02[j+k*dim] = MatrixH02(j, k);
   }
  }
  while (time < timeStop)
  {
   for (int j = 0; j < dim; j++)
   {
    Y[j] = Yout[j];
   }
   for (int j = 0; j < dim*dim; j++)
   {
    H00[j] = H02[j];
   }
   time += 0.5 * timeStep;
   H0.setTime(time);
   H0.getMatrix(MatrixH01);
   for (int k = 0; k < dim; k++)
   {
    for (int j = 0; j < dim; j++)
    {
     H01[j+k*dim] = MatrixH01(j, k);
    }
   }
   time += 0.5 * timeStep;
   H0.setTime(time);
   H0.getMatrix(MatrixH02);
   for (int k = 0; k < dim; k++)
   {
    for (int j = 0; j < dim; j++)
    {
     H02[j+k*dim] = MatrixH02(j, k);
    }
   }
   zrk4_(Y, H00, H01, H02, &n, &h, Yout);
  }
  for (int i = 0; i < dim; i++)
  {
   Vector1[i] = Yout[i];
  }
  delete[] Y;
  delete[] H00; delete[] H01; delete[] H02;
  delete[] Yout;
 }
}

template<class T> TimeEvolution<T>::TimeEvolution() {}

template<class T> TimeEvolution<T>::TimeEvolution(const TimeEvolution<T>& TimeEvolution0) {}

template<class T> TimeEvolution<T>& TimeEvolution<T>::operator=(const TimeEvolution<T>& TimeEvolution0)
{}
