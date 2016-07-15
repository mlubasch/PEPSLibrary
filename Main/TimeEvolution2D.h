/// Singleton template class TimeEvolution2D implements time evolution for PEPS.
/** The singleton template class TimeEvolution2D implements all functions for time evolution with PEPS.
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

template<class T> class TimeEvolution2D
{
 public:

/// Constructs Singleton TimeEvolution2D<T> and returns reference.
  static TimeEvolution2D<T>& reference()
  {
   static TimeEvolution2D<T> TE;
   return TE;
  }

/// Evolves PEPS in time.
/** This function evolves an initial PEPS0 with a Hamiltonian H0 in real or imaginary time.
    It uses a first order Trotter decomposition of the evolution operator written as PEPOs and subsequently approximates the
    product of these time evolution PEPOs with the evolved PEPS.
    For each time evolution PEPO, the tensors are updated sequentially with convergence precision eps within maxNumSweeps sweeps.
    Environment boundary-MPSs are computed with D2Env, epsEnv, and maxNumSweepsEnv.
    The local tensor update can be reduced from the full tensor to the smallest subtensor including the physical bond.
    Either the pseudoinverse or the closest positive norm-matrix is used in the local update algorithm.
    PEPS0 can be a Concatenated PEPS of concatenation level M.
    The evolved state is normalized during evolution.
    \param H0 input: Hamiltonian2D<T>&, the Hamiltonian
    \param RealImaginaryTE input: const string&, must be "real" or "imaginary"
    \param timeStart input: double, the start time
    \param timeStep input: double, the time-step
    \param timeStop input: double, the stop time
    \param eps input: double, the convergence precision per time evolution PEPO
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed per time evolution PEPO
    \param D2Env input: unsigned int, the maximal virtual bond dimension in the environment approximation
    \param epsEnv input: double, the convergence precision in the environment approximation
    \param maxNumSweepsEnv input: unsigned int, the maximal number of sweeps allowed in the environment approximation
    \param UpdateTensor input: const string&, the tensor in the local tensor update, must be "full" or "reduced"
    \param UpdateMode input: const string&, the local update algorithm, must be "Pseudoinverse" for pseudoinverse or
                             "PositiveNormMatrix" for the closest positive norm-matrix
    \param cutoff input: double, the cutoff, singular values or eigenvalues x_{i} with (x_{i}/x_{max} <= cutoff) are put to zero
    \param PEPS0 input/output: PEPS<T>&, on input the initial PEPS, on output the final PEPS
    \param M optional input: unsigned int, if given the concatenation level
    \sa template<class T> void multiplyPEPOPEPS(const PEPO<T>& PEPO0, const PEPS<T>& PEPS0,
                                                double eps, unsigned int maxNumSweeps,
                                                unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv,
                                                const string& UpdateTensor, const string& UpdateMode, double cutoff,
                                                double& epsAchieved, unsigned int& numSweepsDone,
                                                PEPS<T>& PEPS1)
    \sa template<class T> T PEPS<T>::normalize(const string& Direction, unsigned int D2, double eps, unsigned int maxNumSweeps) */
  void timeEvolve(Hamiltonian2D<T>& H0, const string& RealImaginaryTE, double timeStart, double timeStep, double timeStop,
                  double eps, unsigned int maxNumSweeps,
                  unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv,
                  const string& UpdateTensor, const string& UpdateMode, double cutoff,
                  PEPS<T>& PEPS0, unsigned int M = 1) const;

/// Evolves PEPS in time.
/** This function evolves an initial PEPS0 with a Hamiltonian H0 in real or imaginary time, updating one tensor pair after the other.
    For each set of Trotter gates, the tensors are locally updated either independently or sequentially within numSweeps sweeps over the lattice.
    The bTensor can be constructed without or with all Trotter gates in the environment.
    The local tensor update can be reduced from the full tensor to the smallest subtensor including the physical bond.
    Either the pseudoinverse, the closest positive norm-matrix, or the gauge is used in the local update algorithm.
    The evolved state is normalized during evolution using D2=1.
    \param H0 input: Hamiltonian2D<T>&, the Hamiltonian
    \param RealImaginaryTE input: const string&, must be "real" or "imaginary"
    \param timeStart input: double, the start time
    \param timeStep input: double, the time-step
    \param timeStop input: double, the stop time
    \param SweepUpdate: const string&, must be "independent" or "sequential"
    \param numSweeps input: unsigned int, the number of update sweeps per set of Trotter gates
    \param bEnvironment: const string&, without or with the Trotter gates, must be "simplified" or "full"
    \param D2Env input: unsigned int, the maximal virtual bond dimension in the environment approximation
    \param epsEnv input: double, the convergence precision in the environment approximation
    \param maxNumSweepsEnv input: unsigned int, the maximal number of sweeps allowed in the environment approximation
    \param UpdateTensor input: const string&, the tensor in the local tensor update, must be "reduced" or "full"
    \param UpdateMode input: const string&, the local update algorithm, must be "Pseudoinverse" for pseudoinverse,
                             "PositiveNormMatrix" for the closest positive norm-matrix, or "Gauge"
    \param epsLoc input: double, the convergence precision in the local update
    \param maxNumSweepsLoc input: unsigned int, the maximal number of sweeps allowed in the local update
    \param cutoff input: double, the cutoff, singular values or eigenvalues x_{i} with (x_{i}/x_{max} <= cutoff) are put to zero
    \param PEPS0 input/output: PEPS<T>&, on input the initial PEPS, on output the final PEPS
    \sa void TimeEvolution2D<T>::updateTensor(unsigned int positionRow, unsigned int positionCol, const string& UpdateDirection,
                                              MPO<T>& NormMPO, const string& bEnvironment, const MPO<T>& bMPO, const MPO<T>& TEMPO,
                                              const string& UpdateTensor, const string& UpdateMode, double eps, unsigned int maxNumSweeps, double cutoff,
                                              PEPS<T>& PEPS0) const */
  void timeEvolve(Hamiltonian2D<T>& H0, const string& RealImaginaryTE, double timeStart, double timeStep, double timeStop,
                  const string& SweepUpdate, unsigned int numSweeps,
                  const string& bEnvironment, unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv,
                  const string& UpdateTensor, const string& UpdateMode, double epsLoc, unsigned int maxNumSweepsLoc, double cutoff,
                  PEPS<T>& PEPS0) const;

/// Evolves PEPS in time with Cluster-Update.
/** This function evolves an initial PEPS0 with a Hamiltonian H0 in real or imaginary time, updating one tensor pair after the other.
    For each set of Trotter gates, the tensors are locally updated either independently or sequentially within numSweeps sweeps over the lattice.
    The bTensor can be constructed without or with all Trotter gates in the environment.
    D2sEnv specifies the bond dimensions for the environment approximation inside the cluster sorted from closest column or row to farthest,
    D2sEnv.size() determines the cluster-size, and outside the cluster D2Env=1.
    The local tensor update can be reduced from the full tensor to the smallest subtensor including the physical bond.
    Either the pseudoinverse, the closest positive norm-matrix, or the gauge is used in the local update algorithm.
    The evolved state is normalized during evolution using D2=1.
    \param H0 input: Hamiltonian2D<T>&, the Hamiltonian
    \param RealImaginaryTE input: const string&, must be "real" or "imaginary"
    \param timeStart input: double, the start time
    \param timeStep input: double, the time-step
    \param timeStop input: double, the stop time
    \param SweepUpdate: const string&, must be "independent" or "sequential"
    \param numSweeps input: unsigned int, the number of update sweeps per set of Trotter gates
    \param bEnvironment: const string&, without or with the Trotter gates, must be "simplified" or "full"
    \param D2sEnv input: const vector<unsigned int>&, the maximal virtual bond dimensions for the environment approximation in the cluster,
                         sorted from the closest column or row to the farthest
    \param epsEnv input: double, the convergence precision in the environment approximation
    \param maxNumSweepsEnv input: unsigned int, the maximal number of sweeps allowed in the environment approximation
    \param UpdateTensor input: const string&, the tensor in the local tensor update, must be "reduced" or "full"
    \param UpdateMode input: const string&, the local update algorithm, must be "Pseudoinverse" for pseudoinverse,
                             "PositiveNormMatrix" for the closest positive norm-matrix, or "Gauge"
    \param epsLoc input: double, the convergence precision in the local update
    \param maxNumSweepsLoc input: unsigned int, the maximal number of sweeps allowed in the local update
    \param cutoff input: double, the cutoff, singular values or eigenvalues x_{i} with (x_{i}/x_{max} <= cutoff) are put to zero
    \param PEPS0 input/output: PEPS<T>&, on input the initial PEPS, on output the final PEPS
    \sa void TimeEvolution2D<T>::updateTensor(unsigned int positionRow, unsigned int positionCol, const string& UpdateDirection,
                                              MPO<T>& NormMPO, const string& bEnvironment, const MPO<T>& bMPO, const MPO<T>& TEMPO,
                                              const string& UpdateTensor, const string& UpdateMode, double eps, unsigned int maxNumSweeps, double cutoff,
                                              PEPS<T>& PEPS0) const */
  void timeEvolve(Hamiltonian2D<T>& H0, const string& RealImaginaryTE, double timeStart, double timeStep, double timeStop,
                  const string& SweepUpdate, unsigned int numSweeps,
                  const string& bEnvironment, const vector<unsigned int>& D2sEnv, double epsEnv, unsigned int maxNumSweepsEnv,
                  const string& UpdateTensor, const string& UpdateMode, double epsLoc, unsigned int maxNumSweepsLoc, double cutoff,
                  PEPS<T>& PEPS0) const;

/// Evolves PEPS in time with TEBD.
/** This function evolves an initial PEPS0 with a Hamiltonian H0 in real or imaginary time with TEBD, updating one
    tensor pair after the other with
       void TimeEvolution2D<T>::updateTensorTEBD(const string& UpdateTensor,
                                                 unsigned int positionRow, unsigned int positionCol,
                                                 const string& Direction, const MPO<T>& TEMPO,
                                                 PEPS<T>& PEPS0, Matrix< Matrix<T> >& LambdasV,
                                                 Matrix< Matrix<T> >& LambdasH) const   .
    PEPS0 must be in canonical form with lambda matrices LambdasV on vertical bonds and lambda matrices LambdasH on
    horizontal bonds. If both PEPS0 and H0 are uniform, only part of the tensors is processed. The local tensor
    update can be reduced from the full tensor to the smallest subtensor including the physical bond.
    \param Type input: const string&, the type of PEPS0 and H0, must be "uniform" or "general"
    \param UpdateTensor input: const string&, the tensor in the local tensor update, must be "reduced" or "full"
    \param H0 input: Hamiltonian2D<T>&, the Hamiltonian
    \param RealImaginaryTE input: const string&, must be "real" or "imaginary"
    \param timeStart input: double, the start time
    \param timeStep input: double, the time step
    \param timeStop input: double, the stop time
    \param PEPS0 input/output: PEPS<T>&, on input the initial PEPS, on output the final PEPS
    \param LambdasV input/output: Matrix< Matrix<T> >&, lambda matrices on the vertical bonds
    \param LambdasH input/output: Matrix< Matrix<T> >&, lambda matrices on the horizontal bonds
    \sa void TimeEvolution2D<T>::updateTensorTEBD(const string& UpdateTensor,
                                                  unsigned int positionRow, unsigned int positionCol,
                                                  const string& Direction, const MPO<T>& TEMPO,
                                                  PEPS<T>& PEPS0, Matrix< Matrix<T> >& LambdasV,
                                                  Matrix< Matrix<T> >& LambdasH) const */
  void timeEvolveTEBD(const string& Type, const string& UpdateTensor,
                      Hamiltonian2D<T>& H0, const string& RealImaginaryTE,
                      double timeStart, double timeStep, double timeStop,
                      PEPS<T>& PEPS0, Matrix< Matrix<T> >& LambdasV, Matrix< Matrix<T> >& LambdasH) const;

/// Evolves vector in time exactly.
/** This function evolves an initial state vector Vector0 with a Hamiltonian2D H0 in real or imaginary time.
    Initially, if (RealImaginaryTE == imaginary) then Vector0 will be normalized during evolution, and then,
    the resulting vector Vector1 will also be normalized. For a time-independent Hamiltonian, the evolution
    operator is obtained from diagonalization. For a time-dependent Hamiltonian, a Runge-Kutta scheme of
    fourth order is used.
    \param Vector0 input: const vector<T>&, the initial vector,
                          if (RealImaginaryTE == imaginary) then Vector0 will be normalized during evolution
    \param H0 input: Hamiltonian2D<T>&, the Hamiltonian
    \param RealImaginaryTE input: const string&, must be "real" or "imaginary"
    \param timeStart input: double, the start time
    \param timeStep input: double, the time step
    \param timeStop input: double, the stop time
    \param Vector1 output: vector<T>&, the final vector, must have the same size as Vector0,
                           if (RealImaginaryTE == imaginary) then Vector1 will be normalized */
  void timeEvolveExactly(const vector<T>& Vector0, Hamiltonian2D<T>& H0,
                         const string& RealImaginaryTE, double timeStart, double timeStep,
                         double timeStop, vector<T>& Vector1) const;

/// Computes energy expectation value.
/** This function computes the energy expectation value of PEPS0 for a Hamiltonian H0, and returns that value divided by the norm.
    It successively contracts the individual Hamiltonian MPOs, which each comprise all Hamiltonian terms acting on the same two neighbouring sites,
    and makes use of previously obtained environment.
    \param PEPS0 input: const PEPS<T>&, the PEPS
    \param H0 input: const Hamiltonian2D<T>&, the Hamiltonian
    \param D2 input: unsigned int, the maximal virtual bond dimension in the environment approximation
    \param eps input: double, the convergence precision in the environment approximation
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed in the environment approximation
    \return T, the resulting energy expectation value divided by the norm */
  T expectationValue(const PEPS<T>& PEPS0, const Hamiltonian2D<T>& H0, unsigned int D2, double eps, unsigned int maxNumSweeps) const;

/// Computes cluster energy expectation value.
/** This function computes the cluster energy expectation value of PEPS0 for a Hamiltonian H0, and returns that value divided by the cluster norm.
    It successively contracts the individual Hamiltonian MPOs, which each comprise all Hamiltonian terms acting on the same two neighbouring sites,
    and makes use of previously obtained environment.
    D2s specifies the bond dimensions for the environment approximation inside the cluster sorted from closest column or row to farthest,
    D2s.size() determines the cluster-size, and outside the cluster D2=1.
    \param PEPS0 input: const PEPS<T>&, the PEPS
    \param H0 input: const Hamiltonian2D<T>&, the Hamiltonian
    \param D2s input: const vector<unsigned int>&, the maximal virtual bond dimensions for the environment approximation in the cluster,
                      sorted from the closest column or row to the farthest
    \param eps input: double, the convergence precision in the environment approximation
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed in the environment approximation
    \return T, the resulting cluster energy expectation value divided by the cluster norm */
  T expectationValue(const PEPS<T>& PEPS0, const Hamiltonian2D<T>& H0, const vector<unsigned int>& D2s, double eps, unsigned int maxNumSweeps) const;

 private:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  TimeEvolution2D();

/// Standard copy constructor. Not implemented for this Singleton.
/** The standard copy constructor copies the input TimeEvolution2D into this. This constructor must not
    be implemented for this Singleton class.
    \param TimeEvolution0 input: const TimeEvolution2D<T>&, to be copied into this
    \sa TimeEvolution2D<T>& operator=(const TimeEvolution2D<T>& TimeEvolution0) */
  TimeEvolution2D(const TimeEvolution2D<T>& TimeEvolution0);

/// Assigns TimeEvolution2D to this. Not implemented for this Singleton.
/** The operator= allows to assign a TimeEvolution2D to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side TimeEvolution2D. This function must not be
    implemented for this Singleton class.
    \param TimeEvolution0 input: const TimeEvolution2D<T>&, to be copied into this
    \return TimeEvolution2D<T>&, a reference to the new this
    \sa TimeEvolution2D(const TimeEvolution2D<T>& TimeEvolution0) */
  TimeEvolution2D<T>& operator=(const TimeEvolution2D<T>& TimeEvolution0);

/// Returns initial norm-MPSs.
/** This function computes the initial boundary-MPSs for the norm of PEPS0, before a set of Trotter gates is applied.
    If UpdateDirection==vertical, then NormMPSs will have size PEPS0.N[1] and comprise all norm-MPSs for columns colPosition > 0.
    The norm-MPS for column colPosition will be NormMPSs[colPosition] and its physical indices will be pointing left.
    If UpdateDirection==horizontal, then NormMPSs will have size PEPS0.N[0] and comprise all norm-MPSs for rows rowPosition > 0.
    The norm-MPS for row rowPosition will be NormMPSs[rowPosition] and its physical indices will be pointing up.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices. */
  void getInitialNormMPSs(const PEPS<T>& PEPS0, const string& UpdateDirection,
                          unsigned int D2, double eps, unsigned int maxNumSweeps,
                          vector< MPS<T> >& NormMPSs) const;

/// Updates norm-MPSs on boundary.
/** This function updates the boundary-MPSs for the norm of PEPS0, after the update of a boundary column or row.
    If WhichBoundary==left, then the boundary-MPS for the left boundary will be written into NormMPSs[0] with physical indices pointing right.
    If WhichBoundary==top, then the boundary-MPS for the top boundary will be written into NormMPSs[0] with physical indices pointing down.
    If WhichBoundary==right, then the boundary-MPS for the right boundary will be written into NormMPSs[Ncols-1] with physical indices pointing left.
    If WhichBoundary==bottom, then the boundary-MPS for the bottom boundary will be written into NormMPSs[Nrows-1] with physical indices pointing up.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices. */
  void updateNormMPSsBoundary(const PEPS<T>& PEPS0, const string& WhichBoundary,
                              unsigned int D2, double eps, unsigned int maxNumSweeps,
                              vector< MPS<T> >& NormMPSs) const;

/// Updates norm-MPSs in bulk.
/** This function updates the boundary-MPSs for the norm of PEPS0, after the update of a bulk column or row.
    If SweepDirection==right, then the boundary-MPS for column updatePosition will be written into NormMPSs[updatePosition] with physical indices pointing right.
    If SweepDirection==down, then the boundary-MPS for row updatePosition will be written into NormMPSs[updatePosition] with physical indices pointing down.
    If SweepDirection==left, then the boundary-MPS for column updatePosition will be written into NormMPSs[updatePosition] with physical indices pointing left.
    If SweepDirection==up, then the boundary-MPS for row updatePosition will be written into NormMPSs[updatePosition] with physical indices pointing up.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices. */
  void updateNormMPSsBulk(const PEPS<T>& PEPS0, const string& SweepDirection, unsigned int updatePosition,
                          unsigned int D2, double eps, unsigned int maxNumSweeps,
                          vector< MPS<T> >& NormMPSs) const;

/// Returns initial boundary norm-tensors.
/** This function computes the initial norm-tensors of a boundary of the norm of PEPS0, by contracting the boundary with
    the neighbouring norm-MPS. The resulting tensors are the initial norm-tensors for the update of that boundary.
    If WhichBoundary==left or WhichBoundary==right, then the norm-tensors lie above and below the tensor-pair, and they are enumerated like
    the rows in PEPS0.
    If WhichBoundary==top or WhichBoundary==bottom, then the norm-tensors lie left and right the tensor-pair, and they are enumerated like
    the columns in PEPS0.
    Each norm-tensor has the shape required by the norm-MPO. */
  void getInitialBoundaryNormTensors(const PEPS<T>& PEPS0, const string& WhichBoundary,
                                     const MPS<T>& NormMPS,
                                     vector< Tensor<T> >& NormTensors) const;

/// Returns initial bulk norm-tensors.
/** This function computes the initial norm-tensors of the bulk of the norm of PEPS0, by contracting the bulk column or row with
    its neighbouring norm-MPSs. The resulting tensors are the initial norm-tensors for the update of that bulk column or row.
    If UpdateDirection==vertical, then NormMPS1 denotes the left and NormMPS2 the right NormMPS to column updatePosition, and the resulting
    norm-tensors lie above and below the tensor-pair, enumerated like the rows in PEPS0.
    If UpdateDirection==horizontal, then NormMPS1 denotes the lower and NormMPS2 the upper NormMPS to row updatePosition, and the resulting
    norm-tensors lie left and right the tensor-pair, enumerated like the columns in PEPS0.
    Each norm-tensor has the shape required by the norm-MPO. */
  void getInitialBulkNormTensors(const PEPS<T>& PEPS0, const string& UpdateDirection, unsigned int updatePosition,
                                 const MPS<T>& NormMPS1, const MPS<T>& NormMPS2,
                                 vector< Tensor<T> >& NormTensors) const;

/// Updates boundary norm-tensors.
/** This function updates the norm-tensors for the norm of PEPS0, after the update of a tensor-pair on the boundary.
    If WhichBoundary==left, then NormMPS is the boundary-MPS on column 1 with physical indices pointing left, and tensorPosition denotes the
    row position of the first tensor of the the tensor-pair.
    If WhichBoundary==top, then NormMPS is the boundary-MPS on row 1 with physical indices pointing up, and tensorPosition denotes the
    column position of the first tensor of the the tensor-pair.
    If WhichBoundary==right, then NormMPS is the boundary-MPS on column Ncols-2 with physical indices pointing right, and tensorPosition denotes the
    row position of the first tensor of the the tensor-pair.
    If WhichBoundary==bottom, then NormMPS is the boundary-MPS on row Nrows-2 with physical indices pointing down, and tensorPosition denotes the
    column position of the first tensor of the the tensor-pair.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices.
    Each norm-tensor has the shape required by the norm-MPO. */
  void updateBoundaryNormTensors(const PEPS<T>& PEPS0, const string& WhichBoundary,
                                 const MPS<T>& NormMPS, unsigned int tensorPosition,
                                 vector< Tensor<T> >& NormTensors) const;

/// Updates bulk norm-tensors.
/** This function updates the norm-tensors for the norm of PEPS0, after the update of a tensor-pair in the bulk.
    If UpdateDirection==vertical, then NormMPS1 is the left and NormMPS2 the right boundary-MPS to column updatePosition, and tensorPosition
    denotes the row position of the first tensor of the tensor-pair.
    If UpdateDirection==horizontal, then NormMPS1 is the lower and NormMPS2 the upper boundary-MPS to row updatePosition, and tensorPosition
    denotes the column position of the first tensor of the tensor-pair.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices.
    Each norm-tensor has the shape required by the norm-MPO. */
  void updateBulkNormTensors(const PEPS<T>& PEPS0, const string& UpdateDirection, unsigned int updatePosition,
                             const MPS<T>& NormMPS1, const MPS<T>& NormMPS2,
                             unsigned int tensorPosition, vector< Tensor<T> >& NormTensors) const;

/// Returns boundary norm-MPO.
/** This function computes the norm-MPO for the update of a tensor-pair on the boundary. The norm-MPO is the periodic MPO consisting of the 6 tensors
    surrounding the tensor-pair in the norm contraction.
    If WhichBoundary==left, then NormMPS denotes the boundary-MPS on column 1 with physical indices pointing left, and tensorPosition denotes
    the row position of the first tensor of the tensor-pair.
    If WhichBoundary==top, then NormMPS denotes the boundary-MPS on row 1 with physical indices pointing up, and tensorPosition denotes
    the column position of the first tensor of the tensor-pair.
    If WhichBoundary==right, then NormMPS denotes the boundary-MPS on column Ncols-2 with physical indices pointing right, and tensorPosition denotes
    the row position of the first tensor of the tensor-pair.
    If WhichBoundary==bottom, then NormMPS denotes the boundary-MPS on row Nrows-2 with physical indices pointing down, and tensorPosition denotes
    the column position of the first tensor of the tensor-pair.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices.
    NormTensors stores the norm-tensors with required shape for NormMPO.
    NormMPO is enumerated clockwise around the tensor-pair, starting from the top norm-tensor if the tensor-pair is oriented vertically, or starting from
    the left norm-tensor if the tensor-pair is oriented horizontally. */
  void getBoundaryNormMPO(const string& WhichBoundary, const MPS<T>& NormMPS, const vector< Tensor<T> >& NormTensors,
                          unsigned int tensorPosition, MPO<T>& NormMPO) const;

/// Returns bulk norm-MPO.
/** This function computes the norm-MPO for the update of a tensor-pair in the bulk. The norm-MPO is the periodic MPO consisting of the 6 tensors
    surrounding the tensor-pair in the norm contraction.
    NormMPS1 denotes the left or lower norm-MPS and NormMPS2 the right or upper norm-MPS to the updated column or row.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices.
    NormTensors stores the norm-tensors with required shape for NormMPO.
    TensorPosition denotes the row or column position of the first tensor of the tensor-pair.
    NormMPO is enumerated clockwise around the tensor-pair, starting from the top norm-tensor if the tensor-pair is oriented vertically, or starting from
    the left norm-tensor if the tensor-pair is oriented horizontally. */
  void getBulkNormMPO(const MPS<T>& NormMPS1, const MPS<T>& NormMPS2, const vector< Tensor<T> >& NormTensors,
                      unsigned int tensorPosition, MPO<T>& NormMPO) const;

/// MPS0 is initialized as a random positive separable boundary-MPS.
/** This function assumes MPS0 to be a boundary-MPS, i.e. MPS0.d==D*D.
    It initializes MPS0 as a random separable MPS and makes it positive by multiplying each tensor of the corresponding boundary-MPO with its adjoint. */
  static void setRandomPositiveSeparableBoundaryMPS(MPS<T>& MPS0, const vector<unsigned int>& Seed);

/// Returns initial positive separable norm-MPSs for Cluster-Update.
/** This function computes the initial positive separable boundary-MPSs for the norm of PEPS0, which are located outside the cluster of size clusterSize.
    These SeparableNormMPSs are computed initially before a set of Trotter gates is applied to PEPS0.
    If UpdateDirection==vertical, then SeparableNormMPSs will have size PEPS0.N[1] and comprise all separable norm-MPSs for columns colPosition > clusterSize.
    The separable norm-MPS for column colPosition will be SeparableNormMPSs[colPosition] and its physical indices will be pointing left.
    If UpdateDirection==horizontal, then SeparableNormMPSs will have size PEPS0.N[0] and comprise all separable norm-MPSs for rows rowPosition > clusterSize.
    The separable norm-MPS for row rowPosition will be SeparableNormMPSs[rowPosition] and its physical indices will be pointing up.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices. */
  void getInitialSeparableNormMPSs(const PEPS<T>& PEPS0, const string& UpdateDirection,
                                   unsigned int clusterSize, double eps, unsigned int maxNumSweeps,
                                   vector< MPS<T> >& SeparableNormMPSs) const;

/// Updates positive separable norm-MPSs on boundary for Cluster-Update.
/** This function updates the positive separable boundary-MPSs for the norm of PEPS0, which are located outside the cluster of size clusterSize, after the update of
    a boundary column or row, if clusterSize==0.
    If clusterSize!=0, this function does nothing.
    If WhichBoundary==left, then the positive separable boundary-MPS for the left boundary will be written into SeparableNormMPSs[0] with physical indices pointing right.
    If WhichBoundary==top, then the positive separable boundary-MPS for the top boundary will be written into SeparableNormMPSs[0] with physical indices pointing down.
    If WhichBoundary==right, then the positive separable boundary-MPS for the right boundary will be written into SeparableNormMPSs[Ncols-1] with physical indices pointing left.
    If WhichBoundary==bottom, then the positive separable boundary-MPS for the bottom boundary will be written into SeparableNormMPSs[Nrows-1] with physical indices pointing up.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices. */
  void updateSeparableNormMPSsBoundary(const PEPS<T>& PEPS0, const string& WhichBoundary,
                                       unsigned int clusterSize, double eps, unsigned int maxNumSweeps,
                                       vector< MPS<T> >& SeparableNormMPSs) const;

/// Updates positive separable norm-MPSs in bulk for Cluster-Update.
/** This function updates the positive separable boundary-MPSs for the norm of PEPS0, which are located outside the cluster of size clusterSize, after the update of
    a bulk column or row. 
    If SweepDirection==right, then the positive separable boundary-MPS for column updatePosition-clusterSize will be written into
    SeparableNormMPSs[updatePosition-clusterSize] with physical indices pointing right.
    If SweepDirection==down, then the positive separable boundary-MPS for row updatePosition-clusterSize will be written into
    SeparableNormMPSs[updatePosition-clusterSize] with physical indices pointing down.
    If SweepDirection==left, then the positive separable boundary-MPS for column updatePosition+clusterSize will be written into
    SeparableNormMPSs[updatePosition+clusterSize] with physical indices pointing left.
    If SweepDirection==up, then the positive separable boundary-MPS for row updatePosition+clusterSize will be written into
    SeparableNormMPSs[updatePosition+clusterSize] with physical indices pointing up.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices. */
  void updateSeparableNormMPSsBulk(const PEPS<T>& PEPS0, const string& SweepDirection, unsigned int updatePosition,
                                   unsigned int clusterSize, double eps, unsigned int maxNumSweeps,
                                   vector< MPS<T> >& SeparableNormMPSs) const;

/// Returns norm-MPS for boundary of PEPS0 for Cluster-Update.
/** This function computes the norm-MPS for a boundary of the norm of PEPS0 for the Cluster-Update.
    The vector D2s specifies the bond dimensions of the boundary-MPSs in the cluster, ordered from the closest column or row to the farthest.
    If the cluster-size given by D2s.size() is smaller than the system-size, then the norm-MPS is obtained starting from a SeparableNormMPS.
    If WhichBoundary==left, then NormMPS will be the boundary-MPS on column 1 with physical indices pointing left.
    If WhichBoundary==top, then NormMPS will be the boundary-MPS on row 1 with physical indices pointing up.
    If WhichBoundary==right, then NormMPS will be the boundary-MPS on column Ncols-2 with physical indices pointing right.
    If WhichBoundary==bottom, then NormMPS will be the boundary-MPS on row Nrows-2 with physical indices pointing down.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices. */
  void getBoundaryNormMPS(const PEPS<T>& PEPS0, const string& WhichBoundary,
                          const vector<unsigned int>& D2s, double eps, unsigned int maxNumSweeps,
                          const vector< MPS<T> >& SeparableNormMPSs,
                          MPS<T>& NormMPS) const;

/// Returns norm-MPSs for bulk of PEPS0 for Cluster-Update.
/** This function computes the norm-MPSs for the bulk of the norm of PEPS0 for the Cluster-Update.
    The vector D2s specifies the bond dimensions of the boundary-MPSs in the cluster, ordered from the closest column or row to the farthest.
    If the cluster-size given by D2s.size() is smaller than the system-size, then a norm-MPS can be obtained starting from a SeparableNormMPS.
    If UpdateDirection==vertical, then NormMPS1 will be the left and NormMPS2 the right NormMPS to column updatePosition. The physical indices
    of NormMPS1 point right and the physical indices of NormMPS2 point left.
    If UpdateDirection==horizontal, then NormMPS1 will be the lower and NormMPS2 the upper NormMPS to row updatePosition. The physical indices
    of NormMPS1 point up and the physical indices of NormMPS2 point down.
    As always, the tensors in an MPS are enumerated in the natural way, given by the direction of their physical indices. */
  void getBulkNormMPSs(const PEPS<T>& PEPS0, const string& UpdateDirection, unsigned int updatePosition,
                       const vector<unsigned int>& D2s, double eps, unsigned int maxNumSweeps,
                       const vector< MPS<T> >& SeparableNormMPSs,
                       MPS<T>& NormMPS1, MPS<T>& NormMPS2) const;

/// Updates PEPS boundary.
/** This function updates a PEPS boundary in time evolution, starting from tensorPositionStart and stopping before tensorPositionStop after applying
    the Trotter gates TEMPOs.
    When the tensors are updated independently, then the updated tensors are written into PEPS1. When they are updated sequentially, then the updated tensors are written into PEPS0.
    The bTensor can be constructed without or with all Trotter gates in the environment.
    The local tensor update can be reduced from the full tensor to the smallest subtensor including the physical bond.
    Either the pseudoinverse, the closest positive norm-matrix, or the gauge is used in the local update algorithm.
    \param WhichBoundary input: const string&, the boundary, must be "left", "top", "right", or "bottom"
    \param tensorPositionStart input: unsigned int, the start position of the top or left tensor of the first tensor-pair to be updated
    \param tensorPositionStop input: unsigned int, the stop position of the bottom or right tensor of the last tensor-pair to be updated
    \param TEMPOs input: const vector< MPO<T> >&, the Trotter gates
    \param SweepUpdate input: const string&, if SweepUpdate=="independent" then the updated tensors are written into PEPS1,
                              if SweepUpdate=="sequential" then the updated tensors are written into PEPS0
    \param NormMPSs input: const vector< MPS<T> >&, the boundary-MPSs for the norm of PEPS0
    \param bEnvironment: const string&, without or with the Trotter gates, must be "simplified" or "full"
    \param bMPSs input: const vector< MPS<T> >&, the boundary-MPSs for b if bEnvironment=="full"
    \param UpdateTensor input: const string&, the tensor in the local tensor update, must be "reduced" or "full"
    \param UpdateMode input: const string&, the local update algorithm, must be "Pseudoinverse" for pseudoinverse,
                             "PositiveNormMatrix" for the closest positive norm-matrix, or "Gauge"
    \param epsLoc input: double, the convergence precision in the local update
    \param maxNumSweepsLoc input: unsigned int, the maximal number of sweeps allowed in the local update
    \param cutoff input: double, the cutoff, singular values or eigenvalues x_{i} with (x_{i}/x_{max} <= cutoff) are put to zero
    \param PEPS0 input/output: PEPS<T>&, on input the initial PEPS, if SweepUpdate=="sequential" on output the final PEPS
    \param PEPS1 output: PEPS<T>&, if SweepUpdate=="independent" on output the final PEPS */
  void updatePEPSBoundary(const string& WhichBoundary, unsigned int tensorPositionStart, unsigned int tensorPositionStop,
                          const vector< MPO<T> >& TEMPOs, const string& SweepUpdate,
                          const vector< MPS<T> >& NormMPSs, const string& bEnvironment, const vector< MPS<T> >& bMPSs,
                          const string& UpdateTensor, const string& UpdateMode, double epsLoc, unsigned int maxNumSweepsLoc, double cutoff,
                          PEPS<T>& PEPS0, PEPS<T>& PEPS1) const;

/// Updates PEPS bulk.
/** This function updates the bulk of a PEPS in time evolution, starting from tensorPositionStart and stopping before tensorPositionStop after applying
    the Trotter gates TEMPOs.
    When the tensors are updated independently, then the updated tensors are written into PEPS1. When they are updated sequentially, then the updated tensors are written into PEPS0.
    The bTensor can be constructed without or with all Trotter gates in the environment.
    The local tensor update can be reduced from the full tensor to the smallest subtensor including the physical bond.
    Either the pseudoinverse, the closest positive norm-matrix, or the gauge is used in the local update algorithm.
    \param UpdateDirection input: const string&, must be "vertical" or "horizontal"
    \param updatePosition input: unsigned int, the column or row to be updated
    \param tensorPositionStart input: unsigned int, the start position of the top or left tensor of the first tensor-pair to be updated
    \param tensorPositionStop input: unsigned int, the stop position of the bottom or right tensor of the last tensor-pair to be updated
    \param TEMPOs input: const vector< MPO<T> >&, the Trotter gates
    \param SweepUpdate input: const string&, if SweepUpdate=="independent" then the updated tensors are written into PEPS1,
                              if SweepUpdate=="sequential" then the updated tensors are written into PEPS0
    \param NormMPSs input: const vector< MPS<T> >&, the boundary-MPSs for the norm of PEPS0
    \param bEnvironment: const string&, without or with the Trotter gates, must be "simplified" or "full"
    \param bMPSs input: const vector< MPS<T> >&, the boundary-MPSs for b if bEnvironment=="full"
    \param UpdateTensor input: const string&, the tensor in the local tensor update, must be "reduced" or "full"
    \param UpdateMode input: const string&, the local update algorithm, must be "Pseudoinverse" for pseudoinverse,
                             "PositiveNormMatrix" for the closest positive norm-matrix, or "Gauge"
    \param epsLoc input: double, the convergence precision in the local update
    \param maxNumSweepsLoc input: unsigned int, the maximal number of sweeps allowed in the local update
    \param cutoff input: double, the cutoff, singular values or eigenvalues x_{i} with (x_{i}/x_{max} <= cutoff) are put to zero
    \param PEPS0 input/output: PEPS<T>&, on input the initial PEPS, if SweepUpdate=="sequential" on output the final PEPS
    \param PEPS1 output: PEPS<T>&, if SweepUpdate=="independent" on output the final PEPS */
  void updatePEPSBulk(const string& UpdateDirection, unsigned int updatePosition, unsigned int tensorPositionStart, unsigned int tensorPositionStop,
                      const vector< MPO<T> >& TEMPOs, const string& SweepUpdate,
                      const vector< MPS<T> >& NormMPSs, const string& bEnvironment, const vector< MPS<T> >& bMPSs,
                      const string& UpdateTensor, const string& UpdateMode, double epsLoc, unsigned int maxNumSweepsLoc, double cutoff,
                      PEPS<T>& PEPS0, PEPS<T>& PEPS1) const;

/// Updates PEPS boundary with Cluster-Update.
/** This function updates a PEPS boundary in time evolution with the Cluster-Update, starting from tensorPositionStart and stopping before tensorPositionStop after applying
    the Trotter gates TEMPOs.
    When the tensors are updated independently, then the updated tensors are written into PEPS1. When they are updated sequentially, then the updated tensors are written into PEPS0.
    The vector D2sEnv specifies the bond dimensions of the boundary-MPSs in the cluster, ordered from the closest column or row to the farthest.
    If the cluster-size given by D2s.size() is smaller than the system-size, then the norm-MPS is obtained starting from a SeparableNormMPS and the b-MPS is obtained
    starting from a SeparablebMPS.
    The bTensor can be constructed without or with all Trotter gates in the environment.
    The local tensor update can be reduced from the full tensor to the smallest subtensor including the physical bond.
    Either the pseudoinverse, the closest positive norm-matrix, or the gauge is used in the local update algorithm.
    \param WhichBoundary input: const string&, the boundary, must be "left", "top", "right", or "bottom"
    \param tensorPositionStart input: unsigned int, the start position of the top or left tensor of the first tensor-pair to be updated
    \param tensorPositionStop input: unsigned int, the stop position of the bottom or right tensor of the last tensor-pair to be updated
    \param TEMPOs input: const vector< MPO<T> >&, the Trotter gates
    \param SweepUpdate input: const string&, if SweepUpdate=="independent" then the updated tensors are written into PEPS1,
                              if SweepUpdate=="sequential" then the updated tensors are written into PEPS0
    \param D2sEnv input: const vector<unsigned int>&, the maximal virtual bond dimensions for the environment approximation in the cluster,
                         sorted from the closest column or row to the farthest
    \param epsEnv input: double, the convergence precision in the environment approximation
    \param maxNumSweepsEnv input: unsigned int, the maximal number of sweeps allowed in the environment approximation
    \param SeparableNormMPSs input: const vector< MPS<T> >&, the positive separable boundary-MPSs for the norm of PEPS0 outside the cluster
    \param bEnvironment: const string&, without or with the Trotter gates, must be "simplified" or "full"
    \param SeparablebMPSs input: const vector< MPS<T> >&, the positive separable boundary-MPSs for b outside the cluster if bEnvironment=="full"
    \param UpdateTensor input: const string&, the tensor in the local tensor update, must be "reduced" or "full"
    \param UpdateMode input: const string&, the local update algorithm, must be "Pseudoinverse" for pseudoinverse,
                             "PositiveNormMatrix" for the closest positive norm-matrix, or "Gauge"
    \param epsLoc input: double, the convergence precision in the local update
    \param maxNumSweepsLoc input: unsigned int, the maximal number of sweeps allowed in the local update
    \param cutoff input: double, the cutoff, singular values or eigenvalues x_{i} with (x_{i}/x_{max} <= cutoff) are put to zero
    \param PEPS0 input/output: PEPS<T>&, on input the initial PEPS, if SweepUpdate=="sequential" on output the final PEPS
    \param PEPS1 output: PEPS<T>&, if SweepUpdate=="independent" on output the final PEPS */
  void clusterUpdatePEPSBoundary(const string& WhichBoundary, unsigned int tensorPositionStart, unsigned int tensorPositionStop,
                                 const vector< MPO<T> >& TEMPOs, const string& SweepUpdate,
                                 const vector<unsigned int>& D2sEnv, double epsEnv, unsigned int maxNumSweepsEnv,
                                 const vector< MPS<T> >& SeparableNormMPSs, const string& bEnvironment, const vector< MPS<T> >& SeparablebMPSs,
                                 const string& UpdateTensor, const string& UpdateMode, double epsLoc, unsigned int maxNumSweepsLoc, double cutoff,
                                 PEPS<T>& PEPS0, PEPS<T>& PEPS1) const;

/// Updates PEPS bulk with Cluster-Update.
/** This function updates the bulk of a PEPS in time evolution with the Cluster-Update, starting from tensorPositionStart and stopping before tensorPositionStop after applying
    the Trotter gates TEMPOs.
    When the tensors are updated independently, then the updated tensors are written into PEPS1. When they are updated sequentially, then the updated tensors are written into PEPS0.
    The vector D2sEnv specifies the bond dimensions of the boundary-MPSs in the cluster, ordered from the closest column or row to the farthest.
    If the cluster-size given by D2s.size() is smaller than the system-size, then the norm-MPS can be obtained starting from a SeparableNormMPS and the b-MPS can be obtained
    starting from a SeparablebMPS.
    The bTensor can be constructed without or with all Trotter gates in the environment.
    The local tensor update can be reduced from the full tensor to the smallest subtensor including the physical bond.
    Either the pseudoinverse, the closest positive norm-matrix, or the gauge is used in the local update algorithm.
    \param UpdateDirection input: const string&, must be "vertical" or "horizontal"
    \param updatePosition input: unsigned int, the column or row to be updated
    \param tensorPositionStart input: unsigned int, the start position of the top or left tensor of the first tensor-pair to be updated
    \param tensorPositionStop input: unsigned int, the stop position of the bottom or right tensor of the last tensor-pair to be updated
    \param TEMPOs input: const vector< MPO<T> >&, the Trotter gates
    \param SweepUpdate input: const string&, if SweepUpdate=="independent" then the updated tensors are written into PEPS1,
                              if SweepUpdate=="sequential" then the updated tensors are written into PEPS0
    \param D2sEnv input: const vector<unsigned int>&, the maximal virtual bond dimensions for the environment approximation in the cluster,
                         sorted from the closest column or row to the farthest
    \param epsEnv input: double, the convergence precision in the environment approximation
    \param maxNumSweepsEnv input: unsigned int, the maximal number of sweeps allowed in the environment approximation
    \param SeparableNormMPSs input: const vector< MPS<T> >&, the positive separable boundary-MPSs for the norm of PEPS0 outside the cluster
    \param bEnvironment: const string&, without or with the Trotter gates, must be "simplified" or "full"
    \param SeparablebMPSs input: const vector< MPS<T> >&, the positive separable boundary-MPSs for b outside the cluster if bEnvironment=="full"
    \param UpdateTensor input: const string&, the tensor in the local tensor update, must be "reduced" or "full"
    \param UpdateMode input: const string&, the local update algorithm, must be "Pseudoinverse" for pseudoinverse,
                             "PositiveNormMatrix" for the closest positive norm-matrix, or "Gauge"
    \param epsLoc input: double, the convergence precision in the local update
    \param maxNumSweepsLoc input: unsigned int, the maximal number of sweeps allowed in the local update
    \param cutoff input: double, the cutoff, singular values or eigenvalues x_{i} with (x_{i}/x_{max} <= cutoff) are put to zero
    \param PEPS0 input/output: PEPS<T>&, on input the initial PEPS, if SweepUpdate=="sequential" on output the final PEPS
    \param PEPS1 output: PEPS<T>&, if SweepUpdate=="independent" on output the final PEPS */
  void clusterUpdatePEPSBulk(const string& UpdateDirection, unsigned int updatePosition, unsigned int tensorPositionStart, unsigned int tensorPositionStop,
                             const vector< MPO<T> >& TEMPOs, const string& SweepUpdate,
                             const vector<unsigned int>& D2sEnv, double epsEnv, unsigned int maxNumSweepsEnv,
                             const vector< MPS<T> >& SeparableNormMPSs, const string& bEnvironment, const vector< MPS<T> >& SeparablebMPSs,
                             const string& UpdateTensor, const string& UpdateMode, double epsLoc, unsigned int maxNumSweepsLoc, double cutoff,
                             PEPS<T>& PEPS0, PEPS<T>& PEPS1) const;

/// Updates tensor pair.
/** This function updates a tensor pair on which TEMPO acts, with the first tensor at row positionRow and column
    positionCol and the second tensor in UpdateDirection to it.
    \param positionRow input: unsigned int, the row position of the first tensor
    \param positionCol input: unsigned int, the column position of the first tensor
    \param UpdateDirection input: const string&, the direction from the first to the second tensor, must be "vertical" or "horizontal"
    \param NormMPO input: MPO<T>&, the periodic norm-MPO,
                          if (UpdateDirection=="vertical") then the first tensor of NormMPO must be located on top of the tensor pair,
                          if (UpdateDirection=="horizontal") then the first tensor of NormMPO must be located left of the tensor pair,
                          the remaining tensors of NormMPO must be ordered clockwise around the tensor pair
    \param bEnvironment: const string&, without or with the Trotter gates, must be "simplified" or "full"
    \param bMPO input: const MPO<T>&, addressed only if (bEnvironment=="full"), the periodic b-MPO,
                       if (UpdateDirection=="vertical") then the first tensor of bMPO must be located on top of the tensor pair,
                       if (UpdateDirection=="horizontal") then the first tensor of bMPO must be located left of the tensor pair,
                       the remaining tensors of bMPO must be ordered clockwise around the tensor pair
    \param TEMPO input: const MPO<T>&, the time evolution MPO, must have the correct form
    \param UpdateTensor input: const string&, the tensor in the local tensor update, must be "reduced" or "full"
    \param UpdateMode input: const string&, the local update algorithm, must be "Pseudoinverse" for pseudoinverse,
                             "PositiveNormMatrix" for the closest positive norm-matrix, or "Gauge"
    \param eps input: double, the convergence precision
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param cutoff input: double, the cutoff, singular values or eigenvalues x_{i} with (x_{i}/x_{max} <= cutoff) are put to zero
    \param PEPS0 input/output: PEPS<T>&, on input the initial PEPS, on output the final PEPS */
  void updateTensor(unsigned int positionRow, unsigned int positionCol, const string& UpdateDirection,
                    MPO<T>& NormMPO, const string& bEnvironment, const MPO<T>& bMPO, const MPO<T>& TEMPO,
                    const string& UpdateTensor, const string& UpdateMode, double eps, unsigned int maxNumSweeps, double cutoff,
                    PEPS<T>& PEPS0) const;

/// Updates tensor pair with TEBD.
/** This function updates a tensor pair on which TEMPO acts, with the left tensor at row positionRow and column
    positionCol and the right tensor in Direction to it. It performs the simplified update and is used by the
    function
       void TimeEvolution2D<T>::timeEvolveTEBD(const string& Type, const string& UpdateTensor,
                                               Hamiltonian2D<T>& H0, const string& RealImaginaryTE,
                                               double timeStart, double timeStep, double timeStop,
                                               PEPS<T>& PEPS0, Matrix< Matrix<T> >& LambdasV,
                                               Matrix< Matrix<T> >& LambdasH) const   .
    \param UpdateTensor input: const string&, the tensor in the local tensor update, must be "reduced" or "full"
    \param positionRow input: unsigned int, the row position of the left tensor
    \param positionCol input: unsigned int, the column position of the left tensor
    \param Direction input: const string&, the direction to the right tensor, must be "vertical" or "horizontal"
    \param TEMPO input: const MPO<T>&, the time evolution MPO, must have the correct form
    \param PEPS0 input/output: PEPS<T>&, on input the initial PEPS, on output the final PEPS
    \param LambdasV input/output: Matrix< Matrix<T> >&, lambda matrices on the vertical bonds
    \param LambdasH input/output: Matrix< Matrix<T> >&, lambda matrices on the horizontal bonds
    \sa void TimeEvolution2D<T>::timeEvolveTEBD(const string& Type, const string& UpdateTensor,
                                                Hamiltonian2D<T>& H0, const string& RealImaginaryTE,
                                                double timeStart, double timeStep, double timeStop,
                                                PEPS<T>& PEPS0, Matrix< Matrix<T> >& LambdasV,
                                                Matrix< Matrix<T> >& LambdasH) const */
  void updateTensorTEBD(const string& UpdateTensor, unsigned int positionRow, unsigned int positionCol,
                        const string& Direction, const MPO<T>& TEMPO,
                        PEPS<T>& PEPS0, Matrix< Matrix<T> >& LambdasV, Matrix< Matrix<T> >& LambdasH) const;

/// Computes expectation value for PEPS boundary.
/** This function computes the expectation value of HMPOs on a PEPS boundary, starting from tensorPositionStart and stopping before tensorPositionStop,
    and returns that value divided by the norm.
    \param PEPS0 input: const PEPS<T>&, the PEPS
    \param WhichBoundary input: const string&, the boundary, must be "left", "top", "right", or "bottom"
    \param tensorPositionStart input: unsigned int, the start position of the top or left tensor of the first tensor pair
    \param tensorPositionStop input: unsigned int, the stop position of the bottom or right tensor of the last tensor pair
    \param HMPOs input: const vector< MPO<T> >&, the Hamiltonian MPOs
    \param NormMPSs input: const vector< MPS<T> >&, the boundary-MPSs for the norm of PEPS0
    \return T, the resulting expectation value divided by the norm */
  T expectationValuePEPSBoundary(const PEPS<T>& PEPS0, const string& WhichBoundary, unsigned int tensorPositionStart, unsigned int tensorPositionStop,
                                 const vector< MPO<T> >& HMPOs, const vector< MPS<T> >& NormMPSs) const;

/// Computes expectation value for PEPS bulk.
/** This function computes the expectation value of HMPOs in the bulk of a PEPS, starting from tensorPositionStart and stopping before tensorPositionStop,
    and returns that value divided by the norm.
    \param PEPS0 input: const PEPS<T>&, the PEPS
    \param Direction input: const string&, the direction of HMPOs, must be "vertical" or "horizontal"
    \param position input: unsigned int, the column or row position
    \param tensorPositionStart input: unsigned int, the start position of the top or left tensor of the first tensor pair
    \param tensorPositionStop input: unsigned int, the stop position of the bottom or right tensor of the last tensor pair
    \param HMPOs input: const vector< MPO<T> >&, the Hamiltonian MPOs
    \param NormMPSs input: const vector< MPS<T> >&, the boundary-MPSs for the norm of PEPS0
    \return T, the resulting expectation value divided by the norm */
  T expectationValuePEPSBulk(const PEPS<T>& PEPS0, const string& Direction, unsigned int position, unsigned int tensorPositionStart, unsigned int tensorPositionStop,
                             const vector< MPO<T> >& HMPOs, const vector< MPS<T> >& NormMPSs) const;

/// Computes cluster expectation value for PEPS boundary.
/** This function computes the cluster expectation value of HMPOs on a PEPS boundary, starting from tensorPositionStart and stopping before tensorPositionStop,
    and returns that value divided by the cluster norm.
    The vector D2s specifies the bond dimensions of the boundary-MPSs in the cluster, ordered from the closest column or row to the farthest.
    If the cluster-size given by D2s.size() is smaller than the system-size, then the norm-MPS is obtained starting from a SeparableNormMPS.
    \param PEPS0 input: const PEPS<T>&, the PEPS
    \param WhichBoundary input: const string&, the boundary, must be "left", "top", "right", or "bottom"
    \param tensorPositionStart input: unsigned int, the start position of the top or left tensor of the first tensor pair
    \param tensorPositionStop input: unsigned int, the stop position of the bottom or right tensor of the last tensor pair
    \param HMPOs input: const vector< MPO<T> >&, the Hamiltonian MPOs
    \param D2s input: const vector<unsigned int>&, the maximal virtual bond dimensions for the environment approximation in the cluster,
                      sorted from the closest column or row to the farthest
    \param eps input: double, the convergence precision in the environment approximation
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed in the environment approximation
    \param SeparableNormMPSs input: const vector< MPS<T> >&, the positive separable boundary-MPSs for the norm of PEPS0 outside the cluster
    \return T, the resulting cluster expectation value divided by the cluster norm */
  T clusterExpectationValuePEPSBoundary(const PEPS<T>& PEPS0, const string& WhichBoundary, unsigned int tensorPositionStart, unsigned int tensorPositionStop,
                                        const vector< MPO<T> >& HMPOs, 
                                        const vector<unsigned int>& D2s, double eps, unsigned int maxNumSweeps,
                                        const vector< MPS<T> >& SeparableNormMPSs) const;

/// Computes cluster expectation value for PEPS bulk.
/** This function computes the cluster expectation value of HMPOs in the bulk of a PEPS, starting from tensorPositionStart and stopping before tensorPositionStop,
    and returns that value divided by the cluster norm.
    The vector D2s specifies the bond dimensions of the boundary-MPSs in the cluster, ordered from the closest column or row to the farthest.
    If the cluster-size given by D2s.size() is smaller than the system-size, then the norm-MPS is obtained starting from a SeparableNormMPS.
    \param PEPS0 input: const PEPS<T>&, the PEPS
    \param Direction input: const string&, the direction of HMPOs, must be "vertical" or "horizontal"
    \param position input: unsigned int, the column or row position
    \param tensorPositionStart input: unsigned int, the start position of the top or left tensor of the first tensor pair
    \param tensorPositionStop input: unsigned int, the stop position of the bottom or right tensor of the last tensor pair
    \param HMPOs input: const vector< MPO<T> >&, the Hamiltonian MPOs
    \param D2s input: const vector<unsigned int>&, the maximal virtual bond dimensions for the environment approximation in the cluster,
                      sorted from the closest column or row to the farthest
    \param eps input: double, the convergence precision in the environment approximation
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed in the environment approximation
    \param SeparableNormMPSs input: const vector< MPS<T> >&, the positive separable boundary-MPSs for the norm of PEPS0 outside the cluster
    \return T, the resulting cluster expectation value divided by the cluster norm */
  T clusterExpectationValuePEPSBulk(const PEPS<T>& PEPS0, const string& Direction, unsigned int position, unsigned int tensorPositionStart, unsigned int tensorPositionStop,
                                    const vector< MPO<T> >& HMPOs,
                                    const vector<unsigned int>& D2s, double eps, unsigned int maxNumSweeps,
                                    const vector< MPS<T> >& SeparableNormMPSs) const;

/// Computes expectation value.
/** This function computes the expectation value of HMPO at the tensor pair, with the first tensor at row positionRow and column
    positionCol and the second tensor in Direction to it, and returns that value divided by the norm.
    \param PEPS0 input: const PEPS<T>&, the PEPS
    \param positionRow input: unsigned int, the row position of the first tensor
    \param positionCol input: unsigned int, the column position of the first tensor
    \param Direction input: const string&, the direction from the first to the second tensor, must be "vertical" or "horizontal"
    \param HMPO input: const MPO<T>&, the Hamiltonian MPO
    \param NormMPO input: const MPO<T>&, the periodic norm-MPO,
                          if (Direction=="vertical") then the first tensor of NormMPO must be located on top of the tensor pair,
                          if (Direction=="horizontal") then the first tensor of NormMPO must be located left of the tensor pair,
                          the remaining tensors of NormMPO must be ordered clockwise around the tensor pair
    \return T, the resulting expectation value divided by the norm */
  T expectationValue(const PEPS<T>& PEPS0, unsigned int positionRow, unsigned int positionCol, const string& Direction,
                     const MPO<T>& HMPO, const MPO<T>& NormMPO) const;

};

template<class T> void TimeEvolution2D<T>::timeEvolve(Hamiltonian2D<T>& H0, const string& RealImaginaryTE,
                                                      double timeStart, double timeStep, double timeStop,
                                                      double eps, unsigned int maxNumSweeps,
                                                      unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv,
                                                      const string& UpdateTensor, const string& UpdateMode, double cutoff,
                                                      PEPS<T>& PEPS0, unsigned int M) const
{
 unsigned int numSweepsDone;
 double epsAchieved, time;
 string Direction = "horizontal";
 PEPS<T> PEPS1;
 vector< PEPO<T> > TEPEPOs;
// 0. initialize:
// 0 - normalize PEPS0:
 PEPS0.normalize(Direction, D2Env, epsEnv, maxNumSweepsEnv);
// 0 - time-independent H0: get (Concatenated) TEPEPOs
 if (!(H0.isTimeDependent()))
  H0.getTEPEPOs(RealImaginaryTE, timeStep, TEPEPOs);
 if (M > 1)
 {
  for (int i = 0; i < TEPEPOs.size(); i++)
   TEPEPOs[i].setConcatenated(M);
 }
// 1. time-evolve:
 for (time = timeStart; time < timeStop; time += timeStep)
 {
// 1.1 time-dependent H0: get (Concatenated) TEPEPOs:
  if (H0.isTimeDependent())
  {
   H0.setTime(time);
   H0.getTEPEPOs(RealImaginaryTE, timeStep, TEPEPOs);
   if (M > 1)
   {
    for (int i = 0; i < TEPEPOs.size(); i++)
     TEPEPOs[i].setConcatenated(M);
   }
  }
// 1.2 approximate product of TEPEPOs with PEPS0:
  for (int i = 0; i < TEPEPOs.size(); i++)
  {
   PEPS1 = PEPS0;
   multiplyPEPOPEPS(TEPEPOs[i], PEPS1, eps, maxNumSweeps, D2Env, epsEnv, maxNumSweepsEnv,
                    UpdateTensor, UpdateMode, cutoff, epsAchieved, numSweepsDone, PEPS0);
   PEPS0.normalize(Direction, D2Env, epsEnv, maxNumSweepsEnv);
  }
 }
// 2. time-dependent H0: reset time
 if (H0.isTimeDependent())
  H0.setTime(timeStart);
}

template<class T> void TimeEvolution2D<T>::getInitialNormMPSs(const PEPS<T>& PEPS0, const string& UpdateDirection,
                                                              unsigned int D2, double eps, unsigned int maxNumSweeps,
                                                              vector< MPS<T> >& NormMPSs) const
{
 unsigned int D = PEPS0.getD(), Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), numSweepsDone;
 double epsAchieved;
 string BC = "open", Direction, WhichBoundary;
 vector<unsigned int> Seed;
// a) vertical updates:
 if (UpdateDirection == "vertical")
 {
  NormMPSs = vector< MPS<T> >(Ncols);
// a)1. boundary-MPS for boundary:
  NormMPSs[Ncols-1] = MPS<T>(BC, Nrows, D*D, min(D2, D*D));
  Seed = vector<unsigned int>(Nrows);
  for (int i = 0; i < Nrows; i++)
   Seed[i] = time(0)+13*i;
  NormMPSs[Ncols-1].fillRandomly(Seed);
  WhichBoundary = "right";
  PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPSs[Ncols-1]);
// a)2. boundary-MPSs for bulk:
  Direction = "left";
  for (int position = Ncols-2; position > 0; position--)
  {
   NormMPSs[position] = NormMPSs[position+1];
   NormMPSs[position].setD(min(D2, D*D*NormMPSs[position+1].getD()));
   PEPS0.multiplyBulkMPOMPS(Direction, position, NormMPSs[position+1], eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPSs[position]);
  }
 }
// b) horizontal updates:
 else if (UpdateDirection == "horizontal")
 {
  NormMPSs = vector< MPS<T> >(Nrows);
// b)1. boundary-MPS for boundary:
  NormMPSs[Nrows-1] = MPS<T>(BC, Ncols, D*D, min(D2, D*D));
  Seed = vector<unsigned int>(Ncols);
  for (int i = 0; i < Ncols; i++)
   Seed[i] = time(0)+13*i;
  NormMPSs[Nrows-1].fillRandomly(Seed);
  WhichBoundary = "bottom";
  PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPSs[Nrows-1]);
// b)2. boundary-MPSs for bulk:
  Direction = "up";
  for (int position = Nrows-2; position > 0; position--)
  {
   NormMPSs[position] = NormMPSs[position+1];
   NormMPSs[position].setD(min(D2, D*D*NormMPSs[position+1].getD()));
   PEPS0.multiplyBulkMPOMPS(Direction, position, NormMPSs[position+1], eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPSs[position]);
  }
 }
}

template<class T> void TimeEvolution2D<T>::updateNormMPSsBoundary(const PEPS<T>& PEPS0, const string& WhichBoundary,
                                                                  unsigned int D2, double eps, unsigned int maxNumSweeps,
                                                                  vector< MPS<T> >& NormMPSs) const
{
 unsigned int D = PEPS0.getD(), Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), numSweepsDone;
 double epsAchieved;
 string BC = "open";
 vector<unsigned int> Seed;
// a) left boundary:
 if (WhichBoundary == "left")
 {
  NormMPSs[0] = MPS<T>(BC, Nrows, D*D, min(D2, D*D));
  Seed = vector<unsigned int>(Nrows);
  for (int i = 0; i < Nrows; i++)
   Seed[i] = time(0)+13*i;
  NormMPSs[0].fillRandomly(Seed);
  PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPSs[0]);
 }
// b) top boundary:
 else if (WhichBoundary == "top")
 {
  NormMPSs[0] = MPS<T>(BC, Ncols, D*D, min(D2, D*D));
  Seed = vector<unsigned int>(Ncols);
  for (int i = 0; i < Ncols; i++)
   Seed[i] = time(0)+13*i;
  NormMPSs[0].fillRandomly(Seed);
  PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPSs[0]);
 }
// c) right boundary:
 else if (WhichBoundary == "right")
 {
  NormMPSs[Ncols-1] = MPS<T>(BC, Nrows, D*D, min(D2, D*D));
  Seed = vector<unsigned int>(Nrows);
  for (int i = 0; i < Nrows; i++)
   Seed[i] = time(0)+13*i;
  NormMPSs[Ncols-1].fillRandomly(Seed);
  PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPSs[Ncols-1]);
 }
// d) bottom boundary:
 else if (WhichBoundary == "bottom")
 {
  NormMPSs[Nrows-1] = MPS<T>(BC, Ncols, D*D, min(D2, D*D));
  Seed = vector<unsigned int>(Ncols);
  for (int i = 0; i < Ncols; i++)
   Seed[i] = time(0)+13*i;
  NormMPSs[Nrows-1].fillRandomly(Seed);
  PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPSs[Nrows-1]);
 }
}

template<class T> void TimeEvolution2D<T>::updateNormMPSsBulk(const PEPS<T>& PEPS0, const string& SweepDirection, unsigned int updatePosition,
                                                              unsigned int D2, double eps, unsigned int maxNumSweeps,
                                                              vector< MPS<T> >& NormMPSs) const
{
 unsigned int D = PEPS0.getD(), numSweepsDone;
 double epsAchieved;
// a) sweep right:
 if (SweepDirection == "right")
 {
  NormMPSs[updatePosition] = NormMPSs[updatePosition-1];
  NormMPSs[updatePosition].setD(min(D2, D*D*NormMPSs[updatePosition-1].getD()));
  PEPS0.multiplyBulkMPOMPS(SweepDirection, updatePosition, NormMPSs[updatePosition-1], eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPSs[updatePosition]);
 }
// b) sweep down:
 else if (SweepDirection == "down")
 {
  NormMPSs[updatePosition] = NormMPSs[updatePosition-1];
  NormMPSs[updatePosition].setD(min(D2, D*D*NormMPSs[updatePosition-1].getD()));
  PEPS0.multiplyBulkMPOMPS(SweepDirection, updatePosition, NormMPSs[updatePosition-1], eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPSs[updatePosition]);
 }
// c) sweep left:
 else if (SweepDirection == "left")
 {
  NormMPSs[updatePosition] = NormMPSs[updatePosition+1];
  NormMPSs[updatePosition].setD(min(D2, D*D*NormMPSs[updatePosition+1].getD()));
  PEPS0.multiplyBulkMPOMPS(SweepDirection, updatePosition, NormMPSs[updatePosition+1], eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPSs[updatePosition]);
 }
// d) sweep up:
 else if (SweepDirection == "up")
 {
  NormMPSs[updatePosition] = NormMPSs[updatePosition+1];
  NormMPSs[updatePosition].setD(min(D2, D*D*NormMPSs[updatePosition+1].getD()));
  PEPS0.multiplyBulkMPOMPS(SweepDirection, updatePosition, NormMPSs[updatePosition+1], eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPSs[updatePosition]);
 }
}

template<class T> void TimeEvolution2D<T>::getInitialBoundaryNormTensors(const PEPS<T>& PEPS0, const string& WhichBoundary,
                                                                         const MPS<T>& NormMPS,
                                                                         vector< Tensor<T> >& NormTensors) const
{
 unsigned int D = PEPS0.getD(), N, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows();
 vector<unsigned int> Indices0(1), Indices02(2), Indices03(3), Indices1(1), Indices12(2), Indices13(3), Order4(4), Shape3(3), Shape4(4), Shape6(6), Shape8(8);
 Tensor<T> Tensor0, Tensor1;
// a) left boundary column:
 if (WhichBoundary == "left")
 {
  N = Nrows;
  NormTensors = vector< Tensor<T> >(N);
// a)1. position N-1:
  NormMPS.get(N-1, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  PEPS0.get(N-1, 0, Tensor1);
  Indices0[0] = 2; Indices1[0] = 2;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  PEPS0.get(N-1, 0, Tensor1); Tensor1.complexConjugate();
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 2; Indices12[1] = 4;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape8);
  Shape4[0] = Shape8[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormTensors[N-1] = Tensor0;
// a)2. N-1 > position > 1:
  Indices0[0] = 0; Indices1[0] = 1;
  Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 3; Indices12[1] = 2;
  Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 3; Indices13[1] = 2; Indices13[2] = 4;
  Order4[0] = 1; Order4[1] = 0; Order4[2] = 2; Order4[3] = 3;
  for (int position = N-2; position > 1; position--)
  {
   NormMPS.get(position, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(position, 0, Tensor1);
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(position, 0, Tensor1); Tensor1.complexConjugate();
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = 1; Shape4[1] = Shape6[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   Tensor0.permute(Order4);
   NormTensors[position] = Tensor0;
  }
// a)3. position 0:
  NormMPS.get(0, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  PEPS0.get(0, 0, Tensor1);
  Indices0[0] = 2; Indices1[0] = 2;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  PEPS0.get(0, 0, Tensor1); Tensor1.complexConjugate();
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 2; Indices12[1] = 4;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape8);
  Shape4[0] = 1; Shape4[1] = Shape8[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormTensors[0] = Tensor0;
 }
// b) top boundary row:
 else if (WhichBoundary == "top")
 {
  N = Ncols;
  NormTensors = vector< Tensor<T> >(N);
// b)1. position N-1:
  NormMPS.get(0, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  PEPS0.get(0, N-1, Tensor1);
  Indices0[0] = 2; Indices1[0] = 3;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  PEPS0.get(0, N-1, Tensor1); Tensor1.complexConjugate();
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 4;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape8);
  Shape4[0] = 1; Shape4[1] = Shape8[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormTensors[N-1] = Tensor0;
// b)2. N-1 > position > 1:
  Indices0[0] = 1; Indices1[0] = 0;
  Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 2; Indices12[1] = 3;
  Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 2; Indices13[1] = 3; Indices13[2] = 4;
  for (int position = N-2; position > 1; position--)
  {
   NormMPS.get(N-1-position, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(0, position, Tensor1);
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(0, position, Tensor1); Tensor1.complexConjugate();
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = 1; Shape4[1] = Shape6[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[position] = Tensor0;
  }
// b)3. position 0:
  NormMPS.get(N-1, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  PEPS0.get(0, 0, Tensor1);
  Indices0[0] = 2; Indices1[0] = 3;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  PEPS0.get(0, 0, Tensor1); Tensor1.complexConjugate();
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 4;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape8);
  Shape4[0] = Shape8[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormTensors[0] = Tensor0;
 }
// c) right boundary column:
 else if (WhichBoundary == "right")
 {
  N = Nrows;
  NormTensors = vector< Tensor<T> >(N);
// c)1. position N-1:
  NormMPS.get(0, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  PEPS0.get(N-1, Ncols-1, Tensor1);
  Indices0[0] = 2; Indices1[0] = 0;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  PEPS0.get(N-1, Ncols-1, Tensor1); Tensor1.complexConjugate();
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 4;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape8);
  Shape4[0] = 1; Shape4[1] = Shape8[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormTensors[N-1] = Tensor0;
// c)2. N-1 > position > 1:
  Indices0[0] = 1; Indices1[0] = 0;
  Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 3; Indices12[1] = 0;
  Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 3; Indices13[1] = 0; Indices13[2] = 4;
  for (int position = N-2; position > 1; position--)
  {
   NormMPS.get(N-1-position, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(position, Ncols-1, Tensor1);
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(position, Ncols-1, Tensor1); Tensor1.complexConjugate();
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = 1; Shape4[1] = Shape6[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[position] = Tensor0;
  }
// c)3. position 0:
  NormMPS.get(N-1, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  PEPS0.get(0, Ncols-1, Tensor1);
  Indices0[0] = 2; Indices1[0] = 0;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  PEPS0.get(0, Ncols-1, Tensor1); Tensor1.complexConjugate();
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 4;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape8);
  Shape4[0] = Shape8[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormTensors[0] = Tensor0;
 }
// d) NormTensors for bottom boundary row:
 else if (WhichBoundary == "bottom")
 {
  N = Ncols;
  NormTensors = vector< Tensor<T> >(N);
// d)1. position N-1:
  NormMPS.get(N-1, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  PEPS0.get(Nrows-1, N-1, Tensor1);
  Indices0[0] = 2; Indices1[0] = 1;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  PEPS0.get(Nrows-1, N-1, Tensor1); Tensor1.complexConjugate();
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 4;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape8);
  Shape4[0] = Shape8[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormTensors[N-1] = Tensor0;
// d)2. N-1 > position > 1:
  Indices0[0] = 0; Indices1[0] = 1;
  Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 2; Indices12[1] = 1;
  Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 2; Indices13[1] = 1; Indices13[2] = 4;
  Order4[0] = 1; Order4[1] = 0; Order4[2] = 2; Order4[3] = 3;
  for (int position = N-2; position > 1; position--)
  {
   NormMPS.get(position, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(Nrows-1, position, Tensor1);
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(Nrows-1, position, Tensor1); Tensor1.complexConjugate();
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = 1; Shape4[1] = Shape6[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   Tensor0.permute(Order4);
   NormTensors[position] = Tensor0;
  }
// d)3. position 0:
  NormMPS.get(0, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  PEPS0.get(Nrows-1, 0, Tensor1);
  Indices0[0] = 2; Indices1[0] = 1;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  PEPS0.get(Nrows-1, 0, Tensor1); Tensor1.complexConjugate();
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 4;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape8);
  Shape4[0] = 1; Shape4[1] = Shape8[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormTensors[0] = Tensor0;
 }
}

template<class T> void TimeEvolution2D<T>::getInitialBulkNormTensors(const PEPS<T>& PEPS0, const string& UpdateDirection, unsigned int updatePosition,
                                                                     const MPS<T>& NormMPS1, const MPS<T>& NormMPS2,
                                                                     vector< Tensor<T> >& NormTensors) const
{
 unsigned int D = PEPS0.getD(), N, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows();
 vector<unsigned int> Indices0(1), Indices02(2), Indices03(3), Indices1(1), Indices12(2), Indices13(3), Order4(4), Shape3(3), Shape4(4), Shape8(8);
 Tensor<T> Tensor0, Tensor1;
// a) vertical updates:
 if (UpdateDirection == "vertical")
 {
  N = Nrows;
  NormTensors = vector< Tensor<T> >(N);
// a)1. position N-1:
  NormMPS1.get(0, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  PEPS0.get(N-1, updatePosition, Tensor1);
  Indices0[0] = 2; Indices1[0] = 0;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  PEPS0.get(N-1, updatePosition, Tensor1); Tensor1.complexConjugate();
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 4;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  NormMPS2.get(N-1, Tensor1);
  Tensor1.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
  Tensor1.reshape(Shape4);
  Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 2; Indices12[1] = 3;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape8);
  Shape4[0] = Shape8[1]; Shape4[1] = D; Shape4[2] = D; Shape4[3] = Shape8[6];
  Tensor0.reshape(Shape4);
  Order4[0] = 3; Order4[1] = 0; Order4[2] = 1; Order4[3] = 2;
  Tensor0.permute(Order4);
  NormTensors[N-1] = Tensor0;
// a)2. N-1 > position > 1:
  Indices0[0] = 1; Indices1[0] = 0;
  Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 3; Indices12[1] = 0;
  Order4[0] = 3; Order4[1] = 0; Order4[2] = 1; Order4[3] = 2;
  for (int position = N-2; position > 1; position--)
  {
   NormMPS1.get(N-1-position, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(position, updatePosition, Tensor1);
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(position, updatePosition, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 3; Indices13[1] = 0; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   NormMPS2.get(position, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 5; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.permute(Order4);
   NormTensors[position] = Tensor0;
  }
// a)3. position 0:
  NormMPS1.get(N-1, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  PEPS0.get(0, updatePosition, Tensor1);
  Indices0[0] = 2; Indices1[0] = 0;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  PEPS0.get(0, updatePosition, Tensor1); Tensor1.complexConjugate();
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 4;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  NormMPS2.get(0, Tensor1);
  Tensor1.getShape(Shape3);
  Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor1.reshape(Shape4);
  Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 2; Indices12[1] = 3;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape8);
  Shape4[0] = Shape8[0]; Shape4[1] = D; Shape4[2] = D; Shape4[3] = Shape8[7];
  Tensor0.reshape(Shape4);
  Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
  Tensor0.permute(Order4);
  NormTensors[0] = Tensor0;
 }
// b) horizontal updates:
 else if (UpdateDirection == "horizontal")
 {
  N = Ncols;
  NormTensors = vector< Tensor<T> >(N);
// b)1. position N-1:
  NormMPS1.get(0, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  PEPS0.get(updatePosition, N-1, Tensor1);
  Indices0[0] = 2; Indices1[0] = 3;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  PEPS0.get(updatePosition, N-1, Tensor1); Tensor1.complexConjugate();
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 4;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  NormMPS2.get(N-1, Tensor1);
  Tensor1.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
  Tensor1.reshape(Shape4);
  Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 2; Indices12[1] = 3;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape8);
  Shape4[0] = Shape8[1]; Shape4[1] = D; Shape4[2] = D; Shape4[3] = Shape8[6];
  Tensor0.reshape(Shape4);
  Order4[0] = 3; Order4[1] = 0; Order4[2] = 1; Order4[3] = 2;
  Tensor0.permute(Order4);
  NormTensors[N-1] = Tensor0;
// b)2. N-1 > position > 1:
  Indices0[0] = 1; Indices1[0] = 0;
  Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 2; Indices12[1] = 3;
  Order4[0] = 3; Order4[1] = 0; Order4[2] = 1; Order4[3] = 2;
  for (int position = N-2; position > 1; position--)
  {
   NormMPS1.get(N-1-position, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(updatePosition, position, Tensor1);
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(updatePosition, position, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 2; Indices13[1] = 3; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   NormMPS2.get(position, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 5; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.permute(Order4);
   NormTensors[position] = Tensor0;
  }
// b)3. position 0:
  NormMPS1.get(N-1, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  PEPS0.get(updatePosition, 0, Tensor1);
  Indices0[0] = 2; Indices1[0] = 3;
  Tensor0.contract(Indices0, Tensor1, Indices1);
  PEPS0.get(updatePosition, 0, Tensor1); Tensor1.complexConjugate();
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 4;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  NormMPS2.get(0, Tensor1);
  Tensor1.getShape(Shape3);
  Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor1.reshape(Shape4);
  Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 2; Indices12[1] = 3;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  Tensor0.getShape(Shape8);
  Shape4[0] = Shape8[0]; Shape4[1] = D; Shape4[2] = D; Shape4[3] = Shape8[7];
  Tensor0.reshape(Shape4);
  Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
  Tensor0.permute(Order4);
  NormTensors[0] = Tensor0;
 }
}

template<class T> void TimeEvolution2D<T>::updateBoundaryNormTensors(const PEPS<T>& PEPS0, const string& WhichBoundary,
                                                                     const MPS<T>& NormMPS, unsigned int tensorPosition,
                                                                     vector< Tensor<T> >& NormTensors) const
{
 unsigned int D = PEPS0.getD(), N, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows();
 vector<unsigned int> Indices0(1), Indices02(2), Indices03(3), Indices1(1), Indices12(2), Indices13(3), Shape3(3), Shape4(4), Shape6(6), Shape8(8);
 Tensor<T> Tensor0, Tensor1;
// a) left boundary column:
 if (WhichBoundary == "left")
 {
  N = Nrows;
// a)a) tensorPosition == 0:
  if (tensorPosition == 0)
  {
   NormMPS.get(0, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   PEPS0.get(0, 0, Tensor1);
   Indices0[0] = 2; Indices1[0] = 2;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(0, 0, Tensor1); Tensor1.complexConjugate();
   Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 2; Indices12[1] = 4;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape8);
   Shape4[0] = 1; Shape4[1] = Shape8[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[0] = Tensor0;
   NormMPS.get(1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 1; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(1, 0, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(1, 0, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = 1; Shape4[1] = Shape6[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[1] = Tensor0;
  }
// a)b) 0 < tensorPosition < N-3:
  else if ((tensorPosition > 0) && (tensorPosition < N-3))
  {
   Tensor0 = NormTensors[tensorPosition-1];
   NormMPS.get(tensorPosition, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 1; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(tensorPosition, 0, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(tensorPosition, 0, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = 1; Shape4[1] = Shape6[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[tensorPosition] = Tensor0;
   NormMPS.get(tensorPosition+1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 1; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(tensorPosition+1, 0, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(tensorPosition+1, 0, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = 1; Shape4[1] = Shape6[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[tensorPosition+1] = Tensor0;
  }
 }
// b) top boundary row:
 else if (WhichBoundary == "top")
 {
  N = Ncols;
// b)a) tensorPosition == 0:
  if (tensorPosition == 0)
  {
   NormMPS.get(N-1, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   PEPS0.get(0, 0, Tensor1);
   Indices0[0] = 2; Indices1[0] = 3;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(0, 0, Tensor1); Tensor1.complexConjugate();
   Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 4;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape8);
   Shape4[0] = Shape8[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[0] = Tensor0;
   NormMPS.get(N-2, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(0, 1, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(0, 1, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 0; Indices13[1] = 3; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = Shape6[1]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[1] = Tensor0;
  }
// b)b) 0 < tensorPosition < N-3:
  else if ((tensorPosition > 0) && (tensorPosition < N-3))
  {
   Tensor0 = NormTensors[tensorPosition-1];
   NormMPS.get(N-1-tensorPosition, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(0, tensorPosition, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(0, tensorPosition, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 0; Indices13[1] = 3; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = Shape6[1]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[tensorPosition] = Tensor0;
   NormMPS.get(N-1-tensorPosition-1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(0, tensorPosition+1, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(0, tensorPosition+1, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 0; Indices13[1] = 3; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = Shape6[1]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[tensorPosition+1] = Tensor0;
  }
 }
// c) right boundary column:
 else if (WhichBoundary == "right")
 {
  N = Nrows;
// c)a) tensorPosition == 0:
  if (tensorPosition == 0)
  {
   NormMPS.get(N-1, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   PEPS0.get(0, Ncols-1, Tensor1);
   Indices0[0] = 2; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(0, Ncols-1, Tensor1); Tensor1.complexConjugate();
   Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 4;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape8);
   Shape4[0] = Shape8[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[0] = Tensor0;
   NormMPS.get(N-2, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(1, Ncols-1, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 0;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(1, Ncols-1, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 1; Indices13[1] = 0; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = Shape6[1]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[1] = Tensor0;
  }
// c)b) 0 < tensorPosition < N-3:
  else if ((tensorPosition > 0) && (tensorPosition < N-3))
  {
   Tensor0 = NormTensors[tensorPosition-1];
   NormMPS.get(N-1-tensorPosition, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(tensorPosition, Ncols-1, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 0;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(tensorPosition, Ncols-1, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 1; Indices13[1] = 0; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = Shape6[1]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[tensorPosition] = Tensor0;
   NormMPS.get(N-1-tensorPosition-1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(tensorPosition+1, Ncols-1, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 0;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(tensorPosition+1, Ncols-1, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 1; Indices13[1] = 0; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = Shape6[1]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[tensorPosition+1] = Tensor0;
  }
 }
// d) bottom boundary row:
 else if (WhichBoundary == "bottom")
 {
  N = Ncols;
// d)a) tensorPosition == 0:
  if (tensorPosition == 0)
  {
   NormMPS.get(0, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   PEPS0.get(Nrows-1, 0, Tensor1);
   Indices0[0] = 2; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(Nrows-1, 0, Tensor1); Tensor1.complexConjugate();
   Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 4;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape8);
   Shape4[0] = 1; Shape4[1] = Shape8[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[0] = Tensor0;
   NormMPS.get(1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 1; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(Nrows-1, 1, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 1;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(Nrows-1, 1, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = 1; Shape4[1] = Shape6[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[1] = Tensor0;
  }
// d)b) 0 < tensorPosition < N-3:
  else if ((tensorPosition > 0) && (tensorPosition < N-3))
  {
   Tensor0 = NormTensors[tensorPosition-1];
   NormMPS.get(tensorPosition, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 1; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(Nrows-1, tensorPosition, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 1;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(Nrows-1, tensorPosition, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = 1; Shape4[1] = Shape6[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[tensorPosition] = Tensor0;
   NormMPS.get(tensorPosition+1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 1; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(Nrows-1, tensorPosition+1, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 1;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(Nrows-1, tensorPosition+1, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 0; Indices13[1] = 1; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Tensor0.getShape(Shape6);
   Shape4[0] = 1; Shape4[1] = Shape6[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormTensors[tensorPosition+1] = Tensor0;
  }
 }
}

template<class T> void TimeEvolution2D<T>::updateBulkNormTensors(const PEPS<T>& PEPS0, const string& UpdateDirection, unsigned int updatePosition,
                                                                 const MPS<T>& NormMPS1, const MPS<T>& NormMPS2,
                                                                 unsigned int tensorPosition, vector< Tensor<T> >& NormTensors) const
{
 unsigned int D = PEPS0.getD(), N, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows();
 vector<unsigned int> Indices0(1), Indices02(2), Indices03(3), Indices1(1), Indices12(2), Indices13(3), Order4(4), Shape3(3), Shape4(4), Shape8(8);
 Tensor<T> Tensor0, Tensor1;
// a) vertical updates:
 if (UpdateDirection == "vertical")
 {
  N = Nrows;
// a)a) tensorPosition 0:
  if (tensorPosition == 0)
  {
   NormMPS1.get(N-1, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   PEPS0.get(0, updatePosition, Tensor1);
   Indices0[0] = 2; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(0, updatePosition, Tensor1); Tensor1.complexConjugate();
   Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 4;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   NormMPS2.get(0, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 2; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape8);
   Shape4[0] = Shape8[0]; Shape4[1] = D; Shape4[2] = D; Shape4[3] = Shape8[7];
   Tensor0.reshape(Shape4);
   Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
   Tensor0.permute(Order4);
   NormTensors[0] = Tensor0;
   NormMPS1.get(N-2, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(1, updatePosition, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 0;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(1, updatePosition, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 1; Indices13[1] = 0; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   NormMPS2.get(1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices03[0] = 0; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
   Tensor0.permute(Order4);
   NormTensors[1] = Tensor0;
  }
// a)b) 0 < tensorPosition < N-3:
  else if ((tensorPosition > 0) && (tensorPosition < N-3))
  {
   Tensor0 = NormTensors[tensorPosition-1];
   NormMPS1.get(N-1-tensorPosition, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(tensorPosition, updatePosition, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 0;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(tensorPosition, updatePosition, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 1; Indices13[1] = 0; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   NormMPS2.get(tensorPosition, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices03[0] = 0; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
   Tensor0.permute(Order4);
   NormTensors[tensorPosition] = Tensor0;
   NormMPS1.get(N-1-tensorPosition-1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(tensorPosition+1, updatePosition, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 0;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(tensorPosition+1, updatePosition, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 1; Indices13[1] = 0; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   NormMPS2.get(tensorPosition+1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices03[0] = 0; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
   Tensor0.permute(Order4);
   NormTensors[tensorPosition+1] = Tensor0;
  }
 }
// b) horizontal updates:
 else if (UpdateDirection == "horizontal")
 {
  N = Ncols;
// b)a) tensorPosition 0:
  if (tensorPosition == 0)
  {
   NormMPS1.get(N-1, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   PEPS0.get(updatePosition, 0, Tensor1);
   Indices0[0] = 2; Indices1[0] = 3;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(updatePosition, 0, Tensor1); Tensor1.complexConjugate();
   Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 4;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   NormMPS2.get(0, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 2; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   Tensor0.getShape(Shape8);
   Shape4[0] = Shape8[0]; Shape4[1] = D; Shape4[2] = D; Shape4[3] = Shape8[7];
   Tensor0.reshape(Shape4);
   Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
   Tensor0.permute(Order4);
   NormTensors[0] = Tensor0;
   NormMPS1.get(N-2, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(updatePosition, 1, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(updatePosition, 1, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 0; Indices13[1] = 3; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   NormMPS2.get(1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices03[0] = 0; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
   Tensor0.permute(Order4);
   NormTensors[1] = Tensor0;
  }
// b)b) 0 < tensorPosition < N-3:
  else if ((tensorPosition > 0) && (tensorPosition < N-3))
  {
   Tensor0 = NormTensors[tensorPosition-1];
   NormMPS1.get(N-1-tensorPosition, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(updatePosition, tensorPosition, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(updatePosition, tensorPosition, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 0; Indices13[1] = 3; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   NormMPS2.get(tensorPosition, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices03[0] = 0; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
   Tensor0.permute(Order4);
   NormTensors[tensorPosition] = Tensor0;
   NormMPS1.get(N-1-tensorPosition-1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   PEPS0.get(updatePosition, tensorPosition+1, Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 3;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   PEPS0.get(updatePosition, tensorPosition+1, Tensor1); Tensor1.complexConjugate();
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 0; Indices13[1] = 3; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   NormMPS2.get(tensorPosition+1, Tensor1);
   Tensor1.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor1.reshape(Shape4);
   Indices03[0] = 0; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
   Tensor0.permute(Order4);
   NormTensors[tensorPosition+1] = Tensor0;
  }
 }
}

template<class T> void TimeEvolution2D<T>::getBoundaryNormMPO(const string& WhichBoundary, const MPS<T>& NormMPS, const vector< Tensor<T> >& NormTensors,
                                                              unsigned int tensorPosition, MPO<T>& NormMPO) const
{
 unsigned int D = sqrt(NormMPS.getd()), N = NormMPS.getN();
 vector<unsigned int> Shape3(3), Shape4(4);
 Tensor<T> Tensor0;
// OneTensor:
 Shape4[0] = 1; Shape4[1] = 1; Shape4[2] = 1; Shape4[3] = 1;
 Tensor<T> OneTensor(Shape4);
 OneTensor.set(0, 1.0);
// NormMPO:
 NormMPO = MPO<T>("periodic", 6, D, NormMPS.getD());
// a) left boundary column or bottom boundary row:
 if ((WhichBoundary == "left") || (WhichBoundary == "bottom"))
 {
// a)a) tensorPosition == 0:
  if (tensorPosition == 0)
  {
   NormMPO.set(0, OneTensor);
   NormMPS.get(0, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormMPO.set(1, Tensor0);
   NormMPS.get(1, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormMPO.set(2, Tensor0);
   NormMPO.set(3, NormTensors[2]);
   NormMPO.set(4, OneTensor);
   NormMPO.set(5, OneTensor);
  }
// a)b) 0 < tensorPosition < N-2:
  else if ((tensorPosition > 0) && (tensorPosition < N-2))
  {
   NormMPO.set(0, NormTensors[tensorPosition-1]);
   NormMPS.get(tensorPosition, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormMPO.set(1, Tensor0);
   NormMPS.get(tensorPosition+1, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormMPO.set(2, Tensor0);
   NormMPO.set(3, NormTensors[tensorPosition+2]);
   NormMPO.set(4, OneTensor);
   NormMPO.set(5, OneTensor);
  }
// a)c) tensorPosition == N-2:
  else if (tensorPosition == N-2)
  {
   NormMPO.set(0, NormTensors[N-3]);
   NormMPS.get(N-2, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormMPO.set(1, Tensor0);
   NormMPS.get(N-1, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormMPO.set(2, Tensor0);
   NormMPO.set(3, OneTensor);
   NormMPO.set(4, OneTensor);
   NormMPO.set(5, OneTensor);
  }
 }
// b) top boundary row or right boundary column:
 else if ((WhichBoundary == "top") || (WhichBoundary == "right"))
 {
// b)a) tensorPosition == 0:
  if (tensorPosition == 0)
  {
   NormMPO.set(0, OneTensor);
   NormMPO.set(1, OneTensor);
   NormMPO.set(2, OneTensor);
   NormMPO.set(3, NormTensors[2]);
   NormMPS.get(N-2, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormMPO.set(4, Tensor0);
   NormMPS.get(N-1, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormMPO.set(5, Tensor0);
  }
// b)b) 0 < tensorPosition < N-2:
  else if ((tensorPosition > 0) && (tensorPosition < N-2))
  {
   NormMPO.set(0, NormTensors[tensorPosition-1]);
   NormMPO.set(1, OneTensor);
   NormMPO.set(2, OneTensor);
   NormMPO.set(3, NormTensors[tensorPosition+2]);
   NormMPS.get(N-1-tensorPosition-1, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormMPO.set(4, Tensor0);
   NormMPS.get(N-1-tensorPosition, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormMPO.set(5, Tensor0);
  }
// b)c) tensorPosition == N-2:
  else if (tensorPosition == N-2)
  {
   NormMPO.set(0, NormTensors[N-3]);
   NormMPO.set(1, OneTensor);
   NormMPO.set(2, OneTensor);
   NormMPO.set(3, OneTensor);
   NormMPS.get(0, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormMPO.set(4, Tensor0);
   NormMPS.get(1, Tensor0);
   Tensor0.getShape(Shape3);
   Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
   Tensor0.reshape(Shape4);
   NormMPO.set(5, Tensor0);
  }
 }
}

template<class T> void TimeEvolution2D<T>::getBulkNormMPO(const MPS<T>& NormMPS1, const MPS<T>& NormMPS2, const vector< Tensor<T> >& NormTensors,
                                                          unsigned int tensorPosition, MPO<T>& NormMPO) const
{
 unsigned int D = sqrt(NormMPS1.getd()), N = NormMPS1.getN();
 vector<unsigned int> Shape3(3), Shape4(4);
 Tensor<T> Tensor0;
// OneTensor:
 Shape4[0] = 1; Shape4[1] = 1; Shape4[2] = 1; Shape4[3] = 1;
 Tensor<T> OneTensor(Shape4);
 OneTensor.set(0, 1.0);
// NormMPO:
 NormMPO = MPO<T>("periodic", 6, D, max(NormMPS1.getD(), NormMPS2.getD()));
// a) tensorPosition == 0:
 if (tensorPosition == 0)
 {
  NormMPO.set(0, OneTensor);
  NormMPS2.get(0, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormMPO.set(1, Tensor0);
  NormMPS2.get(1, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormMPO.set(2, Tensor0);
  NormMPO.set(3, NormTensors[2]);
  NormMPS1.get(N-2, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormMPO.set(4, Tensor0);
  NormMPS1.get(N-1, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormMPO.set(5, Tensor0);
 }
// b) 0 < tensorPosition < N-2:
 else if ((tensorPosition > 0) && (tensorPosition < N-2))
 {
  NormMPO.set(0, NormTensors[tensorPosition-1]);
  NormMPS2.get(tensorPosition, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormMPO.set(1, Tensor0);
  NormMPS2.get(tensorPosition+1, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormMPO.set(2, Tensor0);
  NormMPO.set(3, NormTensors[tensorPosition+2]);
  NormMPS1.get(N-1-tensorPosition-1, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormMPO.set(4, Tensor0);
  NormMPS1.get(N-1-tensorPosition, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormMPO.set(5, Tensor0);
 }
// c) tensorPosition == N-2:
 else if (tensorPosition == N-2)
 {
  NormMPO.set(0, NormTensors[N-3]);
  NormMPS2.get(N-2, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormMPO.set(1, Tensor0);
  NormMPS2.get(N-1, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormMPO.set(2, Tensor0);
  NormMPO.set(3, OneTensor);
  NormMPS1.get(0, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = 1; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormMPO.set(4, Tensor0);
  NormMPS1.get(1, Tensor0);
  Tensor0.getShape(Shape3);
  Shape4[0] = Shape3[0]; Shape4[1] = Shape3[1]; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  NormMPO.set(5, Tensor0);
 }
}

template<class T> void TimeEvolution2D<T>::setRandomPositiveSeparableBoundaryMPS(MPS<T>& MPS0, const vector<unsigned int>& Seed)
{
 unsigned int D, N = MPS0.getN();
 vector<unsigned int> dMPS0, Indices0(1), Indices1(1), Shape3(3), Shape4(4);
 Tensor<T> Tensor0, Tensor1;
 MPS0.getd(dMPS0);
 MPS0 = MPS<T>("open", N, dMPS0, 1);
 Indices0[0] = 3; Indices1[0] = 3;
 for (int pos = 0; pos < N; pos++)
 {
  MPS0.get(pos, Tensor0);
  Tensor0.fillRandomly(Seed[pos]);
  Tensor0.getShape(Shape3);
  D = sqrt(Shape3[2]);
  Shape4[0] = 1; Shape4[1] = 1; Shape4[2] = D; Shape4[3] = D;
  Tensor0.reshape(Shape4);
  Tensor0.complexConjugate(Tensor1);
  Tensor0.contract(Indices0, Tensor1, Indices1);
  Tensor0.reshape(Shape3);
  MPS0.set(pos, Tensor0);
 }
}

template<class T> void TimeEvolution2D<T>::getInitialSeparableNormMPSs(const PEPS<T>& PEPS0, const string& UpdateDirection,
                                                                       unsigned int clusterSize, double eps, unsigned int maxNumSweeps,
                                                                       vector< MPS<T> >& SeparableNormMPSs) const
{
 if (((UpdateDirection == "vertical") && (clusterSize >= PEPS0.getNcols()-1)) || ((UpdateDirection == "horizontal") && (clusterSize >= PEPS0.getNrows()-1)))
  return;
 unsigned int D = PEPS0.getD(), Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), numSweepsDone;
 double epsAchieved;
 string BC = "open", Direction, WhichBoundary;
 vector<unsigned int> Seed;
 MPS<T> MPS0, MPS1;
// a) vertical updates:
 if (UpdateDirection == "vertical")
 {
  SeparableNormMPSs = vector< MPS<T> >(Ncols);
// a)1. get separable boundary-MPS:
  MPS0 = MPS<T>(BC, Nrows, D*D, 1);
  Seed = vector<unsigned int>(Nrows);
  for (int i = 0; i < Nrows; i++)
   Seed[i] = time(0)+13*i;
  TimeEvolution2D<T>::setRandomPositiveSeparableBoundaryMPS(MPS0, Seed);
  WhichBoundary = "right";
  PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, MPS0);
  SeparableNormMPSs[Ncols-1] = MPS0;
// a)2. get separable bulk-MPSs:
  Direction = "left";
  MPS1 = MPS0;
  for (int position = Ncols-2; position > clusterSize; position--)
  {
   MPS0 = MPS1;
   PEPS0.multiplyBulkMPOMPS(Direction, position, MPS0, eps, maxNumSweeps, epsAchieved, numSweepsDone, MPS1);
   SeparableNormMPSs[position] = MPS1;
  }
 }
// b) horizontal updates:
 else if (UpdateDirection == "horizontal")
 {
  SeparableNormMPSs = vector< MPS<T> >(Nrows);
// b)1. get separable boundary-MPS:
  MPS0 = MPS<T>(BC, Ncols, D*D, 1);
  Seed = vector<unsigned int>(Ncols);
  for (int i = 0; i < Ncols; i++)
   Seed[i] = time(0)+13*i;
  TimeEvolution2D<T>::setRandomPositiveSeparableBoundaryMPS(MPS0, Seed);
  WhichBoundary = "bottom";
  PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, MPS0);
  SeparableNormMPSs[Nrows-1] = MPS0;
// b)2. get separable bulk-MPSs:
  Direction = "up";
  MPS1 = MPS0;
  for (int position = Nrows-2; position > clusterSize; position--)
  {
   MPS0 = MPS1;
   PEPS0.multiplyBulkMPOMPS(Direction, position, MPS0, eps, maxNumSweeps, epsAchieved, numSweepsDone, MPS1);
   SeparableNormMPSs[position] = MPS1;
  }
 }
}

template<class T> void TimeEvolution2D<T>::updateSeparableNormMPSsBoundary(const PEPS<T>& PEPS0, const string& WhichBoundary,
                                                                           unsigned int clusterSize, double eps, unsigned int maxNumSweeps,
                                                                           vector< MPS<T> >& SeparableNormMPSs) const
{
 if (clusterSize > 0)
  return;
 unsigned int D = PEPS0.getD(), Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), numSweepsDone;
 double epsAchieved;
 string BC = "open";
 vector<unsigned int> Seed;
// a) left boundary:
 if (WhichBoundary == "left")
 {
  SeparableNormMPSs[0] = MPS<T>(BC, Nrows, D*D, 1);
  Seed = vector<unsigned int>(Nrows);
  for (int i = 0; i < Nrows; i++)
   Seed[i] = time(0)+13*i;
  TimeEvolution2D<T>::setRandomPositiveSeparableBoundaryMPS(SeparableNormMPSs[0], Seed);
  PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, SeparableNormMPSs[0]);
 }
// b) top boundary:
 else if (WhichBoundary == "top")
 {
  SeparableNormMPSs[0] = MPS<T>(BC, Ncols, D*D, 1);
  Seed = vector<unsigned int>(Ncols);
  for (int i = 0; i < Ncols; i++)
   Seed[i] = time(0)+13*i;
  TimeEvolution2D<T>::setRandomPositiveSeparableBoundaryMPS(SeparableNormMPSs[0], Seed);
  PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, SeparableNormMPSs[0]);
 }
// c) right boundary:
 else if (WhichBoundary == "right")
 {
  SeparableNormMPSs[Ncols-1] = MPS<T>(BC, Nrows, D*D, 1);
  Seed = vector<unsigned int>(Nrows);
  for (int i = 0; i < Nrows; i++)
   Seed[i] = time(0)+13*i;
  TimeEvolution2D<T>::setRandomPositiveSeparableBoundaryMPS(SeparableNormMPSs[Ncols-1], Seed);
  PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, SeparableNormMPSs[Ncols-1]);
 }
// d) bottom boundary:
 else if (WhichBoundary == "bottom")
 {
  SeparableNormMPSs[Nrows-1] = MPS<T>(BC, Ncols, D*D, 1);
  Seed = vector<unsigned int>(Ncols);
  for (int i = 0; i < Ncols; i++)
   Seed[i] = time(0)+13*i;
  TimeEvolution2D<T>::setRandomPositiveSeparableBoundaryMPS(SeparableNormMPSs[Nrows-1], Seed);
  PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, SeparableNormMPSs[Nrows-1]);
 }
}

template<class T> void TimeEvolution2D<T>::updateSeparableNormMPSsBulk(const PEPS<T>& PEPS0, const string& SweepDirection, unsigned int updatePosition,
                                                                       unsigned int clusterSize, double eps, unsigned int maxNumSweeps,
                                                                       vector< MPS<T> >& SeparableNormMPSs) const
{
 unsigned int D = PEPS0.getD(), Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), numSweepsDone;
 double epsAchieved;
 string BC = "open", WhichBoundary;
 vector<unsigned int> Seed;
// a) sweep right:
 if (SweepDirection == "right")
 {
// a)a) get positive separable boundary-MPS for left boundary:
  if (int(updatePosition)-int(clusterSize) == 0)
  {
   SeparableNormMPSs[0] = MPS<T>(BC, Nrows, D*D, 1);
   Seed = vector<unsigned int>(Nrows);
   for (int i = 0; i < Nrows; i++)
    Seed[i] = time(0)+13*i;
   TimeEvolution2D<T>::setRandomPositiveSeparableBoundaryMPS(SeparableNormMPSs[0], Seed);
   WhichBoundary = "left";
   PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, SeparableNormMPSs[0]);
  }
// a)b) get positive separable boundary-MPS in bulk:
  else if (int(updatePosition)-int(clusterSize) > 0)
  {
   SeparableNormMPSs[updatePosition-clusterSize] = SeparableNormMPSs[updatePosition-clusterSize-1];
   PEPS0.multiplyBulkMPOMPS(SweepDirection, updatePosition-clusterSize, SeparableNormMPSs[updatePosition-clusterSize-1],
                            eps, maxNumSweeps, epsAchieved, numSweepsDone, SeparableNormMPSs[updatePosition-clusterSize]);
  }
 }
// b) sweep down:
 else if (SweepDirection == "down")
 {
// b)a) get positive separable boundary-MPS for top boundary:
  if (int(updatePosition)-int(clusterSize) == 0)
  {
   SeparableNormMPSs[0] = MPS<T>(BC, Ncols, D*D, 1);
   Seed = vector<unsigned int>(Ncols);
   for (int i = 0; i < Ncols; i++)
    Seed[i] = time(0)+13*i;
   TimeEvolution2D<T>::setRandomPositiveSeparableBoundaryMPS(SeparableNormMPSs[0], Seed);
   WhichBoundary = "top";
   PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, SeparableNormMPSs[0]);
  }
// b)b) get positive separable boundary-MPS in bulk:
  else if (int(updatePosition)-int(clusterSize) > 0)
  {
   SeparableNormMPSs[updatePosition-clusterSize] = SeparableNormMPSs[updatePosition-clusterSize-1];
   PEPS0.multiplyBulkMPOMPS(SweepDirection, updatePosition-clusterSize, SeparableNormMPSs[updatePosition-clusterSize-1],
                            eps, maxNumSweeps, epsAchieved, numSweepsDone, SeparableNormMPSs[updatePosition-clusterSize]);
  }
 }
// c) sweep left:
 else if (SweepDirection == "left")
 {
// c)a) get positive separable boundary-MPS for right boundary:
  if (updatePosition+clusterSize == Ncols-1)
  {
   SeparableNormMPSs[Ncols-1] = MPS<T>(BC, Nrows, D*D, 1);
   Seed = vector<unsigned int>(Nrows);
   for (int i = 0; i < Nrows; i++)
    Seed[i] = time(0)+13*i;
   TimeEvolution2D<T>::setRandomPositiveSeparableBoundaryMPS(SeparableNormMPSs[Ncols-1], Seed);
   WhichBoundary = "right";
   PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, SeparableNormMPSs[Ncols-1]);
  }
// c)b) get positive separable boundary-MPS in bulk:
  else if (updatePosition+clusterSize < Ncols-1)
  {
   SeparableNormMPSs[updatePosition+clusterSize] = SeparableNormMPSs[updatePosition+clusterSize+1];
   PEPS0.multiplyBulkMPOMPS(SweepDirection, updatePosition+clusterSize, SeparableNormMPSs[updatePosition+clusterSize+1],
                            eps, maxNumSweeps, epsAchieved, numSweepsDone, SeparableNormMPSs[updatePosition+clusterSize]);
  }
 }
// d) sweep up:
 else if (SweepDirection == "up")
 {
// d)a) get positive separable boundary-MPS for bottom boundary:
  if (updatePosition+clusterSize == Nrows-1)
  {
   SeparableNormMPSs[Nrows-1] = MPS<T>(BC, Ncols, D*D, 1);
   Seed = vector<unsigned int>(Ncols);
   for (int i = 0; i < Ncols; i++)
    Seed[i] = time(0)+13*i;
   TimeEvolution2D<T>::setRandomPositiveSeparableBoundaryMPS(SeparableNormMPSs[Nrows-1], Seed);
   WhichBoundary = "bottom";
   PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, SeparableNormMPSs[Nrows-1]);
  }
// d)b) get positive separable boundary-MPS in bulk:
  else if (updatePosition+clusterSize < Nrows-1)
  {
   SeparableNormMPSs[updatePosition+clusterSize] = SeparableNormMPSs[updatePosition+clusterSize+1];
   PEPS0.multiplyBulkMPOMPS(SweepDirection, updatePosition+clusterSize, SeparableNormMPSs[updatePosition+clusterSize+1],
                            eps, maxNumSweeps, epsAchieved, numSweepsDone, SeparableNormMPSs[updatePosition+clusterSize]);
  }
 }
}

template<class T> void TimeEvolution2D<T>::getBoundaryNormMPS(const PEPS<T>& PEPS0, const string& WhichBoundary,
                                                              const vector<unsigned int>& D2s, double eps, unsigned int maxNumSweeps,
                                                              const vector< MPS<T> >& SeparableNormMPSs,
                                                              MPS<T>& NormMPS) const
{
 unsigned int clusterSize = D2s.size(), D = PEPS0.getD(), D2, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), numSweepsDone, positionStart;
 double epsAchieved;
 string BC = "open", Direction, WhichBoundary0;
 vector<unsigned int> Seed;
 MPS<T> MPS0;
// a) NormMPS for left boundary column:
 if (WhichBoundary == "left")
 {
// a)1. get NormMPS from the farthest column:
// a)1. - clusterSize < Ncols-1: get last separable NormMPS right of the cluster:
  if (clusterSize < Ncols-1)
  {
   NormMPS = SeparableNormMPSs[clusterSize+1];   
  }
// a)1. - clusterSize >= Ncols-1: get right boundary:
  else
  {
   D2 = min(D2s[Ncols-2], D*D);
   NormMPS = MPS<T>(BC, Nrows, D*D, D2);
   Seed = vector<unsigned int>(Nrows);
   for (int i = 0; i < Nrows; i++)
    Seed[i] = time(0)+13*i;
   NormMPS.fillRandomly(Seed);
   WhichBoundary0 = "right";
   PEPS0.getBoundaryMPS(WhichBoundary0, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPS);
  }
// a)2. contract the cluster:
  Direction = "left";
  positionStart = min(clusterSize, Ncols-2);
  for (int position = positionStart; position > 0; position--)
  {
   MPS0 = NormMPS;
   D2 = min(D2s[position-1], D*D*MPS0.getD());
   NormMPS.setD(D2);
   PEPS0.multiplyBulkMPOMPS(Direction, position, MPS0, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPS);
  }
 }
// b) NormMPS for top boundary row:
 else if (WhichBoundary == "top")
 {
// b)1. get NormMPS from the farthest row:
// b)1. - clusterSize < Nrows-1: get last separable NormMPS below the cluster:
  if (clusterSize < Nrows-1)
  {
   NormMPS = SeparableNormMPSs[clusterSize+1];   
  }
// b)1. - clusterSize >= Nrows-1: get bottom boundary:
  else
  {
   D2 = min(D2s[Nrows-2], D*D);
   NormMPS = MPS<T>(BC, Ncols, D*D, D2);
   Seed = vector<unsigned int>(Ncols);
   for (int i = 0; i < Ncols; i++)
    Seed[i] = time(0)+13*i;
   NormMPS.fillRandomly(Seed);
   WhichBoundary0 = "bottom";
   PEPS0.getBoundaryMPS(WhichBoundary0, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPS);
  }
// b)2. contract the cluster:
  Direction = "up";
  positionStart = min(clusterSize, Nrows-2);
  for (int position = positionStart; position > 0; position--)
  {
   MPS0 = NormMPS;
   D2 = min(D2s[position-1], D*D*MPS0.getD());
   NormMPS.setD(D2);
   PEPS0.multiplyBulkMPOMPS(Direction, position, MPS0, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPS);
  }
 }
// c) NormMPS for right boundary column:
 else if (WhichBoundary == "right")
 {
// c)1. get NormMPS from the farthest column:
// c)1. - Ncols-1-clusterSize > 0: get last separable NormMPS left of the cluster:
  if (int(Ncols)-1-int(clusterSize) > 0)
  {
   NormMPS = SeparableNormMPSs[Ncols-1-clusterSize-1];   
  }
// c)1. - Ncols-1-clusterSize <= 0: get left boundary:
  else
  {
   D2 = min(D2s[Ncols-2], D*D);
   NormMPS = MPS<T>(BC, Nrows, D*D, D2);
   Seed = vector<unsigned int>(Nrows);
   for (int i = 0; i < Nrows; i++)
    Seed[i] = time(0)+13*i;
   NormMPS.fillRandomly(Seed);
   WhichBoundary0 = "left";
   PEPS0.getBoundaryMPS(WhichBoundary0, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPS);
  }
// c)2. contract the cluster:
  Direction = "right";
  positionStart = max(int(Ncols)-1-int(clusterSize), 1);
  for (int position = positionStart; position < Ncols-1; position++)
  {
   MPS0 = NormMPS;
   D2 = min(D2s[Ncols-1-position-1], D*D*MPS0.getD());
   NormMPS.setD(D2);
   PEPS0.multiplyBulkMPOMPS(Direction, position, MPS0, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPS);
  }
 }
// d) NormMPS for bottom boundary row:
 else if (WhichBoundary == "bottom")
 {
// d)1. get NormMPS from the farthest row:
// d)1. - Nrows-1-clusterSize > 0: get last separable NormMPS above the cluster:
  if (int(Nrows)-1-int(clusterSize) > 0)
  {
   NormMPS = SeparableNormMPSs[Nrows-1-clusterSize-1];   
  }
// d)1. - Nrows-1-clusterSize <= 0: get top boundary:
  else
  {
   D2 = min(D2s[Nrows-2], D*D);
   NormMPS = MPS<T>(BC, Ncols, D*D, D2);
   Seed = vector<unsigned int>(Ncols);
   for (int i = 0; i < Ncols; i++)
    Seed[i] = time(0)+13*i;
   NormMPS.fillRandomly(Seed);
   WhichBoundary0 = "top";
   PEPS0.getBoundaryMPS(WhichBoundary0, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPS);
  }
// d)2. contract the cluster:
  Direction = "down";
  positionStart = max(int(Nrows)-1-int(clusterSize), 1);
  for (int position = positionStart; position < Nrows-1; position++)
  {
   MPS0 = NormMPS;
   D2 = min(D2s[Nrows-1-position-1], D*D*MPS0.getD());
   NormMPS.setD(D2);
   PEPS0.multiplyBulkMPOMPS(Direction, position, MPS0, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPS);
  }
 }
}

template<class T> void TimeEvolution2D<T>::getBulkNormMPSs(const PEPS<T>& PEPS0, const string& UpdateDirection, unsigned int updatePosition,
                                                           const vector<unsigned int>& D2s, double eps, unsigned int maxNumSweeps,
                                                           const vector< MPS<T> >& SeparableNormMPSs,
                                                           MPS<T>& NormMPS1, MPS<T>& NormMPS2) const
{
 unsigned int clusterSize = D2s.size(), D = PEPS0.getD(), D2, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), numSweepsDone, positionStart;
 double epsAchieved;
 string BC = "open", Direction, WhichBoundary;
 vector<unsigned int> Seed;
 MPS<T> MPS0;
// a) UpdateDirection == "vertical":
 if (UpdateDirection == "vertical")
 {
// a)a) left NormMPS for column updatePosition:
// a)a)1. get NormMPS from the farthest column:
// a)a)1. - updatePosition-clusterSize > 0: get last separable NormMPS left of the cluster:
  if (int(updatePosition)-int(clusterSize) > 0)
  {
   NormMPS1 = SeparableNormMPSs[updatePosition-clusterSize-1];   
  }
// a)a)1. - updatePosition-clusterSize <= 0: get left boundary:
  else
  {
   D2 = min(D2s[updatePosition-1], D*D);
   NormMPS1 = MPS<T>(BC, Nrows, D*D, D2);
   Seed = vector<unsigned int>(Nrows);
   for (int i = 0; i < Nrows; i++)
    Seed[i] = time(0)+13*i;
   NormMPS1.fillRandomly(Seed);
   WhichBoundary = "left";
   PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPS1);
  }
// a)a)2. contract the cluster:
  Direction = "right";
  positionStart = max(int(updatePosition)-int(clusterSize), 1);
  for (int position = positionStart; position < updatePosition; position++)
  {
   MPS0 = NormMPS1;
   D2 = min(D2s[updatePosition-position-1], D*D*MPS0.getD());
   NormMPS1.setD(D2);
   PEPS0.multiplyBulkMPOMPS(Direction, position, MPS0, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPS1);
  }
// a)b) right NormMPS for column updatePosition:
// a)b)1. get NormMPS from the farthest column:
// a)b)1. - updatePosition+clusterSize < Ncols-1: get last separable NormMPS right of the cluster:
  if (updatePosition+clusterSize < Ncols-1)
  {
   NormMPS2 = SeparableNormMPSs[updatePosition+clusterSize+1];   
  }
// a)b)1. - updatePosition+clusterSize >= Ncols-1: get right boundary:
  else
  {
   D2 = min(D2s[Ncols-1-updatePosition-1], D*D);
   NormMPS2 = MPS<T>(BC, Nrows, D*D, D2);
   Seed = vector<unsigned int>(Nrows);
   for (int i = 0; i < Nrows; i++)
    Seed[i] = time(0)+13*i;
   NormMPS2.fillRandomly(Seed);
   WhichBoundary = "right";
   PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPS2);
  }
// a)b)2. contract the cluster:
  Direction = "left";
  positionStart = min(updatePosition+clusterSize, Ncols-2);
  for (int position = positionStart; position > updatePosition; position--)
  {
   MPS0 = NormMPS2;
   D2 = min(D2s[position-updatePosition-1], D*D*MPS0.getD());
   NormMPS2.setD(D2);
   PEPS0.multiplyBulkMPOMPS(Direction, position, MPS0, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPS2);
  }
 }
// b) UpdateDirection == "horizontal":
 else if (UpdateDirection == "horizontal")
 {
// b)a) lower NormMPS for row updatePosition:
// b)a)1. get NormMPS from the farthest row:
// b)a)1. - updatePosition+clusterSize < Nrows-1: get last separable NormMPS below the cluster:
  if (updatePosition+clusterSize < Nrows-1)
  {
   NormMPS1 = SeparableNormMPSs[updatePosition+clusterSize+1];   
  }
// b)a)1. - updatePosition+clusterSize >= Nrows-1: get bottom boundary:
  else
  {
   D2 = min(D2s[Nrows-1-updatePosition-1], D*D);
   NormMPS1 = MPS<T>(BC, Ncols, D*D, D2);
   Seed = vector<unsigned int>(Ncols);
   for (int i = 0; i < Ncols; i++)
    Seed[i] = time(0)+13*i;
   NormMPS1.fillRandomly(Seed);
   WhichBoundary = "bottom";
   PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPS1);
  }
// b)a)2. contract the cluster:
  Direction = "up";
  positionStart = min(updatePosition+clusterSize, Nrows-2);
  for (int position = positionStart; position > updatePosition; position--)
  {
   MPS0 = NormMPS1;
   D2 = min(D2s[position-updatePosition-1], D*D*MPS0.getD());
   NormMPS1.setD(D2);
   PEPS0.multiplyBulkMPOMPS(Direction, position, MPS0, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPS1);
  }
// b)b) upper NormMPS for row updatePosition:
// b)b)1. get NormMPS from the farthest row:
// b)b)1. - updatePosition-clusterSize > 0: get last separable NormMPS above the cluster:
  if (int(updatePosition)-int(clusterSize) > 0)
  {
   NormMPS2 = SeparableNormMPSs[updatePosition-clusterSize-1];   
  }
// b)b)1. - updatePosition-clusterSize <= 0: get top boundary:
  else
  {
   D2 = min(D2s[updatePosition-1], D*D);
   NormMPS2 = MPS<T>(BC, Ncols, D*D, D2);
   Seed = vector<unsigned int>(Ncols);
   for (int i = 0; i < Ncols; i++)
    Seed[i] = time(0)+13*i;
   NormMPS2.fillRandomly(Seed);
   WhichBoundary = "top";
   PEPS0.getBoundaryMPS(WhichBoundary, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPS2);
  }
// b)b)2. contract the cluster:
  Direction = "down";
  positionStart = max(int(updatePosition)-int(clusterSize), 1);
  for (int position = positionStart; position < updatePosition; position++)
  {
   MPS0 = NormMPS2;
   D2 = min(D2s[updatePosition-position-1], D*D*MPS0.getD());
   NormMPS2.setD(D2);
   PEPS0.multiplyBulkMPOMPS(Direction, position, MPS0, eps, maxNumSweeps, epsAchieved, numSweepsDone, NormMPS2);
  }
 }
}

template<class T> void TimeEvolution2D<T>::updatePEPSBoundary(const string& WhichBoundary, unsigned int tensorPositionStart, unsigned int tensorPositionStop,
                                                              const vector< MPO<T> >& TEMPOs, const string& SweepUpdate,
                                                              const vector< MPS<T> >& NormMPSs, const string& bEnvironment, const vector< MPS<T> >& bMPSs,
                                                              const string& UpdateTensor, const string& UpdateMode, double epsLoc, unsigned int maxNumSweepsLoc, double cutoff,
                                                              PEPS<T>& PEPS0, PEPS<T>& PEPS1) const
{
 unsigned int col, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), row;
 string UpdateDirection;
 MPO<T> bMPO, NormMPO, TEMPO;
 vector< Tensor<T> > NormTensors;
// a) left boundary column:
 if (WhichBoundary == "left")
 {
  col = 0;
  UpdateDirection = "vertical";
// a)1. get initial boundary norm-tensors:
  this->getInitialBoundaryNormTensors(PEPS0, WhichBoundary, NormMPSs[1], NormTensors);
// a)2. update boundary column:
  for (row = tensorPositionStart; row < tensorPositionStop; row += 2)
  {
// a)2.1. get norm-MPO:
   this->getBoundaryNormMPO(WhichBoundary, NormMPSs[1], NormTensors, row, NormMPO);
// a)2.2. update tensor-pair:
   if (SweepUpdate == "independent")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPOs[0], UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS1);
   else if (SweepUpdate == "sequential")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPOs[0], UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0);
// a)2.3. update norm-tensors:
   this->updateBoundaryNormTensors(PEPS0, WhichBoundary, NormMPSs[1], row, NormTensors);
  }
 }
// b) top boundary row:
 else if (WhichBoundary == "top")
 {
  row = 0;
  UpdateDirection = "horizontal";
// b)1. get initial boundary norm-tensors:
  this->getInitialBoundaryNormTensors(PEPS0, WhichBoundary, NormMPSs[1], NormTensors);
// b)2. update boundary row:
  for (col = tensorPositionStart; col < tensorPositionStop; col += 2)
  {
// b)2.1. get norm-MPO:
   this->getBoundaryNormMPO(WhichBoundary, NormMPSs[1], NormTensors, col, NormMPO);
// b)2.2. update tensor-pair:
// b)2.2. - horizontal-even:
   if (tensorPositionStart == 0)
   {
    if (col == 0)
     TEMPO = TEMPOs[1];
    else if ((col > 0) && (col < Ncols-2))
     TEMPO = TEMPOs[0];
    else if (col == Ncols-2)
     TEMPO = TEMPOs[2];
   }
// b)2.2. - horizontal-odd:
   else if (tensorPositionStart == 1)
   {
    TEMPO = TEMPOs[3];
   }
   if (SweepUpdate == "independent")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPO, UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS1);
   else if (SweepUpdate == "sequential")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPO, UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0);
// b)2.3. update norm-tensors:
   this->updateBoundaryNormTensors(PEPS0, WhichBoundary, NormMPSs[1], col, NormTensors);
  }
 }
// c) right boundary column:
 else if (WhichBoundary == "right")
 {
  col = Ncols-1;
  UpdateDirection = "vertical";
// c)1. get initial boundary norm-tensors:
  this->getInitialBoundaryNormTensors(PEPS0, WhichBoundary, NormMPSs[Ncols-2], NormTensors);
// c)2. update boundary column:
  for (row = tensorPositionStart; row < tensorPositionStop; row += 2)
  {
// c)2.1. get norm-MPO:
   this->getBoundaryNormMPO(WhichBoundary, NormMPSs[Ncols-2], NormTensors, row, NormMPO);
// c)2.2. update tensor-pair:
   if (SweepUpdate == "independent")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPOs[0], UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS1);
   else if (SweepUpdate == "sequential")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPOs[0], UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0);
// c)2.3. update norm-tensors:
   this->updateBoundaryNormTensors(PEPS0, WhichBoundary, NormMPSs[Ncols-2], row, NormTensors);
  }
 }
// d) bottom boundary row:
 else if (WhichBoundary == "bottom")
 {
  row = Nrows-1;
  UpdateDirection = "horizontal";
// d)1. get initial boundary norm-tensors:
  this->getInitialBoundaryNormTensors(PEPS0, WhichBoundary, NormMPSs[Nrows-2], NormTensors);
// d)2. update boundary row:
  for (col = tensorPositionStart; col < tensorPositionStop; col += 2)
  {
// d)2.1. get norm-MPO:
   this->getBoundaryNormMPO(WhichBoundary, NormMPSs[Nrows-2], NormTensors, col, NormMPO);
// d)2.2. update tensor-pair:
// d)2.2. - horizontal-even:
   if (tensorPositionStart == 0)
   {
    if (col == 0)
     TEMPO = TEMPOs[1];
    else if ((col > 0) && (col < Ncols-2))
     TEMPO = TEMPOs[0];
    else if (col == Ncols-2)
     TEMPO = TEMPOs[2];
   }
// d)2.2. - horizontal-odd:
   else if (tensorPositionStart == 1)
   {
    TEMPO = TEMPOs[3];
   }
   if (SweepUpdate == "independent")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPO, UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS1);
   else if (SweepUpdate == "sequential")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPO, UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0);
// d)2.3. update norm-tensors:
   this->updateBoundaryNormTensors(PEPS0, WhichBoundary, NormMPSs[Nrows-2], col, NormTensors);
  }
 }
}

template<class T> void TimeEvolution2D<T>::updatePEPSBulk(const string& UpdateDirection, unsigned int updatePosition, unsigned int tensorPositionStart, unsigned int tensorPositionStop,
                                                          const vector< MPO<T> >& TEMPOs, const string& SweepUpdate,
                                                          const vector< MPS<T> >& NormMPSs, const string& bEnvironment, const vector< MPS<T> >& bMPSs,
                                                          const string& UpdateTensor, const string& UpdateMode, double epsLoc, unsigned int maxNumSweepsLoc, double cutoff,
                                                          PEPS<T>& PEPS0, PEPS<T>& PEPS1) const
{
 unsigned int col, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), row;
 MPO<T> bMPO, NormMPO, TEMPO;
 vector< Tensor<T> > NormTensors;
// a) vertical updates:
 if (UpdateDirection == "vertical")
 {
  col = updatePosition;
// a)1. get initial norm-tensors:
  this->getInitialBulkNormTensors(PEPS0, UpdateDirection, updatePosition, NormMPSs[updatePosition-1], NormMPSs[updatePosition+1], NormTensors);
// a)2. update bulk column:
  for (row = tensorPositionStart; row < tensorPositionStop; row += 2)
  {
// a)2.1. get norm-MPO:
   this->getBulkNormMPO(NormMPSs[updatePosition-1], NormMPSs[updatePosition+1], NormTensors, row, NormMPO);
// a)2.2. update tensor-pair:
   if (SweepUpdate == "independent")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPOs[0], UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS1);
   else if (SweepUpdate == "sequential")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPOs[0], UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0);
// a)2.3. update norm-tensors:
   this->updateBulkNormTensors(PEPS0, UpdateDirection, updatePosition, NormMPSs[updatePosition-1], NormMPSs[updatePosition+1], row, NormTensors);
  }
 }
// b) horizontal updates:
 else if (UpdateDirection == "horizontal")
 {
  row = updatePosition;
// b)1. get initial norm-tensors:
  this->getInitialBulkNormTensors(PEPS0, UpdateDirection, updatePosition, NormMPSs[updatePosition+1], NormMPSs[updatePosition-1], NormTensors);
// b)2. update bulk row:
  for (col = tensorPositionStart; col < tensorPositionStop; col += 2)
  {
// b)2.1. get norm-MPO:
   this->getBulkNormMPO(NormMPSs[updatePosition+1], NormMPSs[updatePosition-1], NormTensors, col, NormMPO);
// b)2.2. update tensor-pair:
// b)2.2. - horizontal-even:
   if (tensorPositionStart == 0)
   {
    if (col == 0)
     TEMPO = TEMPOs[1];
    else if ((col > 0) && (col < Ncols-2))
     TEMPO = TEMPOs[0];
    else if (col == Ncols-2)
     TEMPO = TEMPOs[2];
   }
// b)2.2. - horizontal-odd:
   else if (tensorPositionStart == 1)
   {
    TEMPO = TEMPOs[3];
   }
   if (SweepUpdate == "independent")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPO, UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS1);
   else if (SweepUpdate == "sequential")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPO, UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0);
// b)2.3. update norm-tensors:
   this->updateBulkNormTensors(PEPS0, UpdateDirection, updatePosition, NormMPSs[updatePosition+1], NormMPSs[updatePosition-1], col, NormTensors);
  }
 }
}

template<class T> void TimeEvolution2D<T>::clusterUpdatePEPSBoundary(const string& WhichBoundary, unsigned int tensorPositionStart, unsigned int tensorPositionStop,
                                                                     const vector< MPO<T> >& TEMPOs, const string& SweepUpdate,
                                                                     const vector<unsigned int>& D2sEnv, double epsEnv, unsigned int maxNumSweepsEnv,
                                                                     const vector< MPS<T> >& SeparableNormMPSs, const string& bEnvironment, const vector< MPS<T> >& SeparablebMPSs,
                                                                     const string& UpdateTensor, const string& UpdateMode, double epsLoc, unsigned int maxNumSweepsLoc, double cutoff,
                                                                     PEPS<T>& PEPS0, PEPS<T>& PEPS1) const
{
 unsigned int col, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), row;
 string UpdateDirection;
 MPS<T> NormMPS;
 MPO<T> bMPO, NormMPO, TEMPO;
 vector< Tensor<T> > NormTensors;
// get boundary-MPS for the boundary:
 this->getBoundaryNormMPS(PEPS0, WhichBoundary, D2sEnv, epsEnv, maxNumSweepsEnv, SeparableNormMPSs, NormMPS);
// a) left boundary column:
 if (WhichBoundary == "left")
 {
  col = 0;
  UpdateDirection = "vertical";
// a)1. get initial boundary norm-tensors:
  this->getInitialBoundaryNormTensors(PEPS0, WhichBoundary, NormMPS, NormTensors);
// a)2. update boundary column:
  for (row = tensorPositionStart; row < tensorPositionStop; row += 2)
  {
// a)2.1. get norm-MPO:
   this->getBoundaryNormMPO(WhichBoundary, NormMPS, NormTensors, row, NormMPO);
// a)2.2. update tensor-pair:
   if (SweepUpdate == "independent")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPOs[0], UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS1);
   else if (SweepUpdate == "sequential")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPOs[0], UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0);
// a)2.3. update norm-tensors:
   this->updateBoundaryNormTensors(PEPS0, WhichBoundary, NormMPS, row, NormTensors);
  }
 }
// b) top boundary row:
 else if (WhichBoundary == "top")
 {
  row = 0;
  UpdateDirection = "horizontal";
// b)1. get initial boundary norm-tensors:
  this->getInitialBoundaryNormTensors(PEPS0, WhichBoundary, NormMPS, NormTensors);
// b)2. update boundary row:
  for (col = tensorPositionStart; col < tensorPositionStop; col += 2)
  {
// b)2.1. get norm-MPO:
   this->getBoundaryNormMPO(WhichBoundary, NormMPS, NormTensors, col, NormMPO);
// b)2.2. update tensor-pair:
// b)2.2. - horizontal-even:
   if (tensorPositionStart == 0)
   {
    if (col == 0)
     TEMPO = TEMPOs[1];
    else if ((col > 0) && (col < Ncols-2))
     TEMPO = TEMPOs[0];
    else if (col == Ncols-2)
     TEMPO = TEMPOs[2];
   }
// b)2.2. - horizontal-odd:
   else if (tensorPositionStart == 1)
   {
    TEMPO = TEMPOs[3];
   }
   if (SweepUpdate == "independent")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPO, UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS1);
   else if (SweepUpdate == "sequential")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPO, UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0);
// b)2.3. update norm-tensors:
   this->updateBoundaryNormTensors(PEPS0, WhichBoundary, NormMPS, col, NormTensors);
  }
 }
// c) right boundary column:
 else if (WhichBoundary == "right")
 {
  col = Ncols-1;
  UpdateDirection = "vertical";
// c)1. get initial boundary norm-tensors:
  this->getInitialBoundaryNormTensors(PEPS0, WhichBoundary, NormMPS, NormTensors);
// c)2. update boundary column:
  for (row = tensorPositionStart; row < tensorPositionStop; row += 2)
  {
// c)2.1. get norm-MPO:
   this->getBoundaryNormMPO(WhichBoundary, NormMPS, NormTensors, row, NormMPO);
// c)2.2. update tensor-pair:
   if (SweepUpdate == "independent")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPOs[0], UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS1);
   else if (SweepUpdate == "sequential")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPOs[0], UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0);
// c)2.3. update norm-tensors:
   this->updateBoundaryNormTensors(PEPS0, WhichBoundary, NormMPS, row, NormTensors);
  }
 }
// d) bottom boundary row:
 else if (WhichBoundary == "bottom")
 {
  row = Nrows-1;
  UpdateDirection = "horizontal";
// d)1. get initial boundary norm-tensors:
  this->getInitialBoundaryNormTensors(PEPS0, WhichBoundary, NormMPS, NormTensors);
// d)2. update boundary row:
  for (col = tensorPositionStart; col < tensorPositionStop; col += 2)
  {
// d)2.1. get norm-MPO:
   this->getBoundaryNormMPO(WhichBoundary, NormMPS, NormTensors, col, NormMPO);
// d)2.2. update tensor-pair:
// d)2.2. - horizontal-even:
   if (tensorPositionStart == 0)
   {
    if (col == 0)
     TEMPO = TEMPOs[1];
    else if ((col > 0) && (col < Ncols-2))
     TEMPO = TEMPOs[0];
    else if (col == Ncols-2)
     TEMPO = TEMPOs[2];
   }
// d)2.2. - horizontal-odd:
   else if (tensorPositionStart == 1)
   {
    TEMPO = TEMPOs[3];
   }
   if (SweepUpdate == "independent")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPO, UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS1);
   else if (SweepUpdate == "sequential")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPO, UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0);
// d)2.3. update norm-tensors:
   this->updateBoundaryNormTensors(PEPS0, WhichBoundary, NormMPS, col, NormTensors);
  }
 }
}

template<class T> void TimeEvolution2D<T>::clusterUpdatePEPSBulk(const string& UpdateDirection, unsigned int updatePosition, unsigned int tensorPositionStart, unsigned int tensorPositionStop,
                                                                 const vector< MPO<T> >& TEMPOs, const string& SweepUpdate,
                                                                 const vector<unsigned int>& D2sEnv, double epsEnv, unsigned int maxNumSweepsEnv,
                                                                 const vector< MPS<T> >& SeparableNormMPSs, const string& bEnvironment, const vector< MPS<T> >& SeparablebMPSs,
                                                                 const string& UpdateTensor, const string& UpdateMode, double epsLoc, unsigned int maxNumSweepsLoc, double cutoff,
                                                                 PEPS<T>& PEPS0, PEPS<T>& PEPS1) const
{
 unsigned int col, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), row;
 MPS<T> NormMPS1, NormMPS2;
 MPO<T> bMPO, NormMPO, TEMPO;
 vector< Tensor<T> > NormTensors;
// get boundary-MPSs for the bulk:
 this->getBulkNormMPSs(PEPS0, UpdateDirection, updatePosition, D2sEnv, epsEnv, maxNumSweepsEnv, SeparableNormMPSs, NormMPS1, NormMPS2);
// a) vertical updates:
 if (UpdateDirection == "vertical")
 {
  col = updatePosition;
// a)1. get initial norm-tensors:
  this->getInitialBulkNormTensors(PEPS0, UpdateDirection, updatePosition, NormMPS1, NormMPS2, NormTensors);
// a)2. update bulk column:
  for (row = tensorPositionStart; row < tensorPositionStop; row += 2)
  {
// a)2.1. get norm-MPO:
   this->getBulkNormMPO(NormMPS1, NormMPS2, NormTensors, row, NormMPO);
// a)2.2. update tensor-pair:
   if (SweepUpdate == "independent")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPOs[0], UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS1);
   else if (SweepUpdate == "sequential")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPOs[0], UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0);
// a)2.3. update norm-tensors:
   this->updateBulkNormTensors(PEPS0, UpdateDirection, updatePosition, NormMPS1, NormMPS2, row, NormTensors);
  }
 }
// b) horizontal updates:
 else if (UpdateDirection == "horizontal")
 {
  row = updatePosition;
// b)1. get initial norm-tensors:
  this->getInitialBulkNormTensors(PEPS0, UpdateDirection, updatePosition, NormMPS1, NormMPS2, NormTensors);
// b)2. update bulk row:
  for (col = tensorPositionStart; col < tensorPositionStop; col += 2)
  {
// b)2.1. get norm-MPO:
   this->getBulkNormMPO(NormMPS1, NormMPS2, NormTensors, col, NormMPO);
// b)2.2. update tensor-pair:
// b)2.2. - horizontal-even:
   if (tensorPositionStart == 0)
   {
    if (col == 0)
     TEMPO = TEMPOs[1];
    else if ((col > 0) && (col < Ncols-2))
     TEMPO = TEMPOs[0];
    else if (col == Ncols-2)
     TEMPO = TEMPOs[2];
   }
// b)2.2. - horizontal-odd:
   else if (tensorPositionStart == 1)
   {
    TEMPO = TEMPOs[3];
   }
   if (SweepUpdate == "independent")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPO, UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS1);
   else if (SweepUpdate == "sequential")
    this->updateTensor(row, col, UpdateDirection, NormMPO, bEnvironment, bMPO, TEMPO, UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0);
// b)2.3. update norm-tensors:
   this->updateBulkNormTensors(PEPS0, UpdateDirection, updatePosition, NormMPS1, NormMPS2, col, NormTensors);
  }
 }
}

template<class T> void TimeEvolution2D<T>::updateTensor(unsigned int positionRow, unsigned int positionCol, const string& UpdateDirection,
                                                        MPO<T>& NormMPO, const string& bEnvironment, const MPO<T>& bMPO, const MPO<T>& TEMPO,
                                                        const string& UpdateTensor, const string& UpdateMode, double eps, unsigned int maxNumSweeps, double cutoff,
                                                        PEPS<T>& PEPS0) const
{
 unsigned int D = PEPS0.getD(), dim0, dim0Red, dim1, dim1Red, dimNormTensor, dimPos, index0, index1, sweep;
 double error0, error1, norm;
 T b;
 vector<unsigned int> Indices1(1), Indices2(2), Indices3(3), Indices4(4), Order3(3), Order4(4), Order5(5), Shape2(2), Shape3(3), Shape4(4), Shape5(5);
 vector<unsigned int> Indices01(1), Indices11(1), Indices02(2), Indices12(2), Indices03(3), Indices13(3), Indices04(4), Indices14(4), Indices08(8), Indices18(8);
 vector<unsigned int> Index3(3), Index4(4), Index5(5), Index8(8), Shape0(5), Shape0Red(3), Shape1(5), Shape1Red(3);
 vector<double> W;
 vector<T> bVector, Vector0;
 Tensor<T> bTensor, bTensor0, NormTensor, NormTensor0, SqrtNormTensor, Tensor0Red, Tensor1Red, TensorL, TensorR, TensorPair;
 Tensor<T> EnvTensor, GaugeTensor, InvGaugeTensor, Tensor0, Tensor1, Tensor2, Tensor3, TETensor0, TETensor1;
 Matrix<T> EnvMatrix, EnvMatrixAdj, Lambda, Matrix0, Matrix0Inv, NormTensorMatrix, NormTensorMatrixAdj, Vr;
 Matrix<T> NormMatrix, NormMatrix0, NormMatrixAdj, NormMatrixInv;
 vector< Tensor<T> > GaugeTensors(6), InvGaugeTensors(6), InvGaugeTensors2(2);
 if (bEnvironment == "full")
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void TimeEvolution2D<T>::" <<
          "updateTensor(unsigned int positionRow, unsigned int positionCol, const string& UpdateDirection, " <<
                       "MPO<T>& NormMPO, const string& bEnvironment, const MPO<T>& bMPO, const MPO<T>& TEMPO, " <<
                       "const string& UpdateTensor, const string& UpdateMode, double eps, unsigned int maxNumSweeps, " <<
                       "double cutoff, PEPS<T>& PEPS0) const: " <<
          "(bEnvironment == full) is not implemented." << endl;
  exit(1);
 }
// 0. reshape TETensors:
 TEMPO.get(0, TETensor0); TEMPO.get(1, TETensor1);
 TETensor0.getShape(Shape4);
 Shape3[0] = Shape4[1]; Shape3[1] = Shape4[2]; Shape3[2] = Shape4[3];
 TETensor0.reshape(Shape3);
 TETensor1.getShape(Shape4);
 Shape3[0] = Shape4[0]; Shape3[1] = Shape4[2]; Shape3[2] = Shape4[3];
 TETensor1.reshape(Shape3);
// 1. get tensor pair from PEPS0:
 PEPS0.get(positionRow, positionCol, Tensor0);
 if (UpdateDirection == "vertical")
 {
  PEPS0.get(positionRow+1, positionCol, Tensor1);
// 1. - permute tensor pair to conform it to the order of NormMPO:
  Order5[0] = 1; Order5[1] = 2; Order5[2] = 3; Order5[3] = 0; Order5[4] = 4;
  Tensor0.permute(Order5);
  Tensor1.permute(Order5);
 }
 else if (UpdateDirection == "horizontal")
 {
  PEPS0.get(positionRow, positionCol+1, Tensor1);
 }
// 2. update tensor pair:
// 2.a) separable environment:
 if (NormMPO.getD() == 1)
 {
// 2.a)1. gauge tensors:
  for (int pos = 0; pos < 6; pos++)
  {
   NormMPO.get(pos, EnvTensor);
   EnvTensor.getShape(Shape4);
   if (Shape4[2] == D)
   {
    EnvMatrix = Matrix<T>(D, D);
    for (int i = 0; i < D*D; i++)
     EnvMatrix.set(i, EnvTensor.get(i));
    EnvMatrix.adjoint(EnvMatrixAdj);
    EnvMatrix.add(EnvMatrixAdj);
    EnvMatrix.multiply(0.5);
    EnvMatrix.setType("hermitian");
    W = vector<double>(D); Vr = Matrix<T>(D, D);
    EnvMatrix.eigenDecompose(W, Vr);
    dimPos = 0;
    if (abs(W[D-1]) >= abs(W[0]))
    {
     for (int i = D-1; i >= 0; i--){
      if (W[i] > 0.0)
       dimPos++;
      else
       break;
     }
     Shape2[0] = dimPos; Shape2[1] = D;
     GaugeTensor = Tensor<T>(Shape2);
     for (int j = 0; j < Shape2[1]; j++){
      for (int i = 0; i < Shape2[0]; i++){
       GaugeTensor.set(i+j*Shape2[0], Vr(j, D-1-i)*sqrt(W[D-1-i]));
     }}
     Shape2[0] = D; Shape2[1] = dimPos;
     InvGaugeTensors[pos] = Tensor<T>(Shape2);
     for (int j = 0; j < Shape2[1]; j++){
      for (int i = 0; i < Shape2[0]; i++){
       InvGaugeTensors[pos].set(i+j*Shape2[0], MathAuxiliary::complexConjugate(Vr(i, D-1-j))/sqrt(W[D-1-j]));
     }}
    }
    else
    {
     if (pos == 5)
     {
      cerr << "Program terminated because of error in function " <<
              "template<class T> void TimeEvolution2D<T>::" <<
              "updateTensor(unsigned int positionRow, unsigned int positionCol, const string& UpdateDirection, " <<
                           "MPO<T>& NormMPO, const string& bEnvironment, const MPO<T>& bMPO, const MPO<T>& TEMPO, " <<
                           "const string& UpdateTensor, const string& UpdateMode, double eps, unsigned int maxNumSweeps, " <<
                           "double cutoff, PEPS<T>& PEPS0) const: " <<
              "All eigenvalues of the norm-matrix were negative. Increase the bond dimension of the environment!" << endl;
      exit(1);
     }
     for (int pos2 = pos+1; pos2 < 6; pos2++)
     {
      NormMPO.get(pos2, Tensor2);
      if (Tensor2.getSize() != 1)
      {
       Tensor2.multiply(-1.0);
       NormMPO.set(pos2, Tensor2);
       break;
      }
      if ((pos2 == 5) && (Tensor2.getSize() == 1))
      {
       cerr << "Program terminated because of error in function " <<
               "template<class T> void TimeEvolution2D<T>::" <<
               "updateTensor(unsigned int positionRow, unsigned int positionCol, const string& UpdateDirection, " <<
                            "MPO<T>& NormMPO, const string& bEnvironment, const MPO<T>& bMPO, const MPO<T>& TEMPO, " <<
                            "const string& UpdateTensor, const string& UpdateMode, double eps, unsigned int maxNumSweeps, " <<
                            "double cutoff, PEPS<T>& PEPS0) const: " <<
               "All eigenvalues of the norm-matrix were negative. Increase the bond dimension of the environment!" << endl;
       exit(1);
      }
     }
     for (int i = 0; i < D; i++){
      if (W[i] < 0.0)
       dimPos++;
      else
       break;
     }
     Shape2[0] = dimPos; Shape2[1] = D;
     GaugeTensor = Tensor<T>(Shape2);
     for (int j = 0; j < Shape2[1]; j++){
      for (int i = 0; i < Shape2[0]; i++){
       GaugeTensor.set(i+j*Shape2[0], Vr(j, i)*sqrt(-W[i]));
     }}
     Shape2[0] = D; Shape2[1] = dimPos;
     InvGaugeTensors[pos] = Tensor<T>(Shape2);
     for (int j = 0; j < Shape2[1]; j++){
      for (int i = 0; i < Shape2[0]; i++){
       InvGaugeTensors[pos].set(i+j*Shape2[0], MathAuxiliary::complexConjugate(Vr(i, j))/sqrt(-W[j]));
     }}
    }
// 2.a)1. - multiply gauges:
    if (pos == 0)
    {
     Indices01[0] = 1; Indices11[0] = 0;
     GaugeTensor.contract(Indices01, Tensor0, Indices11);
     Tensor0 = GaugeTensor;
    }
    else if (pos == 1)
    {
     Indices01[0] = 1; Indices11[0] = 1;
     Tensor0.contract(Indices01, GaugeTensor, Indices11);
     Order5[0] = 0; Order5[1] = 4; Order5[2] = 1; Order5[3] = 2; Order5[4] = 3;
     Tensor0.permute(Order5);
    }
    else if (pos == 2)
    {
     Indices01[0] = 1; Indices11[0] = 1;
     Tensor1.contract(Indices01, GaugeTensor, Indices11);
     Order5[0] = 0; Order5[1] = 4; Order5[2] = 1; Order5[3] = 2; Order5[4] = 3;
     Tensor1.permute(Order5);
    }
    else if (pos == 3)
    {
     Indices01[0] = 2; Indices11[0] = 1;
     Tensor1.contract(Indices01, GaugeTensor, Indices11);
     Order5[0] = 0; Order5[1] = 1; Order5[2] = 4; Order5[3] = 2; Order5[4] = 3;
     Tensor1.permute(Order5);
    }
    else if (pos == 4)
    {
     Indices01[0] = 3; Indices11[0] = 1;
     Tensor1.contract(Indices01, GaugeTensor, Indices11);
     Order5[0] = 0; Order5[1] = 1; Order5[2] = 2; Order5[3] = 4; Order5[4] = 3;
     Tensor1.permute(Order5);
    }
    else if (pos == 5)
    {
     Indices01[0] = 3; Indices11[0] = 1;
     Tensor0.contract(Indices01, GaugeTensor, Indices11);
     Order5[0] = 0; Order5[1] = 1; Order5[2] = 2; Order5[3] = 4; Order5[4] = 3;
     Tensor0.permute(Order5);
    }
   }
  }
// 2.a)2. split tensors and compute reduced tensors Tensor0Red and Tensor1Red:
  Indices3[0] = 0; Indices3[1] = 1; Indices3[2] = 3;
  Indices2[0] = 2; Indices2[1] = 4;
  Tensor0.QRDecompose(Indices3, Indices2, Tensor0Red);
  TensorL = Tensor0;
  Indices2[0] = 0; Indices2[1] = 4;
  Indices3[0] = 1; Indices3[1] = 2; Indices3[2] = 3;
  Tensor1.LQDecompose(Indices2, Indices3, Tensor1Red);
  TensorR = Tensor1;
  Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
  Tensor1Red.permute(Order3);
// 2.a)3. multiply inverse gauges:
  for (int pos = 0; pos < 6; pos++)
  {
   NormMPO.get(pos, EnvTensor);
   EnvTensor.getShape(Shape4);
   if (Shape4[2] == D)
   {
    if (pos == 0)
    {
     Indices01[0] = 1; Indices11[0] = 0;
     InvGaugeTensors[0].contract(Indices01, TensorL, Indices11);
     TensorL = InvGaugeTensors[0];
    }
    else if (pos == 1)
    {
     Indices01[0] = 1; Indices11[0] = 1;
     TensorL.contract(Indices01, InvGaugeTensors[1], Indices11);
     Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
     TensorL.permute(Order4);
    }
    else if (pos == 2)
    {
     Indices01[0] = 1; Indices11[0] = 1;
     TensorR.contract(Indices01, InvGaugeTensors[2], Indices11);
     Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
     TensorR.permute(Order4);
    }
    else if (pos == 3)
    {
     Indices01[0] = 2; Indices11[0] = 1;
     TensorR.contract(Indices01, InvGaugeTensors[3], Indices11);
     Order4[0] = 0; Order4[1] = 1; Order4[2] = 3; Order4[3] = 2;
     TensorR.permute(Order4);
    }
    else if (pos == 4)
    {
     Indices01[0] = 3; Indices11[0] = 1;
     TensorR.contract(Indices01, InvGaugeTensors[4], Indices11);
    }
    else if (pos == 5)
    {
     Indices01[0] = 2; Indices11[0] = 1;
     TensorL.contract(Indices01, InvGaugeTensors[5], Indices11);
     Order4[0] = 0; Order4[1] = 1; Order4[2] = 3; Order4[3] = 2;
     TensorL.permute(Order4);
    }
   }
  }
// 2.a)4. minimize error:
  TensorPair = Tensor0Red;
  Indices01[0] = 1; Indices11[0] = 0;
  TensorPair.contract(Indices01, Tensor1Red, Indices11);
  Indices01[0] = 0;
  TETensor0.contract(Indices01, TETensor1, Indices11);
  Indices02[0] = 1; Indices02[1] = 3; Indices12[0] = 0; Indices12[1] = 2;
  TensorPair.contract(Indices02, TETensor0, Indices12);
  Indices02[0] = 0; Indices02[1] = 2; Indices12[0] = 1; Indices12[1] = 3;
  TensorPair.singularValueDecompose(Indices02, Indices12, D, Tensor0Red, Lambda, Tensor1Red);
  for (int i = 0; i < D; i++)
   Lambda(i, i) = sqrt(Lambda(i, i));
  Tensor1 = Lambda;
  Indices01[0] = 2; Indices11[0] = 0;
  Tensor0Red.contract(Indices01, Tensor1, Indices11);
  Tensor0 = Lambda;
  Indices01[0] = 1; Indices11[0] = 0;
  Tensor0.contract(Indices01, Tensor1Red, Indices11);
  Tensor1Red = Tensor0;
// 2.a)5. recover Tensor0 and Tensor1:
  Tensor0 = TensorL;
  Indices01[0] = 3; Indices11[0] = 0;
  Tensor0.contract(Indices01, Tensor0Red, Indices11);
  Order5[0] = 0; Order5[1] = 1; Order5[2] = 4; Order5[3] = 2; Order5[4] = 3;
  Tensor0.permute(Order5);
  Tensor1 = Tensor1Red;
  Indices01[0] = 1; Indices11[0] = 0;
  Tensor1.contract(Indices01, TensorR, Indices11);
  Order5[0] = 0; Order5[1] = 2; Order5[2] = 3; Order5[3] = 4; Order5[4] = 1;
  Tensor1.permute(Order5);
// 2.a)6. UpdateMode == "Gauge": put tensors on same footing
  if (UpdateMode == "Gauge")
  {
   Indices4[0] = 0; Indices4[1] = 1; Indices4[2] = 3; Indices4[3] = 4;
   Indices1[0] = 2;
   Tensor0.QRDecompose(Indices4, Indices1, TensorR);
   Indices1[0] = 0;
   Indices4[0] = 1; Indices4[1] = 2; Indices4[2] = 3; Indices4[3] = 4;
   Tensor1.LQDecompose(Indices1, Indices4, TensorL);
   Indices01[0] = 1; Indices11[0] = 0;
   TensorR.contract(Indices01, TensorL, Indices11);
   Indices01[0] = 0; Indices11[0] = 1;
   TensorR.singularValueDecompose(Indices01, Indices11, D, Tensor0Red, Lambda, Tensor1Red);
   for (int i = 0; i < D; i++)
    Lambda(i, i) = sqrt(Lambda(i, i));
   Tensor2 = Lambda;
   Indices01[0] = 1; Indices11[0] = 0;
   Tensor0Red.contract(Indices01, Tensor2, Indices11);
   Tensor2 = Lambda;
   Tensor2.contract(Indices01, Tensor1Red, Indices11);
   Indices01[0] = 4; Indices11[0] = 0;
   Tensor0.contract(Indices01, Tensor0Red, Indices11);
   Order5[0] = 0; Order5[1] = 1; Order5[2] = 4; Order5[3] = 2; Order5[4] = 3;
   Tensor0.permute(Order5);
   Indices01[0] = 1; Indices11[0] = 0;
   Tensor2.contract(Indices01, Tensor1, Indices11);
   Tensor1 = Tensor2;
  }
 }
// 2.b) nonseparable environment and UpdateTensor == "reduced":
 else if ((NormMPO.getD() != 1) && (UpdateTensor == "reduced"))
 {
// 2.b)1. split tensors and compute reduced tensors Tensor0Red and Tensor1Red:
  Indices3[0] = 0; Indices3[1] = 1; Indices3[2] = 3;
  Indices2[0] = 2; Indices2[1] = 4;
  Tensor0.QRDecompose(Indices3, Indices2, Tensor0Red);
  TensorL = Tensor0;
  Indices2[0] = 0; Indices2[1] = 4;
  Indices3[0] = 1; Indices3[1] = 2; Indices3[2] = 3;
  Tensor1.LQDecompose(Indices2, Indices3, Tensor1Red);
  TensorR = Tensor1;
  Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
  Tensor1Red.permute(Order3);
// 2.b)2. compute NormTensor:
// 2.b)2. - left half:
  NormMPO.get(5, NormTensor); NormMPO.get(0, Tensor1);
  Indices01[0] = 1; Indices11[0] = 0;
  NormTensor.contract(Indices01, Tensor1, Indices11);
  Tensor1 = TensorL;
  Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 2; Indices12[1] = 0;
  NormTensor.contract(Indices02, Tensor1, Indices12);
  TensorL.complexConjugate(Tensor1);
  Indices02[1] = 3;
  NormTensor.contract(Indices02, Tensor1, Indices12);
  NormMPO.get(1, Tensor1);
  Indices03[0] = 1; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
  NormTensor.contract(Indices03, Tensor1, Indices13);
// 2.b)2. - right half:
  NormMPO.get(2, Tensor0); NormMPO.get(3, Tensor1);
  Indices01[0] = 1; Indices11[0] = 0;
  Tensor0.contract(Indices01, Tensor1, Indices11);
  Tensor1 = TensorR;
  Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  TensorR.complexConjugate(Tensor1);
  Indices02[1] = 3;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  NormMPO.get(4, Tensor1);
  Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 5; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
  Tensor0.contract(Indices03, Tensor1, Indices13);
// 2.b)2. - contract left with right half and permute NormTensor:
  Indices02[0] = 0; Indices02[1] = 3; Indices12[0] = 3; Indices12[1] = 0;
  NormTensor.contract(Indices02, Tensor0, Indices12);
  Order4[0] = 1; Order4[1] = 3; Order4[2] = 0; Order4[3] = 2;
  NormTensor.permute(Order4);
// 2.b)3. make NormTensor positive:
  if ((UpdateMode == "PositiveNormMatrix") || (UpdateMode == "Gauge"))
  {
   NormTensor.getShape(Shape4);
   dimNormTensor = Shape4[0]*Shape4[1];
   NormTensorMatrix = Matrix<T>(dimNormTensor, dimNormTensor);
   for (int i = 0; i < dimNormTensor*dimNormTensor; i++)
    NormTensorMatrix.set(i, NormTensor.get(i));
   if (UpdateMode == "Gauge")
    NormTensorMatrix.transpose();
   NormTensorMatrix.adjoint(NormTensorMatrixAdj);
   NormTensorMatrix.add(NormTensorMatrixAdj);
   NormTensorMatrix.multiply(0.5);
   NormTensorMatrix.setType("hermitian");
   W = vector<double>(dimNormTensor); Vr = Matrix<T>(dimNormTensor, dimNormTensor);
   NormTensorMatrix.eigenDecompose(W, Vr);
   dimPos = 0;
   for (int i = dimNormTensor-1; i >= 0; i--){
    if (W[i] > 0.0)
     dimPos++;
    else
     break;
   }
   if (dimPos == 0)
   {
    cerr << "Program terminated because of error in function " <<
            "template<class T> void TimeEvolution2D<T>::" <<
            "updateTensor(unsigned int positionRow, unsigned int positionCol, const string& UpdateDirection, " <<
                         "MPO<T>& NormMPO, const string& bEnvironment, const MPO<T>& bMPO, const MPO<T>& TEMPO, " <<
                         "const string& UpdateTensor, const string& UpdateMode, double eps, unsigned int maxNumSweeps, " <<
                         "double cutoff, PEPS<T>& PEPS0) const: " <<
            "All eigenvalues of the norm-matrix were negative. Increase the bond dimension of the environment!" << endl;
    exit(1);
   }
// 2.b)3.a) UpdateMode == "PositiveNormMatrix": only make positive
   if (UpdateMode == "PositiveNormMatrix")
   {
    NormTensorMatrix.fillZeroes();
    for (int k = dimNormTensor-1; k >= 0; k--){
     if (W[k] > 0.0){
      for (int j = 0; j < dimNormTensor; j++){
       for (int i = 0; i < dimNormTensor; i++){
        NormTensorMatrix(i, j) += Vr(i, k)*W[k]*MathAuxiliary::complexConjugate(Vr(j, k));
      }}
     }
     else
      break;
    }
    for (int i = 0; i < dimNormTensor*dimNormTensor; i++)
     NormTensor.set(i, NormTensorMatrix.get(i));
   }
// 2.b)3.b) UpdateMode == "Gauge": make positive and gauge
   else if (UpdateMode == "Gauge")
   {
    Shape2[0] = dimPos; Shape2[1] = dimNormTensor;
    SqrtNormTensor = Tensor<T>(Shape2);
    for (int j = 0; j < Shape2[1]; j++){
     for (int i = 0; i < Shape2[0]; i++){
      SqrtNormTensor.set(i+j*Shape2[0], Vr(j, dimNormTensor-1-i)*sqrt(W[dimNormTensor-1-i]));
    }}
    Shape3[0] = dimPos; Shape3[1] = Shape4[2]; Shape3[2] = Shape4[3];
    SqrtNormTensor.reshape(Shape3);
    Tensor0 = SqrtNormTensor;
    Indices1[0] = 1; Indices2[0] = 0; Indices2[1] = 2;
    Tensor0.LQDecompose(Indices1, Indices2, Tensor1);
    Tensor1.getShape(Shape2);
    Matrix0 = Matrix<T>(Shape2[0], Shape2[1]);
    for (int i = 0; i < Shape2[0]*Shape2[1]; i++)
     Matrix0.set(i, Tensor1.get(i));
    Matrix0.pseudoinvert(cutoff, Matrix0Inv);
    InvGaugeTensors2[0] = Matrix0Inv;
    Indices01[0] = 0; Indices11[0] = 0;
    Tensor1.contract(Indices01, Tensor0Red, Indices11);
    Tensor0Red = Tensor1;
    Indices2[0] = 0; Indices2[1] = 1; Indices1[0] = 2;
    SqrtNormTensor.QRDecompose(Indices2, Indices1, Tensor0);
    Tensor0.getShape(Shape2);
    Matrix0 = Matrix<T>(Shape2[0], Shape2[1]);
    for (int i = 0; i < Shape2[0]*Shape2[1]; i++)
     Matrix0.set(i, Tensor0.get(i));
    Matrix0.pseudoinvert(cutoff, Matrix0Inv);
    InvGaugeTensors2[1] = Matrix0Inv;
    Indices01[0] = 1; Indices11[0] = 1;
    Tensor1Red.contract(Indices01, Tensor0, Indices11);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor1Red.permute(Order3);
    Tensor0 = InvGaugeTensors2[0];
    Tensor0.contract(Indices01, SqrtNormTensor, Indices11);
    NormTensor = Tensor0;
    NormTensor.complexConjugate();
    NormTensor.contract(Indices01, Tensor0, Indices11);
   }
  }
// 2.b)4. compute bTensor:
  bTensor = NormTensor;
  Tensor0 = Tensor0Red; Tensor1 = Tensor1Red;
  Indices01[0] = 1; Indices11[0] = 0;
  Tensor0.contract(Indices01, Tensor1, Indices11);
  Indices02[0] = 2; Indices02[1] = 3; Indices12[0] = 0; Indices12[1] = 2;
  bTensor.contract(Indices02, Tensor0, Indices12);
  Tensor0 = TETensor0; Tensor1 = TETensor1;
  Indices01[0] = 0; Indices11[0] = 0;
  Tensor0.contract(Indices01, Tensor1, Indices11);
  bTensor.contract(Indices02, Tensor0, Indices12);
// 2.b)5. minimize error:
// 2.b)5.0. compute initial error:
  error0 = 1.0e10;
  TensorPair = Tensor0Red; Tensor1 = Tensor1Red;
  Indices01[0] = 1; Indices11[0] = 0;
  TensorPair.contract(Indices01, Tensor1, Indices11);
  Tensor0 = NormTensor; Tensor1 = TensorPair;
  Indices02[0] = 2; Indices02[1] = 3; Indices12[0] = 0; Indices12[1] = 2;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  TensorPair.complexConjugate();
  Tensor1 = TensorPair;
  Indices04[0] = 0; Indices04[1] = 1; Indices04[2] = 2; Indices04[3] = 3; Indices14[0] = 0; Indices14[1] = 2; Indices14[2] = 1; Indices14[3] = 3;
  Tensor0.contract(Indices04, Tensor1, Indices14);
  norm = MathAuxiliary::convertToDouble(Tensor0.get(0));
  Tensor0 = bTensor;
  Tensor0.contract(Indices04, TensorPair, Indices14);
  b = Tensor0.get(0);
  error1 = norm - 2.0*MathAuxiliary::convertToDouble(b);
// 2.b)5.0. - UpdateMode == "Gauge": initialize tensor pair via SVD
  if (UpdateMode == "Gauge")
  {
   TensorPair = Tensor0Red;
   Indices01[0] = 1; Indices11[0] = 0;
   TensorPair.contract(Indices01, Tensor1Red, Indices11);
   Indices01[0] = 0;
   TETensor0.contract(Indices01, TETensor1, Indices11);
   Indices02[0] = 1; Indices02[1] = 3; Indices12[0] = 0; Indices12[1] = 2;
   TensorPair.contract(Indices02, TETensor0, Indices12);
   Indices02[0] = 0; Indices02[1] = 2; Indices12[0] = 1; Indices12[1] = 3;
   TensorPair.singularValueDecompose(Indices02, Indices12, D, Tensor0Red, Lambda, Tensor1Red);
   for (int i = 0; i < D; i++)
    Lambda(i, i) = sqrt(Lambda(i, i));
   Tensor1 = Lambda;
   Indices01[0] = 2; Indices11[0] = 0;
   Tensor0Red.contract(Indices01, Tensor1, Indices11);
   Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
   Tensor0Red.permute(Order3);
   Tensor0 = Lambda;
   Indices01[0] = 1; Indices11[0] = 0;
   Tensor0.contract(Indices01, Tensor1Red, Indices11);
   Tensor1Red = Tensor0;
   TensorPair = Tensor0Red; Tensor1 = Tensor1Red;
   Indices01[0] = 1; Indices11[0] = 0;
   TensorPair.contract(Indices01, Tensor1, Indices11);
   Tensor0 = NormTensor; Tensor1 = TensorPair;
   Indices02[0] = 2; Indices02[1] = 3; Indices12[0] = 0; Indices12[1] = 2;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   TensorPair.complexConjugate();
   Tensor1 = TensorPair;
   Indices04[0] = 0; Indices04[1] = 1; Indices04[2] = 2; Indices04[3] = 3; Indices14[0] = 0; Indices14[1] = 2; Indices14[2] = 1; Indices14[3] = 3;
   Tensor0.contract(Indices04, Tensor1, Indices14);
   norm = MathAuxiliary::convertToDouble(Tensor0.get(0));
   Tensor0 = bTensor;
   Tensor0.contract(Indices04, TensorPair, Indices14);
   b = Tensor0.get(0);
   error0 = error1;
   error1 = norm - 2.0*MathAuxiliary::convertToDouble(b);
  }
// 2.b)5. sweep over tensor pair:
  dim0Red = Tensor0Red.getSize(); Tensor0Red.getShape(Shape0Red);
  dim1Red = Tensor1Red.getSize(); Tensor1Red.getShape(Shape1Red);
  sweep = 0;
  while((abs((error1-error0)/error1) > eps) && (sweep < maxNumSweeps))
  {
// 2.b)5.1. update Tensor0Red:
   if (UpdateMode == "Gauge")
   {
    Indices1[0] = 0; Indices2[0] = 1; Indices2[1] = 2;
    Tensor1Red.LQDecompose(Indices1, Indices2, Tensor1);
   }
   NormTensor0 = NormTensor; Tensor1 = Tensor1Red;
   Indices01[0] = 3; Indices11[0] = 1;
   NormTensor0.contract(Indices01, Tensor1, Indices11);
   Tensor1Red.complexConjugate(Tensor1);
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
   NormTensor0.contract(Indices02, Tensor1, Indices12);
   NormMatrix = Matrix<T>(dim0Red, dim0Red);
   NormMatrix.fillZeroes();
   for (int k = 0; k < Shape0Red[2]; k++){
    for (int j1 = 0; j1 < Shape0Red[1]; j1++){
     Index4[2] = j1;
    for (int j0 = 0; j0 < Shape0Red[0]; j0++){
     Index4[1] = j0;
     index1 = j0 + j1*Shape0Red[0] + k*Shape0Red[0]*Shape0Red[1];
     for (int i1 = 0; i1 < Shape0Red[1]; i1++){
      Index4[3] = i1;
     for (int i0 = 0; i0 < Shape0Red[0]; i0++){
      Index4[0] = i0;
      index0 = i0 + i1*Shape0Red[0] + k*Shape0Red[0]*Shape0Red[1];
      NormMatrix(index0, index1) = NormTensor0.get(Index4);
     }}
    }}
   }
   bTensor0 = bTensor; Tensor1Red.complexConjugate(Tensor1);
   Indices02[0] = 1; Indices02[1] = 3; Indices12[0] = 1; Indices12[1] = 2;
   bTensor0.contract(Indices02, Tensor1, Indices12);
   bVector = vector<T>(dim0Red, 0.0);
   for (int i2 = 0; i2 < Shape0Red[2]; i2++){
    Index3[1] = i2;
   for (int i1 = 0; i1 < Shape0Red[1]; i1++){
    Index3[2] = i1;
   for (int i0 = 0; i0 < Shape0Red[0]; i0++){
    Index3[0] = i0;
    index0 = i0 + i1*Shape0Red[0] + i2*Shape0Red[0]*Shape0Red[1];
    bVector[index0] = bTensor0.get(Index3);
   }}}
   if (UpdateMode == "Pseudoinverse")
   {
    NormMatrix.pseudoinvert(cutoff, NormMatrixInv);
   }
   else if ((UpdateMode == "PositiveNormMatrix") || (UpdateMode == "Gauge"))
   {
    NormMatrix.adjoint(NormMatrixAdj);
    NormMatrix.add(NormMatrixAdj);
    NormMatrix.multiply(0.5);
    NormMatrix.setType("hermitian");
    W = vector<double>(dim0Red); Vr = Matrix<T>(dim0Red, dim0Red);
    NormMatrix.eigenDecompose(W, Vr);
    dimPos = 0;
    for (int i = dim0Red-1; i >= 0; i--){
     if (W[i] > 0.0)
      dimPos++;
     else
      break;
    }
    if (dimPos == 0)
    {
     cerr << "Program terminated because of error in function " <<
             "template<class T> void TimeEvolution2D<T>::" <<
             "updateTensor(unsigned int positionRow, unsigned int positionCol, const string& UpdateDirection, " <<
                          "MPO<T>& NormMPO, const string& bEnvironment, const MPO<T>& bMPO, const MPO<T>& TEMPO, " <<
                          "const string& UpdateTensor, const string& UpdateMode, double eps, unsigned int maxNumSweeps, " <<
                          "double cutoff, PEPS<T>& PEPS0) const: " <<
             "All eigenvalues of the norm-matrix were negative. Increase the bond dimension of the environment!" << endl;
     exit(1);
    }
    NormMatrixInv = Matrix<T>(dim0Red, dim0Red);
    NormMatrixInv.fillZeroes();
    for (int k = dim0Red-1; k >= 0; k--){
     if (W[k]/abs(W[dim0Red-1]) > cutoff){
      for (int j = 0; j < dim0Red; j++){
       for (int i = 0; i < dim0Red; i++){
        NormMatrixInv(i, j) += Vr(i, k)*MathAuxiliary::complexConjugate(Vr(j, k))/W[k];
      }}
     }
     else
      break;
    }
   }
   Vector0 = vector<T>(dim0Red);
   NormMatrixInv.multiply(bVector, Vector0);
   for (int i = 0; i < dim0Red; i++)
    Tensor0Red.set(i, Vector0[i]);
// 2.b)5.2. update Tensor1Red:
   if (UpdateMode == "Gauge")
   {
    Indices2[0] = 0; Indices2[1] = 2; Indices1[0] = 1;
    Tensor0Red.QRDecompose(Indices2, Indices1, Tensor1);
    Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
    Tensor0Red.permute(Order3);
   }
   NormTensor0 = NormTensor; Tensor1 = Tensor0Red;
   Indices01[0] = 2; Indices11[0] = 0;
   NormTensor0.contract(Indices01, Tensor1, Indices11);
   Tensor0Red.complexConjugate(Tensor1);
   Indices02[0] = 0; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 2;
   NormTensor0.contract(Indices02, Tensor1, Indices12);
   NormMatrix = Matrix<T>(dim1Red, dim1Red);
   NormMatrix.fillZeroes();
   for (int k = 0; k < Shape1Red[2]; k++){
    for (int j1 = 0; j1 < Shape1Red[1]; j1++){
     Index4[1] = j1;
    for (int j0 = 0; j0 < Shape1Red[0]; j0++){
     Index4[2] = j0;
     index1 = j0 + j1*Shape1Red[0] + k*Shape1Red[0]*Shape1Red[1];
     for (int i1 = 0; i1 < Shape1Red[1]; i1++){
      Index4[0] = i1;
     for (int i0 = 0; i0 < Shape1Red[0]; i0++){
      Index4[3] = i0;
      index0 = i0 + i1*Shape1Red[0] + k*Shape1Red[0]*Shape1Red[1];
      NormMatrix(index0, index1) = NormTensor0.get(Index4);
     }}
    }}
   }
   bTensor0 = bTensor; Tensor0Red.complexConjugate(Tensor1);
   Indices02[0] = 0; Indices02[1] = 2; Indices12[0] = 0; Indices12[1] = 2;
   bTensor0.contract(Indices02, Tensor1, Indices12);
   bVector = vector<T>(dim1Red, 0.0);
   for (int i2 = 0; i2 < Shape1Red[2]; i2++){
    Index3[1] = i2;
   for (int i1 = 0; i1 < Shape1Red[1]; i1++){
    Index3[0] = i1;
   for (int i0 = 0; i0 < Shape1Red[0]; i0++){
    Index3[2] = i0;
    index0 = i0 + i1*Shape1Red[0] + i2*Shape1Red[0]*Shape1Red[1];
    bVector[index0] = bTensor0.get(Index3);
   }}}
   if (UpdateMode == "Pseudoinverse")
   {
    NormMatrix.pseudoinvert(cutoff, NormMatrixInv);
   }
   else if ((UpdateMode == "PositiveNormMatrix") || (UpdateMode == "Gauge"))
   {
    NormMatrix.adjoint(NormMatrixAdj);
    NormMatrix.add(NormMatrixAdj);
    NormMatrix.multiply(0.5);
    NormMatrix.setType("hermitian");
    W = vector<double>(dim1Red); Vr = Matrix<T>(dim1Red, dim1Red);
    NormMatrix.eigenDecompose(W, Vr);
    dimPos = 0;
    for (int i = dim1Red-1; i >= 0; i--){
     if (W[i] > 0.0)
      dimPos++;
     else
      break;
    }
    if (dimPos == 0)
    {
     cerr << "Program terminated because of error in function " <<
             "template<class T> void TimeEvolution2D<T>::" <<
             "updateTensor(unsigned int positionRow, unsigned int positionCol, const string& UpdateDirection, " <<
                          "MPO<T>& NormMPO, const string& bEnvironment, const MPO<T>& bMPO, const MPO<T>& TEMPO, " <<
                          "const string& UpdateTensor, const string& UpdateMode, double eps, unsigned int maxNumSweeps, " <<
                          "double cutoff, PEPS<T>& PEPS0) const: " <<
             "All eigenvalues of the norm-matrix were negative. Increase the bond dimension of the environment!" << endl;
     exit(1);
    }
    NormMatrixInv = Matrix<T>(dim1Red, dim1Red);
    NormMatrixInv.fillZeroes();
    for (int k = dim1Red-1; k >= 0; k--){
     if (W[k]/abs(W[dim1Red-1]) > cutoff){
      for (int j = 0; j < dim1Red; j++){
       for (int i = 0; i < dim1Red; i++){
        NormMatrixInv(i, j) += Vr(i, k)*MathAuxiliary::complexConjugate(Vr(j, k))/W[k];
      }}
     }
     else
      break;
    }
   }
   Vector0 = vector<T>(dim1Red);
   NormMatrixInv.multiply(bVector, Vector0);
   for (int i = 0; i < dim1Red; i++)
    Tensor1Red.set(i, Vector0[i]);
// 2.b)5.3. compute error:
   Tensor1 = Tensor1Red;
   Indices02[0] = 1; Indices02[1] = 2; Indices12[0] = 1; Indices12[1] = 0;
   NormTensor0.contract(Indices02, Tensor1, Indices12);
   Tensor1Red.complexConjugate(Tensor1);
   Indices03[0] = 0; Indices03[1] = 1; Indices03[2] = 2; Indices13[0] = 1; Indices13[1] = 0; Indices13[2] = 2;
   NormTensor0.contract(Indices03, Tensor1, Indices13);
   norm = MathAuxiliary::convertToDouble(NormTensor0.get(0));
   Tensor1Red.complexConjugate(Tensor1);
   Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 0;
   bTensor0.contract(Indices03, Tensor1, Indices13);
   b = bTensor0.get(0);
   error0 = error1;
   error1 = norm - 2.0*MathAuxiliary::convertToDouble(b);
   sweep++;
  }
// 2.b)6. UpdateMode == "Gauge": multiply inverse gauges and put tensors on same footing
  if (UpdateMode == "Gauge")
  {
   Indices01[0] = 0; Indices11[0] = 0;
   Tensor0Red.contract(Indices01, InvGaugeTensors2[0], Indices11);
   Indices01[0] = 1; Indices11[0] = 1;
   Tensor1Red.contract(Indices01, InvGaugeTensors2[1], Indices11);
   Indices2[0] = 1; Indices2[1] = 2;
   Indices1[0] = 0;
   Tensor0Red.QRDecompose(Indices2, Indices1, Tensor0);
   Indices1[0] = 0;
   Indices2[0] = 1; Indices2[1] = 2;
   Tensor1Red.LQDecompose(Indices1, Indices2, Tensor1);
   Indices01[0] = 1; Indices11[0] = 0;
   Tensor0.contract(Indices01, Tensor1, Indices11);
   Indices01[0] = 0; Indices11[0] = 1;
   Tensor0.singularValueDecompose(Indices01, Indices11, D, Tensor1, Lambda, Tensor2);
   for (int i = 0; i < D; i++)
    Lambda(i, i) = sqrt(Lambda(i, i));
   Tensor3 = Lambda;
   Indices01[0] = 1; Indices11[0] = 0;
   Tensor1.contract(Indices01, Tensor3, Indices11);
   Tensor3 = Lambda;
   Tensor3.contract(Indices01, Tensor2, Indices11);
   Indices01[0] = 2; Indices11[0] = 0;
   Tensor0Red.contract(Indices01, Tensor1, Indices11);
   Order3[0] = 1; Order3[1] = 2; Order3[2] = 0;
   Tensor0Red.permute(Order3);
   Indices01[0] = 1;
   Tensor3.contract(Indices01, Tensor1Red, Indices11);
   Order3[0] = 0; Order3[2] = 1;
   Tensor3.permute(Order3);
   Tensor1Red = Tensor3;
  }
// 2.b)7. recover Tensor0 and Tensor1:
  Tensor0 = TensorL;
  Indices01[0] = 3; Indices11[0] = 0;
  Tensor0.contract(Indices01, Tensor0Red, Indices11);
  Order5[0] = 0; Order5[1] = 1; Order5[2] = 3; Order5[3] = 2; Order5[4] = 4;
  Tensor0.permute(Order5);
  Tensor1 = Tensor1Red;
  Indices01[0] = 1; Indices11[0] = 0;
  Tensor1.contract(Indices01, TensorR, Indices11);
  Order5[0] = 0; Order5[1] = 2; Order5[2] = 3; Order5[3] = 4; Order5[4] = 1;
  Tensor1.permute(Order5);
 }
// 2.c) nonseparable environment and UpdateTensor == "full":
 else if ((NormMPO.getD() != 1) && (UpdateTensor == "full"))
 {
  TensorL = Tensor0;
  TensorR = Tensor1;
  TensorL.getShape(Shape0);
  TensorR.getShape(Shape1);
// 2.c)0. UpdateMode == "Gauge": gauge tensors
  if (UpdateMode == "Gauge")
  {
// 2.c)0.1. get gauges for TensorL:
   NormMPO.get(0, Tensor0); NormMPO.get(1, Tensor1);
   Indices01[0] = 1; Indices11[0] = 0;
   Tensor0.contract(Indices01, Tensor1, Indices11);
   NormMPO.get(2, Tensor1); NormMPO.get(3, Tensor2);
   Tensor1.contract(Indices01, Tensor2, Indices11);
   Tensor2 = TensorR;
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
   Tensor1.contract(Indices02, Tensor2, Indices12);
   TensorR.complexConjugate(Tensor2);
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 4;
   Tensor1.contract(Indices03, Tensor2, Indices13);
   NormMPO.get(4, Tensor2);
   Indices03[2] = 5; Indices13[0] = 0; Indices13[2] = 3;
   Tensor1.contract(Indices03, Tensor2, Indices13);
   NormMPO.get(5, Tensor2);
   Indices01[0] = 3; Indices11[0] = 0;
   Tensor1.contract(Indices01, Tensor2, Indices11);
   Indices02[0] = 0; Indices02[1] = 3; Indices12[0] = 3; Indices12[1] = 0;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   dimNormTensor = Shape0[0]*Shape0[1]*Shape0[2]*Shape0[3];
   NormTensorMatrix = Matrix<T>(dimNormTensor, dimNormTensor);
   NormTensorMatrix.fillZeroes();
   for (int j3 = 0; j3 < Shape0[3]; j3++){
    Index8[6] = j3;
   for (int j2 = 0; j2 < Shape0[2]; j2++){
    Index8[4] = j2;
   for (int j1 = 0; j1 < Shape0[1]; j1++){
    Index8[2] = j1;
   for (int j0 = 0; j0 < Shape0[0]; j0++){
    Index8[0] = j0;
    index1 = j0 + j1*Shape0[0] + j2*Shape0[0]*Shape0[1] + j3*Shape0[0]*Shape0[1]*Shape0[2];
    for (int i3 = 0; i3 < Shape0[3]; i3++){
     Index8[7] = i3;
    for (int i2 = 0; i2 < Shape0[2]; i2++){
     Index8[5] = i2;
    for (int i1 = 0; i1 < Shape0[1]; i1++){
     Index8[3] = i1;
    for (int i0 = 0; i0 < Shape0[0]; i0++){
     Index8[1] = i0;
     index0 = i0 + i1*Shape0[0] + i2*Shape0[0]*Shape0[1] + i3*Shape0[0]*Shape0[1]*Shape0[2];
     NormTensorMatrix(index1, index0) = Tensor0.get(Index8);
    }}}}
   }}}}
   NormTensorMatrix.adjoint(NormTensorMatrixAdj);
   NormTensorMatrix.add(NormTensorMatrixAdj);
   NormTensorMatrix.multiply(0.5);
   NormTensorMatrix.setType("hermitian");
   W = vector<double>(dimNormTensor); Vr = Matrix<T>(dimNormTensor, dimNormTensor);
   NormTensorMatrix.eigenDecompose(W, Vr);
   dimPos = 0;
   for (int i = dimNormTensor-1; i >= 0; i--){
    if (W[i] > 0.0)
     dimPos++;
    else
     break;
   }
   if (dimPos == 0)
   {
    cerr << "Program terminated because of error in function " <<
            "template<class T> void TimeEvolution2D<T>::" <<
            "updateTensor(unsigned int positionRow, unsigned int positionCol, const string& UpdateDirection, " <<
                         "MPO<T>& NormMPO, const string& bEnvironment, const MPO<T>& bMPO, const MPO<T>& TEMPO, " <<
                         "const string& UpdateTensor, const string& UpdateMode, double eps, unsigned int maxNumSweeps, " <<
                         "double cutoff, PEPS<T>& PEPS0) const: " <<
            "All eigenvalues of the norm-matrix were negative. Increase the bond dimension of the environment!" << endl;
    exit(1);
   }
   Shape2[0] = dimPos; Shape2[1] = dimNormTensor;
   SqrtNormTensor = Tensor<T>(Shape2);
   for (int j = 0; j < Shape2[1]; j++){
    for (int i = 0; i < Shape2[0]; i++){
     SqrtNormTensor.set(i+j*Shape2[0], Vr(j, dimNormTensor-1-i)*sqrt(W[dimNormTensor-1-i]));
   }}
   Shape5[0] = dimPos; Shape5[1] = Shape0[0]; Shape5[2] = Shape0[1]; Shape5[3] = Shape0[2]; Shape5[4] = Shape0[3];
   SqrtNormTensor.reshape(Shape5);
   if (Shape0[0] == D)
   {
    Tensor0 = SqrtNormTensor;
    Indices4[0] = 0; Indices4[1] = 2; Indices4[2] = 3; Indices4[3] = 4;
    Indices1[0] = 1;
    Tensor0.QRDecompose(Indices4, Indices1, GaugeTensors[0]);
    GaugeTensors[0].getShape(Shape2);
    Matrix0 = Matrix<T>(Shape2[0], Shape2[1]);
    for (int i = 0; i < Shape2[0]*Shape2[1]; i++)
     Matrix0.set(i, GaugeTensors[0].get(i));
    Matrix0.pseudoinvert(cutoff, Matrix0Inv);
    InvGaugeTensors[0] = Matrix0Inv;
   }
   if (Shape0[1] == D)
   {
    Tensor0 = SqrtNormTensor;
    Indices4[0] = 0; Indices4[1] = 1; Indices4[2] = 3; Indices4[3] = 4;
    Indices1[0] = 2;
    Tensor0.QRDecompose(Indices4, Indices1, GaugeTensors[1]);
    GaugeTensors[1].getShape(Shape2);
    Matrix0 = Matrix<T>(Shape2[0], Shape2[1]);
    for (int i = 0; i < Shape2[0]*Shape2[1]; i++)
     Matrix0.set(i, GaugeTensors[1].get(i));
    Matrix0.pseudoinvert(cutoff, Matrix0Inv);
    InvGaugeTensors[1] = Matrix0Inv;
   }
   if (Shape0[3] == D)
   {
    Indices4[0] = 0; Indices4[1] = 1; Indices4[2] = 2; Indices4[3] = 3;
    Indices1[0] = 4;
    SqrtNormTensor.QRDecompose(Indices4, Indices1, GaugeTensors[5]);
    GaugeTensors[5].getShape(Shape2);
    Matrix0 = Matrix<T>(Shape2[0], Shape2[1]);
    for (int i = 0; i < Shape2[0]*Shape2[1]; i++)
     Matrix0.set(i, GaugeTensors[5].get(i));
    Matrix0.pseudoinvert(cutoff, Matrix0Inv);
    InvGaugeTensors[5] = Matrix0Inv;
   }
// 2.c)0.2. get gauges for TensorR:
   NormMPO.get(5, Tensor0); NormMPO.get(0, Tensor1);
   Indices01[0] = 1; Indices11[0] = 0;
   Tensor0.contract(Indices01, Tensor1, Indices11);
   Tensor1 = TensorL;
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 3; Indices12[1] = 0;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   TensorL.complexConjugate(Tensor1);
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 3; Indices13[1] = 0; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   NormMPO.get(1, Tensor1);
   Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   NormMPO.get(2, Tensor1);
   Indices01[0] = 3; Indices11[0] = 0;
   Tensor0.contract(Indices01, Tensor1, Indices11);
   NormMPO.get(3, Tensor1); NormMPO.get(4, Tensor2);
   Indices01[0] = 1;
   Tensor1.contract(Indices01, Tensor2, Indices11);
   Indices02[0] = 0; Indices02[1] = 3; Indices12[0] = 3; Indices12[1] = 0;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   dimNormTensor = Shape1[0]*Shape1[1]*Shape1[2]*Shape1[3];
   NormTensorMatrix = Matrix<T>(dimNormTensor, dimNormTensor);
   NormTensorMatrix.fillZeroes();
   for (int j3 = 0; j3 < Shape1[3]; j3++){
    Index8[6] = j3;
   for (int j2 = 0; j2 < Shape1[2]; j2++){
    Index8[4] = j2;
   for (int j1 = 0; j1 < Shape1[1]; j1++){
    Index8[2] = j1;
   for (int j0 = 0; j0 < Shape1[0]; j0++){
    Index8[0] = j0;
    index1 = j0 + j1*Shape1[0] + j2*Shape1[0]*Shape1[1] + j3*Shape1[0]*Shape1[1]*Shape1[2];
    for (int i3 = 0; i3 < Shape1[3]; i3++){
     Index8[7] = i3;
    for (int i2 = 0; i2 < Shape1[2]; i2++){
     Index8[5] = i2;
    for (int i1 = 0; i1 < Shape1[1]; i1++){
     Index8[3] = i1;
    for (int i0 = 0; i0 < Shape1[0]; i0++){
     Index8[1] = i0;
     index0 = i0 + i1*Shape1[0] + i2*Shape1[0]*Shape1[1] + i3*Shape1[0]*Shape1[1]*Shape1[2];
     NormTensorMatrix(index1, index0) = Tensor0.get(Index8);
    }}}}
   }}}}
   NormTensorMatrix.adjoint(NormTensorMatrixAdj);
   NormTensorMatrix.add(NormTensorMatrixAdj);
   NormTensorMatrix.multiply(0.5);
   NormTensorMatrix.setType("hermitian");
   W = vector<double>(dimNormTensor); Vr = Matrix<T>(dimNormTensor, dimNormTensor);
   NormTensorMatrix.eigenDecompose(W, Vr);
   dimPos = 0;
   for (int i = dimNormTensor-1; i >= 0; i--){
    if (W[i] > 0.0)
     dimPos++;
    else
     break;
   }
   if (dimPos == 0)
   {
    cerr << "Program terminated because of error in function " <<
            "template<class T> void TimeEvolution2D<T>::" <<
            "updateTensor(unsigned int positionRow, unsigned int positionCol, const string& UpdateDirection, " <<
                         "MPO<T>& NormMPO, const string& bEnvironment, const MPO<T>& bMPO, const MPO<T>& TEMPO, " <<
                         "const string& UpdateTensor, const string& UpdateMode, double eps, unsigned int maxNumSweeps, " <<
                         "double cutoff, PEPS<T>& PEPS0) const: " <<
            "All eigenvalues of the norm-matrix were negative. Increase the bond dimension of the environment!" << endl;
    exit(1);
   }
   Shape2[0] = dimPos; Shape2[1] = dimNormTensor;
   SqrtNormTensor = Tensor<T>(Shape2);
   for (int j = 0; j < Shape2[1]; j++){
    for (int i = 0; i < Shape2[0]; i++){
     SqrtNormTensor.set(i+j*Shape2[0], Vr(j, dimNormTensor-1-i)*sqrt(W[dimNormTensor-1-i]));
   }}
   Shape5[0] = dimPos; Shape5[1] = Shape1[0]; Shape5[2] = Shape1[1]; Shape5[3] = Shape1[2]; Shape5[4] = Shape1[3];
   SqrtNormTensor.reshape(Shape5);
   if (Shape1[1] == D)
   {
    Tensor0 = SqrtNormTensor;
    Indices4[0] = 0; Indices4[1] = 1; Indices4[2] = 3; Indices4[3] = 4;
    Indices1[0] = 2;
    Tensor0.QRDecompose(Indices4, Indices1, GaugeTensors[2]);
    GaugeTensors[2].getShape(Shape2);
    Matrix0 = Matrix<T>(Shape2[0], Shape2[1]);
    for (int i = 0; i < Shape2[0]*Shape2[1]; i++)
     Matrix0.set(i, GaugeTensors[2].get(i));
    Matrix0.pseudoinvert(cutoff, Matrix0Inv);
    InvGaugeTensors[2] = Matrix0Inv;
   }
   if (Shape1[2] == D)
   {
    Tensor0 = SqrtNormTensor;
    Indices4[0] = 0; Indices4[1] = 1; Indices4[2] = 2; Indices4[3] = 4;
    Indices1[0] = 3;
    Tensor0.QRDecompose(Indices4, Indices1, GaugeTensors[3]);
    GaugeTensors[3].getShape(Shape2);
    Matrix0 = Matrix<T>(Shape2[0], Shape2[1]);
    for (int i = 0; i < Shape2[0]*Shape2[1]; i++)
     Matrix0.set(i, GaugeTensors[3].get(i));
    Matrix0.pseudoinvert(cutoff, Matrix0Inv);
    InvGaugeTensors[3] = Matrix0Inv;
   }
   if (Shape1[3] == D)
   {
    Indices4[0] = 0; Indices4[1] = 1; Indices4[2] = 2; Indices4[3] = 3;
    Indices1[0] = 4;
    SqrtNormTensor.QRDecompose(Indices4, Indices1, GaugeTensors[4]);
    GaugeTensors[4].getShape(Shape2);
    Matrix0 = Matrix<T>(Shape2[0], Shape2[1]);
    for (int i = 0; i < Shape2[0]*Shape2[1]; i++)
     Matrix0.set(i, GaugeTensors[4].get(i));
    Matrix0.pseudoinvert(cutoff, Matrix0Inv);
    InvGaugeTensors[4] = Matrix0Inv;
   }
// 2.c)0.3. multiply gauges and inverse gauges:
   if (Shape0[0] == D)
   {
    Indices01[0] = 1; Indices11[0] = 0;
    GaugeTensors[0].contract(Indices01, TensorL, Indices11);
    TensorL = GaugeTensors[0];
    NormMPO.get(0, Tensor0); InvGaugeTensor = InvGaugeTensors[0];
    Indices01[0] = 2;
    Tensor0.contract(Indices01, InvGaugeTensor, Indices11);
    InvGaugeTensors[0].complexConjugate(InvGaugeTensor);
    Tensor0.contract(Indices01, InvGaugeTensor, Indices11);
    NormMPO.set(0, Tensor0);
   }
   if (Shape0[1] == D)
   {
    Indices01[0] = 1; Indices11[0] = 1;
    TensorL.contract(Indices01, GaugeTensors[1], Indices11);
    Order5[0] = 0; Order5[1] = 4; Order5[2] = 1; Order5[3] = 2; Order5[4] = 3;
    TensorL.permute(Order5);
    NormMPO.get(1, Tensor0); InvGaugeTensor = InvGaugeTensors[1];
    Indices01[0] = 2; Indices11[0] = 0;
    Tensor0.contract(Indices01, InvGaugeTensor, Indices11);
    InvGaugeTensors[1].complexConjugate(InvGaugeTensor);
    Tensor0.contract(Indices01, InvGaugeTensor, Indices11);
    NormMPO.set(1, Tensor0);
   }
   if (Shape0[3] == D)
   {
    Indices01[0] = 3; Indices11[0] = 1;
    TensorL.contract(Indices01, GaugeTensors[5], Indices11);
    Order5[0] = 0; Order5[1] = 1; Order5[2] = 2; Order5[3] = 4; Order5[4] = 3;
    TensorL.permute(Order5);
    NormMPO.get(5, Tensor0); InvGaugeTensor = InvGaugeTensors[5];
    Indices01[0] = 2; Indices11[0] = 0;
    Tensor0.contract(Indices01, InvGaugeTensor, Indices11);
    InvGaugeTensors[5].complexConjugate(InvGaugeTensor);
    Tensor0.contract(Indices01, InvGaugeTensor, Indices11);
    NormMPO.set(5, Tensor0);
   }
   if (Shape1[1] == D)
   {
    Indices01[0] = 1; Indices11[0] = 1;
    TensorR.contract(Indices01, GaugeTensors[2], Indices11);
    Order5[0] = 0; Order5[1] = 4; Order5[2] = 1; Order5[3] = 2; Order5[4] = 3;
    TensorR.permute(Order5);
    NormMPO.get(2, Tensor0); InvGaugeTensor = InvGaugeTensors[2];
    Indices01[0] = 2; Indices11[0] = 0;
    Tensor0.contract(Indices01, InvGaugeTensor, Indices11);
    InvGaugeTensors[2].complexConjugate(InvGaugeTensor);
    Tensor0.contract(Indices01, InvGaugeTensor, Indices11);
    NormMPO.set(2, Tensor0);
   }
   if (Shape1[2] == D)
   {
    Indices01[0] = 2; Indices11[0] = 1;
    TensorR.contract(Indices01, GaugeTensors[3], Indices11);
    Order5[0] = 0; Order5[1] = 1; Order5[2] = 4; Order5[3] = 2; Order5[4] = 3;
    TensorR.permute(Order5);
    NormMPO.get(3, Tensor0); InvGaugeTensor = InvGaugeTensors[3];
    Indices01[0] = 2; Indices11[0] = 0;
    Tensor0.contract(Indices01, InvGaugeTensor, Indices11);
    InvGaugeTensors[3].complexConjugate(InvGaugeTensor);
    Tensor0.contract(Indices01, InvGaugeTensor, Indices11);
    NormMPO.set(3, Tensor0);
   }
   if (Shape1[3] == D)
   {
    Indices01[0] = 3; Indices11[0] = 1;
    TensorR.contract(Indices01, GaugeTensors[4], Indices11);
    Order5[0] = 0; Order5[1] = 1; Order5[2] = 2; Order5[3] = 4; Order5[4] = 3;
    TensorR.permute(Order5);
    NormMPO.get(4, Tensor0); InvGaugeTensor = InvGaugeTensors[4];
    Indices01[0] = 2; Indices11[0] = 0;
    Tensor0.contract(Indices01, InvGaugeTensor, Indices11);
    InvGaugeTensors[4].complexConjugate(InvGaugeTensor);
    Tensor0.contract(Indices01, InvGaugeTensor, Indices11);
    NormMPO.set(4, Tensor0);
   }
  }
// 2.c)1. compute bTensor:
// 2.c)1. - left half:
  NormMPO.get(5, bTensor); NormMPO.get(0, Tensor1);
  Indices01[0] = 1; Indices11[0] = 0;
  bTensor.contract(Indices01, Tensor1, Indices11);
  Tensor1 = TensorL;
  Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 3; Indices12[1] = 0;
  bTensor.contract(Indices02, Tensor1, Indices12);
  NormMPO.get(1, Tensor1);
  Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 2;
  bTensor.contract(Indices02, Tensor1, Indices12);
// 2.c)1. - right half:
  NormMPO.get(2, Tensor0); NormMPO.get(3, Tensor1);
  Indices01[0] = 1; Indices11[0] = 0;
  Tensor0.contract(Indices01, Tensor1, Indices11);
  Tensor1 = TensorR;
  Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  NormMPO.get(4, Tensor1);
  Indices02[0] = 2; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 2;
  Tensor0.contract(Indices02, Tensor1, Indices12);
// 2.c)1. - contract left with right half:
  Indices03[0] = 0; Indices03[1] = 3; Indices03[2] = 5; Indices13[0] = 5; Indices13[1] = 3; Indices13[2] = 0;
  bTensor.contract(Indices03, Tensor0, Indices13);
// 2.c)1. - contract with TEMPO:
  Tensor0 = TETensor0; Tensor1 = TETensor1;
  Indices01[0] = 0; Indices11[0] = 0;
  Tensor0.contract(Indices01, Tensor1, Indices11);
  Indices02[0] = 2; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
  bTensor.contract(Indices02, Tensor0, Indices12);
// 2.c)2. minimize error:
// 2.c)2.0. compute initial error:
  error0 = 1.0e10;
  NormMPO.get(5, Tensor0); NormMPO.get(0, Tensor1);
  Indices01[0] = 1; Indices11[0] = 0;
  Tensor0.contract(Indices01, Tensor1, Indices11);
  Tensor1 = TensorL;
  Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 3; Indices12[1] = 0;
  Tensor0.contract(Indices02, Tensor1, Indices12);
  TensorL.complexConjugate(Tensor1);
  Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 3; Indices13[1] = 0; Indices13[2] = 4;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  NormMPO.get(1, Tensor1);
  Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
  Tensor0.contract(Indices03, Tensor1, Indices13);
  NormMPO.get(2, Tensor1); NormMPO.get(3, Tensor2);
  Indices01[0] = 1; Indices11[0] = 0;
  Tensor1.contract(Indices01, Tensor2, Indices11);
  Tensor2 = TensorR;
  Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
  Tensor1.contract(Indices02, Tensor2, Indices12);
  TensorR.complexConjugate(Tensor2);
  Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 4;
  Tensor1.contract(Indices03, Tensor2, Indices13);
  NormMPO.get(4, Tensor2);
  Indices03[2] = 5; Indices13[0] = 0; Indices13[2] = 3;
  Tensor1.contract(Indices03, Tensor2, Indices13);
  Indices04[0] = 0; Indices04[1] = 1; Indices04[2] = 2; Indices04[3] = 3; Indices14[0] = 3; Indices14[1] = 1; Indices14[2] = 2; Indices14[3] = 0;
  Tensor0.contract(Indices04, Tensor1, Indices14);
  norm = MathAuxiliary::convertToDouble(Tensor0.get(0));
  TensorPair = TensorL; Tensor1 = TensorR;
  Indices01[0] = 2; Indices11[0] = 0;
  TensorPair.contract(Indices01, Tensor1, Indices11);
  TensorPair.complexConjugate();
  Tensor0 = bTensor;
  Indices08[0] = 0; Indices08[1] = 1; Indices08[2] = 2; Indices08[3] = 3; Indices08[4] = 4; Indices08[5] = 5; Indices08[6] = 6; Indices08[7] = 7;
  Indices18[0] = 2; Indices18[1] = 0; Indices18[2] = 1; Indices18[3] = 4; Indices18[4] = 5; Indices18[5] = 6; Indices18[6] = 3; Indices18[7] = 7;
  Tensor0.contract(Indices08, TensorPair, Indices18);
  b = Tensor0.get(0);
  error1 = norm - 2.0*MathAuxiliary::convertToDouble(b);
// 2.c)2.0. UpdateMode == "Gauge": initialize tensor pair via SVD
  if (UpdateMode == "Gauge")
  {
// 2.c)2.0. - compute reduced tensors Tensor0Red and Tensor1Red:
   Indices3[0] = 0; Indices3[1] = 1; Indices3[2] = 3;
   Indices2[0] = 2; Indices2[1] = 4;
   TensorL.QRDecompose(Indices3, Indices2, Tensor0Red);
   Indices2[0] = 0; Indices2[1] = 4;
   Indices3[0] = 1; Indices3[1] = 2; Indices3[2] = 3;
   TensorR.LQDecompose(Indices2, Indices3, Tensor1Red);
// 2.c)2.0. - multiply reduced tensor-pair with TEMPO and perform SVD:
   TensorPair = Tensor0Red;
   Indices01[0] = 1; Indices11[0] = 0;
   TensorPair.contract(Indices01, Tensor1Red, Indices11);
   Indices01[0] = 0;
   TETensor0.contract(Indices01, TETensor1, Indices11);
   Indices02[0] = 1; Indices02[1] = 2; Indices12[0] = 0; Indices12[1] = 2;
   TensorPair.contract(Indices02, TETensor0, Indices12);
   Indices02[0] = 0; Indices12[0] = 1; Indices12[1] = 3;
   TensorPair.singularValueDecompose(Indices02, Indices12, D, Tensor0Red, Lambda, Tensor1Red);
   for (int i = 0; i < D; i++)
    Lambda(i, i) = sqrt(Lambda(i, i));
   Tensor1 = Lambda;
   Indices01[0] = 2; Indices11[0] = 0;
   Tensor0Red.contract(Indices01, Tensor1, Indices11);
   Tensor1 = Lambda;
   Indices01[0] = 0; Indices11[0] = 1;
   Tensor1Red.contract(Indices01, Tensor1, Indices11);
// 2.c)2.0 - recover TensorL and TensorR:
   Indices01[0] = 3; Indices11[0] = 0;
   TensorL.contract(Indices01, Tensor0Red, Indices11);
   Order5[0] = 0; Order5[1] = 1; Order5[2] = 4; Order5[3] = 2; Order5[4] = 3;
   TensorL.permute(Order5);
   Indices01[0] = 0; Indices11[0] = 0;
   TensorR.contract(Indices01, Tensor1Red, Indices11);
   Order5[0] = 4; Order5[1] = 0; Order5[2] = 1; Order5[3] = 2; Order5[4] = 3;
   TensorR.permute(Order5);
// 2.c)2.0. - compute error:
   NormMPO.get(5, Tensor0); NormMPO.get(0, Tensor1);
   Indices01[0] = 1; Indices11[0] = 0;
   Tensor0.contract(Indices01, Tensor1, Indices11);
   Tensor1 = TensorL;
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 3; Indices12[1] = 0;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   TensorL.complexConjugate(Tensor1);
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 3; Indices13[1] = 0; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   NormMPO.get(1, Tensor1);
   Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   NormMPO.get(2, Tensor1); NormMPO.get(3, Tensor2);
   Indices01[0] = 1; Indices11[0] = 0;
   Tensor1.contract(Indices01, Tensor2, Indices11);
   Tensor2 = TensorR;
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
   Tensor1.contract(Indices02, Tensor2, Indices12);
   TensorR.complexConjugate(Tensor2);
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 4;
   Tensor1.contract(Indices03, Tensor2, Indices13);
   NormMPO.get(4, Tensor2);
   Indices03[2] = 5; Indices13[0] = 0; Indices13[2] = 3;
   Tensor1.contract(Indices03, Tensor2, Indices13);
   Indices04[0] = 0; Indices04[1] = 1; Indices04[2] = 2; Indices04[3] = 3; Indices14[0] = 3; Indices14[1] = 1; Indices14[2] = 2; Indices14[3] = 0;
   Tensor0.contract(Indices04, Tensor1, Indices14);
   norm = MathAuxiliary::convertToDouble(Tensor0.get(0));
   TensorPair = TensorL; Tensor1 = TensorR;
   Indices01[0] = 2; Indices11[0] = 0;
   TensorPair.contract(Indices01, Tensor1, Indices11);
   TensorPair.complexConjugate();
   Tensor0 = bTensor;
   Indices08[0] = 0; Indices08[1] = 1; Indices08[2] = 2; Indices08[3] = 3; Indices08[4] = 4; Indices08[5] = 5; Indices08[6] = 6; Indices08[7] = 7;
   Indices18[0] = 2; Indices18[1] = 0; Indices18[2] = 1; Indices18[3] = 4; Indices18[4] = 5; Indices18[5] = 6; Indices18[6] = 3; Indices18[7] = 7;
   Tensor0.contract(Indices08, TensorPair, Indices18);
   b = Tensor0.get(0);
   error0 = error1;
   error1 = norm - 2.0*MathAuxiliary::convertToDouble(b);
  }
// 2.c)2. sweep over tensor pair:
  dim0 = TensorL.getSize();
  dim1 = TensorR.getSize();
  sweep = 0;
  while((abs((error1-error0)/error1) > eps) && (sweep < maxNumSweeps))
  {
// 2.c)2.1. update TensorL:
   if (UpdateMode == "Gauge")
   {
    Indices1[0] = 0;
    Indices4[0] = 1; Indices4[1] = 2; Indices4[2] = 3; Indices4[3] = 4;
    TensorR.LQDecompose(Indices1, Indices4, Tensor1);
   }
   NormMPO.get(0, Tensor0); NormMPO.get(1, Tensor1);
   Indices01[0] = 1; Indices11[0] = 0;
   Tensor0.contract(Indices01, Tensor1, Indices11);
   NormMPO.get(2, Tensor1); NormMPO.get(3, Tensor2);
   Tensor1.contract(Indices01, Tensor2, Indices11);
   Tensor2 = TensorR;
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
   Tensor1.contract(Indices02, Tensor2, Indices12);
   TensorR.complexConjugate(Tensor2);
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 1; Indices13[1] = 2; Indices13[2] = 4;
   Tensor1.contract(Indices03, Tensor2, Indices13);
   NormMPO.get(4, Tensor2);
   Indices03[2] = 5; Indices13[0] = 0; Indices13[2] = 3;
   Tensor1.contract(Indices03, Tensor2, Indices13);
   NormMPO.get(5, Tensor2);
   Indices01[0] = 3; Indices11[0] = 0;
   Tensor1.contract(Indices01, Tensor2, Indices11);
   Indices02[0] = 0; Indices02[1] = 3; Indices12[0] = 3; Indices12[1] = 0;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   NormMatrix = Matrix<T>(dim0, dim0);
   NormMatrix.fillZeroes();
   for (int k = 0; k < Shape0[4]; k++){
    for (int j3 = 0; j3 < Shape0[3]; j3++){
     Index8[6] = j3;
    for (int j2 = 0; j2 < Shape0[2]; j2++){
     Index8[4] = j2;
    for (int j1 = 0; j1 < Shape0[1]; j1++){
     Index8[2] = j1;
    for (int j0 = 0; j0 < Shape0[0]; j0++){
     Index8[0] = j0;
     index1 = j0 + j1*Shape0[0] + j2*Shape0[0]*Shape0[1] + j3*Shape0[0]*Shape0[1]*Shape0[2] + k*Shape0[0]*Shape0[1]*Shape0[2]*Shape0[3];
     for (int i3 = 0; i3 < Shape0[3]; i3++){
      Index8[7] = i3;
     for (int i2 = 0; i2 < Shape0[2]; i2++){
      Index8[5] = i2;
     for (int i1 = 0; i1 < Shape0[1]; i1++){
      Index8[3] = i1;
     for (int i0 = 0; i0 < Shape0[0]; i0++){
      Index8[1] = i0;
      index0 = i0 + i1*Shape0[0] + i2*Shape0[0]*Shape0[1] + i3*Shape0[0]*Shape0[1]*Shape0[2] + k*Shape0[0]*Shape0[1]*Shape0[2]*Shape0[3];
      NormMatrix(index0, index1) = Tensor0.get(Index8);
     }}}}
    }}}}
   }
   Tensor0 = bTensor; TensorR.complexConjugate(Tensor1);
   Indices04[0] = 3; Indices04[1] = 4; Indices04[2] = 5; Indices04[3] = 7; Indices14[0] = 1; Indices14[1] = 2; Indices14[2] = 3; Indices14[3] = 4;
   Tensor0.contract(Indices04, Tensor1, Indices14);
   bVector = vector<T>(dim0, 0.0);
   for (int k = 0; k < Shape0[4]; k++){
    Index5[3] = k;
    for (int i3 = 0; i3 < Shape0[3]; i3++){
     Index5[0] = i3;
    for (int i2 = 0; i2 < Shape0[2]; i2++){
     Index5[4] = i2;
    for (int i1 = 0; i1 < Shape0[1]; i1++){
     Index5[2] = i1;
    for (int i0 = 0; i0 < Shape0[0]; i0++){
     Index5[1] = i0;
     index0 = i0 + i1*Shape0[0] + i2*Shape0[0]*Shape0[1] + i3*Shape0[0]*Shape0[1]*Shape0[2] + k*Shape0[0]*Shape0[1]*Shape0[2]*Shape0[3];
     bVector[index0] = Tensor0.get(Index5);
    }}}}
   }
   if (UpdateMode == "Pseudoinverse")
   {
    NormMatrix.pseudoinvert(cutoff, NormMatrixInv);
   }
   else if ((UpdateMode == "PositiveNormMatrix") || (UpdateMode == "Gauge"))
   {
    NormMatrix.adjoint(NormMatrixAdj);
    NormMatrix.add(NormMatrixAdj);
    NormMatrix.multiply(0.5);
    NormMatrix.setType("hermitian");
    W = vector<double>(dim0); Vr = Matrix<T>(dim0, dim0);
    NormMatrix.eigenDecompose(W, Vr);
    dimPos = 0;
    for (int i = dim0-1; i >= 0; i--){
     if (W[i] > 0.0)
      dimPos++;
     else
      break;
    }
    if (dimPos == 0)
    {
     cerr << "Program terminated because of error in function " <<
             "template<class T> void TimeEvolution2D<T>::" <<
             "updateTensor(unsigned int positionRow, unsigned int positionCol, const string& UpdateDirection, " <<
                          "MPO<T>& NormMPO, const string& bEnvironment, const MPO<T>& bMPO, const MPO<T>& TEMPO, " <<
                          "const string& UpdateTensor, const string& UpdateMode, double eps, unsigned int maxNumSweeps, " <<
                          "double cutoff, PEPS<T>& PEPS0) const: " <<
             "All eigenvalues of the norm-matrix were negative. Increase the bond dimension of the environment!" << endl;
     exit(1);
    }
    NormMatrixInv = Matrix<T>(dim0, dim0);
    NormMatrixInv.fillZeroes();
    for (int k = dim0-1; k >= 0; k--){
     if (W[k]/abs(W[dim0-1]) > cutoff){
      for (int j = 0; j < dim0; j++){
       for (int i = 0; i < dim0; i++){
        NormMatrixInv(i, j) += Vr(i, k)*MathAuxiliary::complexConjugate(Vr(j, k))/W[k];
      }}
     }
     else
      break;
    }
   }
   Vector0 = vector<T>(dim0);
   NormMatrixInv.multiply(bVector, Vector0);
   for (int i = 0; i < dim0; i++)
    TensorL.set(i, Vector0[i]);
// 2.c)2.2. update TensorR:
   if (UpdateMode == "Gauge")
   {
    Indices4[0] = 0; Indices4[1] = 1; Indices4[2] = 3; Indices4[3] = 4;
    Indices1[0] = 2;
    TensorL.QRDecompose(Indices4, Indices1, Tensor1);
    Order5[0] = 0; Order5[1] = 1; Order5[2] = 4; Order5[3] = 2; Order5[4] = 3;
    TensorL.permute(Order5);
   }
   NormMPO.get(5, Tensor0); NormMPO.get(0, Tensor1);
   Indices01[0] = 1; Indices11[0] = 0;
   Tensor0.contract(Indices01, Tensor1, Indices11);
   Tensor1 = TensorL;
   Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 3; Indices12[1] = 0;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   TensorL.complexConjugate(Tensor1);
   Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 6; Indices13[0] = 3; Indices13[1] = 0; Indices13[2] = 4;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   NormMPO.get(1, Tensor1);
   Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
   Tensor0.contract(Indices03, Tensor1, Indices13);
   NormMPO.get(2, Tensor1);
   Indices01[0] = 3; Indices11[0] = 0;
   Tensor0.contract(Indices01, Tensor1, Indices11);
   NormMPO.get(3, Tensor1); NormMPO.get(4, Tensor2);
   Indices01[0] = 1;
   Tensor1.contract(Indices01, Tensor2, Indices11);
   Indices02[0] = 0; Indices02[1] = 3; Indices12[0] = 3; Indices12[1] = 0;
   Tensor0.contract(Indices02, Tensor1, Indices12);
   NormMatrix = Matrix<T>(dim1, dim1);
   NormMatrix.fillZeroes();
   for (int k = 0; k < Shape1[4]; k++){
    for (int j3 = 0; j3 < Shape1[3]; j3++){
     Index8[6] = j3;
    for (int j2 = 0; j2 < Shape1[2]; j2++){
     Index8[4] = j2;
    for (int j1 = 0; j1 < Shape1[1]; j1++){
     Index8[2] = j1;
    for (int j0 = 0; j0 < Shape1[0]; j0++){
     Index8[0] = j0;
     index1 = j0 + j1*Shape1[0] + j2*Shape1[0]*Shape1[1] + j3*Shape1[0]*Shape1[1]*Shape1[2] + k*Shape1[0]*Shape1[1]*Shape1[2]*Shape1[3];
     for (int i3 = 0; i3 < Shape1[3]; i3++){
      Index8[7] = i3;
     for (int i2 = 0; i2 < Shape1[2]; i2++){
      Index8[5] = i2;
     for (int i1 = 0; i1 < Shape1[1]; i1++){
      Index8[3] = i1;
     for (int i0 = 0; i0 < Shape1[0]; i0++){
      Index8[1] = i0;
      index0 = i0 + i1*Shape1[0] + i2*Shape1[0]*Shape1[1] + i3*Shape1[0]*Shape1[1]*Shape1[2] + k*Shape1[0]*Shape1[1]*Shape1[2]*Shape1[3];
      NormMatrix(index0, index1) = Tensor0.get(Index8);
     }}}}
    }}}}
   }
   NormMatrix0 = NormMatrix;
   Tensor0 = bTensor; TensorL.complexConjugate(Tensor1);
   Indices04[0] = 0; Indices04[1] = 1; Indices04[2] = 2; Indices04[3] = 6; Indices14[0] = 3; Indices14[1] = 0; Indices14[2] = 1; Indices14[3] = 4;
   Tensor0.contract(Indices04, Tensor1, Indices14);
   bVector = vector<T>(dim1, 0.0);
   for (int k = 0; k < Shape1[4]; k++){
    Index5[3] = k;
    for (int i3 = 0; i3 < Shape1[3]; i3++){
     Index5[2] = i3;
    for (int i2 = 0; i2 < Shape1[2]; i2++){
     Index5[1] = i2;
    for (int i1 = 0; i1 < Shape1[1]; i1++){
     Index5[0] = i1;
    for (int i0 = 0; i0 < Shape1[0]; i0++){
     Index5[4] = i0;
     index0 = i0 + i1*Shape1[0] + i2*Shape1[0]*Shape1[1] + i3*Shape1[0]*Shape1[1]*Shape1[2] + k*Shape1[0]*Shape1[1]*Shape1[2]*Shape1[3];
     bVector[index0] = Tensor0.get(Index5);
    }}}}
   }
   if (UpdateMode == "Pseudoinverse")
   {
    NormMatrix.pseudoinvert(cutoff, NormMatrixInv);
   }
   else if ((UpdateMode == "PositiveNormMatrix") || (UpdateMode == "Gauge"))
   {
    NormMatrix.adjoint(NormMatrixAdj);
    NormMatrix.add(NormMatrixAdj);
    NormMatrix.multiply(0.5);
    NormMatrix.setType("hermitian");
    W = vector<double>(dim1); Vr = Matrix<T>(dim1, dim1);
    NormMatrix.eigenDecompose(W, Vr);
    dimPos = 0;
    for (int i = dim1-1; i >= 0; i--){
     if (W[i] > 0.0)
      dimPos++;
     else
      break;
    }
    if (dimPos == 0)
    {
     cerr << "Program terminated because of error in function " <<
             "template<class T> void TimeEvolution2D<T>::" <<
             "updateTensor(unsigned int positionRow, unsigned int positionCol, const string& UpdateDirection, " <<
                          "MPO<T>& NormMPO, const string& bEnvironment, const MPO<T>& bMPO, const MPO<T>& TEMPO, " <<
                          "const string& UpdateTensor, const string& UpdateMode, double eps, unsigned int maxNumSweeps, " <<
                          "double cutoff, PEPS<T>& PEPS0) const: " <<
             "All eigenvalues of the norm-matrix were negative. Increase the bond dimension of the environment!" << endl;
     exit(1);
    }
    NormMatrixInv = Matrix<T>(dim1, dim1);
    NormMatrixInv.fillZeroes();
    for (int k = dim1-1; k >= 0; k--){
     if (W[k]/abs(W[dim1-1]) > cutoff){
      for (int j = 0; j < dim1; j++){
       for (int i = 0; i < dim1; i++){
        NormMatrixInv(i, j) += Vr(i, k)*MathAuxiliary::complexConjugate(Vr(j, k))/W[k];
      }}
     }
     else
      break;
    }
   }
   Vector0 = vector<T>(dim1);
   NormMatrixInv.multiply(bVector, Vector0);
   for (int i = 0; i < dim1; i++)
    TensorR.set(i, Vector0[i]);
// 2.c)2.3. compute error:
   norm = MathAuxiliary::convertToDouble(expectationValueMatrix(Vector0, NormMatrix0));
   b = scalarProductVector(Vector0, bVector);
   error0 = error1;
   error1 = norm - 2.0*MathAuxiliary::convertToDouble(b);
   sweep++;
  }
// 2.c)3. UpdateMode == "Gauge": multiply inverse gauges and put tensors on same footing
  if (UpdateMode == "Gauge")
  {
   if (Shape0[0] == D)
   {
    Indices01[0] = 1; Indices11[0] = 0;
    InvGaugeTensors[0].contract(Indices01, TensorL, Indices11);
    TensorL = InvGaugeTensors[0];
   }
   if (Shape0[1] == D)
   {
    Indices01[0] = 1; Indices11[0] = 1;
    TensorL.contract(Indices01, InvGaugeTensors[1], Indices11);
    Order5[0] = 0; Order5[1] = 4; Order5[2] = 1; Order5[3] = 2; Order5[4] = 3;
    TensorL.permute(Order5);
   }
   if (Shape0[3] == D)
   {
    Indices01[0] = 3; Indices11[0] = 1;
    TensorL.contract(Indices01, InvGaugeTensors[5], Indices11);
    Order5[0] = 0; Order5[1] = 1; Order5[2] = 2; Order5[3] = 4; Order5[4] = 3;
    TensorL.permute(Order5);
   }
   if (Shape1[1] == D)
   {
    Indices01[0] = 1; Indices11[0] = 1;
    TensorR.contract(Indices01, InvGaugeTensors[2], Indices11);
    Order5[0] = 0; Order5[1] = 4; Order5[2] = 1; Order5[3] = 2; Order5[4] = 3;
    TensorR.permute(Order5);
   }
   if (Shape1[2] == D)
   {
    Indices01[0] = 2; Indices11[0] = 1;
    TensorR.contract(Indices01, InvGaugeTensors[3], Indices11);
    Order5[0] = 0; Order5[1] = 1; Order5[2] = 4; Order5[3] = 2; Order5[4] = 3;
    TensorR.permute(Order5);
   }
   if (Shape1[3] == D)
   {
    Indices01[0] = 3; Indices11[0] = 1;
    TensorR.contract(Indices01, InvGaugeTensors[4], Indices11);
    Order5[0] = 0; Order5[1] = 1; Order5[2] = 2; Order5[3] = 4; Order5[4] = 3;
    TensorR.permute(Order5);
   }
   Indices4[0] = 0; Indices4[1] = 1; Indices4[2] = 3; Indices4[3] = 4;
   Indices1[0] = 2;
   TensorL.QRDecompose(Indices4, Indices1, Tensor0);
   Indices1[0] = 0;
   Indices4[0] = 1; Indices4[1] = 2; Indices4[2] = 3; Indices4[3] = 4;
   TensorR.LQDecompose(Indices1, Indices4, Tensor1);
   Indices01[0] = 1; Indices11[0] = 0;
   Tensor0.contract(Indices01, Tensor1, Indices11);
   Indices01[0] = 0; Indices11[0] = 1;
   Tensor0.singularValueDecompose(Indices01, Indices11, D, Tensor0Red, Lambda, Tensor1Red);
   for (int i = 0; i < D; i++)
    Lambda(i, i) = sqrt(Lambda(i, i));
   Tensor1 = Lambda;
   Indices01[0] = 1; Indices11[0] = 0;
   Tensor0Red.contract(Indices01, Tensor1, Indices11);
   Tensor1 = Lambda;
   Tensor1.contract(Indices01, Tensor1Red, Indices11);
   Indices01[0] = 4; Indices11[0] = 0;
   TensorL.contract(Indices01, Tensor0Red, Indices11);
   Order5[0] = 0; Order5[1] = 1; Order5[2] = 4; Order5[3] = 2; Order5[4] = 3;
   TensorL.permute(Order5);
   Indices01[0] = 1; Indices11[0] = 0;
   Tensor1.contract(Indices01, TensorR, Indices11);
   TensorR = Tensor1;
  }
// 2.c)4. recover Tensor0 and Tensor1:
  Tensor0 = TensorL;
  Tensor1 = TensorR;
 }
// 3. set tensor pair in PEPS0:
 if (UpdateDirection == "vertical")
 {
// 3. - permute tensor pair to conform it to the order of PEPS0:
  Order5[0] = 3; Order5[1] = 0; Order5[2] = 1; Order5[3] = 2; Order5[4] = 4;
  Tensor0.permute(Order5);
  Tensor1.permute(Order5);
  PEPS0.set(positionRow, positionCol, Tensor0);
  PEPS0.set(positionRow+1, positionCol, Tensor1);
 }
 else if (UpdateDirection == "horizontal")
 {
  PEPS0.set(positionRow, positionCol, Tensor0);
  PEPS0.set(positionRow, positionCol+1, Tensor1);
 }
}

template<class T> void TimeEvolution2D<T>::timeEvolve(Hamiltonian2D<T>& H0, const string& RealImaginaryTE, double timeStart, double timeStep, double timeStop,
                                                      const string& SweepUpdate, unsigned int numSweeps,
                                                      const string& bEnvironment, unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv,
                                                      const string& UpdateTensor, const string& UpdateMode, double epsLoc, unsigned int maxNumSweepsLoc, double cutoff,
                                                      PEPS<T>& PEPS0) const
{
 unsigned int col, colStart, colStop, Nrows = PEPS0.getNrows(), Ncols = PEPS0.getNcols(), row, rowStart, rowStop, sweep;
 double time;
 string SweepDirection, UpdateDirection, WhichBoundary;
 PEPS<T> PEPS1;
 vector< MPS<T> > bMPSs, NormMPSs;
 vector< MPO<T> > TEMPOs;
 colStop = Ncols-1;
 rowStop = Nrows-1;
 PEPS0.normalize("horizontal", 1, epsEnv, maxNumSweepsEnv);
 if (SweepUpdate == "independent")
  PEPS1 = PEPS0;
// start time-evolution:
 if (!(H0.isTimeDependent()))
  H0.getTEMPOs(RealImaginaryTE, timeStep, TEMPOs);
 for (time = timeStart; time < timeStop; time += timeStep)
 {
  if (H0.isTimeDependent())
  {
   H0.setTime(time);
   H0.getTEMPOs(RealImaginaryTE, timeStep, TEMPOs);
  }
// 1. vertical updates:
  UpdateDirection = "vertical";
  SweepDirection = "right";
// 1.a) vertical-even updates:
  rowStart = 0;
  for (sweep = 0; sweep < numSweeps; sweep++)
  {
   this->getInitialNormMPSs(PEPS0, UpdateDirection, D2Env, epsEnv, maxNumSweepsEnv, NormMPSs);
// 1.a)I) update column 0:
   col = 0;
   WhichBoundary = "left";
   this->updatePEPSBoundary(WhichBoundary, rowStart, rowStop, TEMPOs, SweepUpdate, NormMPSs, bEnvironment, bMPSs,
                            UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
   this->updateNormMPSsBoundary(PEPS0, WhichBoundary, D2Env, epsEnv, maxNumSweepsEnv, NormMPSs);
// 1.a)II) update columns 1 to colStop-2:
   for (col = 1; col < colStop; col++)
   {
    this->updatePEPSBulk(UpdateDirection, col, rowStart, rowStop, TEMPOs, SweepUpdate, NormMPSs, bEnvironment, bMPSs,
                         UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
    this->updateNormMPSsBulk(PEPS0, SweepDirection, col, D2Env, epsEnv, maxNumSweepsEnv, NormMPSs);
   }
// 1.a)III) update column Ncols-1:
   col = Ncols-1;
   WhichBoundary = "right";
   this->updatePEPSBoundary(WhichBoundary, rowStart, rowStop, TEMPOs, SweepUpdate, NormMPSs, bEnvironment, bMPSs,
                            UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
// 1.a)IV) SweepUpdate == "independent":
   if (SweepUpdate == "independent")
    PEPS0 = PEPS1;
  }
  PEPS0.normalize("horizontal", 1, epsEnv, maxNumSweepsEnv);
// 1.b) vertical-odd updates:
  rowStart = 1;
  for (sweep = 0; sweep < numSweeps; sweep++)
  {
   this->getInitialNormMPSs(PEPS0, UpdateDirection, D2Env, epsEnv, maxNumSweepsEnv, NormMPSs);
// 1.b)I) update column 0:
   col = 0;
   WhichBoundary = "left";
   this->updatePEPSBoundary(WhichBoundary, rowStart, rowStop, TEMPOs, SweepUpdate, NormMPSs, bEnvironment, bMPSs,
                            UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
   this->updateNormMPSsBoundary(PEPS0, WhichBoundary, D2Env, epsEnv, maxNumSweepsEnv, NormMPSs);
// 1.b)II) update columns 1 to colStop-2:
   for (col = 1; col < colStop; col++)
   {
    this->updatePEPSBulk(UpdateDirection, col, rowStart, rowStop, TEMPOs, SweepUpdate, NormMPSs, bEnvironment, bMPSs,
                         UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
    this->updateNormMPSsBulk(PEPS0, SweepDirection, col, D2Env, epsEnv, maxNumSweepsEnv, NormMPSs);
   }
// 1.b)III) update column Ncols-1:
   col = Ncols-1;
   WhichBoundary = "right";
   this->updatePEPSBoundary(WhichBoundary, rowStart, rowStop, TEMPOs, SweepUpdate, NormMPSs, bEnvironment, bMPSs,
                            UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
// 1.b)IV) SweepUpdate == "independent":
   if (SweepUpdate == "independent")
    PEPS0 = PEPS1;
  }
  PEPS0.normalize("horizontal", 1, epsEnv, maxNumSweepsEnv);
// 2. horizontal updates:
  UpdateDirection = "horizontal";
  SweepDirection = "down";
// 2.a) horizontal-even updates:
  colStart = 0;
  for (sweep = 0; sweep < numSweeps; sweep++)
  {
   this->getInitialNormMPSs(PEPS0, UpdateDirection, D2Env, epsEnv, maxNumSweepsEnv, NormMPSs);
// 2.a)I) update row 0:
   row = 0;
   WhichBoundary = "top";
   this->updatePEPSBoundary(WhichBoundary, colStart, colStop, TEMPOs, SweepUpdate, NormMPSs, bEnvironment, bMPSs,
                            UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
   this->updateNormMPSsBoundary(PEPS0, WhichBoundary, D2Env, epsEnv, maxNumSweepsEnv, NormMPSs);
// 2.a)II) update rows 1 to rowStop-2:
   for (row = 1; row < rowStop; row++)
   {
    this->updatePEPSBulk(UpdateDirection, row, colStart, colStop, TEMPOs, SweepUpdate, NormMPSs, bEnvironment, bMPSs,
                         UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
    this->updateNormMPSsBulk(PEPS0, SweepDirection, row, D2Env, epsEnv, maxNumSweepsEnv, NormMPSs);
   }
// 2.a)III) update row Nrows-1:
   row = Nrows-1;
   WhichBoundary = "bottom";
   this->updatePEPSBoundary(WhichBoundary, colStart, colStop, TEMPOs, SweepUpdate, NormMPSs, bEnvironment, bMPSs,
                            UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
// 2.a)IV) SweepUpdate == "independent":
   if (SweepUpdate == "independent")
    PEPS0 = PEPS1;
  }
  PEPS0.normalize("vertical", 1, epsEnv, maxNumSweepsEnv);
// 2.b) horizontal-odd updates:
  colStart = 1;
  for (sweep = 0; sweep < numSweeps; sweep++)
  {
   this->getInitialNormMPSs(PEPS0, UpdateDirection, D2Env, epsEnv, maxNumSweepsEnv, NormMPSs);
// 2.b)I) update row 0:
   row = 0;
   WhichBoundary = "top";
   this->updatePEPSBoundary(WhichBoundary, colStart, colStop, TEMPOs, SweepUpdate, NormMPSs, bEnvironment, bMPSs,
                            UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
   this->updateNormMPSsBoundary(PEPS0, WhichBoundary, D2Env, epsEnv, maxNumSweepsEnv, NormMPSs);
// 2.b)II) update rows 1 to rowStop-2:
   for (row = 1; row < rowStop; row++)
   {
    this->updatePEPSBulk(UpdateDirection, row, colStart, colStop, TEMPOs, SweepUpdate, NormMPSs, bEnvironment, bMPSs,
                         UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
    this->updateNormMPSsBulk(PEPS0, SweepDirection, row, D2Env, epsEnv, maxNumSweepsEnv, NormMPSs);
   }
// 2.b)III) update row Nrows-1:
   row = Nrows-1;
   WhichBoundary = "bottom";
   this->updatePEPSBoundary(WhichBoundary, colStart, colStop, TEMPOs, SweepUpdate, NormMPSs, bEnvironment, bMPSs,
                            UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
// 2.b)IV) SweepUpdate == "independent":
   if (SweepUpdate == "independent")
    PEPS0 = PEPS1;
  }
  PEPS0.normalize("vertical", 1, epsEnv, maxNumSweepsEnv);
 }// loop over time from timeStart to timeStop in steps of timeStep
 if (H0.isTimeDependent())
  H0.setTime(timeStart);
}

template<class T> void TimeEvolution2D<T>::timeEvolve(Hamiltonian2D<T>& H0, const string& RealImaginaryTE, double timeStart, double timeStep, double timeStop,
                                                      const string& SweepUpdate, unsigned int numSweeps,
                                                      const string& bEnvironment, const vector<unsigned int>& D2sEnv, double epsEnv, unsigned int maxNumSweepsEnv,
                                                      const string& UpdateTensor, const string& UpdateMode, double epsLoc, unsigned int maxNumSweepsLoc, double cutoff,
                                                      PEPS<T>& PEPS0) const
{
 unsigned int clusterSize = D2sEnv.size(), col, colStart, colStop, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), row, rowStart, rowStop, sweep;
 double time;
 string SweepDirection, UpdateDirection, WhichBoundary;
 PEPS<T> PEPS1;
 vector< MPS<T> > SeparablebMPSs, SeparableNormMPSs;
 vector< MPO<T> > TEMPOs;
 colStop = Ncols-1;
 rowStop = Nrows-1;
 PEPS0.normalize("horizontal", 1, epsEnv, maxNumSweepsEnv);
 if (SweepUpdate == "independent")
  PEPS1 = PEPS0;
// start time-evolution:
 if (!(H0.isTimeDependent()))
  H0.getTEMPOs(RealImaginaryTE, timeStep, TEMPOs);
 for (time = timeStart; time < timeStop; time += timeStep)
 {
  if (H0.isTimeDependent())
  {
   H0.setTime(time);
   H0.getTEMPOs(RealImaginaryTE, timeStep, TEMPOs);
  }
// 1. vertical updates:
  UpdateDirection = "vertical";
  SweepDirection = "right";
// 1.a) vertical-even updates:
  rowStart = 0;
  for (sweep = 0; sweep < numSweeps; sweep++)
  {
   this->getInitialSeparableNormMPSs(PEPS0, UpdateDirection, clusterSize, epsEnv, maxNumSweepsEnv, SeparableNormMPSs);
// 1.a)I) update column 0:
   col = 0;
   WhichBoundary = "left";
   this->clusterUpdatePEPSBoundary(WhichBoundary, rowStart, rowStop, TEMPOs, SweepUpdate, D2sEnv, epsEnv, maxNumSweepsEnv, SeparableNormMPSs, bEnvironment, SeparablebMPSs,
                                   UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
   this->updateSeparableNormMPSsBoundary(PEPS0, WhichBoundary, clusterSize, epsEnv, maxNumSweepsEnv, SeparableNormMPSs);
// 1.a)II) update columns 1 to colStop-2:
   for (col = 1; col < colStop; col++)
   {
    this->clusterUpdatePEPSBulk(UpdateDirection, col, rowStart, rowStop, TEMPOs, SweepUpdate, D2sEnv, epsEnv, maxNumSweepsEnv, SeparableNormMPSs, bEnvironment, SeparablebMPSs,
                                UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
    this->updateSeparableNormMPSsBulk(PEPS0, SweepDirection, col, clusterSize, epsEnv, maxNumSweepsEnv, SeparableNormMPSs);
   }
// 1.a)III) update column Ncols-1:
   col = Ncols-1;
   WhichBoundary = "right";
   this->clusterUpdatePEPSBoundary(WhichBoundary, rowStart, rowStop, TEMPOs, SweepUpdate, D2sEnv, epsEnv, maxNumSweepsEnv, SeparableNormMPSs, bEnvironment, SeparablebMPSs,
                                   UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
// 1.a)IV) SweepUpdate == "independent":
   if (SweepUpdate == "independent")
    PEPS0 = PEPS1;
  }
  PEPS0.normalize("horizontal", 1, epsEnv, maxNumSweepsEnv);
// 1.b) vertical-odd updates:
  rowStart = 1;
  for (sweep = 0; sweep < numSweeps; sweep++)
  {
   this->getInitialSeparableNormMPSs(PEPS0, UpdateDirection, clusterSize, epsEnv, maxNumSweepsEnv, SeparableNormMPSs);
// 1.b)I) update column 0:
   col = 0;
   WhichBoundary = "left";
   this->clusterUpdatePEPSBoundary(WhichBoundary, rowStart, rowStop, TEMPOs, SweepUpdate, D2sEnv, epsEnv, maxNumSweepsEnv, SeparableNormMPSs, bEnvironment, SeparablebMPSs,
                                   UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
   this->updateSeparableNormMPSsBoundary(PEPS0, WhichBoundary, clusterSize, epsEnv, maxNumSweepsEnv, SeparableNormMPSs);
// 1.b)II) update columns 1 to colStop-2:
   for (col = 1; col < colStop; col++)
   {
    this->clusterUpdatePEPSBulk(UpdateDirection, col, rowStart, rowStop, TEMPOs, SweepUpdate, D2sEnv, epsEnv, maxNumSweepsEnv, SeparableNormMPSs, bEnvironment, SeparablebMPSs,
                                UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
    this->updateSeparableNormMPSsBulk(PEPS0, SweepDirection, col, clusterSize, epsEnv, maxNumSweepsEnv, SeparableNormMPSs);
   }
// 1.b)III) update column Ncols-1:
   col = Ncols-1;
   WhichBoundary = "right";
   this->clusterUpdatePEPSBoundary(WhichBoundary, rowStart, rowStop, TEMPOs, SweepUpdate, D2sEnv, epsEnv, maxNumSweepsEnv, SeparableNormMPSs, bEnvironment, SeparablebMPSs,
                                   UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
// 1.b)IV) SweepUpdate == "independent":
   if (SweepUpdate == "independent")
    PEPS0 = PEPS1;
  }
  PEPS0.normalize("horizontal", 1, epsEnv, maxNumSweepsEnv);
// 2. horizontal updates:
  UpdateDirection = "horizontal";
  SweepDirection = "down";
// 2.a) horizontal-even updates:
  colStart = 0;
  for (sweep = 0; sweep < numSweeps; sweep++)
  {
   this->getInitialSeparableNormMPSs(PEPS0, UpdateDirection, clusterSize, epsEnv, maxNumSweepsEnv, SeparableNormMPSs);
// 2.a)I) update row 0:
   row = 0;
   WhichBoundary = "top";
   this->clusterUpdatePEPSBoundary(WhichBoundary, colStart, colStop, TEMPOs, SweepUpdate, D2sEnv, epsEnv, maxNumSweepsEnv, SeparableNormMPSs, bEnvironment, SeparablebMPSs,
                                   UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
   this->updateSeparableNormMPSsBoundary(PEPS0, WhichBoundary, clusterSize, epsEnv, maxNumSweepsEnv, SeparableNormMPSs);
// 2.a)II) update rows 1 to rowStop-2:
   for (row = 1; row < rowStop; row++)
   {
    this->clusterUpdatePEPSBulk(UpdateDirection, row, colStart, colStop, TEMPOs, SweepUpdate, D2sEnv, epsEnv, maxNumSweepsEnv, SeparableNormMPSs, bEnvironment, SeparablebMPSs,
                                UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
    this->updateSeparableNormMPSsBulk(PEPS0, SweepDirection, row, clusterSize, epsEnv, maxNumSweepsEnv, SeparableNormMPSs);
   }
// 2.a)III) update row Nrows-1:
   row = Nrows-1;
   WhichBoundary = "bottom";
   this->clusterUpdatePEPSBoundary(WhichBoundary, colStart, colStop, TEMPOs, SweepUpdate, D2sEnv, epsEnv, maxNumSweepsEnv, SeparableNormMPSs, bEnvironment, SeparablebMPSs,
                                   UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
// 2.a)IV) SweepUpdate == "independent":
   if (SweepUpdate == "independent")
    PEPS0 = PEPS1;
  }
  PEPS0.normalize("vertical", 1, epsEnv, maxNumSweepsEnv);
// 2.b) horizontal-odd updates:
  colStart = 1;
  for (sweep = 0; sweep < numSweeps; sweep++)
  {
   this->getInitialSeparableNormMPSs(PEPS0, UpdateDirection, clusterSize, epsEnv, maxNumSweepsEnv, SeparableNormMPSs);
// 2.b)I) update row 0:
   row = 0;
   WhichBoundary = "top";
   this->clusterUpdatePEPSBoundary(WhichBoundary, colStart, colStop, TEMPOs, SweepUpdate, D2sEnv, epsEnv, maxNumSweepsEnv, SeparableNormMPSs, bEnvironment, SeparablebMPSs,
                                   UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
   this->updateSeparableNormMPSsBoundary(PEPS0, WhichBoundary, clusterSize, epsEnv, maxNumSweepsEnv, SeparableNormMPSs);
// 2.b)II) update rows 1 to rowStop-2:
   for (row = 1; row < rowStop; row++)
   {
    this->clusterUpdatePEPSBulk(UpdateDirection, row, colStart, colStop, TEMPOs, SweepUpdate, D2sEnv, epsEnv, maxNumSweepsEnv, SeparableNormMPSs, bEnvironment, SeparablebMPSs,
                                UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
    this->updateSeparableNormMPSsBulk(PEPS0, SweepDirection, row, clusterSize, epsEnv, maxNumSweepsEnv, SeparableNormMPSs);
   }
// 2.b)III) update row Nrows-1:
   row = Nrows-1;
   WhichBoundary = "bottom";
   this->clusterUpdatePEPSBoundary(WhichBoundary, colStart, colStop, TEMPOs, SweepUpdate, D2sEnv, epsEnv, maxNumSweepsEnv, SeparableNormMPSs, bEnvironment, SeparablebMPSs,
                                   UpdateTensor, UpdateMode, epsLoc, maxNumSweepsLoc, cutoff, PEPS0, PEPS1);
// 2.b)IV) SweepUpdate == "independent":
   if (SweepUpdate == "independent")
    PEPS0 = PEPS1;
  }
  PEPS0.normalize("vertical", 1, epsEnv, maxNumSweepsEnv);
 }// loop over time from timeStart to timeStop in steps of timeStep
 if (H0.isTimeDependent())
  H0.setTime(timeStart);
}

template<class T> void TimeEvolution2D<T>::timeEvolveTEBD(const string& Type, const string& UpdateTensor,
                                                          Hamiltonian2D<T>& H0, const string& RealImaginaryTE,
                                                          double timeStart, double timeStep, double timeStop,
                                                          PEPS<T>& PEPS0, Matrix< Matrix<T> >& LambdasV,
                                                          Matrix< Matrix<T> >& LambdasH) const
{
 unsigned int Nrows = PEPS0.getNrows(), Ncols = PEPS0.getNcols();
 vector< MPO<T> > TEMPOs;
 double time = timeStart;
 unsigned int row, col;
 if (Type == "uniform")
 {
  if (!(H0.isTimeDependent()))
   H0.getTEMPOs(RealImaginaryTE, timeStep, TEMPOs);
  while (time < timeStop)
  {
   if (H0.isTimeDependent())
   {
    H0.setTime(time);
    H0.getTEMPOs(RealImaginaryTE, timeStep, TEMPOs);
   }
// perform TEBD approximation:
// - vertical-even:
   for (col = 0; col < Ncols; col++)
   {
    for (row = 0; row < Nrows-1; row += 2)
    {
     this->updateTensorTEBD(UpdateTensor, row, col, "vertical", TEMPOs[0], PEPS0, LambdasV, LambdasH);
    }
   }
// - vertical-odd:
   for (col = 0; col < Ncols; col++)
   {
    for (row = 1; row < Nrows-1; row += 2)
    {
     this->updateTensorTEBD(UpdateTensor, row, col, "vertical", TEMPOs[0], PEPS0, LambdasV, LambdasH);
    }
   }
// - horizontal-even:
   col = 0;
   for (row = 0; row < Nrows; row++)
   {
    this->updateTensorTEBD(UpdateTensor, row, col, "horizontal", TEMPOs[1], PEPS0, LambdasV, LambdasH);
   }
   for (col = 2; col < Ncols-2; col += 2)
   {
    for (row = 0; row < Nrows; row++)
    {
     this->updateTensorTEBD(UpdateTensor, row, col, "horizontal", TEMPOs[0], PEPS0, LambdasV, LambdasH);
    }
   }
   if (Ncols%2 == 0)
   {
    col = Ncols-2;
    for (row = 0; row < Nrows; row++)
    {
     this->updateTensorTEBD(UpdateTensor, row, col, "horizontal", TEMPOs[2], PEPS0, LambdasV, LambdasH);
    }
   }
// - horizontal-odd:
   for (col = 1; col < Ncols-1; col += 2)
   {
    for (row = 0; row < Nrows; row++)
    {
     this->updateTensorTEBD(UpdateTensor, row, col, "horizontal", TEMPOs[3], PEPS0, LambdasV, LambdasH);
    }
   }
   time += timeStep;
  }
  if (H0.isTimeDependent())
   H0.setTime(timeStart);
 }
 else if (Type == "general")
 {}
}

template<class T> void TimeEvolution2D<T>::timeEvolveExactly(const vector<T>& Vector0, Hamiltonian2D<T>& H0,
                                                             const string& RealImaginaryTE,
                                                             double timeStart, double timeStep,
                                                             double timeStop,
                                                             vector<T>& Vector1) const
{
 Matrix<unsigned int> d; H0.getd(d);
 unsigned int dim = 1;
 unsigned int Nrows = H0.getNrows(), Ncols = H0.getNcols();
 for (int col = 0; col < Ncols; col++)
 {
  for (int row = 0; row < Nrows; row++)
  {
   dim *= d(row, col);
  }
 }
#ifdef DEBUG
 if ((Vector0.size() != dim) || (Nrows == 0) ||
     ((RealImaginaryTE != "real") && (RealImaginaryTE != "imaginary")) || (Vector1.size() != dim))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void TimeEvolution2D<T>::" <<
          "timeEvolveExactly(const vector<T>& Vector0, const Hamiltonian2D<T>& H0, " <<
                            "const string& RealImaginaryTE, double timeStart, double timeStep, " <<
                            "double timeStop, vector<T>& Vector1) const: " <<
          "((Vector0.size() != dim) || (Nrows == 0) || " <<
           "((RealImaginaryTE != real) && (RealImaginaryTE != imaginary)) || " <<
           "(Vector1.size() != dim))." << endl;
  exit(1);
 }
#endif
 double time = timeStart;
 if (!(H0.isTimeDependent()))
// if H0 is time-independent, determine the exact time evolution operator by diagonalizing H0:
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
// if H0 is time-dependent, do time evolution with a Runge-Kutta algorithm of fourth order:
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

template<class T> T TimeEvolution2D<T>::expectationValue(const PEPS<T>& PEPS0, const Hamiltonian2D<T>& H0, unsigned int D2, double eps, unsigned int maxNumSweeps) const
{
 unsigned int col, colStart, colStop, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), row, rowStart, rowStop;
 string Direction, SweepDirection, WhichBoundary;
 T result = 0.0;
 vector< MPS<T> > NormMPSs;
 vector< MPO<T> > HMPOs;
 colStop = Ncols-1;
 rowStop = Nrows-1;
// 0. get HMPOs:
 H0.getHMPOs(HMPOs);
// 1. vertical expectation values:
 Direction = "vertical";
 SweepDirection = "right";
// 1.a) vertical-even expectation values:
 rowStart = 0;
 this->getInitialNormMPSs(PEPS0, Direction, D2, eps, maxNumSweeps, NormMPSs);
// 1.a)I) contract column 0:
 col = 0;
 WhichBoundary = "left";
 result += this->expectationValuePEPSBoundary(PEPS0, WhichBoundary, rowStart, rowStop, HMPOs, NormMPSs);
 this->updateNormMPSsBoundary(PEPS0, WhichBoundary, D2, eps, maxNumSweeps, NormMPSs);
// 1.a)II) contract columns 1 to colStop-2:
 for (col = 1; col < colStop; col++)
 {
  result += this->expectationValuePEPSBulk(PEPS0, Direction, col, rowStart, rowStop, HMPOs, NormMPSs);
  this->updateNormMPSsBulk(PEPS0, SweepDirection, col, D2, eps, maxNumSweeps, NormMPSs);
 }
// 1.a)III) contract column Ncols-1:
 col = Ncols-1;
 WhichBoundary = "right";
 result += this->expectationValuePEPSBoundary(PEPS0, WhichBoundary, rowStart, rowStop, HMPOs, NormMPSs);
// 1.b) vertical-odd expectation values:
 rowStart = 1;
 this->getInitialNormMPSs(PEPS0, Direction, D2, eps, maxNumSweeps, NormMPSs);
// 1.b)I) contract column 0:
 col = 0;
 WhichBoundary = "left";
 result += this->expectationValuePEPSBoundary(PEPS0, WhichBoundary, rowStart, rowStop, HMPOs, NormMPSs);
 this->updateNormMPSsBoundary(PEPS0, WhichBoundary, D2, eps, maxNumSweeps, NormMPSs);
// 1.b)II) contract columns 1 to colStop-2:
 for (col = 1; col < colStop; col++)
 {
  result += this->expectationValuePEPSBulk(PEPS0, Direction, col, rowStart, rowStop, HMPOs, NormMPSs);
  this->updateNormMPSsBulk(PEPS0, SweepDirection, col, D2, eps, maxNumSweeps, NormMPSs);
 }
// 1.b)III) contract column Ncols-1:
 col = Ncols-1;
 WhichBoundary = "right";
 result += this->expectationValuePEPSBoundary(PEPS0, WhichBoundary, rowStart, rowStop, HMPOs, NormMPSs);
// 2. horizontal expectation values:
 Direction = "horizontal";
 SweepDirection = "down";
// 2.a) horizontal-even expectation values:
 colStart = 0;
 this->getInitialNormMPSs(PEPS0, Direction, D2, eps, maxNumSweeps, NormMPSs);
// 2.a)I) contract row 0:
 row = 0;
 WhichBoundary = "top";
 result += this->expectationValuePEPSBoundary(PEPS0, WhichBoundary, colStart, colStop, HMPOs, NormMPSs);
 this->updateNormMPSsBoundary(PEPS0, WhichBoundary, D2, eps, maxNumSweeps, NormMPSs);
// 2.a)II) contract rows 1 to rowStop-2:
 for (row = 1; row < rowStop; row++)
 {
  result += this->expectationValuePEPSBulk(PEPS0, Direction, row, colStart, colStop, HMPOs, NormMPSs);
  this->updateNormMPSsBulk(PEPS0, SweepDirection, row, D2, eps, maxNumSweeps, NormMPSs);
 }
// 2.a)III) contract row Nrows-1:
 row = Nrows-1;
 WhichBoundary = "bottom";
 result += this->expectationValuePEPSBoundary(PEPS0, WhichBoundary, colStart, colStop, HMPOs, NormMPSs);
// 2.b) horizontal-odd expectation values:
 colStart = 1;
 this->getInitialNormMPSs(PEPS0, Direction, D2, eps, maxNumSweeps, NormMPSs);
// 2.b)I) contract row 0:
 row = 0;
 WhichBoundary = "top";
 result += this->expectationValuePEPSBoundary(PEPS0, WhichBoundary, colStart, colStop, HMPOs, NormMPSs);
 this->updateNormMPSsBoundary(PEPS0, WhichBoundary, D2, eps, maxNumSweeps, NormMPSs);
// 2.b)II) contract rows 1 to rowStop-2:
 for (row = 1; row < rowStop; row++)
 {
  result += this->expectationValuePEPSBulk(PEPS0, Direction, row, colStart, colStop, HMPOs, NormMPSs);
  this->updateNormMPSsBulk(PEPS0, SweepDirection, row, D2, eps, maxNumSweeps, NormMPSs);
 }
// 2.b)III) contract row Nrows-1:
 row = Nrows-1;
 WhichBoundary = "bottom";
 result += this->expectationValuePEPSBoundary(PEPS0, WhichBoundary, colStart, colStop, HMPOs, NormMPSs);
// 3. return result:
 return result;
}

template<class T> T TimeEvolution2D<T>::expectationValue(const PEPS<T>& PEPS0, const Hamiltonian2D<T>& H0, const vector<unsigned int>& D2s, double eps, unsigned int maxNumSweeps) const
{
 unsigned int clusterSize = D2s.size(), col, colStart, colStop, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), row, rowStart, rowStop;
 string Direction, SweepDirection, WhichBoundary;
 T result = 0.0;
 vector< MPS<T> > SeparableNormMPSs;
 vector< MPO<T> > HMPOs;
 colStop = Ncols-1;
 rowStop = Nrows-1;
// 0. get HMPOs:
 H0.getHMPOs(HMPOs);
// 1. vertical expectation values:
 Direction = "vertical";
 SweepDirection = "right";
// 1.a) vertical-even expectation values:
 rowStart = 0;
 this->getInitialSeparableNormMPSs(PEPS0, Direction, clusterSize, eps, maxNumSweeps, SeparableNormMPSs);
// 1.a)I) contract column 0:
 col = 0;
 WhichBoundary = "left";
 result += this->clusterExpectationValuePEPSBoundary(PEPS0, WhichBoundary, rowStart, rowStop, HMPOs, D2s, eps, maxNumSweeps, SeparableNormMPSs);
 this->updateSeparableNormMPSsBoundary(PEPS0, WhichBoundary, clusterSize, eps, maxNumSweeps, SeparableNormMPSs);
// 1.a)II) contract columns 1 to Ncols-2:
 for (col = 1; col < colStop; col++)
 {
  result += this->clusterExpectationValuePEPSBulk(PEPS0, Direction, col, rowStart, rowStop, HMPOs, D2s, eps, maxNumSweeps, SeparableNormMPSs);
  this->updateSeparableNormMPSsBulk(PEPS0, SweepDirection, col, clusterSize, eps, maxNumSweeps, SeparableNormMPSs);
 }
// 1.a)III) contract column Ncols-1:
 col = Ncols-1;
 WhichBoundary = "right";
 result += this->clusterExpectationValuePEPSBoundary(PEPS0, WhichBoundary, rowStart, rowStop, HMPOs, D2s, eps, maxNumSweeps, SeparableNormMPSs);
// 1.b) vertical-odd expectation values:
 rowStart = 1;
 this->getInitialSeparableNormMPSs(PEPS0, Direction, clusterSize, eps, maxNumSweeps, SeparableNormMPSs);
// 1.b)I) contract column 0:
 col = 0;
 WhichBoundary = "left";
 result += this->clusterExpectationValuePEPSBoundary(PEPS0, WhichBoundary, rowStart, rowStop, HMPOs, D2s, eps, maxNumSweeps, SeparableNormMPSs);
 this->updateSeparableNormMPSsBoundary(PEPS0, WhichBoundary, clusterSize, eps, maxNumSweeps, SeparableNormMPSs);
// 1.b)II) contract columns 1 to Ncols-2:
 for (col = 1; col < colStop; col++)
 {
  result += this->clusterExpectationValuePEPSBulk(PEPS0, Direction, col, rowStart, rowStop, HMPOs, D2s, eps, maxNumSweeps, SeparableNormMPSs);
  this->updateSeparableNormMPSsBulk(PEPS0, SweepDirection, col, clusterSize, eps, maxNumSweeps, SeparableNormMPSs);
 }
// 1.b)III) contract column Ncols-1:
 col = Ncols-1;
 WhichBoundary = "right";
 result += this->clusterExpectationValuePEPSBoundary(PEPS0, WhichBoundary, rowStart, rowStop, HMPOs, D2s, eps, maxNumSweeps, SeparableNormMPSs);
// 2. horizontal expectation values:
 Direction = "horizontal";
 SweepDirection = "down";
// 2.a) horizontal-even expectation values:
 colStart = 0;
 this->getInitialSeparableNormMPSs(PEPS0, Direction, clusterSize, eps, maxNumSweeps, SeparableNormMPSs);
// 2.a)I) contract row 0:
 row = 0;
 WhichBoundary = "top";
 result += this->clusterExpectationValuePEPSBoundary(PEPS0, WhichBoundary, colStart, colStop, HMPOs, D2s, eps, maxNumSweeps, SeparableNormMPSs);
 this->updateSeparableNormMPSsBoundary(PEPS0, WhichBoundary, clusterSize, eps, maxNumSweeps, SeparableNormMPSs);
// 2.a)II) contract rows 1 to Nrows-2:
 for (row = 1; row < rowStop; row++)
 {
  result += this->clusterExpectationValuePEPSBulk(PEPS0, Direction, row, colStart, colStop, HMPOs, D2s, eps, maxNumSweeps, SeparableNormMPSs);
  this->updateSeparableNormMPSsBulk(PEPS0, SweepDirection, row, clusterSize, eps, maxNumSweeps, SeparableNormMPSs);
 }
// 2.a)III) contract row Nrows-1:
 row = Nrows-1;
 WhichBoundary = "bottom";
 result += this->clusterExpectationValuePEPSBoundary(PEPS0, WhichBoundary, colStart, colStop, HMPOs, D2s, eps, maxNumSweeps, SeparableNormMPSs);
// 2.b) horizontal-odd expectation values:
 colStart = 1;
 this->getInitialSeparableNormMPSs(PEPS0, Direction, clusterSize, eps, maxNumSweeps, SeparableNormMPSs);
// 2.b)I) contract row 0:
 row = 0;
 WhichBoundary = "top";
 result += this->clusterExpectationValuePEPSBoundary(PEPS0, WhichBoundary, colStart, colStop, HMPOs, D2s, eps, maxNumSweeps, SeparableNormMPSs);
 this->updateSeparableNormMPSsBoundary(PEPS0, WhichBoundary, clusterSize, eps, maxNumSweeps, SeparableNormMPSs);
// 2.b)II) contract rows 1 to Nrows-2:
 for (row = 1; row < rowStop; row++)
 {
  result += this->clusterExpectationValuePEPSBulk(PEPS0, Direction, row, colStart, colStop, HMPOs, D2s, eps, maxNumSweeps, SeparableNormMPSs);
  this->updateSeparableNormMPSsBulk(PEPS0, SweepDirection, row, clusterSize, eps, maxNumSweeps, SeparableNormMPSs);
 }
// 2.b)III) contract row Nrows-1:
 row = Nrows-1;
 WhichBoundary = "bottom";
 result += this->clusterExpectationValuePEPSBoundary(PEPS0, WhichBoundary, colStart, colStop, HMPOs, D2s, eps, maxNumSweeps, SeparableNormMPSs);
// 3. return result:
 return result;
}

template<class T> TimeEvolution2D<T>::TimeEvolution2D() {}

template<class T> TimeEvolution2D<T>::TimeEvolution2D(const TimeEvolution2D<T>& TimeEvolution0) {}

template<class T> TimeEvolution2D<T>& TimeEvolution2D<T>::operator=(const TimeEvolution2D<T>& TimeEvolution0)
{}

template<class T> void TimeEvolution2D<T>::updateTensorTEBD(const string& UpdateTensor,
                                                            unsigned int positionRow, unsigned int positionCol,
                                                            const string& Direction, const MPO<T>& TEMPO,
                                                            PEPS<T>& PEPS0, Matrix< Matrix<T> >& LambdasV,
                                                            Matrix< Matrix<T> >& LambdasH) const
{
 unsigned int Nrows = PEPS0.getNrows(), Ncols = PEPS0.getNcols(), D = PEPS0.getD();
 vector<unsigned int> Indices0(1), Indices1(1), Order(5);
 vector<unsigned int> Shape3(3), Shape4(4), Indices02(2), Indices12(2), Indices04(4), Indices14(4);
 Tensor<T> Tensor0, Tensor1, TETensor0, TETensor1, Tensor2;
 Matrix<T> Lambda, LambdaInv(D, D);
// reshape TETensors:
 TEMPO.get(0, TETensor0); TEMPO.get(1, TETensor1);
 TETensor0.getShape(Shape4);
 Shape3[0] = Shape4[1]; Shape3[1] = Shape4[2]; Shape3[2] = Shape4[3];
 TETensor0.reshape(Shape3);
 TETensor1.getShape(Shape4);
 Shape3[0] = Shape4[0]; Shape3[1] = Shape4[2]; Shape3[2] = Shape4[3];
 TETensor1.reshape(Shape3);
 if (Direction == "vertical")
 {
// obtain upper tensor Tensor0 and lower tensor Tensor1:
  PEPS0.get(positionRow, positionCol, Tensor0);
  PEPS0.get(positionRow+1, positionCol, Tensor1);
// multiply lambda matrices:
  if (positionCol != 0)
  {
   Indices0[0] = 0; Indices1[0] = 1;
   Lambda = LambdasH(positionRow, positionCol-1);
   Tensor0.contract(Indices0, Lambda, Indices1);
   Order[0] = 4; Order[1] = 0; Order[2] = 1; Order[3] = 2; Order[4] = 3;
   Tensor0.permute(Order);
   Lambda = LambdasH(positionRow+1, positionCol-1);
   Tensor1.contract(Indices0, Lambda, Indices1);
   Tensor1.permute(Order);
  }
  if (positionRow != 0)
  {
   Indices0[0] = 1; Indices1[0] = 1;
   Lambda = LambdasV(positionRow-1, positionCol);
   Tensor0.contract(Indices0, Lambda, Indices1);
   Order[0] = 0; Order[1] = 4; Order[2] = 1; Order[3] = 2; Order[4] = 3;
   Tensor0.permute(Order);
  }
  if (positionCol != Ncols-1)
  {
   Indices0[0] = 2; Indices1[0] = 0;
   Lambda = LambdasH(positionRow, positionCol);
   Tensor0.contract(Indices0, Lambda, Indices1);
   Order[0] = 0; Order[1] = 1; Order[2] = 4; Order[3] = 2; Order[4] = 3;
   Tensor0.permute(Order);
   Lambda = LambdasH(positionRow+1, positionCol);
   Tensor1.contract(Indices0, Lambda, Indices1);
   Tensor1.permute(Order);
  }
  if (positionRow+1 != Nrows-1)
  {
   Indices0[0] = 3; Indices1[0] = 0;
   Lambda = LambdasV(positionRow+1, positionCol);
   Tensor1.contract(Indices0, Lambda, Indices1);
   Order[0] = 0; Order[1] = 1; Order[2] = 2; Order[3] = 4; Order[4] = 3;
   Tensor1.permute(Order);
  }
// multiply middle lambda matrix:
  Indices0[0] = 3; Indices1[0] = 0;
  Lambda = LambdasV(positionRow, positionCol);
  Tensor0.contract(Indices0, Lambda, Indices1);
// perform SVD:
  if (UpdateTensor == "full")
  {
   Indices1[0] = 1;
   Tensor0.contract(Indices0, TETensor0, Indices1);
   Indices0[0] = 4;
   Tensor1.contract(Indices0, TETensor1, Indices1);
   Indices02[0] = 3; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 4;
   Tensor2 = Tensor0;
   Tensor2.contract(Indices02, Tensor1, Indices12);
   Indices04[0] = 0; Indices04[1] = 1; Indices04[2] = 2; Indices04[3] = 3;
   Indices14[0] = 4; Indices14[1] = 5; Indices14[2] = 6; Indices14[3] = 7;
   Tensor2.singularValueDecompose(Indices04, Indices14, D, Tensor0, Lambda, Tensor1);
   Order[0] = 0; Order[1] = 1; Order[2] = 2; Order[3] = 4; Order[4] = 3;
   Tensor0.permute(Order);
   Order[0] = 1; Order[1] = 0; Order[2] = 2; Order[3] = 3; Order[4] = 4;
   Tensor1.permute(Order);
  }
  else if (UpdateTensor == "reduced")
  {
   vector<unsigned int> Indices3(3), Indices2(2), Shape0, Shape1;
   Tensor0.getShape(Shape0); Tensor1.getShape(Shape1);
   Tensor<T> Tensor0U, Tensor0D, Tensor1U, Tensor1D;
   Matrix<T> Lambda0, Lambda1;
   Indices3[0] = 0; Indices3[1] = 1; Indices3[2] = 2;
   Indices2[0] = 3; Indices2[1] = 4;
   Tensor0.singularValueDecompose(Indices3, Indices2, min(Shape0[3]*Shape0[4], Shape0[0]*Shape0[1]*Shape0[2]),
                                  Tensor0U, Lambda0, Tensor0D);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0D.contract(Indices0, Lambda0, Indices1);
   Indices2[0] = 1; Indices2[1] = 4;
   Indices3[0] = 0; Indices3[1] = 2; Indices3[2] = 3;
   Tensor1.singularValueDecompose(Indices2, Indices3, min(Shape1[1]*Shape1[4], Shape1[0]*Shape1[2]*Shape1[3]),
                                  Tensor1U, Lambda1, Tensor1D);
   Indices0[0] = 2; Indices1[0] = 0;
   Tensor1U.contract(Indices0, Lambda1, Indices1);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0D.contract(Indices0, TETensor0, Indices1);
   Indices0[0] = 1;
   Tensor1U.contract(Indices0, TETensor1, Indices1);
   Indices02[0] = 0; Indices02[1] = 2; Indices12[0] = 0; Indices12[1] = 2;
   Tensor2 = Tensor0D;
   Tensor2.contract(Indices02, Tensor1U, Indices12);
   Indices02[0] = 0; Indices02[1] = 1; Indices12[0] = 2; Indices12[1] = 3;
   Tensor2.singularValueDecompose(Indices02, Indices12, D, Tensor0D, Lambda, Tensor1U);
   Tensor0 = Tensor0U;
   Indices0[0] = 3; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor0D, Indices1);
   Order[0] = 0; Order[1] = 1; Order[2] = 2; Order[3] = 4; Order[4] = 3;
   Tensor0.permute(Order);
   Tensor1 = Tensor1U;
   Indices0[0] = 1; Indices1[0] = 0;
   Tensor1.contract(Indices0, Tensor1D, Indices1);
   Order[0] = 2; Order[1] = 0; Order[2] = 3; Order[3] = 4; Order[4] = 1;
   Tensor1.permute(Order);
  }
// set new lambda matrix:
  LambdasV(positionRow, positionCol).fillZeroes();
  for (int i = 0; i < D; i++)
   LambdasV(positionRow, positionCol)(i, i) = Lambda(i, i)/Lambda(0, 0);
// multiply inverse lambda matrices to Tensor0 and Tensor1:
  if (positionCol != 0)
  {
   Indices0[0] = 0; Indices1[0] = 1;
   LambdaInv.fillZeroes();
   for (int i = 0; i < D; i++)
    LambdaInv(i, i) = 1.0/LambdasH(positionRow, positionCol-1)(i, i);
   Tensor0.contract(Indices0, LambdaInv, Indices1);
   Order[0] = 4; Order[1] = 0; Order[2] = 1; Order[3] = 2; Order[4] = 3;
   Tensor0.permute(Order);
   LambdaInv.fillZeroes();
   for (int i = 0; i < D; i++)
    LambdaInv(i, i) = 1.0/LambdasH(positionRow+1, positionCol-1)(i, i);
   Tensor1.contract(Indices0, LambdaInv, Indices1);
   Tensor1.permute(Order);
  }
  if (positionRow != 0)
  {
   Indices0[0] = 1; Indices1[0] = 1;
   LambdaInv.fillZeroes();
   for (int i = 0; i < D; i++)
    LambdaInv(i, i) = 1.0/LambdasV(positionRow-1, positionCol)(i, i);
   Tensor0.contract(Indices0, LambdaInv, Indices1);
   Order[0] = 0; Order[1] = 4; Order[2] = 1; Order[3] = 2; Order[4] = 3;
   Tensor0.permute(Order);
  }
  if (positionCol != Ncols-1)
  {
   Indices0[0] = 2; Indices1[0] = 0;
   LambdaInv.fillZeroes();
   for (int i = 0; i < D; i++)
    LambdaInv(i, i) = 1.0/LambdasH(positionRow, positionCol)(i, i);
   Tensor0.contract(Indices0, LambdaInv, Indices1);
   Order[0] = 0; Order[1] = 1; Order[2] = 4; Order[3] = 2; Order[4] = 3;
   Tensor0.permute(Order);
   LambdaInv.fillZeroes();
   for (int i = 0; i < D; i++)
    LambdaInv(i, i) = 1.0/LambdasH(positionRow+1, positionCol)(i, i);
   Tensor1.contract(Indices0, LambdaInv, Indices1);
   Tensor1.permute(Order);
  }
  if (positionRow+1 != Nrows-1)
  {
   Indices0[0] = 3; Indices1[0] = 0;
   LambdaInv.fillZeroes();
   for (int i = 0; i < D; i++)
    LambdaInv(i, i) = 1.0/LambdasV(positionRow+1, positionCol)(i, i);
   Tensor1.contract(Indices0, LambdaInv, Indices1);
   Order[0] = 0; Order[1] = 1; Order[2] = 2; Order[3] = 4; Order[4] = 3;
   Tensor1.permute(Order);
  }
// set new Tensor0 and Tensor1:
  Tensor0.normalize(); Tensor1.normalize();
  PEPS0.set(positionRow, positionCol, Tensor0);
  PEPS0.set(positionRow+1, positionCol, Tensor1);
 }
 else if (Direction == "horizontal")
 {
// obtain left tensor Tensor0 and right tensor Tensor1:
  PEPS0.get(positionRow, positionCol, Tensor0);
  PEPS0.get(positionRow, positionCol+1, Tensor1);
// multiply lambda matrices:
  if (positionCol != 0)
  {
   Indices0[0] = 0; Indices1[0] = 1;
   Lambda = LambdasH(positionRow, positionCol-1);
   Tensor0.contract(Indices0, Lambda, Indices1);
   Order[0] = 4; Order[1] = 0; Order[2] = 1; Order[3] = 2; Order[4] = 3;
   Tensor0.permute(Order);
  }
  if (positionRow != 0)
  {
   Indices0[0] = 1; Indices1[0] = 1;
   Lambda = LambdasV(positionRow-1, positionCol);
   Tensor0.contract(Indices0, Lambda, Indices1);
   Order[0] = 0; Order[1] = 4; Order[2] = 1; Order[3] = 2; Order[4] = 3;
   Tensor0.permute(Order);
   Lambda = LambdasV(positionRow-1, positionCol+1);
   Tensor1.contract(Indices0, Lambda, Indices1);
   Tensor1.permute(Order);
  }
  if (positionCol+1 != Ncols-1)
  {
   Indices0[0] = 2; Indices1[0] = 0;
   Lambda = LambdasH(positionRow, positionCol+1);
   Tensor1.contract(Indices0, Lambda, Indices1);
   Order[0] = 0; Order[1] = 1; Order[2] = 4; Order[3] = 2; Order[4] = 3;
   Tensor1.permute(Order);
  }
  if (positionRow != Nrows-1)
  {
   Indices0[0] = 3; Indices1[0] = 0;
   Lambda = LambdasV(positionRow, positionCol);
   Tensor0.contract(Indices0, Lambda, Indices1);
   Order[0] = 0; Order[1] = 1; Order[2] = 2; Order[3] = 4; Order[4] = 3;
   Tensor0.permute(Order);
   Lambda = LambdasV(positionRow, positionCol+1);
   Tensor1.contract(Indices0, Lambda, Indices1);
   Tensor1.permute(Order);
  }
// multiply middle lambda matrix:
  Indices0[0] = 2; Indices1[0] = 0;
  Lambda = LambdasH(positionRow, positionCol);
  Tensor0.contract(Indices0, Lambda, Indices1);
// perform SVD:
  if (UpdateTensor == "full")
  {
   Indices0[0] = 3; Indices1[0] = 1;
   Tensor0.contract(Indices0, TETensor0, Indices1);
   Indices0[0] = 4;
   Tensor1.contract(Indices0, TETensor1, Indices1);
   Indices02[0] = 3; Indices02[1] = 4; Indices12[0] = 0; Indices12[1] = 4;
   Tensor2 = Tensor0;
   Tensor2.contract(Indices02, Tensor1, Indices12);
   Indices04[0] = 0; Indices04[1] = 1; Indices04[2] = 2; Indices04[3] = 3;
   Indices14[0] = 4; Indices14[1] = 5; Indices14[2] = 6; Indices14[3] = 7;
   Tensor2.singularValueDecompose(Indices04, Indices14, D, Tensor0, Lambda, Tensor1);
   Order[0] = 0; Order[1] = 1; Order[2] = 4; Order[3] = 2; Order[4] = 3;
   Tensor0.permute(Order);
  }
  else if (UpdateTensor == "reduced")
  {
   vector<unsigned int> Indices3(3), Indices2(2), Shape0, Shape1;
   Tensor0.getShape(Shape0); Tensor1.getShape(Shape1);
   Tensor<T> Tensor0L, Tensor0R, Tensor1L, Tensor1R;
   Matrix<T> Lambda0, Lambda1;
   Indices3[0] = 0; Indices3[1] = 1; Indices3[2] = 2;
   Indices2[0] = 3; Indices2[1] = 4;
   Tensor0.singularValueDecompose(Indices3, Indices2, min(Shape0[3]*Shape0[4], Shape0[0]*Shape0[1]*Shape0[2]),
                                  Tensor0L, Lambda0, Tensor0R);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0R.contract(Indices0, Lambda0, Indices1);
   Indices2[0] = 0; Indices2[1] = 4;
   Indices3[0] = 1; Indices3[1] = 2; Indices3[2] = 3;
   Tensor1.singularValueDecompose(Indices2, Indices3, min(Shape1[0]*Shape1[4], Shape1[1]*Shape1[2]*Shape1[3]),
                                  Tensor1L, Lambda1, Tensor1R);
   Indices0[0] = 2; Indices1[0] = 0;
   Tensor1L.contract(Indices0, Lambda1, Indices1);
   Indices0[0] = 0; Indices1[0] = 1;
   Tensor0R.contract(Indices0, TETensor0, Indices1);
   Indices0[0] = 1;
   Tensor1L.contract(Indices0, TETensor1, Indices1);
   Indices02[0] = 0; Indices02[1] = 2; Indices12[0] = 0; Indices12[1] = 2;
   Tensor2 = Tensor0R;
   Tensor2.contract(Indices02, Tensor1L, Indices12);
   Indices02[0] = 0; Indices02[1] = 1; Indices12[0] = 2; Indices12[1] = 3;
   Tensor2.singularValueDecompose(Indices02, Indices12, D, Tensor0R, Lambda, Tensor1L);
   Tensor0 = Tensor0L;
   Indices0[0] = 3; Indices1[0] = 0;
   Tensor0.contract(Indices0, Tensor0R, Indices1);
   Order[0] = 0; Order[1] = 1; Order[2] = 4; Order[3] = 2; Order[4] = 3;
   Tensor0.permute(Order);
   Tensor1 = Tensor1L;
   Indices0[0] = 1; Indices1[0] = 0;
   Tensor1.contract(Indices0, Tensor1R, Indices1);
   Order[0] = 0; Order[1] = 2; Order[2] = 3; Order[3] = 4; Order[4] = 1;
   Tensor1.permute(Order);
  }
// set new lambda matrix:
  LambdasH(positionRow, positionCol).fillZeroes();
  for (int i = 0; i < D; i++)
   LambdasH(positionRow, positionCol)(i, i) = Lambda(i, i)/Lambda(0, 0);
// multiply inverse lambda matrices to Tensor0 and Tensor1:
  if (positionCol != 0)
  {
   Indices0[0] = 0; Indices1[0] = 1;
   LambdaInv.fillZeroes();
   for (int i = 0; i < D; i++)
    LambdaInv(i, i) = 1.0/LambdasH(positionRow, positionCol-1)(i, i);
   Tensor0.contract(Indices0, LambdaInv, Indices1);
   Order[0] = 4; Order[1] = 0; Order[2] = 1; Order[3] = 2; Order[4] = 3;
   Tensor0.permute(Order);
  }
  if (positionRow != 0)
  {
   Indices0[0] = 1; Indices1[0] = 1;
   LambdaInv.fillZeroes();
   for (int i = 0; i < D; i++)
    LambdaInv(i, i) = 1.0/LambdasV(positionRow-1, positionCol)(i, i);
   Tensor0.contract(Indices0, LambdaInv, Indices1);
   Order[0] = 0; Order[1] = 4; Order[2] = 1; Order[3] = 2; Order[4] = 3;
   Tensor0.permute(Order);
   LambdaInv.fillZeroes();
   for (int i = 0; i < D; i++)
    LambdaInv(i, i) = 1.0/LambdasV(positionRow-1, positionCol+1)(i, i);
   Tensor1.contract(Indices0, LambdaInv, Indices1);
   Tensor1.permute(Order);
  }
  if (positionCol+1 != Ncols-1)
  {
   Indices0[0] = 2; Indices1[0] = 0;
   LambdaInv.fillZeroes();
   for (int i = 0; i < D; i++)
    LambdaInv(i, i) = 1.0/LambdasH(positionRow, positionCol+1)(i, i);
   Tensor1.contract(Indices0, LambdaInv, Indices1);
   Order[0] = 0; Order[1] = 1; Order[2] = 4; Order[3] = 2; Order[4] = 3;
   Tensor1.permute(Order);
  }
  if (positionRow != Nrows-1)
  {
   Indices0[0] = 3; Indices1[0] = 0;
   LambdaInv.fillZeroes();
   for (int i = 0; i < D; i++)
    LambdaInv(i, i) = 1.0/LambdasV(positionRow, positionCol)(i, i);
   Tensor0.contract(Indices0, LambdaInv, Indices1);
   Order[0] = 0; Order[1] = 1; Order[2] = 2; Order[3] = 4; Order[4] = 3;
   Tensor0.permute(Order);
   LambdaInv.fillZeroes();
   for (int i = 0; i < D; i++)
    LambdaInv(i, i) = 1.0/LambdasV(positionRow, positionCol+1)(i, i);
   Tensor1.contract(Indices0, LambdaInv, Indices1);
   Tensor1.permute(Order);
  }
// set new Tensor0 and Tensor1:
  Tensor0.normalize(); Tensor1.normalize();
  PEPS0.set(positionRow, positionCol, Tensor0);
  PEPS0.set(positionRow, positionCol+1, Tensor1);
 }
}

template<class T> T TimeEvolution2D<T>::expectationValuePEPSBoundary(const PEPS<T>& PEPS0, const string& WhichBoundary,
                                                                     unsigned int tensorPositionStart, unsigned int tensorPositionStop,
                                                                     const vector< MPO<T> >& HMPOs, const vector< MPS<T> >& NormMPSs) const
{
 unsigned int col, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), row;
 string Direction;
 T result = 0.0;
 MPO<T> HMPO, NormMPO;
 vector< Tensor<T> > NormTensors;
// a) left boundary column:
 if (WhichBoundary == "left")
 {
  col = 0;
  Direction = "vertical";
// a)1. get initial boundary norm-tensors:
  this->getInitialBoundaryNormTensors(PEPS0, WhichBoundary, NormMPSs[1], NormTensors);
// a)2. contract boundary column:
  for (row = tensorPositionStart; row < tensorPositionStop; row += 2)
  {
// a)2.1. get norm-MPO:
   this->getBoundaryNormMPO(WhichBoundary, NormMPSs[1], NormTensors, row, NormMPO);
// a)2.2. compute expectation value:
   result += this->expectationValue(PEPS0, row, col, Direction, HMPOs[0], NormMPO);
// a)2.3. update norm-tensors:
   this->updateBoundaryNormTensors(PEPS0, WhichBoundary, NormMPSs[1], row, NormTensors);
  }
 }
// b) top boundary row:
 else if (WhichBoundary == "top")
 {
  row = 0;
  Direction = "horizontal";
// b)1. get initial boundary norm-tensors:
  this->getInitialBoundaryNormTensors(PEPS0, WhichBoundary, NormMPSs[1], NormTensors);
// b)2. contract boundary row:
  for (col = tensorPositionStart; col < tensorPositionStop; col += 2)
  {
// b)2.1. get norm-MPO:
   this->getBoundaryNormMPO(WhichBoundary, NormMPSs[1], NormTensors, col, NormMPO);
// b)2.2. compute expectation value:
// b)2.2. - horizontal-even:
   if (tensorPositionStart == 0)
   {
    if (col == 0)
     HMPO = HMPOs[1];
    else if ((col > 0) && (col < Ncols-2))
     HMPO = HMPOs[0];
    else if (col == Ncols-2)
     HMPO = HMPOs[2];
   }
// b)2.2. - horizontal-odd:
   else if (tensorPositionStart == 1)
   {
    HMPO = HMPOs[3];
   }
   result += this->expectationValue(PEPS0, row, col, Direction, HMPO, NormMPO);
// b)2.3. update norm-tensors:
   this->updateBoundaryNormTensors(PEPS0, WhichBoundary, NormMPSs[1], col, NormTensors);
  }
 }
// c) right boundary column:
 else if (WhichBoundary == "right")
 {
  col = Ncols-1;
  Direction = "vertical";
// c)1. get initial boundary norm-tensors:
  this->getInitialBoundaryNormTensors(PEPS0, WhichBoundary, NormMPSs[Ncols-2], NormTensors);
// c)2. contract boundary column:
  for (row = tensorPositionStart; row < tensorPositionStop; row += 2)
  {
// c)2.1. get norm-MPO:
   this->getBoundaryNormMPO(WhichBoundary, NormMPSs[Ncols-2], NormTensors, row, NormMPO);
// c)2.2. compute expectation value:
   result += this->expectationValue(PEPS0, row, col, Direction, HMPOs[0], NormMPO);
// c)2.3. update norm-tensors:
   this->updateBoundaryNormTensors(PEPS0, WhichBoundary, NormMPSs[Ncols-2], row, NormTensors);
  }
 }
// d) bottom boundary row:
 else if (WhichBoundary == "bottom")
 {
  row = Nrows-1;
  Direction = "horizontal";
// d)1. get initial boundary norm-tensors:
  this->getInitialBoundaryNormTensors(PEPS0, WhichBoundary, NormMPSs[Nrows-2], NormTensors);
// d)2. contract boundary row:
  for (col = tensorPositionStart; col < tensorPositionStop; col += 2)
  {
// d)2.1. get norm-MPO:
   this->getBoundaryNormMPO(WhichBoundary, NormMPSs[Nrows-2], NormTensors, col, NormMPO);
// d)2.2. compute expectation value:
// d)2.2. - horizontal-even:
   if (tensorPositionStart == 0)
   {
    if (col == 0)
     HMPO = HMPOs[1];
    else if ((col > 0) && (col < Ncols-2))
     HMPO = HMPOs[0];
    else if (col == Ncols-2)
     HMPO = HMPOs[2];
   }
// d)2.2. - horizontal-odd:
   else if (tensorPositionStart == 1)
   {
    HMPO = HMPOs[3];
   }
   result += this->expectationValue(PEPS0, row, col, Direction, HMPO, NormMPO);
// d)2.3. update norm-tensors:
   this->updateBoundaryNormTensors(PEPS0, WhichBoundary, NormMPSs[Nrows-2], col, NormTensors);
  }
 }
 return result;
}

template<class T> T TimeEvolution2D<T>::expectationValuePEPSBulk(const PEPS<T>& PEPS0, const string& Direction, unsigned int position,
                                                                 unsigned int tensorPositionStart, unsigned int tensorPositionStop,
                                                                 const vector< MPO<T> >& HMPOs, const vector< MPS<T> >& NormMPSs) const
{
 unsigned int col, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), row;
 T result = 0.0;
 MPO<T> HMPO, NormMPO;
 vector< Tensor<T> > NormTensors;
// a) vertical expectation value:
 if (Direction == "vertical")
 {
  col = position;
// a)1. get initial norm-tensors:
  this->getInitialBulkNormTensors(PEPS0, Direction, position, NormMPSs[position-1], NormMPSs[position+1], NormTensors);
// a)2. contract bulk column:
  for (row = tensorPositionStart; row < tensorPositionStop; row += 2)
  {
// a)2.1. get norm-MPO:
   this->getBulkNormMPO(NormMPSs[position-1], NormMPSs[position+1], NormTensors, row, NormMPO);
// a)2.2. compute expectation value:
   result += this->expectationValue(PEPS0, row, col, Direction, HMPOs[0], NormMPO);
// a)2.3. update norm-tensors:
   this->updateBulkNormTensors(PEPS0, Direction, position, NormMPSs[position-1], NormMPSs[position+1], row, NormTensors);
  }
 }
// b) horizontal expectation value:
 else if (Direction == "horizontal")
 {
  row = position;
// b)1. get initial norm-tensors:
  this->getInitialBulkNormTensors(PEPS0, Direction, position, NormMPSs[position+1], NormMPSs[position-1], NormTensors);
// b)2. contract bulk row:
  for (col = tensorPositionStart; col < tensorPositionStop; col += 2)
  {
// b)2.1. get norm-MPO:
   this->getBulkNormMPO(NormMPSs[position+1], NormMPSs[position-1], NormTensors, col, NormMPO);
// b)2.2. compute expectation value:
// b)2.2. - horizontal-even:
   if (tensorPositionStart == 0)
   {
    if (col == 0)
     HMPO = HMPOs[1];
    else if ((col > 0) && (col < Ncols-2))
     HMPO = HMPOs[0];
    else if (col == Ncols-2)
     HMPO = HMPOs[2];
   }
// b)2.2. - horizontal-odd:
   else if (tensorPositionStart == 1)
   {
    HMPO = HMPOs[3];
   }
   result += this->expectationValue(PEPS0, row, col, Direction, HMPO, NormMPO);
// b)2.3. update norm-tensors:
   this->updateBulkNormTensors(PEPS0, Direction, position, NormMPSs[position+1], NormMPSs[position-1], col, NormTensors);
  }
 }
 return result;
}

template<class T> T TimeEvolution2D<T>::clusterExpectationValuePEPSBoundary(const PEPS<T>& PEPS0, const string& WhichBoundary,
                                                                            unsigned int tensorPositionStart, unsigned int tensorPositionStop,
                                                                            const vector< MPO<T> >& HMPOs,
                                                                            const vector<unsigned int>& D2s, double eps, unsigned int maxNumSweeps,
                                                                            const vector< MPS<T> >& SeparableNormMPSs) const
{
 unsigned int col, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), row;
 string Direction;
 T result = 0.0;
 MPS<T> NormMPS;
 MPO<T> HMPO, NormMPO;
 vector< Tensor<T> > NormTensors;
// get boundary-MPS for the boundary and get initial boundary norm-tensors:
 this->getBoundaryNormMPS(PEPS0, WhichBoundary, D2s, eps, maxNumSweeps, SeparableNormMPSs, NormMPS);
 this->getInitialBoundaryNormTensors(PEPS0, WhichBoundary, NormMPS, NormTensors);
// a) left boundary column:
 if (WhichBoundary == "left")
 {
  col = 0;
  Direction = "vertical";
// a) contract boundary column:
  for (row = tensorPositionStart; row < tensorPositionStop; row += 2)
  {
// a)1. get norm-MPO:
   this->getBoundaryNormMPO(WhichBoundary, NormMPS, NormTensors, row, NormMPO);
// a)2. compute expectation value:
   result += this->expectationValue(PEPS0, row, col, Direction, HMPOs[0], NormMPO);
// a)3. update norm-tensors:
   this->updateBoundaryNormTensors(PEPS0, WhichBoundary, NormMPS, row, NormTensors);
  }
 }
// b) top boundary row:
 else if (WhichBoundary == "top")
 {
  row = 0;
  Direction = "horizontal";
// b) contract boundary row:
  for (col = tensorPositionStart; col < tensorPositionStop; col += 2)
  {
// b)1. get norm-MPO:
   this->getBoundaryNormMPO(WhichBoundary, NormMPS, NormTensors, col, NormMPO);
// b)2. compute expectation value:
// b)2. - horizontal-even:
   if (tensorPositionStart == 0)
   {
    if (col == 0)
     HMPO = HMPOs[1];
    else if ((col > 0) && (col < Ncols-2))
     HMPO = HMPOs[0];
    else if (col == Ncols-2)
     HMPO = HMPOs[2];
   }
// b)2. - horizontal-odd:
   else if (tensorPositionStart == 1)
   {
    HMPO = HMPOs[3];
   }
   result += this->expectationValue(PEPS0, row, col, Direction, HMPO, NormMPO);
// b)3. update norm-tensors:
   this->updateBoundaryNormTensors(PEPS0, WhichBoundary, NormMPS, col, NormTensors);
  }
 }
// c) right boundary column:
 else if (WhichBoundary == "right")
 {
  col = Ncols-1;
  Direction = "vertical";
// c) contract boundary column:
  for (row = tensorPositionStart; row < tensorPositionStop; row += 2)
  {
// c)1. get norm-MPO:
   this->getBoundaryNormMPO(WhichBoundary, NormMPS, NormTensors, row, NormMPO);
// c)2. compute expectation value:
   result += this->expectationValue(PEPS0, row, col, Direction, HMPOs[0], NormMPO);
// c)3. update norm-tensors:
   this->updateBoundaryNormTensors(PEPS0, WhichBoundary, NormMPS, row, NormTensors);
  }
 }
// d) bottom boundary row:
 else if (WhichBoundary == "bottom")
 {
  row = Nrows-1;
  Direction = "horizontal";
// d) contract boundary row:
  for (col = tensorPositionStart; col < tensorPositionStop; col += 2)
  {
// d)1. get norm-MPO:
   this->getBoundaryNormMPO(WhichBoundary, NormMPS, NormTensors, col, NormMPO);
// d)2. compute expectation value:
// d)2. - horizontal-even:
   if (tensorPositionStart == 0)
   {
    if (col == 0)
     HMPO = HMPOs[1];
    else if ((col > 0) && (col < Ncols-2))
     HMPO = HMPOs[0];
    else if (col == Ncols-2)
     HMPO = HMPOs[2];
   }
// d)2. - horizontal-odd:
   else if (tensorPositionStart == 1)
   {
    HMPO = HMPOs[3];
   }
   result += this->expectationValue(PEPS0, row, col, Direction, HMPO, NormMPO);
// d)3. update norm-tensors:
   this->updateBoundaryNormTensors(PEPS0, WhichBoundary, NormMPS, col, NormTensors);
  }
 }
 return result;
}

template<class T> T TimeEvolution2D<T>::clusterExpectationValuePEPSBulk(const PEPS<T>& PEPS0, const string& Direction, unsigned int position,
                                                                        unsigned int tensorPositionStart, unsigned int tensorPositionStop,
                                                                        const vector< MPO<T> >& HMPOs,
                                                                        const vector<unsigned int>& D2s, double eps, unsigned int maxNumSweeps,
                                                                        const vector< MPS<T> >& SeparableNormMPSs) const
{
 unsigned int col, Ncols = PEPS0.getNcols(), Nrows = PEPS0.getNrows(), row;
 T result = 0.0;
 MPS<T> NormMPS1, NormMPS2;
 MPO<T> HMPO, NormMPO;
 vector< Tensor<T> > NormTensors;
// get boundary-MPSs for the bulk and get initial bulk norm-tensors:
 this->getBulkNormMPSs(PEPS0, Direction, position, D2s, eps, maxNumSweeps, SeparableNormMPSs, NormMPS1, NormMPS2);
 this->getInitialBulkNormTensors(PEPS0, Direction, position, NormMPS1, NormMPS2, NormTensors);
// a) vertical expectation value:
 if (Direction == "vertical")
 {
  col = position;
// a) contract bulk column:
  for (row = tensorPositionStart; row < tensorPositionStop; row += 2)
  {
// a)1. get norm-MPO:
   this->getBulkNormMPO(NormMPS1, NormMPS2, NormTensors, row, NormMPO);
// a)2. compute expectation value:
   result += this->expectationValue(PEPS0, row, col, Direction, HMPOs[0], NormMPO);
// a)3. update norm-tensors:
   this->updateBulkNormTensors(PEPS0, Direction, position, NormMPS1, NormMPS2, row, NormTensors);
  }
 }
// b) horizontal expectation value:
 else if (Direction == "horizontal")
 {
  row = position;
// b) contract bulk row:
  for (col = tensorPositionStart; col < tensorPositionStop; col += 2)
  {
// b)1. get norm-MPO:
   this->getBulkNormMPO(NormMPS1, NormMPS2, NormTensors, col, NormMPO);
// b)2. compute expectation value:
// b)2. - horizontal-even:
   if (tensorPositionStart == 0)
   {
    if (col == 0)
     HMPO = HMPOs[1];
    else if ((col > 0) && (col < Ncols-2))
     HMPO = HMPOs[0];
    else if (col == Ncols-2)
     HMPO = HMPOs[2];
   }
// b)2. - horizontal-odd:
   else if (tensorPositionStart == 1)
   {
    HMPO = HMPOs[3];
   }
   result += this->expectationValue(PEPS0, row, col, Direction, HMPO, NormMPO);
// b)3. update norm-tensors:
   this->updateBulkNormTensors(PEPS0, Direction, position, NormMPS1, NormMPS2, col, NormTensors);
  }
 }
 return result;
}

template<class T> T TimeEvolution2D<T>::expectationValue(const PEPS<T>& PEPS0, unsigned int positionRow, unsigned int positionCol, const string& Direction,
                                                         const MPO<T>& HMPO, const MPO<T>& NormMPO) const
{
 T exp, norm, result;
 vector<unsigned int> Indices01(1), Indices11(1), Indices02(2), Indices12(2), Indices03(3), Indices13(3), Indices04(4), Indices14(4);
 vector<unsigned int> Indices2(2), Indices3(3), Order5(5), Shape3(3), Shape4(4);
 Tensor<T> NormTensor, Tensor0, Tensor0Red, Tensor1, Tensor1Red, TensorL, TensorPair, TensorR;
// 1. get tensor pair from PEPS0:
 PEPS0.get(positionRow, positionCol, Tensor0);
 if (Direction == "vertical")
 {
  PEPS0.get(positionRow+1, positionCol, Tensor1);
// 1. - permute tensor pair to conform it to the order of NormMPO:
  Order5[0] = 1; Order5[1] = 2; Order5[2] = 3; Order5[3] = 0; Order5[4] = 4;
  Tensor0.permute(Order5);
  Tensor1.permute(Order5);
 }
 else if (Direction == "horizontal")
 {
  PEPS0.get(positionRow, positionCol+1, Tensor1);
 }
// 2. split tensors and compute reduced tensors Tensor0Red and Tensor1Red:
 Indices3[0] = 0; Indices3[1] = 1; Indices3[2] = 3;
 Indices2[0] = 2; Indices2[1] = 4;
 Tensor0.QRDecompose(Indices3, Indices2, Tensor0Red);
 TensorL = Tensor0;
 Indices2[0] = 0; Indices2[1] = 4;
 Indices3[0] = 1; Indices3[1] = 2; Indices3[2] = 3;
 Tensor1.LQDecompose(Indices2, Indices3, Tensor1Red);
 TensorR = Tensor1;
// 3. compute NormTensor:
// 3.1. left half:
 NormMPO.get(5, NormTensor); NormMPO.get(0, Tensor1);
 Indices01[0] = 1; Indices11[0] = 0;
 NormTensor.contract(Indices01, Tensor1, Indices11);
 Tensor1 = TensorL;
 Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 2; Indices12[1] = 0;
 NormTensor.contract(Indices02, Tensor1, Indices12);
 TensorL.complexConjugate(Tensor1);
 Indices02[1] = 3;
 NormTensor.contract(Indices02, Tensor1, Indices12);
 NormMPO.get(1, Tensor1);
 Indices03[0] = 1; Indices03[1] = 2; Indices03[2] = 4; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
 NormTensor.contract(Indices03, Tensor1, Indices13);
// 3.2. right half:
 NormMPO.get(2, Tensor0); NormMPO.get(3, Tensor1);
 Indices01[0] = 1; Indices11[0] = 0;
 Tensor0.contract(Indices01, Tensor1, Indices11);
 Tensor1 = TensorR;
 Indices02[0] = 1; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 TensorR.complexConjugate(Tensor1);
 Indices02[1] = 3;
 Tensor0.contract(Indices02, Tensor1, Indices12);
 NormMPO.get(4, Tensor1);
 Indices03[0] = 1; Indices03[1] = 3; Indices03[2] = 5; Indices13[0] = 0; Indices13[1] = 2; Indices13[2] = 3;
 Tensor0.contract(Indices03, Tensor1, Indices13);
// 3.3. contract left with right half:
 Indices02[0] = 0; Indices02[1] = 3; Indices12[0] = 3; Indices12[1] = 0;
 NormTensor.contract(Indices02, Tensor0, Indices12);
// 4. compute norm:
 TensorPair = Tensor0Red; Tensor1 = Tensor1Red;
 Indices01[0] = 1; Indices11[0] = 0;
 TensorPair.contract(Indices01, Tensor1, Indices11);
 Tensor1 = TensorPair;
 Indices02[0] = 0; Indices02[1] = 2; Indices12[0] = 0; Indices12[1] = 3;
 NormTensor.contract(Indices02, Tensor1, Indices12);
 Tensor0 = NormTensor; TensorPair.complexConjugate(Tensor1);
 Indices04[0] = 0; Indices04[1] = 1; Indices04[2] = 2; Indices04[3] = 3; Indices14[0] = 0; Indices14[1] = 3; Indices14[2] = 1; Indices14[3] = 2;
 Tensor0.contract(Indices04, Tensor1, Indices14);
 norm = Tensor0.get(0);
// 5. compute expectation value:
 HMPO.get(0, Tensor0); HMPO.get(1, Tensor1);
 Tensor0.getShape(Shape4);
 Shape3[0] = Shape4[1]; Shape3[1] = Shape4[2]; Shape3[2] = Shape4[3];
 Tensor0.reshape(Shape3);
 Tensor1.getShape(Shape4);
 Shape3[0] = Shape4[0]; Shape3[1] = Shape4[2]; Shape3[2] = Shape4[3];
 Tensor1.reshape(Shape3);
 Indices01[0] = 0; Indices11[0] = 0;
 Tensor0.contract(Indices01, Tensor1, Indices11);
 Indices02[0] = 2; Indices02[1] = 3; Indices12[0] = 0; Indices12[1] = 2;
 NormTensor.contract(Indices02, Tensor0, Indices12);
 TensorPair.complexConjugate();
 NormTensor.contract(Indices04, TensorPair, Indices14);
 exp = NormTensor.get(0);
// 6. return result:
 result = exp/norm;
 return result;
}
