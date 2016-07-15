/// Singleton template class GroundState2D implements ground state computation for PEPS.
/** The singleton template class GroundState2D implements all functions for ground state computation with
    PEPS.
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

template<class T> class GroundState2D
{
 public:

/// Constructs Singleton GroundState2D<T> and returns reference.
  static GroundState2D<T>& reference()
  {
   static GroundState2D<T> GS;
   return GS;
  }

/// Computes ground state.
/** This function computes the ground state of a given Hamiltonian2D H0 by directly minimizing
       <PEPS0|H0|PEPS0>/<PEPS0|PEPS0>
    updating one tensor after the other with
       void GroundState2D<T>::updateTensor(Tensor<T>& NormTensor, Tensor<T>& H0Tensor,
                                           Tensor<T>& H0CenterTensor, Tensor<T>& Tensor0,
                                           double cutoff = 1.0e-4, unsigned int mode = 0) const   .
    It starts with PEPS0 and sweeps through the tensors going from one column to the next, each time
    minimizing the energy. The number of local sweeps per column are given by numSweepsLocal.
    The energy is computed after each global sweep over all columns and if the relative change in energy from one
    global sweep to the next is below eps, the function returns. The function always returns if the number
    of global sweeps done exceeds maxNumSweeps.
    The environment and Hamiltonian MPSs of the updated column are computed with virtual bond dimension
    D2Env, convergence precision epsEnv and maximal numbers of sweeps maxNumSweepsEnv using
       void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps,
                           double& errorAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1)   .
    \param H0 input: const Hamiltonian2D<T>&, the Hamiltonian
    \param D2Env input: unsigned int, the maximal virtual bond dimension of the environment MPSs,
                        must fulfill (D2Env > 0)
    \param epsEnv input: double, the convergence precision for the environment MPSs,
                         must fulfill (epsEnv >= 0.0)
    \param maxNumSweepsEnv input: unsigned int, the maximal numbers of sweeps allowed in
                                  the environment MPS approximation
    \param numSweepsLocal input: unsigned int, the number of local sweeps per column
    \param eps input: double, the global convergence precision,
                      i.e. convergence in global sweep n if eps >= abs(energy(n) - energy(n-1))/abs(energy(n)),
                      must fulfill (eps >= 0.0)
    \param maxNumSweeps input: unsigned int, the maximal number of global sweeps allowed
    \param Energies output: vector<double>&, the energies obtained after each global sweep,
                            must fulfill (Energies.size() == maxNumSweeps+1)
    \param numSweepsDone output: unsigned int&, the final number of global sweeps done
    \param PEPS0 input/output: PEPS<T>&, on input the initial PEPS, on output the final PEPS,
                               must have the correct form
    \param cutoff optional input: double, the cutoff used in mode 0 of void GroundState2D<T>::updateTensor(...),
                                  must fulfill (cutoff >= 0.0), has default value 1.0e-4
    \param mode optional input: unsigned int, the tensor update mode used in
                                void GroundState2D<T>::updateTensor(...), has default value 0
    \sa void GroundState2D<T>::updateTensor(Tensor<T>& NormTensor, Tensor<T>& H0Tensor,
                                            Tensor<T>& H0CenterTensor, Tensor<T>& Tensor0,
                                            double cutoff = 1.0e-4, unsigned int mode = 0) const
    \sa void multiplyMPOMPS(const MPO<T>& MPO0, const MPS<T>& MPS0, double eps, unsigned int maxNumSweeps,
                            double& errorAchieved, unsigned int& numSweepsDone, MPS<T>& MPS1) */
  void computeGroundState(const Hamiltonian2D<T>& H0,
                          unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv,
                          unsigned int numSweepsLocal,
                          double eps, unsigned int maxNumSweeps,
                          vector<double>& Energies, unsigned int& numSweepsDone,
                          PEPS<T>& PEPS0,
                          double cutoff = 1.0e-4, unsigned int mode = 0) const;

/// Computes ground state exactly.
/** This function computes the ground state of a given Hamiltonian2D H0 exactly by using its matrix
    representation and diagonalizing this matrix with LAPACK's routines for hermitian matrices. The
    underlying diagonalisation routine is a QR factorization.
    \param H0 input: const Hamiltonian2D<T>&, the Hamiltonian
    \param energy output: double&, the ground state energy
    \param State output: vector<T>&, the ground state */
  void computeGroundStateExactly(const Hamiltonian2D<T>& H0, double& energy, vector<T>& State) const;

 private:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  GroundState2D();

/// Standard copy constructor. Not implemented for this Singleton.
/** The standard copy constructor copies the input GroundState2D into this. This constructor must not
    be implemented for this Singleton class.
    \param GroundState0 input: const GroundState2D<T>&, to be copied into this
    \sa GroundState2D<T>& operator=(const GroundState2D<T>& GroundState0) */
  GroundState2D(const GroundState2D<T>& GroundState0);

/// Assigns GroundState2D to this. Not implemented for this Singleton.
/** The operator= allows to assign a GroundState2D to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side GroundState2D. This function must not be
    implemented for this Singleton class.
    \param GroundState0 input: const GroundState2D<T>&, to be copied into this
    \return GroundState2D<T>&, a reference to the new this
    \sa GroundState2D(const GroundState2D<T>& GroundState0) */
  GroundState2D<T>& operator=(const GroundState2D<T>& GroundState0);

/// Updates tensor.
/** This function updates a tensor Tensor0 given an environment NormTensor, a Hamiltonian H0Tensor and
    a Hamiltonian center tensor H0CenterTensor.
    It finds the tensor that minimizes
       <Tensor0|H0Tensor|Tensor0>/<Tensor0|NormTensor|Tensor0>
    and is used by the function
       void GroundState2D<T>::computeGroundState(const Hamiltonian2D<T>& H0,
                                                 unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv,
                                                 unsigned int numSweepsLocal,
                                                 double eps, unsigned int maxNumSweeps,
                                                 vector<double>& Energies, unsigned int& numSweepsDone,
                                                 PEPS<T>& PEPS0,
                                                 double cutoff = 1.0e-4, unsigned int mode = 0) const   .
    In mode 0, the norm matrix is made positive definite and then all its eigenvalues smaller than cutoff are
    discarded. LAPACK's XSYEV and XHEEV solve the standard hermitian eigenproblem. This is the default mode.
    In mode 1, LAPACK's XGGEVX solve the generalized nonsymmetric eigenproblem.
    In mode 2, a steepest descent method starts from Tensor0 and successively minimizes the energy.
    \param NormTensor input: Tensor<T>&, the norm tensor, must have the correct shape
    \param H0Tensor input: Tensor<T>&, the Hamiltonian tensor, must have the correct shape
    \param H0CenterTensor input: Tensor<T>&, the Hamiltonian center tensor, must have the correct shape
    \param Tensor0 input/output: Tensor<T>&, on input the initial tensor, on output the final tensor,
                                 must have the correct shape
    \param cutoff optional input: double, the cutoff used in mode 0, must fulfill (cutoff >= 0.0),
                                  has default value 1.0e-4
    \param mode optional input: unsigned int, the mode, if (mode == 0) then a standard hermitian eigenproblem
                                is solved, if (mode == 1) then a generalized nonsymmetric eigenproblem is
                                solved, if (mode == 2) then a steepest descent method is applied,
                                has default value 0
    \sa void GroundState2D<T>::computeGroundState(const Hamiltonian2D<T>& H0,
                                                  unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv,
                                                  unsigned int numSweepsLocal,
                                                  double eps, unsigned int maxNumSweeps,
                                                  vector<double>& Energies, unsigned int& numSweepsDone,
                                                  PEPS<T>& PEPS0,
                                                  double cutoff = 1.0e-4, unsigned int mode = 0) const */
  void updateTensor(Tensor<T>& NormTensor, Tensor<T>& H0Tensor, Tensor<T>& H0CenterTensor,
                    Tensor<T>& Tensor0, double cutoff = 1.0e-4, unsigned int mode = 0) const;
};

template<class T> void GroundState2D<T>::computeGroundState(const Hamiltonian2D<T>& H0,
                                          unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv,
                                          unsigned int numSweepsLocal,
                                          double eps, unsigned int maxNumSweeps,
                                          vector<double>& Energies, unsigned int& numSweepsDone,
                                          PEPS<T>& PEPS0,
                                          double cutoff, unsigned int mode) const
{
 string BC; H0.getBC(BC);
 unsigned int Nrows = H0.getNrows(), Ncols = H0.getNcols();
 Matrix<unsigned int> d(Nrows, Ncols); H0.getd(d);
#ifdef DEBUG
 string BCPEPS0; PEPS0.getBC(BCPEPS0);
 Matrix<unsigned int> dPEPS0; PEPS0.getd(dPEPS0);
 if ((Nrows == 0) || (D2Env == 0) || (epsEnv < 0.0) || (eps < 0.0) || (Energies.size() != maxNumSweeps+1) ||
     (BC != BCPEPS0) || (Nrows != PEPS0.getNrows()) || (Ncols != PEPS0.getNcols()) || (d != dPEPS0) ||
     (cutoff < 0.0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void GroundState2D<T>::" <<
          "computeGroundState(const Hamiltonian2D<T>& H0, " <<
                             "unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv, " <<
                             "unsigned int numSweepsLocal, " <<
                             "double eps, unsigned int maxNumSweeps, " <<
                             "vector<double>& Energies, unsigned int& numSweepsDone, " <<
                             "PEPS<T>& PEPS0, " <<
                             "double cutoff = 1.0e-4, unsigned int mode = 0) const: " <<
          "((Nrows == 0) || (D2Env == 0) || (epsEnv < 0.0) || (eps < 0.0) || " <<
           "(Energies.size() != maxNumSweeps+1) || (BC != BCPEPS0) || (Nrows != PEPS0.getNrows()) || " <<
           "(Ncols != PEPS0.getNcols()) || (d != dPEPS0) || (cutoff < 0.0))." << endl;
  exit(1);
 }
#endif
 MPS<T> MPS0, MPS1; MPO<T> MPO0;
 unsigned int D0, D1, D, numSweepsDoneEnv, numSweepsDoneLocal, row, col;
 double energy, errorAchievedEnv, normSquaredPEPS0, expectationValueH0, newEnergy;
 string RepresentationH0, Direction, DirectionLocal;
 Tensor<T> Tensor0, Tensor1, Tensor2, NormTensor, H0Tensor, H0CenterTensor;
 vector<unsigned int> Indices01(1), Indices11(1), Indices02(2), Indices12(2);
 D = PEPS0.getD();
 H0.getRepresentation(RepresentationH0);
 energy = 1.0e300;
 numSweepsDone = 0;
 if (RepresentationH0 == "PEPO")
 {
// Get the PEPO representation of H0:
  PEPO<T> PEPOH0; H0.getPEPO(PEPOH0);
  unsigned int DH0 = PEPOH0.getD();
  vector< MPS<T> > NormMPSs(Ncols), H0MPSs(Ncols);
  vector< Tensor<T> > NormTensors(Nrows), H0Tensors(Nrows);
  vector<unsigned int> Shape3(3), Shape4(4), Order4(4), Order6(6), Order8(8), Order12(12), Order28(8);
  Order8[0] = 0; Order8[1] = 4; Order8[2] = 1; Order8[3] = 5;
  Order8[4] = 2; Order8[5] = 6; Order8[6] = 3; Order8[7] = 7;
  Order12[0] = 0; Order12[1] = 4; Order12[2] = 8; Order12[3] = 1;
  Order12[4] = 5; Order12[5] = 9; Order12[6] = 2; Order12[7] = 6;
  Order12[8] = 10; Order12[9] = 3; Order12[10] = 7; Order12[11] = 11;
  Order28[0] = 1; Order28[1] = 7; Order28[2] = 5; Order28[3] = 0;
  Order28[4] = 2; Order28[5] = 3; Order28[6] = 4; Order28[7] = 6;
  if (BC == "open")
  {
// Compute the norm MPSs for the norm to each column:
   PEPS0.getMPS("left", MPS1);
   NormMPSs[Ncols-1] = MPS1;
   for (col = Ncols-2; col > 0; col--)
   {
    MPS0 = MPS1;
    PEPS0.getMPO("left", col, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2Env);
    MPS1.setD(D1);
    multiplyMPOMPS(MPO0, MPS0, epsEnv, maxNumSweepsEnv, errorAchievedEnv, numSweepsDoneEnv, MPS1);
    NormMPSs[col] = MPS1;
   }
// Compute the Hamiltonian MPSs for H0 to each column:
   PEPS0.getMPS(PEPOH0, "left", MPS1);
   H0MPSs[Ncols-1] = MPS1;
   for (col = Ncols-2; col > 0; col--)
   {
    MPS0 = MPS1;
    PEPS0.getMPO(PEPOH0, "left", col, MPO0);
    D1 = min(MPS0.getD()*MPO0.getD(), D2Env);
    MPS1.setD(D1);
    multiplyMPOMPS(MPO0, MPS0, epsEnv, maxNumSweepsEnv, errorAchievedEnv, numSweepsDoneEnv, MPS1);
    H0MPSs[col] = MPS1;
   }
// sweep:
   string Direction = "right";
   while (numSweepsDone < maxNumSweeps)
   {
    if (Direction == "right")
    {
// compute energy:
     PEPS0.getMPS(Direction, MPS0);
     MPS1 = NormMPSs[1];
     normSquaredPEPS0 = MPS0.contractReverse(MPS1);
     PEPS0.getMPS(PEPOH0, Direction, MPS0);
     MPS1 = H0MPSs[1];
     expectationValueH0 = MPS0.contractReverse(MPS1);
     newEnergy = expectationValueH0 / normSquaredPEPS0;
     if (abs(newEnergy - energy)/abs(newEnergy) <= eps)
     {
      Energies[numSweepsDone] = newEnergy;
      return;
     }
     energy = newEnergy;
     Energies[numSweepsDone] = newEnergy;
// sweep from left to right:
     for (col = 0; col < Ncols-1; col++)
     {
// sweep from left to right starting with column 0:
      if (col == 0)
      {
// compute norm tensors:
       MPS1 = NormMPSs[col+1];
       PEPS0.getMPS(Direction, MPS0);
       MPS0.get(0, Tensor0); MPS1.get(Nrows-1, Tensor1);
       Indices01[0] = 2; Indices11[0] = 2;
       Tensor0.contract(Indices01, Tensor1, Indices11);
       Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
       Tensor0.permute(Order4);
       NormTensors[Nrows-1] = Tensor0;
       Indices01[0] = 2; Indices11[0] = 0;
       Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
       for (row = Nrows-2; row > 0; row--)
       {
        MPS0.get(Nrows-row-1, Tensor1);
        Tensor0.contract(Indices01, Tensor1, Indices11);
        MPS1.get(row, Tensor1);
        Tensor0.contract(Indices02, Tensor1, Indices12);
        NormTensors[row] = Tensor0;
       }
// compute H0 tensors:
       MPS1 = H0MPSs[col+1];
       PEPS0.getMPS(PEPOH0, Direction, MPS0);
       MPS0.get(0, Tensor0); MPS1.get(Nrows-1, Tensor1);
       Indices01[0] = 2; Indices11[0] = 2;
       Tensor0.contract(Indices01, Tensor1, Indices11);
       Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
       Tensor0.permute(Order4);
       H0Tensors[Nrows-1] = Tensor0;
       Indices01[0] = 2; Indices11[0] = 0;
       Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
       for (row = Nrows-2; row > 0; row--)
       {
        MPS0.get(Nrows-row-1, Tensor1);
        Tensor0.contract(Indices01, Tensor1, Indices11);
        MPS1.get(row, Tensor1);
        Tensor0.contract(Indices02, Tensor1, Indices12);
        H0Tensors[row] = Tensor0;
       }
       numSweepsDoneLocal = 0;
       DirectionLocal = "down";
       while (numSweepsDoneLocal < numSweepsLocal)
       {
        if (DirectionLocal == "down")
        {
// sweep from up to down starting with row 0:
         row = 0;
// update tensor:
         NormMPSs[col+1].get(row, Tensor0); Tensor1 = NormTensors[row+1];
         Indices01[0] = 1; Indices11[0] = 3;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         NormTensor = Tensor0;
         H0MPSs[col+1].get(row, Tensor0); Tensor1 = H0Tensors[row+1];
         Tensor0.contract(Indices01, Tensor1, Indices11);
         H0Tensor = Tensor0;
         PEPOH0.get(row, col, H0CenterTensor);
         PEPS0.get(row, col, Tensor0);
         this->updateTensor(NormTensor, H0Tensor, H0CenterTensor, Tensor0, cutoff, mode);
         PEPS0.set(row, col, Tensor0);
// compute norm tensor:
         Tensor0.complexConjugate(Tensor1);
         Indices01[0] = 4; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order8);
         Shape3[0] = 1; Shape3[1] = D*D; Shape3[2] = D*D;
         Tensor0.reshape(Shape3);
         NormMPSs[col+1].get(row, Tensor1);
         Indices01[0] = 1; Indices11[0] = 2;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Order4[0] = 0; Order4[1] = 2; Order4[2] = 1; Order4[3] = 3;
         Tensor0.permute(Order4);
         NormTensors[row] = Tensor0;
// compute H0 tensor:
         PEPS0.get(row, col, Tensor0); PEPOH0.get(row, col, Tensor1);
         Indices01[0] = 4; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         PEPS0.get(row, col, Tensor1); Tensor1.complexConjugate();
         Indices01[0] = 8; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order12);
         Shape3[0] = 1; Shape3[1] = D*DH0*D; Shape3[2] = D*DH0*D;
         Tensor0.reshape(Shape3);
         H0MPSs[col+1].get(row, Tensor1);
         Indices01[0] = 1; Indices11[0] = 2;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Order4[0] = 0; Order4[1] = 2; Order4[2] = 1; Order4[3] = 3;
         Tensor0.permute(Order4);
         H0Tensors[row] = Tensor0;
// sweep from up to down from row 1 to Nrows-2:
         for (row = 1; row < Nrows-1; row++)
         {
// update tensor:
          Tensor0 = NormTensors[row-1]; NormMPSs[col+1].get(row, Tensor1);
          Indices01[0] = 3; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor1 = NormTensors[row+1];
          Indices01[0] = 3; Indices11[0] = 3;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          NormTensor = Tensor0;
          Tensor0 = H0Tensors[row-1]; H0MPSs[col+1].get(row, Tensor1);
          Indices01[0] = 3; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor1 = H0Tensors[row+1];
          Indices01[0] = 3; Indices11[0] = 3;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          H0Tensor = Tensor0;
          PEPOH0.get(row, col, H0CenterTensor);
          PEPS0.get(row, col, Tensor0);
          this->updateTensor(NormTensor, H0Tensor, H0CenterTensor, Tensor0, cutoff, mode);
          PEPS0.set(row, col, Tensor0);
// compute norm tensor:
          Tensor0.complexConjugate(Tensor1);
          Indices01[0] = 4; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor0.permute(Order8);
          Shape3[0] = D*D; Shape3[1] = D*D; Shape3[2] = D*D;
          Tensor0.reshape(Shape3);
          Tensor1 = Tensor0; Tensor0 = NormTensors[row-1];
          Indices01[0] = 2; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          NormMPSs[col+1].get(row, Tensor1);
          Indices02[0] = 2; Indices02[1] = 3; Indices12[0] = 0; Indices12[1] = 2;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          NormTensors[row] = Tensor0;
// compute H0 tensor:
          PEPS0.get(row, col, Tensor0); PEPOH0.get(row, col, Tensor1);
          Indices01[0] = 4; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          PEPS0.get(row, col, Tensor1); Tensor1.complexConjugate();
          Indices01[0] = 8; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor0.permute(Order12);
          Shape3[0] = D*DH0*D; Shape3[1] = D*DH0*D; Shape3[2] = D*DH0*D;
          Tensor0.reshape(Shape3);
          Tensor1 = Tensor0; Tensor0 = H0Tensors[row-1];
          Indices01[0] = 2; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          H0MPSs[col+1].get(row, Tensor1);
          Indices02[0] = 2; Indices02[1] = 3; Indices12[0] = 0; Indices12[1] = 2;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          H0Tensors[row] = Tensor0;
         }
         DirectionLocal = "up";
        }
        else if (DirectionLocal == "up")
        {
// sweep from down to up starting with row Nrows-1:
         row = Nrows-1;
// update tensor:
         Tensor0 = NormTensors[row-1]; NormMPSs[col+1].get(row, Tensor1);
         Indices01[0] = 3; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         NormTensor = Tensor0;
         Tensor0 = H0Tensors[row-1]; H0MPSs[col+1].get(row, Tensor1);
         Tensor0.contract(Indices01, Tensor1, Indices11);
         H0Tensor = Tensor0;
         PEPOH0.get(row, col, H0CenterTensor);
         PEPS0.get(row, col, Tensor0);
         this->updateTensor(NormTensor, H0Tensor, H0CenterTensor, Tensor0, cutoff, mode);
         PEPS0.set(row, col, Tensor0);
// compute norm tensor:
         Tensor0.complexConjugate(Tensor1);
         Indices01[0] = 4; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order8);
         Shape3[0] = D*D; Shape3[1] = D*D; Shape3[2] = 1;
         Tensor0.reshape(Shape3);
         NormMPSs[col+1].get(row, Tensor1);
         Indices01[0] = 1; Indices11[0] = 2;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Order4[0] = 1; Order4[1] = 3; Order4[2] = 0; Order4[3] = 2;
         Tensor0.permute(Order4);
         NormTensors[row] = Tensor0;
// compute H0 tensor:
         PEPS0.get(row, col, Tensor0); PEPOH0.get(row, col, Tensor1);
         Indices01[0] = 4; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         PEPS0.get(row, col, Tensor1); Tensor1.complexConjugate();
         Indices01[0] = 8; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order12);
         Shape3[0] = D*DH0*D; Shape3[1] = D*DH0*D; Shape3[2] = 1;
         Tensor0.reshape(Shape3);
         H0MPSs[col+1].get(row, Tensor1);
         Indices01[0] = 1; Indices11[0] = 2;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Order4[0] = 1; Order4[1] = 3; Order4[2] = 0; Order4[3] = 2;
         Tensor0.permute(Order4);
         H0Tensors[row] = Tensor0;
// sweep from down to up from row Nrows-2 to 1:
         for (row = Nrows-2; row > 0; row--)
         {
// update tensor:
          Tensor0 = NormTensors[row-1]; NormMPSs[col+1].get(row, Tensor1);
          Indices01[0] = 3; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor1 = NormTensors[row+1];
          Indices01[0] = 3; Indices11[0] = 3;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          NormTensor = Tensor0;
          Tensor0 = H0Tensors[row-1]; H0MPSs[col+1].get(row, Tensor1);
          Indices01[0] = 3; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor1 = H0Tensors[row+1];
          Indices01[0] = 3; Indices11[0] = 3;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          H0Tensor = Tensor0;
          PEPOH0.get(row, col, H0CenterTensor);
          PEPS0.get(row, col, Tensor0);
          this->updateTensor(NormTensor, H0Tensor, H0CenterTensor, Tensor0, cutoff, mode);
          PEPS0.set(row, col, Tensor0);
// compute norm tensor:
          Tensor0.complexConjugate(Tensor1);
          Indices01[0] = 4; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor0.permute(Order8);
          Shape3[0] = D*D; Shape3[1] = D*D; Shape3[2] = D*D;
          Tensor0.reshape(Shape3);
          Tensor1 = Tensor0; Tensor0 = NormTensors[row+1];
          Indices01[0] = 2; Indices11[0] = 2;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          NormMPSs[col+1].get(row, Tensor1);
          Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          NormTensors[row] = Tensor0;
// compute H0 tensor:
          PEPS0.get(row, col, Tensor0); PEPOH0.get(row, col, Tensor1);
          Indices01[0] = 4; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          PEPS0.get(row, col, Tensor1); Tensor1.complexConjugate();
          Indices01[0] = 8; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor0.permute(Order12);
          Shape3[0] = D*DH0*D; Shape3[1] = D*DH0*D; Shape3[2] = D*DH0*D;
          Tensor0.reshape(Shape3);
          Tensor1 = Tensor0; Tensor0 = H0Tensors[row+1];
          Indices01[0] = 2; Indices11[0] = 2;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          H0MPSs[col+1].get(row, Tensor1);
          Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          H0Tensors[row] = Tensor0;
         }
         DirectionLocal = "down";
        }
        numSweepsDoneLocal++;
       }
       PEPS0.getMPS(Direction, MPS0);
       NormMPSs[col] = MPS0;
       PEPS0.getMPS(PEPOH0, Direction, MPS0);
       H0MPSs[col] = MPS0;
      }
// sweep from left to right from column 1 to Ncols-2:
      else
      {
// compute norm tensors:
       MPS0 = NormMPSs[col-1]; PEPS0.getMPO(Direction, col, MPO0); MPS1 = NormMPSs[col+1];
       MPS0.get(0, Tensor0); MPO0.get(0, Tensor1);
       Indices01[0] = 2; Indices11[0] = 2;
       Tensor0.contract(Indices01, Tensor1, Indices11);
       MPS1.get(Nrows-1, Tensor1);
       Indices01[0] = 4; Indices11[0] = 2;
       Tensor0.contract(Indices01, Tensor1, Indices11);
       Order6[0] = 0; Order6[1] = 2; Order6[2] = 5; Order6[3] = 1; Order6[4] = 3; Order6[5] = 4;
       Tensor0.permute(Order6);
       NormTensors[Nrows-1] = Tensor0;
       Indices01[0] = 3; Indices11[0] = 0;
       Indices02[0] = 3; Indices02[1] = 6;
       for (row = Nrows-2; row > 0; row--)
       {
        MPS0.get(Nrows-row-1, Tensor1);
        Tensor0.contract(Indices01, Tensor1, Indices11);
        MPO0.get(Nrows-row-1, Tensor1);
        Indices12[0] = 0; Indices12[1] = 2;
        Tensor0.contract(Indices02, Tensor1, Indices12);
        MPS1.get(row, Tensor1);
        Indices12[0] = 1; Indices12[1] = 2;
        Tensor0.contract(Indices02, Tensor1, Indices12);
        NormTensors[row] = Tensor0;
       }
// compute H0 tensors:
       MPS0 = H0MPSs[col-1]; PEPS0.getMPO(PEPOH0, Direction, col, MPO0); MPS1 = H0MPSs[col+1];
       MPS0.get(0, Tensor0); MPO0.get(0, Tensor1);
       Indices01[0] = 2; Indices11[0] = 2;
       Tensor0.contract(Indices01, Tensor1, Indices11);
       MPS1.get(Nrows-1, Tensor1);
       Indices01[0] = 4; Indices11[0] = 2;
       Tensor0.contract(Indices01, Tensor1, Indices11);
       Tensor0.permute(Order6);
       H0Tensors[Nrows-1] = Tensor0;
       Indices01[0] = 3; Indices11[0] = 0;
       Indices02[0] = 3; Indices02[1] = 6;
       for (row = Nrows-2; row > 0; row--)
       {
        MPS0.get(Nrows-row-1, Tensor1);
        Tensor0.contract(Indices01, Tensor1, Indices11);
        MPO0.get(Nrows-row-1, Tensor1);
        Indices12[0] = 0; Indices12[1] = 2;
        Tensor0.contract(Indices02, Tensor1, Indices12);
        MPS1.get(row, Tensor1);
        Indices12[0] = 1; Indices12[1] = 2;
        Tensor0.contract(Indices02, Tensor1, Indices12);
        H0Tensors[row] = Tensor0;
       }
       numSweepsDoneLocal = 0;
       DirectionLocal = "down";
       while (numSweepsDoneLocal < numSweepsLocal)
       {
        if (DirectionLocal == "down")
        {
// sweep from up to down starting with row 0:
         row = 0;
// update tensor:
         NormMPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = NormTensors[row+1];
         Indices01[0] = 0; Indices11[0] = 3;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         NormMPSs[col+1].get(row, Tensor1);
         Indices01[0] = 6; Indices11[0] = 1;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order28);
         NormTensor = Tensor0;
         H0MPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = H0Tensors[row+1];
         Indices01[0] = 0; Indices11[0] = 3;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         H0MPSs[col+1].get(row, Tensor1);
         Indices01[0] = 6; Indices11[0] = 1;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order28);
         H0Tensor = Tensor0;
         PEPOH0.get(row, col, H0CenterTensor);
         PEPS0.get(row, col, Tensor0);
         this->updateTensor(NormTensor, H0Tensor, H0CenterTensor, Tensor0, cutoff, mode);
         PEPS0.set(row, col, Tensor0);
// compute norm tensor:
         Tensor0.complexConjugate(Tensor1);
         Indices01[0] = 4; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order8);
         Shape4[0] = D*D; Shape4[1] = 1; Shape4[2] = D*D; Shape4[3] = D*D;
         Tensor0.reshape(Shape4);
         Tensor1 = Tensor0; NormMPSs[col-1].get(Nrows-row-1, Tensor0);
         Indices01[0] = 2; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         NormMPSs[col+1].get(row, Tensor1);
         Indices01[0] = 3; Indices11[0] = 2;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Order6[0] = 1; Order6[1] = 2; Order6[2] = 4; Order6[3] = 0; Order6[4] = 3; Order6[5] = 5;
         Tensor0.permute(Order6);
         NormTensors[row] = Tensor0;
// compute H0 tensor:
         PEPS0.get(row, col, Tensor0); PEPOH0.get(row, col, Tensor1);
         Indices01[0] = 4; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         PEPS0.get(row, col, Tensor1); Tensor1.complexConjugate();
         Indices01[0] = 8; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order12);
         Shape4[0] = D*DH0*D; Shape4[1] = 1; Shape4[2] = D*DH0*D; Shape4[3] = D*DH0*D;
         Tensor0.reshape(Shape4);
         Tensor1 = Tensor0; H0MPSs[col-1].get(Nrows-row-1, Tensor0);
         Indices01[0] = 2; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         H0MPSs[col+1].get(row, Tensor1);
         Indices01[0] = 3; Indices11[0] = 2;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order6);
         H0Tensors[row] = Tensor0;
// sweep from up to down from row 1 to Nrows-2:
         for (row = 1; row < Nrows-1; row++)
         {
// update tensor:
          NormMPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = NormTensors[row-1];
          Indices01[0] = 1; Indices11[0] = 3;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          NormMPSs[col+1].get(row, Tensor1);
          Indices01[0] = 6; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor1 = NormTensors[row+1];
          Indices02[0] = 0; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 5;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          NormTensor = Tensor0;
          H0MPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = H0Tensors[row-1];
          Indices01[0] = 1; Indices11[0] = 3;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          H0MPSs[col+1].get(row, Tensor1);
          Indices01[0] = 6; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor1 = H0Tensors[row+1];
          Indices02[0] = 0; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 5;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          H0Tensor = Tensor0;
          PEPOH0.get(row, col, H0CenterTensor);
          PEPS0.get(row, col, Tensor0);
          this->updateTensor(NormTensor, H0Tensor, H0CenterTensor, Tensor0, cutoff, mode);
          PEPS0.set(row, col, Tensor0);
// compute norm tensor:
          Tensor0.complexConjugate(Tensor1);
          Indices01[0] = 4; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor0.permute(Order8);
          Shape4[0] = D*D; Shape4[1] = D*D; Shape4[2] = D*D; Shape4[3] = D*D;
          Tensor0.reshape(Shape4);
          Tensor2 = Tensor0;
          Tensor0 = NormTensors[row-1]; NormMPSs[col-1].get(Nrows-row-1, Tensor1);
          Indices01[0] = 3; Indices11[0] = 1;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 0;
          Tensor0.contract(Indices02, Tensor2, Indices12);
          NormMPSs[col+1].get(row, Tensor1);
          Indices02[0] = 3; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 2;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          NormTensors[row] = Tensor0;
// compute H0 tensor:
          PEPS0.get(row, col, Tensor0); PEPOH0.get(row, col, Tensor1);
          Indices01[0] = 4; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          PEPS0.get(row, col, Tensor1); Tensor1.complexConjugate();
          Indices01[0] = 8; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor0.permute(Order12);
          Shape4[0] = D*DH0*D; Shape4[1] = D*DH0*D; Shape4[2] = D*DH0*D; Shape4[3] = D*DH0*D;
          Tensor0.reshape(Shape4);
          Tensor2 = Tensor0;
          Tensor0 = H0Tensors[row-1]; H0MPSs[col-1].get(Nrows-row-1, Tensor1);
          Indices01[0] = 3; Indices11[0] = 1;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 0;
          Tensor0.contract(Indices02, Tensor2, Indices12);
          H0MPSs[col+1].get(row, Tensor1);
          Indices02[0] = 3; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 2;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          H0Tensors[row] = Tensor0;
         }
         DirectionLocal = "up";
        }
        else if (DirectionLocal == "up")
        {
// sweep from down to up starting with row Nrows-1:
         row = Nrows-1;
// update tensor:
         NormMPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = NormTensors[row-1];
         Indices01[0] = 1; Indices11[0] = 3;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         NormMPSs[col+1].get(row, Tensor1);
         Indices01[0] = 6; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         NormTensor = Tensor0;
         H0MPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = H0Tensors[row-1];
         Indices01[0] = 1; Indices11[0] = 3;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         H0MPSs[col+1].get(row, Tensor1);
         Indices01[0] = 6; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         H0Tensor = Tensor0;
         PEPOH0.get(row, col, H0CenterTensor);
         PEPS0.get(row, col, Tensor0);
         this->updateTensor(NormTensor, H0Tensor, H0CenterTensor, Tensor0, cutoff, mode);
         PEPS0.set(row, col, Tensor0);
// compute norm tensor:
         Tensor0.complexConjugate(Tensor1);
         Indices01[0] = 4; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order8);
         Shape4[0] = D*D; Shape4[1] = D*D; Shape4[2] = D*D; Shape4[3] = 1;
         Tensor0.reshape(Shape4);
         Tensor1 = Tensor0; NormMPSs[col-1].get(Nrows-row-1, Tensor0);
         Indices01[0] = 2; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         NormMPSs[col+1].get(row, Tensor1);
         Indices01[0] = 3; Indices11[0] = 2;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Order6[0] = 0; Order6[1] = 3; Order6[2] = 5; Order6[3] = 1; Order6[4] = 2; Order6[5] = 4;
         Tensor0.permute(Order6);
         NormTensors[row] = Tensor0;
// compute H0 tensor:
         PEPS0.get(row, col, Tensor0); PEPOH0.get(row, col, Tensor1);
         Indices01[0] = 4; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         PEPS0.get(row, col, Tensor1); Tensor1.complexConjugate();
         Indices01[0] = 8; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order12);
         Shape4[0] = D*DH0*D; Shape4[1] = D*DH0*D; Shape4[2] = D*DH0*D; Shape4[3] = 1;
         Tensor0.reshape(Shape4);
         Tensor1 = Tensor0; H0MPSs[col-1].get(Nrows-row-1, Tensor0);
         Indices01[0] = 2; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         H0MPSs[col+1].get(row, Tensor1);
         Indices01[0] = 3; Indices11[0] = 2;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order6);
         H0Tensors[row] = Tensor0;
// sweep from down to up from row Nrows-2 to 1:
         for (row = Nrows-2; row > 0; row--)
         {
// update tensor:
          NormMPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = NormTensors[row-1];
          Indices01[0] = 1; Indices11[0] = 3;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          NormMPSs[col+1].get(row, Tensor1);
          Indices01[0] = 6; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor1 = NormTensors[row+1];
          Indices02[0] = 0; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 5;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          NormTensor = Tensor0;
          H0MPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = H0Tensors[row-1];
          Indices01[0] = 1; Indices11[0] = 3;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          H0MPSs[col+1].get(row, Tensor1);
          Indices01[0] = 6; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor1 = H0Tensors[row+1];
          Indices02[0] = 0; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 5;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          H0Tensor = Tensor0;
          PEPOH0.get(row, col, H0CenterTensor);
          PEPS0.get(row, col, Tensor0);
          this->updateTensor(NormTensor, H0Tensor, H0CenterTensor, Tensor0, cutoff, mode);
          PEPS0.set(row, col, Tensor0);
// compute norm tensor:
          Tensor0.complexConjugate(Tensor1);
          Indices01[0] = 4; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor0.permute(Order8);
          Shape4[0] = D*D; Shape4[1] = D*D; Shape4[2] = D*D; Shape4[3] = D*D;
          Tensor0.reshape(Shape4);
          Tensor2 = Tensor0;
          Tensor0 = NormTensors[row+1]; NormMPSs[col-1].get(Nrows-row-1, Tensor1);
          Indices01[0] = 3; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 0;
          Tensor0.contract(Indices02, Tensor2, Indices12);
          NormMPSs[col+1].get(row, Tensor1);
          Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          NormTensors[row] = Tensor0;
// compute H0 tensor:
          PEPS0.get(row, col, Tensor0); PEPOH0.get(row, col, Tensor1);
          Indices01[0] = 4; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          PEPS0.get(row, col, Tensor1); Tensor1.complexConjugate();
          Indices01[0] = 8; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor0.permute(Order12);
          Shape4[0] = D*DH0*D; Shape4[1] = D*DH0*D; Shape4[2] = D*DH0*D; Shape4[3] = D*DH0*D;
          Tensor0.reshape(Shape4);
          Tensor2 = Tensor0;
          Tensor0 = H0Tensors[row+1]; H0MPSs[col-1].get(Nrows-row-1, Tensor1);
          Indices01[0] = 3; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 0;
          Tensor0.contract(Indices02, Tensor2, Indices12);
          H0MPSs[col+1].get(row, Tensor1);
          Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          H0Tensors[row] = Tensor0;
         }
         DirectionLocal = "down";
        }
        numSweepsDoneLocal++;
       }
       MPS0 = NormMPSs[col-1]; PEPS0.getMPO(Direction, col, MPO0); MPS1 = MPS0;
       D1 = min(MPS0.getD()*MPO0.getD(), D2Env);
       MPS1.setD(D1);
       multiplyMPOMPS(MPO0, MPS0, epsEnv, maxNumSweepsEnv, errorAchievedEnv, numSweepsDoneEnv, MPS1);
       NormMPSs[col] = MPS1;
       MPS0 = H0MPSs[col-1]; PEPS0.getMPO(PEPOH0, Direction, col, MPO0); MPS1 = MPS0;
       D1 = min(MPS0.getD()*MPO0.getD(), D2Env);
       MPS1.setD(D1);
       multiplyMPOMPS(MPO0, MPS0, epsEnv, maxNumSweepsEnv, errorAchievedEnv, numSweepsDoneEnv, MPS1);
       H0MPSs[col] = MPS1;
      }
     }
     Direction = "left";
    }
    else if (Direction == "left")
    {
// compute energy:
     PEPS0.getMPS(Direction, MPS1);
     MPS0 = NormMPSs[Ncols-2];
     normSquaredPEPS0 = MPS0.contractReverse(MPS1);
     PEPS0.getMPS(PEPOH0, Direction, MPS1);
     MPS0 = H0MPSs[Ncols-2];
     expectationValueH0 = MPS0.contractReverse(MPS1);
     newEnergy = expectationValueH0 / normSquaredPEPS0;
     if (abs(newEnergy - energy)/abs(newEnergy) <= eps)
     {
      Energies[numSweepsDone] = newEnergy;
      return;
     }
     energy = newEnergy;
     Energies[numSweepsDone] = newEnergy;
// sweep from right to left:
     for (col = Ncols-1; col > 0; col--)
     {
// sweep from right to left starting with column Ncols-1:
      if (col == Ncols-1)
      {
// compute norm tensors:
       MPS0 = NormMPSs[col-1]; PEPS0.getMPS(Direction, MPS1);
       MPS0.get(0, Tensor0); MPS1.get(Nrows-1, Tensor1);
       Indices01[0] = 2; Indices11[0] = 2;
       Tensor0.contract(Indices01, Tensor1, Indices11);
       Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
       Tensor0.permute(Order4);
       NormTensors[Nrows-1] = Tensor0;
       Indices01[0] = 2; Indices11[0] = 0;
       Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
       for (row = Nrows-2; row > 0; row--)
       {
        MPS0.get(Nrows-row-1, Tensor1);
        Tensor0.contract(Indices01, Tensor1, Indices11);
        MPS1.get(row, Tensor1);
        Tensor0.contract(Indices02, Tensor1, Indices12);
        NormTensors[row] = Tensor0;
       }
// compute H0 tensors:
       MPS0 = H0MPSs[col-1]; PEPS0.getMPS(PEPOH0, Direction, MPS1);
       MPS0.get(0, Tensor0); MPS1.get(Nrows-1, Tensor1);
       Indices01[0] = 2; Indices11[0] = 2;
       Tensor0.contract(Indices01, Tensor1, Indices11);
       Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
       Tensor0.permute(Order4);
       H0Tensors[Nrows-1] = Tensor0;
       Indices01[0] = 2; Indices11[0] = 0;
       Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 2;
       for (row = Nrows-2; row > 0; row--)
       {
        MPS0.get(Nrows-row-1, Tensor1);
        Tensor0.contract(Indices01, Tensor1, Indices11);
        MPS1.get(row, Tensor1);
        Tensor0.contract(Indices02, Tensor1, Indices12);
        H0Tensors[row] = Tensor0;
       }
       numSweepsDoneLocal = 0;
       DirectionLocal = "down";
       while (numSweepsDoneLocal < numSweepsLocal)
       {
        if (DirectionLocal == "down")
        {
// sweep from up to down starting with row 0:
         row = 0;
// update tensor:
         NormMPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = NormTensors[row+1];
         Indices01[0] = 0; Indices11[0] = 2;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         NormTensor = Tensor0;
         H0MPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = H0Tensors[row+1];
         Tensor0.contract(Indices01, Tensor1, Indices11);
         H0Tensor = Tensor0;
         PEPOH0.get(row, col, H0CenterTensor);
         PEPS0.get(row, col, Tensor0);
         this->updateTensor(NormTensor, H0Tensor, H0CenterTensor, Tensor0, cutoff, mode);
         PEPS0.set(row, col, Tensor0);
// compute norm tensor:
         Tensor0.complexConjugate(Tensor1);
         Indices01[0] = 4; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order8);
         Shape3[0] = D*D; Shape3[1] = 1; Shape3[2] = D*D;
         Tensor0.reshape(Shape3);
         Tensor1 = Tensor0; NormMPSs[col-1].get(Nrows-row-1, Tensor0);
         Indices01[0] = 2; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Order4[0] = 1; Order4[1] = 2; Order4[2] = 0; Order4[3] = 3;
         Tensor0.permute(Order4);
         NormTensors[row] = Tensor0;
// compute H0 tensor:
         PEPS0.get(row, col, Tensor0); PEPOH0.get(row, col, Tensor1);
         Indices01[0] = 4; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         PEPS0.get(row, col, Tensor1); Tensor1.complexConjugate();
         Indices01[0] = 8; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order12);
         Shape3[0] = D*DH0*D; Shape3[1] = 1; Shape3[2] = D*DH0*D;
         Tensor0.reshape(Shape3);
         Tensor1 = Tensor0; H0MPSs[col-1].get(Nrows-row-1, Tensor0);
         Indices01[0] = 2; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Order4[0] = 1; Order4[1] = 2; Order4[2] = 0; Order4[3] = 3;
         Tensor0.permute(Order4);
         H0Tensors[row] = Tensor0;
// sweep from up to down from row 1 to Nrows-2:
         for (row = 1; row < Nrows-1; row++)
         {
// update tensor:
          NormMPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = NormTensors[row-1];
          Indices01[0] = 1; Indices11[0] = 2;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor1 = NormTensors[row+1];
          Indices01[0] = 0; Indices11[0] = 2;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          NormTensor = Tensor0;
          H0MPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = H0Tensors[row-1];
          Indices01[0] = 1; Indices11[0] = 2;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor1 = H0Tensors[row+1];
          Indices01[0] = 0; Indices11[0] = 2;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          H0Tensor = Tensor0;
          PEPOH0.get(row, col, H0CenterTensor);
          PEPS0.get(row, col, Tensor0);
          this->updateTensor(NormTensor, H0Tensor, H0CenterTensor, Tensor0, cutoff, mode);
          PEPS0.set(row, col, Tensor0);
// compute norm tensor:
          Tensor0.complexConjugate(Tensor1);
          Indices01[0] = 4; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor0.permute(Order8);
          Shape3[0] = D*D; Shape3[1] = D*D; Shape3[2] = D*D;
          Tensor0.reshape(Shape3);
          Tensor2 = Tensor0;
          Tensor0 = NormTensors[row-1]; NormMPSs[col-1].get(Nrows-row-1, Tensor1);
          Indices01[0] = 2; Indices11[0] = 1;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 0;
          Tensor0.contract(Indices02, Tensor2, Indices12);
          NormTensors[row] = Tensor0;
// compute H0 tensor:
          PEPS0.get(row, col, Tensor0); PEPOH0.get(row, col, Tensor1);
          Indices01[0] = 4; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          PEPS0.get(row, col, Tensor1); Tensor1.complexConjugate();
          Indices01[0] = 8; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor0.permute(Order12);
          Shape3[0] = D*DH0*D; Shape3[1] = D*DH0*D; Shape3[2] = D*DH0*D;
          Tensor0.reshape(Shape3);
          Tensor2 = Tensor0;
          Tensor0 = H0Tensors[row-1]; H0MPSs[col-1].get(Nrows-row-1, Tensor1);
          Indices01[0] = 2; Indices11[0] = 1;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 1; Indices12[1] = 0;
          Tensor0.contract(Indices02, Tensor2, Indices12);
          H0Tensors[row] = Tensor0;
         }
         DirectionLocal = "up";
        }
        else if (DirectionLocal == "up")
        {
// sweep from down to up starting with row Nrows-1:
         row = Nrows-1;
// update tensor:
         NormMPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = NormTensors[row-1];
         Indices01[0] = 1; Indices11[0] = 2;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         NormTensor = Tensor0;
         H0MPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = H0Tensors[row-1];
         Tensor0.contract(Indices01, Tensor1, Indices11);
         H0Tensor = Tensor0;
         PEPOH0.get(row, col, H0CenterTensor);
         PEPS0.get(row, col, Tensor0);
         this->updateTensor(NormTensor, H0Tensor, H0CenterTensor, Tensor0, cutoff, mode);
         PEPS0.set(row, col, Tensor0);
// compute norm tensor:
         Tensor0.complexConjugate(Tensor1);
         Indices01[0] = 4; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order8);
         Shape3[0] = D*D; Shape3[1] = D*D; Shape3[2] = 1;
         Tensor0.reshape(Shape3);
         Tensor1 = Tensor0; NormMPSs[col-1].get(Nrows-row-1, Tensor0);
         Indices01[0] = 2; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
         Tensor0.permute(Order4);
         NormTensors[row] = Tensor0;
// compute H0 tensor:
         PEPS0.get(row, col, Tensor0); PEPOH0.get(row, col, Tensor1);
         Indices01[0] = 4; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         PEPS0.get(row, col, Tensor1); Tensor1.complexConjugate();
         Indices01[0] = 8; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order12);
         Shape3[0] = D*DH0*D; Shape3[1] = D*DH0*D; Shape3[2] = 1;
         Tensor0.reshape(Shape3);
         Tensor1 = Tensor0; H0MPSs[col-1].get(Nrows-row-1, Tensor0);
         Indices01[0] = 2; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Order4[0] = 0; Order4[1] = 3; Order4[2] = 1; Order4[3] = 2;
         Tensor0.permute(Order4);
         H0Tensors[row] = Tensor0;
// sweep from down to up from row Nrows-2 to 1:
         for (row = Nrows-2; row > 0; row--)
         {
// update tensor:
          NormMPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = NormTensors[row-1];
          Indices01[0] = 1; Indices11[0] = 2;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor1 = NormTensors[row+1];
          Indices01[0] = 0; Indices11[0] = 2;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          NormTensor = Tensor0;
          H0MPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = H0Tensors[row-1];
          Indices01[0] = 1; Indices11[0] = 2;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor1 = H0Tensors[row+1];
          Indices01[0] = 0; Indices11[0] = 2;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          H0Tensor = Tensor0;
          PEPOH0.get(row, col, H0CenterTensor);
          PEPS0.get(row, col, Tensor0);
          this->updateTensor(NormTensor, H0Tensor, H0CenterTensor, Tensor0, cutoff, mode);
          PEPS0.set(row, col, Tensor0);
// compute norm tensor:
          Tensor0.complexConjugate(Tensor1);
          Indices01[0] = 4; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor0.permute(Order8);
          Shape3[0] = D*D; Shape3[1] = D*D; Shape3[2] = D*D;
          Tensor0.reshape(Shape3);
          Tensor2 = Tensor0;
          Tensor0 = NormTensors[row+1]; NormMPSs[col-1].get(Nrows-row-1, Tensor1);
          Indices01[0] = 2; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 2; Indices12[1] = 0;
          Tensor0.contract(Indices02, Tensor2, Indices12);
          NormTensors[row] = Tensor0;
// compute H0 tensor:
          PEPS0.get(row, col, Tensor0); PEPOH0.get(row, col, Tensor1);
          Indices01[0] = 4; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          PEPS0.get(row, col, Tensor1); Tensor1.complexConjugate();
          Indices01[0] = 8; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor0.permute(Order12);
          Shape3[0] = D*DH0*D; Shape3[1] = D*DH0*D; Shape3[2] = D*DH0*D;
          Tensor0.reshape(Shape3);
          Tensor2 = Tensor0;
          Tensor0 = H0Tensors[row+1]; H0MPSs[col-1].get(Nrows-row-1, Tensor1);
          Indices01[0] = 2; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Indices02[0] = 2; Indices02[1] = 4; Indices12[0] = 2; Indices12[1] = 0;
          Tensor0.contract(Indices02, Tensor2, Indices12);
          H0Tensors[row] = Tensor0;
         }
         DirectionLocal = "down";
        }
        numSweepsDoneLocal++;
       }
       PEPS0.getMPS(Direction, MPS0);
       NormMPSs[col] = MPS0;
       PEPS0.getMPS(PEPOH0, Direction, MPS0);
       H0MPSs[col] = MPS0;
      }
// sweep from right to left from column Ncols-2 to 1:
      else
      {
// compute norm tensors:
       MPS0 = NormMPSs[col-1]; PEPS0.getMPO("right", col, MPO0); MPS1 = NormMPSs[col+1];
       MPS0.get(0, Tensor0); MPO0.get(0, Tensor1);
       Indices01[0] = 2; Indices11[0] = 2;
       Tensor0.contract(Indices01, Tensor1, Indices11);
       MPS1.get(Nrows-1, Tensor1);
       Indices01[0] = 4; Indices11[0] = 2;
       Tensor0.contract(Indices01, Tensor1, Indices11);
       Order6[0] = 0; Order6[1] = 2; Order6[2] = 5; Order6[3] = 1; Order6[4] = 3; Order6[5] = 4;
       Tensor0.permute(Order6);
       NormTensors[Nrows-1] = Tensor0;
       Indices01[0] = 3; Indices11[0] = 0;
       Indices02[0] = 3; Indices02[1] = 6;
       for (row = Nrows-2; row > 0; row--)
       {
        MPS0.get(Nrows-row-1, Tensor1);
        Tensor0.contract(Indices01, Tensor1, Indices11);
        MPO0.get(Nrows-row-1, Tensor1);
        Indices12[0] = 0; Indices12[1] = 2;
        Tensor0.contract(Indices02, Tensor1, Indices12);
        MPS1.get(row, Tensor1);
        Indices12[0] = 1; Indices12[1] = 2;
        Tensor0.contract(Indices02, Tensor1, Indices12);
        NormTensors[row] = Tensor0;
       }
// compute H0 tensors:
       MPS0 = H0MPSs[col-1]; PEPS0.getMPO(PEPOH0, "right", col, MPO0); MPS1 = H0MPSs[col+1];
       MPS0.get(0, Tensor0); MPO0.get(0, Tensor1);
       Indices01[0] = 2; Indices11[0] = 2;
       Tensor0.contract(Indices01, Tensor1, Indices11);
       MPS1.get(Nrows-1, Tensor1);
       Indices01[0] = 4; Indices11[0] = 2;
       Tensor0.contract(Indices01, Tensor1, Indices11);
       Tensor0.permute(Order6);
       H0Tensors[Nrows-1] = Tensor0;
       Indices01[0] = 3; Indices11[0] = 0;
       Indices02[0] = 3; Indices02[1] = 6;
       for (row = Nrows-2; row > 0; row--)
       {
        MPS0.get(Nrows-row-1, Tensor1);
        Tensor0.contract(Indices01, Tensor1, Indices11);
        MPO0.get(Nrows-row-1, Tensor1);
        Indices12[0] = 0; Indices12[1] = 2;
        Tensor0.contract(Indices02, Tensor1, Indices12);
        MPS1.get(row, Tensor1);
        Indices12[0] = 1; Indices12[1] = 2;
        Tensor0.contract(Indices02, Tensor1, Indices12);
        H0Tensors[row] = Tensor0;
       }
       numSweepsDoneLocal = 0;
       DirectionLocal = "down";
       while (numSweepsDoneLocal < numSweepsLocal)
       {
        if (DirectionLocal == "down")
        {
// sweep from up to down starting with row 0:
         row = 0;
// update tensor:
         NormMPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = NormTensors[row+1];
         Indices01[0] = 0; Indices11[0] = 3;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         NormMPSs[col+1].get(row, Tensor1);
         Indices01[0] = 6; Indices11[0] = 1;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order28);
         NormTensor = Tensor0;
         H0MPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = H0Tensors[row+1];
         Indices01[0] = 0; Indices11[0] = 3;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         H0MPSs[col+1].get(row, Tensor1);
         Indices01[0] = 6; Indices11[0] = 1;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order28);
         H0Tensor = Tensor0;
         PEPOH0.get(row, col, H0CenterTensor);
         PEPS0.get(row, col, Tensor0);
         this->updateTensor(NormTensor, H0Tensor, H0CenterTensor, Tensor0, cutoff, mode);
         PEPS0.set(row, col, Tensor0);
// compute norm tensor:
         Tensor0.complexConjugate(Tensor1);
         Indices01[0] = 4; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order8);
         Shape4[0] = D*D; Shape4[1] = 1; Shape4[2] = D*D; Shape4[3] = D*D;
         Tensor0.reshape(Shape4);
         Tensor1 = Tensor0; NormMPSs[col-1].get(Nrows-row-1, Tensor0);
         Indices01[0] = 2; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         NormMPSs[col+1].get(row, Tensor1);
         Indices01[0] = 3; Indices11[0] = 2;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Order6[0] = 1; Order6[1] = 2; Order6[2] = 4; Order6[3] = 0; Order6[4] = 3; Order6[5] = 5;
         Tensor0.permute(Order6);
         NormTensors[row] = Tensor0;
// compute H0 tensor:
         PEPS0.get(row, col, Tensor0); PEPOH0.get(row, col, Tensor1);
         Indices01[0] = 4; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         PEPS0.get(row, col, Tensor1); Tensor1.complexConjugate();
         Indices01[0] = 8; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order12);
         Shape4[0] = D*DH0*D; Shape4[1] = 1; Shape4[2] = D*DH0*D; Shape4[3] = D*DH0*D;
         Tensor0.reshape(Shape4);
         Tensor1 = Tensor0; H0MPSs[col-1].get(Nrows-row-1, Tensor0);
         Indices01[0] = 2; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         H0MPSs[col+1].get(row, Tensor1);
         Indices01[0] = 3; Indices11[0] = 2;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order6);
         H0Tensors[row] = Tensor0;
// sweep from up to down from row 1 to Nrows-2:
         for (row = 1; row < Nrows-1; row++)
         {
// update tensor:
          NormMPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = NormTensors[row-1];
          Indices01[0] = 1; Indices11[0] = 3;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          NormMPSs[col+1].get(row, Tensor1);
          Indices01[0] = 6; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor1 = NormTensors[row+1];
          Indices02[0] = 0; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 5;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          NormTensor = Tensor0;
          H0MPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = H0Tensors[row-1];
          Indices01[0] = 1; Indices11[0] = 3;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          H0MPSs[col+1].get(row, Tensor1);
          Indices01[0] = 6; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor1 = H0Tensors[row+1];
          Indices02[0] = 0; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 5;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          H0Tensor = Tensor0;
          PEPOH0.get(row, col, H0CenterTensor);
          PEPS0.get(row, col, Tensor0);
          this->updateTensor(NormTensor, H0Tensor, H0CenterTensor, Tensor0, cutoff, mode);
          PEPS0.set(row, col, Tensor0);
// compute norm tensor:
          Tensor0.complexConjugate(Tensor1);
          Indices01[0] = 4; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor0.permute(Order8);
          Shape4[0] = D*D; Shape4[1] = D*D; Shape4[2] = D*D; Shape4[3] = D*D;
          Tensor0.reshape(Shape4);
          Tensor2 = Tensor0;
          Tensor0 = NormTensors[row-1]; NormMPSs[col-1].get(Nrows-row-1, Tensor1);
          Indices01[0] = 3; Indices11[0] = 1;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 0;
          Tensor0.contract(Indices02, Tensor2, Indices12);
          NormMPSs[col+1].get(row, Tensor1);
          Indices02[0] = 3; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 2;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          NormTensors[row] = Tensor0;
// compute H0 tensor:
          PEPS0.get(row, col, Tensor0); PEPOH0.get(row, col, Tensor1);
          Indices01[0] = 4; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          PEPS0.get(row, col, Tensor1); Tensor1.complexConjugate();
          Indices01[0] = 8; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor0.permute(Order12);
          Shape4[0] = D*DH0*D; Shape4[1] = D*DH0*D; Shape4[2] = D*DH0*D; Shape4[3] = D*DH0*D;
          Tensor0.reshape(Shape4);
          Tensor2 = Tensor0;
          Tensor0 = H0Tensors[row-1]; H0MPSs[col-1].get(Nrows-row-1, Tensor1);
          Indices01[0] = 3; Indices11[0] = 1;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 0;
          Tensor0.contract(Indices02, Tensor2, Indices12);
          H0MPSs[col+1].get(row, Tensor1);
          Indices02[0] = 3; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 2;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          H0Tensors[row] = Tensor0;
         }
         DirectionLocal = "up";
        }
        else if (DirectionLocal == "up")
        {
// sweep from down to up starting with row Nrows-1:
         row = Nrows-1;
// update tensor:
         NormMPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = NormTensors[row-1];
         Indices01[0] = 1; Indices11[0] = 3;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         NormMPSs[col+1].get(row, Tensor1);
         Indices01[0] = 6; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         NormTensor = Tensor0;
         H0MPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = H0Tensors[row-1];
         Indices01[0] = 1; Indices11[0] = 3;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         H0MPSs[col+1].get(row, Tensor1);
         Indices01[0] = 6; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         H0Tensor = Tensor0;
         PEPOH0.get(row, col, H0CenterTensor);
         PEPS0.get(row, col, Tensor0);
         this->updateTensor(NormTensor, H0Tensor, H0CenterTensor, Tensor0, cutoff, mode);
         PEPS0.set(row, col, Tensor0);
// compute norm tensor:
         Tensor0.complexConjugate(Tensor1);
         Indices01[0] = 4; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order8);
         Shape4[0] = D*D; Shape4[1] = D*D; Shape4[2] = D*D; Shape4[3] = 1;
         Tensor0.reshape(Shape4);
         Tensor1 = Tensor0; NormMPSs[col-1].get(Nrows-row-1, Tensor0);
         Indices01[0] = 2; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         NormMPSs[col+1].get(row, Tensor1);
         Indices01[0] = 3; Indices11[0] = 2;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Order6[0] = 0; Order6[1] = 3; Order6[2] = 5; Order6[3] = 1; Order6[4] = 2; Order6[5] = 4;
         Tensor0.permute(Order6);
         NormTensors[row] = Tensor0;
// compute H0 tensor:
         PEPS0.get(row, col, Tensor0); PEPOH0.get(row, col, Tensor1);
         Indices01[0] = 4; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         PEPS0.get(row, col, Tensor1); Tensor1.complexConjugate();
         Indices01[0] = 8; Indices11[0] = 4;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order12);
         Shape4[0] = D*DH0*D; Shape4[1] = D*DH0*D; Shape4[2] = D*DH0*D; Shape4[3] = 1;
         Tensor0.reshape(Shape4);
         Tensor1 = Tensor0; H0MPSs[col-1].get(Nrows-row-1, Tensor0);
         Indices01[0] = 2; Indices11[0] = 0;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         H0MPSs[col+1].get(row, Tensor1);
         Indices01[0] = 3; Indices11[0] = 2;
         Tensor0.contract(Indices01, Tensor1, Indices11);
         Tensor0.permute(Order6);
         H0Tensors[row] = Tensor0;
// sweep from down to up from row Nrows-2 to 1:
         for (row = Nrows-2; row > 0; row--)
         {
// update tensor:
          NormMPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = NormTensors[row-1];
          Indices01[0] = 1; Indices11[0] = 3;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          NormMPSs[col+1].get(row, Tensor1);
          Indices01[0] = 6; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor1 = NormTensors[row+1];
          Indices02[0] = 0; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 5;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          NormTensor = Tensor0;
          H0MPSs[col-1].get(Nrows-row-1, Tensor0); Tensor1 = H0Tensors[row-1];
          Indices01[0] = 1; Indices11[0] = 3;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          H0MPSs[col+1].get(row, Tensor1);
          Indices01[0] = 6; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor1 = H0Tensors[row+1];
          Indices02[0] = 0; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 5;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          H0Tensor = Tensor0;
          PEPOH0.get(row, col, H0CenterTensor);
          PEPS0.get(row, col, Tensor0);
          this->updateTensor(NormTensor, H0Tensor, H0CenterTensor, Tensor0, cutoff, mode);
          PEPS0.set(row, col, Tensor0);
// compute norm tensor:
          Tensor0.complexConjugate(Tensor1);
          Indices01[0] = 4; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor0.permute(Order8);
          Shape4[0] = D*D; Shape4[1] = D*D; Shape4[2] = D*D; Shape4[3] = D*D;
          Tensor0.reshape(Shape4);
          Tensor2 = Tensor0;
          Tensor0 = NormTensors[row+1]; NormMPSs[col-1].get(Nrows-row-1, Tensor1);
          Indices01[0] = 3; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 0;
          Tensor0.contract(Indices02, Tensor2, Indices12);
          NormMPSs[col+1].get(row, Tensor1);
          Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          NormTensors[row] = Tensor0;
// compute H0 tensor:
          PEPS0.get(row, col, Tensor0); PEPOH0.get(row, col, Tensor1);
          Indices01[0] = 4; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          PEPS0.get(row, col, Tensor1); Tensor1.complexConjugate();
          Indices01[0] = 8; Indices11[0] = 4;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Tensor0.permute(Order12);
          Shape4[0] = D*DH0*D; Shape4[1] = D*DH0*D; Shape4[2] = D*DH0*D; Shape4[3] = D*DH0*D;
          Tensor0.reshape(Shape4);
          Tensor2 = Tensor0;
          Tensor0 = H0Tensors[row+1]; H0MPSs[col-1].get(Nrows-row-1, Tensor1);
          Indices01[0] = 3; Indices11[0] = 0;
          Tensor0.contract(Indices01, Tensor1, Indices11);
          Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 3; Indices12[1] = 0;
          Tensor0.contract(Indices02, Tensor2, Indices12);
          H0MPSs[col+1].get(row, Tensor1);
          Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
          Tensor0.contract(Indices02, Tensor1, Indices12);
          H0Tensors[row] = Tensor0;
         }
         DirectionLocal = "down";
        }
        numSweepsDoneLocal++;
       }
       MPS0 = NormMPSs[col+1]; PEPS0.getMPO(Direction, col, MPO0); MPS1 = MPS0;
       D1 = min(MPS0.getD()*MPO0.getD(), D2Env);
       MPS1.setD(D1);
       multiplyMPOMPS(MPO0, MPS0, epsEnv, maxNumSweepsEnv, errorAchievedEnv, numSweepsDoneEnv, MPS1);
       NormMPSs[col] = MPS1;
       MPS0 = H0MPSs[col+1]; PEPS0.getMPO(PEPOH0, Direction, col, MPO0); MPS1 = MPS0;
       D1 = min(MPS0.getD()*MPO0.getD(), D2Env);
       MPS1.setD(D1);
       multiplyMPOMPS(MPO0, MPS0, epsEnv, maxNumSweepsEnv, errorAchievedEnv, numSweepsDoneEnv, MPS1);
       H0MPSs[col] = MPS1;
      }
     }
     Direction = "right";
    }
    numSweepsDone++;
   }
// compute energy after maxNumSweeps global sweeps:
   if (Direction == "left")
   {
    PEPS0.getMPS(Direction, MPS1);
    MPS0 = NormMPSs[Ncols-2];
    normSquaredPEPS0 = MPS0.contractReverse(MPS1);
    PEPS0.getMPS(PEPOH0, Direction, MPS1);
    MPS0 = H0MPSs[Ncols-2];
    expectationValueH0 = MPS0.contractReverse(MPS1);
    Energies[numSweepsDone] = expectationValueH0 / normSquaredPEPS0;
   }
   else if (Direction == "right")
   {
    PEPS0.getMPS(Direction, MPS0);
    MPS1 = NormMPSs[1];
    normSquaredPEPS0 = MPS0.contractReverse(MPS1);
    PEPS0.getMPS(PEPOH0, Direction, MPS0);
    MPS1 = H0MPSs[1];
    expectationValueH0 = MPS0.contractReverse(MPS1);
    Energies[numSweepsDone] = expectationValueH0 / normSquaredPEPS0;
   }
  }
  else if (BC == "periodic")
  {
   cout << "The following function is not implemented for periodic boundary conditions: " <<
           "template<class T> void GroundState2D<T>::" <<
           "computeGroundState(const Hamiltonian2D<T>& H0, " <<
                              "unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv, " <<
                              "unsigned int numSweepsLocal, " <<
                              "double eps, unsigned int maxNumSweeps, " <<
                              "vector<double>& Energies, unsigned int& numSweepsDone, " <<
                              "PEPS<T>& PEPS0, " <<
                              "double cutoff = 1.0e-4, unsigned int mode = 0) const." << endl;
  }
 }
 else if (RepresentationH0 == "interactions")
 {
  cout << "The following function is not implemented for interactions representation: " <<
          "template<class T> void GroundState2D<T>::" <<
          "computeGroundState(const Hamiltonian2D<T>& H0, " <<
                             "unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv, " <<
                             "unsigned int numSweepsLocal, " <<
                             "double eps, unsigned int maxNumSweeps, " <<
                             "vector<double>& Energies, unsigned int& numSweepsDone, " <<
                             "PEPS<T>& PEPS0, " <<
                             "double cutoff = 1.0e-4, unsigned int mode = 0) const." << endl;
 }
 else if (RepresentationH0 == "matrix")
 {
  cout << "The following function is not implemented for matrix representation: " <<
          "template<class T> void GroundState2D<T>::" <<
          "computeGroundState(const Hamiltonian2D<T>& H0, " <<
                             "unsigned int D2Env, double epsEnv, unsigned int maxNumSweepsEnv, " <<
                             "unsigned int numSweepsLocal, " <<
                             "double eps, unsigned int maxNumSweeps, " <<
                             "vector<double>& Energies, unsigned int& numSweepsDone, " <<
                             "PEPS<T>& PEPS0, " <<
                             "double cutoff = 1.0e-4, unsigned int mode = 0) const." << endl;
 }
}

template<class T> void GroundState2D<T>::computeGroundStateExactly(const Hamiltonian2D<T>& H0,
                                                                   double& energy,
                                                                   vector<T>& State) const
{
 unsigned int Nrows = H0.getNrows(), Ncols = H0.getNcols();
 Matrix<unsigned int> d(Nrows, Ncols); H0.getd(d);
 unsigned int dim = 1;
 for (int j = 0; j < Ncols; j++)
 {
  for (int i = 0; i < Nrows; i++)
  {
   dim *= d(i, j);
  }
 }
#ifdef DEBUG
 if ((Nrows == 0) || (State.size() != dim))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void GroundState2D<T>::" <<
          "computeGroundStateExactly(const Hamiltonian2D<T>& H0, double& energy, " <<
                                    "vector<T>& State) const: " <<
          "((H0.getNrows() == 0) || (State.size() != d^{Nrows*Ncols}))." << endl;
  exit(1);
 }
#endif
 Matrix<T> MatrixH0(dim, dim); vector<double> W(dim); Matrix<T> Vr(dim, dim);
 H0.getMatrix(MatrixH0);
 MatrixH0.eigenDecompose(W, Vr);
 energy = W[0];
 for (int i = 0; i < dim; i++)
 {
  State[i] = Vr(i, 0);
 }
}

template<class T> GroundState2D<T>::GroundState2D() {}

template<class T> GroundState2D<T>::GroundState2D(const GroundState2D<T>& GroundState0) {}

template<class T> GroundState2D<T>& GroundState2D<T>::operator=(const GroundState2D<T>& GroundState0)
{}

template<class T> void GroundState2D<T>::updateTensor(Tensor<T>& NormTensor, Tensor<T>& H0Tensor,
                                                      Tensor<T>& H0CenterTensor, Tensor<T>& Tensor0,
                                                      double cutoff, unsigned int mode) const
{
 vector<unsigned int> H0CenterShape, Shape0;
 H0CenterTensor.getShape(H0CenterShape); Tensor0.getShape(Shape0);
#ifdef DEBUG
 if ((H0CenterShape.size() != 6) || (Shape0.size() != 5) || (H0CenterShape[4] != Shape0[4]) ||
     (H0CenterShape[5] != Shape0[4]) || (cutoff < 0.0) || ((mode != 0) && (mode != 1) && (mode != 2)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void GroundState2D<T>::" <<
          "updateTensor(Tensor<T>& NormTensor, Tensor<T>& H0Tensor, " <<
                       "Tensor<T>& H0CenterTensor, Tensor<T>& Tensor0, " <<
                       "double cutoff = 1.0e-4, unsigned int mode = 0) const: " <<
          "H0CenterTensor or Tensor0 has incorrect shape or " <<
          "((cutoff < 0.0) || ((mode != 0) && (mode != 1) && (mode != 2)))." << endl;
  exit(1);
 }
#endif
 vector<unsigned int> Shape8(8), Index8(8), Shape12(12), Indices04(4), Indices14(4), Index10(10);
 unsigned int dim = Shape0[0]*Shape0[1]*Shape0[2]*Shape0[3]*Shape0[4], i, j;
 Matrix<T> NormMatrix(dim, dim), OneMatrix(Shape0[4], Shape0[4]), H0Matrix(dim, dim);
// compute NormMatrix and H0Matrix:
// - reshape NormTensor:
 Shape8[0] = Shape0[0]; Shape8[1] = Shape0[0]; Shape8[2] = Shape0[1]; Shape8[3] = Shape0[1];
 Shape8[4] = Shape0[2]; Shape8[5] = Shape0[2]; Shape8[6] = Shape0[3]; Shape8[7] = Shape0[3];
 NormTensor.reshape(Shape8);
// - reshape H0Tensor:
 Shape12[0] = Shape0[0]; Shape12[1] = H0CenterShape[0]; Shape12[2] = Shape0[0];
 Shape12[3] = Shape0[1]; Shape12[4] = H0CenterShape[1]; Shape12[5] = Shape0[1];
 Shape12[6] = Shape0[2]; Shape12[7] = H0CenterShape[2]; Shape12[8] = Shape0[2];
 Shape12[9] = Shape0[3]; Shape12[10] = H0CenterShape[3]; Shape12[11] = Shape0[3];
 H0Tensor.reshape(Shape12);
// - define OneMatrix on physical indices:
 OneMatrix.fillZeroes();
 for (i = 0; i < Shape0[4]; i++)
  OneMatrix(i, i) = 1.0;
// - contract H0Tensor with H0CenterTensor:
 Indices04[0] = 1; Indices04[1] = 4; Indices04[2] = 7; Indices04[3] = 10;
 Indices14[0] = 0; Indices14[1] = 1; Indices14[2] = 2; Indices14[3] = 3;
 H0Tensor.contract(Indices04, H0CenterTensor, Indices14);
// - fill NormMatrix and H0Matrix:
 NormMatrix.fillZeroes(); H0Matrix.fillZeroes();
// row number i and column number j:
 for (int i4 = 0; i4 < Shape0[4]; i4++)
 {
  Index10[9] = i4;
 for (int j4 = 0; j4 < Shape0[4]; j4++)
 {
  Index10[8] = j4;
  for (int i3 = 0; i3 < Shape0[3]; i3++)
  {
   Index8[7] = i3;
   Index10[7] = i3;
  for (int j3 = 0; j3 < Shape0[3]; j3++)
  {
   Index8[6] = j3;
   Index10[6] = j3;
   for (int i2 = 0; i2 < Shape0[2]; i2++)
   {
    Index8[5] = i2;
    Index10[5] = i2;
   for (int j2 = 0; j2 < Shape0[2]; j2++)
   {
    Index8[4] = j2;
    Index10[4] = j2;
    for (int i1 = 0; i1 < Shape0[1]; i1++)
    {
     Index8[3] = i1;
     Index10[3] = i1;
    for (int j1 = 0; j1 < Shape0[1]; j1++)
    {
     Index8[2] = j1;
     Index10[2] = j1;
     for (int i0 = 0; i0 < Shape0[0]; i0++)
     {
      Index8[1] = i0;
      Index10[1] = i0;
     for (int j0 = 0; j0 < Shape0[0]; j0++)
     {
      Index8[0] = j0;
      Index10[0] = j0;
      i = i0 + i1*Shape0[0] + i2*Shape0[0]*Shape0[1] + i3*Shape0[0]*Shape0[1]*Shape0[2] +
          i4*Shape0[0]*Shape0[1]*Shape0[2]*Shape0[3];
      j = j0 + j1*Shape0[0] + j2*Shape0[0]*Shape0[1] + j3*Shape0[0]*Shape0[1]*Shape0[2] +
          j4*Shape0[0]*Shape0[1]*Shape0[2]*Shape0[3];
      NormMatrix(i, j) = NormTensor.get(Index8)*OneMatrix(i4, j4);
      H0Matrix(i, j) = H0Tensor.get(Index10);
     }
     }
    }
    }
   }
   }
  }
  }
 }
 }
 if (mode == 0)
 {
  unsigned int row, col;
  T element;
// diagonalize NormMatrixH*NormMatrix and determine NormMatrix=F*FH and pseudoinverses FInv and FHInv:
  Matrix<T> NormMatrixH(NormMatrix);
  NormMatrixH.transpose(); NormMatrixH.complexConjugate();
  NormMatrixH.multiply(NormMatrix);
  NormMatrixH.setType("hermitian");
  vector<double> W(dim); Matrix<T> Vr(dim, dim);
  NormMatrixH.eigenDecompose(W, Vr);
  unsigned int dimRed = 0;
  for (i = 0; i < dim; i++)
  {
   if (sqrt(W[i]) > cutoff)
    dimRed++;
  }
  Matrix<T> FInv(dimRed, dim), FHInv(dim, dimRed);
  col = 0;
  for (i = 0; i < dim; i++)
  {
   if (sqrt(W[i]) > cutoff)
   {
    for (row = 0; row < dim; row++)
    {
     element = Vr(row, i)/sqrt(sqrt(W[i]));
     FInv(col, row) = MathAuxiliary::complexConjugate(element);
     FHInv(row, col) = element;
    }
    col++;
   }
  }
  Matrix<T> H0MatrixRed(FInv);
  H0MatrixRed.multiply(H0Matrix); H0MatrixRed.multiply(FHInv);
  H0MatrixRed.setType("hermitian");
  vector<double> WRed(dimRed); Matrix<T> VrRed(dimRed, dimRed);
  H0MatrixRed.eigenDecompose(WRed, VrRed);
  vector<T> Vector0(dimRed), Vector1(dim);
  for (i = 0; i < dimRed; i++)
   Vector0[i] = VrRed(i, 0);
  FHInv.multiply(Vector0, Vector1);
  for (i = 0; i < dim; i++)
   Tensor0.set(i, Vector1[i]);
  Tensor0.normalize();
 }
 else if (mode == 1)
 {
  char balanc = 'B', jobvl = 'V', jobvr = 'V', sense = 'B';
  vector< complex<double> > Alpha(dim), Beta(dim);
  Matrix<T> Vl(dim, dim), Vr(dim, dim);
  int ilo, ihi; double abnrm, bbnrm;
  vector<double> Lscale(dim), Rscale(dim), Rconde(dim), Rcondv(dim);
  H0Matrix.eigenDecompose(NormMatrix, balanc, jobvl, jobvr, sense, Alpha, Beta, Vl, Vr, ilo, ihi, Lscale, Rscale,
                          abnrm, bbnrm, Rconde, Rcondv);
  unsigned int evalMinNum = 0; complex<double> evalMin = Alpha[0]/Beta[0], eval;
  for (i = 1; i < dim; i++)
  {
   eval = Alpha[i]/Beta[i];
   if (eval.real() < evalMin.real())
   {
    evalMinNum = i;
    evalMin = eval;
   }
  }
  for (i = 0; i < dim; i++)
   Tensor0.set(i, Vr(i, evalMinNum));
  Tensor0.normalize();
 }
 else if (mode == 2)
 {
  T norm, energy, element, result, resultOld, change, changeOld;
  unsigned int numSteps;
  vector<T> Vector0(dim), H0Vector(dim), NormVector(dim);
// set stepwidth delta:
  double delta, deltaStart = 0.0001, deltaStep = 0.1, deltaStop = 0.00001, eps = 1.0e-8, changeFac = 1.0;
  unsigned int numStepsTotal = 10000, numStepsLocal = 1000;
// determine Vector0=Tensor0:
  for (i = 0; i < dim; i++)
   Vector0[i] = Tensor0.get(i);
  vector<T> Vector0C(Vector0);
  for (delta = deltaStart; delta >= deltaStop; delta *= deltaStep)
  {
   norm = expectationValueMatrix(Vector0, NormMatrix);
   energy = expectationValueMatrix(Vector0, H0Matrix);
   resultOld = energy/norm;
   numSteps = 0;
// one step against the gradient:
   for (i = 0; i < numStepsLocal; i++)
   {
// H0Vector=H0Matrix*Vector0 and NormVector=NormMatrix*Vector0:
    H0Matrix.multiply(Vector0, H0Vector);
    NormMatrix.multiply(Vector0, NormVector);
// the remainder:
    element = 1.0/norm;
    multiply(element, H0Vector);
    element = -energy/(norm*norm);
    multiply(element, NormVector);
    add(H0Vector, NormVector);
    element = -delta;
    multiply(element, H0Vector);
    add(Vector0, H0Vector);
// compute norm=Vector0*NormMatrix*Vector0 and energy=Vector0*H0Matrix*Vector0:
    norm = expectationValueMatrix(Vector0, NormMatrix);
    energy = expectationValueMatrix(Vector0, H0Matrix);
   }
   result = energy/norm;
   change = abs(result-resultOld)/abs(result);
   if ((result < resultOld) && (change > eps))
   {
    do
    {
     numSteps += numStepsLocal;
     Vector0C = Vector0;
     resultOld = result;
     changeOld = change;
     for (i = 0; i < numStepsLocal; i++)
     {
// H0Vector=H0Matrix*Vector0 and NormVector=NormMatrix*Vector0:
      H0Matrix.multiply(Vector0, H0Vector);
      NormMatrix.multiply(Vector0, NormVector);
// the remainder:
      element = 1.0/norm;
      multiply(element, H0Vector);
      element = -energy/(norm*norm);
      multiply(element, NormVector);
      add(H0Vector, NormVector);
      element = -delta;
      multiply(element, H0Vector);
      add(Vector0, H0Vector);
// compute norm=Vector0*NormMatrix*Vector0 and energy=Vector0*H0Matrix*Vector0:
      norm = expectationValueMatrix(Vector0, NormMatrix);
      energy = expectationValueMatrix(Vector0, H0Matrix);
     }
     result = energy/norm;
     change = abs(result-resultOld)/abs(result);
    }
    while ((result < resultOld) && (abs(change) < changeFac*abs(changeOld)) && (change > eps) &&
           (numSteps < numStepsTotal));
   }
   Vector0 = Vector0C;
  }
  for (i = 0; i < dim; i++)
   Tensor0.set(i, Vector0[i]);
  Tensor0.normalize();
 }
}
