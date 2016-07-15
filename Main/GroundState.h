/// Singleton template class GroundState implements ground state computation for MPS.
/** The singleton template class GroundState implements all functions for ground state computation of
    MPS.
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

template<class T> class GroundState
{
 public:

/// Constructs Singleton GroundState<T> and returns reference.
  static GroundState<T>& reference()
  {
   static GroundState<T> GS;
   return GS;
  }

/// Computes ground state.
/** This function computes the ground state of a given Hamiltonian H0 by minimizing
    <MPS0|H0|MPS0>/<MPS0|MPS0> directly. It starts with MPS0 and sweeps through the tensors, each time
    minimizing the energy. The energy is computed after each sweep and if the change in energy from one
    sweep to the next is below eps, the function returns. The function always returns if the number of
    sweeps done exceeds maxNumSweeps. This function is only implemented for a real template parameter T.
    \param H0 input: const Hamiltonian<T>&, the Hamiltonian
    \param eps input: double, the convergence precision, i.e. convergence in sweep n if
                      eps > abs(energy(n) - energy(n-1))
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param energy output: double&, the final energy
    \param numSweepsDone output: unsigned int&, the final number of sweeps done
    \param MPS0 input/output: MPS<T>&, on input the initial MPS, on output the final MPS, must have
                              the correct form */
  void computeGroundState(const Hamiltonian<T>& H0, double eps, unsigned int maxNumSweeps,
                          double& energy, unsigned int& numSweepsDone, MPS<T>& MPS0) const;

/// Computes low eigenstate.
/** This function computes a low lying eigenstate specified by which of a given Hamiltonian H0 by minimizing
    <MPS0|H0|MPS0>/<MPS0|MPS0> directly:
    which==0 implies the ground state, which==1 the first excited, and so on.
    It starts with MPS0 and sweeps through the tensors, each time minimizing the energy.
    The energy is computed after each sweep and if the change in energy from one sweep to the next is below
    eps, the function returns. The function always returns if the number of sweeps done exceeds
    maxNumSweeps.
    This function uses static void Matrix<T>::computeArnoldi to find the optimal tensors.
    T must be real.
    \param H0 input: const Hamiltonian<T>&, the Hamiltonian
    \param which input: unsigned int, specifies the desired eigenstate,
                        which==0 implies the ground state, which==1 the first excited, and so on
    \param eps input: double, the convergence precision,
                      i.e. convergence in sweep n if eps > abs(energy(n) - energy(n-1))
    \param maxNumSweeps input: unsigned int, the maximal number of sweeps allowed
    \param energy output: double&, the final energy
    \param numSweepsDone output: unsigned int&, the final number of sweeps done
    \param MPS0 input/output: MPS<T>&, on input the initial MPS, on output the final MPS,
                              must have the correct form
    \sa static void Matrix<T>::computeArnoldi(const Tensor<T>& TensorLeft, const Tensor<T>& TensorMiddle,
                                              const Tensor<T>& TensorRight, unsigned int which,
                                              Tensor<T>& Tensor0) */
  void computeLowEigenstate(const Hamiltonian<T>& H0, unsigned int which, double eps,
                            unsigned int maxNumSweeps, double& energy, unsigned int& numSweepsDone,
                            MPS<T>& MPS0) const;

/// Computes ground state and low-lying excitations with Arnoldi.
/** This function computes the ground state and low-lying excitations of a given Hamiltonian H0 with the
    Arnoldi method.
    It considers the Hamiltonian H0-sigma and starting from MPS0 builds up a Krylov basis.
    The Hamiltonian in the Krylov basis is diagonalized after 5 basis vectors are generated and then
    whenever a new basis vector is added.
    The five lowest lying energies are returned in Energies.
    \param H0 input: const Hamiltonian<T>&, the Hamiltonian
    \param sigma input: double, the shift of the energies
    \param MPS0 input: const MPS<T>&, the initial MPS starting the Krylov construction */
  void computeGroundStateArnoldi(const Hamiltonian<T>& H0, double sigma, MPS<T>& MPS0) const;

/// Computes ground state exactly.
/** This function computes the ground state of a given Hamiltonian H0 exactly by using its matrix
    representation and diagonalizing this matrix with LAPACK's routines for hermitian matrices. The
    underlying diagonalisation routine is a QR factorization.
    \param H0 input: const Hamiltonian<T>&, the Hamiltonian
    \param energy output: double&, the ground state energy
    \param State output: vector<T>&, the ground state */
  void computeGroundStateExactly(const Hamiltonian<T>& H0, double& energy, vector<T>& State) const;

 private:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  GroundState();

/// Standard copy constructor. Not implemented for this Singleton.
/** The standard copy constructor copies the input GroundState into this. This constructor must not
    be implemented for this Singleton class.
    \param GroundState0 input: const GroundState<T>&, to be copied into this
    \sa GroundState<T>& operator=(const GroundState<T>& GroundState0) */
  GroundState(const GroundState<T>& GroundState0);

/// Assigns GroundState to this. Not implemented for this Singleton.
/** The operator= allows to assign a GroundState to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side GroundState. This function must not be
    implemented for this Singleton class.
    \param GroundState0 input: const GroundState<T>&, to be copied into this
    \return GroundState<T>&, a reference to the new this
    \sa GroundState(const GroundState<T>& GroundState0) */
  GroundState<T>& operator=(const GroundState<T>& GroundState0);
};

template<class T> void GroundState<T>::computeGroundState(const Hamiltonian<T>& H0, double eps,
                                                          unsigned int maxNumSweeps,
                                                          double& energy, unsigned int& numSweepsDone,
                                                          MPS<T>& MPS0) const
{
 string BC; H0.getBC(BC);
 unsigned int N = H0.getN();
 vector<unsigned int> d(N); H0.getd(d);
#ifdef DEBUG
 string BCMPS0; MPS0.getBC(BCMPS0);
 vector<unsigned int> dMPS0; MPS0.getd(dMPS0);
 if ((BC != BCMPS0) || (N != MPS0.getN()) || (N == 0) || (d != dMPS0) || (eps <= 0.0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void GroundState<T>::" <<
          "computeGroundState(const Hamiltonian<T>& H0, double eps, unsigned int maxNumSweeps, " <<
                             "double& energy, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "H0 and MPS0 are not of the same form or (N == 0) or (eps <= 0.0)." << endl;
  exit(1);
 }
 else if ((typeid(T) != typeid(float)) && (typeid(T) != typeid(double)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void GroundState<T>::" <<
          "computeGroundState(const Hamiltonian<T>& H0, double eps, unsigned int maxNumSweeps, " <<
                             "double& energy, unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "This function is only implemented for a real template parameter T." << endl;
  exit(1);
 }
#endif
 Matrix<T> MatrixH0, X;
 T e; vector<T> V;
 string RepresentationH0; H0.getRepresentation(RepresentationH0);
 double expectationValueMPS0MPOH0MPS0, normSquaredMPS0, newEnergy;
 unsigned int dim, rest;
 energy = 1.0e6;
 numSweepsDone = 0;
 if (RepresentationH0 == "interactions")
 {
  cout << "The following function has not been implemented yet for interactions representation: " <<
          "computeGroundState(const Hamiltonian<T>& H0, double eps, unsigned int maxNumSweeps, " <<
                             "double& energy, unsigned int& numSweepsDone, MPS<T>& MPS0) " <<
                             "const." << endl;
 }
 else if (RepresentationH0 == "MPO")
 {
// We get the MPO representation of H0:
  MPO<T> MPOH0; H0.getMPO(MPOH0);
  if (BC == "open")
  {
// First we bring MPS0 into normal form from right to left:
    MPS0.bringIntoNormalForm();
// Then we construct the H for each site:
    Tensor<T> TensorsH[N];
    Tensor<T> Tensor0, Tensor1, Tensor2;
    vector<unsigned int> Indices0(1), Indices1(1), Order(6), Indices02(2), Indices12(2);
    MPS0.get(N-1, Tensor0); MPOH0.get(N-1, Tensor1);
    Indices0[0] = 2; Indices1[0] = 2;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    MPS0.get(N-1, Tensor1); Tensor1.complexConjugate();
    Indices0[0] = 4; Indices1[0] = 2;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    Order[0] = 1; Order[1] = 3; Order[2] = 5; Order[3] = 0; Order[4] = 2; Order[5] = 4;
    Tensor0.permute(Order);
    TensorsH[N-1] = Tensor0;
    for (int i = N-2; i > 0; i--)
    {
     MPS0.get(i, Tensor1);
     Indices0[0] = 3; Indices1[0] = 1;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     MPOH0.get(i, Tensor1);
     Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     MPS0.get(i, Tensor1); Tensor1.complexConjugate();
     Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     TensorsH[i] = Tensor0;
    }
   string direction = "right";
   vector<unsigned int> Indices06(6), Index02(2), Index3(3), Shape2(2), Shape3(3), Index06(6),
   Order0(3);
   while (numSweepsDone < maxNumSweeps)
   {
    if (direction == "right")
    {
     Tensor0 = TensorsH[1]; MPOH0.get(0, Tensor1);
     Indices0[0] = 4; Indices1[0] = 1;
     Tensor0.contract(Indices0, Tensor1, Indices1);
// we compute the energy:
     Tensor1 = Tensor0;
     MPS0.get(0, Tensor2);
     Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
     Tensor1.contract(Indices02, Tensor2, Indices12);
     MPS0.get(0, Tensor2); Tensor2.complexConjugate();
     Indices02[0] = 3; Indices02[1] = 5; Indices12[0] = 1; Indices12[1] = 2;
     Tensor1.contract(Indices02, Tensor2, Indices12);
     Index06[0] = 0; Index06[1] = 0; Index06[2] = 0; Index06[3] = 0; Index06[4] = 0; Index06[5] = 0;
     expectationValueMPS0MPOH0MPS0 = MathAuxiliary::convertToDouble(Tensor1.get(Index06));
     MPS0.get(0, Tensor1); Tensor1.complexConjugate(Tensor2);
     Indices02[0] = 1; Indices02[1] = 2; Indices12[0] = 1; Indices12[1] = 2;
     Tensor1.contract(Indices02, Tensor2, Indices12);
     Indices02[0] = 0; Indices02[1] = 0;
     normSquaredMPS0 = MathAuxiliary::convertToDouble(Tensor1.get(Indices02));
     newEnergy = expectationValueMPS0MPOH0MPS0 / normSquaredMPS0;
     if (abs(newEnergy - energy) <= eps)
     {
      energy = newEnergy;
      return;
     }
     energy = newEnergy;
// we prepare the sweep to the right:
// we compute the new tensor at position 0:
     Indices06[0] = 2; Indices06[1] = 4; Indices06[2] = 7; Indices06[3] = 0; Indices06[4] = 3;
     Indices06[5] = 6;
     Index02[0] = 0; Index02[1] = 0;
     Tensor0.getSubtensor(Indices06, Index02, Tensor1);
     MPS0.get(0, Tensor0); Tensor0.getShape(Shape3);
     dim = Shape3[0]*Shape3[1]*Shape3[2];
     Shape2[0] = dim; Shape2[1] = dim;
     Tensor1.reshape(Shape2);
     MatrixH0 = Matrix<T>(dim, dim);
     for (int i = 0; i < dim; i++)
     {
      Index02[0] = i;
      for (int j = 0; j < dim; j++)
      {
       Index02[1] = j;
       MatrixH0(i, j) = Tensor1.get(Index02);
      }
     }
     MatrixH0.setType("hermitian");
     V = vector<T>(dim);
     MatrixH0.computeGroundState(e, V);
     Tensor0 = Tensor<T>(Shape3);
     for (int i = 0; i < dim; i++)
     {
      Index3[2] = i / (Shape3[0]*Shape3[1]); rest = i % (Shape3[0]*Shape3[1]);
      Index3[1] = rest / Shape3[0]; rest = rest % Shape3[0];
      Index3[0] = rest;
      Tensor0.set(Index3, V[i]);
     }
     MPS0.set(0, Tensor0);
     MPS0.normalizeTensor(0, direction, X);
// we compute TensorsH[0]:
     MPS0.get(0, Tensor0); MPOH0.get(0, Tensor1);
     Indices0[0] = 2; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     MPS0.get(0, Tensor1); Tensor1.complexConjugate();
     Indices0[0] = 4; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Order[0] = 0; Order[1] = 2; Order[2] = 4; Order[3] = 1; Order[4] = 3; Order[5] = 5;
     Tensor0.permute(Order);
     TensorsH[0] = Tensor0;
// we do a sweep from left to right:
     for (int i = 1; i < N-1; i++)
     {
// we compute the new tensor at position i:
      MPOH0.get(i, Tensor1);
      Indices0[0] = 4; Indices1[0] = 0;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      Tensor1 = TensorsH[i+1];
      Indices0[0] = 5; Indices1[0] = 4;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      Indices06[0] = 4; Indices06[1] = 11; Indices06[2] = 6; Indices06[3] = 3; Indices06[4] = 10;
      Indices06[5] = 5;
      Index06[0] = 0; Index06[1] = 0; Index06[2] = 0; Index06[3] = 0; Index06[4] = 0; Index06[5] = 0;
      Tensor0.getSubtensor(Indices06, Index06, Tensor1);
      MPS0.get(i, Tensor0); Tensor0.getShape(Shape3);
      dim = Shape3[0]*Shape3[1]*Shape3[2];
      Shape2[0] = dim; Shape2[1] = dim;
      Tensor1.reshape(Shape2);
      MatrixH0 = Matrix<T>(dim, dim);
      for (int j = 0; j < dim; j++)
      {
       Index02[0] = j;
       for (int k = 0; k < dim; k++)
       {
        Index02[1] = k;
        MatrixH0(j, k) = Tensor1.get(Index02);
       }
      }
      MatrixH0.setType("hermitian");
      V = vector<T>(dim);
      MatrixH0.computeGroundState(e, V);
      Tensor0 = Tensor<T>(Shape3);
      for (int j = 0; j < dim; j++)
      {
       Index3[2] = j / (Shape3[0]*Shape3[1]); rest = j % (Shape3[0]*Shape3[1]);
       Index3[1] = rest / Shape3[0]; rest = rest % Shape3[0];
       Index3[0] = rest;
       Tensor0.set(Index3, V[j]);
      }
      MPS0.set(i, Tensor0);
      MPS0.normalizeTensor(i, direction, X);
// we compute TensorsH[i]:
      Tensor0 = TensorsH[i-1]; MPS0.get(i, Tensor1);
      Indices0[0] = 3; Indices1[0] = 0;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      MPOH0.get(i, Tensor1);
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      MPS0.get(i, Tensor1); Tensor1.complexConjugate();
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      TensorsH[i] = Tensor0;
     }
     MPS0.get(N-1, Tensor1);
     Indices0[0] = 1; Indices1[0] = 0;
     X.contract(Indices0, Tensor1, Indices1);
     MPS0.set(N-1, X);
     direction = "left";
    }
    else if (direction == "left")
    {
     Tensor0 = TensorsH[N-2]; MPOH0.get(N-1, Tensor1);
     Indices0[0] = 4; Indices1[0] = 0;
     Tensor0.contract(Indices0, Tensor1, Indices1);
// we compute the energy:
     Tensor1 = Tensor0;
     MPS0.get(N-1, Tensor2);
     Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
     Tensor1.contract(Indices02, Tensor2, Indices12);
     MPS0.get(N-1, Tensor2); Tensor2.complexConjugate();
     Indices02[0] = 3; Indices02[1] = 5; Indices12[0] = 0; Indices12[1] = 2;
     Tensor1.contract(Indices02, Tensor2, Indices12);
     Index06[0] = 0; Index06[1] = 0; Index06[2] = 0; Index06[3] = 0; Index06[4] = 0; Index06[5] = 0;
     expectationValueMPS0MPOH0MPS0 = MathAuxiliary::convertToDouble(Tensor1.get(Index06));
     MPS0.get(N-1, Tensor1); Tensor1.complexConjugate(Tensor2);
     Indices02[0] = 0; Indices02[1] = 2; Indices12[0] = 0; Indices12[1] = 2;
     Tensor1.contract(Indices02, Tensor2, Indices12);
     Indices02[0] = 0; Indices02[1] = 0;
     normSquaredMPS0 = MathAuxiliary::convertToDouble(Tensor1.get(Indices02));
     newEnergy = expectationValueMPS0MPOH0MPS0 / normSquaredMPS0;
     if (abs(newEnergy - energy) <= eps)
     {
      energy = newEnergy;
      return;
     }
     energy = newEnergy;
// we prepare the sweep to the left:
// we compute the new tensor at position N-1:
     Indices06[0] = 4; Indices06[1] = 2; Indices06[2] = 7; Indices06[3] = 3; Indices06[4] = 0;
     Indices06[5] = 6;
     Index02[0] = 0; Index02[1] = 0;
     Tensor0.getSubtensor(Indices06, Index02, Tensor1);
     MPS0.get(N-1, Tensor0); Tensor0.getShape(Shape3);
     dim = Shape3[0]*Shape3[1]*Shape3[2];
     Shape2[0] = dim; Shape2[1] = dim;
     Tensor1.reshape(Shape2);
     MatrixH0 = Matrix<T>(dim, dim);
     for (int i = 0; i < dim; i++)
     {
      Index02[0] = i;
      for (int j = 0; j < dim; j++)
      {
       Index02[1] = j;
       MatrixH0(i, j) = Tensor1.get(Index02);
      }
     }
     MatrixH0.setType("hermitian");
     V = vector<T>(dim);
     MatrixH0.computeGroundState(e, V);
     Tensor0 = Tensor<T>(Shape3);
     for (int i = 0; i < dim; i++)
     {
      Index3[2] = i / (Shape3[0]*Shape3[1]); rest = i % (Shape3[0]*Shape3[1]);
      Index3[1] = rest / Shape3[0]; rest = rest % Shape3[0];
      Index3[0] = rest;
      Tensor0.set(Index3, V[i]);
     }
     MPS0.set(N-1, Tensor0);
     MPS0.normalizeTensor(N-1, direction, X);
// we compute TensorsH[N-1]:
     MPS0.get(N-1, Tensor0); MPOH0.get(N-1, Tensor1);
     Indices0[0] = 2; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     MPS0.get(N-1, Tensor1); Tensor1.complexConjugate();
     Indices0[0] = 4; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Order[0] = 1; Order[1] = 3; Order[2] = 5; Order[3] = 0; Order[4] = 2; Order[5] = 4;
     Tensor0.permute(Order);
     TensorsH[N-1] = Tensor0;
// we do a sweep from right to left:
     for (int i = N-2; i >= 1; i--)
     {
// we compute the new tensor at position i:
      MPOH0.get(i, Tensor1);
      Indices0[0] = 4; Indices1[0] = 1;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      Tensor1 = TensorsH[i-1];
      Indices0[0] = 5; Indices1[0] = 4;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      Indices06[0] = 11; Indices06[1] = 4; Indices06[2] = 6; Indices06[3] = 10; Indices06[4] = 3;
      Indices06[5] = 5;
      Index06[0] = 0; Index06[1] = 0; Index06[2] = 0; Index06[3] = 0; Index06[4] = 0; Index06[5] = 0;
      Tensor0.getSubtensor(Indices06, Index06, Tensor1);
      MPS0.get(i, Tensor0); Tensor0.getShape(Shape3);
      dim = Shape3[0]*Shape3[1]*Shape3[2];
      Shape2[0] = dim; Shape2[1] = dim;
      Tensor1.reshape(Shape2);
      MatrixH0 = Matrix<T>(dim, dim);
      for (int j = 0; j < dim; j++)
      {
       Index02[0] = j;
       for (int k = 0; k < dim; k++)
       {
        Index02[1] = k;
        MatrixH0(j, k) = Tensor1.get(Index02);
       }
      }
      MatrixH0.setType("hermitian");
      V = vector<T>(dim);
      MatrixH0.computeGroundState(e, V);
      Tensor0 = Tensor<T>(Shape3);
      for (int j = 0; j < dim; j++)
      {
       Index3[2] = j / (Shape3[0]*Shape3[1]); rest = j % (Shape3[0]*Shape3[1]);
       Index3[1] = rest / Shape3[0]; rest = rest % Shape3[0];
       Index3[0] = rest;
       Tensor0.set(Index3, V[j]);
      }
      MPS0.set(i, Tensor0);
      MPS0.normalizeTensor(i, direction, X);
// we compute TensorsH[i]:
      Tensor0 = TensorsH[i+1]; MPS0.get(i, Tensor1);
      Indices0[0] = 3; Indices1[0] = 1;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      MPOH0.get(i, Tensor1);
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      MPS0.get(i, Tensor1); Tensor1.complexConjugate();
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      TensorsH[i] = Tensor0;
     }
     MPS0.get(0, Tensor0);
     Indices0[0] = 1; Indices1[0] = 0;
     Tensor0.contract(Indices0, X, Indices1);
     Order0[0] = 0; Order0[1] = 2; Order0[2] = 1;
     Tensor0.permute(Order0);
     MPS0.set(0, Tensor0);
     direction = "right";
    }
    numSweepsDone++;
   }
  }
  else if (BC == "periodic")
  {
   cout << "The following function has not been implemented yet for periodic boundary conditions: " <<
           "computeGroundState(const Hamiltonian<T>& H0, double eps, unsigned int maxNumSweeps, " <<
                              "double& energy, unsigned int& numSweepsDone, MPS<T>& MPS0) " <<
                              "const." << endl;
  }
 }
}

template<class T> void GroundState<T>::computeLowEigenstate(const Hamiltonian<T>& H0, unsigned int which,
                                                            double eps, unsigned int maxNumSweeps,
                                                            double& energy, unsigned int& numSweepsDone,
                                                            MPS<T>& MPS0) const
{
 string BC; H0.getBC(BC);
 unsigned int N = H0.getN();
 vector<unsigned int> d(N); H0.getd(d);
#ifdef DEBUG
 string BCMPS0; MPS0.getBC(BCMPS0);
 vector<unsigned int> dMPS0; MPS0.getd(dMPS0);
 if ((BC != BCMPS0) || (N != MPS0.getN()) || (N == 0) || (d != dMPS0) || (eps <= 0.0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void GroundState<T>::" <<
          "computeLowEigenstate(const Hamiltonian<T>& H0, unsigned int which, double eps, " <<
                               "unsigned int maxNumSweeps, double& energy, " <<
                               "unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "H0 and MPS0 are not of the same form or (N == 0) or (eps <= 0.0)." << endl;
  exit(1);
 }
 else if ((typeid(T) != typeid(float)) && (typeid(T) != typeid(double)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void GroundState<T>::" <<
          "computeLowEigenstate(const Hamiltonian<T>& H0, unsigned int which, double eps, " <<
                               "unsigned int maxNumSweeps, double& energy, " <<
                               "unsigned int& numSweepsDone, MPS<T>& MPS0) const: " <<
          "This function is only implemented for a real template parameter T." << endl;
  exit(1);
 }
#endif
 string RepresentationH0; H0.getRepresentation(RepresentationH0);
 double expectationValueMPS0MPOH0MPS0, normSquaredMPS0, newEnergy;
 energy = 1.0e6;
 numSweepsDone = 0;
 if (RepresentationH0 == "interactions")
 {
  cout << "The following function has not been implemented yet for interactions representation: " <<
          "template<class T> void GroundState<T>::" <<
          "computeLowEigenstate(const Hamiltonian<T>& H0, unsigned int which, double eps, " <<
                               "unsigned int maxNumSweeps, double& energy, " <<
                               "unsigned int& numSweepsDone, MPS<T>& MPS0) const." << endl;
 }
 else if (RepresentationH0 == "MPO")
 {
// We get the MPO representation of H0:
  MPO<T> MPOH0; H0.getMPO(MPOH0);
  if (BC == "open")
  {
   Matrix<T> X;
// First we construct the unit tensor:
   vector<unsigned int> Shape(6), Index(6);
   Shape[0] = 1; Shape[1] = 1; Shape[2] = 1; Shape[3] = 1; Shape[4] = 1; Shape[5] = 1;
   Tensor<T> UnitTensor(Shape);
   Index[0] = 0; Index[1] = 0; Index[2] = 0; Index[3] = 0; Index[4] = 0; Index[5] = 0;
   UnitTensor.set(Index, 1.0);
// We bring MPS0 into normal form from right to left:
   MPS0.bringIntoNormalForm();
// Then we construct the H for each site:
   vector< Tensor<T> > TensorsH(N);
   Tensor<T> Tensor0, Tensor1;
   vector<unsigned int> Indices0(1), Indices1(1), Indices02(2), Indices12(2), Order(6), Order3(3);
   MPS0.get(N-1, Tensor0); MPOH0.get(N-1, Tensor1);
   Indices0[0] = 2; Indices1[0] = 2;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   MPS0.get(N-1, Tensor1); Tensor1.complexConjugate();
   Indices0[0] = 4; Indices1[0] = 2;
   Tensor0.contract(Indices0, Tensor1, Indices1);
   Order[0] = 1; Order[1] = 3; Order[2] = 5; Order[3] = 0; Order[4] = 2; Order[5] = 4;
   Tensor0.permute(Order);
   TensorsH[N-1] = Tensor0;
   for (int i = N-2; i > 0; i--)
   {
    MPS0.get(i, Tensor1);
    Indices0[0] = 3; Indices1[0] = 1;
    Tensor0.contract(Indices0, Tensor1, Indices1);
    MPOH0.get(i, Tensor1);
    Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    MPS0.get(i, Tensor1); Tensor1.complexConjugate();
    Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
    Tensor0.contract(Indices02, Tensor1, Indices12);
    TensorsH[i] = Tensor0;
   }
   string direction = "right";
   while (numSweepsDone < maxNumSweeps)
   {
    if (direction == "right")
    {
// we compute the energy:
     Tensor0 = TensorsH[1]; MPS0.get(0, Tensor1);
     Indices0[0] = 3; Indices1[0] = 1;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     MPOH0.get(0, Tensor1);
     Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     MPS0.get(0, Tensor1); Tensor1.complexConjugate();
     Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     expectationValueMPS0MPOH0MPS0 = MathAuxiliary::convertToDouble(Tensor0.get(Index));
     MPS0.get(0, Tensor0); Tensor0.complexConjugate(Tensor1);
     Indices02[0] = 1; Indices02[1] = 2; Indices12[0] = 1; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     Indices02[0] = 0; Indices02[1] = 0;
     normSquaredMPS0 = MathAuxiliary::convertToDouble(Tensor0.get(Indices02));
     newEnergy = expectationValueMPS0MPOH0MPS0 / normSquaredMPS0;
     if (abs(newEnergy - energy) <= eps)
     {
      energy = newEnergy;
      return;
     }
     energy = newEnergy;
// we prepare the sweep to the right:
// we compute the new tensor at position 0:
     MPOH0.get(0, Tensor0);
     Matrix<T>::computeArnoldi(UnitTensor, Tensor0, TensorsH[1], which, Tensor1);
     MPS0.set(0, Tensor1);
     MPS0.normalizeTensor(0, direction, X);
// we compute TensorsH[0]:
     MPS0.get(0, Tensor0); MPOH0.get(0, Tensor1);
     Indices0[0] = 2; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     MPS0.get(0, Tensor1); Tensor1.complexConjugate();
     Indices0[0] = 4; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Order[0] = 0; Order[1] = 2; Order[2] = 4; Order[3] = 1; Order[4] = 3; Order[5] = 5;
     Tensor0.permute(Order);
     TensorsH[0] = Tensor0;
// we do a sweep from left to right:
     for (int i = 1; i < N-1; i++)
     {
// we compute the new tensor at position i:
      MPOH0.get(i, Tensor0);
      Matrix<T>::computeArnoldi(TensorsH[i-1], Tensor0, TensorsH[i+1], which, Tensor1);
      MPS0.set(i, Tensor1);
      MPS0.normalizeTensor(i, direction, X);
// we compute TensorsH[i]:
      Tensor0 = TensorsH[i-1]; MPS0.get(i, Tensor1);
      Indices0[0] = 3; Indices1[0] = 0;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      MPOH0.get(i, Tensor1);
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      MPS0.get(i, Tensor1); Tensor1.complexConjugate();
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      TensorsH[i] = Tensor0;
     }
     MPS0.get(N-1, Tensor1);
     Indices0[0] = 1; Indices1[0] = 0;
     X.contract(Indices0, Tensor1, Indices1);
     MPS0.set(N-1, X);
     direction = "left";
    }
    else if (direction == "left")
    {
// we compute the energy:
     Tensor0 = TensorsH[N-2]; MPS0.get(N-1, Tensor1);
     Indices0[0] = 3; Indices1[0] = 0;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     MPOH0.get(N-1, Tensor1);
     Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     MPS0.get(N-1, Tensor1); Tensor1.complexConjugate();
     Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 0; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     expectationValueMPS0MPOH0MPS0 = MathAuxiliary::convertToDouble(Tensor0.get(Index));
     MPS0.get(N-1, Tensor0); Tensor0.complexConjugate(Tensor1);
     Indices02[0] = 0; Indices02[1] = 2; Indices12[0] = 0; Indices12[1] = 2;
     Tensor0.contract(Indices02, Tensor1, Indices12);
     Indices02[0] = 0; Indices02[1] = 0;
     normSquaredMPS0 = MathAuxiliary::convertToDouble(Tensor0.get(Indices02));
     newEnergy = expectationValueMPS0MPOH0MPS0 / normSquaredMPS0;
     if (abs(newEnergy - energy) <= eps)
     {
      energy = newEnergy;
      return;
     }
     energy = newEnergy;
// we prepare the sweep to the left:
// we compute the new tensor at position N-1:
     MPOH0.get(N-1, Tensor0);
     Matrix<T>::computeArnoldi(TensorsH[N-2], Tensor0, UnitTensor, which, Tensor1);
     MPS0.set(N-1, Tensor1);
     MPS0.normalizeTensor(N-1, direction, X);
// we compute TensorsH[N-1]:
     MPS0.get(N-1, Tensor0); MPOH0.get(N-1, Tensor1);
     Indices0[0] = 2; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     MPS0.get(N-1, Tensor1); Tensor1.complexConjugate();
     Indices0[0] = 4; Indices1[0] = 2;
     Tensor0.contract(Indices0, Tensor1, Indices1);
     Order[0] = 1; Order[1] = 3; Order[2] = 5; Order[3] = 0; Order[4] = 2; Order[5] = 4;
     Tensor0.permute(Order);
     TensorsH[N-1] = Tensor0;
// we do a sweep from right to left:
     for (int i = N-2; i > 0; i--)
     {
// we compute the new tensor at position i:
      MPOH0.get(i, Tensor0);
      Matrix<T>::computeArnoldi(TensorsH[i-1], Tensor0, TensorsH[i+1], which, Tensor1);
      MPS0.set(i, Tensor1);
      MPS0.normalizeTensor(i, direction, X);
// we compute TensorsH[i]:
      Tensor0 = TensorsH[i+1]; MPS0.get(i, Tensor1);
      Indices0[0] = 3; Indices1[0] = 1;
      Tensor0.contract(Indices0, Tensor1, Indices1);
      MPOH0.get(i, Tensor1);
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      MPS0.get(i, Tensor1); Tensor1.complexConjugate();
      Indices02[0] = 3; Indices02[1] = 6; Indices12[0] = 1; Indices12[1] = 2;
      Tensor0.contract(Indices02, Tensor1, Indices12);
      TensorsH[i] = Tensor0;
     }
     MPS0.get(0, Tensor0);
     Indices0[0] = 1; Indices1[0] = 0;
     Tensor0.contract(Indices0, X, Indices1);
     Order3[0] = 0; Order3[1] = 2; Order3[2] = 1;
     Tensor0.permute(Order3);
     MPS0.set(0, Tensor0);
     direction = "right";
    }
    numSweepsDone++;
   }
  }
  else if (BC == "periodic")
  {
   cout << "The following function has not been implemented yet for periodic boundary conditions: " <<
           "template<class T> void GroundState<T>::" <<
           "computeLowEigenstate(const Hamiltonian<T>& H0, unsigned int which, double eps, " <<
                                "unsigned int maxNumSweeps, double& energy, " <<
                                "unsigned int& numSweepsDone, MPS<T>& MPS0) const." << endl;
  }
 }
}

template<class T> void GroundState<T>::computeGroundStateArnoldi(const Hamiltonian<T>& H0,
                                                                 double sigma,
                                                                 MPS<T>& MPS0) const
{
 unsigned int numEvalues = 5;
 unsigned int maxDim = 20;
 double eps = 1.0e-9; unsigned int maxNumSweeps = 100;
 double errorAchieved; unsigned int numSweepsDone;
 string BC; MPS0.getBC(BC);
 unsigned int N = MPS0.getN();
 vector<unsigned int> d(N); MPS0.getd(d);
 unsigned int D0 = MPS0.getD();
// MPS0:
 MPS0.normalize();
// Hamiltonian MPO:
 MPO<T> MPOH0; H0.getMPO(MPOH0);
 unsigned int DH = MPOH0.getD();
// Matrix K, H:
 Matrix<T> K(maxDim, maxDim), H(maxDim, maxDim);
 K.fillZeroes(); H.fillZeroes();
 K(0, 0) = MPS0.scalarProduct();
 H(0, 0) = expectationValueMPO(MPS0, MPOH0);
// Matrix Hamiltonian:
 Matrix<T> Hamiltonian0(1, 1), Hamiltonian1;
 Hamiltonian0(0, 0) = H(0, 0);
// Vector Coefficients:
 vector<T> Coefficients;
// Norm:
 T norm;
// Basis matrix:
 Matrix<T> Basis(maxDim, maxDim);
 Basis.fillZeroes();
 Basis(0, 0) = 1.0;
// Initializing:
 vector< MPS<T> > Krylov;
 Krylov.push_back(MPS0);
 MPS<T> MPS1, HMPS0, SigmaMPS0;
 vector< MPS<T> > MPSs(2);
 unsigned int D1 = 10;
 vector<unsigned int> Seed(N);
 for (int i = 0; i < N; i++)
  Seed[i] = 13*(i+1);
 T result;
 vector<double> W; Matrix<T> Vr;
 for (int i = 1; i < maxDim; i++)
 {
// 1. We compute the new Krylov vector |\phi_{i}>=(H-\sigma)|\phi_{i-1}> and normalize it:
  MPS0 = Krylov[i-1];
  D0 = MPS0.getD();
  HMPS0 = MPS<T>(BC, N, d, DH*D0);
  HMPS0.fillRandomly(Seed);
  HMPS0.normalize();
  multiplyMPOMPS(MPOH0, MPS0, eps, maxNumSweeps, errorAchieved, numSweepsDone, HMPS0);
  SigmaMPS0 = MPS0;
  SigmaMPS0.multiply(-sigma);
  MPSs[0] = HMPS0; MPSs[1] = SigmaMPS0;
  MPS1 = MPS<T>(BC, N, d, D1);
  D1 += 2;
  MPS1.fillRandomly(Seed);
  MPS1.normalize();
  MPS1.approximate(MPSs, eps, maxNumSweeps, errorAchieved, numSweepsDone);
  while (errorAchieved > 1.0e-4)
  {
   D1 += 5;
   MPS1 = MPS<T>(BC, N, d, D1);
   MPS1.fillRandomly(Seed);
   MPS1.normalize();
   MPS1.approximate(MPSs, eps, maxNumSweeps, errorAchieved, numSweepsDone);
  }
  MPS1.normalize();
  Krylov.push_back(MPS1);
// 2. We construct the new column and row i of K and of H:
  for (int j = 0; j < i; j++)
  {
   MPS0 = Krylov[j];
   K(j, i) = MPS0.scalarProduct(MPS1);
   K(i, j) = MathAuxiliary::complexConjugate(K(j, i));
   H(j, i) = scalarProductMPOMPS(MPS0, MPOH0, MPS1);
   H(i, j) = MathAuxiliary::complexConjugate(H(j, i));
  }
  MPS0 = Krylov[i];
  K(i, i) = MPS0.scalarProduct();
  H(i, i) = expectationValueMPO(MPS0, MPOH0);
// 3. We build the vector of Coefficients:
  Coefficients = vector<T>(i);
  for (int j = 0; j < i; j++)
  {
   Coefficients[j] = 0.0;
   for (int k = 0; k <= j; k++)
   {
    Coefficients[j] += MathAuxiliary::complexConjugate(Basis(k, j)) * K(k, i);
   }
  }
// 4. We compute the norm:
  norm = 1.0;
  for (int j = 0; j < i; j++)
  {
   norm -= Coefficients[j]*MathAuxiliary::complexConjugate(Coefficients[j]);
  }
// 5. We construct the new basis column i:
  for (int j = 0; j < i; j++)
  {
   result = 0.0;
   for (int k = j; k < i; k++)
   {
    result += Coefficients[k]*Basis(j, k);
   }
   Basis(j, i) = -result/sqrt(norm);
  }
  Basis(i, i) = 1.0/sqrt(norm);
// 6. We build the Hamiltonian and diagonalize it:
  Hamiltonian1 = Matrix<T>(i+1, i+1);
  for (int j = 0; j < i; j++)
  {
   for (int k = 0; k < i; k++)
   {
    Hamiltonian1(j, k) = Hamiltonian0(j, k);
   }
  }
  for (int j = 0; j < i; j++)
  {
   Hamiltonian1(j, i) = 0.0;
   for (int k = 0; k <= j; k++)
   {
    for (int l = 0; l <= i; l++)
    {
     Hamiltonian1(j, i) += MathAuxiliary::complexConjugate(Basis(k, j))*Basis(l, i)*H(k, l);
    }
   }
   Hamiltonian1(i, j) = MathAuxiliary::complexConjugate(Hamiltonian1(j, i));
  }
  Hamiltonian1(i, i) = 0.0;
  for (int k = 0; k <= i; k++)
  {
   for (int l = 0; l <= i; l++)
   {
    Hamiltonian1(i, i) += MathAuxiliary::complexConjugate(Basis(k, i))*Basis(l, i)*H(k, l);
   }
  }
  Hamiltonian0 = Hamiltonian1;
  Hamiltonian1.setType("hermitian");
  W = vector<double>(i+1);
  Vr = Matrix<T>(i+1, i+1);
  Hamiltonian1.eigenDecompose(W, Vr);
  cout << "Eigenvalue in iteration " << i << ": " << W[0] << endl;
 }
}

template<class T> void GroundState<T>::computeGroundStateExactly(const Hamiltonian<T>& H0,
                                                                 double& energy,
                                                                 vector<T>& State) const
{
 unsigned int N = H0.getN();
 vector<unsigned int> d(N); H0.getd(d);
 unsigned int dim = 1;
 for (int i = 0; i < N; i++)
 {
  dim *= d[i];
 }
#ifdef DEBUG
 if ((N == 0) || (State.size() != dim))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void GroundState<T>::" <<
          "computeGroundStateExactly(const Hamiltonian<T>& H0, double& energy, " <<
                                    "vector<T>& State) const: " <<
          "((N == 0) || (State.size() != d^{N}))." << endl;
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

template<class T> GroundState<T>::GroundState() {}

template<class T> GroundState<T>::GroundState(const GroundState<T>& GroundState0) {}

template<class T> GroundState<T>& GroundState<T>::operator=(const GroundState<T>& GroundState0)
{}
