/// Abstract template class Operator implements interface for operators.
/** The abstract template class Operator implements the interface for operators in 1D.
    Every operator is represented by interactions, by a MPO or by a matrix. The interactions
    are implemented as a vector of positions and a vector of matrices where each vector entry
    corresponds to one term of the sum making up the operator. This class is
    parent class to the abstract template class Hamiltonian and to the abstract template class
    ProductOperator.
    \param Representation string, the representation, must be "interactions", "MPO" or "matrix"
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

template<class T> class Operator
{
 public:

/// Sets representation of Operator.
/** This function sets the representation of this Operator.
    \param Representation0 input: const string&, the representation of this Operator,
                                  must be "interactions", "MPO" or "matrix" */
  inline void setRepresentation(const string& Representation0);

/// Returns representation of Operator.
/** The returned representation is either "interactions", "MPO" or "matrix".
    \param Representation0 output: string&, the representation of this Operator, must
                                   be either "interactions", "MPO" or "matrix" */
  inline void getRepresentation(string& Representation0) const;

/// Returns Interactions representation of Operator.
/** This pure virtual function returns the interactions of this Operator. The interactions
    are implemented as a vector of positions and a vector of matrices where each vector entry
    corresponds to one term of the sum making up the operator.
    \param Positions output: vector< vector<unsigned int> >&, the positions of the matrices
    \param Matrices output: vector< vector< Matrix<T> > >&, the matrices */
  virtual void getInteractions(vector< vector<unsigned int> >& Positions,
                               vector< vector< Matrix<T> > >& Matrices) const = 0;

/// Returns MPO representation of Operator.
/** This pure virtual function returns this Operator as a MPO.
    \param MPO0 output: MPO<T>&, the MPO representing this Operator */
  virtual void getMPO(MPO<T>& MPO0) const = 0;

/// Returns matrix representation of Operator.
/** This pure virtual function returns this Operator as a matrix.
    \param Matrix0 output: Matrix<T>&, the matrix representing this Operator */
  virtual void getMatrix(Matrix<T>& Matrix0) const = 0;

 protected:

/// Representation.
  string Representation;
};

template<class T> inline void Operator<T>::setRepresentation(const string& Representation0)
{
#ifdef DEBUG
 if ((Representation0 != "interactions") && (Representation0 != "MPO"))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline void Operator<T>::" <<
          "setRepresentation(const string& Representation0)" << endl;
  exit(1);
 }
#endif
 this->Representation = Representation0;
}

template<class T> inline void Operator<T>::getRepresentation(string& Representation0) const
{
 Representation0 = this->Representation;
}
