/// Abstract template class Operator2D implements interface for operators in 2D.
/** The abstract template class Operator2D implements the interface for operators in 2D.
    Every operator is represented by Interactions, by a PEPO or by a Matrix. The Interactions
    are implemented as a vector of row positions, a vector of column positions, and a vector of matrices, where each
    vector entry corresponds to one Interaction, i.e. one term of the sum making up the operator.
    This class is parent class to the abstract template class Hamiltonian2D and to the template class
    ProductOperator2D.
    \param Representation string, the representation, is "Interactions", "PEPO" or "Matrix"
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

template<class T> class Operator2D
{
 public:

/// Sets representation.
/** This function sets the representation of this Operator2D.
    \param Representation0 input: const string&, the representation,
                                  must be "Interactions", "PEPO" or "Matrix" */
  void setRepresentation(const string& Representation0);

/// Returns representation.
/** This function returns the representation of this Operator2D.
    \param Representation0 output: string&, the representation of this Operator2D */
  void getRepresentation(string& Representation0) const { Representation0 = this->Representation; }

/// Returns Interactions representation.
/** This pure virtual function returns the Interactions of this Operator2D.
    The Interactions are implemented as a vector of row positions, a vector of column positions, and a vector of
    matrices, where each vector entry corresponds to one Interaction, i.e. one term of the sum making up the
    operator.
    \param PositionsRow output: vector< vector<unsigned int> >&, the row positions of the Interactions
    \param PositionsCol output: vector< vector<unsigned int> >&, the column positions of the Interactions
    \param Interactions output: vector< vector< Matrix<T> > >&, the Interactions */
  virtual void getInteractions(vector< vector<unsigned int> >& PositionsRow,
                               vector< vector<unsigned int> >& PositionsCol,
                               vector< vector< Matrix<T> > >& Interactions) const = 0;

/// Returns PEPO representation.
/** This pure virtual function returns this Operator2D as a PEPO.
    \param PEPO0 output: PEPO<T>&, the PEPO representing this Operator2D */
  virtual void getPEPO(PEPO<T>& PEPO0) const = 0;

/// Returns Matrix representation.
/** This pure virtual function returns this Operator2D as a matrix.
    \param Matrix0 output: Matrix<T>&, the matrix representing this Operator2D */
  virtual void getMatrix(Matrix<T>& Matrix0) const = 0;

 protected:

/// Representation.
/** The representation is either "Interactions", "PEPO", or "Matrix". */
  string Representation;

};

template<class T> void Operator2D<T>::setRepresentation(const string& Representation0)
{
#ifdef DEBUG
 if ((Representation0 != "Interactions") && (Representation0 != "PEPO") &&
     (Representation0 != "Matrix"))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Operator2D<T>::" <<
          "setRepresentation(const string& Representation0): " <<
          "((Representation0 != Interactions) && (Representation0 != PEPO) && " <<
           "(Representation0 != Matrix))." << endl;
  exit(1);
 }
#endif
 this->Representation = Representation0;
}
