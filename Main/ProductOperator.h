/// Abstract template class ProductOperator implements interface for product operators.
/** The abstract template class ProductOperator implements the interface for product operators in 1D.
    Every product operator is represented by interactions which are implemented as a vector of
    positions and a vector of matrices. This class is child class to the abstract template class
    Operator.
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

template<class T> class ProductOperator
{
 public:

/// Returns Interactions representation of ProductOperator.
/** This function returns the interactions of this ProductOperator. The interactions
    are implemented as a vector of positions and a vector of matrices.
    \param Positions output: vector< vector<unsigned int> >&, the positions of the matrices, the outer
                             vector must have only one entry
    \param Matrices output: vector< vector< Matrix<T> > >&, the matrices, the outer vector must have
                            only one entry */
  void getInteractions(vector< vector<unsigned int> >& Positions,
                       vector< vector< Matrix<T> > >& Matrices);

/// Returns MPO representation of ProductOperator.
/** This pure virtual function returns this ProductOperator as a MPO.
    \param MPO output: MPO<T>&, the MPO representing this ProductOperator */
  void getMPO(MPO<T>& MPO);
};
