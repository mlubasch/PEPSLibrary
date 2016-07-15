/// Functions related to C++'s vector class.
/** Vectors are used as linear arrays, e.g. to store physical states in exact calculations.
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

// Declaration of template classes:
template<class T> class Matrix;

/// Fills Vector0 with zeroes.
/** Vector0 is filled with zeroes. */
template<class T> void fillZeroes(vector<T>& Vector0);

/// Multiplies Vector0 with element.
/** This function multiplies Vector0 with a scalar element.
    \param element input: T, the scalar
    \param Vector0 input/output: vector<T>&, the vector */
template<class T> void multiply(T element, vector<T>& Vector0);

/// Adds Vector0 and Vector1.
/** This function adds Vector0 and Vector1 and stores the result in Vector0.
    \param Vector0 input/output: vector<T>&, on input one summand, on output the final sum
    \param Vector1 input: const vector<T>&, the other summand */
template<class T> void add(vector<T>& Vector0, const vector<T>& Vector1);

/// Takes scalar product of Vector0 with itself.
/** This function computes the scalar product of Vector0 with itself.
    \return double, the resulting value for the scalar product
    \sa template<class T> T normalize(vector<T>& Vector0)
    \sa template<class T> T scalarProduct(const vector<T>& Vector0, const vector<T>& Vector1) */
template<class T> double scalarProduct(const vector<T>& Vector0);

/// Normalizes vector.
/** Vector0 is normalized by computing the norm with
    template<class T> double scalarProduct(const vector<T>& Vector0)
    and then dividing each entry by norm. The resulting norm of Vector0 is returned.
    \return double, the norm
    \sa template<class T> double scalarProduct(const vector<T>& Vector0) */
template<class T> double normalize(vector<T>& Vector0);

/// Takes scalar product of Vector with another Vector.
/** This function computes the scalar product of Vector0 as bra with Vector1 as ket.
    \param Vector0 input: const vector<T>&, the bra vector
    \param Vector1 input: const vector<T>&, the ket vector
    \return T, the resulting value for the scalar product
    \sa template<class T> double scalarProduct(const vector<T>& Vector0) */
template<class T> T scalarProductVector(const vector<T>& Vector0, const vector<T>& Vector1);

/// Multiplies Vector0 with Vector1 to right as direct product and stores result in Vector2.
/** Vector0 is multiplied with Vector1 to the right as a direct product and the result is stored in
    Vector2.
    \param Vector0 input: const vector<T>&, the left vector of the direct product
    \param Vector1 input: const vector<T>&, the right vector of the direct product
    \param Vector2 output: vector<T>&, the resulting vector, on input it does not have to have a
                           specific shape, on output it will have a new shape and new elements */
template<class T> void multiplyDirectProduct(const vector<T>& Vector0, const vector<T>& Vector1,
                                             vector<T>& Vector2);

/// Computes expectation value of Matrix with Vector.
/** This function computes the expectation value of the input Matrix with the input Vector. It is friend
    function of template class Matrix.
    \param Vector0 input: const vector<T>&, the vector
    \param Matrix0 input: const Matrix<T>&, the Matrix
    \return T, the resulting expectation value
    \sa template<class T> T scalarProductMatrixVector(const vector<T>& Vector0,
                                                      const Matrix<T>& Matrix0,
                                                      const vector<T>& Vector1) */
template<class T> T expectationValueMatrix(const vector<T>& Vector0, const Matrix<T>& Matrix0);

/// Computes scalar product between Vector and Matrix applied to another Vector.
/** This function computes the scalar product of Vector0 as bra with Matrix0 applied to Vector1 as ket.
    It is friend function of template class Matrix.
    \param Vector0 input: const vector<T>&, the bra vector, must fulfill
                          Vector0.size()==Matrix0.getDim0()
    \param Matrix0 input: const Matrix<T>&, the Matrix
    \param Vector1 input: const vector<T>&, the ket vector, must fulfill
                          Vector1.size()==Matrix0.getDim1()
    \return T, the resulting scalar product
    \sa template<class T> T expectationValueMatrix(const vector<T>& Vector0,
                                                   const Matrix<T>& Matrix0) */
template<class T> T scalarProductMatrixVector(const vector<T>& Vector0, const Matrix<T>& Matrix0,
                                              const vector<T>& Vector1);

template<class T> void fillZeroes(vector<T>& Vector0)
{
#ifdef DEBUG
 if (Vector0.size() == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void " <<
          "fillZeroes(vector<T>& Vector0): " <<
          "(Vector0.size() == 0)." << endl;
  exit(1);
 }
#endif
 for (int i = 0; i < Vector0.size(); i++)
  Vector0[i] = 0.0;
}

template<class T> void multiply(T element, vector<T>& Vector0)
{
#ifdef DEBUG
 if (Vector0.size() == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void " <<
          "multiply(T element, vector<T>& Vector0): " <<
          "(Vector0.size() == 0)." << endl;
  exit(1);
 }
#endif
 for (int i = 0; i < Vector0.size(); i++)
  Vector0[i] *= element;
}

template<class T> void add(vector<T>& Vector0, const vector<T>& Vector1)
{
#ifdef DEBUG
 if ((Vector0.size() == 0) || (Vector0.size() != Vector1.size()))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void " <<
          "add(vector<T>& Vector0, vector<T>& Vector1): " <<
          "((Vector0.size() == 0) || (Vector0.size() != Vector1.size()))." << endl;
  exit(1);
 }
#endif
 for (int i = 0; i < Vector0.size(); i++)
  Vector0[i] += Vector1[i];
}

template<class T> double scalarProduct(const vector<T>& Vector0)
{
#ifdef DEBUG
 if (Vector0.size() == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> double " <<
          "scalarProduct(const vector<T>& Vector0): " <<
          "(Vector0.size() == 0)." << endl;
  exit(1);
 }
#endif
 T result = 0.0;
 for (int i = 0; i < Vector0.size(); i++)
  result += MathAuxiliary::complexConjugate(Vector0[i]) * Vector0[i];
 return MathAuxiliary::convertToDouble(result);
}

template<class T> double normalize(vector<T>& Vector0)
{
#ifdef DEBUG
 if (Vector0.size() == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> double " <<
          "normalize(vector<T>& Vector0): " <<
          "(Vector0.size() == 0)." << endl;
  exit(1);
 }
#endif
 double norm = sqrt(scalarProduct(Vector0));
 for (int i = 0; i < Vector0.size(); i++)
  Vector0[i] /= norm;
 return norm;
}

template<class T> T scalarProductVector(const vector<T>& Vector0, const vector<T>& Vector1)
{
#ifdef DEBUG
 if ((Vector0.size() == 0) || (Vector0.size() != Vector1.size()))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T " <<
          "scalarProductVector(const vector<T>& Vector0, const vector<T>& Vector1): " <<
          "((Vector0.size() == 0) || (Vector0.size() != Vector1.size()))." << endl;
  exit(1);
 }
#endif
 T result = 0.0;
 for (int i = 0; i < Vector0.size(); i++)
  result += MathAuxiliary::complexConjugate(Vector0[i]) * Vector1[i];
 return result;
}

template<class T> void multiplyDirectProduct(const vector<T>& Vector0, const vector<T>& Vector1,
                                             vector<T>& Vector2)
{
 unsigned int dim0 = Vector0.size(), dim1 = Vector1.size();
#ifdef DEBUG
 if ((dim0 == 0) || (dim1 == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void " <<
          "multiplyDirectProduct(const vector<T>& Vector0, const vector<T>& Vector1, " <<
                                "vector<T>& Vector2): " <<
          "((Vector0.size() == 0) || (Vector1.size() == 0))." << endl;
  exit(1);
 }
#endif
 unsigned int dim2 = dim0*dim1;
 Vector2 = vector<T>(dim2);
 int i2;
 for (int i0 = 0; i0 < dim0; i0++)
 {
  for (int i1 = 0; i1 < dim1; i1++)
  {
   i2 = i0*dim1+i1;
   Vector2[i2] = Vector0[i0]*Vector1[i1];
  }
 }
}

template<class T> T expectationValueMatrix(const vector<T>& Vector0, const Matrix<T>& Matrix0)
{
#ifdef DEBUG
 if ((Vector0.size() == 0) || (Vector0.size() != Matrix0.getDim0()) ||
     (Vector0.size() != Matrix0.getDim1()))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T " <<
          "expectationValueMatrix(const vector<T>& Vector0, const Matrix<T>& Matrix0): " <<
          "((Vector0.size() == 0) || (Vector0.size() != Matrix0.getDim0()) || " <<
           "(Vector0.size() != Matrix0.getDim1()))." << endl;
  exit(1);
 }
#endif
 unsigned int dim = Vector0.size();
 vector<T> Vector1(dim);
 Matrix0.multiply(Vector0, Vector1);
 T result = 0.0;
 for (int i = 0; i < dim; i++)
  result += MathAuxiliary::complexConjugate(Vector0[i]) * Vector1[i];
 return result;
}

template<class T> T scalarProductMatrixVector(const vector<T>& Vector0, const Matrix<T>& Matrix0,
                                              const vector<T>& Vector1)
{
 int m = Matrix0.Shape[0]; int n = Matrix0.Shape[1];
#ifdef DEBUG
 if ((m != Vector0.size()) || (n != Vector1.size()))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T " <<
          "scalarProductMatrixVector(const vector<T>& Vector0, const Matrix<T>& Matrix0, " <<
                                    "const vector<T>& Vector1)." << endl;
  exit(1);
 }
#endif
 vector<T> Vector2(m);
 Matrix0.multiply(Vector1, Vector2);
 T result = 0.0;
 for (int i = 0; i < m; i++)
  result += MathAuxiliary::complexConjugate(Vector0[i]) * Vector2[i];
 return result;
}
