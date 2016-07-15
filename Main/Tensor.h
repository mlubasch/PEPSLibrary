/// Template class Tensor implements tensors.
/** The template class Tensor implements tensors as building blocks of MPS and PEPS.
    \param rank unsigned int, the number of indices of the Tensor
    \param Shape vector<unsigned int>, the dimension of the indices of the Tensor
    \param size unsigned int, the number of elements of the Tensor
    \param Elements T*, the elements of the Tensor
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

// Declaration of template classes:
template<class T> class Tensor;
template<class T> class Matrix;
template<class T> class MPS;

// Declaration of friend functions:
template<class T> void getMPS(Tensor<T>& Tensor0, MPS<T>& MPS0);
template<class T> void getTensor(const MPS<T>& MPS0, Tensor<T>& Tensor0);

template<class T> class Tensor
{
 public:

/// Standard constructor.
/** The standard constructor sets all member attributes to 0. */
  Tensor();

/// Constructor for Tensor with specific Shape.
/** This constructor initializes a Tensor with a specific Shape.
    \param Shape0 input: const vector<unsigned int>&, contains the desired Shape */
  Tensor(const vector<unsigned int>& Shape0);

/// Standard copy constructor.
/** The standard copy constructor copies the input Tensor into this.
    \param Tensor0 input: const Tensor<T>&, to be copied into this
    \sa Tensor<T>& operator=(const Tensor<T>& Tensor0) */
  Tensor(const Tensor<T>& Tensor0);

/// Standard destructor.
/** The standard destructor deletes the elements of the Tensor. */
  ~Tensor();

/// Assigns Tensor to this.
/** The operator= allows to assign a Tensor to this. Hereby this is destroyed and newly
    constructed to be a copy of the right-hand side Tensor.
    \param Tensor0 input: const Tensor<T>&, to be copied into this
    \return Tensor<T>&, a reference to the new this
    \sa Tensor(const Tensor<T>& Tensor0) */
  Tensor<T>& operator=(const Tensor<T>& Tensor0);

/// Returns rank of Tensor.
/** The returned rank is the number of indices of the Tensor.
    \return unsigned int, the rank of the Tensor */
  unsigned int getRank() const { return this->rank; }

/// Sets Shape.
/** This function changes the Shape of this Tensor. This rank remains the same, but this
    size may change.
    If an index dimension is decreased, Shape0[index] < this->Shape[index], then the
    dispensable elements are simply discarded.
    If an index dimension is increased, Shape0[index] > this->Shape[index], then the new elements
    are set as random numbers multiplied by element.
    Element has the default value 0.0, in which case the new elements are zeroes.
    \param Shape0 input: const vector<unsigned int>&, the new Shape, must fulfill
                         Shape0.size()==this->rank
    \param element optional input: T, multiplies all random numbers
    \param seed optional input: unsigned int, if given it is the seed value for srand, else
                                srand is initialized with time(0) */
  void setShape(const vector<unsigned int>& Shape0, T element = 0.0, unsigned int seed = 0);

/// Returns Shape of Tensor.
/** The returned Shape is a vector that holds the extent of each index.
    \param Shape0 output: vector<unsigned int>&, a copy of the Shape of this */
  void getShape(vector<unsigned int>& Shape0) const { Shape0 = this->Shape; }

/// Returns size of Tensor.
/** The returned size is the number of elements of the Tensor.
    \return unsigned int, the size of this */
  unsigned int getSize() const { return this->size; }

/// Checks equivalence of Tensors.
/** This function returns true, if this Tensor and Tensor0 have the same Shape and if they have the
    same elements at each Index. Otherwise it returns false.
    \param Tensor0 input: const Tensor<T>&, the Tensor with which this Tensor is compared
    \return bool, is true if this Tensor is equal to Tensor0, otherwise is false
    \sa bool operator!=(const Tensor<T>& Tensor0) const */
  bool operator==(const Tensor<T>& Tensor0) const;

/// Checks inequivalence of Tensors.
/** This function returns true, if this Tensor and Tensor0 have different Shape or if they have
    different elements at the same Index. Otherwise it returns false.
    \param Tensor0 input: const Tensor<T>&, the Tensor with which this Tensor is compared
    \return bool, is true if this Tensor is not equal to Tensor0, otherwise is false
    \sa bool operator==(const Tensor<T>& Tensor0) const */
  bool operator!=(const Tensor<T>& Tensor0) const;

/// Sets element at Index.
/** Given an Index the corresponding element in Elements is set.
    \param Index input: const vector<unsigned int>&, the desired Index
    \param element input: T, to be written at Index */
  inline void set(const vector<unsigned int>& Index, T element);

/// Sets element at position.
/** Given a position the corresponding element in Elements is set.
    \param position input: unsigned int, the position
    \param element input: T, to be written at position */
  inline void set(unsigned int position, T element);

/// Returns element at Index.
/** Given an Index the corresponding element in Elements is returned.
    \param Index input: const vector<unsigned int>&, the desired Index
    \return T, a copy of the element at Index */
  inline T get(const vector<unsigned int>& Index) const;

/// Returns element at position.
/** Given a position the corresponding element in Elements is returned.
    \param position input: unsigned int, the position
    \return T, a copy of the element at position */
  inline T get(unsigned int position) const;

/// Writes this Tensor to binary file.
/** Given a file name FileName, a new binary file is constructed into which this Tensor is
    written. A Tensor is represented in a binary file by:
    {rank, Shape[0], ..., Shape[rank-1], size, Elements[0], ..., Elements[size-1]}.
    \param FileName input: const string&, the name for the new binary file to which this Tensor is
                           written
    \sa void read(const string& FileName) */
  void write(const string& FileName) const;

/// Reads Tensor from binary file.
/** Given a binary file called FileName, this Tensor is replaced by the Tensor in FileName.
    A Tensor is represented in a binary file by:
    {rank, Shape[0], ..., Shape[rank-1], size, Elements[0], ..., Elements[size-1]}.
    \param FileName input: const string&, the binary file from which this Tensor is read
    \sa void write(const string& FileName) const */
  void read(const string& FileName);

/// Normalizes this tensor.
/** This function normalizes this tensor by finding its largest element in absolute value
    and dividing all elements by this absolute value.
    \return T, the absolute value of the largest element in absolute value */
  T normalize();

/// Normalizes this tensor in Frobenius norm.
/** This function normalizes this tensor in Frobenius norm by computing
       frobeniusNorm = \sqrt{\sum_{i} |this->Elements[i]|^{2}}
    and dividing all elements by frobeniusNorm.
    \return T, the Frobenius norm */
  T frobeniusNormalize();

/// Converts float Tensor to complex<float>.
/** This function converts this Tensor of type float to complex<float>. This Tensor must have T==float.
    \param Tensor0 output: Tensor< complex<float> >&, the complex equivalent to this Tensor */
  void convertToComplex(Tensor< complex<float> >& Tensor0) const;

/// Converts double Tensor to complex<double>.
/** This function converts this Tensor of type double to complex<double>. This Tensor must have
    T==double.
    \param Tensor0 output: Tensor< complex<double> >&, the complex equivalent to this Tensor */
  void convertToComplex(Tensor< complex<double> >& Tensor0) const;

/// Returns subtensor of this Tensor.
/** This function returns the subtensor specified by Indices where the values of the remaining indices
    are fixed to be Index.
    \param Indices input: const vector<unsigned int>&, the indices of the subtensor
    \param Index input: const vector<unsigned int>&, the values for the remaining indices
    \param Tensor0 output: Tensor<T>&, the resulting Tensor */
  void getSubtensor(const vector<unsigned int>& Indices, const vector<unsigned int>& Index,
                    Tensor<T>& Tensor0) const;

/// Permutes Tensor according to Order.
/** Given a specific Order this Tensor is permuted according to it.
    \param Order input: const vector<unsigned int>&, the desired order */
  void permute(const vector<unsigned int>& Order);

/// Reshapes Tensor according to Shape.
/** Given a specific Shape Tensor is reshaped according to it.
    \param Shape0 input: const vector<unsigned int>&, the new Shape of this */
  void reshape(const vector<unsigned int>& Shape0);

/// Multiplies this Tensor with element.
/** This Tensor is multiplied by element.
    \param element input: T, the scalar with which this Tensor is multiplied */
  void multiply(T element);

/// Adds this with second Tensor and stores result in this.
/** This is added with Tensor0 and the result is stored in this. This and Tensor0 must have the same
    Shape.
    \param Tensor0 input: const Tensor<T>&, the second Tensor for addition, must have the same Shape as
                          this Tensor */
  void add(const Tensor<T>& Tensor0);

/// Takes scalar product of this tensor with itself.
/** This function computes the scalar product of this tensor with itself.
    \return T, the resulting value for the scalar product */
  T scalarProduct() const;

/// Contracts Tensor via index0 and index1.
/** Given index0 and index1 this is contracted.
    \param index0 input: unsigned int, the first index over which this is contracted
    \param index1 input: unsigned int, the second index over which this is contracted
    \sa void contract(const vector<unsigned int>& Indices, const Tensor<T>& Tensor0,
                      const vector<unsigned int>& Indices0)
    \sa void contract(const vector<unsigned int>& Indices, const Tensor<T>& Tensor0,
                      const vector<unsigned int>& Indices0, Tensor<T>& Tensor1) const */
  void contract(unsigned int index0, unsigned int index1);

/// Contracts this with second Tensor over Indices and stores result in this.
/** This is contracted with Tensor0 over given Indices and the result is stored in this. LAPACK's XGEMM
    is used.
    \param Indices input: const vector<unsigned int>&, the indices for contraction of this
    \param Tensor0 input: Tensor<T>&, the second Tensor for contraction, is altered on output
    \param Indices0 input: const vector<unsigned int>&, the indices for contraction of
                    Tensor0
    \sa void contract(unsigned int index0, unsigned int index1)
    \sa void contract(const vector<unsigned int>& Indices, const Tensor<T>& Tensor0,
                      const vector<unsigned int>& Indices0, Tensor<T>& Tensor1) const */
  void contract(const vector<unsigned int>& Indices, Tensor<T>& Tensor0,
                const vector<unsigned int>& Indices0);

/// Contracts this with second Tensor over Indices and stores result in third tensor.
/** This is contracted with Tensor0 over given Indices and the result is stored in Tensor1. LAPACK's
    XGEMM is used.
    \param Indices input: const vector<unsigned int>&, the indices for contraction of this
    \param Tensor0 input: const Tensor<T>&, the second Tensor for contraction
    \param Indices0 input: const vector<unsigned int>&, the indices for contraction of
                    Tensor0
    \param Tensor1 output: Tensor<T>&, the resulting Tensor
    \sa void contract(unsigned int index0, unsigned int index1)
    \sa void contract(const vector<unsigned int>& Indices, Tensor<T>& Tensor0,
                      const vector<unsigned int>& Indices0) */
  void contract(const vector<unsigned int>& Indices, const Tensor<T>& Tensor0,
                const vector<unsigned int>& Indices0, Tensor<T>& Tensor1) const;

/// This is filled with random entries.
/** This Tensor is filled with uniformly distributed random entries from [-1,1] which are multiplied
    by element. LAPACK's XLARNV is used.
    \param seed optional input: unsigned int, if given it is the seed value for srand, else srand is
                                initialized with time(0)
    \param element optional input: T, multiplies all random numbers
    \sa void Tensor<T>::fillZeroes() */
  void fillRandomly(unsigned int seed = 0, T element = 1.0);

/// Fills this Tensor with zeroes.
/** This Tensor is filled with zeroes.
    \sa void Tensor<T>::fillRandomly() */
  void fillZeroes();

/// Complex conjugates this Tensor.
/** The complex conjugate of this Tensor is returned in this Tensor. */
  void complexConjugate();

/// Constructs complex conjugate of this Tensor.
/** The complex conjugate of this Tensor is stored in the input Tensor.
    \param Tensor0 output: Tensor<T>&, the complex conjugate of this Tensor, is overwritten */
  void complexConjugate(Tensor<T>& Tensor0) const;

/// Returns SVD of this Tensor.
/** This function performs a SVD on this Tensor=Tensor0*Sigma*Tensor1 via the matrix having
    Indices0 as row and Indices1 as column multiindex. The elements in Indices0 and Indices1
    must be sorted in ascending order. On output, the new bond has range
       min(Dcut, min(dim0, dim1))
    and is addressed as the last index of Tensor0 and as the first index of Tensor1. As
    Tensor0 and Tensor1 will be isometric, Sigma will be a diagonal matrix storing the singular values.
    \param Indices0 input: const vector<unsigned int>&, the row multiindex,
                           must be sorted in ascending order
    \param Indices1 input: const vector<unsigned int>&, the column multiindex,
                           must be sorted in ascending order
    \param Dcut input: unsigned int, min(Dcut, min(dim0, dim1)) defines the number of singular
                       values kept and the range of the new bond
    \param Tensor0 output: Tensor<T>&, resulting isometric left tensor with new bond as last index
    \param Sigma output: Matrix<T>&, resulting diagonal matrix with min(Dcut, min(dim0, dim1))
                         singular values
    \param Tensor1 output: Tensor<T>&, resulting isometric right tensor with new bond as first index */
  void singularValueDecompose(const vector<unsigned int>& Indices0,
                              const vector<unsigned int>& Indices1,
                              unsigned int Dcut, Tensor<T>& Tensor0, Matrix<T>& Sigma,
                              Tensor<T>& Tensor1);

/// Returns QR decomposition of this Tensor.
/** This function performs a QR decomposition on
       this Tensor = TensorQ*TensorR
    via the matrix having Indices0 as row and Indices1 as column multiindex.
    The elements in Indices0 and Indices1 must be sorted in ascending order. On output, the new bond
    has range min(dim0, dim1) and is addressed as the last index of TensorQ and as the first index
    of TensorR. This Tensor is replaced by TensorQ.
    \param Indices0 input: const vector<unsigned int>&, the row multiindex,
                           must be sorted in ascending order
    \param Indices1 input: const vector<unsigned int>&, the column multiindex,
                           must be sorted in ascending order
    \param TensorR output: Tensor<T>&, resulting right tensor with new bond as first index */
  void QRDecompose(const vector<unsigned int>& Indices0, const vector<unsigned int>& Indices1,
                   Tensor<T>& TensorR);

/// Returns LQ decomposition of this Tensor.
/** This function performs a LQ decomposition on
       this Tensor = TensorL*TensorQ
    via the matrix having Indices0 as row and Indices1 as column multiindex.
    The elements in Indices0 and Indices1 must be sorted in ascending order. On output, the new bond
    has range min(dim0, dim1) and is addressed as the last index of TensorL and as the first index
    of TensorQ. This Tensor is replaced by TensorQ.
    \param Indices0 input: const vector<unsigned int>&, the row multiindex,
                           must be sorted in ascending order
    \param Indices1 input: const vector<unsigned int>&, the column multiindex,
                           must be sorted in ascending order
    \param TensorL output: Tensor<T>&, resulting left tensor with new bond as last index */
  void LQDecompose(const vector<unsigned int>& Indices0, const vector<unsigned int>& Indices1,
                   Tensor<T>& TensorL);

/// Returns SVD of this two body operator.
/** This function performs a SVD on this Tensor=U*Sigma*Vt and returns two tensors
    TensorLeft=U*sqrt(Sigma) and TensorRight=sqrt(Sigma)*Vt representing the SVD. This Tensor must
    be a two body operator, i.e. this->rank==4 and this->Shape[0]==this->Shape[2] and
    this->Shape[1]==this->Shape[3]. On output, this Tensor is destroyed.
    \param TensorLeft output: Tensor<T>&, the left tensor of the SVD, must have the correct form
    \param TensorRight output: Tensor<T>&, the right tensor of the SVD, must have the correct form
    \sa void Matrix<T>::twoBodyExponentSVD(T delta, unsigned int d0, unsigned int d1,
                                           Tensor<T>& TensorLeft, Tensor<T>& TensorRight) const */
  void twoBodyOperatorSVD(Tensor<T>& TensorLeft, Tensor<T>& TensorRight);

/// Returns NLDA.
/** This function approximates this Tensor by the NonLocal Density Approximation (NLDA)
       G(n_{0}, n_{1}, ..., n_{L-1}) = \sum_{x=0}^{range} \sum_{l=0}^{L-x-1} g_{x,l}(n_{l}, n_{l+x})   ,
    where range denotes the maximal range of the two-site terms.
    Each non-zero entry of this Tensor constitutes one training density. Given a training set of
    m densities and given n variables in the NLDA, we solve a linear least squares problem.
    If mode==0 then
       void linearLeastSquares(double rcond, Matrix<T>& b, unsigned int& rank, vector<double>& S)
    solves the linear least squares problem directly for the mxn problem matrix, which is reduced
    to a nxn problem matrix if mode==1.
    If mode==2 then an alternating least squares method is used instead.
    \param range input: unsigned int, the maximal range of the two-site terms
    \param rank output: unsigned int&, if mode==0 or mode==1 the rank of the problem matrix
    \param Sigma output: vector<double>&, if mode==0 or mode==1 the singular values of the problem matrix
    \param NLDA output: vector< Tensor<T> >&, the desired NLDA sorted in the order of the above sum
    \param cutoff optional input: the cutoff used for the pseudoinverse,
                                  default value is 1.0e-14
    \param mode optional input: the mode,
                                if mode==0 the least squares problem is solved directly via the mxn
                                problem matrix, whereas if mode==1 via the nxn problem matrix,
                                if mode==2 then alternating least squares is used,
                                default value is 0
    \return double, the final relative distance between this Tensor and the NLDA */
  double getNLDA(unsigned int range, unsigned int& rank, vector<double>& Sigma,
                 vector< Tensor<T> >& NLDA, double cutoff = 1.0e-14, unsigned int mode = 0) const;

/// Returns NLDA as polynomial interpolation.
/** This function assumes this Tensor to store the FHK values to the densities in nTensors.
    It approximates this Tensor by the NonLocal Density Approximation (NLDA)
       G(n_{0}, n_{1}, ..., n_{L-1}) = \sum_{x=0}^{range} \sum_{l=0}^{L-x-1} g_{x,l}(n_{l}, n_{l+x})   ,
    where range denotes the maximal range of the two-site terms and
       g_{x==0,l}(n_{l}) = \sum_{s_{l}=0}^{deg0} g_{x=0,l}^{s_{l}} n_{l}^{s_{l}}
    and
       g_{x!=0,l}(n_{l}, n_{l+x}) = \sum_{s_{l}=0, s_{l+x}=0}^{deg1} g_{x,l}^{s_{l}, s_{l+x}} n_{l}^{s_{l}} n_{l+x}^{s_{l+x}}
    are polynomial interpolations of degree deg0 and deg1.
    Given a training set of m densities and given n variables in the NLDA, we solve a linear least squares problem with
       void linearLeastSquares(double rcond, Matrix<T>& b, unsigned int& rank, vector<double>& S)   .
    \param nTensors input: const vector< Tensor<T> >&, the densities
    \param range input: unsigned int, the maximal range of the two-site terms
    \param deg0 input: unsigned int, the degree of the polynomial interpolation for the local terms to x==0
    \param deg1 input: unsigned int, the degree of the polynomial interpolation for the nonlocal terms to x!=0
    \param Sigma output: vector<double>&, the singular values of the problem matrix
    \param NLDA output: vector< Tensor<T> >&, the desired NLDA sorted in the order of the above sums
    \param cutoff optional input: the cutoff used for the pseudoinverse, default value is 1.0e-14
    \return double, the final relative distance between this Tensor and the NLDA */
  double getPolyNLDA(const vector< Tensor<T> >& nTensors, unsigned int range, unsigned int deg0,
                     unsigned int deg1, vector<double>& Sigma, vector< Tensor<T> >& NLDA, double cutoff = 1.0e-14) const;

/// Returns NLDA as polynomial interpolation from Taylor expansion around LDA.
/** This function assumes this Tensor to store the FHK values to the densities in nTensors.
    It approximates this Tensor by the NonLocal Density Approximation (NLDA)
       G(n_{0}, n_{1}, ..., n_{L-1}) = \sum_{x=0}^{range} \sum_{l=0}^{L-x-1} g_{x,l}(n_{l}, n_{l+x})   ,
    where range denotes the maximal range of the two-site terms and
       g_{x==0,l}(n_{l}) = \sum_{s_{l}=0}^{deg0} g_{x=0,l}^{s_{l}} (n_{l}-1)^{2s_{l}}
    and
       g_{x!=0,l}(n_{l}, n_{l+x}) = \sum_{s_{l}=1, s_{l+x}=1}^{deg1} g_{x,l}^{s_{l}, s_{l+x}} (n_{l}-1)^{s_{l}} (n_{l+x}-1)^{s_{l+x}}
    are polynomial interpolations of degree deg0 and deg1.
    Given a training set of m densities and given n variables in the NLDA, we solve a linear least squares problem with
       void linearLeastSquares(double rcond, Matrix<T>& b, unsigned int& rank, vector<double>& S)   .
    \param nTensors input: const vector< Tensor<T> >&, the densities
    \param range input: unsigned int, the maximal range of the two-site terms
    \param deg0 input: unsigned int, the degree of the polynomial interpolation for the local terms to x==0
    \param deg1 input: unsigned int, the degree of the polynomial interpolation for the nonlocal terms to x!=0
    \param Sigma output: vector<double>&, the singular values of the problem matrix
    \param NLDA output: vector< Tensor<T> >&, the desired NLDA sorted in the order of the above sums
    \param cutoff optional input: the cutoff used for the pseudoinverse, default value is 1.0e-14
    \return double, the final relative distance between this Tensor and the NLDA */
  double getTaylorNLDA(const vector< Tensor<T> >& nTensors, unsigned int range, unsigned int deg0,
                       unsigned int deg1, vector<double>& Sigma, vector< Tensor<T> >& NLDA, double cutoff = 1.0e-14) const;

/// Returns NLDA as polynomial interpolation from Taylor expansion around LDA.
/** This function assumes this Tensor to store the FHK values to the densities in nTensors.
    It approximates this Tensor by the NonLocal Density Approximation (NLDA)
       G(n_{0}, n_{1}, ..., n_{L-1}) = \sum_{x=0}^{range} \sum_{l=0}^{L-x-1} g_{x,l}(n_{l}, n_{l+x})   ,
    where range denotes the maximal range of the two-site terms and
       g_{x==0,l}(n_{l}) = \sum_{s_{l}=0}^{deg0} g_{x=0,l}^{s_{l}} (n_{l}-1)^{2s_{l}}
    and
       g_{x!=0,l}(n_{l}, n_{l+x}) = \sum_{s_{l}=1, s_{l+x}=1}^{deg1} g_{x,l}^{s_{l}, s_{l+x}} (n_{l}-1)^{s_{l}} (n_{l+x}-1)^{s_{l+x}}
    are polynomial interpolations of degree deg0 and deg1, and we include only the coefficients
    g_{x,l}^{s_{l}, s_{l+x}} with (s_{l}+s_{l+x})%2 == 0   .
    Given a training set of m densities and given n variables in the NLDA, we solve a linear least squares problem with
       void linearLeastSquares(double rcond, Matrix<T>& b, unsigned int& rank, vector<double>& S)   .
    \param nTensors input: const vector< Tensor<T> >&, the densities
    \param range input: unsigned int, the maximal range of the two-site terms
    \param deg0 input: unsigned int, the degree of the polynomial interpolation for the local terms to x==0
    \param deg1 input: unsigned int, the degree of the polynomial interpolation for the nonlocal terms to x!=0
    \param Sigma output: vector<double>&, the singular values of the problem matrix
    \param NLDA output: vector< Tensor<T> >&, the desired NLDA sorted in the order of the above sums
    \param cutoff optional input: the cutoff used for the pseudoinverse, default value is 1.0e-14
    \return double, the final relative distance between this Tensor and the NLDA */
  double getTaylorNLDA2(const vector< Tensor<T> >& nTensors, unsigned int range, unsigned int deg0,
                        unsigned int deg1, vector<double>& Sigma, vector< Tensor<T> >& NLDA, double cutoff = 1.0e-14) const;

/// Returns MPS representation.
/** This friend function computes the MPS MPS0 with open BC representing the Tensor Tensor0.
    It performs successive singular value decompositions on Tensor0 with
       void singularValueDecompose(const vector<unsigned int>& Indices0,
                                   const vector<unsigned int>& Indices1,
                                   unsigned int Dcut, Tensor<T>& Tensor0, Matrix<T>& Sigma,
                                   Tensor<T>& Tensor1)   .
    The exact MPS representation is returned.
    \param Tensor0 input: Tensor<T>&, the Tensor
    \param MPS0 output: MPS<T>&, the MPS with open BC representing Tensor0
    \sa void singularValueDecompose(const vector<unsigned int>& Indices0,
                                    const vector<unsigned int>& Indices1,
                                    unsigned int Dcut, Tensor<T>& Tensor0, Matrix<T>& Sigma,
                                    Tensor<T>& Tensor1)
    \sa void getTensor(const MPS<T>& MPS0, Tensor<T>& Tensor0) */
  friend void getMPS<>(Tensor<T>& Tensor0, MPS<T>& MPS0);

/// Returns Tensor representation.
/** This friend function computes the Tensor Tensor0 representing the MPS MPS0.
    The range 1 boundary indices are removed in the case of open boundary conditions.
    \param MPS0 input: const MPS<T>&, the MPS
    \param Tensor0 output: Tensor<T>&, the Tensor representing MPS0
    \sa void getMPS(Tensor<T>& Tensor0, MPS<T>& MPS0) */
  friend void getTensor<>(const MPS<T>& MPS0, Tensor<T>& Tensor0);

/// Returns position to Index.
/** Given an Index the corresponding position in Elements is returned.
    \param Index input: const vector<unsigned int>&, the desired Index
    \return unsigned int, the desired position in Elements */
  inline unsigned int getPosition(const vector<unsigned int>& Index) const;

/// Returns position to Index in specific Shape.
/** Assuming a Tensor of a specific Shape and given an Index the corresponding position in
    Elements is returned.
    \param Shape0 input: const vector<unsigned int>&, the assumed Shape
    \param Index0 input: const vector<unsigned int>&, the desired Index
    \return unsigned int, the desired position */
  inline static unsigned int getPosition(const vector<unsigned int>& Shape0,
                                         const vector<unsigned int>& Index0);

/// Returns Index to position.
/** Given a position in Elements the corresponding Index is returned.
    \param position input: unsigned int, the desired position
    \param Index output: vector<unsigned int>&, the resulting Index */
  inline void getIndex(unsigned int position, vector<unsigned int>& Index) const;

/// Returns Index to position in specific Shape.
/** Assuming a Tensor of specific Shape and given a position in Elements the corresponding Index
    is returned.
    \param Shape0 input: const vector<unsigned int>&, the assumed Shape
    \param position0 input: unsigned int, the desired position
    \param Index0 output: vector<unsigned int>&, the resulting Index */
  inline static void getIndex(const vector<unsigned int>& Shape0, unsigned int position0,
                              vector<unsigned int>& Index0);

 protected:

/// Rank of Tensor.
  unsigned int rank;

/// Shape of Tensor.
  vector<unsigned int> Shape;

/// Number of elements of Tensor.
  unsigned int size;

/// Elements of Tensor.
  T* Elements;

};

template<class T> Tensor<T>::Tensor()
{
 this->rank = 0;
 this->size = 0;
 this->Elements = 0;
}

template<class T> Tensor<T>::Tensor(const vector<unsigned int>& Shape0)
{
 this->rank = Shape0.size();
 this->Shape = Shape0;
 this->size = 1;
 for (int i = 0; i < this->rank; i++)
 {
  this->size *= this->Shape[i];
 }
 this->Elements = new T[this->size];
}

template<class T> Tensor<T>::Tensor(const Tensor<T>& Tensor0)
{
 this->rank = Tensor0.rank;
 this->Shape = Tensor0.Shape;
 this->size = Tensor0.size;
 this->Elements = new T[this->size];
 for (int i = 0; i < this->size; i++)
 {
  this->Elements[i] = Tensor0.Elements[i];
 }
}

template<class T> Tensor<T>::~Tensor()
{
 delete[] this->Elements;
}

template<class T> Tensor<T>& Tensor<T>::operator=(const Tensor<T>& Tensor0)
{
 if (this != &Tensor0)
 {
  this->rank = Tensor0.rank;
  this->Shape = Tensor0.Shape;
  this->size = Tensor0.size;
  delete[] this->Elements;
  this->Elements = new T[this->size];
  for (int i = 0; i < this->size; i++)
  {
   this->Elements[i] = Tensor0.Elements[i];
  }
 }
 return *this;
}

template<class T> void Tensor<T>::setShape(const vector<unsigned int>& Shape0, T element, unsigned int seed)
{
#ifdef DEBUG
 if (this->size == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Tensor<T>::" <<
          "setShape(const vector<unsigned int>& Shape0, T element, unsigned int seed): " <<
          "(this->size == 0)." << endl;
  exit(1);
 }
 if (Shape0.size() != this->rank)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Tensor<T>::" <<
          "setShape(const vector<unsigned int>& Shape0, T element, unsigned int seed): " <<
          "(Shape0.size() != this->rank)." << endl;
  exit(1);
 }
#endif
 if (Shape0 == this->Shape)
 {
  return;
 }
 bool IndexIsInside;
 vector<unsigned int> Index(this->rank);
 Tensor<T> Tensor0(Shape0);
// 1. fill Tensor0:
 if (element == 0.0)
 {
  Tensor0.fillZeroes();
 }
 else
 {
  Tensor0.fillRandomly(seed, element);
 }
// 2. write this Tensor into Tensor0:
 for (int pos = 0; pos < this->size; pos++)
 {
  this->getIndex(pos, Index);
  IndexIsInside = true;
  for (int i = 0; i < this->rank; i++)
  {
   if (Index[i] >= Shape0[i])
   {
    IndexIsInside = false;
   }
  }
  if (IndexIsInside)
  {
   Tensor0.set(Index, this->get(Index));
  }
 }
// 3. write Tensor0 into this Tensor:
 this->Shape = Shape0;
 this->size = Tensor0.size;
 delete[] this->Elements;
 this->Elements = new T[this->size];
 for (int i = 0; i < this->size; i++)
 {
  this->Elements[i] = Tensor0.Elements[i];
 }
}

template<class T> bool Tensor<T>::operator==(const Tensor<T>& Tensor0) const
{
 if (this->Shape != Tensor0.Shape)
 {
  return false;
 }
 else if (this->Shape == Tensor0.Shape)
 {
  for (int i = 0; i < this->size; i++)
  {
   if (this->Elements[i] != Tensor0.Elements[i])
    return false;
  }
  return true;
 }
}

template<class T> bool Tensor<T>::operator!=(const Tensor<T>& Tensor0) const
{
 if (this->Shape != Tensor0.Shape)
 {
  return true;
 }
 else if (this->Shape == Tensor0.Shape)
 {
  for (int i = 0; i < this->size; i++)
  {
   if (this->Elements[i] != Tensor0.Elements[i])
    return true;
  }
  return false;
 }
}

template<class T> inline void Tensor<T>::set(const vector<unsigned int>& Index, T element)
{
#ifdef DEBUG
 if (Index.size() != this->rank)
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline void Tensor<T>::" <<
          "set(const vector<unsigned int>& Index, T element)" << endl;
  exit(1);
 }
 for (int i = 0; i < this->rank; i++)
 {
  if (Index[i] >= this->Shape[i])
  {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline void Tensor<T>::" <<
          "set(const vector<unsigned int>& Index, T element)" << endl;
  exit(1);
  }
 }
#endif
 unsigned int position = this->getPosition(Index);
 this->Elements[position] = element;
}

template<class T> inline void Tensor<T>::set(unsigned int position, T element)
{
#ifdef DEBUG
 if (position > this->size-1)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> inline void Tensor<T>::" <<
          "set(unsigned int position, T element): " <<
          "(position > this->size-1)." << endl;
  exit(1);
 }
#endif
 this->Elements[position] = element;
}

template<class T> inline T Tensor<T>::get(const vector<unsigned int>& Index) const
{
#ifdef DEBUG
 if (Index.size() != this->rank)
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline T Tensor<T>::" <<
          "get(const vector<unsigned int>& Index)" << endl;
  exit(1);
 }
 for (int i = 0; i < this->rank; i++)
 {
  if (Index[i] >= this->Shape[i])
  {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline T Tensor<T>::" <<
          "get(const vector<unsigned int>& Index)" << endl;
  exit(1);
  }
 }
#endif
 unsigned int position = this->getPosition(Index);
 return this->Elements[position];
}

template<class T> inline T Tensor<T>::get(unsigned int position) const
{
#ifdef DEBUG
 if (position > this->size-1)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> inline T Tensor<T>::" <<
          "get(unsigned int position) const: " <<
          "(position > this->size-1)." << endl;
  exit(1);
 }
#endif
 return this->Elements[position];
}

template<class T> void Tensor<T>::write(const string& FileName) const
{
 ofstream File(FileName.c_str(), ios::out | ios::binary);
 if (File.is_open())
 {
  File.write((char*)&(this->rank), sizeof(this->rank));
  for (int i = 0; i < this->rank; i++)
  {
   File.write((char*)&(this->Shape[i]), sizeof(this->Shape[i]));
  }
  File.write((char*)&(this->size), sizeof(this->size));
  for (int i = 0; i < this->size; i++)
  {
   File.write((char*)&(this->Elements[i]), sizeof(this->Elements[i]));
  }
  File.close();
 }
 else
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Tensor<T>::" <<
          "write(const string& FileName) const: " <<
          "Binary output file " << FileName << " could not be opened." << endl;
  exit(1);
 }
}

template<class T> void Tensor<T>::read(const string& FileName)
{
 ifstream File(FileName.c_str(), ios::in | ios::binary);
 if (File.is_open())
 {
  File.read((char*)&(this->rank), sizeof(this->rank));
  this->Shape = vector<unsigned int>(this->rank);
  for (int i = 0; i < this->rank; i++)
  {
   File.read((char*)&(this->Shape[i]), sizeof(this->Shape[i]));
  }
  File.read((char*)&(this->size), sizeof(this->size));
  delete[] this->Elements;
  this->Elements = new T[this->size];
  for (int i = 0; i < this->size; i++)
  {
   File.read((char*)&(this->Elements[i]), sizeof(this->Elements[i]));
  }
  File.close();
 }
 else
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Tensor<T>::" <<
          "read(const string& FileName): " <<
          "Binary input file " << FileName << " could not be opened." << endl;
  exit(1);
 }
}

template<class T> T Tensor<T>::normalize()
{
#ifdef DEBUG
 if (this->size == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T Tensor<T>::" <<
          "normalize(): " <<
          "(this->size == 0)." << endl;
  exit(1);
 }
#endif
 T maxElement = 0.0;
 for (int position = 0; position < this->size; position++)
 {
  if (abs(this->Elements[position]) > abs(maxElement))
   maxElement = this->Elements[position];
 }
 for (int position = 0; position < this->size; position++)
  this->Elements[position] /= abs(maxElement);
 return abs(maxElement);
}

template<class T> T Tensor<T>::frobeniusNormalize()
{
#ifdef DEBUG
 if (this->size == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T Tensor<T>::" <<
          "frobeniusNormalize(): " <<
          "(this->size == 0)." << endl;
  exit(1);
 }
#endif
 T frobeniusNorm = 0.0;
 for (int position = 0; position < this->size; position++)
  frobeniusNorm += this->Elements[position]*MathAuxiliary::complexConjugate(this->Elements[position]);
 frobeniusNorm = sqrt(frobeniusNorm);
 for (int position = 0; position < this->size; position++)
  this->Elements[position] /= frobeniusNorm;
 return frobeniusNorm;
}

template<class T> void Tensor<T>::convertToComplex(Tensor< complex<float> >& Tensor0) const
{
#ifdef DEBUG
 if (typeid(T) != typeid(float))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Tensor<T>::" <<
          "convertToComplex(Tensor< complex<float> >& Tensor0) const: " <<
          "(typeid(T) != typeid(float))." << endl;
  exit(1);
 }
#endif
 Tensor0 = Tensor< complex<float> >(this->Shape);
 for (int position = 0; position < this->size; position++)
 {
  Tensor0.set(position, complex<float>(this->Elements[position]));
 }
}

template<class T> void Tensor<T>::convertToComplex(Tensor< complex<double> >& Tensor0) const
{
#ifdef DEBUG
 if (typeid(T) != typeid(double))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Tensor<T>::" <<
          "convertToComplex(Tensor< complex<double> >& Tensor0) const: " <<
          "(typeid(T) != typeid(double))." << endl;
  exit(1);
 }
#endif
 Tensor0 = Tensor< complex<double> >(this->Shape);
 for (int position = 0; position < this->size; position++)
 {
  Tensor0.set(position, complex<double>(this->Elements[position]));
 }
}

template<class T> void Tensor<T>::getSubtensor(const vector<unsigned int>& Indices,
                                               const vector<unsigned int>& Index,
                                               Tensor<T>& Tensor0) const
{
#ifdef DEBUG
 if (Indices.size()+Index.size() != this->rank)
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void Tensor<T>::" <<
          "getSubtensor(const vector<unsigned int>& Indices, " <<
                       "const vector<unsigned int>& Index, " <<
                       "Tensor<T>& Tensor0) const" << endl;
  exit(1);
 }
 for (int i = 0; i < Indices.size(); i++)
 {
  if (Indices[i] >= this->rank)
  {
   cerr << "Program terminated because of error in function: " <<
           "template<class T> void Tensor<T>::" <<
           "getSubtensor(const vector<unsigned int>& Indices, " <<
                        "const vector<unsigned int>& Index, " <<
                        "Tensor<T>& Tensor0) const" << endl;
   exit(1);
  }
  for (int j = i+1; j < Indices.size(); j++)
  {
   if (Indices[i] == Indices[j])
   {
    cerr << "Program terminated because of error in function: " <<
            "template<class T> void Tensor<T>::" <<
            "getSubtensor(const vector<unsigned int>& Indices, " <<
                         "const vector<unsigned int>& Index, " <<
                         "Tensor<T>& Tensor0) const" << endl;
    exit(1);
   }
  }
 }
#endif
 unsigned int rank0 = Indices.size();
 vector<unsigned int> Shape0(rank0);
 for (int i = 0; i < rank0; i++)
 {
  Shape0[i] = this->Shape[Indices[i]];
 }
 Tensor0 = Tensor<T>(Shape0);
 vector<unsigned int> Index0(this->rank), Index1(Tensor0.rank);
 bool flag; unsigned int k = 0;
 for (int i = 0; i < this->rank; i++)
 {
  flag = false;
  for (int j = 0; j < Indices.size(); j++)
  {
   if (i == Indices[j])
   {
    flag = true;
   }
  }
  if (flag == false)
  {
   Index0[i] = Index[k++];
  }
 }
 for (int i = 0; i < Tensor0.size; i++)
 {
  Tensor0.getIndex(i, Index1);
  for (int j = 0; j < Indices.size(); j++)
  {
   Index0[Indices[j]] = Index1[j];
  }
  Tensor0.Elements[i] = this->get(Index0);
 }
}

template<class T> void Tensor<T>::permute(const vector<unsigned int>& Order)
{
#ifdef DEBUG
 if (Order.size() != this->rank)
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void Tensor<T>::" <<
          "permute(const vector<unsigned int>& Order)" << endl;
  exit(1);
 }
 for (int i = 0; i < this->rank; i++)
 {
  if (Order[i] >= this->rank)
  {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void Tensor<T>::" <<
          "permute(const vector<unsigned int>& Order)" << endl;
  exit(1);
  }
  for (int j = i+1; j < this->rank; j++)
  {
   if (Order[i] == Order[j])
   {
   cerr << "Program terminated because of error in function: " <<
           "template<class T> void Tensor<T>::" <<
           "permute(const vector<unsigned int>& Order)" << endl;
   exit(1);
   }
  }
 }
#endif
 vector<unsigned int> NewShape(this->rank);
 for (int i = 0; i < this->rank; i++)
 {
  NewShape[i] = this->Shape[Order[i]];
 }
 T* NewElements = new T[this->size];
 unsigned int position, newPosition;
 vector<unsigned int> Index(this->rank), NewIndex(this->rank);
 for (int i = 0; i < this->size; i++)
 {
  position = i;
  this->getIndex(position, Index);
  for (int j = 0; j < this->rank; j++)
  {
   NewIndex[j] = Index[Order[j]];
  }
  newPosition = Tensor<T>::getPosition(NewShape, NewIndex);
  NewElements[newPosition] = this->Elements[position];
 }
 this->Shape = NewShape;
 delete[] this->Elements;
 this->Elements = NewElements;
}

template<class T> void Tensor<T>::reshape(const vector<unsigned int>& Shape0)
{
 unsigned int rank0 = Shape0.size();
#ifdef DEBUG
 unsigned int size0 = 1;
 for (int i = 0; i < rank0; i++)
 {
  size0 *= Shape0[i];
 }
 if (size0 != this->size)
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void Tensor<T>::" <<
          "reshape(const vector<unsigned int>& Shape0)" << endl;
  exit(1);
 }
#endif
 this->rank = rank0;
 this->Shape = Shape0;
}

template<class T> void Tensor<T>::multiply(T element)
{
#ifdef DEBUG
 if (this->size == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Tensor<T>::" <<
          "multiply(T element): " <<
          "(this->size == 0)." << endl;
  exit(1);
 }
#endif
 for (int i = 0; i < this->size; i++)
 {
  this->Elements[i] *= element;
 }
}

template<class T> void Tensor<T>::add(const Tensor<T>& Tensor0)
{
#ifdef DEBUG
 if (this->Shape != Tensor0.Shape)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Tensor<T>::" <<
          "add(const Tensor<T>& Tensor0): " <<
          "(this->Shape != Tensor0.Shape)." << endl;
  exit(1);
 }
#endif
 for (int i = 0; i < this->size; i++)
 {
  this->Elements[i] += Tensor0.Elements[i];
 }
}

template<class T> T Tensor<T>::scalarProduct() const
{
#ifdef DEBUG
 if (this->size == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> T Tensor<T>::" <<
          "scalarProduct() const: " <<
          "(this->size == 0)." << endl;
  exit(1);
 }
#endif
 T scalarProduct = 0.0;
 for (int position = 0; position < this->size; position++)
  scalarProduct += this->Elements[position]*MathAuxiliary::complexConjugate(this->Elements[position]);
 return scalarProduct;
}

template<class T> void Tensor<T>::contract(unsigned int index0, unsigned int index1)
{
#ifdef DEBUG
 if ((index0 >= this->rank) || (index1 >= this->rank) || (index0 == index1) ||
     (this->Shape[index0] != this->Shape[index1]))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void Tensor<T>::" <<
          "contract(unsigned int index0, unsigned int index1)" << endl;
  exit(1);
 }
#endif
 unsigned int newRank = this->rank - 2;
 vector<unsigned int> NewShape(newRank);
 unsigned int j = 0;
 for (int i = 0; i < this->rank; i++)
 {
  if ((i != index0) && (i != index1))
  {
   NewShape[j++] = this->Shape[i];
  }
 }
 unsigned int newSize = 1;
 for (int i = 0; i < newRank; i++)
 {
  newSize *= NewShape[i];
 }
 T* NewElements = new T[newSize];
 unsigned int newPosition, k;
 vector<unsigned int> NewIndex(newRank), Index(this->rank);
 T element;
 for (int i = 0; i < newSize; i++)
 {
  newPosition = i;
  Tensor<T>::getIndex(NewShape, newPosition, NewIndex);
  k = 0;
  for (int j = 0; j < this->rank; j++)
  {
   if ((j != index0) && (j != index1))
   {
    Index[j] = NewIndex[k++];
   }
  }
  element = 0.0;
  for (int j = 0; j < this->Shape[index0]; j++)
  {
   Index[index0] = j;
   Index[index1] = j;
   element += this->get(Index);
  }
  NewElements[newPosition] = element;
 }
 this->rank = newRank;
 this->Shape = NewShape;
 this->size = newSize;
 delete[] this->Elements;
 this->Elements = NewElements;
}

template<class T> void Tensor<T>::contract(const vector<unsigned int>& Indices,
                                           Tensor<T>& Tensor0,
                                           const vector<unsigned int>& Indices0)
{
 unsigned int numIndices = Indices.size();
 unsigned int numIndices0 = Indices0.size();
 unsigned int rank0 = Tensor0.getRank();
 vector<unsigned int> Shape0;
 Tensor0.getShape(Shape0);
#ifdef DEBUG
 unsigned int numElements = 1, numElements0 = 1;
 for (int i = 0; i < numIndices; i++)
 {
  numElements *= this->Shape[Indices[i]];
 }
 for (int i = 0; i < numIndices0; i++)
 {
  numElements0 *= Shape0[Indices0[i]];
 }
 if ((numIndices != numIndices0) || (numElements != numElements0))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void Tensor<T>::" <<
          "contract(const vector<unsigned int>& Indices, Tensor<T>& Tensor0, " <<
          "const vector<unsigned int>& Indices0)" << endl;
  exit(1);
 }
 for (int i = 0; i < numIndices; i++)
 {
  if ((Indices[i] >= this->rank) || (Indices0[i] >= rank0))
  {
   cerr << "Program terminated because of error in function: " <<
           "template<class T> void Tensor<T>::" <<
           "contract(const vector<unsigned int>& Indices, Tensor<T>& Tensor0, " <<
           "const vector<unsigned int>& Indices0)" << endl;
   exit(1);
  }
  for (int j = i+1; j < numIndices; j++)
  {
   if ((Indices[i] == Indices[j]) || (Indices0[i] == Indices0[j]))
   {
    cerr << "Program terminated because of error in function: " <<
            "template<class T> void Tensor<T>::" <<
            "contract(const vector<unsigned int>& Indices, Tensor<T>& Tensor0, " <<
            "const vector<unsigned int>& Indices0)" << endl;
    exit(1);
   }
  }
 }
#endif
 unsigned int restIndices = this->rank - numIndices;
 unsigned int restIndices0 = rank0 - numIndices;
 unsigned int finalRank = restIndices + restIndices0;
 vector<unsigned int> FinalShape(finalRank);
 vector<unsigned int> Order(this->rank);
 int k = 0, l = restIndices, m = 1, n = 1, o, p = 0;
 for (int i = 0; i < this->rank; i++)
 {
  o = 0;
  for (int j = 0; j < numIndices; j++)
  {
   if (Indices[j] == i)
    o = 1;
  }
  if (o == 0)
  {
   Order[k++] = i;
   m *= this->Shape[i];
   FinalShape[p++] = this->Shape[i];
  }
 }
 for (int i = 0; i < numIndices; i++)
 {
  Order[l++] = Indices[i];
  n *= this->Shape[Indices[i]];
 }
 this->permute(Order);
 unsigned int newRank = 2;
 vector<unsigned int> NewShape(newRank);
 NewShape[0] = m;
 NewShape[1] = n;
 this->reshape(NewShape);
 vector<unsigned int> Order0(rank0);
 int k0 = 0, l0 = numIndices, m0 = 1, n0 = 1;
 for (int i = 0; i < rank0; i++)
 {
  o = 0;
  for (int j = 0; j < numIndices; j++)
  {
   if (Indices0[j] == i)
   {
    o = 1;
   }
  }
  if (o == 0)
  {
   Order0[l0++] = i;
   n0 *= Shape0[i];
   FinalShape[p++] = Shape0[i];
  }
 }
 for (int i = 0; i < numIndices; i++)
 {
  Order0[k0++] = Indices0[i];
  m0 *= Shape0[Indices0[i]];
 }
 Tensor0.permute(Order0);
 unsigned int newRank0 = 2;
 vector<unsigned int> newShape0(newRank0);
 newShape0[0] = m0;
 newShape0[1] = n0;
 Tensor0.reshape(newShape0);
 char transa = 'N', transb = 'N';
 unsigned int finalSize = m * n0;
 T* FinalElements = new T[finalSize];
 if (typeid(T) == typeid(float))
 {
  float alpha = 1.0, beta = 0.0;
  sgemm_(&transa, &transb, &m, &n0, &n, &alpha, (float*)this->Elements, &m,
         (float*)Tensor0.Elements, &n, &beta, (float*)FinalElements, &m);
 }
 else if (typeid(T) == typeid(double))
 {
  double alpha = 1.0, beta = 0.0;
  dgemm_(&transa, &transb, &m, &n0, &n, &alpha, (double*)this->Elements, &m,
         (double*)Tensor0.Elements, &n, &beta, (double*)FinalElements, &m);
 }
 else if (typeid(T) == typeid(complex<float>))
 {
  complex<float> alpha = 1.0, beta = 0.0;
  cgemm_(&transa, &transb, &m, &n0, &n, &alpha, (complex<float>*)this->Elements, &m,
         (complex<float>*)Tensor0.Elements, &n, &beta, (complex<float>*)FinalElements, &m);
 }
 else if (typeid(T) == typeid(complex<double>))
 {
  complex<double> alpha = 1.0, beta = 0.0;
  zgemm_(&transa, &transb, &m, &n0, &n, &alpha, (complex<double>*)this->Elements, &m,
         (complex<double>*)Tensor0.Elements, &n, &beta, (complex<double>*)FinalElements, &m);
 }
 this->rank = finalRank;
 this->Shape = FinalShape;
 this->size = finalSize;
 delete[] this->Elements;
 this->Elements = FinalElements;
}

template<class T> void Tensor<T>::contract(const vector<unsigned int>& Indices,
                                           const Tensor<T>& Tensor0,
                                           const vector<unsigned int>& Indices0,
                                           Tensor<T>& Tensor1) const
{
 Tensor<T> thisC(*this), Tensor0C(Tensor0);
 unsigned int numIndices = Indices.size();
 unsigned int numIndices0 = Indices0.size();
 unsigned int rankC = thisC.getRank();
 unsigned int rank0C = Tensor0C.getRank();
 vector<unsigned int> ShapeC, Shape0C;
 thisC.getShape(ShapeC);
 Tensor0C.getShape(Shape0C);
#ifdef DEBUG
 unsigned int numElements = 1, numElements0 = 1;
 for (int i = 0; i < numIndices; i++)
 {
  numElements *= ShapeC[Indices[i]];
 }
 for (int i = 0; i < numIndices0; i++)
 {
  numElements0 *= Shape0C[Indices0[i]];
 }
 if ((numIndices != numIndices0) || (numElements != numElements0))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void Tensor<T>::" <<
          "contract(const vector<unsigned int>& Indices, const Tensor<T>& Tensor0, " <<
          "const vector<unsigned int>& Indices0, Tensor<T>& Tensor1) const" << endl;
  exit(1);
 }
 for (int i = 0; i < numIndices; i++)
 {
  if ((Indices[i] >= rankC) || (Indices0[i] >= rank0C))
  {
   cerr << "Program terminated because of error in function: " <<
           "template<class T> void Tensor<T>::" <<
           "contract(const vector<unsigned int>& Indices, const Tensor<T>& Tensor0, " <<
           "const vector<unsigned int>& Indices0, Tensor<T>& Tensor1) const" << endl;
   exit(1);
  }
  for (int j = i+1; j < numIndices; j++)
  {
   if ((Indices[i] == Indices[j]) || (Indices0[i] == Indices0[j]))
   {
    cerr << "Program terminated because of error in function: " <<
            "template<class T> void Tensor<T>::" <<
            "contract(const vector<unsigned int>& Indices, const Tensor<T>& Tensor0, " <<
            "const vector<unsigned int>& Indices0, Tensor<T>& Tensor1) const" << endl;
    exit(1);
   }
  }
 }
#endif
 unsigned int restIndices = rankC - numIndices;
 unsigned int restIndices0 = rank0C - numIndices;
 unsigned int rank1 = restIndices + restIndices0;
 vector<unsigned int> Shape1(rank1);
 vector<unsigned int> Order(rankC);
 int k = 0, l = restIndices, m = 1, n = 1, o, p = 0;
 for (int i = 0; i < rankC; i++)
 {
  o = 0;
  for (int j = 0; j < numIndices; j++)
  {
   if (Indices[j] == i)
    o = 1;
  }
  if (o == 0)
  {
   Order[k++] = i;
   m *= ShapeC[i];
   Shape1[p++] = ShapeC[i];
  }
 }
 for (int i = 0; i < numIndices; i++)
 {
  Order[l++] = Indices[i];
  n *= ShapeC[Indices[i]];
 }
 thisC.permute(Order);
 unsigned int newRankC = 2;
 vector<unsigned int> NewShapeC(newRankC);
 NewShapeC[0] = m;
 NewShapeC[1] = n;
 thisC.reshape(NewShapeC);
 vector<unsigned int> Order0(rank0C);
 int k0 = 0, l0 = numIndices, m0 = 1, n0 = 1;
 for (int i = 0; i < rank0C; i++)
 {
  o = 0;
  for (int j = 0; j < numIndices; j++)
  {
   if (Indices0[j] == i)
   {
    o = 1;
   }
  }
  if (o == 0)
  {
   Order0[l0++] = i;
   n0 *= Shape0C[i];
   Shape1[p++] = Shape0C[i];
  }
 }
 for (int i = 0; i < numIndices; i++)
 {
  Order0[k0++] = Indices0[i];
  m0 *= Shape0C[Indices0[i]];
 }
 Tensor0C.permute(Order0);
 unsigned int newRank0C = 2;
 vector<unsigned int> newShape0C(newRank0C);
 newShape0C[0] = m0;
 newShape0C[1] = n0;
 Tensor0C.reshape(newShape0C);
 char transa = 'N', transb = 'N';
 unsigned int size1 = m * n0;
 T* Elements1 = new T[size1];
 if (typeid(T) == typeid(float))
 {
  float alpha = 1.0, beta = 0.0;
  sgemm_(&transa, &transb, &m, &n0, &n, &alpha, (float*)thisC.Elements, &m,
         (float*)Tensor0C.Elements, &n, &beta, (float*)Elements1, &m);
 }
 else if (typeid(T) == typeid(double))
 {
  double alpha = 1.0, beta = 0.0;
  dgemm_(&transa, &transb, &m, &n0, &n, &alpha, (double*)thisC.Elements, &m,
         (double*)Tensor0C.Elements, &n, &beta, (double*)Elements1, &m);
 }
 else if (typeid(T) == typeid(complex<float>))
 {
  complex<float> alpha = 1.0, beta = 0.0;
  cgemm_(&transa, &transb, &m, &n0, &n, &alpha, (complex<float>*)thisC.Elements, &m,
         (complex<float>*)Tensor0C.Elements, &n, &beta, (complex<float>*)Elements1, &m);
 }
 else if (typeid(T) == typeid(complex<double>))
 {
  complex<double> alpha = 1.0, beta = 0.0;
  zgemm_(&transa, &transb, &m, &n0, &n, &alpha, (complex<double>*)thisC.Elements, &m,
         (complex<double>*)Tensor0C.Elements, &n, &beta, (complex<double>*)Elements1, &m);
 }
 Tensor1.rank = rank1;
 Tensor1.Shape = Shape1;
 Tensor1.size = size1;
 delete[] Tensor1.Elements;
 Tensor1.Elements = Elements1;
}

template<class T> void Tensor<T>::fillRandomly(unsigned int seed, T element)
{
#ifdef DEBUG
 if (this->size == 0)
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Tensor<T>::" <<
          "fillRandomly(unsigned int seed, T element): " <<
          "(this->size == 0)." << endl;
  exit(1);
 }
#endif
// random numbers uniformly distributed between [-1, 1]:
 int idist = 2;
 if (seed == 0)
  srand(time(0));
 else
  srand(seed);
// randomNums have to lie in [0, 4095] and randomNum3 has to be odd:
 int randomNum0 = rand() % 4096;
 int randomNum1 = rand() % 4096;
 int randomNum2 = rand() % 4096;
 int randomNum3 = ((rand() % 2048) * 2) + 1;
 int Iseed[] = {randomNum0, randomNum1, randomNum2, randomNum3};
 int n = this->size;
 if (typeid(T) == typeid(float))
  slarnv_(&idist, Iseed, &n, (float*)this->Elements);
 else if (typeid(T) == typeid(double))
  dlarnv_(&idist, Iseed, &n, (double*)this->Elements);
 else if (typeid(T) == typeid(complex<float>))
  clarnv_(&idist, Iseed, &n, (complex<float>*)this->Elements);
 else if (typeid(T) == typeid(complex<double>))
  zlarnv_(&idist, Iseed, &n, (complex<double>*)this->Elements);
 for (int i = 0; i < this->size; i++)
  this->Elements[i] *= element;
}

template<class T> void Tensor<T>::fillZeroes()
{
#ifdef DEBUG
 if (this->size == 0)
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void Tensor<T>::" <<
          "fillZeroes()" << endl;
  exit(1);
 }
#endif
 for (int i = 0; i < this->size; i++)
 {
  this->Elements[i] = 0.0;
 }
}

template<class T> void Tensor<T>::complexConjugate()
{
#ifdef DEBUG
 if (this->size == 0)
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void Tensor<T>::" <<
          "complexConjugate()" << endl;
  exit(1);
 }
#endif
 for (int i = 0; i < this->size; i++)
 {
  this->Elements[i] = MathAuxiliary::complexConjugate(this->Elements[i]);
 }
}

template<class T> void Tensor<T>::complexConjugate(Tensor<T>& Tensor0) const
{
#ifdef DEBUG
 if (this->size == 0)
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> void Tensor<T>::" <<
          "complexConjugate(Tensor<T>& Tensor0) const" << endl;
  exit(1);
 }
#endif
 Tensor0.rank = this->rank;
 Tensor0.Shape = this->Shape;
 Tensor0.size = this->size;
 delete[] Tensor0.Elements;
 Tensor0.Elements = new T[this->size];
 for (int i = 0; i < this->size; i++)
 {
  Tensor0.Elements[i] = MathAuxiliary::complexConjugate(this->Elements[i]);
 }
}

template<class T> void Tensor<T>::singularValueDecompose(const vector<unsigned int>& Indices0,
                                                         const vector<unsigned int>& Indices1,
                                                         unsigned int Dcut,
                                                         Tensor<T>& Tensor0, Matrix<T>& Sigma,
                                                         Tensor<T>& Tensor1)
{
 unsigned int rank0 = Indices0.size(), rank1 = Indices1.size();
#ifdef DEBUG
 if ((this->size == 0) || (this->rank != rank0+rank1))
 {
  cerr << "Program terminated because of error in function " <<
          "singularValueDecompose(const vector<unsigned int>& Indices0, " <<
                                 "const vector<unsigned int>& Indices1, " <<
                                 "unsigned int Dcut, " <<
                                 "Tensor<T>& Tensor0, Matrix<T>& Sigma, " <<
                                 "Tensor<T>& Tensor1): " <<
          "((this->size == 0) || (this->rank != Indices0.size()+Indices1.size()))." << endl;
  exit(1);
 }
 for (int i = 0; i < rank0; i++)
 {
  if (Indices0[i] >= this->rank)
  {
   cerr << "Program terminated because of error in function " <<
           "singularValueDecompose(const vector<unsigned int>& Indices0, " <<
                                  "const vector<unsigned int>& Indices1, " <<
                                  "unsigned int Dcut, " <<
                                  "Tensor<T>& Tensor0, Matrix<T>& Sigma, " <<
                                  "Tensor<T>& Tensor1): " <<
           "(Indices0[" << i << "] >= this->rank)." << endl;
   exit(1);
  }
  for (int j = i+1; j < rank0; j++)
  {
   if (Indices0[j] <= Indices0[i])
   {
    cerr << "Program terminated because of error in function " <<
            "singularValueDecompose(const vector<unsigned int>& Indices0, " <<
                                   "const vector<unsigned int>& Indices1, " <<
                                   "unsigned int Dcut, " <<
                                   "Tensor<T>& Tensor0, Matrix<T>& Sigma, " <<
                                   "Tensor<T>& Tensor1): " <<
            "(Indices0[" << j << "] <= Indices0[" << i << "])." << endl;
    exit(1);
   }
  }
  for (int j = 0; j < rank1; j++)
  {
   if (Indices1[j] == Indices0[i])
   {
    cerr << "Program terminated because of error in function " <<
            "singularValueDecompose(const vector<unsigned int>& Indices0, " <<
                                   "const vector<unsigned int>& Indices1, " <<
                                   "unsigned int Dcut, " <<
                                   "Tensor<T>& Tensor0, Matrix<T>& Sigma, " <<
                                   "Tensor<T>& Tensor1): " <<
            "(Indices1[" << j << "] == Indices0[" << i << "])." << endl;
    exit(1);
   }
  }
 }
 for (int i = 0; i < rank1; i++)
 {
  if (Indices1[i] >= this->rank)
  {
   cerr << "Program terminated because of error in function " <<
           "singularValueDecompose(const vector<unsigned int>& Indices0, " <<
                                  "const vector<unsigned int>& Indices1, " <<
                                  "unsigned int Dcut, " <<
                                  "Tensor<T>& Tensor0, Matrix<T>& Sigma, " <<
                                  "Tensor<T>& Tensor1): " <<
           "(Indices1[" << i << "] >= this->rank)." << endl;
   exit(1);
  }
  for (int j = i+1; j < rank1; j++)
  {
   if (Indices1[j] <= Indices1[i])
   {
    cerr << "Program terminated because of error in function " <<
            "singularValueDecompose(const vector<unsigned int>& Indices0, " <<
                                   "const vector<unsigned int>& Indices1, " <<
                                   "unsigned int Dcut, " <<
                                   "Tensor<T>& Tensor0, Matrix<T>& Sigma, " <<
                                   "Tensor<T>& Tensor1): " <<
            "(Indices1[" << j << "] <= Indices1[" << i << "])." << endl;
    exit(1);
   }
  }
 }
#endif
// get matrix dimensions for the SVD:
 unsigned int dim0 = 1, dim1 = 1;
 for (int i = 0; i < rank0; i++)
  dim0 *= this->Shape[Indices0[i]];
 for (int i = 0; i < rank1; i++)
  dim1 *= this->Shape[Indices1[i]];
// compute shape of Tensor0 and Tensor1:
 unsigned int dimRed = min(Dcut, min(dim0, dim1));
 vector<unsigned int> Shape0(rank0+1);
 for (int i = 0; i < rank0; i++)
  Shape0[i] = this->Shape[Indices0[i]];
 Shape0[rank0] = dimRed;
 vector<unsigned int> Shape1(rank1+1);
 Shape1[0] = dimRed;
 for (int i = 1; i < rank1+1; i++)
  Shape1[i] = this->Shape[Indices1[i-1]];
// permute and reshape this Tensor:
 vector<unsigned int> Order(this->rank);
 unsigned int position = 0;
 for (int i = 0; i < rank0; i++)
 {
  Order[position] = Indices0[i];
  position++;
 }
 for (int i = 0; i < rank1; i++)
 {
  Order[position] = Indices1[i];
  position++;
 }
 this->permute(Order);
 vector<unsigned int> Shape(2);
 Shape[0] = dim0; Shape[1] = dim1;
 this->reshape(Shape);
// construct pointer thisMatrixP and cast this Tensor to Matrix:
 Matrix<T>* thisMatrixP;
 thisMatrixP = static_cast<Matrix<T>*>(this);
// singular value decompose this Matrix:
 Matrix<T> U(dim0, dim0), Sigma0(dim0, dim1), Vt(dim1, dim1);
 thisMatrixP->singularValueDecompose(U, Sigma0, Vt);
// build Sigma, Tensor0 and Tensor1:
 Sigma = Matrix<T>(dimRed, dimRed);
 Sigma.fillZeroes();
 for (int i = 0; i < dimRed; i++)
  Sigma(i, i) = Sigma0(i, i);
 Shape[0] = dim0; Shape[1] = dimRed;
 Tensor0 = Tensor<T>(Shape);
 vector<unsigned int> Index(2);
 for (int j = 0; j < dimRed; j++)
 {
  Index[1] = j;
  for (int i = 0; i < dim0; i++)
  {
   Index[0] = i;
   Tensor0.set(Index, U(i, j));
  }
 }
 Tensor0.reshape(Shape0);
 Shape[0] = dimRed; Shape[1] = dim1;
 Tensor1 = Tensor<T>(Shape);
 for (int j = 0; j < dim1; j++)
 {
  Index[1] = j;
  for (int i = 0; i < dimRed; i++)
  {
   Index[0] = i;
   Tensor1.set(Index, Vt(i, j));
  }
 }
 Tensor1.reshape(Shape1);
}

template<class T> void Tensor<T>::QRDecompose(const vector<unsigned int>& Indices0,
                                              const vector<unsigned int>& Indices1,
                                              Tensor<T>& TensorR)
{
 unsigned int rank0 = Indices0.size(), rank1 = Indices1.size();
#ifdef DEBUG
 if ((this->size == 0) || (this->rank != rank0+rank1))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Tensor<T>::" <<
          "QRDecompose(const vector<unsigned int>& Indices0, " <<
                      "const vector<unsigned int>& Indices1, " <<
                      "Tensor<T>& TensorR): " <<
          "((this->size == 0) || (this->rank != Indices0.size()+Indices1.size()))." << endl;
  exit(1);
 }
 for (int i = 0; i < rank0; i++)
 {
  if (Indices0[i] >= this->rank)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Tensor<T>::" <<
           "QRDecompose(const vector<unsigned int>& Indices0, " <<
                       "const vector<unsigned int>& Indices1, " <<
                       "Tensor<T>& TensorR): " <<
           "(Indices0[" << i << "] >= this->rank)." << endl;
   exit(1);
  }
  for (int j = i+1; j < rank0; j++)
  {
   if (Indices0[j] <= Indices0[i])
   {
    cerr << "Program terminated because of error in function " <<
            "template<class T> void Tensor<T>::" <<
            "QRDecompose(const vector<unsigned int>& Indices0, " <<
                        "const vector<unsigned int>& Indices1, " <<
                        "Tensor<T>& TensorR): " <<
            "(Indices0[" << j << "] <= Indices0[" << i << "])." << endl;
    exit(1);
   }
  }
  for (int j = 0; j < rank1; j++)
  {
   if (Indices1[j] == Indices0[i])
   {
    cerr << "Program terminated because of error in function " <<
            "template<class T> void Tensor<T>::" <<
            "QRDecompose(const vector<unsigned int>& Indices0, " <<
                        "const vector<unsigned int>& Indices1, " <<
                        "Tensor<T>& TensorR): " <<
            "(Indices1[" << j << "] == Indices0[" << i << "])." << endl;
    exit(1);
   }
  }
 }
 for (int i = 0; i < rank1; i++)
 {
  if (Indices1[i] >= this->rank)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Tensor<T>::" <<
           "QRDecompose(const vector<unsigned int>& Indices0, " <<
                       "const vector<unsigned int>& Indices1, " <<
                       "Tensor<T>& TensorR): " <<
           "(Indices1[" << i << "] >= this->rank)." << endl;
   exit(1);
  }
  for (int j = i+1; j < rank1; j++)
  {
   if (Indices1[j] <= Indices1[i])
   {
    cerr << "Program terminated because of error in function " <<
            "template<class T> void Tensor<T>::" <<
            "QRDecompose(const vector<unsigned int>& Indices0, " <<
                        "const vector<unsigned int>& Indices1, " <<
                        "Tensor<T>& TensorR): " <<
            "(Indices1[" << j << "] <= Indices1[" << i << "])." << endl;
    exit(1);
   }
  }
 }
#endif
// get matrix dimensions:
 unsigned int dim0 = 1, dim1 = 1;
 for (int i = 0; i < rank0; i++)
  dim0 *= this->Shape[Indices0[i]];
 for (int i = 0; i < rank1; i++)
  dim1 *= this->Shape[Indices1[i]];
// compute new shape of this Tensor and TensorR:
 unsigned int dimRed = min(dim0, dim1);
 vector<unsigned int> Shape0(rank0+1);
 for (int i = 0; i < rank0; i++)
  Shape0[i] = this->Shape[Indices0[i]];
 Shape0[rank0] = dimRed;
 vector<unsigned int> Shape1(rank1+1);
 Shape1[0] = dimRed;
 for (int i = 1; i < rank1+1; i++)
  Shape1[i] = this->Shape[Indices1[i-1]];
// permute and reshape this Tensor:
 vector<unsigned int> Order(this->rank);
 unsigned int position = 0;
 for (int i = 0; i < rank0; i++)
 {
  Order[position] = Indices0[i];
  position++;
 }
 for (int i = 0; i < rank1; i++)
 {
  Order[position] = Indices1[i];
  position++;
 }
 this->permute(Order);
 vector<unsigned int> Shape(2);
 Shape[0] = dim0; Shape[1] = dim1;
 this->reshape(Shape);
// construct Matrix to this Tensor:
 Matrix<T> ThisMatrix(dim0, dim1);
 for (int i = 0; i < this->size; i++)
  ThisMatrix.Elements[i] = this->Elements[i];
// QR decompose this Matrix:
 Matrix<T> MatrixR;
 ThisMatrix.QRDecompose(MatrixR);
// reshape this Tensor and TensorR:
 *this = ThisMatrix;
 this->reshape(Shape0);
 TensorR = MatrixR;
 TensorR.reshape(Shape1);
}

template<class T> void Tensor<T>::LQDecompose(const vector<unsigned int>& Indices0,
                                              const vector<unsigned int>& Indices1,
                                              Tensor<T>& TensorL)
{
 unsigned int rank0 = Indices0.size(), rank1 = Indices1.size();
#ifdef DEBUG
 if ((this->size == 0) || (this->rank != rank0+rank1))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Tensor<T>::" <<
          "LQDecompose(const vector<unsigned int>& Indices0, " <<
                      "const vector<unsigned int>& Indices1, " <<
                      "Tensor<T>& TensorL): " <<
          "((this->size == 0) || (this->rank != Indices0.size()+Indices1.size()))." << endl;
  exit(1);
 }
 for (int i = 0; i < rank0; i++)
 {
  if (Indices0[i] >= this->rank)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Tensor<T>::" <<
           "LQDecompose(const vector<unsigned int>& Indices0, " <<
                       "const vector<unsigned int>& Indices1, " <<
                       "Tensor<T>& TensorL): " <<
           "(Indices0[" << i << "] >= this->rank)." << endl;
   exit(1);
  }
  for (int j = i+1; j < rank0; j++)
  {
   if (Indices0[j] <= Indices0[i])
   {
    cerr << "Program terminated because of error in function " <<
            "template<class T> void Tensor<T>::" <<
            "LQDecompose(const vector<unsigned int>& Indices0, " <<
                        "const vector<unsigned int>& Indices1, " <<
                        "Tensor<T>& TensorL): " <<
            "(Indices0[" << j << "] <= Indices0[" << i << "])." << endl;
    exit(1);
   }
  }
  for (int j = 0; j < rank1; j++)
  {
   if (Indices1[j] == Indices0[i])
   {
    cerr << "Program terminated because of error in function " <<
            "template<class T> void Tensor<T>::" <<
            "LQDecompose(const vector<unsigned int>& Indices0, " <<
                        "const vector<unsigned int>& Indices1, " <<
                        "Tensor<T>& TensorL): " <<
            "(Indices1[" << j << "] == Indices0[" << i << "])." << endl;
    exit(1);
   }
  }
 }
 for (int i = 0; i < rank1; i++)
 {
  if (Indices1[i] >= this->rank)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> void Tensor<T>::" <<
           "LQDecompose(const vector<unsigned int>& Indices0, " <<
                       "const vector<unsigned int>& Indices1, " <<
                       "Tensor<T>& TensorL): " <<
           "(Indices1[" << i << "] >= this->rank)." << endl;
   exit(1);
  }
  for (int j = i+1; j < rank1; j++)
  {
   if (Indices1[j] <= Indices1[i])
   {
    cerr << "Program terminated because of error in function " <<
            "template<class T> void Tensor<T>::" <<
            "LQDecompose(const vector<unsigned int>& Indices0, " <<
                        "const vector<unsigned int>& Indices1, " <<
                        "Tensor<T>& TensorL): " <<
            "(Indices1[" << j << "] <= Indices1[" << i << "])." << endl;
    exit(1);
   }
  }
 }
#endif
// get matrix dimensions:
 unsigned int dim0 = 1, dim1 = 1;
 for (int i = 0; i < rank0; i++)
  dim0 *= this->Shape[Indices0[i]];
 for (int i = 0; i < rank1; i++)
  dim1 *= this->Shape[Indices1[i]];
// compute new shape of this Tensor and TensorL:
 unsigned int dimRed = min(dim0, dim1);
 vector<unsigned int> Shape0(rank0+1);
 for (int i = 0; i < rank0; i++)
  Shape0[i] = this->Shape[Indices0[i]];
 Shape0[rank0] = dimRed;
 vector<unsigned int> Shape1(rank1+1);
 Shape1[0] = dimRed;
 for (int i = 1; i < rank1+1; i++)
  Shape1[i] = this->Shape[Indices1[i-1]];
// transpose, permute and reshape this Tensor:
 vector<unsigned int> Order(this->rank);
 unsigned int position = 0;
 for (int i = 0; i < rank1; i++)
 {
  Order[position] = Indices1[i];
  position++;
 }
 for (int i = 0; i < rank0; i++)
 {
  Order[position] = Indices0[i];
  position++;
 }
 this->permute(Order);
 vector<unsigned int> Shape(2);
 Shape[0] = dim1; Shape[1] = dim0;
 this->reshape(Shape);
// construct Matrix to this Tensor:
 Matrix<T> ThisMatrix(dim1, dim0);
 for (int i = 0; i < this->size; i++)
  ThisMatrix.Elements[i] = this->Elements[i];
// QR decompose this Matrix:
 Matrix<T> MatrixR;
 ThisMatrix.QRDecompose(MatrixR);
// transpose again, and reshape this Tensor and TensorL:
 MatrixR.transpose();
 TensorL = MatrixR;
 TensorL.reshape(Shape0);
 ThisMatrix.transpose();
 *this = ThisMatrix;
 this->reshape(Shape1);
}

template<class T> void Tensor<T>::twoBodyOperatorSVD(Tensor<T>& TensorLeft, Tensor<T>& TensorRight)
{
 unsigned int d0 = this->Shape[0], d1 = this->Shape[1];
 unsigned int numSV = min(d0*d0, d1*d1);
#ifdef DEBUG
 if ((this->rank != 4) || (this->Shape[0] != this->Shape[2]) || (this->Shape[1] != this->Shape[3]) ||
     (TensorLeft.Shape[0] != 1) || (TensorLeft.Shape[1] != numSV) || (TensorLeft.Shape[2] != d0) ||
     (TensorLeft.Shape[3] != d0) || (TensorRight.Shape[0] != numSV) || (TensorRight.Shape[1] != 1) ||
     (TensorRight.Shape[2] != d1) || (TensorRight.Shape[3] != d1))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> void Tensor<T>::" <<
          "twoBodyOperatorSVD(Tensor<T>& TensorLeft, Tensor<T>& TensorRight): " <<
          "This Tensor or the input arguments have an incorrect form." << endl;
  exit(1);
 }
#endif
 Matrix<T> ThisMatrix(d0*d0, d1*d1), U(d0*d0, d0*d0), Sigma(d0*d0, d1*d1), Vt(d1*d1, d1*d1);
 vector<unsigned int> Order(4);
 Order[0] = 0; Order[1] = 2; Order[2] = 1; Order[3] = 3;
 this->permute(Order);
 vector<unsigned int> Shape0(2);
 Shape0[0] = d0*d0; Shape0[1] = d1*d1;
 this->reshape(Shape0);
 vector<unsigned int> Index(2);
 for (int i = 0; i < d0*d0; i++)
 {
  Index[0] = i;
  for (int j = 0; j < d1*d1; j++)
  {
   Index[1] = j;
   ThisMatrix(i, j) = this->get(Index);
  }
 }
 ThisMatrix.singularValueDecompose(U, Sigma, Vt);
 vector<unsigned int> Index0(4);
 T element;
 Index0[0] = 0;
 for (int i = 0; i < d0*d0; i++)
 {
  Index0[2] = i % d0;
  Index0[3] = i / d0;
  for (int j = 0; j < numSV; j++)
  {
   Index0[1] = j;
   element = U(i, j) * sqrt(Sigma(j, j));
   TensorLeft.set(Index0, element);
  }
 }
 Index0[1] = 0;
 for (int i = 0; i < numSV; i++)
 {
  Index0[0] = i;
  for (int j = 0; j < d1*d1; j++)
  {
   Index0[2] = j % d1;
   Index0[3] = j / d1;
   element = sqrt(Sigma(i, i)) * Vt(i, j);
   TensorRight.set(Index0, element);
  }
 }
}

template<class T> double Tensor<T>::getNLDA(unsigned int range, unsigned int& rank, vector<double>& Sigma,
                                            vector< Tensor<T> >& NLDA, double cutoff, unsigned int mode) const
{
 unsigned int L = this->rank, S = this->Shape[0];
#ifdef DEBUG
 if ((this->size == 0) || ((mode != 0) && (mode != 1) && (mode != 2)))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> double Tensor<T>::" <<
          "getNLDA(unsigned int range, unsigned int& rank, vector<double>& Sigma, " <<
                  "vector< Tensor<T> >& NLDA, double cutoff, unsigned int mode) const: " <<
          "((this->size == 0) || ((mode != 0) && (mode != 1) && (mode != 2)))." << endl;
  exit(1);
 }
 for (int i = 0; i < L; i++)
 {
  if (this->Shape[i] != S)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> double Tensor<T>::" <<
           "getNLDA(unsigned int range, unsigned int& rank, vector<double>& Sigma, " <<
                   "vector< Tensor<T> >& NLDA, double cutoff, unsigned int mode) const: " <<
           "(this->Shape[" << i << "] != R+1)." << endl;
   exit(1);
  }
 }
#endif
 range = min(range, L-1);
 if (mode == 0)
 {
// number of training densities = number of non-zero entries in this Tensor =: m:
 unsigned int m = 0;
 for (int i = 0; i < this->size; i++)
 {
  if (this->Elements[i] != 0.0)
   m++;
 }
// number of NLDA parameters =: n:
 unsigned int n = L*S;
 for (int x = 1; x <= range; x++)
  n += (L-x)*S*S;
// define linear least squares problem:
 vector<T> F(m);
 Matrix<T> A(m, n); A.fillZeroes();
 unsigned int row = 0, xStart;
 vector<unsigned int> Index(L);
 for (int i = 0; i < this->size; i++)
 {
  if (this->Elements[i] != 0.0)
  {
// set b:
   F[row] = this->Elements[i];
// set A:
   this->getIndex(i, Index);
// - terms for range 0:
   for (int j = 0; j < L; j++)
    A(row, j*S+Index[j]) = 1.0;
// - terms for range x > 0:
   xStart = L*S;
   for (int x = 1; x <= range; x++)
   {
    for (int j = 0; j < L-x; j++)
    {
     A(row, xStart+j*S*S+Index[j]+S*Index[j+x]) = 1.0;
    }
    xStart += (L-x)*S*S;
   }
   row++;
  }
 }
// solve linear least squares:
 Matrix<T> b(m, 1);
 for (int i = 0; i < m; i++)
  b(i, 0) = F[i];
 Matrix<T> AC(A);
 AC.linearLeastSquares(cutoff, b, rank, Sigma);
 vector<T> Sol(n);
 for (int i = 0; i < n; i++)
  Sol[i] = b(i, 0);
// write NLDA:
 unsigned int sizeNLDA = L;
 for (int x = 1; x <= range; x++)
  sizeNLDA += (L-x);
 NLDA = vector< Tensor<T> >(sizeNLDA);
 unsigned int posNLDA = 0;
// - Tensors for range 0:
 vector<unsigned int> Shape0(1); Shape0[0] = S;
 Tensor<T> Tensor0(Shape0);
 for (int i = 0; i < L; i++)
 {
  for (int j = 0; j < S; j++)
   Tensor0.set(j, b(i*S+j, 0));
  NLDA[posNLDA] = Tensor0;
  posNLDA++;
 }
// - Tensors for range x > 0:
 Shape0 = vector<unsigned int>(2); Shape0[0] = S; Shape0[1] = S;
 Tensor0 = Tensor<T>(Shape0);
 xStart = L*S;
 for (int x = 1; x <= range; x++)
 {
  for (int i = 0; i < L-x; i++)
  {
   for (int j = 0; j < S*S; j++)
    Tensor0.set(j, b(xStart+i*S*S+j, 0));
   NLDA[posNLDA] = Tensor0;
   posNLDA++;
  }
  xStart += (L-x)*S*S;
 }
// compute distance:
 vector<T> ASol(m);
 A.multiply(Sol, ASol);
 double distance = 0.0;
 for (int i = 0; i < m; i++)
  distance += pow(ASol[i]-F[i], 2);
 distance = sqrt(abs(distance));
 double mean = 0.0;
 for (int i = 0; i < m; i++)
  mean += pow(F[i], 2);
 mean = sqrt(abs(mean));
 return distance/mean;
 }
 else if (mode == 1)
 {}
 else if (mode == 2)
 {}
}

template<class T> double Tensor<T>::getPolyNLDA(const vector< Tensor<T> >& nTensors, unsigned int range, unsigned int deg0,
                                                unsigned int deg1, vector<double>& Sigma, vector< Tensor<T> >& NLDA, double cutoff) const
{
// number of training densities = number of entries in this Tensor =: m:
 unsigned int m = this->size;
// system size L:
 unsigned int L = nTensors.size();
// number of NLDA parameters =: n:
 unsigned int n = L*(deg0+1);
 for (int x = 1; x <= range; x++)
  n += (L-x)*(deg1+1)*(deg1+1);
 vector<unsigned int> Shape0(1); Shape0[0] = m;
#ifdef DEBUG
 if ((m == 0) || (this->Shape != Shape0) || (m < n) || (L == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> double Tensor<T>::" <<
          "getPolyNLDA(const vector< Tensor<T> >& nTensors, unsigned int range, unsigned int deg0, " <<
                      "unsigned int deg1, vector<double>& Sigma, vector< Tensor<T> >& NLDA, double cutoff) const: " <<
          "((this->size == 0) || (this Tensor has wrong Shape) || (this->size < number of NLDA parameters) || (nTensors.size() == 0))." << endl;
  exit(1);
 }
 for (int i = 0; i < L; i++)
 {
  if (nTensors[i].Shape != Shape0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> double Tensor<T>::" <<
           "getPolyNLDA(const vector< Tensor<T> >& nTensors, unsigned int range, unsigned int deg0, " <<
                       "unsigned int deg1, vector<double>& Sigma, vector< Tensor<T> >& NLDA, double cutoff) const: " <<
           "(nTensors[" << i << "] has wrong Shape)." << endl;
   exit(1);
  }
 }
#endif
 range = min(range, L-1);
// define linear least squares problem:
 vector<T> F(m);
 Matrix<T> A(m, n); A.fillZeroes();
 unsigned int xStart;
 for (int i = 0; i < m; i++)
 {
// set b:
  F[i] = this->get(i);
// set A:
// - terms for range 0:
  for (int j = 0; j < L; j++)
  {
   for (int k = 0; k <= deg0; k++)
   {
    A(i, j*(deg0+1)+k) = pow(nTensors[j].get(i), k);
   }
  }
// - terms for range x > 0:
  xStart = L*(deg0+1);
  for (int x = 1; x <= range; x++)
  {
   for (int j = 0; j < L-x; j++)
   {
    for (int k2 = 0; k2 <= deg1; k2++){
    for (int k1 = 0; k1 <= deg1; k1++){
     A(i, xStart+j*(deg1+1)*(deg1+1)+k2*(deg1+1)+k1) = pow(nTensors[j].get(i), k1)*pow(nTensors[j+x].get(i), k2);
    }
    }
   }
   xStart += (L-x)*(deg1+1)*(deg1+1);
  }
 }
// solve linear least squares:
 Matrix<T> b(m, 1);
 for (int i = 0; i < m; i++)
  b(i, 0) = F[i];
 Matrix<T> AC(A);
 unsigned int rank;
 AC.linearLeastSquares(cutoff, b, rank, Sigma);
 vector<T> Sol(n);
 for (int i = 0; i < n; i++)
  Sol[i] = b(i, 0);
// write NLDA:
 unsigned int sizeNLDA = L;
 for (int x = 1; x <= range; x++)
  sizeNLDA += (L-x);
 NLDA = vector< Tensor<T> >(sizeNLDA);
 unsigned int posNLDA = 0;
// - Tensors for range x == 0:
 Shape0[0] = deg0+1;
 Tensor<T> Tensor0(Shape0);
 for (int i = 0; i < L; i++)
 {
  for (int j = 0; j <= deg0; j++)
   Tensor0.set(j, b(i*(deg0+1)+j, 0));
  NLDA[posNLDA] = Tensor0;
  posNLDA++;
 }
// - Tensors for range x > 0:
 Shape0 = vector<unsigned int>(2); Shape0[0] = deg1+1; Shape0[1] = deg1+1;
 Tensor0 = Tensor<T>(Shape0);
 xStart = L*(deg0+1);
 for (int x = 1; x <= range; x++)
 {
  for (int i = 0; i < L-x; i++)
  {
   for (int j = 0; j < (deg1+1)*(deg1+1); j++)
    Tensor0.set(j, b(xStart+i*(deg1+1)*(deg1+1)+j, 0));
   NLDA[posNLDA] = Tensor0;
   posNLDA++;
  }
  xStart += (L-x)*(deg1+1)*(deg1+1);
 }
// compute distance:
 vector<T> ASol(m);
 A.multiply(Sol, ASol);
 double distance = 0.0;
 for (int i = 0; i < m; i++)
  distance += pow(ASol[i]-F[i], 2);
 distance = sqrt(abs(distance));
 double mean = 0.0;
 for (int i = 0; i < m; i++)
  mean += pow(F[i], 2);
 mean = sqrt(abs(mean));
 return distance/mean;
}

template<class T> double Tensor<T>::getTaylorNLDA(const vector< Tensor<T> >& nTensors, unsigned int range,
                                                  unsigned int deg0, unsigned int deg1, vector<double>& Sigma,
                                                  vector< Tensor<T> >& NLDA, double cutoff) const
{
// number of training densities = number of entries in this Tensor =: m:
 unsigned int m = this->size;
// system size L:
 unsigned int L = nTensors.size();
// number of NLDA parameters =: n:
 unsigned int n = L*(deg0+1);
 for (int x = 1; x <= range; x++)
  n += (L-x)*deg1*deg1;
 vector<unsigned int> Shape0(1); Shape0[0] = m;
#ifdef DEBUG
 if ((m == 0) || (this->Shape != Shape0) || (m < n) || (L == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> double Tensor<T>::" <<
          "getTaylorNLDA(const vector< Tensor<T> >& nTensors, unsigned int range, unsigned int deg0, " <<
                        "unsigned int deg1, vector<double>& Sigma, vector< Tensor<T> >& NLDA, double cutoff) const: " <<
          "((this->size == 0) || (this Tensor has wrong Shape) || (this->size < number of NLDA parameters) || (nTensors.size() == 0))." << endl;
  exit(1);
 }
 for (int i = 0; i < L; i++)
 {
  if (nTensors[i].Shape != Shape0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> double Tensor<T>::" <<
           "getTaylorNLDA(const vector< Tensor<T> >& nTensors, unsigned int range, unsigned int deg0, " <<
                         "unsigned int deg1, vector<double>& Sigma, vector< Tensor<T> >& NLDA, double cutoff) const: " <<
           "(nTensors[" << i << "] has wrong Shape)." << endl;
   exit(1);
  }
 }
#endif
 range = min(range, L-1);
// define linear least squares problem:
 vector<T> F(m);
 Matrix<T> A(m, n); A.fillZeroes();
 unsigned int xStart;
 for (int i = 0; i < m; i++)
 {
// set b:
  F[i] = this->get(i);
// set A:
// - terms for range 0:
  for (int j = 0; j < L; j++)
  {
   for (int k = 0; k <= deg0; k++)
   {
    A(i, j*(deg0+1)+k) = pow(nTensors[j].get(i)-1.0, 2*k);
   }
  }
// - terms for range x > 0:
  xStart = L*(deg0+1);
  for (int x = 1; x <= range; x++)
  {
   for (int j = 0; j < L-x; j++)
   {
    for (int k2 = 1; k2 <= deg1; k2++){
    for (int k1 = 1; k1 <= deg1; k1++){
     A(i, xStart+j*deg1*deg1+(k2-1)*deg1+k1-1) = pow(nTensors[j].get(i)-1.0, k1)*pow(nTensors[j+x].get(i)-1.0, k2);
    }
    }
   }
   xStart += (L-x)*deg1*deg1;
  }
 }
// solve linear least squares:
 Matrix<T> b(m, 1);
 for (int i = 0; i < m; i++)
  b(i, 0) = F[i];
 Matrix<T> AC(A);
 unsigned int rank;
 AC.linearLeastSquares(cutoff, b, rank, Sigma);
 vector<T> Sol(n);
 for (int i = 0; i < n; i++)
  Sol[i] = b(i, 0);
// write NLDA:
 unsigned int sizeNLDA = L;
 for (int x = 1; x <= range; x++)
  sizeNLDA += (L-x);
 NLDA = vector< Tensor<T> >(sizeNLDA);
 unsigned int posNLDA = 0;
// - Tensors for range x == 0:
 Shape0[0] = deg0+1;
 Tensor<T> Tensor0(Shape0);
 for (int i = 0; i < L; i++)
 {
  for (int j = 0; j <= deg0; j++)
   Tensor0.set(j, b(i*(deg0+1)+j, 0));
  NLDA[posNLDA] = Tensor0;
  posNLDA++;
 }
// - Tensors for range x > 0:
 Shape0 = vector<unsigned int>(2); Shape0[0] = deg1; Shape0[1] = deg1;
 Tensor0 = Tensor<T>(Shape0);
 xStart = L*(deg0+1);
 for (int x = 1; x <= range; x++)
 {
  for (int i = 0; i < L-x; i++)
  {
   for (int j = 0; j < deg1*deg1; j++)
    Tensor0.set(j, b(xStart+i*deg1*deg1+j, 0));
   NLDA[posNLDA] = Tensor0;
   posNLDA++;
  }
  xStart += (L-x)*deg1*deg1;
 }
// compute distance:
 vector<T> ASol(m);
 A.multiply(Sol, ASol);
 double distance = 0.0;
 for (int i = 0; i < m; i++)
  distance += pow(ASol[i]-F[i], 2);
 distance = sqrt(abs(distance));
 double mean = 0.0;
 for (int i = 0; i < m; i++)
  mean += pow(F[i], 2);
 mean = sqrt(abs(mean));
 return distance/mean;
}

template<class T> double Tensor<T>::getTaylorNLDA2(const vector< Tensor<T> >& nTensors, unsigned int range,
                                                   unsigned int deg0, unsigned int deg1, vector<double>& Sigma,
                                                   vector< Tensor<T> >& NLDA, double cutoff) const
{
// number of training densities = number of entries in this Tensor =: m:
 unsigned int m = this->size;
// system size L:
 unsigned int L = nTensors.size();
// number of NLDA parameters =: n:
 unsigned int n = L*(deg0+1);
 unsigned int numParam = 0;
 for (int j = 1; j <= deg1; j++){
  for (int i = 1; i <= deg1; i++){
   if ((i+j)%2 == 0)
    numParam++;
  }
 }
 for (int x = 1; x <= range; x++)
  n += (L-x)*numParam;
 vector<unsigned int> Shape0(1); Shape0[0] = m;
#ifdef DEBUG
 if ((m == 0) || (this->Shape != Shape0) || (m < n) || (L == 0))
 {
  cerr << "Program terminated because of error in function " <<
          "template<class T> double Tensor<T>::" <<
          "getTaylorNLDA2(const vector< Tensor<T> >& nTensors, unsigned int range, unsigned int deg0, " <<
                         "unsigned int deg1, vector<double>& Sigma, vector< Tensor<T> >& NLDA, double cutoff) const: " <<
          "((this->size == 0) || (this Tensor has wrong Shape) || (this->size < number of NLDA parameters) || (nTensors.size() == 0))." << endl;
  exit(1);
 }
 for (int i = 0; i < L; i++)
 {
  if (nTensors[i].Shape != Shape0)
  {
   cerr << "Program terminated because of error in function " <<
           "template<class T> double Tensor<T>::" <<
           "getTaylorNLDA2(const vector< Tensor<T> >& nTensors, unsigned int range, unsigned int deg0, " <<
                          "unsigned int deg1, vector<double>& Sigma, vector< Tensor<T> >& NLDA, double cutoff) const: " <<
           "(nTensors[" << i << "] has wrong Shape)." << endl;
   exit(1);
  }
 }
#endif
 range = min(range, L-1);
// define linear least squares problem:
 vector<T> F(m);
 Matrix<T> A(m, n); A.fillZeroes();
 unsigned int xStart, count;
 for (int i = 0; i < m; i++)
 {
// set b:
  F[i] = this->get(i);
// set A:
// - terms for range 0:
  for (int j = 0; j < L; j++){
   for (int k = 0; k <= deg0; k++){
    A(i, j*(deg0+1)+k) = pow(nTensors[j].get(i)-1.0, 2*k);
   }
  }
// - terms for range x > 0:
  xStart = L*(deg0+1);
  for (int x = 1; x <= range; x++){
   for (int j = 0; j < L-x; j++){
    count = 0;
    for (int k2 = 1; k2 <= deg1; k2++){
    for (int k1 = 1; k1 <= deg1; k1++){
     if ((k1+k2)%2 == 0){
      A(i, xStart+j*numParam+count) = pow(nTensors[j].get(i)-1.0, k1)*pow(nTensors[j+x].get(i)-1.0, k2);
       count++;
     }
    }
    }
   }
   xStart += (L-x)*numParam;
  }
 }
// solve linear least squares:
 Matrix<T> b(m, 1);
 for (int i = 0; i < m; i++)
  b(i, 0) = F[i];
 Matrix<T> AC(A);
 unsigned int rank;
 AC.linearLeastSquares(cutoff, b, rank, Sigma);
 vector<T> Sol(n);
 for (int i = 0; i < n; i++)
  Sol[i] = b(i, 0);
// write NLDA:
 unsigned int sizeNLDA = L;
 for (int x = 1; x <= range; x++)
  sizeNLDA += (L-x);
 NLDA = vector< Tensor<T> >(sizeNLDA);
 unsigned int posNLDA = 0;
// - Tensors for range x == 0:
 Shape0[0] = deg0+1;
 Tensor<T> Tensor0(Shape0);
 for (int i = 0; i < L; i++)
 {
  for (int j = 0; j <= deg0; j++)
   Tensor0.set(j, b(i*(deg0+1)+j, 0));
  NLDA[posNLDA] = Tensor0;
  posNLDA++;
 }
// - Tensors for range x > 0:
 Shape0 = vector<unsigned int>(2); Shape0[0] = deg1; Shape0[1] = deg1;
 Tensor0 = Tensor<T>(Shape0);
 xStart = L*(deg0+1);
 for (int x = 1; x <= range; x++){
  for (int i = 0; i < L-x; i++){
   count = 0;
   for (int k2 = 1; k2 <= deg1; k2++){
   for (int k1 = 1; k1 <= deg1; k1++){
    if ((k1+k2)%2 == 0){
     Tensor0.set((k2-1)*deg1+k1-1, b(xStart+i*numParam+count, 0));
     count++;
    }
    else
     Tensor0.set((k2-1)*deg1+k1-1, 0.0);
   }
   }
   NLDA[posNLDA] = Tensor0;
   posNLDA++;
  }
  xStart += (L-x)*numParam;
 }
// compute distance:
 vector<T> ASol(m);
 A.multiply(Sol, ASol);
 double distance = 0.0;
 for (int i = 0; i < m; i++)
  distance += pow(ASol[i]-F[i], 2);
 distance = sqrt(abs(distance));
 double mean = 0.0;
 for (int i = 0; i < m; i++)
  mean += pow(F[i], 2);
 mean = sqrt(abs(mean));
 return distance/mean;
}

template<class T> void getMPS(Tensor<T>& Tensor0, MPS<T>& MPS0)
{
#ifdef DEBUG
 if (Tensor0.size == 0)
 {
  cerr << "Program terminated because of error in friend function " <<
          "template<class T> void " <<
          "getMPS(Tensor<T>& Tensor0, MPS<T>& MPS0): " <<
          "(Tensor0.size == 0)." << endl;
  exit(1);
 }
#endif
// add one index of range 1 on each boundary of Tensor0:
 unsigned int rank0 = Tensor0.rank+2;
 vector<unsigned int> Shape0(rank0);
 Shape0[0] = 1; Shape0[rank0-1] = 1;
 for (int l = 1; l < rank0-1; l++)
  Shape0[l] = Tensor0.Shape[l-1];
 Tensor0.rank = rank0;
 Tensor0.Shape = Shape0;
// define MPS0:
 string BC0 = "open";
 unsigned int N0 = rank0-2;
 unsigned int d0 = 0;
 for (int l = 1; l < rank0-1; l++)
 {
  if (Shape0[l] > d0)
   d0 = Shape0[l];
 }
 unsigned int dL = Shape0[1], dR = 1;
 for (int l = 2; l < rank0-1; l++)
  dR *= Shape0[l];
 unsigned int D0 = min(dL, dR);
 for (int l = 2; l < rank0-2; l++)
 {
  dL *= Shape0[l];
  dR /= Shape0[l];
  if (min(dL, dR) > D0)
   D0 = min(dL, dR);
 }
 MPS0 = MPS<T>(BC0, N0, d0, D0);
// perform successive SVDs:
 unsigned int Dcut = D0;
 vector<unsigned int> Indices1(2), Indices2, Order(3), IndexR(1), IndexL(1);
 Indices1[0] = 0; Indices1[1] = 1;
 Order[0] = 0; Order[1] = 2; Order[2] = 1;
 IndexR[0] = 1; IndexL[0] = 0;
 Tensor<T> Tensor1, Tensor2;
 Matrix<T> Sigma;
 for (int l = 2; l < rank0-1; l++)
 {
  Indices2 = vector<unsigned int>(rank0-l);
  for (int m = 0; m < rank0-l; m++)
   Indices2[m] = 2+m;
  Tensor0.singularValueDecompose(Indices1, Indices2, Dcut, Tensor1, Sigma, Tensor2);
  Tensor1.permute(Order);
  MPS0.Tensors[l-2] = Tensor1;
  Sigma.contract(IndexR, Tensor2, IndexL);
  Tensor0 = Sigma;
 }
 Tensor0.permute(Order);
 MPS0.set(rank0-3, Tensor0);
}

template<class T> inline unsigned int Tensor<T>::getPosition(const vector<unsigned int>& Index)
                                                 const
{
#ifdef DEBUG
 if (Index.size() != this->rank)
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline unsigned int Tensor<T>::" <<
          "getPosition(const vector<unsigned int>& Index) const" << endl;
  exit(1);
 }
 for (int i = 0; i < this->rank; i++)
 {
  if (Index[i] >= this->Shape[i])
  {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline unsigned int Tensor<T>::" <<
          "getPosition(const vector<unsigned int>& Index) const" << endl;
  exit(1);
  }
 }
#endif
 unsigned int position = 0, i = 1;
 for (int j = 0; j < this->rank; j++)
 {
  position += Index[j] * i;
  i *= this->Shape[j];
 }
 return position;
}

template<class T> inline unsigned int Tensor<T>::getPosition(const vector<unsigned int>& Shape0,
                                                             const vector<unsigned int>& Index0)
{
 unsigned int rank0 = Shape0.size();
#ifdef DEBUG
 if (rank0 != Index0.size())
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline unsigned int Tensor<T>::" <<
          "getPosition(const vector<unsigned int>& Shape0, " <<
          "const vector<unsigned int>& Index0)" << endl;
  exit(1);
 }
 for (int i = 0; i < rank0; i++)
 {
  if (Shape0[i] <= Index0[i])
  {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline unsigned int Tensor<T>::" <<
          "getPosition(const vector<unsigned int>& Shape0, " <<
          "const vector<unsigned int>& Index0)" << endl;
  exit(1);
  }
 }
#endif
 unsigned int position = 0, i = 1;
 for (int j = 0; j < rank0; j++)
 {
  position += Index0[j] * i;
  i *= Shape0[j];
 }
 return position;
}

template<class T> inline void Tensor<T>::getIndex(unsigned int position,
                                                  vector<unsigned int>& Index) const
{
#ifdef DEBUG
 if ((position >= this->size) || (Index.size() != this->rank))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline void Tensor<T>::" <<
          "getIndex(unsigned int position, vector<unsigned int>& Index) " <<
          "const" << endl;
  exit(1);
 }
#endif
 unsigned int ShapeProd[this->rank];
 ShapeProd[0] = 1;
 for (int i = 1; i < this->rank; i++)
 {
  ShapeProd[i] = ShapeProd[i-1] * this->Shape[i-1];
 }
 unsigned int rest = position;
 for (int i = this->rank-1; i >= 0; i--)
 {
  Index[i] = rest / ShapeProd[i];
  rest %= ShapeProd[i];
 }
}

template<class T> inline void Tensor<T>::getIndex(const vector<unsigned int>& Shape0,
                                                  unsigned int position0,
                                                  vector<unsigned int>& Index0)
{
 unsigned int rank0 = Shape0.size();
#ifdef DEBUG
 unsigned int size0 = 1;
 for (int i = 0; i < rank0; i++)
 {
  size0 *= Shape0[i];
 }
 if ((rank0 != Index0.size()) || (position0 >= size0))
 {
  cerr << "Program terminated because of error in function: " <<
          "template<class T> inline void Tensor<T>::" <<
          "getIndex(const vector<unsigned int>& Shape0, unsigned int position0, " <<
          "vector<unsigned int>& Index0)" << endl;
  exit(1);
 }
#endif
 unsigned int ShapeProd[rank0];
 ShapeProd[0] = 1;
 for (int i = 1; i < rank0; i++)
 {
  ShapeProd[i] = ShapeProd[i-1] * Shape0[i-1];
 }
 unsigned int rest = position0;
 for (int i = rank0-1; i >= 0; i--)
 {
  Index0[i] = rest / ShapeProd[i];
  rest %= ShapeProd[i];
 }
}
