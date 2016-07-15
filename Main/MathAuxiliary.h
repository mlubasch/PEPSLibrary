/// The required auxiliary math functions.
/** This header file contains the required auxiliary math functions.
    \remark This software is part of the PEPS Library.
    \remark This software was programmed by Michael Lubasch when he was a PhD student in the Theory Division of
            the Max Planck Institute of Quantum Optics where his work was supervised by Prof. Dr. J. Ignacio Cirac
            and Dr. Mari-Carmen Ba\~{n}uls.
    \copyright This software is distributed under the PEPS Library License v1.0
               (see accompanying file LICENSE.txt). */

using namespace std;

namespace MathAuxiliary
{
 float complexConjugate(float element) { return element; }

 double complexConjugate(double element) { return element; }

 complex<float> complexConjugate(complex<float> element) { return conj(element); }

 complex<double> complexConjugate(complex<double> element) { return conj(element); }

 void convertComplex(const complex<float>& c, float& r0) { r0 = c.real(); }

 void convertComplex(const complex<float>& c, double& r0) { r0 = c.real(); }

 void convertComplex(const complex<float>& c, complex<float>& c0) { c0 = c; }

 void convertComplex(const complex<float>& c, complex<double>& c0) { c0 = c; }

 void convertComplex(const complex<double>& c, float& r0) { r0 = c.real(); }

 void convertComplex(const complex<double>& c, double& r0) { r0 = c.real(); }

 void convertComplex(const complex<double>& c, complex<float>& c0) { c0 = c; }

 void convertComplex(const complex<double>& c, complex<double>& c0) { c0 = c; }

 double convertToDouble(float f) { return f; }

 double convertToDouble(double d) { return d; }

 double convertToDouble(complex<float> cf) { return cf.real(); }

 double convertToDouble(complex<double> cd) { return cd.real(); }
};
