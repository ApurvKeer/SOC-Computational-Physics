#ifndef BLACKHOLERAYTRACER_Matrix3x3_H
#define BLACKHOLERAYTRACER_Matrix3x3_H

#include "Vector3D.h"

#include <iosfwd>

  class Matrix3x3 {

  public:


    // The default constructor.
    Matrix3x3(void) { }

    // Constructor for row major form data.
    // Transposes to the internal column major form.
    // REQUIRES: data should be of size 9.
    Matrix3x3(double * data)
    {
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
        {
          // Transpostion happens within the () query.
          (*this)(i,j) = data[i*3 + j];
        }

    }


    /**
     * Sets all elements to val.
     */
    void zero(double val = 0.0){
      entries[0] =
      entries[1] =
      entries[2] = Vector3D( val, val, val );
    }

    // /**
    //  * Returns the determinant of A.
    //  */
    // double det( void ) const;

    /**
     * Returns the Frobenius norm of A.
     */
    double norm( void ) const{
      return sqrt( entries[0].norm2() +
                 entries[1].norm2() +
                 entries[2].norm2());
    }

    /**
     * Returns a fresh 3x3 identity matrix.
     */
    static Matrix3x3 identity( void ){
      Matrix3x3 B;

    B(0,0) = 1.; B(0,1) = 0.; B(0,2) = 0.;
    B(1,0) = 0.; B(1,1) = 1.; B(1,2) = 0.;
    B(2,0) = 0.; B(2,1) = 0.; B(2,2) = 1.;

    return B;
    }

    // No Cross products for 3 by 3 matrix.

    /**
     * Returns the ith column.
     */
    Vector3D& column( int i ){return entries[i];}
    const Vector3D& column( int i ) const{return entries[i];}

    /**
     * Returns the transpose of A.
     */
    Matrix3x3 T( void ) const{
      const Matrix3x3& A( *this );
    Matrix3x3 B;

    for( int i = 0; i < 3; i++ )
      for( int j = 0; j < 3; j++ )
      {
        B(i,j) = A(j,i);
      }

    return B;
    }

    // /**
    //  * Returns the inverse of A.
    //  */
    // Matrix3x3 inv( void ) const;

    // accesses element (i,j) of A using 0-based indexing
    // where (i, j) is (row, column).
    double& operator()( int i, int j ){return entries[j][i];}
    const double& operator()( int i, int j ) const{return entries[j][i];}

    // accesses the ith column of A
    Vector3D& operator[]( int i ){return entries[i];}
    const Vector3D& operator[]( int i ) const{return entries[i];}

    // increments by B
    void operator+=( const Matrix3x3& B ){
       Matrix3x3& A( *this );
    double* Aij = (double*) &A;
    const double* Bij = (const double*) &B;

    // Add the 9 contigous vector packed double values.
    *Aij++ += *Bij++;//0
    *Aij++ += *Bij++;
    *Aij++ += *Bij++;
    *Aij++ += *Bij++;
    *Aij++ += *Bij++;//4
    *Aij++ += *Bij++;
    *Aij++ += *Bij++;
    *Aij++ += *Bij++;
    *Aij++ += *Bij++;//8
    *Aij++ += *Bij++;
    //9.
    }

    // returns -A
    Matrix3x3 operator-( void ) const{
      // returns -A (Negation).
    const Matrix3x3& A( *this );
    Matrix3x3 B;

    B(0,0) = -A(0,0); B(0,1) = -A(0,1); B(0,2) = -A(0,2);
    B(1,0) = -A(1,0); B(1,1) = -A(1,1); B(1,2) = -A(1,2);
    B(2,0) = -A(2,0); B(2,1) = -A(2,1); B(2,2) = -A(2,2);

    return B;
    }

    // returns A-B
    Matrix3x3 operator-( const Matrix3x3& B ) const{
      const Matrix3x3& A( *this );
    Matrix3x3 C;

    for( int i = 0; i < 3; i++ )
      for( int j = 0; j < 3; j++ )
      {
        C(i,j) = A(i,j) - B(i,j);
      }

    return C;
    }

    // returns c*A
    Matrix3x3 operator*( double c ) const{
      const Matrix3x3& A( *this );
    Matrix3x3 B;

    for( int i = 0; i < 3; i++ )
      for( int j = 0; j < 3; j++ )
      {
        B(i,j) = c*A(i,j);
      }

    return B;
    }

    // returns A*B
    Matrix3x3 operator*( const Matrix3x3& B ) const{
      const Matrix3x3& A( *this );
    Matrix3x3 C;

    for( int i = 0; i < 3; i++ )
      for( int j = 0; j < 3; j++ )
      {
        C(i,j) = 0.;

        for( int k = 0; k < 3; k++ )
        {
          C(i,j) += A(i,k)*B(k,j);
        }
      }

    return C;
    }

    // returns A*x
    Vector3D operator*( const Vector3D& x ) const{
      return x[0]*entries[0] + // Add up products for each matrix column.
           x[1]*entries[1] +
           x[2]*entries[2];
    }

    // divides each element by x
    void operator/=( double x ){
      Matrix3x3& A( *this );
    double rx = 1./x;

    for( int i = 0; i < 3; i++ )
      for( int j = 0; j < 3; j++ )
      {
        A( i, j ) *= rx;
      }
    }

  protected:

    // 3 by 3 matrices are represented by an array of 3 column vectors.
    Vector3D entries[3];

  }; // class Matrix3x3

// returns the outer product of u and v.
  Matrix3x3 outer( const Vector3D& u, const Vector3D& v ){
    Matrix3x3 B;

    // Opposite of an inner product.
    for( int i = 0; i < 3; i++ )
      for( int j = 0; j < 3; j++ )
      {
        B( i, j ) = u[i]*v[j];
      }

    return B;
  }

// returns c*A
  Matrix3x3 operator*( double c, const Matrix3x3& A ){
    Matrix3x3 cA;
    const double* Aij = (const double*) &A;
    double* cAij = (double*) &cA;

    *cAij++ = c * (*Aij++);//0
    *cAij++ = c * (*Aij++);
    *cAij++ = c * (*Aij++);
    *cAij++ = c * (*Aij++);
    *cAij++ = c * (*Aij++);//4
    *cAij++ = c * (*Aij++);
    *cAij++ = c * (*Aij++);
    *cAij++ = c * (*Aij++);
    *cAij++ = c * (*Aij++);//8
    *cAij++ = c * (*Aij++);
    //9
    return cA;
  }

// prints entries
  std::ostream& operator<<( std::ostream& os, const Matrix3x3& A ){
    for( int i = 0; i < 3; i++ )
    {
      os << "[ ";

      for( int j = 0; j < 3; j++ )
      {
        os << A(i,j) << " ";
      }

      os << "]" << std::endl;
    }

    return os;
  }




#endif //BLACKHOLERAYTRACER_Matrix3x3_H
