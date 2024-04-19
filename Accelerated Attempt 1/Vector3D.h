#include <ostream>
#include <cmath>



/**
 * Defines 3D vectors.
 */
  class Vector3D {
  public:

    // components
    double x, y, z;

    /**
     * Constructor.
     * Initializes tp vector (0,0,0).
     */
    __host__ __device__ Vector3D() : x( 0.0f ), y( 0.0f ), z( 0.0f ) { }

    /**
     * Constructor.
     * Initializes to vector (x,y,z).
     */
    __host__ __device__ Vector3D( double x, double y, double z) : x( x ), y( y ), z( z ) { }

    /**
     * Constructor.
     * Initializes to vector (c,c,c)
     */
    __host__ __device__ Vector3D( double c ) : x( c ), y( c ), z( c ) { }

    /**
     * Constructor.
     * Initializes from existing vector
     */
    __host__ __device__ Vector3D( const Vector3D& v ) : x( v.x ), y( v.y ), z( v.z ) { }

    // returns reference to the specified component (0-based indexing: x, y, z)
    __host__ __device__ inline double& operator[] ( const int& index ) {
      return ( &x )[ index ];
    }

    // returns const reference to the specified component (0-based indexing: x, y, z)
    __host__ __device__ inline const double& operator[] ( const int& index ) const {
      return ( &x )[ index ];
    }

    __host__ __device__ inline bool operator==( const Vector3D& v) const {
      return v.x == x && v.y == y && v.z == z;
    }

    // negation
    __host__ __device__ inline Vector3D operator-( void ) const {
      return Vector3D( -x, -y, -z );
    }

    // addition
    __host__ __device__ inline Vector3D operator+( const Vector3D& v ) const {
      return Vector3D( x + v.x, y + v.y, z + v.z );
    }

    // subtraction
    __host__ __device__ inline Vector3D operator-( const Vector3D& v ) const {
      return Vector3D( x - v.x, y - v.y, z - v.z );
    }

    // right scalar multiplication
    __host__ __device__ inline Vector3D operator*( const double& c ) const {
      return Vector3D( x * c, y * c, z * c );
    }

    // scalar division
    __host__ __device__ inline Vector3D operator/( const double& c ) const {
      const double rc = 1.0f/c;
      return Vector3D( rc * x, rc * y, rc * z );
    }

    // addition / assignment
    __host__ __device__ inline void operator+=( const Vector3D& v ) {
      x += v.x; y += v.y; z += v.z;
    }

    // subtraction / assignment
    __host__ __device__ inline void operator-=( const Vector3D& v ) {
      x -= v.x; y -= v.y; z -= v.z;
    }

    // scalar multiplication / assignment
    __host__ __device__ inline void operator*=( const double& c ) {
      x *= c; y *= c; z *= c;
    }

    // scalar division / assignment
    __host__ __device__ inline void operator/=( const double& c ) {
      (*this) *= ( 1.f/c );
    }

    /**
     * Returns Euclidean length.
     */
    __host__ __device__ inline double norm( void ) const {
      return sqrt( x*x + y*y + z*z );
    }

    /**
     * Returns Euclidean length squared.
     */
    __host__ __device__ inline double norm2( void ) const {
      return x*x + y*y + z*z;
    }

    /**
     * Returns unit vector.
     */
    __host__ __device__ inline Vector3D unit( void ) const {
      double rNorm = 1.f / sqrt( x*x + y*y + z*z );
      return Vector3D( rNorm*x, rNorm*y, rNorm*z );
    }

    /**
     * Divides by Euclidean length.
     */
    __host__ __device__ inline void normalize( void ) {
      (*this) /= norm();
    }

  }; // class Vector3D

// left scalar multiplication
  __host__ __device__ inline Vector3D operator* ( const double& c, const Vector3D& v ) {
    return Vector3D( c * v.x, c * v.y, c * v.z );
  }

// dot product (a.k.a. inner or scalar product)
  __host__ __device__ inline double dot( const Vector3D& u, const Vector3D& v ) {
    return u.x*v.x + u.y*v.y + u.z*v.z ;
  }

// cross product
  __host__ __device__ inline Vector3D cross( const Vector3D& u, const Vector3D& v ) {
    return Vector3D( u.y*v.z - u.z*v.y,
                     u.z*v.x - u.x*v.z,
                     u.x*v.y - u.y*v.x );
  }

// prints components
  std::ostream& operator<<( std::ostream& os, const Vector3D& v );


