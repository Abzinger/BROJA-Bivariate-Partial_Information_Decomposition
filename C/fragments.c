// Desparate' man's way to avoid mixing up x,y,z indices:

typedef  struct { unsigned X, Y, Z; } N_t;

typedef  struct { unsigned x; }  X_t;
typedef  struct { unsigned y; }  Y_t;
typedef  struct { unsigned z; }  Z_t;

#define InitX ( X_t{0} )
#define InitY ( Y_t{0} )
#define InitZ ( Z_t{0} )
static inline X_t incX(X_t x) { ++x.x; return x; }
static inline Y_t incY(Y_t y) { ++y.y; return y; }
static inline Z_t incZ(Z_t z) { ++z.z; return z; }
static inline bool lessX(X_t x, N_t n) { return x.x<n.X; }
static inline bool lessY(Y_t y, N_t n) { return y.y<n.Y; }
static inline bool lessZ(Z_t z, N_t n) { return z.z<n.Z; }

static inline unsigned atXY(X_t x, Y_t y, N_t n)         { return x.x + n.X*y.y; }
static inline unsigned atXZ(X_t x, Z_t z, N_t n)         { return x.x + n.X*z.z; }
static inline unsigned atYZ(Y_t y, Z_t z, N_t n)         { return y.y + n.Y*z.z; }
static inline unsigned atXYZ(X_t x, Y_t y, Z_t z, N_t n) { return x.x + n.X*atYZ(y,z,n); }

typedef  struct { unsigned x,y; }  XY_t;
typedef  struct { unsigned x,z; }  XZ_t;

struct Marginals
{
     XY_t        * xy;
     XZ_t        * xz;
     long double * b_xy;
     long double * b_xz;
};


#define MAX(x,y)                               \
     ({                                                  \
          typeof(x) _x = x;                              \
          typeof(y) _y = y;                              \
          _x >= _y ? _x : _y;                            \
     })

#define MAX0(x)                                            \
     ({                                                             \
          typeof(x) _x   = x;                                       \
          const typeof(x) zero = 0;                                 \
          _x >= zero ? _x : zero;                                   \
     })


