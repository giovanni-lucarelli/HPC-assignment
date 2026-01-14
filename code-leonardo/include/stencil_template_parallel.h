#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <math.h>
#include <float.h>

#include <omp.h>
#include <mpi.h>

#define NORTH 0
#define SOUTH 1
#define EAST  2
#define WEST  3

#define SEND 0
#define RECV 1

#define OLD 0
#define NEW 1

#define _x_ 0
#define _y_ 1

typedef unsigned int uint;

typedef uint vec2_t[2]; // array of 2 uns. int
typedef double *restrict buffers_t[4];

typedef struct {
    double   *restrict data;
    vec2_t     size;
} plane_t;


/* ------------------------------- Declarations ------------------------------- */
int initialize ( MPI_Comm *, int, int, int, char **, vec2_t *, vec2_t *, int *, int *, uint *, int *, int *, int *, vec2_t **, double *, plane_t *, buffers_t * );

int memory_release ( buffers_t *buffers, plane_t *planes );

extern int inject_energy ( const int, const int, const vec2_t  *, const double, plane_t *, const vec2_t );

extern int update_plane ( const int, const vec2_t, const plane_t *, plane_t * );

extern int get_total_energy( plane_t *, double  * );

int output_energy_stat ( int, plane_t *, double, int, MPI_Comm * );

uint simple_factorization( uint, int *, uint ** );

int initialize_sources( int, int , MPI_Comm  *, uint [2], int, int *, vec2_t ** );

int memory_allocate ( const uint *, const vec2_t, buffers_t *, plane_t * );

int dump ( const double *, const uint [2], const char *, double *, double * );

double* gather_global_plane(const plane_t *,
                            int , int,
                            const vec2_t ,   // process grid: N[_x_], N[_y_]
                            const vec2_t ,   // global size: S[_x_], S[_y_]
                            MPI_Comm);

int dump_global(const double *, int , int , const char *, double *, double *);

/* ---------------------------- Inline functions ---------------------------- */

inline int inject_energy ( 
    const int periodic,
    const int Nsources,
    const vec2_t *Sources,
    const double energy,
    plane_t *plane,
    const vec2_t Nproc   // griglia di processi: Nproc[_x_], Nproc[_y_]
)
{
    const int Nx = (int)plane->size[_x_];
    const int Ny = (int)plane->size[_y_];
    const int pitch = Nx + 2;
    double * restrict data = plane->data;

    #define IDX(i,j) ((j)*pitch + (i))

    for (int s = 0; s < Nsources; s++) {
        int x = (int)Sources[s][_x_]; // 1..Nx
        int y = (int)Sources[s][_y_]; // 1..Ny

        data[ IDX(x,y) ] += energy;

        if (periodic) {
            // copia periodica solo se tutta la direzione sta nello stesso rank
            if (Nproc[_x_] == 1) {
                if (x == 1)   data[ IDX(Nx+1, y) ] += energy; // halo destro
                if (x == Nx)  data[ IDX(0,    y) ] += energy; // halo sinistro
            }
            if (Nproc[_y_] == 1) {
                if (y == 1)   data[ IDX(x, Ny+1) ] += energy; // halo alto
                if (y == Ny)  data[ IDX(x, 0   ) ] += energy; // halo basso
            }
        }
    }

    #undef IDX
    return 0;
}

    
inline int update_plane ( 
    const int periodic,
    const vec2_t N,         // the grid of MPI tasks
    const plane_t *oldplane,
    plane_t *newplane
)
{
    uint register fxsize = oldplane->size[_x_]+2;
    // uint register fysize = oldplane->size[_y_]+2;
    uint register xsize = oldplane->size[_x_];
    uint register ysize = oldplane->size[_y_];
    
    #define IDX( i, j ) ( (j)*fxsize + (i) )

    double * restrict old = oldplane->data;
    double * restrict new = newplane->data;

    // OpenMP parallelization
    #pragma omp parallel for schedule(static)
    for (uint j = 1; j <= ysize; j++)
        for ( uint i = 1; i <= xsize; i++)
        {

            new[ IDX(i,j) ] =
                old[ IDX(i,j) ] / 2.0 + ( old[IDX(i-1, j)] + old[IDX(i+1, j)] +
                                            old[IDX(i, j-1)] + old[IDX(i, j+1)] ) /4.0 / 2.0;
            
        }

    if (periodic) {
        if (N[_y_] == 1) { // periodicità verticale tutta dentro lo stesso rank
            for (int i = 1; i <= xsize; ++i) {
                new[IDX(i, 0)]        = new[IDX(i, ysize)];
                new[IDX(i, ysize+1)]  = new[IDX(i, 1)];
            }
        }
        if (N[_x_] == 1) { // periodicità orizzontale tutta dentro lo stesso rank
            for (int j = 1; j <= ysize; ++j) {
                new[IDX(0, j)]        = new[IDX(xsize, j)];
                new[IDX(xsize+1, j)]  = new[IDX(1, j)];
            }
        }
    }
    
    #undef IDX
    return 0;
}

inline int update_internal ( 
    const vec2_t N,         // the grid of MPI tasks
    const plane_t *oldplane,
    plane_t *newplane
)
{
    uint register fxsize = oldplane->size[_x_]+2;
    uint register xsize = oldplane->size[_x_];
    uint register ysize = oldplane->size[_y_];
    
    #define IDX( i, j ) ( (j)*fxsize + (i) )

    double * restrict old = oldplane->data;
    double * restrict new = newplane->data;

    // OpenMP parallelization
    #pragma omp parallel for schedule(static)
    for (uint j = 2; j <= ysize-1; j++)
        for ( uint i = 2; i <= xsize-1; i++)
        {

            new[ IDX(i,j) ] =
                old[ IDX(i,j) ] / 2.0 + ( old[IDX(i-1, j)] + old[IDX(i+1, j)] +
                                            old[IDX(i, j-1)] + old[IDX(i, j+1)] ) /4.0 / 2.0;
            
        }
    
    #undef IDX
    return 0;
}

inline int update_boundary ( 
    const int periodic,
    const vec2_t N,         // the grid of MPI tasks
    const plane_t *oldplane,
    plane_t *newplane
)
{
    uint register fxsize = oldplane->size[_x_]+2;
    uint register xsize = oldplane->size[_x_];
    uint register ysize = oldplane->size[_y_];
    
    #define IDX( i, j ) ( (j)*fxsize + (i) )

    double * restrict old = oldplane->data;
    double * restrict new = newplane->data;

    const double alpha = 0.5;
    const double constant =  (1-alpha) / 4.0;

    double center, neighbors;

    // update top and bottom edges
    // do not update the corners here
    #pragma omp parallel for schedule(static)
    for (uint i = 2; i <= xsize-1; i++){
        // top edge
        double center = old[ IDX(i,1) ];
        double neighbors = old[IDX(i-1, 1)] + old[IDX(i+1, 1)] +
                    old[IDX(i, 0)] + old[IDX(i, 2)];
        new[ IDX(i,1) ] = center * alpha + neighbors * constant;

        // bottom edge
        center = old[ IDX(i,ysize) ];
        neighbors = old[IDX(i-1, ysize)] + old[IDX(i+1, ysize)] + old[IDX(i, ysize-1)] + old[IDX(i, ysize+1)];
        new[ IDX(i,ysize) ] = center * alpha + neighbors * constant;
    }

    // update left and right edges
    // update the corners here
    #pragma omp parallel for schedule(static)
    for (uint j = 1; j <= ysize; j++){
        // left edge
        double center = old[ IDX(1,j) ];
        double neighbors = old[IDX(0, j)] + old[IDX(2, j)] + old[IDX(1, j-1)] + old[IDX(1, j+1)];
        new[ IDX(1,j) ] = center * alpha + neighbors * constant;
        // right edge
        center = old[ IDX(xsize,j) ];
        neighbors = old[IDX(xsize-1, j)] + old[IDX(xsize+1, j)] + old[IDX(xsize, j-1)] + old[IDX(xsize, j+1)];
        new[ IDX(xsize,j) ] = center * alpha + neighbors * constant;

    }

    if (periodic) {
        if (N[_y_] == 1) { // periodicità verticale tutta dentro lo stesso rank
            for (int i = 1; i <= xsize; ++i) {
                new[IDX(i, 0)]        = new[IDX(i, ysize)];
                new[IDX(i, ysize+1)]  = new[IDX(i, 1)];
            }
        }
        if (N[_x_] == 1) { // periodicità orizzontale tutta dentro lo stesso rank
            for (int j = 1; j <= ysize; ++j) {
                new[IDX(0, j)]        = new[IDX(xsize, j)];
                new[IDX(xsize+1, j)]  = new[IDX(1, j)];
            }
        }
    }
    
    #undef IDX
    return 0;
}

inline int get_total_energy( plane_t *plane, double *energy )
{

    const int register xsize = plane->size[_x_];
    const int register ysize = plane->size[_y_];
    const int register fsize = xsize+2;

    double * restrict data = plane->data;
    
    #define IDX( i, j ) ( (j)*fsize + (i) )

    #if defined(LONG_ACCURACY)    
    long double totenergy = 0;
    #else
    double totenergy = 0;    
    #endif

    #pragma omp parallel for reduction(+:totenergy) schedule(static)
    for ( int j = 1; j <= ysize; j++ )
        for ( int i = 1; i <= xsize; i++ )
            totenergy += data[ IDX(i, j) ];

    #undef IDX

    *energy = (double)totenergy;
    return 0;
}