/*
	mysizex   :   local x-extension of your patch
	mysizey   :   local y-extension of your patch
 */

#include "stencil_template_parallel.h"

static inline double *cell(double *base, int pitch, int i, int j) {
	/* Return a pointer to cell (i,j) in a 2D array stored in row-major order. */
    return base + j*pitch + i;
}

static inline void pack_columns_send(const plane_t *pl, buffers_t *bufs) {
    double *base = pl->data;
    const int Nx = (int)pl->size[_x_];
    const int Ny = (int)pl->size[_y_];
    const int pitch = Nx + 2;

    for (int j = 1; j <= Ny; ++j) bufs[SEND][WEST][j-1] = *(base + j*pitch + 1);
    for (int j = 1; j <= Ny; ++j) bufs[SEND][EAST][j-1] = *(base + j*pitch + Nx);
}

static inline void unpack_columns_recv(plane_t *pl, const buffers_t *bufs) {
    double *base = pl->data;
    const int Nx = (int)pl->size[_x_];
    const int Ny = (int)pl->size[_y_];
    const int pitch = Nx + 2;

    for (int j = 1; j <= Ny; ++j) *(base + j*pitch + 0)     = bufs[RECV][WEST][j-1];
    for (int j = 1; j <= Ny; ++j) *(base + j*pitch + Nx+1 ) = bufs[RECV][EAST][j-1];
}

int main(int argc, char **argv){

	MPI_Comm myCOMM_WORLD;
	int  Rank, Ntasks;
	uint neighbours[4];

	int  Niterations;
	int  periodic;
	vec2_t S, N;
	
	int      Nsources;
	int      Nsources_local;
	vec2_t    *Sources_local;
	double   energy_per_source;

	plane_t   planes[2];  
	buffers_t buffers[2];
	
	int output_energy_stat_perstep;

  	// initialize MPI environment
  	{
		int level_obtained;
		
		// NOTE: change MPI_FUNNELED if appropriate
		MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &level_obtained );

		// Check if the desired thread level is supported
		if ( level_obtained < MPI_THREAD_FUNNELED ) {
			printf("MPI_thread level obtained is %d instead of %d\n",
				level_obtained, MPI_THREAD_FUNNELED );
			MPI_Finalize();
			exit(1); 
		}
		
		MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
		MPI_Comm_size(MPI_COMM_WORLD, &Ntasks);
		MPI_Comm_dup (MPI_COMM_WORLD, &myCOMM_WORLD);
  	}
  
	/* argument checking and setting */
	int ret = initialize ( &myCOMM_WORLD, Rank, Ntasks, argc, argv, 
		&S, // the size of the global plane 
		&N, // the size of the MPI tasks grid
		&periodic, &output_energy_stat_perstep,
		neighbours, &Niterations, &Nsources, &Nsources_local, &Sources_local, &energy_per_source, 
		&planes[0], // planes in the current task (rank)
		           // planes[0] is the "old" plane, planes[1] is the "new" plane
		&buffers[0] );

	if ( ret ){
		printf("task %d is opting out with termination code %d\n",
			Rank, ret );
		
		MPI_Finalize();
		return 0;
	}
	
	
	int current = OLD;
	double t1 = MPI_Wtime();   /* take wall-clock time */
	
	for (int iter = 0; iter < Niterations; ++iter){

		MPI_Request reqs[8]; // 4 sends + 4 recvs (for each process)

		// inject energy from local sources, only once per iteration
		inject_energy( periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N );

		// [A] fill the buffers, and/or make the buffers' pointers pointing to the correct position

		// plane base and geometry
		double *base = planes[current].data;
		// const int Nx = (int)N[_x_], Ny = (int)N[_y_];

		const int Nx = (int)planes[current].size[_x_];
		const int Ny = (int)planes[current].size[_y_];
		const int pitch = Nx + 2;

		// ROW buffers: point directly into plane (no copy)
		buffers[SEND][NORTH] = cell(base, pitch, 1, 1);        // first interior row
		buffers[RECV][NORTH] = cell(base, pitch, 1, 0);        // north halo

		buffers[SEND][SOUTH] = cell(base, pitch, 1, Ny);       // last interior row
		buffers[RECV][SOUTH] = cell(base, pitch, 1, Ny+1);     // south halo

		// COLUMN buffers: pack into allocated contiguous arrays
		pack_columns_send(&planes[current], &buffers);

		// [B] perform the halo communications
		//     (1) use Send / Recv
		//     (2) use Isend / Irecv
		//         --> can you overlap communication and compution in this way?

		int r = 0;
		// Post receives FIRST
		// receive from north, send to south, west to east, east to west
		MPI_Irecv(buffers[RECV][NORTH], Nx, MPI_DOUBLE, neighbours[NORTH], 100, myCOMM_WORLD, &reqs[r++]);
		MPI_Irecv(buffers[RECV][SOUTH], Nx, MPI_DOUBLE, neighbours[SOUTH], 101, myCOMM_WORLD, &reqs[r++]);
		MPI_Irecv(buffers[RECV][WEST],  Ny, MPI_DOUBLE, neighbours[WEST],  102, myCOMM_WORLD, &reqs[r++]);
		MPI_Irecv(buffers[RECV][EAST],  Ny, MPI_DOUBLE, neighbours[EAST],  103, myCOMM_WORLD, &reqs[r++]);

		// Then sends
		MPI_Isend(buffers[SEND][NORTH], Nx, MPI_DOUBLE, neighbours[NORTH], 101, myCOMM_WORLD, &reqs[r++]);
		MPI_Isend(buffers[SEND][SOUTH], Nx, MPI_DOUBLE, neighbours[SOUTH], 100, myCOMM_WORLD, &reqs[r++]);
		MPI_Isend(buffers[SEND][WEST],  Ny, MPI_DOUBLE, neighbours[WEST],  103, myCOMM_WORLD, &reqs[r++]);
		MPI_Isend(buffers[SEND][EAST],  Ny, MPI_DOUBLE, neighbours[EAST],  102, myCOMM_WORLD, &reqs[r++]);

		// Optional overlap: here we could update only interior cells (i=2..Nx-1, j=2..Ny-1)

		// Wait all comms done
		MPI_Waitall(r, reqs, MPI_STATUSES_IGNORE);
		
		// [C] copy the haloes data

		// Columns received into temporary arrays -> write into halo columns
		unpack_columns_recv(&planes[current], &buffers);

		// update grid points
		update_plane( periodic, N, &planes[current], &planes[!current] );

		// output if needed
		if ( output_energy_stat_perstep )
		output_energy_stat ( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD );

		// char fname[64];
		// sprintf(fname, "rank%03d_step%05d.bin", Rank, iter);

		// dump(planes[!current].data, planes[!current].size, fname, NULL, NULL);

		double *G = gather_global_plane(&planes[!current], Rank, Ntasks, N, S, myCOMM_WORLD);

		if (Rank == 0) {
			char fname[64];
			sprintf(fname, "global_step_%05d.bin", iter);
			dump_global(G, (int)S[_x_], (int)S[_y_], fname, NULL, NULL);
			free(G);
		}

		/* swap plane indexes for the new iteration */
		current = !current;
	}
	
	t1 = MPI_Wtime() - t1;

	output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	
	memory_release( buffers, planes );
	
	MPI_Finalize();
	return 0;
}

int initialize ( 
	MPI_Comm *Comm,
	int      Me,                  // the rank of the calling process
	int      Ntasks,              // the total number of MPI ranks
	int      argc,                // the argc from command line
	char   **argv,                // the argv from command line
	vec2_t  *S,                   // the size of the plane
	vec2_t  *N,                   // two-uint array defining the MPI tasks' grid
	int     *periodic,            // periodic-boundary tag
	int     *output_energy_stat,
	int     *neighbours,          // four-int array that gives back the neighbours of the calling task
	int     *Niterations,         // how many iterations
	int     *Nsources,            // how many heat sources
	int     *Nsources_local,
	vec2_t **Sources_local,
	double  *energy_per_source,   // how much heat per source
	plane_t *planes,
	buffers_t *buffers
)
{
	int halt = 0;
	int ret;
	int verbose = 0;
	
	// set deffault values

	(*S)[_x_]         = 10000;
	(*S)[_y_]         = 10000;
	*periodic         = 0;
	*Nsources         = 4;
	*Nsources_local   = 0;
	*Sources_local    = NULL;
	*Niterations      = 1000;
	*energy_per_source = 1.0;

	if ( planes == NULL ) {
		// manage the situation
		return -1;
	}

	planes[OLD].size[0] = planes[OLD].size[0] = 0;
	planes[NEW].size[0] = planes[NEW].size[0] = 0;
	
	for ( int i = 0; i < 4; i++ )
		neighbours[i] = MPI_PROC_NULL;

	for ( int b = 0; b < 2; b++ )
		for ( int d = 0; d < 4; d++ )
		buffers[b][d] = NULL;
	
	// CLI parsing
	while ( 1 )
	{
		int opt;
		while((opt = getopt(argc, argv, ":hx:y:e:E:n:o:p:v:")) != -1)
		{
		switch( opt )
		{
		case 'x': (*S)[_x_] = (uint)atoi(optarg);
			break;

		case 'y': (*S)[_y_] = (uint)atoi(optarg);
			break;

		case 'e': *Nsources = atoi(optarg);
			break;

		case 'E': *energy_per_source = atof(optarg);
			break;

		case 'n': *Niterations = atoi(optarg);
			break;

		case 'o': *output_energy_stat = (atoi(optarg) > 0);
			break;

		case 'p': *periodic = (atoi(optarg) > 0);
			break;

		case 'v': verbose = atoi(optarg);
			break;

		case 'h': {
			if ( Me == 0 )
			printf( "\nvalid options are ( values btw [] are the default values ):\n"
				"-x    x size of the plate [10000]\n"
				"-y    y size of the plate [10000]\n"
				"-e    how many energy sources on the plate [4]\n"
				"-E    how many energy sources on the plate [1.0]\n"
				"-n    how many iterations [1000]\n"
				"-p    whether periodic boundaries applies  [0 = false]\n\n"
				);
			halt = 1; }
			break;
			
			
		case ':': printf( "option -%c requires an argument\n", optopt);
			break;
			
		case '?': printf(" -------- help unavailable ----------\n");
			break;
		}
		}

		if ( opt == -1 )
		break;
	}

	if ( halt )
		return 1;
	
	// check for all the parms being meaningful
	if ( (*S)[_x_] < 1 || (*S)[_y_] < 1 ){
		fprintf(stderr, "Error: invalid grid size %u x %u\n", (*S)[_x_], (*S)[_y_] );
		return 1;
	}

	if ( *Nsources < 1 ){
		fprintf(stderr, "Error: invalid number of sources %d\n", *Nsources );
		return 2;
	}

	if ( *Niterations < 1 ){
		fprintf(stderr, "Error: invalid number of iterations %d\n", *Niterations );
		return 3;
	}

	if ( *energy_per_source < 0 ){
		fprintf(stderr, "Error: invalid energy per source %g\n", *energy_per_source );
		return 4;
	}

	
	/*
	* find a suitable domain decomposition
	* very simple algorithm, you may want to
	* substitute it with a better one
	*
	* the plane Sx x Sy will be solved with a grid
	* of Nx x Ny MPI tasks
	*/

	vec2_t Grid;
	double formfactor = ((*S)[_x_] >= (*S)[_y_] ? (double)(*S)[_x_]/(*S)[_y_] : (double)(*S)[_y_]/(*S)[_x_] );
	int    dimensions = 2 - (Ntasks <= ((int)formfactor+1) );

	
	if ( dimensions == 1 ){
		if ( (*S)[_x_] >= (*S)[_y_] )
		Grid[_x_] = Ntasks, Grid[_y_] = 1;
		else
		Grid[_x_] = 1, Grid[_y_] = Ntasks;
	} else {
		int   Nf;
		uint *factors;
		uint  first = 1;
		ret = simple_factorization( Ntasks, &Nf, &factors );
		
		for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ )
		first *= factors[i];

		if ( (*S)[_x_] > (*S)[_y_] )
		Grid[_x_] = Ntasks/first, Grid[_y_] = first;
		else
		Grid[_x_] = first, Grid[_y_] = Ntasks/first;
	}

	(*N)[_x_] = Grid[_x_];
	(*N)[_y_] = Grid[_y_];
	

	// my coordinates in the grid of processors
	int X = Me % Grid[_x_];
	int Y = Me / Grid[_x_];

	// find my neighbours
	if ( Grid[_x_] > 1 ){  
		if ( *periodic ) {       
		neighbours[EAST]  = Y*Grid[_x_] + (Me + 1 ) % Grid[_x_];
		neighbours[WEST]  = (X%Grid[_x_] > 0 ? Me-1 : (Y+1)*Grid[_x_]-1); }
		
		else {
		neighbours[EAST]  = ( X < Grid[_x_]-1 ? Me+1 : MPI_PROC_NULL );
		neighbours[WEST]  = ( X > 0 ? (Me-1)%Ntasks : MPI_PROC_NULL ); }  
	}

	if ( Grid[_y_] > 1 ){
		if ( *periodic ) {      
		neighbours[NORTH] = (Ntasks + Me - Grid[_x_]) % Ntasks;
		neighbours[SOUTH] = (Ntasks + Me + Grid[_x_]) % Ntasks; }

		else {    
		neighbours[NORTH] = ( Y > 0 ? Me - Grid[_x_]: MPI_PROC_NULL );
		neighbours[SOUTH] = ( Y < Grid[_y_]-1 ? Me + Grid[_x_] : MPI_PROC_NULL ); }
	}

	/*
	* every MPI task determines the size sx x sy of its own domain
	* REMIND: the computational domain will be embedded into a frame
	*         that is (sx+2) x (sy+2)
	*         the outern frame will be used for halo communication or
	*/
	
	vec2_t mysize;
	// compute the size of the patch
	uint s = (*S)[_x_] / Grid[_x_];
	uint r = (*S)[_x_] % Grid[_x_];
	mysize[_x_] = s + (X < r);
	s = (*S)[_y_] / Grid[_y_];
	r = (*S)[_y_] % Grid[_y_];
	mysize[_y_] = s + (Y < r);

	planes[OLD].size[0] = mysize[0];
	planes[OLD].size[1] = mysize[1];
	planes[NEW].size[0] = mysize[0];
	planes[NEW].size[1] = mysize[1];
	

	if ( verbose > 0 ){
		if ( Me == 0 ){
			printf("Tasks are decomposed in a grid %d x %d\n\n",
				Grid[_x_], Grid[_y_] );
			fflush(stdout);
		}

		MPI_Barrier(*Comm);
		
		for ( int t = 0; t < Ntasks; t++ ){
			if ( t == Me ){
				printf("Task %4d :: "
					"\tgrid coordinates : %3d, %3d\n"
					"\tneighbours: N %4d    E %4d    S %4d    W %4d\n",
					Me, X, Y,
					neighbours[NORTH], neighbours[EAST],
					neighbours[SOUTH], neighbours[WEST] );
				fflush(stdout);
			}

			MPI_Barrier(*Comm);
		}
	}

	
	// allocate the needed memory
	ret = memory_allocate( neighbours, mysize, buffers, planes );

	// allocate the heat sources
	ret = initialize_sources( Me, Ntasks, Comm, mysize, *Nsources, Nsources_local, Sources_local );
	
	return 0;  
}


uint simple_factorization( uint A, int *Nfactors, uint **factors )
/*
 * rought factorization;
 * assumes that A is small, of the order of <~ 10^5 max,
 * since it represents the number of tasks
 #
 */
{
	int N = 0;
	int f = 2;
	uint _A_ = A;

	while ( f < A )
		{
		while( _A_ % f == 0 ) {
		N++;
		_A_ /= f; }

		f++;
		}

	*Nfactors = N;
	uint *_factors_ = (uint*)malloc( N * sizeof(uint) );

	N   = 0;
	f   = 2;
	_A_ = A;

	while ( f < A )
		{
		while( _A_ % f == 0 ) {
		_factors_[N++] = f;
		_A_ /= f; }
		f++;
		}

	*factors = _factors_;
	return 0;
}

int initialize_sources( 
	int       Me,
	int       Ntasks,
	MPI_Comm *Comm,
	vec2_t    mysize,
	int       Nsources,
	int      *Nsources_local,
	vec2_t **Sources )
{

	srand48(time(NULL) ^ Me);
	int *tasks_with_sources = (int*)malloc( Nsources * sizeof(int) );
	
	if ( Me == 0 )
		{
		for ( int i = 0; i < Nsources; i++ )
		tasks_with_sources[i] = (int)lrand48() % Ntasks;
		}
	
	MPI_Bcast( tasks_with_sources, Nsources, MPI_INT, 0, *Comm );

	int nlocal = 0;
	for ( int i = 0; i < Nsources; i++ )
		nlocal += (tasks_with_sources[i] == Me);
	*Nsources_local = nlocal;
	
	if ( nlocal > 0 )
		{
		vec2_t * restrict helper = (vec2_t*)malloc( nlocal * sizeof(vec2_t) );      
		for ( int s = 0; s < nlocal; s++ )
		{
		helper[s][_x_] = 1 + lrand48() % mysize[_x_];
		helper[s][_y_] = 1 + lrand48() % mysize[_y_];
		}

		*Sources = helper;
		}
	
	free( tasks_with_sources );

	return 0;
}

// int memory_allocate ( const int *neighbours, const vec2_t N, buffers_t *buffers_ptr, plane_t *planes_ptr){
// 	/*
// 	here you allocate the memory buffers that you need to
// 	(i)  hold the results of your computation
// 	(ii) communicate with your neighbours

// 	The memory layout that I propose to you is as follows:

// 	(i) --- calculations
// 	you need 2 memory regions: the "OLD" one that contains the
// 	results for the step (i-1)th, and the "NEW" one that will contain
// 	the updated results from the step ith.

// 	Then, the "NEW" will be treated as "OLD" and viceversa.

// 	These two memory regions are indexed by *plate_ptr:

// 	planew_ptr[0] ==> the "OLD" region
// 	plamew_ptr[1] ==> the "NEW" region


// 	(ii) --- communications

// 	you may need two buffers (one for sending and one for receiving)
// 	for each one of your neighnours, that are at most 4:
// 	north, south, east amd west.      

// 	To them you need to communicate at most mysizex or mysizey
// 	double data.

// 	These buffers are indexed by the buffer_ptr pointer so
// 	that

// 	(*buffers_ptr)[SEND][ {NORTH,...,WEST} ] = .. some memory regions
// 	(*buffers_ptr)[RECV][ {NORTH,...,WEST} ] = .. some memory regions
	
// 	--->> Of course you can change this layout as you prefer
	
// 	*/

// 	if (planes_ptr == NULL )
// 		// an invalid pointer has been passed
// 		// manage the situation
// 		;


// 	if (buffers_ptr == NULL )
// 		// an invalid pointer has been passed
// 		// manage the situation
// 		;
		
// 	// ··················································
// 	// allocate memory for data
// 	// we allocate the space needed for the plane plus a contour frame
// 	// that will contains data form neighbouring MPI tasks
// 	unsigned int frame_size = (planes_ptr[OLD].size[_x_]+2) * (planes_ptr[OLD].size[_y_]+2);

// 	planes_ptr[OLD].data = (double*)malloc( frame_size * sizeof(double) );
// 	if ( planes_ptr[OLD].data == NULL )
// 		// manage the malloc fail
// 		;
// 	memset ( planes_ptr[OLD].data, 0, frame_size * sizeof(double) );

// 	planes_ptr[NEW].data = (double*)malloc( frame_size * sizeof(double) );
// 	if ( planes_ptr[NEW].data == NULL )
// 		// manage the malloc fail
// 		;
// 	memset ( planes_ptr[NEW].data, 0, frame_size * sizeof(double) );


// 	// ··················································
// 	// buffers for north and south communication 
// 	// are not really needed
// 	//
// 	// in fact, they are already contiguous, just the
// 	// first and last line of every rank's plane
// 	//
// 	// you may just make some pointers pointing to the
// 	// correct positions

// 	// or, if you prefer, just go on and allocate buffers
// 	// also for north and south communications

// 	// ··················································
// 	// allocate buffers
// 	// ··················································
	
// 	return 0;
// }


int memory_allocate ( const int *neighbours, const vec2_t N, buffers_t *buffers_ptr, plane_t *planes_ptr )
{
    if (!buffers_ptr || !planes_ptr) return 1;

    // Allocate planes (OLD/NEW), each (Nx+2)*(Ny+2)
    const uint Nx = planes_ptr[OLD].size[_x_];
    const uint Ny = planes_ptr[OLD].size[_y_];
    const size_t frame_elems = (size_t)(Nx + 2) * (Ny + 2);

    planes_ptr[OLD].data = (double*) calloc(frame_elems, sizeof(double));
    planes_ptr[NEW].data = (double*) calloc(frame_elems, sizeof(double));
    if (!planes_ptr[OLD].data || !planes_ptr[NEW].data) return 2;

    // Initialize buffer pointers to NULL
    for (int d=0; d<4; ++d) {
        (*buffers_ptr)[d] = NULL;            // SEND side will be set each step for rows; cols allocated below
    }

    // We have buffers[2] in main → buffers_ptr refers to buffers[SEND] only here.
    // Convention we’ll use:
    //   - We will later pass BOTH buffers arrays to this function twice,
    //     or (simpler) call this alloc twice from initialize: once for SEND, once for RECV.
    // To keep your current signature and usage, do both allocations here by assuming
    // buffers_ptr points to the FIRST of a pair [SEND, RECV] laid out contiguously.

    buffers_t *buf_pair = buffers_ptr; // points to buffers[SEND]
    buffers_t *buf_pair_recv = buffers_ptr + 1; // buffers[RECV] lives at +1 (because in main you have buffers[2])

    // Allocate column buffers (WEST/EAST) for SEND & RECV (size Ny)
    (*buf_pair)[WEST]  = (double*) malloc(Ny * sizeof(double));
    (*buf_pair)[EAST]  = (double*) malloc(Ny * sizeof(double));
    (*buf_pair_recv)[WEST] = (double*) malloc(Ny * sizeof(double));
    (*buf_pair_recv)[EAST] = (double*) malloc(Ny * sizeof(double));

    if (!(*buf_pair)[WEST] || !(*buf_pair)[EAST] || !(*buf_pair_recv)[WEST] || !(*buf_pair_recv)[EAST])
        return 3;

    // NOTE:
    // For NORTH/SOUTH:
    //   SEND pointers will point directly to plane rows each iteration.
    //   RECV pointers will point directly to plane halo rows each iteration.
    //   (No allocation needed.)

    return 0;
}



int memory_release ( buffers_t *buffers, plane_t *planes ){

	/*
		Deallocate all memory previously allocated in memory_allocate()
	*/

    if (planes) {
        free(planes[OLD].data);
        free(planes[NEW].data);
    }
    if (buffers) {
        // Free column buffers we allocated (WEST/EAST for SEND & RECV)
        if (buffers[SEND][WEST]) free(buffers[SEND][WEST]);
        if (buffers[SEND][EAST]) free(buffers[SEND][EAST]);
        if (buffers[RECV][WEST]) free(buffers[RECV][WEST]);
        if (buffers[RECV][EAST]) free(buffers[RECV][EAST]);
        // NORTH/SOUTH were never malloc’d
    }
    return 0;
}

int output_energy_stat ( int step, plane_t *plane, double budget, int Me, MPI_Comm *Comm ){

	/*
		Compute the total energy in the system and print it out.	
	*/

	double system_energy = 0;
	double tot_system_energy = 0;
	get_total_energy ( plane, &system_energy );
	
	MPI_Reduce ( &system_energy, &tot_system_energy, 1, MPI_DOUBLE, MPI_SUM, 0, *Comm );
	
	if ( Me == 0 )
		{
		if ( step >= 0 )
		printf(" [ step %4d ] ", step ); fflush(stdout);

		
		printf( "total injected energy is %g, "
			"system energy is %g "
			"( in avg %g per grid point)\n",
			budget,
			tot_system_energy,
			tot_system_energy / (plane->size[_x_]*plane->size[_y_]) );
		}
	
	return 0;
}


int dump ( const double *data, const uint size[2], const char *filename, double *min, double *max ){
    if (!filename || filename[0] == '\0') return 1;

    FILE *outfile = fopen(filename, "wb");
    if (!outfile) return 2;

    uint Sx = size[0], Sy = size[1];
    uint pitch = Sx + 2;                  // include gli halo
    float *row = (float*)malloc(Sx * sizeof(float));

    double _min_ = DBL_MAX, _max_ = -DBL_MAX;

    for (uint j = 1; j <= Sy; ++j)        // solo righe interne
    {
        const double *line = data + j*pitch + 1;  // salta l'halo sinistro
        for (uint i = 0; i < Sx; ++i) {
            float v = (float)line[i];
            row[i] = v;
            if (v < _min_) _min_ = v;
            if (v > _max_) _max_ = v;
        }
        fwrite(row, sizeof(float), Sx, outfile);
    }

    free(row);
    fclose(outfile);

    if (min) *min = _min_;
    if (max) *max = _max_;
    return 0;
}

static inline int block_size_1d(int S, int G, int X) {
    int base = S / G, rem = S % G;
    return base + (X < rem);
}
static inline int block_offset_1d(int S, int G, int X) {
    int base = S / G, rem = S % G;
    return X * base + (X < rem ? X : rem);
}

double* gather_global_plane(const plane_t *plane,
                            int Rank, int Ntasks,
                            const vec2_t N,   // process grid: N[_x_], N[_y_]
                            const vec2_t S,   // global size: S[_x_], S[_y_]
                            MPI_Comm comm)
{
    const int Nx = (int)plane->size[_x_];
    const int Ny = (int)plane->size[_y_];
    const int pitch = Nx + 2;

    // 1) Pack local interior (no halos) into contiguous buffer
    size_t local_size = (size_t)Nx * Ny;
    double *localbuf = (double*) malloc(local_size * sizeof(double));
    if (!localbuf) return NULL;

    for (int j = 0; j < Ny; ++j) {
        const double *row = plane->data + (j+1)*pitch + 1;
        for (int i = 0; i < Nx; ++i) {
            localbuf[(size_t)j * Nx + i] = row[i];
        }
    }

    // Trivial case: single rank => local buffer is the global grid
    if (Ntasks == 1) {
        return localbuf;  // caller must free
    }

    // 2) Only root needs recvcounts/displs and the big 'gathered' buffer
    double *gathered = NULL;
    int *recvcounts = NULL, *displs = NULL;

    int Gx = (int)N[_x_], Gy = (int)N[_y_];
    int Sx = (int)S[_x_], Sy = (int)S[_y_];

    if (Rank == 0) {
        recvcounts = (int*) malloc((size_t)Ntasks * sizeof(int));
        displs     = (int*) malloc((size_t)Ntasks * sizeof(int));
        if (!recvcounts || !displs) {
            free(localbuf);
            free(recvcounts); free(displs);
            return NULL;
        }
        // compute recvcounts and displacements in the packed array
        int offset = 0;
        for (int r = 0; r < Ntasks; ++r) {
            int gx = r % Gx;
            int gy = r / Gx;
            int sx = block_size_1d(Sx, Gx, gx);
            int sy = block_size_1d(Sy, Gy, gy);
            recvcounts[r] = sx * sy;   // how many doubles from rank r
            displs[r]     = offset;    // where r's chunk starts in 'gathered'
            offset       += recvcounts[r];
        }
        gathered = (double*) malloc((size_t)Sx * Sy * sizeof(double));
        if (!gathered) {
            free(localbuf);
            free(recvcounts); free(displs);
            return NULL;
        }
    }

    // 3) Gather the packed chunks to rank 0
    MPI_Gatherv(localbuf, (int)local_size, MPI_DOUBLE,
                gathered, recvcounts, displs, MPI_DOUBLE,
                0, comm);
    free(localbuf);

    // 4) On rank 0, remap chunks into the correct 2D positions
    if (Rank == 0) {
        double *global_grid = (double*) malloc((size_t)Sx * Sy * sizeof(double));
        if (!global_grid) {
            free(gathered); free(recvcounts); free(displs);
            return NULL;
        }

        for (int r = 0; r < Ntasks; ++r) {
            int gx = r % Gx;
            int gy = r / Gx;
            int sx = block_size_1d(Sx, Gx, gx);
            int sy = block_size_1d(Sy, Gy, gy);

            int px = block_offset_1d(Sx, Gx, gx);  // x offset of rank r patch
            int py = block_offset_1d(Sy, Gy, gy);  // y offset

            int base = displs[r]; // start of r's chunk inside 'gathered'

            for (int j = 0; j < sy; ++j) {
                for (int i = 0; i < sx; ++i) {
                    global_grid[(size_t)(py + j) * Sx + (px + i)] =
                        gathered[(size_t)base + j * sx + i];
                }
            }
        }

        free(gathered);
        free(recvcounts);
        free(displs);
        return global_grid; // caller frees
    }

    // Non-root: nothing to return
    return NULL;
}

int dump_global(const double *data, int Sx, int Sy,
                const char *filename, double *min, double *max)
{
    if (!filename || !filename[0]) return 1;
    FILE *f = fopen(filename, "wb");
    if (!f) return 2;

    float *row = (float*) malloc((size_t)Sx * sizeof(float));
    if (!row) { fclose(f); return 3; }

    double lo =  DBL_MAX, hi = -DBL_MAX;

    for (int j = 0; j < Sy; ++j) {
        const double *line = data + (size_t)j * Sx;
        for (int i = 0; i < Sx; ++i) {
            float v = (float)line[i];
            row[i] = v;
            if (v < lo) lo = v;
            if (v > hi) hi = v;
        }
        fwrite(row, sizeof(float), (size_t)Sx, f);
    }

    free(row);
    fclose(f);
    if (min) *min = lo;
    if (max) *max = hi;
    return 0;
}
