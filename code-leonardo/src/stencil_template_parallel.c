#include "stencil_template_parallel.h"
#include <sys/stat.h>   // mkdir
#include <stdio.h>

static const char* get_env_or_na(const char* key) {
    const char* v = getenv(key);
    return (v && v[0]) ? v : "n/a";
}

int main(int argc, char **argv){

	MPI_Comm myCOMM_WORLD;
	int Rank, Ntasks;
	uint neighbours[4];

	int Niterations;
	int periodic;
	vec2_t S, N;

	int Nsources;
	int Nsources_local;
	vec2_t *Sources_local;
	double energy_per_source;

	plane_t planes[2];
	buffers_t buffers[2];
	
	int output_energy_stat_perstep;

	// declaring times
	double init_t = 0.0;
	init_t -= MPI_Wtime();

  	// initialize MPI environment
  	{
		int level_obtained;

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
	int ret = initialize ( &myCOMM_WORLD, 
		Rank, // each rank will have its own subplane
		Ntasks, argc, argv, 
		&S, // the size of the global plane 
		&N, // the size of the MPI tasks grid
		&periodic, &output_energy_stat_perstep,
		neighbours, &Niterations, &Nsources, &Nsources_local, &Sources_local, &energy_per_source, 
		&planes[0], // planes in the current task (rank)
		&buffers[0] );

	if ( ret ){
		printf("task %d is opting out with termination code %d\n",
			Rank, ret );

		MPI_Finalize();
		return 0;
	}

	init_t += MPI_Wtime();

	double t_comm = 0.0, t_comp = 0.0, t_total;

	t_total = -MPI_Wtime();

	int current = OLD;
	for (int iter = 0; iter < Niterations; ++iter){

		MPI_Request reqs[8]; // 4 sends + 4 recvs (for each process)		

		double t0 = MPI_Wtime();
		// inject energy from local sources, only once per iteration
		inject_energy( periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N );
		t_comp += MPI_Wtime() - t0;

		/* -------------------------------------------------------------------------- */
		/*                             Filling the buffers                            */
		/* -------------------------------------------------------------------------- */
		
		// [A] fill the buffers, and/or make the buffers' pointers pointing to the correct position
		
		// plane base and geometry
		double *base = planes[current].data;
		const int Nx = (int)planes[current].size[_x_];
		const int Ny = (int)planes[current].size[_y_];
		const int pitch = Nx + 2;

		// start timing communication
		double t0_comm = MPI_Wtime();

		// ROW buffers: point directly into plane (no copy)
		// NORTH
		if (neighbours[NORTH] != MPI_PROC_NULL) {
			buffers[SEND][NORTH] = base + 1*pitch + 1;   // riga interna j=1
			buffers[RECV][NORTH] = base + 0*pitch + 1;   // halo j=0
		}
		// SOUTH
		if (neighbours[SOUTH] != MPI_PROC_NULL) {
			buffers[SEND][SOUTH] = base + Ny*pitch + 1;      // riga interna j=Ny
			buffers[RECV][SOUTH] = base + (Ny+1)*pitch + 1;  // halo j=Ny+1
		}
		
		if (neighbours[WEST] != MPI_PROC_NULL) {
			for (int j = 1; j <= Ny; ++j)
				buffers[SEND][WEST][j-1] = *(base + j*pitch + 1);
		}
		// EAST: prendo colonna i=Nx -> impacchetto in buffers[SEND][EAST][0..Ny-1]
		if (neighbours[EAST] != MPI_PROC_NULL) {
			for (int j = 1; j <= Ny; ++j)
				buffers[SEND][EAST][j-1] = *(base + j*pitch + Nx);
		}
	
		/* -------------------------------------------------------------------------- */
		/*                                Comunication                                */
		/* -------------------------------------------------------------------------- */
		
		// [B] perform the halo communications
		//     (1) use Send / Recv
		//     (2) use Isend / Irecv
		//         --> can you overlap communication and compution in this way?
		
		#define TAG_N_TO_S  10
		#define TAG_S_TO_N  11
		#define TAG_W_TO_E  12
		#define TAG_E_TO_W  13

		int r = 0;

		// --- NORTH: ricevo dalla mia NORTH halo ciò che il vicino NORTH invia con S_TO_N
		if (neighbours[NORTH] != MPI_PROC_NULL) {
			if (neighbours[NORTH] == Rank) {
				for (int i=0;i<Nx;i++) buffers[RECV][NORTH][i] = buffers[SEND][NORTH][i];
			} else {
				MPI_Irecv(buffers[RECV][NORTH], Nx, MPI_DOUBLE, neighbours[NORTH], TAG_S_TO_N, myCOMM_WORLD, &reqs[r++]);
				MPI_Isend(buffers[SEND][NORTH], Nx, MPI_DOUBLE, neighbours[NORTH], TAG_N_TO_S, myCOMM_WORLD, &reqs[r++]);
			}
		}

		// --- SOUTH: ricevo dal vicino SOUTH ciò che lui manda con N_TO_S
		if (neighbours[SOUTH] != MPI_PROC_NULL) {
			if (neighbours[SOUTH] == Rank) {
				for (int i=0;i<Nx;i++) buffers[RECV][SOUTH][i] = buffers[SEND][SOUTH][i];
			} else {
				MPI_Irecv(buffers[RECV][SOUTH], Nx, MPI_DOUBLE, neighbours[SOUTH], TAG_N_TO_S, myCOMM_WORLD, &reqs[r++]);
				MPI_Isend(buffers[SEND][SOUTH], Nx, MPI_DOUBLE, neighbours[SOUTH], TAG_S_TO_N, myCOMM_WORLD, &reqs[r++]);
			}
		}

		// --- WEST: ricevo dal WEST ciò che lui manda con E_TO_W; io mando al WEST con W_TO_E
		if (neighbours[WEST] != MPI_PROC_NULL) {
			if (neighbours[WEST] == Rank) {
				for (int j=0;j<Ny;j++) buffers[RECV][WEST][j] = buffers[SEND][WEST][j];
			} else {
				MPI_Irecv(buffers[RECV][WEST], Ny, MPI_DOUBLE, neighbours[WEST], TAG_E_TO_W, myCOMM_WORLD, &reqs[r++]);
				MPI_Isend(buffers[SEND][WEST], Ny, MPI_DOUBLE, neighbours[WEST], TAG_W_TO_E, myCOMM_WORLD, &reqs[r++]);
			}
		}

		// --- EAST: ricevo dall'EAST ciò che lui manda con W_TO_E; io mando all'EAST con E_TO_W
		if (neighbours[EAST] != MPI_PROC_NULL) {
			if (neighbours[EAST] == Rank) {
				for (int j=0;j<Ny;j++) buffers[RECV][EAST][j] = buffers[SEND][EAST][j];
			} else {
				MPI_Irecv(buffers[RECV][EAST], Ny, MPI_DOUBLE, neighbours[EAST], TAG_W_TO_E, myCOMM_WORLD, &reqs[r++]);
				MPI_Isend(buffers[SEND][EAST], Ny, MPI_DOUBLE, neighbours[EAST], TAG_E_TO_W, myCOMM_WORLD, &reqs[r++]);
			}
		}

		t_comm += MPI_Wtime() - t0_comm;

		// update internals 

		double t0_comp = MPI_Wtime();
		update_internal(N, &planes[current], &planes[!current]);
		t_comp += MPI_Wtime() - t0_comp;

		// Wait all comms done
		t0_comm = MPI_Wtime();
		MPI_Waitall(r, reqs, MPI_STATUSES_IGNORE);

		/* -------------------------------------------------------------------------- */
		/*                                  Unpacking                                 */
		/* -------------------------------------------------------------------------- */

		// [C] copy the haloes data

		// Columns received into temporary arrays -> write into halo columns
		if (neighbours[WEST] != MPI_PROC_NULL) {
			for (int j = 1; j <= Ny; ++j)
				*(base + j*pitch + 0) = buffers[RECV][WEST][j-1];      // halo i=0
		}
		if (neighbours[EAST] != MPI_PROC_NULL) {
			for (int j = 1; j <= Ny; ++j)
				*(base + j*pitch + (Nx+1)) = buffers[RECV][EAST][j-1];  // halo i=Nx+1
		}

		t_comm += MPI_Wtime() - t0_comm;

		// update borders 
		t0_comp = MPI_Wtime();
		update_boundary( periodic, N, &planes[current], &planes[!current] );
		t_comp += MPI_Wtime() - t0_comp;

		// output if needed
		#ifdef ENABLE_OUTPUT
		if ( output_energy_stat_perstep )
		output_energy_stat( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD );

		if (output_energy_stat){
			double *G = gather_global_plane(&planes[!current], Rank, Ntasks, N, S, myCOMM_WORLD);
			if (Rank == 0) {
				char fname[64];
				sprintf(fname, "global_step_%05d.bin", iter);
				dump_global(G, (int)S[_x_], (int)S[_y_], fname, NULL, NULL);
				free(G);
			}
		}
		#endif

		/* swap plane indexes for the new iteration */
		current = !current;
	}

	t_total += MPI_Wtime();

	output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	
	memory_release( buffers, planes );

	// reduce the times among the processes

	double t_tot_max = 0;
	MPI_Reduce ( &t_total, &t_tot_max, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD );

	double t_comp_max = 0;
	MPI_Reduce ( &t_comp, &t_comp_max, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD );

	double t_comm_max = 0;
	MPI_Reduce ( &t_comm, &t_comm_max, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD );

	// store the data from each experiment into a CSV

	if ( Rank == 0 ){
		// print global pct
		// printf("avg communication time percentage: %f %%\n", (t_comm_avg/t_tot_avg)*100.0);
		// printf("avg computation time percentage: %f %%\n", (t_comp_avg/t_tot_avg)*100.0);
		// printf("avg overhead time percentage: %f %%\n", ((t_tot_avg - t_comm_avg - t_comp_avg)/t_tot_avg)*100.0);

		/* ensure directory exists */
		mkdir("data", 0755); /* ignore error if it already exists */

		/* grab SLURM env (strings) */
		const char* slurm_job_id        = get_env_or_na("SLURM_JOB_ID");
		const char* slurm_job_name      = get_env_or_na("SLURM_JOB_NAME");
		const char* slurm_partition     = get_env_or_na("SLURM_PARTITION");
		const char* slurm_nnodes        = get_env_or_na("SLURM_NNODES");          // --nodes
		const char* slurm_cpus_per_task = get_env_or_na("SLURM_CPUS_PER_TASK");   // --cpus-per-task
		const char* slurm_tasks_per_node= get_env_or_na("SLURM_TASKS_PER_NODE");  // --ntasks-per-node (string)
		const char* slurm_ntasks        = get_env_or_na("SLURM_NTASKS");          // total MPI tasks
		const char* omp_threads = getenv("OMP_NUM_THREADS");


		FILE *fp = fopen("data/results.csv", "a");

		/* write header if new file (or if you prefer, track with a separate flag) */
		if (fp && ftell(fp) == 0) {
			fprintf(fp,
				"size_x,size_y,Ntasks,Niterations,"
				"total_time,comm_time,comp_time,"
				"slurm_job_id,slurm_job_name,partition,"
				"nodes,cpus_per_task,tasks_per_node,ntasks,omp_threads\n");
		}

		if (fp) {
			fprintf(fp,
				"%d,%d,%d,%d,"
				"%f,%f,%f,"
				"%s,%s,%s,"
				"%s,%s,%s,%s,%s\n",
				S[_x_], S[_y_], Ntasks, Niterations,
				t_tot_max, t_comm_max, t_comp_max,
				slurm_job_id, slurm_job_name, slurm_partition,
				slurm_nnodes, slurm_cpus_per_task, slurm_tasks_per_node, slurm_ntasks, omp_threads
			);
			fclose(fp);
		}
	}

	MPI_Finalize();
	return 0;
}

int initialize ( 
	MPI_Comm *Comm,
	int      Me,                  // the rank of the calling process
	int      Ntasks,              // the total number of MPI ranks
	int      argc,                // for the CLI
	char   **argv,                // 
	vec2_t  *S,                   // the size of the plane
	vec2_t  *N,                   // two-uint array defining the MPI tasks' grid
	int     *periodic,            // periodic-boundary tag
	int     *output_energy_stat,
	uint     *neighbours,          // four-int array that gives back the neighbors of the calling task
	int     *Niterations,         // numbero of iterations
	int     *Nsources,            // number of heat sources
	int     *Nsources_local,
	vec2_t **Sources_local,
	double  *energy_per_source,   // how much heat per source
	plane_t *planes,
	buffers_t *buffers
){
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
		return -1;
	}

	planes[OLD].size[_x_] = 0;
	planes[OLD].size[_y_] = 0;
	planes[NEW].size[_x_] = 0;
	planes[NEW].size[_y_] = 0;

	
	for ( int i = 0; i < 4; i++ )
		neighbours[i] = MPI_PROC_NULL;

	for ( int b = 0; b < 2; b++ )
		for ( int d = 0; d < 4; d++ )
		buffers[b][d] = NULL;
	
	/* -------------------------------------------------------------------------- */
	/*                                 CLI PARSING                                */
	/* -------------------------------------------------------------------------- */

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

	/* -------------------------------------------------------------------------- */
	/*                            DOMAIN DECOMPOSITION                            */
	/* -------------------------------------------------------------------------- */
	
	/*
	* find a suitable domain decomposition
	* very simple algorithm, you may want to
	* substitute it with a better one

	* the plane Sx x Sy will be solved with a grid
	* of Nx x Ny MPI tasks
	*/

	vec2_t Grid;
	double formfactor = ((*S)[_x_] >= (*S)[_y_] ? (double)(*S)[_x_]/(*S)[_y_] : (double)(*S)[_y_]/(*S)[_x_] );

	int    dimensions = 2 - (Ntasks <= ((int)formfactor+1) );
	// this is needed to decide if we want a 1d grid or 2d grid
	// dimensions can be either 1 or 2. 
	// 1 if Ntasks smaller than the floor of form factor
	// 2 if larger


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

		if ( ret ){
			fprintf(stderr, "Error: cannot factorize %d\n", Ntasks );
			return 5;
		}
		
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

int memory_allocate ( const uint *neighbours, const vec2_t N, buffers_t *buffers_ptr, plane_t *planes_ptr ){
    if (!buffers_ptr || !planes_ptr) return 1;

    // Allocate planes (OLD/NEW), each (Nx+2)*(Ny+2)
    const uint Nx = planes_ptr[OLD].size[_x_];
    const uint Ny = planes_ptr[OLD].size[_y_];
    const size_t frame_elems = (size_t)(Nx + 2) * (Ny + 2);

    // planes_ptr[OLD].data = (double*) calloc(frame_elems, sizeof(double));
    // planes_ptr[NEW].data = (double*) calloc(frame_elems, sizeof(double));

	planes_ptr[OLD].data = (double*) malloc(frame_elems * sizeof(double));
	planes_ptr[NEW].data = (double*) malloc(frame_elems * sizeof(double));

	if (!planes_ptr[OLD].data || !planes_ptr[NEW].data) return 2;
	
	#pragma omp parallel for schedule(static)
	for (int j = 0; j < Ny + 2; ++j){
		for (int i = 0; i < Nx + 2; ++i) {
			size_t idx = (size_t)j * (Nx + 2) + i;
			planes_ptr[OLD].data[idx] = 0.0;
			planes_ptr[NEW].data[idx] = 0.0;
		}
	}

    // Initialize buffer pointers to NULL
    for (int d=0; d<4; ++d) {
        (*buffers_ptr)[d] = NULL;            // SEND side will be set each step for rows; cols allocated below
    }

    buffers_t *buf_pair = buffers_ptr; // points to buffers[SEND]
    buffers_t *buf_pair_recv = buffers_ptr + 1; // buffers[RECV] lives at +1 (because in main we have buffers[2])

    // Allocate column buffers (WEST/EAST) for SEND & RECV (size Ny)
    (*buf_pair)[WEST]  = (double*) malloc(Ny * sizeof(double));
    (*buf_pair)[EAST]  = (double*) malloc(Ny * sizeof(double));
    (*buf_pair_recv)[WEST] = (double*) malloc(Ny * sizeof(double));
    (*buf_pair_recv)[EAST] = (double*) malloc(Ny * sizeof(double));

    if (!(*buf_pair)[WEST] || !(*buf_pair)[EAST] || !(*buf_pair_recv)[WEST] || !(*buf_pair_recv)[EAST])
        return 3;

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
        // NORTH/SOUTH were never malloc'd
    }
    return 0;
}

int output_energy_stat( int step, plane_t *plane, double budget, int Me, MPI_Comm *Comm ){

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
		printf(" [ step %4d ] ", step ); 
		fflush(stdout);
		
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
// S: total number of elements
// G: number of tasks
// X: rank of the task
// returns: number of elements that belongs to process X
    int base = S / G, rem = S % G;
    return base + (X < rem);
}
static inline int block_offset_1d(int S, int G, int X) {
// S: total number of elements
// G: number of tasks
// X: rank of the task
// returns: index of the first element that belongs to process X
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
        // recvcounts[r] number of elements root task received from task r
        displs     = (int*) malloc((size_t)Ntasks * sizeof(int));
        // displs[r] displacement in the buffer for the data received from task r
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
    MPI_Gatherv(localbuf, // send buffer
        (int)local_size, // number of elements to be gathered
        MPI_DOUBLE, // send datatype
        gathered, // recv buffer
        recvcounts, // array[r] that tells the number of elems to be received from process r 
        displs, // array[r] that tells the starting index in recv buffer for data of process r
        MPI_DOUBLE, // recv datatype
        0, // root
        comm
        );
        
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
