/*
 *  mysizex   :   local x-extendion of your patch
 *  mysizey   :   local y-extension of your patch
 */

#include "stencil_template_serial.h"

int main(int argc, char **argv){

	int  Niterations; // discrete time steps of dynamics
	int  periodic;
	int  S[2];  // the size of the plate in x and y directions

	int     Nsources; // number of heat sources
	int    *Sources;  // positions of the heat sources in (x,y) pairs
	double  energy_per_source;  // how much energy each source injects at each injection event

	double *planes[2]; // planes[0] is the "old" plane, planes[1] is the "new" plane
	
	double  injected_heat = 0; // total amount of injected heat in the system

	int injection_frequency; // how often the energy is injected, i.e. after how many iterations
	int output_energy_at_steps = 0; 
   
  /* argument checking and setting */

  // this function checks and sets all the needed parameters
  // and allocates the needed memory, defines the CLI
	initialize ( argc, argv, &S[0], &periodic, &Niterations,
	    &Nsources, &Sources, &energy_per_source, &planes[0],
	    &output_energy_at_steps, &injection_frequency );

	int current = OLD;

	// first energy injection
	if ( injection_frequency > 1 ){
		inject_energy( periodic, Nsources, Sources, energy_per_source, S, planes[current] );
	}

	for (int iter = 0; iter < Niterations; iter++){
		
		/* new energy from sources */
		if ( iter % injection_frequency == 0 ){
			inject_energy( periodic, Nsources, Sources, energy_per_source, S, planes[current] );
			injected_heat += Nsources*energy_per_source;
		}
					
		/* update grid points after injection */
		update_plane(periodic, S, planes[current], planes[!current] );

		if ( output_energy_at_steps ){
			double system_heat;
			get_total_energy( S, planes[!current], &system_heat);
					
			printf("step %d :: injected energy is %g, updated system energy is %g\n", iter, 
				injected_heat, system_heat );

			char filename[100];
			sprintf( filename, "plane_%05d.bin", iter );
			dump( planes[!current], S, filename, NULL, NULL );		
		}

		/* swap planes for the new iteration */
		current = !current;
		
	}
	
	/* get final heat in the system */
	double system_heat;
	get_total_energy( S, planes[current], &system_heat);

	printf("injected energy is %g, system energy is %g\n",
		injected_heat, system_heat );

	// deallocate memory
	memory_release( planes[OLD], Sources );
	return 0;
}


int initialize ( 
	int      argc,                // the argc from command line
	char   **argv,                // the argv from command line
	int     *S,                   // two-uint array defining the x,y dimensions of the grid
	int     *periodic,            // periodic-boundary tag
	int     *Niterations,         // how many iterations
	int     *Nsources,            // how many heat sources
	int    **Sources,
	double  *energy_per_source,   // how much heat per source
	double **planes,
	int     *output_energy_at_steps,
	int     *injection_frequency
)
{
	int ret;

	// set default values

	S[_x_]            = 1000;
	S[_y_]            = 1000;
	*periodic         = 0;
	*Nsources         = 1;
	*Niterations      = 99;
	*output_energy_at_steps = 0;
	*energy_per_source = 1.0;
	*injection_frequency = *Niterations;

	double freq = 0;
	
	// process the commadn line arguments

	while ( 1 ){

		int opt;

		while((opt = getopt(argc, argv, ":x:y:e:E:f:n:p:o:")) != -1){

			switch( opt ){
				case 'x': S[_x_] = (uint)atoi(optarg);
					break;

				case 'y': S[_y_] = (uint)atoi(optarg);
					break;

				case 'e': *Nsources = atoi(optarg);
					break;

				case 'E': *energy_per_source = atof(optarg);
					break;

				case 'n': *Niterations = atoi(optarg);
					break;

				case 'p': *periodic = (atoi(optarg) > 0);
					break;

				case 'o': *output_energy_at_steps = (atoi(optarg) > 0);
					break;

				case 'f': freq = atof(optarg);
					break;
					
				case 'h': printf( "valid options are ( values btw [] are the default values ):\n"
							"-x    x size of the plate [1000]\n"
							"-y    y size of the plate [1000]\n"
							"-e    how many energy sources on the plate [1]\n"
							"-E    how many energy sources on the plate [1.0]\n"
							"-f    the frequency of energy injection [0.0]\n"
							"-n    how many iterations [100]\n"
							"-p    whether periodic boundaries applies  [0 = false]\n"
							"-o    whether to print the energy budgest at every step [0 = false]\n"
							);
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

	if ( freq == 0 )
		*injection_frequency = 1;
	else
		{
		freq = (freq > 1.0 ? 1.0 : freq );
		*injection_frequency = freq * *Niterations;
		}


	// check for all the parms being meaningful
	if ( S[_x_] < 1 || S[_y_] < 1 ){
		fprintf(stderr, "Error: invalid grid size %d x %d\n", S[_x_], S[_y_] );
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

	if ( *injection_frequency < 1 || *injection_frequency > *Niterations ){
		fprintf(stderr, "Error: invalid injection frequency %d\n", *injection_frequency );
		return 5;
	}

	// allocate the needed memory
	ret = memory_allocate( S, planes );
	if ( ret != 0 ){
		fprintf(stderr, "Error: memory allocation failed\n");
		return 10;
	}
	
	// randomly spread the heat sources
	ret = initialize_sources( S, *Nsources, Sources );
	if ( ret != 0 ){
		fprintf(stderr, "Error: heat source initialization failed\n");
		return 11;
	}
	
	return 0;  
}


int memory_allocate ( const int size[2], double **planes_ptr )
	/*
	* allocate the memory for the planes
	* we need 2 planes: the first contains the
	* current data, the second the updated data
	*
	* in the integration loop then the roles are
	* swapped at every iteration
	*/
{
	if (planes_ptr == NULL ){
		fprintf(stderr, "Error: invalid pointer in memory_allocate\n");
		return 1; // TODO: is this enough?
	}

	// define the number of bytes to allocate for each plane
	// we add 2 to each dimension to account for the ghost cells
	// (i.e. the boundary conditions)
	unsigned int bytes = (size[_x_]+2)*(size[_y_]+2);

	planes_ptr[OLD] = (double*)malloc( 2*bytes*sizeof(double) );
	// fill with zeros
	memset( planes_ptr[OLD], 0, 2*bytes*sizeof(double) );
	// smart trick: planes_ptr[NEW] points to the memory
	// right after planes_ptr[OLD]
	// so we do not need to call malloc again
	// and we are sure that both planes are contiguous in memory
	planes_ptr[NEW] = planes_ptr[OLD] + bytes;
		
	return 0;
}


int initialize_sources( uint size[2], int Nsources, int **Sources )
	/*
	* randomly spread heat sources
	*/
{
	*Sources = (int*)malloc( Nsources * 2 *sizeof(uint) );
	for ( int s = 0; s < Nsources; s++ )
	{
		(*Sources)[s*2] = 1+ lrand48() % size[_x_];
		(*Sources)[s*2+1] = 1+ lrand48() % size[_y_];
	}

	return 0;
}

int memory_release ( double *data, int *sources )
{
	if( data != NULL )
		free( data );

	if( sources != NULL )
		free( sources );

	return 0;
}


// int dump ( const double *data, const uint size[2], const char *filename, double *min, double *max )
// {
// 	if ( (filename != NULL) && (filename[0] != '\0') )
// 	{
// 		FILE *outfile = fopen( filename, "w" );
// 		if ( outfile == NULL )
// 		return 2;
		
// 		float *array = (float*)malloc( size[0] * sizeof(float) );
		
// 		double _min_ = DBL_MAX;
// 		double _max_ = 0;

// 		for ( int j = 0; j < size[1]; j++ )
// 		{
// 			/*
// 			float y = (float)j / size[1];
// 			fwrite ( &y, sizeof(float), 1, outfile );
// 			*/	
			
// 			const double * restrict line = data + j*size[0];

// 			for ( int i = 0; i < size[0]; i++ ) 
// 			{
// 				array[i] = (float)line[i];
// 				_min_ = ( line[i] < _min_? line[i] : _min_ );
// 				_max_ = ( line[i] > _max_? line[i] : _max_ ); 
// 			}
			
// 			fwrite( array, sizeof(float), size[0], outfile );
// 		}
		
// 		free( array );
		
// 		fclose( outfile );

// 		if ( min != NULL ){
// 			*min = _min_;
// 		}
// 		if ( max != NULL ){
// 			*max = _max_;
// 		}

// 	}

// 	else return 1;
// }


int dump ( const double *data, const uint size[2], const char *filename,
           double *min, double *max )
{
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


