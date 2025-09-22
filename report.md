# Theoretical Background

Fourier's heat equation:

let $u(x,y,t)$ be the internal energy at point $(x,y)$ and time $t$, then

$$
\frac{\partial u}{\partial t} = \alpha \nabla^2 u
$$

where $\alpha$ is the thermal diffusivity constant s.t. $\alpha > 0$, and $\nabla^2$ is the Laplacian operator over the spatial dimensions $x,y,z$.

In 2D, the heat equation is:

$$
\frac{\partial u}{\partial t} = \alpha \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right)
$$

We can discretize the spatial domain into a grid of points, and approximate the second derivatives using finite differences. For a grid point $(i,j)$, the discrete Laplacian can be approximated as:

$$
\nabla^2 u_{i,j} \approx u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - 4u_{i,j}
$$	

Using an explicit time-stepping scheme, we can update the temperature at each grid point as follows:

$$
u_{i,j}^{n+1} = u_{i,j}^n + \Delta t \cdot \alpha \left( u_{i+1,j}^n + u_{i-1,j}^n + u_{i,j+1}^n + u_{i,j-1}^n - 4u_{i,j}^n \right)
$$

where $\Delta t$ is the time step size, and $n$ denotes the current time step.




# Correctness of the serial code

## in `stencil_template_serial.c`

* in `memory_allocate` function, handled the case of invalid pointer:

```c
if (planes_ptr == NULL ){
        fprintf(stderr, "Error: invalid pointer in memory_allocate\n");
		return 1; // TODO: is this enough?
	}
```

* added `typedef unsigned int uint;` at the top of the file, after includes, to define the missing `uint` type

* param initialization handling for bad input.

```c
	if ( S[_x_] < 1 || S[_y_] < 1 ){
		fprintf(stderr, "Error: invalid grid size %d x %d\n", S[_x_], S[_y_] );
		return 1;
	}
```


## in header file
* errore che introduceva anisotropia

```c

for ( int i = 1; i <= xsize; i++ )
                {
                    new[ i ] = new[ IDX(i, ysize) ];
                    new[ IDX(i, ysize+1) ] = new[ i ];
                }

```


* funzione `dump` era buggata

# Parallelization strategy