# What the code does (big picture)

It simulates 2-D heat diffusion on a rectangular grid using the classic **5-point stencil**. You can optionally “inject” heat from a set of point sources every few steps, choose periodic or non-periodic boundaries, and (optionally) print/save status along the way.

Internally it keeps **two arrays** (planes) with a **1-cell halo** in each direction (so allocated size is `(Nx+2)×(Ny+2)`), and it **ping-pongs** between them:

* `planes[OLD]` → read (current state)
* `planes[NEW]` → write (next state)
  then swap.

# Main moving parts

## Command-line interface (CLI)

Handled in `initialize(...)` via `getopt`. Options (with defaults):

* `-x <int>`: grid size in x (default `1000`)
* `-y <int>`: grid size in y (default `1000`)
* `-e <int>`: number of energy sources (default `1`)
* `-E <double>`: energy injected **per source** when injection happens (default `1.0`)
* `-f <double>`: **frequency** of energy injection as a *fraction* of total iterations (default `0.0`, i.e., never; if you set 1.0 it injects every step; in general the code turns this into an integer “every k steps”)
* `-n <int>`: number of iterations (default `100`)
* `-p <0|1>`: periodic boundaries (default `0` = off)
* `-o <0|1>`: print the energy budget at every step (default `0` = off)

You’ll also be asked (or the code will produce) a list of **source positions** (array `Sources` of length `2*Nsources`, storing pairs `(x, y)` in *interior* indices).

## Core kernels (declared in the header)

* `update_plane(periodic, size, old, new)`

  * The 5-point stencil: for each interior cell $(i,j)$,

    $$
    \texttt{new}[i,j] = \texttt{old}[i,j] + \alpha \big( \texttt{old}[i-1,j]+\texttt{old}[i+1,j]+\texttt{old}[i,j-1]+\texttt{old}[i,j+1]-4\,\texttt{old}[i,j] \big)
    $$

    where `alpha` is a constant in the header (there’s a `#define alpha_guess 0.5` comment there to “mimic heat diffusivity”).
    The function knows the array has **halo cells** (hence `+2` stride in macros) and it honors `periodic` for neighbors at the border (wrap-around vs. fixed edge behavior).
* `inject_energy(periodic, Nsources, Sources, energy_per_source, S, plane)` (inline in the header)

  * For each source `(x,y)` it does `plane[x,y] += energy_per_source`. (With periodic boundaries, source coordinates that would fall on halos are wrapped to interior.)
* `get_total_energy(S, plane, &E)`

  * Sums the interior `S[0]×S[1]` cells (skipping the halos) and returns the total as a `double`. (There’s a comment about possible loop unrolling.)

## Program flow (`main`)

1. **Parse options & allocate** the two planes with halos, plus the sources (via `initialize`).
2. If `injection_frequency > 1`, **inject initial energy** once before the loop.
3. For `iter = 0 .. Niterations-1`:

   * If `iter % injection_frequency == 0`, **inject energy** from all sources.
   * **Stencil update**: `update_plane(periodic, S, planes[OLD], planes[NEW])`.
   * Optionally compute and print energies:

     * `get_total_energy(S, current_plane, &system_energy)`
     * keep track of `injected_heat`
     * if `-o 1`, print: `step k :: injected energy is ... updated system energy is ...`
   * **Swap** `OLD`/`NEW`.
4. After the loop, print the **final energy budget** (`injected` vs `system`).
5. Free memory (`memory_release`).

## Output / dump

There’s a `dump(const double *data, const uint size[2], const char *filename, double *min, double *max)` helper that writes each row as **float32** to a binary file named like `plane_%05d.bin`. (It also tracks min/max while writing.) If you enable dumping (you can easily add a call inside the loop), you can post-process the raw grid.

# How to build

Plain GCC (needs `-lm` for `math.h`):

```bash
gcc -O3 stencil_template_serial.c -o stencil_serial -lm
```

(If you want more aggressive vectorization, you can try `-O3 -march=native -ffast-math` on your own machine.)

# How to run (examples)

* 256×256 grid, 1 source at default place, 100 steps, inject every step, periodic boundaries, print energy every step:

  ```bash
  ./stencil_serial -x 256 -y 256 -n 100 -f 1.0 -p 1 -o 1
  ```

* 1024×1024, 4 sources, inject twice over the whole simulation (i.e., at steps 0 and 50 if `-n 100` → `-f 0.5`):

  ```bash
  ./stencil_serial -x 1024 -y 1024 -e 4 -n 100 -f 0.5 -E 0.25 -p 0 -o 1
  ```

> Notes
> • `-f` is a *fraction* of the iteration count; the code converts it to an integer frequency.
> • `-E` is the energy added **per source** at each injection time.
> • Source coordinates are stored in the `Sources` array; you can set them in `initialize_sources(...)` or by modifying how that function picks positions.

# Reading the binary dumps (optional)

The `dump` writes each row as `float32` (width = `S[0]`, height = `S[1]`). In Python:

```python
import numpy as np
nx, ny = 256, 256
A = np.fromfile('plane_00050.bin', dtype=np.float32).reshape(ny, nx)
# Now visualize
import matplotlib.pyplot as plt
plt.imshow(A); plt.colorbar(); plt.show()
```

# Key implementation details (so you can tweak safely)

* **Indexing & halos**
  Arrays are allocated as `(S[_x_] + 2) * (S[_y_] + 2)`. The header defines an `IDX(i,j)` macro like:

  ```c
  #define IDX(i,j) ( (j) * (Sx+2) + (i) )
  ```

  and loops typically run `i = 1..Sx` and `j = 1..Sy` for the interior. That’s why you’ll see `+2` in strides.

* **Periodic boundaries**
  When `periodic=1`, neighbor lookups on the edges wrap to the opposite side. Otherwise, edges effectively use the existing halo (which you can keep equal to edge values or zeros depending on how you initialize).

* **Energy accounting**
  The code tracks `injected_heat` (sum of all `E` added through sources) and compares it to the **system energy** (sum of grid values) to help you sanity-check conservation/dissipation.

* **Performance knobs**
  The inner loops are “plain C” and marked with comments about unrolling. If you want to experiment:

  ```c
  #pragma GCC unroll 4
  for (int j = 1; j <= Sy; ++j)
    for (int i = 1; i <= Sx; ++i)
      /* stencil */
  ```

  You can also try tiling (block the `i,j` loops) for cache locality on large grids.

# Typical modifications you might want

* **Change the diffusivity $\alpha$**
  There’s a define in the header (e.g., `alpha_guess 0.5`). Lower it to be more stable if you see blow-ups with large `-E` or very frequent injections; raise it to diffuse faster.

* **Different boundary conditions**
  If you want fixed-temperature edges (Dirichlet), set halo cells to the boundary value each iteration before the stencil. For insulating (Neumann), mirror interior values into halos.

* **Custom source locations**
  Edit `initialize_sources(...)` to set `Sources[2*s] = xi; Sources[2*s+1] = yj;` for your preferred coordinates (1..Sx, 1..Sy).

---

If you want, paste the exact `initialize_sources(...)` body you’re using and I’ll annotate how it positions the sources and where to tweak it. Or tell me a run you care about (sizes, frequency, periodic or not), and I’ll give you a one-liner command plus a quick plotting snippet.
