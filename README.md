<h1>cuPoisson</h1>
CuPoisson is a GPU implementation of the 2D fast Poisson solver using CUDA. The method solves the discrete Poisson equation on a rectangular grid, assuming zero Dirichlet boundary conditions.

This code is the result of a master's thesis written by Folkert Bleichrodt at Utrecht University under the supervision of Henk Dijkstra and Rob Bisseling. A scientific paper has been published, discussing the Poisson solver as part of a method for solving a 2D PDE ocean model. The paper can be downloaded at: http://dx.doi.org/10.1016/j.ocemod.2011.10.001

<h3>How to acknowledge?</h3>
If you use this code for your research, we would appreciate it if you would refer to the following paper:

F. Bleichrodt, R. H. Bisseling, and H. A. Dijkstra. "Accelerating a barotropic ocean model using a GPU." Ocean Modelling, Volume 41 <http://www.sciencedirect.com/science/journal/14635003/41/supp/C>, 2012, Pages 16–21

<h3>Feedback</h3>
I would be grateful for any feedback/comments on the code, or questions regarding the documentation.
Please sent me an email at F.Bleichrodt [a| cwi.nl

<h3>Move to github/news</h3>
As most of you probably know by now, Google code is being discontinued. I have moved the code to Github, which also should make it easier to allow pull-requests. I'm planning to update the code which now includes a Matlab interface.
<h2>1. Running the program</h2>
A simple driver file main.c has been provided to show how you can use cupoisson in your own C code. Currently the code is not provided as a library, since the codebase is quiet small. 

After compilation, go to the build folder and run

$ ./cupoisson

to execute the sample file.

The solution is stored as a binary (or text) file. To make a contour plot, a matlab m-file plotSolutionSingle.m has been provided (or plotSolutionDouble.m if you are using double precision).

<h3>1.1 Compiling</h3>
Go to the build folder:

$ cd build

run compiler:

$ make

This will put the program cupoisson inside the build folder. To run the program you should execute

$ ./cupoisson

For cleaning up object files, run

$ make clean

<h3>1.2 Using double precision</h3>
Default, the Poisson solver is compiled in single precision format. On the GPU, the highest performance is gained when using single precision computations. If this is not accurate enough for your application, the program has to be compiled to support double precision. Please follow these steps exactly:

Edit precision.h and uncomment the line:
//#define DOUBLE_PRECISION
this will tell the compiler to use double in place of float.
2. If the code has been compiled using single precision previously,
clean up the compilation:
$ make clean
3. Edit the Makefile and uncomment the line:
NVCCPARMS += -arch sm_13
this will tell the nvcc compiler to utilize all features of devices of compute capability 1.3 (including double precision)
4. Compile:
$ make
Note that a performance penalty of approximately a factor 2 is payed by using double precision. The difference in performance might be different for your hardware. If the accuracy is high enough, always consider using single precision first!

To compile again for single precision accuracy, revert all steps 1-3.

<h3>1.3 Trouble shooting</h3>
If compilation fails, make sure that:

CUDA is installed
The CUDA library path is in your $LD_LIBRARY_PATH
Check the include directory in CFLAGS of the makefile
<h2>2. Source code contents</h2>
I will give a short description of the contents of the source code.

<h3>2.1 main.c</h3>
A driver/example file for using the Poisson solver from your C/C++ code.

<h3>2.2 poisson.cu</h3>
The poisson solver using CUFFT. This file also has the implemation of the realFFT. Currently only square domains are supported. For performance benefits, grids of size 2^n+1 are best. In this case the FFT transform of length 2(2^n+1-1) is computed (again a power of 2) for which the FFT has the highest performance.

Currently, only the grids interior is computed since zero Dirichlet boundary conditions are assumed. When writing out the grid as a text-file, zeros are padded to the boundaries. If another type of boundary condition is needed, this needs to be adjusted in the code.

<h3>2.3 utils.cu</h3>
Some utilities for error checking, as well as writing out the solution to a file.

<h3>2.4 precision.h</h3>
Here you can specify if you want to use double precision instead of single precision GPU code. Please refer to section 1.2 for all the details.

<h2>3. Additional remarks</h2>
The current code is optimized to minimize data transfer from CPU to GPU to shared mem (cache). This is why there is some amount of code duplication. At certain points in the program, data is available on the shared memory of the multiprocessors. This data has to be copied to the GPU global memory. If the result must be transposed, it would not be optimal to first copy the data to global memory and then use another kernel to do the transpose operation (which will probably use shared memory again). As a result, our CUDA kernels are quiet big and have several tasks instead of only one.

The main details of the algorithm are described in a research paper: http://dx.doi.org/10.1016/j.ocemod.2011.10.001

The fast Poisson solver exploits an eigen decomposition of the discrete Poisson operator. Discrete sine transforms, computed using a (real)FFT, form the main body of the code. We suggest the user to read sections 3.3 and 4.1 of the paper for implementation details.

Note that the code is oblivious to actual domain size. The spacing h between gridpoints is important for the solution. Therefore, multiply the solution by h*h (h^2) since the discrete Poisson operator follows from the Poisson equation using a central finite difference scheme.
