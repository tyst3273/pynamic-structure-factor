# pynamic-structure-factor 
python code to calculate inelastic-neutron-scattering dynamic structure factor, S(**Q**,*w*), from molecular dynamics trajectories using MPI parallelism over **Q**-points.  

## note:
- this code assumes the classical approx. so that everything commutes and we can directly evaluate S(**Q**,*w*) by time and space FTs. this violates detailed balance (i.e. there is no distinction between energy gain/loss of the neutron and some other symmetry requirements aren't fulfilled). nevertheless, this approximation makes it possible to do this calculation. if you can write a code that calculates it from purely quantum trajectories, please share that with me :)
- also, none of the FTs are properly scaled. i.e. it should be true that for some square-integrable function f, Integral[f^2] = Integral[FT(f^2)]. this is enforced by properly normalizing the FTs in a code, but since the scattered intensity is proportional to S(**Q**,*w*)/flux, and we don't have any notion of flux in our calculations, i didnt care about scaling the FT appropriately. 

## contents
- modules
  - all the classes and objects used in the calculation
- example_inputs
  - example input_params file read by the main code
  - main python script to run the calculation
- utils
  - some python scripts for plotting
  - making input Q lists
  - compressing txt MD outputs into hdf5 databases to be read by the code. 

## intstructions
- to be done in detail eventually... until then, if you really want to use, please email me:
  - ty.sterling@colorado.edu
- also can use input files as examples and try to figure it out. i will write better docs when i am done with the project i am using this for. 

## new in v1.0
- parallelized using mpi4py
  - can now run in parallel using mpi4py. install this with conda and use the mpiexec executable that comes with it (i recommend a new environment to not interfere with your local MPI installation). 
  - only parallelized over Q-points, but could also split over 'blocks' of trajectory data.
  - if you dont want to use mpi4py, you can run the 'serial' version that is in the example_inputs directory. it uses the same methods without parallelization. this way, can avoid installing mpi4py.
  - can also run the mpi version with -np 1 or without calling the mpi executable. the results are identical. but why? 
  - it is mandatory when running with mpi to pass the '-m mpi4py' argument to the python interpreter when calling. if not, exceptions or errors will only kill the calling process and the code will enter deadlock and run forever. if you pass the module argument, the code will know to kill all processes if an exception is raised on any process. see the mpi4py docs for more info, cause i dont understand anymore than what i wrote here.
  - example syntax: mpiexec -np 4 python -m mpi4py PSF.py > log
  - the number of Q points will be split as evenly as possible over all procs with the remainder going on rank 0.

## TODO
- use pyFFTW with threading instead of scipy or numpy ffts, as these are way faster. 
  - can i convert my space-FT's into FFT? i don't think so, since we are summing over Q in other BZ's.

