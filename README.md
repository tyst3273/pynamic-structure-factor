# pynamic-structure-factor

### Most recent updates:
- June 20, 2024: adding symmetry to reduce number of **Q**-points to irreducible set. *in-progress, not ready!*
- May 21, 2024: add command line flag `-n` to set the number of threads. some other improvements.

## Description
Python code to calculate inelastic neutron and x-ray scattering dynamic structure factors, S(**Q**,*w*), from molecular dynamics trajectories using parallelism over **Q**-points.

This code is pretty fast as written, but could be a lot faster if redone in `c`, `c++`, etc. I'm happy to help you rewrite it ;)

Look at the manual (doc/manual.pdf) and the examples for how to use the code.

## Contact
Please contact me at:
    ty.sterling@colorado.edu
if you have any questions, bug reports, requests for features/functionality, etc.

## Citing
If you use this code in any publications, please consider citing my github repository and the following 
publication for which this code and method was developed: 

repo:
> @misc{psf2023, 
  author = {Sterling, T. C.}, 
  title = {pynamic-structure-factor},
  year = {2023},
  publisher = {GitHub},
  journal = {GitHub repository},
  url = {https://github.com/tyst3273/pynamic-structure-factor }}

paper:
> @article{weadock2023,
  title = {The nature of dynamic local order in CH3NH3PbI3 and CH3NH3PbBr3},
  journal = {Joule},
  year = {2023},
  issn = {2542-4351},
  doi = {https://doi.org/10.1016/j.joule.2023.03.017 },
  url = {https://www.sciencedirect.com/science/article/pii/S2542435123001290 },
  author = {Nicholas J. Weadock and Tyler C. Sterling and Julian A. Vigil and 
           Aryeh Gold-Parker and Ian C. Smith and Ballal Ahammed and Matthew J. Krogstad 
           and Feng Ye and David Voneshen and Peter M. Gehring and Andrew M. Rappe and 
           Hans-Georg Steinrück and Elif Ertekin and Hemamala I. Karunadasa and 
           Dmitry Reznik and Michael F. Toney}}

## Updates (that aren't in the manual yet)
- Using parallelism over **Q**-points with multiprocessing may clash with OMP parallelism on the backend (`numpy`/`scipy` are built against `openBLAS` or `MKL` which may be threaded.) Neither my code nor the backend check if the other is parallized, so if you request all procs and the backend uses all threads, your computer will be very over utilized and performance will be poor. You can set the number of threads with the enviroment variable `OMP_NUM_THREADS` or you can set it at run-time by passing a flag `-n NUMBER_OF_THREADS` when you run `PSF.py`. e.g. `python PSF.py -i input.py -n 1`. Note, I havent noticed any gain with hybrid parallelism in my implementation, so I usually request 1 thread and let multiprocessing handle everything with `num_Qpoint_procs`. I will update the manual when I have time! 
- new arguments `use_symmetry`, `symm_basis_pos`, `symm_basis_types`, and `symm_magmoms` allow you to reduce a set of Q-points into an irreducible set to speed up calculations. This isn't necesarilly sensible for disorderd systems, but can be used to quickly create colorplots, volumetric plots, etc. and anyway, `mantid` symmetrizes the real experimental data this way. Note, this depends on the python API to `spglib` , so be sure to install it. `use_symmetry` switches on this functionality. `symm_basis_pos` are the idealized, fractional basis positions for the unitcell specified by `lattice_vectors`. it is an Nx3 list. these don't *necessarily* have to be the same as the positions of the basis in the simulation, but for consistency they should be. `symm_basis_types` is a list of N integers that designates the types of atoms in the unitcell, e.g. [1,2,3,3,3] for a cubic perovskite. `symm_basis_magmoms` is either collinear (scalar) or non-collinear (vector) magnetic moments on the atoms. A better option than giving these args explicitly may be just to give a space group number (which I think `spglib` supports), but there is a nuance about how the unitcell is oriented wrt. to standard setting `spglib` expects for a given spacegroup. the current implementation is fool-proof. -- WARNING - I am still testing and cant promise it works correctly, so use at your own risk. i need to check of the symops spglib returns are on the standardized primitive cell or the conventional cell the user gives. Anyway, `use_symmetry = False` works just like before.

## To-do:
Note, if you want any of these features or any others, contact me (contact info above). I am busy, so will only implement these when I am bored/I need them... I may be motivated to do it sooner by someone asking me nicely!
- Add *plugins* to read trajectories depending on the user specified format. Should take 'block inds' as arg and return types, positions, etc. 
- Look into using `multiprocessing.shared_memory` to place the trajectory into shared memory that can later be accessed by forked processes
- Consider a low-memory algorithm to calculate x-ray atomic form factors; the form factors f(Q) depend on |**Q**|, so I pre-calculate them outside of looping over trajectory data. The fastest way is to calculate an array with shape [num_atoms, num_Qpts], but this array is huge for big simulations and big **Q**-point grids (easily runs out of ram). Current plan: just calculate f(Q) for each atom type, then map to all atoms in the loop over **Q**-points. This will be a little slower, but will be memory efficient. 
- Convert big arrays to `C` order in memory ... I stupidly used `F` order (impliclty) before knowing it made a difference. This may be a big speed up!
- apply quantum correction to the distribution function. Check my notes on force constant from greens functions.
- apply quantum correction to the stationary phase approximation (prove it is stationary phase approx first.) see Tuckermans book and my notes. 
- currently, lattice vectors are constant during the simulation. it should be possible to read the lattice vectors from the file for each time step, but recalculating in cartesian coords will be slow...
- calculate incoherent neutron scattering. currently, what is calculated is all contributions using the coherent scattering length. i.e. summing over ij with no restriction in i=j

## Bugs
- it has been brought to my attention that the parallel part of my code only works on linux. the problem is with `multiprocessing`. i am looking for a solution to get it to run in parallel on windows and mac ... please bear with me ... or switch to linux! you can still run in serial (just set `num_Qpoint_procs = 1`) on windows and mac. UPDATE: I think `multiprocess`, which is a fork of `multiprocessing`, fixes this, but I
haven't investigated in detail.

*Thanks, Tyler*

![image](https://user-images.githubusercontent.com/35535765/220440178-00a59db5-2dae-4774-9e0d-2f3de4752dfd.png)

