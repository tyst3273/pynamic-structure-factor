# pynamic-structure-factor 
python code to calculate inelastic-neutron-scattering dynamic structure factor, S(**Q**,*w*).  

## note:
- this code assumes the classical approx. so that everything commutes and we can directly evaluate S(**Q**,*w*) by time and space FTs. this violates detailed balance (i.e. there is no distinction between energy gain/loss of the neutron and some other symmetry requirements aren't fulfilled). nevertheless, this approximation makes it possible to do this calculation. if you can write a code that calculates it from purely quantum trajectories, please share that with me :)
- also, none of the FTs are properly scaled. i.e. it should be true that for some square-integrable function f, Integral[f^2] = Integral[FT(f^2)]. this is enforced by properly normalizing the FTs in a code, but since the scattered intensity is proportional to S(**Q**,*w*)/flux, and we don't have any notion of flux in our calculations, i didnt care about scaling the FT appropriately. 

## contents
- modules
  - all the classes and objects used in the calculation
- example_inputs
  - some example input scripts for scans along an automatically generated BZ path
  - or on an abritrary list of Q points read from a file
- utils
  some python scripts for plotting, making input Q lists, or compresses txt file MD outputs into hdf5 databases to be read by the code. 

## intstructions
- to be done in detail eventually... until then, if you really want to use, please email me:
  - ty.sterling@colorado.edu
- also can use input files as examples and try to figure it out. i will write better docs when i am done with the project i am using this for. 

