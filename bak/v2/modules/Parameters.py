"""
This file contains all the parameters used in the MD simulation, etc.
The scattering length b is coherect. Might ought to convert total cross section
to lenght since I am trying to include coherent and incoherent scattering.
"""
import numpy as np

class params:

    def __init__(self,traj_file='lammps/pos.dat',dt=1e-15,stride=32,total_steps=2**21,
            num_atoms=512,num_types=2,b={'1':4.1491e-5},supercell=[16,2,2],
            lattice_vectors=[[1,0,0],[0,1,0],[0,0,1]],Qmin=[1,0,1],Qmax=[2,0,1],
            Qsteps=10,blocks=10,debug=False,minimum_image=True):

        self.traj_file = traj_file # lammps trajectory file
        self.dt = dt # seconds
        self.stride = stride # printing interval in steps
        self.total_steps = total_steps

        self.num_atoms = num_atoms
        self.num_types = num_types
        self.b = b # angstrom
        self.lattice_vectors = np.array(lattice_vectors) # angstrom

        self.Qmin = Qmin # rlu
        self.Qmax = Qmax # rlu
        self.Qsteps = Qsteps
        self.blocks = blocks
        self.debug = debug
        self.minimum_image = minimum_image

        self.log_file = 'log'
        self.log_handle = open(self.log_file,'w')
        self.traj_handle = open(self.traj_file,'r')

        self.log_handle.write('\n** WARNING **\n')
        self.log_handle.write(' This might not work for non-orthogonal lattice vectors\n')
        self.log_handle.flush()

        self._make_Qparams()
        self._make_frequency_grid()

    def _make_Qparams(self):

        # compute reciprocal lattice vectors
        self.r_lattice_vectors = np.zeros((3,3))
        self.cell_vol = self.lattice_vectors[0,:].dot(np.cross(self.lattice_vectors[1,:],
            self.lattice_vectors[2,:]))
        self.r_lattice_vectors[0,:] = 2*np.pi*np.cross(self.lattice_vectors[1,:],
                self.lattice_vectors[2,:])/self.cell_vol
        self.r_lattice_vectors[1,:] = 2*np.pi*np.cross(self.lattice_vectors[2,:],
                self.lattice_vectors[0,:])/self.cell_vol
        self.r_lattice_vectors[2,:] = 2*np.pi*np.cross(self.lattice_vectors[0,:],
                self.lattice_vectors[1,:])/self.cell_vol

        # populate BZ slice
        self.Qpoints = np.zeros((self.Qsteps,3))
        reduced_Q = np.zeros((self.Qsteps,3))
        reduced_Q[:,0] = np.linspace(self.Qmin[0],self.Qmax[0],self.Qsteps)
        reduced_Q[:,1] = np.linspace(self.Qmin[1],self.Qmax[1],self.Qsteps)
        reduced_Q[:,2] = np.linspace(self.Qmin[2],self.Qmax[2],self.Qsteps)

        self.log_handle.write('\n** Brillouin zone path **\n'
                f' Qpoints: {self.Qsteps}:\n')
        self.log_handle.flush()
        for Q in range(self.Qsteps):
            self.Qpoints[Q,:] = (self.r_lattice_vectors[0,:]*reduced_Q[Q,0]+
                    self.r_lattice_vectors[1,:]*reduced_Q[Q,1]+
                    self.r_lattice_vectors[2,:]*reduced_Q[Q,2])
            self.log_handle.write(f' {Q+1}\t{self.Qpoints[Q,0]:2.3f} '
                    f'{self.Qpoints[Q,1]:2.3f} {self.Qpoints[Q,2]:2.3f}\n')
            self.log_handle.flush()

    def _make_frequency_grid(self):
        
        self.block_steps = (self.total_steps//self.stride)//self.blocks
        self.max_freq = 1e-12/self.dt/self.stride/2*4.13567 # FFT is symmetrical
        self.meV = np.linspace(0,self.max_freq*2,self.block_steps*2) # meV
        self.num_freq = self.meV.shape[0] # same as block_steps
        self.df = self.meV[1]-self.meV[0]

        self.log_handle.write('\n** Frequency Grid **\n')
        self.log_handle.write(f' Max freq: {self.max_freq:2.3f} meV\n')
        self.log_handle.write(f' Number of freq.: {self.num_freq}\n')
        self.log_handle.write(f' Resolution: {self.df:2.3e} meV\n')
        self.log_handle.flush()














