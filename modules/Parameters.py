
import numpy as np
import h5py 

class params:

    def __init__(self,traj_file='lammps/pos.hdf5',unwrap_pos=False,dt=1e-15,stride=32,
            total_steps=2**21,num_atoms=512,num_types=2,b={'1':4.1491e-5},supercell=[8,8,8],
            lattice_vectors=[[5.431,0,0],[0,5.431,0],[0,0,5.431]],file_format='hdf5',
            log_file='log',Qpoints_file=False,Qmin=[0,0,0],Qmax=[2,0,0],Qsteps=10,
            num_blocks=1,blocks=[0],debug=False,recalculate_cell_lengths=False):

        """
        initialize input options here. everything has a 'default' 
        but is probably not sensible in most cases. the options are
        breifly explained here. will do it better in a *.pdf 
        instruction file in the future.
        """

        # ============== logging info during calc  =================
        self.log_file = log_file                    # file to write status and other info
        self.log_handle = open(self.log_file,'w')   # open the file

        # ================== input trajectory ======================
        self.traj_file = traj_file                          # trajectory file
        self.traj_handle = h5py.File(self.traj_file,'r')    # open the database
        self.file_format = 'hdf5'                           # only hdf5 is implemented now

        # ================= MD simulation params ===================
        self.dt = dt                                            # MD step in seconds
        self.stride = stride                                    # trajectory printing interval in steps
        self.total_steps = total_steps                          # length of MD simulation (in the file!) in steps
        self.unwrap_pos = unwrap_pos                            # unimpose minimum-image convention
        self.recalculate_cell_lengths = recalculate_cell_lengths  # compute cell lengths if e.g. NPT calc

        self.num_atoms = num_atoms      # total number of atoms
        self.num_types = num_types      # total numer of distinct types
        self.b = b                      # INS coherent scattering lengths in Angstrom

        self.lattice_vectors = lattice_vectors  # angstroms, orthogonal only for now
        self.supercell = supercell              # [nx,ny,nz] integers specifying supercell size

        self.log_handle.write('\n** WARNING **\n')  # print some info to log file
        self.log_handle.write(' This might not work for non-orthogonal lattice vectors\n')
        self.log_handle.flush()

        # ============ Q points on which to calc S(Q,w) =============
        self.Qpoints_file = Qpoints_file    # optional Q-points file. overrides 2d scan if given
        self.Qmin = Qmin                    # [h,k,l] in rlu, 2d scan
        self.Qmax = Qmax                    # [h,k,l] in rlu, 2d scan 
        self.Qsteps = Qsteps                # number of steps in 2d Q scan

        # ================ other calculation options ================
        self.num_blocks = num_blocks    # number of 'blocks' to split data in file, average all blocks at the end
        self.blocks = blocks            # list containing indicies of which block to include 

        # =============================================================
        # need Q in 1/A (not rlu) to compute space FT. if using cell 
        # lengths from input, do it once and for all. if getting cell 
        # lengths from traj file (done for each block)... wait to c
        # convert Q to 1/A once cell lenghts read from file.
        # but we need Qsteps from the file before then, so get it now
        # =============================================================

        if self.recalculate_cell_lengths != True: 
            self.gen_Qpoints()
        else:
            if self.Qpoints_file != False: # if not using Qpoints file, Qsteps is from input
                self._get_Qsteps()
            

        self._make_frequency_grid() # compute frequencies corresponding to time FFT




    # =============================================================
    # ****************    public methods    ***********************
    # =============================================================
    
    def gen_Qpoints(self):

        """ 
        generate Q points, either a 2d scan or read list from file.
        """

        if self.Qpoints_file == False:
            self._Qpoints_from_2d_scan() # gen 2d scan from inputs, convert Q to 1/A
        else:
            self._Qpoints_from_list()    # gen from input file, convert Q to 1/A


    def clean_up(self):

        """
        best practice to close files when done with them, especially hdf5
        """

        self.traj_handle.close()
        self.log_handle.close()
    

    # =============================================================
    # ****************    private methods    **********************
    # =============================================================

    def _Qpoints_from_2d_scan(self): 

        """
        generate Q-points between Qmin and Qmax for 2d BZ scan. 
        there will be Qsteps number of them. generate list in rlu
        and then convert to 1/A based on lattice_vectors
        designed for orthogonal lattice vecs, but set up to be 
        easily exteneded to non-orthogonal case. 
        """

        # ========= generate BZ slice from input options ===========
        self.Qpoints = np.zeros((self.Qsteps,3))
        reduced_Q = np.zeros((self.Qsteps,3))
        reduced_Q[:,0] = np.linspace(self.Qmin[0],self.Qmax[0],self.Qsteps)
        reduced_Q[:,1] = np.linspace(self.Qmin[1],self.Qmax[1],self.Qsteps)
        reduced_Q[:,2] = np.linspace(self.Qmin[2],self.Qmax[2],self.Qsteps)
        self.reduced_Q = np.copy(reduced_Q)

        self._compute_reciprocal_latt()
        self._convert_Q_to_1_over_A()

    
    def _Qpoints_from_list(self):

        """
        generate Qpoint in 1/A. read from file Qpoints_file,  
        Give a csv file of Qpoints. 1 per line, each coord seperated by spaces. Overwrites other 
        definitions for Q slices
        """

        # ============== read Q point list from file =================
        try:
            self.reduced_Q = np.loadtxt(self.Qpoints_file)
        except: 
            print('\n\tERROR: Qpoint list seems wrong\n')
            exit()
        if len(self.reduced_Q.shape) == 1:
            self.reduced_Q = self.reduced_Q.reshape((1,self.reduced_Q.shape[0]))
        if self.reduced_Q.shape[1] != 3:
            print('\n\tERROR: Qpoint list seems wrong\n')
            exit()

        self.Qsteps = self.reduced_Q.shape[0]
        self.Qpoints = np.zeros((self.Qsteps,3)) # in this case, just total number of Q points

        self._compute_reciprocal_latt() # compute reciprocal lattice vectors
        self._convert_Q_to_1_over_A()   # convert from rlu to 1/A


    def _compute_reciprocal_latt(self):

        """
        compute reciprocal lattice vectors from real lattice
        """

        self.lattice_vectors = np.array(self.lattice_vectors)
        self.r_lattice_vectors = np.zeros((3,3))
        self.cell_vol = self.lattice_vectors[0,:].dot(np.cross(self.lattice_vectors[1,:],
            self.lattice_vectors[2,:]))
        self.r_lattice_vectors[0,:] = 2*np.pi*np.cross(self.lattice_vectors[1,:],
                self.lattice_vectors[2,:])/self.cell_vol
        self.r_lattice_vectors[1,:] = 2*np.pi*np.cross(self.lattice_vectors[2,:],
                self.lattice_vectors[0,:])/self.cell_vol
        self.r_lattice_vectors[2,:] = 2*np.pi*np.cross(self.lattice_vectors[0,:],
                self.lattice_vectors[1,:])/self.cell_vol

    def _convert_Q_to_1_over_A(self):

        """
        convert Q points on rlu grid to 1/A
        """

        self.log_handle.write('\n** Brillouin zone path **\n' 
                f' Qpoints: {self.Qsteps}:\n')
        self.log_handle.flush()

        for Q in range(self.Qsteps):
            self.Qpoints[Q,:] = (self.r_lattice_vectors[0,:]*self.reduced_Q[Q,0]+    
                    self.r_lattice_vectors[1,:]*self.reduced_Q[Q,1]+
                    self.r_lattice_vectors[2,:]*self.reduced_Q[Q,2])
            self.log_handle.write(f' {Q+1}\t{self.Qpoints[Q,0]:2.3f} '          
                    f'{self.Qpoints[Q,1]:2.3f} {self.Qpoints[Q,2]:2.3f} 1/Angst\t'
                    f'({self.reduced_Q[Q,0]:2.3f} {self.reduced_Q[Q,1]:2.3f} '
                    f'{self.reduced_Q[Q,2]:2.3f} r.l.u.) \n')
            self.log_handle.flush()


    def _get_Qsteps(self):

        """
        open the file and get the number of Qsteps in it. we need it to set up arrays
        """

        self.Qsteps = np.loadtxt(self.Qpoints_file).shape[0]


    def _make_frequency_grid(self):

        """
        generate frequencies corresponding to time FT. not used in 
        calc, but needed for specifying energy axis in plots, etc. 
        convert to 1/s in THz and convert that to meV
        """

        self.block_steps = (self.total_steps//self.stride)//self.num_blocks # dt in array
        self.max_freq = 1e-12/self.dt/self.stride/2*4.13567             # FFT of real fn
        self.meV = np.linspace(0,self.max_freq*2,self.block_steps*2)    # convert to meV
        self.num_freq = self.meV.shape[0]                               # same as block_steps
        self.df = self.meV[1]-self.meV[0]                               # frequency resolution

        self.log_handle.write('\n** Frequency Grid **\n')               # print to log file
        self.log_handle.write(f' Max freq: {self.max_freq:2.3f} meV\n')
        self.log_handle.write(f' Number of freq.: {self.num_freq}\n')
        self.log_handle.write(f' Resolution: {self.df:2.3e} meV\n')
        self.log_handle.flush()


        











