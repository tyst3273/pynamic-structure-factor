
import numpy as np
import h5py 
from FileIO import print_log, print_herald

class params:

    def __init__(self,parser):

        """
        i probably should be inhereting from Parser, but i havent mastered 
        this without side-effects yet and just want this to work for now. 

        for now, this is for serial execution. later, i will prob. refine this 
        to be initialized on each rank and only pass the Q points and blocks I want this rank
        to work on. this is why i dont want to use inheritance. i only want to copy the attributes
        that each rank should have.

        note on logging: i dont want each rank writing to the log file. I am going to put a new 
        method in FileIO that will accept a message and file handle as as an argument and only 
        write to the file if rank == 0.

        """

        # ============ copy attributes from parser ==============

        # I/O
        self.pos_dir = parser.pos_dir
        self.traj_file = parser.traj_file
        self.output_dir = parser.output_dir
        self.outfile_prefix = parser.outfile_prefix

        # MD params
        self.dt = parser.dt
        self.stride = parser.stride
        self.total_steps = parser.total_steps
        self.num_atoms = parser.num_atoms
        self.supercell = parser.supercell
        self.lattice_vectors = parser.lattice_vectors

        # calc options
        self.unwrap_pos = parser.unwrap_pos
        self.recalculate_cell_lengths = parser.recalculate_cell_lengths
        self.b = parser.b
        self.log_file = parser.log_file
        self.Qpoints_file = parser.Qpoints_file
        self.Qmin = parser.Qmin
        self.Qmax = parser.Qmax
        self.Qsteps = parser.Qsteps
        self.num_blocks = parser.num_blocks
        self.blocks = parser.blocks


        # ==========================================================
        # ------ initialize the rest of the stuff we need ----------
        # ==========================================================
    
        # open files
        self.traj_handle = h5py.File(self.traj_file,'r')
        self.log_handle = open(self.log_file,'w')

        print_herald(self.log_handle)

        message = 'this wont work for non-orthogonal lattice vectors (yet)\n'
        print_log(self.log_handle,message,msg_type='WARNING')

        # build arrays
        self.gen_Qpoints() # read Q points from file or build slice from input
        self._make_frequency_grid() # compute frequencies corresponding to time FFT

        # print input cell lengths
        message = (f'cell lengths from input: {self.lattice_vectors[0][0]} '
                   f'{self.lattice_vectors[1][1]} '
                   f'{self.lattice_vectors[2][2]} Angstrom\n')
        print_log(self.log_handle,message,msg_type='NOTE')

        # do once and for all if not recalculating
        if not self.recalculate_cell_lengths:
            message = 'using cell lenghts from input\n'
            print_log(self.log_handle,message,msg_type='NOTE')

            self.convert_Q_to_1_over_A()

        # other wise, print that will be done later
        else:
            message = 'using cell lenghts from hdf5 trajectory file\n'
            print_log(self.log_handle,message,msg_type='NOTE')

            # will convert Q to 1/A each time traj file is read


    # =============================================================
    # ****************    public methods    ***********************
    # =============================================================
    
    def gen_Qpoints(self):

        """ 
        generate Q points, either a 2d scan or read list from file.
        """

        if self.Qpoints_file == False:
            self._Qpoints_from_2d_scan() # gen 2d scan from inputs
        else:
            self._Qpoints_from_list()    # gen from input file

        message = f'Qpoints: {self.Qsteps}'
        print_log(self.log_handle,message,msg_type='Brillouin zone path')

        for Q in range(self.Qsteps):
            message = (f'{Q+1}\t{self.reduced_Q[Q,0]:2.3f} {self.reduced_Q[Q,1]:2.3f} '
                    f'{self.reduced_Q[Q,2]:2.3f} r.l.u.')
            print_log(self.log_handle,message)


    def convert_Q_to_1_over_A(self):

        """
        convert Q points on rlu grid to 1/A
        """

        self._compute_reciprocal_latt()

        self.Qpoints = np.zeros((self.Qsteps,3))
        for Q in range(self.Qsteps):
            self.Qpoints[Q,:] = (self.r_lattice_vectors[0,:]*self.reduced_Q[Q,0]+
                    self.r_lattice_vectors[1,:]*self.reduced_Q[Q,1]+
                    self.r_lattice_vectors[2,:]*self.reduced_Q[Q,2])


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
        self.reduced_Q = np.zeros((self.Qsteps,3))
        self.reduced_Q[:,0] = np.linspace(self.Qmin[0],self.Qmax[0],self.Qsteps)
        self.reduced_Q[:,1] = np.linspace(self.Qmin[1],self.Qmax[1],self.Qsteps)
        self.reduced_Q[:,2] = np.linspace(self.Qmin[2],self.Qmax[2],self.Qsteps)

    
    def _Qpoints_from_list(self):

        """
        generate Qpoint in 1/A. read from file Qpoints_file,  
        Give a csv file of Qpoints. 1 per line, each coord seperated by spaces. Overwrites other 
        definitions for Q slices
        """

        # ============== read Q point list from file =================
        self.reduced_Q = np.loadtxt(self.Qpoints_file)

        if len(self.reduced_Q.shape) == 1:
            self.reduced_Q = self.reduced_Q.reshape((1,self.reduced_Q.shape[0]))

        self.Qsteps = self.reduced_Q.shape[0]


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

        message = (f'real space lattice (Angstrom):\n'
                f'  {self.lattice_vectors[0,0]:2.3f} {self.lattice_vectors[0,1]:2.3f}'
                f' {self.lattice_vectors[0,2]:2.3f}\n  {self.lattice_vectors[1,0]:2.3f}'
                f' {self.lattice_vectors[1,1]:2.3f} {self.lattice_vectors[1,2]:2.3f}\n'
                f'  {self.lattice_vectors[2,0]:2.3f} {self.lattice_vectors[2,1]:2.3f}'
                f' {self.lattice_vectors[2,2]:2.3f}\n')
        print_log(self.log_handle,message,msg_type='NOTE')

        message = (f'reciprocal space lattice (1/Angstrom):\n'
                f'  {self.r_lattice_vectors[0,0]:2.3f} {self.r_lattice_vectors[0,1]:2.3f}'
                f' {self.r_lattice_vectors[0,2]:2.3f}\n  {self.r_lattice_vectors[1,0]:2.3f}'
                f' {self.r_lattice_vectors[1,1]:2.3f} {self.r_lattice_vectors[1,2]:2.3f}\n'
                f'  {self.r_lattice_vectors[2,0]:2.3f} {self.r_lattice_vectors[2,1]:2.3f}'
                f' {self.r_lattice_vectors[2,2]:2.3f}\n')
        print_log(self.log_handle,message)

    
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

        message = (f'max freq: {self.max_freq:2.3f} meV\n number of freq.: {self.num_freq}\n'
                   f' resolution: {self.df:2.3e} meV\n')
        print_log(self.log_handle,message,msg_type='frequency Grid')


        











