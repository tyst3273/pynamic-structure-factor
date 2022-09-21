
import numpy as np
from psf.m_error import crash


class c_trajectory:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,comm,timers):

        """
        template
        """

        # copy refs to stuff
        self.comm = comm
        self.config = config
        self.timers = timers

        # the trajectory file
        self.trajectory_file = self.config.trajectory_file

        # set up the 'blocks' of indices for calculating on
        self._get_block_inds()
       
     # ----------------------------------------------------------------------------------------------

    def _get_block_inds(self):

        """
        setup indices for block averaging
        """

        _md_steps = self.config.md_num_steps
        _num_blocks = self.config.num_trajectory_blocks

        # averging is done in frequency space; having blocks of data that are on
        # different frequency grids will require interpolation to average together
        # i want to avoid doing that so i crop the dataset so that all blocks
        # have same frequency
        if _md_steps % _num_blocks != 0:
            msg = '\n*** warning! ***\n'
            msg += 'md_num_steps must be divisible by num_trajectory_blocks\n'
            msg += 'discarding the extra timesteps\n'
            print(msg)

        _blocks = self.config.trajectory_blocks
        
        _rem = _md_steps % _num_blocks
        self.num_block_avg = self.config.num_block_avg
        self.num_steps = _md_steps-_rem
        self.block_steps = self.num_steps // _num_blocks

        # lo/hi inds to index the blocks from the trajectory file 
        self.block_inds = np.zeros((self.num_block_avg,2),dtype=int)
        self.block_inds[:,0] = _blocks*self.num_steps
        self.block_inds[:,1] = (_blocks+1)*self.num_steps

        msg = f'\n*** blocks ***\nnum_block_avg: {self.num_block_avg}\n'
        msg += f'block_steps: {self.block_steps}\n'
        print(msg)

    # ----------------------------------------------------------------------------------------------

