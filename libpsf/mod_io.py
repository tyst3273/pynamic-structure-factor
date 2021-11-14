#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   !                                                                           !
#   ! Copyright 2021 by Tyler C. Sterling and Dmitry Reznik,                    !
#   ! University of Colorado Boulder                                            !
#   !                                                                           !
#   ! This file is part of the pynamic-structure-factor (PSF) software.         !
#   ! PSF is free software: you can redistribute it and/or modify it under      !
#   ! the terms of the GNU General Public License as published by the           !
#   ! Free software Foundation, either version 3 of the License, or             !
#   ! (at your option) any later version. PSF is distributed in the hope        !
#   ! that it will be useful, but WITHOUT ANY WARRANTY; without even the        !
#   ! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  !
#   ! See the GNU General Public License for more details.                      !
#   !                                                                           !
#   ! A copy of the GNU General Public License should be available              !
#   ! alongside this source in a file named gpl-3.0.txt. If not see             !
#   ! <http://www.gnu.org/licenses/>.                                           !
#   !                                                                           !
#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import numpy as np
import os
import h5py
from libpsf.mod_utils import print_stdout, PSF_exception

class traj_file:

    """
    have to open traj_file on each rank seperately. python complained when i tried to pass open file
    with mpi communicator. 
    """

    # -----------------------------------------------------------------------------------------

    def __init__(self,invars):

        """
        open the hdf5 file. had to open it on each process rather than open and copy. mpi4py 
        complains when trying to pass open files, atleast using the 'pickle' versions of send
        & recv. 
        """

        try:
            self.handle = h5py.File(invars.traj_file,'r')
        except:
            message = 'file \'{invars.traj_file}\' seems borked'
            raise PSF_exception(message)

    # -----------------------------------------------------------------------------------------

    def parse_trajectory(self,invars,sqw):

        """
        get the trajectories, atom_types, and box sizes from the hdf5 files. can add methods to 
        get the data from differnt file formats, but ultimately the pos, atom_types, and box_lengths
        arrays that are returend need to be consistent with my code. 
        pos has shape [number of time steps in block, number of atoms, 3-dimensions (i.e. x,y,z)]
        atom_types has shape [number of time steps in block, number of atoms]
        note that ids here really means TYPES, but i am too lazy to go change the rest of my code.
        for now, my code assumes that the ids are sorted at each time step so that they are 
        identical at each step. if they arent sorted, the ids(=TYPES) at each step can be used to 
        assign the correct scattering lenght to the atoms.
        box_bounds has shape [3], 1 for each cartesian direction.
        to use the lammps read, lammps should dump using a command like:

        dump            pos all h5md ${dt_dump} pos.h5 position species
        dump_modify     pos sort id

        """

        message = 'now reading positions'
        print_stdout(message,msg_type='NOTE')

        # the indicies to be sliced.
        inds = [sqw.block_index*sqw.block_steps,(sqw.block_index+1)*sqw.block_steps]

        if invars.parse_custom: # get from my old format made using the 'compressor' tool
            sqw.box_lengths[0] = np.mean(self.handle['box_bounds'][inds[0]:inds[1],1]-
                    self.handle['box_bounds'][inds[0]:inds[1],0])
            sqw.box_lengths[1] = np.mean(self.handle['box_bounds'][inds[0]:inds[1],3]-
                    self.handle['box_bounds'][inds[0]:inds[1],2])
            sqw.box_lengths[2] = np.mean(self.handle['box_bounds'][inds[0]:inds[1],5]-
                    self.handle['box_bounds'][inds[0]:inds[1],4])
            sqw.pos[:,:,:] = self.handle['pos_data'][inds[0]:inds[1],:,:]   # get the positins
            sqw.atom_types[0,:] = self.handle['atom_types'][:]                # get the atom TYPES
        else: # read from lammps output.
            sqw.box_lengths[0] = np.mean(self.handle['particles']['all']['box']['edges']['value']
                        [inds[0]:inds[1],0],axis=0)
            sqw.box_lengths[1] = np.mean(self.handle['particles']['all']['box']['edges']['value']
                        [inds[0]:inds[1],1],axis=0)
            sqw.box_lengths[2] = np.mean(self.handle['particles']['all']['box']['edges']['value']
                        [inds[0]:inds[1],2],axis=0)
            sqw.pos[:,:,:] = self.handle['particles']['all']['position']['value'][inds[0]:inds[1],:,:]
            sqw.atom_types[:,:] = self.handle['particles']['all']['species']['value'][inds[0]:inds[1],:]

        # optionally unimpose minimum image convention
        if invars.unwrap_pos:

            message = 'unwrapping positions'
            print_stdout(message,msg_type='NOTE')

            self._unwrap_positions(invars,sqw)

    # -----------------------------------------------------------------------------------------

    def close(self):

        """
        close the hdf5 file
        """

        self.handle.close()

    # =======================================================================================
    # ------------------------------ private methods ----------------------------------------
    # =======================================================================================

    def _unwrap_positions(self,invars,sqw):

        """
        un-apply minimum image convention so that positions dont jump discontinuously
        by a full box length. this wont be memory effecient but should be fast... ish

        the issue is that, if the positions are written out 'wrapped' (e.g. lammps
        x y z instead of xu yu zu), then atoms that cross the box boundary are wrapped
        back to the other side of the box. in that case, atoms near the box boundary will
        occasionally jump discontinously by a full-box lenght. i did some checking and this
        screws up the debye waller factor at low Q (i.e. wave-length ~ the box). it was
        only a minor effect, but this takes care of it and isn't super expensive.

        the solution is to 'un-impose' the minimum image convention. i.e. treat every
        atom as the center of the cell at t=0 and at all other times t', if the atom has
        deviated by half a box length (i.e. outside the cell since the atom is at the
        center), shift it back. See e.g. Allen: "Computer Simulation of Liquids".
        """

        lx = sqw.box_lengths[0] 
        ly = sqw.box_lengths[1]
        lz = sqw.box_lengths[2]

        # build an array to be added to the positions to shift them 
        shift = sqw.pos-np.tile(sqw.pos[0,:,:].reshape((1,invars.num_atoms,3)),
                reps=[sqw.block_steps,1,1])

        # check whether to shift atoms in the x direction     
        dr = -lx*(shift[:,:,0] > lx/2).astype(int)
        dr = dr+lx*(shift[:,:,0] <= -lx/2).astype(int)
        shift[:,:,0] = np.copy(dr)

        # check whether to shift atoms in the y direction  
        dr = -ly*(shift[:,:,1] > ly/2).astype(int)
        dr = dr+ly*(shift[:,:,1] <= -ly/2).astype(int)
        shift[:,:,1] = np.copy(dr)

        # check whether to shift atoms in the z direction   
        dr = -lz*(shift[:,:,2] > lz/2).astype(int)
        dr = dr+lz*(shift[:,:,2] <= -lz/2).astype(int)
        shift[:,:,2] = np.copy(dr)

        # apply the shift    
        sqw.pos = sqw.pos+shift
       
        del shift, dr, lx, ly, lz

# ============================================================================================
# --------------------------------------------------------------------------------------------
# ============================================================================================

def save_sqw(invars,Qpts,energy,sqw,f_name='SQW.hdf5'):

    """
    write the output to an hdf5 file. use the read_sqw method to read it 
    """

    num_e = energy.shape[0]
    num_Q = Qpts.shape[0]
    f_name = os.path.join(invars.output_dir,f_name)
    with h5py.File(f_name,'w') as db:
        db_energy = db.create_dataset('energy_meV',[num_e])
        db_Qpts = db.create_dataset('Qpts_rlu',[num_Q,3])
        db_sqw = db.create_dataset('sqw',[num_e,num_Q])
        db_energy[:] = energy[:]
        db_Qpts[:,:] = Qpts[:,:]
        db_sqw[:,:] = sqw[:,:]

# --------------------------------------------------------------------------------------------

def read_sqw(f_name):

    """
    read the sqw data written with save_sqw
    """

    with h5py.File(f_name,'r') as db:
        energy = db['energy_meV'][:]
        Qpts = db['Qpts_rlu'][:,:]
        sqw = db['sqw'][:,:]

    return energy, Qpts, sqw

# ---------------------------------------------------------------------------------------------

def save_bragg(invars,Qpts,bragg,f_name='BRAGG.hdf5'):

    """
    write the output to an hdf5 file. use the read_bragg method to read it 
    """

    num_Q = Qpts.shape[0]
    f_name = os.path.join(invars.output_dir,f_name)
    with h5py.File(f_name,'w') as db:
        db_Qpts = db.create_dataset('Qpts_rlu',[num_Q,3])
        db_bragg = db.create_dataset('bragg',[num_Q])
        db_Qpts[:,:] = Qpts[:,:]
        db_bragg[:] = bragg[:]

# --------------------------------------------------------------------------------------------

def read_bragg(f_name):

    """
    read the bragg data written with save_bragg
    """

    with h5py.File(f_name,'r') as db:
        Qpts = db['Qpts_rlu'][:,:]
        bragg = db['bragg'][:]

    return  Qpts, bragg

# ---------------------------------------------------------------------------------------------

def save_timeavg(invars,Qpts,timeavg,f_name='TIMEAVG.hdf5'):

    """
    write the output to an hdf5 file. use the read_timeavg method to read it 
    """

    num_Q = Qpts.shape[0]
    f_name = os.path.join(invars.output_dir,f_name)
    with h5py.File(f_name,'w') as db:
        db_Qpts = db.create_dataset('Qpts_rlu',[num_Q,3])
        db_timeavg = db.create_dataset('time_averaged',[num_Q])
        db_Qpts[:,:] = Qpts[:,:]
        db_timeavg[:] = timeavg[:]

# --------------------------------------------------------------------------------------------

def read_timeavg(f_name):

    """
    read the time averaged intensity written with save_timeavg
    """

    with h5py.File(f_name,'r') as db:
        Qpts = db['Qpts_rlu'][:,:]
        timeavg = db['time_averaged'][:]

    return  Qpts, timeavg

# ---------------------------------------------------------------------------------------------


