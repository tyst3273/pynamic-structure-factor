
import numpy as np
import h5py
import sys

# --------------------------------------------------------------------------------------------------

class c_convert_to_reduced_coords:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,input_file,output_file,num_reps):

        """
        convert a file containing MD trajectories in cartesian coordinates to reduced coordinates.
        
        NOTE: the input file MUST be written with my converter tool; the data MUST have the same 
            order at every time step and MUST be unwrapped (lammps ~ xu, yu, zu)
        """

        self.input_file = input_file
        self.output_file = output_file
        self.num_reps = num_reps

        print('')
        print('number of unitcells along a, b, c:',self.num_reps)
        print('')
        print('reading from file:',self.input_file)

        with h5py.File(self.input_file,'r') as idb, h5py.File(self.output_file,'w') as odb:

            self.box_vectors = idb['box_vectors'][...]
            types = idb['types'][...]

            self.cartesian_pos = idb['cartesian_pos']
            shape = self.cartesian_pos.shape

            self.num_steps = shape[0]
            self.num_atoms = shape[1]
            print('num_steps:',self.num_steps)
            print('num_atoms:',self.num_atoms)
            print('')
            
            odb.create_dataset('reduced_pos',shape=shape,dtype=float)
            odb.create_dataset('box_vectors',data=self.box_vectors)
            odb.create_dataset('types',data=types)

            self.reduced_pos = np.zeros((self.num_atoms,3),dtype=float,order='F')
            for step in range(self.num_steps):
                
                if step % 100 == 0:
                    print('now on step',step,'out of',self.num_steps)

                self._get_reduced_pos(step)

                if step == 0:
                    self.reduced_pos_0 = self.reduced_pos.copy(order='F')

                self._unwrap()

                odb['reduced_pos'][step,...] = self.reduced_pos

    # ----------------------------------------------------------------------------------------------

    def _get_reduced_pos(self,step):

        """
        convert positions of all atoms at a single timestep to reduced coordinates per unitcell.
        need supercell reps (num_reps) to convert box vectors to lattice vectors
        """

        latvecs = self.box_vectors[step,...].copy()
        latvecs[0,:] /= self.num_reps[0]
        latvecs[1,:] /= self.num_reps[1]
        latvecs[2,:] /= self.num_reps[2]

        inv_latvecs = np.linalg.inv(latvecs)
        inv_latvecs = np.tile(inv_latvecs.reshape(3,3,1),reps=(1,1,self.num_atoms))

        cart_pos = self.cartesian_pos[step,...].copy(order='F')

        self.reduced_pos[:,0] = inv_latvecs[0,0,:] * cart_pos[:,0] + \
                                inv_latvecs[0,1,:] * cart_pos[:,1] + \
                                inv_latvecs[0,2,:] * cart_pos[:,2] 
        
        self.reduced_pos[:,1] = inv_latvecs[1,0,:] * cart_pos[:,0] + \
                                inv_latvecs[1,1,:] * cart_pos[:,1] + \
                                inv_latvecs[1,2,:] * cart_pos[:,2] 
        
        self.reduced_pos[:,2] = inv_latvecs[2,0,:] * cart_pos[:,0] + \
                                inv_latvecs[2,1,:] * cart_pos[:,1] + \
                                inv_latvecs[2,2,:] * cart_pos[:,2] 
    
    # ----------------------------------------------------------------------------------------------

    def _unwrap(self):

        """
        if data at step 'step' have moved by more than 1/2 a box vector, shift back
        """

        _x = self.num_reps[0]
        _y = self.num_reps[1]
        _z = self.num_reps[2]

        delta = self.reduced_pos-self.reduced_pos_0
        dx = -_x * (delta[:,0] >= _x/2).astype(float) 
        dx += _x * (delta[:,0] < -_x/2).astype(float) 
        dy = -_y * (delta[:,1] >= _y/2).astype(float) 
        dy += _y * (delta[:,1] < -_y/2).astype(float) 
        dz = -_z * (delta[:,2] >= _z/2).astype(float) 
        dz += _z * (delta[:,2] < -_z/2).astype(float) 

        self.reduced_pos[:,0] += dx
        self.reduced_pos[:,1] += dy
        self.reduced_pos[:,2] += dz

    # ----------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    if len(sys.argv) != 6:
        msg = 'must give input hdf5 filename with cartesian coords,\n'
        msg += 'target hdf5 filename to write reduced coords,\n'
        msg += 'and supercell size along a, b, and c as args'
        exit(msg)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    num_reps = [ int(_) for _ in sys.argv[3:] ]

    converter = c_convert_to_reduced_coords(input_file,output_file,num_reps)

# --------------------------------------------------------------------------------------------------
