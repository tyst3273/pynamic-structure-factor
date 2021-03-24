import numpy as np
import h5py


class compressor:

    """
    takes the files in 'files' and writes an hdf5 database named 'output_file'.

    contains:
    pos_data = [num_time_steps x num_atoms x xyz]
    time_step = [num_time_steps]
    atom_types = [num_atoms]
    box_bounds = [num_atoms x xlo,xhi,ylo,yhi,zlo,zhi]

    """

    def __init__(self,files,output_file,num_steps=0,num_atoms=0,debug=False):
        
        self.debug = debug
        self.num_steps = num_steps
        self.num_atoms = num_atoms
        self.num_files = len(files)

        # open output hdf5 files, for debug==True, write steps only
        if self.debug == True:
            self.output_file = h5py.File('debug.hdf5','w')
        else:
            self.output_file = h5py.File(output_file,'w')    

        # open all of the input fils
        self.files = []
        for f in files:
            self.files.append(open(f,'r'))

        # get number of atoms if not specified
        if self.num_atoms == 0:
            with open(files[0],'r') as tmp:
                for i in range(3):
                    tmp.readline()
                self.num_atoms = int(tmp.readline().strip())
        print('\n\t there are {:g} atoms\n'.format(self.num_atoms))

        # get number of steps if not specified
        if self.num_steps == 0:
            with open(files[0],'r') as tmp:
                self.num_steps = sum(1 for line in tmp)
            self.num_steps = self.num_steps//(self.num_atoms+9)
        self.num_steps = int(self.num_steps*self.num_files)
        print('\n\t there are {:g} steps\n'.format(self.num_steps))
        
        # can write new methods but this one is good
        self._parse_files()

        # close all of the input fils
        for f in self.files:
            self.files.close()
        
        # close output file
        self.output_file.close()

    
    def _parse_files(self):
        
        # initialize datasets. For debug==True, all but time_steps are place holders
        if self.debug == True:
            self.pos_data = self.output_file.create_dataset('pos_data',[1])
            self.atom_types = self.output_file.create_dataset('atom_types',[1])
            self.box_bounds = self.output_file.create_dataset('box_bounds',[1])
            self.time_steps = self.output_file.create_dataset('time_steps',[self.num_steps])
        else:
            self.pos_data = self.output_file.create_dataset('pos_data',
                                        [self.num_steps,self.num_atoms,3])
            self.atom_types = self.output_file.create_dataset('atom_types',
                                        [self.num_atoms])
            self.time_steps = self.output_file.create_dataset('time_steps',
                                        [self.num_steps])
            self.box_bounds = self.output_file.create_dataset('box_bounds',
                                        [self.num_steps,6])

        step = 0 # count the time steps across all files
        for f in self.files:
            print("\n\tnow reading positions ...\n")
            for i in range(self.num_steps//self.num_files):
                if step%100 == 0: # progress report
                    print('\t\tstep: {:g} out of {:g}\n'.format(step,self.num_steps))

                # this part skips headers and get the timestep
                f.readline()
                self.time_steps[step] = int(f.readline().strip())
                f.readline()
                f.readline()
                f.readline()

                # if debug==True, skip writing the data
                if self.debug == True:
                    for i in range(self.num_atoms+4):
                        f.readline()

                # otherwise, read the data and write it to the hdf5 database
                else:
                    self.box_bounds[step,[0,1]] = np.array(
                            f.readline().strip().split()).astype(float)
                    self.box_bounds[step,[2,3]] = np.array(
                            f.readline().strip().split()).astype(float)
                    self.box_bounds[step,[4,5]] = np.array(
                            f.readline().strip().split()).astype(float)
                    f.readline()
                
                    pos = np.zeros((self.num_atoms,3))
                    ids = np.zeros(self.num_atoms).astype(int)
                    types = np.zeros(self.num_atoms).astype(int)
                    if self.atom_types[0] == 0:
                        for atom in range(self.num_atoms):
                            tmp = f.readline().strip().split()
                            ids[atom] = int(tmp[0])
                            types[atom] = int(tmp[1])
                            pos[atom,:] = np.array(tmp[2:5])
                        self.atom_types[:] = np.copy(types[np.argsort(ids)[:]])
                        self.pos_data[step,:,:] = np.copy(pos[np.argsort(ids)[:],:])
                    else:
                        for atom in range(self.num_atoms):
                            tmp = f.readline().strip().split()
                            ids[atom] = int(tmp[0])
                            pos[atom,:] = np.array(tmp[2:5]).astype(float)
                        self.pos_data[step,:,:] = np.copy(pos[np.argsort(ids)[:],:])
                
                # update step count
                step = step+1



                 














