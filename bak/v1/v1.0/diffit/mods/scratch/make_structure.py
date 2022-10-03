import numpy as np



class rutile_frenkel:

    def __init__(self):

        """
        class for making rutile supercell with frenkel defects
        """

        self.prim_lens = np.array([4.6235537460367500,  4.6235537460367500,  2.9786839831313254])
        self.prim_pos = np.array([[0.5000000000000000,  0.5000000000000000,  0.5000000000000000],
                                  [0.0000000000000000,  0.0000000000000000,  0.0000000000000000], 
                                  [0.1953400114833092,  0.8046599885166907,  0.5000000000000000],  
                                  [0.8046599885166907,  0.1953400114833092,  0.5000000000000000],  
                                  [0.3046599885166907,  0.3046599885166907,  0.0000000000000000],   
                                  [0.6953400114833093,  0.6953400114833093,  0.0000000000000000]]) 
        self.prim_types = np.array([1,1,2,2,2,2],dtype=int)
        self.prim_charges = np.array([2.196,-1.098])
        self.num_prim_atoms = 6

    def make_sc(self,sc_reps=[1,1,1]):

        """
        get num_reps copies of primitive pos (sc_pos) and array of 
        vectors to add to copies to create supercell (sc_shift)
        """

        self.sc_reps = np.array(sc_reps,dtype=int)
        self.sc_lens = self.prim_lens*self.sc_reps
        self.num_reps = np.prod(self.sc_reps)
        self.num_sc_atoms = self.num_prim_atoms*self.num_reps
        self.sc_pos = np.tile(self.prim_pos.reshape(1,self.num_prim_atoms,3),reps=(self.num_reps,1,1))
        self.sc_types = np.tile(self.prim_types.reshape(1,self.num_prim_atoms),reps=(self.num_reps,1))

        _x, _y, _z = np.meshgrid(np.arange(self.sc_reps[0]),
                                 np.arange(self.sc_reps[1]),
                                 np.arange(self.sc_reps[2]),indexing='ij')
        _x = _x.flatten(); _y = _y.flatten(); _z = _z.flatten()

        self.sc_shift = np.array((_x,_y,_z),dtype=int).T
        self.sc_shift.shape = [self.num_reps,1,3]
        self.sc_shift = np.tile(self.sc_shift,reps=(1,self.num_prim_atoms,1))

        self.sc_pos += self.sc_shift
        self.sc_pos[:,:,0] /= self.sc_reps[0]
        self.sc_pos[:,:,1] /= self.sc_reps[1]
        self.sc_pos[:,:,2] /= self.sc_reps[2]

    def make_random_oxygen_frenkels(self,defect_fraction=0.05,num_defects=None):
        
        """
        'defect_fraction' is number of UNITCELLS with defects

        self.prim_pos = np.array([[0.5000000000000000,  0.5000000000000000,  0.5000000000000000],
                                  [0.0000000000000000,  0.0000000000000000,  0.0000000000000000],
                                  [0.1953400114833092,  0.8046599885166907,  0.5000000000000000],  # O1
                                  [0.8046599885166907,  0.1953400114833092,  0.5000000000000000],  # O2
                                  [0.3046599885166907,  0.3046599885166907,  0.0000000000000000],  
                                  [0.6953400114833093,  0.6953400114833093,  0.0000000000000000]]) 

                                # O1 => (0.00,0.25,0.50) or (0.75,1.00,0.50)
                                # O2 => (0.25,0.00,0.50) or (1.00,0.75,0.50)

        """

        if num_defects is None:
            if defect_fraction <= 0 or defect_fraction > 1:
                exit('fuck!')
            num_defects = int(self.num_reps*defect_fraction)
        else:
            if num_defects <= 0 or num_defects > self.num_reps:
                exit('fuck!')

        _uc, _s = np.meshgrid(np.arange(self.num_reps),np.arange(2,6),indexing='ij')
        _uc = _uc.flatten(); _s = _s.flatten()
        sites = np.array((_uc,_s),dtype=int).T
        num_sites = sites.shape[0]

        defects = np.arange(num_sites)
        np.random.shuffle(defects)
        defects = defects[:num_defects]

        d1 = 0.19534001; d2 = 0.47465999

        disp_O1 = np.array([[-d1,-d2, 0],
                            [ d2, d1, 0]])
        disp_O1[:,0] /= self.sc_reps[0]
        disp_O1[:,1] /= self.sc_reps[1]
        disp_O1[:,2] /= self.sc_reps[2]

        disp_O2 = np.array([[ d1, d2, 0],
                            [-d2,-d1, 0]])
        disp_O2[:,0] /= self.sc_reps[0]
        disp_O2[:,1] /= self.sc_reps[1]
        disp_O2[:,2] /= self.sc_reps[2]

        disp_O3 = np.array([[ d1,-d2, 0],
                            [-d2, d1, 0]])
        disp_O3[:,0] /= self.sc_reps[0]
        disp_O3[:,1] /= self.sc_reps[1]
        disp_O3[:,2] /= self.sc_reps[2]

        disp_O4 = np.array([[-d1, d2, 0],
                            [ d2,-d1, 0]])
        disp_O4[:,0] /= self.sc_reps[0]
        disp_O4[:,1] /= self.sc_reps[1]
        disp_O4[:,2] /= self.sc_reps[2]

        for ii in range(num_defects):

            ind = defects[ii]

            disp = np.arange(2)
            np.random.shuffle(disp)

            if sites[ind,1] == 2:
                disp = disp_O1[disp[0],:]
            elif sites[ind,1] == 3:
                disp = disp_O2[disp[0],:]
            elif sites[ind,1] == 4:
                disp = disp_O3[disp[0],:]
            elif sites[ind,1] == 5:
                disp = disp_O4[disp[0],:]

            self.sc_pos[sites[ind,0],sites[ind,1],:] += disp
       
    def make_random_titanium_frenkels(self,defect_fraction=0.05,num_defects=None):

        """
        'defect_fraction' is number of Ti sites that should be made into Frenkels
        """

        if num_defects is None:
            if defect_fraction <= 0 or defect_fraction > 1:
                exit('fuck!')
            num_defects = int(self.num_reps*defect_fraction)
        else:
            if num_defects <= 0 or num_defects > self.num_reps:      
                exit('fuck!')

        _uc, _s = np.meshgrid(np.arange(self.num_reps),np.arange(2),indexing='ij')
        _uc = _uc.flatten(); _s = _s.flatten()
        sites = np.array((_uc,_s),dtype=int).T
        num_sites = sites.shape[0]

        defects = np.arange(num_sites)
        np.random.shuffle(defects)
        defects = defects[:num_defects]

        displacement = np.array([[-0.5, 0.0, 0.0],
                                 [ 0.0,-0.5, 0.0],
                                 [ 0.5, 0.0, 0.0],
                                 [ 0.0, 0.5, 0.0]])
        displacement[:,0] /= self.sc_reps[0]
        displacement[:,1] /= self.sc_reps[1]
        displacement[:,2] /= self.sc_reps[2]

        for ii in range(num_defects):

            ind = defects[ii]

            disp = np.arange(4)
            np.random.shuffle(disp)
            disp = displacement[disp[0],:]
            self.sc_pos[sites[ind,0],sites[ind,1],:] += disp

    def write_poscar(self,file_name='POSCAR'):
            
        """
        write a VASP 'POSCAR' file for calculating/visualizing with VESTA
        """

        pos = self.sc_pos.reshape(self.num_sc_atoms,3)
        types = self.sc_types.reshape(self.num_sc_atoms)
        inds = np.argsort(types)
        
        num_ti = np.count_nonzero(types == 1)
        num_o = np.count_nonzero(types == 2)

        types = types[inds]
        pos = pos[inds]

        with open(file_name,'w') as f_out:
            _= 0.0
            f_out.write(f'auto generated\n 1.0\n')
            f_out.write(f'  {self.sc_lens[0]:10.7f}  {_:10.7f}  {_:10.7f}\n')
            f_out.write(f'  {_:10.7f}  {self.sc_lens[1]:10.7f}  {_:10.7f}\n')
            f_out.write(f'  {_:10.7f}  {_:10.7f}  {self.sc_lens[2]:10.7f}\n')
            f_out.write(f' Ti O \n')
            f_out.write(f'   {num_ti:g}  {num_o:g}\nDirect\n')
            for ii in range(self.num_sc_atoms):
                f_out.write(f' {pos[ii,0]:10.9f}  {pos[ii,1]:10.9f}  {pos[ii,2]:10.9f}\n')

    def write_lammps(self,file_name='lammps.pos'):

        """
        write a lammps input file
        """

        pos = self.sc_pos.reshape(self.num_sc_atoms,3)
        types = self.sc_types.reshape(self.num_sc_atoms)
        inds = np.argsort(types)

        num_ti = np.count_nonzero(types == 1)
        num_o = np.count_nonzero(types == 2)

        types = types[inds]
        pos = pos[inds]

        pos[:,0] *= self.sc_lens[0]
        pos[:,1] *= self.sc_lens[1]
        pos[:,2] *= self.sc_lens[2]

        with open(file_name,'w') as f_out:
            f_out.write('# auto generated\n\n')
            f_out.write(f'{self.num_sc_atoms} atoms\n')
            f_out.write('2 atom types\n\n')
            f_out.write(f'0.0  {self.sc_lens[0]:10.7f} xlo xhi\n')
            f_out.write(f'0.0  {self.sc_lens[1]:10.7f} ylo yhi\n')
            f_out.write(f'0.0  {self.sc_lens[2]:10.7f} zlo zhi\n')
            f_out.write('\nAtoms \n')
            for ii in range(self.num_sc_atoms):
                f_out.write(f'\n {ii+1:5g} {types[ii]:2g} {self.prim_charges[types[ii]-1]: 6.4f}' \
                    f' {pos[ii,0]:12.9f} {pos[ii,1]:12.9f} {pos[ii,2]:12.9f}')



# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    struct = rutile_frenkel()

    struct.make_sc(sc_reps=[10,10,12])
    struct.make_random_titanium_frenkels(defect_fraction=0.10)
    struct.write_lammps('ti_frenkel_0.10.pos')
    struct.write_poscar()

    struct.make_sc(sc_reps=[10,10,12])
    struct.make_random_oxygen_frenkels(defect_fraction=0.05)
    struct.write_lammps('o_frenkel_0.05.pos')
#    struct.write_poscar()



