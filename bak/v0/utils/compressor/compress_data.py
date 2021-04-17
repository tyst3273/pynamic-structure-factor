from Compressor import compressor


debug = False

    
# name of the output files
output_file = 'pos.hdf5'
    
# get the input files 
files = ['posvel_196.dat','posvel_197.dat','posvel_198.dat',
                                    'posvel_199.dat','posvel_200.dat']

# print the names to check
print(files)
print(output_file)

# run my compression code
comp = compressor(files=files,output_file=output_file,
                num_steps=1000,num_atoms=48000,debug=debug)

# will exit after 1st ns if debug == True
if debug == True:
    exit()



