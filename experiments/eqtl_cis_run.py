# IMPORT GENERAL STUFF
import os
import scipy as SP
import h5py
import sys
import pdb

n_folds    = int(sys.argv[1])
min_fold   = int(sys.argv[2])
max_fold   = int(sys.argv[3])

# CREATION OF THE GROUP
group_name = 'drug_eqtl_cis'
print ""
command = "bgadd /%s" % group_name
print command
os.system(command)
n_jobs = max_fold
command = "bgmod -L %d /%s" % (n_jobs,group_name)
print command
os.system(command)
print ""

#Create temp dir
temp_folder   = '/Users/fatemehzamanzad/Desktop/limix/out/temp/%s'%group_name
if not os.path.exists(temp_folder):
    os.makedirs(temp_folder)

pdb.set_trace()

# GO! GO! GO!
for j in range(min_fold,max_fold):
    stdout_file = os.path.join(temp_folder,'stdout_%d_%d.txt'%(n_folds,j))
    stderr_file = os.path.join(temp_folder,'stderr_%d_%d.txt'%(n_folds,j))
    command  = "bsub -R \"select[gpfs]\" -g /%s " % group_name 
    command += "-o %s " % stdout_file
    command += "-e %s " % stderr_file
    command += "python eqtl_cis.py %d %d" % (n_folds,j)
    print command
    os.system(command)