# calculate Boltzmann weighted correlation by implementing parallelization

from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


import numpy as np
import time
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import XTC, GRO
import MDAnalysis.analysis.rms
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.analysis import align
import glob


#####-------100ns------#####
alpha_start = 107
alpha_end = 125
actloop_start = 208
actloop_end = 233
#############################
pwth = '/scratch/rradhak1/patilk/MEK/e203k/traj_e203k/'
start = time.time()
active_ref1 = '/scratch/rradhak1/patilk/MEK/e203k/traj_e203k/mek_active_CA.pdb'
inactive_ref1 = '/scratch/rradhak1/patilk/MEK/e203k/traj_e203k/mek_inactive_CA.pdb'
active_ref = mda.Universe(active_ref1)
inactive_ref = mda.Universe(inactive_ref1)

var_xtc = []
i = 9
print(i)
myPath = pwth #'/Volumes/3KSPAT/from_2KSPAT/metad_convergence/wt/trajectory_files/traj_curr/traj_wt/'
tifCounter = len(glob.glob1(myPath + 'trajectories_'+str(i),"*.xtc"))
for j in range(0,tifCounter):
    var_xtc.append( myPath + 'trajectories_'+str(i) +'/' +glob.glob1(myPath + 'trajectories_'+str(i),"*.xtc")[j])

    
################################################################################################################ 
################################################################################################################
var_gro =  '/scratch/rradhak1/patilk/MEK/e203k/traj_e203k/mek_e203k_protein.gro' 
u = mda.Universe(var_gro, var_xtc)

print(len(u.trajectory))
N_residue = 360


num = np.loadtxt('num_e203k.txt')
partition_func = np.loadtxt('partition_func_e203k.txt')

print("Computing delta r's")
start_time = time.time()
# compute the delta_r's
nume = np.zeros((len(u.trajectory),N_residue,3))
for i in range(0,len(u.trajectory)):
    u.trajectory[i]
    align.alignto(u, inactive_ref, select="protein and not(resid {}:{}) and not(resid {}:{})  and name CA".format(alpha_start,alpha_end,actloop_start,actloop_end), weights="mass")

    a2helixloopi   = u.select_atoms("protein and name CA".format(alpha_start,alpha_end,actloop_start,actloop_end))
    a2CAhelixloopi  = a2helixloopi.positions
    
    nume[i] = a2CAhelixloopi - num

print("Done with delta r's")

comm.Barrier()

num_per_rank = len(u.trajectory) // size

lower_bound = rank * num_per_rank
upper_bound = (rank+1) * num_per_rank

print("This is core", rank, "I am working on trajectories from", lower_bound,"to",upper_bound-1,flush=True)

C = np.zeros((N_residue, N_residue))
combo = np.zeros((len(u.trajectory),N_residue, N_residue,3))

for i in range(0,N_residue):
    for j in range(0,N_residue):
        numer = np.zeros(len(u.trajectory))
        deno1 = np.zeros(len(u.trajectory))
        deno2 = np.zeros(len(u.trajectory))
        
        
        for k in range(lower_bound, upper_bound):
            numer[k] = np.dot(nume[k][i], nume[k][j]) * partition_func[k]
            deno1[k] = np.dot(nume[k][i], nume[k][i]) * partition_func[k]
            deno2[k] = np.dot(nume[k][j], nume[k][j]) * partition_func[k]
            combo[k][i][j][0] = numer[k] 
            combo[k][i][j][1] = deno1[k] 
            combo[k][i][j][2] = deno2[k] 
    print(i) 
  

      
combo_sum = np.sum(combo, axis=0)

np.savetxt(str(comm.rank) + 'combo_sum_numer.txt' , combo_sum[:,:,0])
np.savetxt(str(comm.rank) + 'combo_sum_deno1.txt' , combo_sum[:,:,1])
np.savetxt(str(comm.rank) + 'combo_sum_deno2.txt' , combo_sum[:,:,2])

end_time = time.time()
t = end_time - start_time
print("------", t)             
    
