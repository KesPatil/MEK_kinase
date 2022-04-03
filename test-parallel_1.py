# Parallelized Python code to compute mean-difference location of the residues 
# over big data trajectories


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


start_time = time.time()
active_ref1 = '/Volumes/3KSPAT/from_2KSPAT/metad_convergence/wt/free_ener_wt/8us/mek_active_CA.pdb'
inactive_ref1 = '/Volumes/3KSPAT/from_2KSPAT/metad_convergence/wt/free_ener_wt/8us/mek_inactive_CA.pdb'
active_ref = mda.Universe(active_ref1)
inactive_ref = mda.Universe(inactive_ref1)

var_xtc = []
for i in range(0,10): # number of HILLS files
    myPath = '/Volumes/3KSPAT/from_2KSPAT/metad_convergence/wt/trajectory_files/traj_curr/traj_wt/'
    tifCounter = len(glob.glob1(myPath + 'trajectories_'+str(i),"*.xtc"))
    for j in range(0,tifCounter):
        var_xtc.append( myPath + 'trajectories_'+str(i) +'/' +glob.glob1(myPath + 'trajectories_'+str(i),"*.xtc")[j])

    
################################################################################################################ 
################################################################################################################
var_gro = '/Volumes/3KSPAT/from_2KSPAT/metad_convergence/wt/trajectory_files/active_mek_wt_protein.gro' 
u = mda.Universe(var_gro, var_xtc)


N_residue = 360


dnum = np.loadtxt('num.txt') 

comm.Barrier()
# compute the delta_r's
bume = np.zeros((len(u.trajectory),N_residue,3))

num_per_rank = len(u.trajectory) // size

lower_bound = rank * num_per_rank
upper_bound = (rank+1) * num_per_rank

print("This is core", rank, "I am working on trajectories from", lower_bound,"to",upper_bound-1,flush=True)


for i in range(lower_bound, upper_bound):
    u.trajectory[i]
    align.alignto(u, inactive_ref, select="protein and not(resid {}:{}) and not(resid {}:{})  and name CA".format(alpha_start,alpha_end,actloop_start,actloop_end), weights="mass")

    a2helixloopi   = u.select_atoms("protein and name CA".format(alpha_start,alpha_end,actloop_start,actloop_end))
    a2CAhelixloopi  = a2helixloopi.positions
    
    bume[i] = a2CAhelixloopi - dnum
comm.Barrier()

# the sum of all cores will be held in the rank 0 core
if comm.rank==0:
    bume_mid = np.zeros_like(bume)
else:
    bume_mid = None

# use MPI to get the totals 
comm.Reduce(
    [bume, MPI.DOUBLE],
    [bume_mid, MPI.DOUBLE],
    op = MPI.SUM,
    root = 0
)
comm.Barrier()

if comm.rank==0:
    bume_reshaped = bume_mid.reshape(bume_mid.shape[0],-1)
    np.savetxt("bume.txt",bume_reshaped)
    print("--ms---", int((time.time() - start_time) * 1000),flush=True)


