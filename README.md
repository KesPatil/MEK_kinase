

# MEK_kinase
This GitHub repo has folders required for setting up the molecular dynamics simulations, metadynamics simulations, and  post-processing codes for studies on MEK kinase and mutations


<p align="justify">
README:  There are three folders F_unbiased that includes the files to setup and perform Unbiased Molecular Dynamics, F_metad contains files to setup and perform metadynamics, F_post_processing_codes contains files to perform post processing <br />
</p>

## Folders

***F_unbiased***: This folder has the input gromacs .mdp scripts required to perform Unbiased Molecular Dynamics and the active and the inactive structures of ALK <br />
1. The topology files .top are generated using Gromacs function : "pdb2gmx" and we implement charmm27 forcefield <br />
2. The outputs (.gro, .cpt, .top, .tpr) of the unbiased MD simualtion are inputs to Metadynamics. The active and inactive MEK structures at 101 ns of MD simulation were chosen as reference in the Metadynamics run. These structures are included in F_metad as  active_mek_wt_protein.pdb and  inactive_mek_wt_protein.pdb <br />

***F_metad***: This folder has the input scripts required to perform Metadynamics and also the output free energy file that is obtained by summing over HILLS files <br />
1. plumed.dat - a plumed script that needs to be run for performing Metadynamcs <br />
2. mek_active_CA.pdb - the reference active structure, one of the two collective variables is RMSD from this structure <br />
3. mek_inactive_CA.pdb - the reference inactive structure, one of the two collective variables is RMSD from this structure <br />

PS note: Additionally, you would need the usual gromacs file:  .mdp file, .top file, .cpt file, .tpr file which are the outputs of unbiased MD <br />


4. fes_8us_wt_mintozero.dat - The free energy file which has free energy at every coordinate on CV space obtained by summing over HILLS files. <br />

***F_post_processing_codes***: This folder contain the python scripts used for post-processing trajectories data to calculate boltzmann weighted correlations matrix, sasa, hydrogen bond occupancies and plot free energy landscapes, extract structures from zones and check convergence of metadynamics free energy zones  </br> 



## Forcefield implemented

Our MD and metadynamics simulations use charmm27.ff from GROMACS 5.0.7

## Procedure
### 1. MD simulation in GROMACS <br />
<p align="justify">
The MD simulations are performed in GROMACS. The paper mentions Biophyscode which actually is a wrapper on Gromacs. Please feel free to try out Biophyscode, a creation from Radhakrishnan lab. https://biophyscode.github.io. Please note: Following is the walk-through on how to  run a simulation using Biophyscode:

"https://biophyscode.github.io/molecular_dynamics_lab/". <br />

All mutations are introduced using BioPhysCode Automacs routine based on MODELLER. Please see MODELLER here: <br />

https://salilab.org/modeller/ <br />



All the requisite files needed to setup and run the simulations are  included in the F_unbiased. We followed Bevan Lab Tutorials: "Lysozyme in water" example for equilibration and production <br />
http://www.mdtutorials.com/gmx/lysozyme/index.html
</p>

### 2. Enhanced sampling simulation of Metadynamics using PLUMED <br />
<p align="justify">
We patch GROMACS 5.0.7 with PLUMED 2.3.5 and use multiple walker (num_of_walkers=10) metadynamics. For installation procedure of PLUMED and the subsequent patching with GROMACS, see here: https://www.plumed.org/doc-v2.6/user-doc/html/_installation.html <br />
</p>

<p align="justify">
In the folder F_metad we have the scripts required to run the metadynamics. Please note: Metadynamics has to be run using the output .gro and .cpt files of the unbiased MD  simulations. The folder contains the .mdp file, plumed.dat and the reference active and inactive MEK alpha Carbons structure files <br />
</p>

<p align="justify">
plumed.dat metadynamics script is run in GROMACS 5.0.7 patched with PLUMED 2.3.5 using multiple walker (num_of_walkers=10) metadynamics. We set in the PLUMED script, the energy in kcal/mol and length in Å. The parameters used in this study to perform Well-Tempered multiple walker metadynamics (WTMD) are:bias factor γ = T +∆T/T = 20, height = 0.6 and pace = 500. <br />
<br />
</p>

<p align="justify">
Aggregate of 8us metadynamics run in GROMACS 5.0.7 patched with PLUMED 2.3.5 using multiple walker (num_of_walkers=10) metadynamics until the convergence criterion is met for the zones of interest (See section xxx) <br />
</p>

<p align="justify">
The output of the metadynamics run are the HILLS files. We will have ten of them since we used ten walkers. To get the free energy file from the HILLS use <b>plumed sum_hills</b> action of PLUMED (https://www.plumed.org/doc-v2.5/user-doc/html/sum_hills.html) <br />
</p>

### 3. Description of the post-processing codes <br />

Analysis of the trajectories was done using python and mostly MDanalysis package in python: https://www.mdanalysis.org  <br />

 The folder F_post_processing contains eleven codes: <br />

 + `Aloop_helicity.ipynb` is the code to compute hbonds that contribute to making the partial helix in the activation loop of kinase <br />
 + `Analysis_weighted_correlation.ipynb` is the code that computes boltzmann weighted SASA, Hbond and their log ratio plot  <br />
 + `Boltz_corr_plotter.ipynb` is the code to plot the boltzmann weighted correlation matrices  <br />
 + `Dihedral_compute.ipynb` is the code to compute the distribution of the dihedral angle representing the DFG flip  <br />
 + `FELandscape_plotter.ipynb` is the code to plot the free energy landscape  <br />
 + `FELandscape_to_zone_structure_extractor.ipynb` is the code to extract structures from chosen zones on the free energy landscape   <br />
 + `Hierarchical_cluster.ipynb` is the code to plot the hierarchical cluster differentiating dynamics of mutated MEK systems  <br />
 + `KE_salt_bridge.ipynb` is the code to compute the distribution of the KE salt bridge distance  <br />
 + `boltz_corr1.py` is the parallelized code to compute boltzmann weighted correlations  <br />
 + `service.sh` is bash script to extract protein trajectory from system trajectory  <br />
 + `zone_highlighter.ipynb` is the code to identify specific zones on the free energy landscape  <br />


## Citations

If you have any suggestions or queries please feel free to reach out at : patilk@seas.upenn.edu  <br />
If you found the above scripts and/or codes helpful in your work, please cite: <br />
1. Jordan, E. Joseph, et al. "Computational algorithms for in silico profiling of activating mutations in cancer." Cellular and Molecular Life Sciences 76.14 (2019): 2663-2679.
2. Patil, Keshav, et al. "Computational studies of anaplastic lymphoma kinase mutations reveal common mechanisms of oncogenic activation." Proceedings of the National Academy of Sciences 118.10 (2021).
3.
</p>

