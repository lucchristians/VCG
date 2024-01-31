# Overview

This is a repository for scripts used in the publication "" found at

## Software used:

1. Python library information is found in the "py_info.txt" file. Important libraries to include are numpy, networkx, mdtraj, openmscg, and matplotlib.

2. Gromacs was used to simulation atomistic simulations. Our versions of gromacs were compiled with plumed enabled.

3. Lammps was used to run coarse grained simulations. Included packages are shown in the "lammps_info.txt" file.

4. The grossfield WHAM code was used to generate PMFs from umbrella sampling data.

## Methods by directories

0_setup: Setting up and simulating atomistic data for protein Q. The scripts are used in the following order:
1. bash AAsim_setup.bash
2. bash submit_sims.bash #set up for use on a SLURM enabled HPC

1_processing: AA processing shows how RMSD was used to check the progress of the simulation towards equilibrating and the set up of the CA-mapped CG represnetation used to train the REM model. In the "AA_processing" directory the scripts are used in teh following order:
1. bash ca_setup_and_rmsd.bash
2. vmd # convert concatenated_CA.xtc to lammpstrj format.
3. python make_vs_pdb.py CA.lammpstrj vs_frame.lammpstrj
4. python lammps2pdb.py vs_frame.lammpstrj

The "CG_processing" subdirectory houses the scripts used to generate the lammps simulation input files for the CG-mapped version of the atomistic system. This is all completed in the "framework.bash" file. Further processing was done withe the "replace_with_r0.py" script to better align virtual site bonds r0 distance to the peak position of the coarse-grained mapped atomistic to allow for better convergence.

2_rem_optimization: Using the CG model produced in the "CG_processing directory, relative entropy minimization can be used to optimize the gaussian and bonded potentials for the virtual sites. This is accomplished using the "rem_iterate.bash" script, which houses the whole REM process. this can also be submitted to a SLURM enabled HPC using the "sim_submit.sbatch" script. remaining scripts were used to check the progress and efficacy of convergence.

3_assembly: Scripts used to set up, run, and analyze the assembly of Q protein. First a CG-mapped monomer with virtual sites needs to be isolatd in a box large enough to represent a 10 mM system (an example 'outfile.data is included to showwhat such a input file should look like. the is replicated using the 'replicate' command in lammps to make a roughly square-shaped environment. In the case of protein Q, the HP sites can be added using the "define_hp_vs.py" script. HP site parametrization was accomplished using the "invbolt_henm_fit.py" script. The simulation it self can be run using a script like the "assemb_submit.sbatch" script. Clusters were graphically defined for the assembled structures using the find_cluster.py script. the process for cluster analysis was completed in the following order:
1. bash multicluster.bash
2. multiconcat.bash
3. bash plots.bash 

4_cg_pmf_calc: scripts used to generate a pmf for the CG model. typical inputs are included in the PMF. "run_setup.bash" was used as the framework for creating the files needed to conduct the sampling over a number of windows in a given range and in a particular axial direction. "multirun.bash" was used to run the sampling simulation. The "intra_excl.py" and the "intra_plot.py" are used ot generate the CG pmfs from these simulations and compare them to an atomistic PMF produced in the '1.5_us' directory.

1.5_us: Scripts used to run the atomistic umbrella sampling simulations. windows were setup using pull sims. the frame work for setting up and conducting these sims is in the "pull_sim_framework.bash" script. windows were set up using the "us_framework.bash" script. by modifying the path in the "us_sim_mps_base.sbatch" script different windows can be run. WHAM was accomplished in these steps:
1. python make_2dmetafile.py {rxn_coord1} {rxn_coord2} {temp}
2. bash run2dwham.bash
3. python mfep_2dfes.py
4. python min_plot_fig.py
This would take the 2D restraint and represent them in 1D along hte minimum energy pathway.
