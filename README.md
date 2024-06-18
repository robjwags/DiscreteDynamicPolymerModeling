# DiscreteDynamicPolymerModeling
Discrete dynamic polymer modeling

# Main README for "A foundational framework for the mesoscale modeling of dynamic elastomers and gels"
Robert J. Wagner
June 12th, 2024


# DESCRIPTION

 This folder contains the necessary codes written in MATLAB2022a required for the following tasks in relation to this manuscript:
 1. Generating initial topology files (.txt) and input files (.in) for LAMMPS
 2. Running LAMMPS using the <sytem(command)> line in MATLAB2022a and .in/.txt files from step 1
 3. Compiling the raw .dump files output from LAMMPS and storing them as .m structures for post-processing
 4. Post-processing and analyzing the compiled data for each section of the manuscript

# REQUIRED SOFTWARE

 This source code requires install of MATLAB2022a or later
 The following LAMMPS packages must also be built to run these studies using cpus (see INSTALL instructions below for install guide):
	A) BROWNIAN
	B) TNT
	C) MOLECULE
 Following might be potentially useful in future:
	D) DPD-BASIC
	E) EXTRA-PAIR
	F) GRANULAR
	G) RIGID

# DIRECTORY ORGANIZATION & STUDY PURPOSES

 Several distinct studies with unique particle topologies and bonding conditions were conducted in the scope of this work. Each study has its own folder. Distinct studies are orgranized as listed below. Visit each subdirectory to find a study-specific README instructing code usage:

 1. SingleChain/ (Section 3.1 of manuscript)
	A. ForceStudy/  (VALIDATION STUDY)
		i.	 Description: A single bead-spring chain pulled apart slowly and held at specified end-to-end lengths for some long duration. 
		ii.	 Measured: Total chain force end-to-end length
		iii. Goal: To validate bead-spring force against the analytical relation for Langevin chains, and calibrate the finitely extensible, nonlinear elastic potential used for Kuhn segments
	
	B. DiffusionStudy/  
		i. 	 Description: A single chain with one end tethered (i.e., fixed in space) and the other (and intermediate portions of the chain, for bead-spring models) subject to Brownian diffusion.
		ii.	 Measured: End-to-end distributions and MSD MSD (net, radial, and tangential) of an ensemble of many such chains with respect to their anchoring node's position.
		iii. Goal: to confirm agreement between mesoscale, bead-spring, and Rouse models in the exploratory sub-diffusion of open end groups attached to tethered chains.

 2. SingleBond/ (Relates to Appendix E of manuscript)
	A. Check_ka/ (VALIDATION STUDY)
		i. 	 Description: Pairs of adjacent, fixed stickers within attachment cutoff distance, b, of each other that are allowed to bind and unbind in time.
		ii.	 Measured: ka and kd
		iii. Goal: To confirm that the attachment and detachment rates set a priori match those that are measured as outputs from the code.

	B. SingleKinetics/ (Relates to Section 3.2 of manuscript)
		i. 	 Description: A single tethered chain whose fixed end is stable and whose free end is subject to Brownian diffusion. A fixed sticker node is also placed distance, d, from the fixed stable end of the tethered chain such that it may bind/unbind with the open sticker of the chain.
		ii. Measured: ka as a function of number of Kuhn segments (N), separation distance (d), and activation energy (epsilon_a)
		iii. Goal: To confirm agreement between the pairwise rates of attachment and detachment between the mesoscale and bead-spring models. Also to develop scaling theory for pairwise attachment kinetics that might be useable for implicit form of attachment probability.

	C. PairwiseKinetics/ (Relates to Appendix F of manuscript)
		i. Description: Two tethered chains whose fixed ends are stable nodes separated by a distance, d, and whose free end are stickers subjected to Brownian diffusion. The free ends of the chains may bind/unbind.
		ii. Measured: ka as a function of number of Kuhn segments (N), separation distance (d), and activation energy (epsilon_a)
		iii. Goal: To confirm agreement between the pairwise rates of attachment and detachment between the mesoscale and bead-spring models and confirm applicability of scaling theory to two-chain (as opposed to one-chain) scenario.  

 3. BondExchange/ (Relates to Section 3.3 of manuscript) 
		i. 	 Description: An ensemble of chains, each with one end fixed in space, are allowed to undergo tethered diffusion and bond/unbond at their distal ends. The fixed ends are arranged in a 3D grid with prescribed spacing separation, d.
		ii.  Measured: ka, kexc, krpt, kd, tau_a, tau_d, tau_exc, tau_rpt, tau_rnm, fa, and fd as functions of fixed node separation distance, d. 
		iii. Goal: To evaluate agreement in the models' bond kinetics and partner exchange kinetics as a function of chain concentration and chain length.

 4. NetworkScale/StressRlx/2024_03_01_Main (Relates to Sections 4 and 5 of manuscript)
	See the dedicated README within this directory for a detailed explanation of contents

# INSTALL INSTRUCTIONS

1) If using Windows, install WSL, update everything and open WSL

2) Navigate to directory where desire to build lammps
	Ex. <cd /mnt/c/Users/rjwag/Documents/>

3) Clone dynamics branch 
	<git clone --branch dynamic https://github.com/slamont1/lammps.git>

4) Go into <cd lammps/>

5) Navigate to the source code directory and TNT package
	<cd src/TNT/>

6) Copy bond_pade.cpp and bond_pade.h to the MOLECULES package
	<cp bond_pade.cpp ../MOLECULE>
	<cp bond_pade.h ../MOLECULE>

7) Make a build directory in lammps folder and navigate to it
	<cd ..>
	<mkdir build>
	<cd build/>

8) Build lammps with apropriate packages and with MPI capability
	<cmake ../cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran -D PKG_MOLECULE=yes -D PKG_EXTRA-MOLECULE=yes -D PKG_TNT=yes -D PKG_FEP=yes -D PKG_BROWNIAN=yes>

9) Make lammps
	<make -j N> where N is number of processecors for make
