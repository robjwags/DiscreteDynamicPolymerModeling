function DefineNumbersOfTethersAndParticles

global Np Nt Ns N_Kuhn...
    N N_bead N_meso 

N = Np*Nt;      % is the number of total tethers for equilibration

% Define numbers of particles
N_backbone_beads = Nt+(Nt-1)*(N_Kuhn-1);
N_sidechain_beads = Nt*N_Kuhn*Ns;
N_bead = Np*(N_backbone_beads+N_sidechain_beads);
N_meso = Np*(Nt*Ns+Nt);

end