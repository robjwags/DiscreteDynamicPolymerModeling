% Branched newtork topolgoy generator for input to LAMMPS via read_data cmd
% Robert J. Wagner
% 2022-11-29

% This version generates the x-links first and then connects

function GenerateTopology(Package,OverrideTop,OverrideIn,...
    CF,OF,TD,UnitType,LC,DC,DTF)

global dims ShowFiguresDebug Np...
    Ns N OldOrNew PlotStickers MakeMovie...
    TurnOnDynamics TopologyFileName...
    SIOrNormalized InputFileName ConstantsFileName...
    LengthConversion DamperConversion CurrentFolder OutputFolder dtFact

%% Set Basic Toggles
OldOrNew = 1;           %Different looping options for two diff. connectivity algorithms
ShowFiguresDebug = 1;   %1 to plot intermediate figures for debugging
TurnOnDynamics = TD;
SIOrNormalized = UnitType;
LengthConversion = LC;
DamperConversion = DC;
% BeadSpringOrMeso = BSOM;
CurrentFolder = CF;
OutputFolder = OF;
dtFact = DTF;

N_Kuhn_out = [];

%% Unpack Swept Input Parameters
Sample = Package(1);    %Sample index
Np = Package(2);        %Number of molecules
Ns = Package(3);        %Number of stickeres per tether site
N = Np*Ns;              %Total numbber of permanent crosslinks (i.e., tether sites)
T = Package(4);         %Temperature
EdgeFactor = Package(5);    %Contour length of chains
ka = Package(6);        %Activation energy of association
kd = Package(7);        %Activation energy of dissociation
f0 = Package(8);       %Force sensitivity to dissociation
dt = Package(9);       %timestep size
damp = Package(10);     %damping coefficient in units [mass/time]
% DiffCoeff = Package(11);
N_Kuhn = Package(12);
b = Package(13);
Separation = Package(14);

%% Callout Input Script
% Initializes non-sweeping parameters
dims = 3;   %Dimensionality of system
InputScript(EdgeFactor,Sample,Np,N,ka,kd,f0,dt,damp,N_Kuhn,b,Separation);

%% Make diriectories and set filenames
SetDirAndFileNames;

%% Convert all to correct units
kb =1;
Tab = table(Np,Ns,N_Kuhn,b,ka,kd,f0,kb,T,dt,damp);
writetable(Tab,ConstantsFileName);

%%%%%%%%%%%%%%%%%%%%%%%% CHANGE LOOPING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
if ~isfile(TopologyFileName) || OverrideTop==1

    %% Size the domain
    % Determines correct domain size based on no. xls and nominal density
    SizeTheDomain(Np,N_Kuhn,b,Separation);

    %% Initialize the Nodes
    InitializeNodes;

    figure(1); clf; hold on
    PlotStickers=1;
    if ShowFiguresDebug==1
        PlotNetwork
    end

    %% Output Topological Data
    SaveMatlabTopologyData;
    mass = 0.05;    %arbitrary
    SaveTopologyDataLAMMPS(N_Kuhn,kb,T,b,mass);

    MakeMovie=0;
    if MakeMovie==1
        MakeMovies(N_Kuhn);
    end

end

if ~isfile(InputFileName) || OverrideIn==1
    if isempty(N_Kuhn_out)
        dat = readtable(ConstantsFileName);
        dat = table2array(dat);
        Np = dat(1); Ns = dat(2);
        N_Kuhn_out = dat(3); b_out = dat(4);
        ka_out = dat(5); kd_out = dat(6);
        f0_out = dat(7); T_out = dat(9);
        dt_out = dat(10); 
    end
    PrintInputFiles(N_Kuhn_out,b_out,ka_out,kd_out,f0_out,T_out,...
        dt_out,damp,CurrentFolder);
end

% toc
close all
clear
clc
clear global

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PrintInputFiles(N,b,ka0,kd0,f0,T,dt,damp,CurrentFolder)

global TopologyFileName OutputAtom_loc OutputBond_loc InputFileName...
    TurnOnDynamics SIOrNormalized dtFact

TurnOnLJ = 0;       % 0 for no volume exclusion (do not need repulsive interactions
damp_eff = damp;
Ldyn = b;               % Attachment cutoff length - smaller for mesoscale

% Initialize string
Assembled = [];

% Define units
L = '####################################'; Assembled = [Assembled,L];
L = '\n# DEFINE UNITS'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
if SIOrNormalized==0
    L = '\nunits        si'; Assembled = [Assembled,L];
elseif SIOrNormalized==1
    L = '\nunits        lj'; Assembled = [Assembled,L];
else
    L = '\nunits        nano'; Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

kb = 1.38e-23;      % Note, this doesn't actually set anything

% Define lengthscales
Nb = N*b;           %chain length
Nb2 = Nb*b;         %N*b^2
r0 = sqrt(N)*b;     %mean random walk chain length
Kstiff = 3*kb*T/Nb2;     %chain stiffness (1/2 multiple to capture effective damper)
ftol = Kstiff*b/10000;  %equilibration force tolerance
Lmax = 2*b;      %max LJ cutoff
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# DEFINE LENGTH & FORCE SCALES'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Defining length scales'; Assembled = [Assembled,L];
L = ['\nvariable    b       equal ',num2str(b,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    N       equal ',num2str(N,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    Nb      equal ',num2str(Nb,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    Nb2     equal ',num2str(Nb2,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    r0      equal ',num2str(r0,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    Kstiff  equal ',num2str(Kstiff,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    ftol    equal ',num2str(ftol,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    Lmax    equal ',num2str(Lmax,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Define Thermodynamic properties
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# DEFINE LENGTH & FORCE SCALES'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
L = ['\nvariable    T       equal ',num2str(T,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    damp    equal ',num2str(damp_eff,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Define repulsive properties, if any properties
sig_r0 = 0.25*r0;           % backbone LJ length scale
eps_r0 = r0*Kstiff/10;      % backbone LJ energy parameter
lambda_r0 = 0.5;            % for soft core LJ potential
% sig_b = 2*b;                % backbone LJ length scale
eps_b = b*Kstiff/10;        % backbone LJ energy parameter
lambda_b = lambda_r0;       % for soft core LJ potential
rc = 2*sig_r0;              % repulsive cutoff lengthscale
if TurnOnLJ==0
    eps_r0 = 0;
    eps_b = 0;
end
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# DEFINE REPULSION (IF ANY)'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
L = ['\nvariable    sig_r0  equal ',num2str(sig_r0,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    eps_r0  equal ',num2str(eps_r0,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    lambda_r0 equal ',num2str(lambda_r0,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    sig_b   equal ',num2str(2*Ldyn,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    eps_b   equal ',num2str(eps_b,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    lambda_b equal ',num2str(lambda_b,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    r_c     equal ',num2str(rc,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Define dynamic bond properties (ONLY USED IF DYNAMICS TURNED ON)
fsens = f0*3*kb*T/(b);      % Detachment force-sensitivity parameter
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# DYNAMIC BOND PARAMETERS'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Kinetic rates'; Assembled = [Assembled,L];
L = ['\nvariable    ka  equal ',num2str(ka0,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    kd  equal ',num2str(kd0,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Attachment cutoff distance'; Assembled = [Assembled,L];
L = ['\nvariable    Ldyn    equal ',num2str(Ldyn,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Force sensitivity'; Assembled = [Assembled,L];
L = ['\nvariable    fsens  equal ',num2str(fsens,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Define runtime scales and number of steps for each phase
taukd = 1/kd0;
tauka = 1/ka0;
tau_heal = 1/(kd0+ka0);
theal = 100*dt*dtFact;
tload = 0;
trelax = 0;
nheal = ceil(theal/dt);
nload = ceil(tload/dt);
nrelax = ceil(trelax/dt);
nmax = nheal+nload+nrelax;

N_extra_bonds = 2;
N_extra_special = 1000;
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# TIMESCALES & LOADING RATES'; Assembled = [Assembled,L];
L = '\n##########################`##########'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Timescales'; Assembled = [Assembled,L];
L = ['\nvariable    taukd   equal ',num2str(taukd,'%.1f')]; Assembled = [Assembled,L];
L = ['\nvariable    tauka   equal ',num2str(tauka,'%.1f')]; Assembled = [Assembled,L];
L = ['\nvariable    tau_heal equal ',num2str(tau_heal,'%.1f')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Timestep'; Assembled = [Assembled,L];
L = ['\nvariable    dt      equal ',num2str(dt,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Time for each portion of the simulation'; Assembled = [Assembled,L];
L = ['\nvariable    theal   equal ',num2str(theal,'%.1f')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Number of timesteps for each sequence'; Assembled = [Assembled,L];
L = ['\nvariable    nheal   equal ',num2str(nheal,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    nmax    equal ',num2str(nmax,'%.0f')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Running parameters
tdata = dt*dtFact*1/10;  %Basically measures with frequency matching 1/10x the diffusion timescale
ttot = nmax*dt;
NData = ceil(ttot/tdata);
iout = ceil(nmax/NData); 
ithermo = ceil(iout/2);
N_dyn = 1;  %check dynamics every N_dyn steps 
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# SIMULATION SPECIFICS'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Output frequencies'; Assembled = [Assembled,L];
L = ['\nvariable    iout    equal ',num2str(iout,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    ithermo equal ',num2str(ithermo,'%.0f')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# For inter-processor communication'; Assembled = [Assembled,L];
L = ['\nneighbor    ',num2str(Lmax,'%.2e'),' bin']; Assembled = [Assembled,L];
L = ['\ncomm_modify cutoff  ',num2str(Lmax,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Initialize system
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# INITIALIZE THE SYSTEM'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Define atom style'; Assembled = [Assembled,L];
L = '\natom_style   bond'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Define bond type'; Assembled = [Assembled,L];
L = '\nbond_style    harmonic'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% InputFolder1 = ['/mnt/c/Users/rjwag/Documents/Research/NSF/',CurrentFolder,'&'];
wls_folder = strrep(CurrentFolder,'C:/','/mnt/c/');
InputFolder1 = [wls_folder,'&'];
InputFolder2 = TopologyFileName;
L = '\n# Import the topology'; Assembled = [Assembled,L];
L = ['\nread_data ',InputFolder1]; Assembled = [Assembled,L];
L = ['\n',InputFolder2,'&']; Assembled = [Assembled,L];
L = ['\n extra/bond/per/atom ',num2str(N_extra_bonds,'%.0f'),...
    ' extra/special/per/atom ',num2str(N_extra_special,'%.0f')]; 
Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Define repulsive potentials between each atom'; Assembled = [Assembled,L];
L = '\npair_style  lj/cut/soft $(2.0) $(1.0) ${lambda_b}'; Assembled = [Assembled,L];
L = '\npair_coeff  1 1	${eps_b} ${sig_b} ${lambda_b}'; Assembled = [Assembled,L];

L = '\n'; Assembled = [Assembled,L];

L = '\n# Label the groups'; Assembled = [Assembled,L];
L = '\ngroup	stationarydynamic type 1'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
L = '\nlabelmap	atom	1 stckr'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Initial Equilibration
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# REACH STEADY STATE'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% NO DYNAMICS FOR THIS STUDY
if TurnOnDynamics==1
    L = ['\nfix		dynamic stationarydynamic bond/dynamic ',num2str(N_dyn,'%.0f'),...
        ' 1 1 ${ka} ${kd} ${Ldyn} maxbond 1'];% bell ${fsens}'];
    Assembled = [Assembled,L];
end

L = '\n# Compute bond information'; Assembled = [Assembled,L];
L = '\ncompute 	Pairs all property/local btype batom1 batom2'; Assembled = [Assembled,L];
L = '\ncompute 	Bonds all bond/local dx dy dz dist force'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Outputs files'; Assembled = [Assembled,L];
L = ['\ndump mydump all custom ${iout} ',OutputAtom_loc,' id x y z vx vy vz type mol']; 
Assembled = [Assembled,L];
L = ['\ndump bondsdump all local ${iout} ',OutputBond_loc,' index c_Pairs[*] c_Bonds[*]']; 
Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Outputs in command window'; Assembled = [Assembled,L];
L = '\nthermo_style    custom step atoms etotal ke pe epair ebond'; Assembled = [Assembled,L];
L = '\nthermo          ${ithermo}'; Assembled = [Assembled,L];
L = '\nthermo_modify	flush yes'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Set timestep and minimization conditions'; Assembled = [Assembled,L];
L = '\ntimestep ${dt}'; Assembled = [Assembled,L];

L = '\n# Compute virial stress on each atom and sum throughout system';
Assembled = [Assembled,L];
L = '\ncompute  sigma all stress/atom NULL bond'; Assembled = [Assembled,L];
L = '\ncompute  sig11 all reduce sum c_sigma[1]'; Assembled = [Assembled,L];
L = '\ncompute  sig22 all reduce sum c_sigma[2]'; Assembled = [Assembled,L];
L = '\ncompute  sig33 all reduce sum c_sigma[3]'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Output equilibration data'; Assembled = [Assembled,L];
L = '\nvariable	S1 equal c_sig11'; Assembled = [Assembled,L];
L = '\nvariable	S2 equal c_sig22'; Assembled = [Assembled,L];
L = '\nvariable	S3 equal c_sig33'; Assembled = [Assembled,L];
L = '\nvariable	teq equal ${dt}*elapsed'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Run initial equilibration'; Assembled = [Assembled,L];
L = '\nrun  ${nheal}'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

fid = fopen(InputFileName,'wt');
fprintf(fid,Assembled);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeMovies(N)

global dims Np Nt Ns

% Plot periodic and unwrapped domains
figure(4)
clf
PlotNetwork

% figure(5)
% clf
% PlotNetworkUnwrapped;

if dims==2
    if ~isfolder('Charts')
        mkdir('Charts')
    end
    figure(4)
    PlotName = ['Charts/Np_',num2str(Np),'.Nt_',num2str(Nt),...
        '.Nb_',num2str(N),'.z_',num2str(Ns),'.png'];
    saveas(gcf,PlotName);

    figure(5)
    PlotName = ['Charts/Unwrapped.Np_',num2str(Np),'.Nt_',num2str(Nt),...
        '.Nb_',num2str(N),'.z_',num2str(Ns),'.png'];
    saveas(gcf,PlotName);
end

if dims==3
    figure(4)
    Frames = 120;
    ang0 = 45;
    anglerange = (linspace(ang0,ang0+360,Frames))';
    for i=1:Frames
        ang = anglerange(i);
        set(gca,'CameraViewAngleMode','Manual')
        axis equal
        view(ang,30)
        F1 = getframe(gcf);
        %         writeVideo(v,F1)
        mov(i) = F1;
    end
    if ~isfolder('Gifs')
        mkdir('Gifs')
    end
    GifName = ['Gifs/Np_',num2str(Np),'.Nt_',num2str(Nt),...
        '.Nb_',num2str(N),'.z_',num2str(Ns),'.gif'];
    movie2gif(mov, GifName, 'DelayTime', 0,'LoopCount',Inf)
end

% figure(5)
% %         if ~isfolder('Movies')
% %             mkdir('Movies')
% %         end
% %     MovieName = ['Movies/Unwrapped.Np_',num2str(Np),'.Nt_',num2str(N),...
% %        '.Nb_',num2str(Nt),'.z_',num2str(Ns)];
% %     if MakeMovie==1
% %         v = VideoWriter(MovieName);
% %         v.FrameRate = 10;
% %         open(v)
% %     end
% ang0 = 45;
% anglerange = (linspace(ang0,ang0+360,Frames))';
% for i=1:Frames
%     ang = anglerange(i);
%     set(gca,'CameraViewAngleMode','Manual')
%     axis equal
%     view(ang,30)
%     F1 = getframe(gcf);
%     %         writeVideo(v,F1)
%     mov(i) = F1;
% end
% %     close(v)
% if ~isfolder('Gifs')
%     mkdir('Gifs')
% end
% GifName = ['Gifs/Unwrapped.Np_',num2str(Np),'.Nt_',num2str(Nt),...
%     '.Nb_',num2str(N),'.z_',num2str(Ns),'.gif'];
% movie2gif(mov, GifName, 'DelayTime', 0,'LoopCount',Inf)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveTopologyDataLAMMPS(N_Kuhn,kb,T,b,mass_mer)

global FileTag Positions dims Lx Ly Lz ConnectionsSP...
    N_atoms N_bonds TopologyFileName Ns 

% Define N_atoms and N_bonds
N_atoms = size(Positions,1);
Froms = zeros(N_atoms*(size(ConnectionsSP,2)-1),1);
Tos = Froms;
for i=1:size(ConnectionsSP,2)-1
    rng = ((N_atoms*i-N_atoms+1):(N_atoms*i))';
    Froms(rng) = ConnectionsSP(:,i);
    Tos(rng) = ConnectionsSP(:,end);
end
Pairs = [Froms Tos];
UniquePairs = sort(Pairs,2);
UniquePairs(UniquePairs(:,1)==0,:) = [];
UniquePairs = unique(UniquePairs,'rows');
N_bonds = size(UniquePairs,1);

% Formats output data as follows for LAMMPS:

% LAMMPS Description           (1st line of file)
% 
% 100 atoms         (this must be the 3rd line, 1st 2 lines are ignored)
% 95 bonds                (# of bonds to be simulated)
% 50 angles               (include these lines even if number = 0)
% 30 dihedrals
% 20 impropers
% 
% 5 atom types           (# of nonbond atom types)
% 10 bond types          (# of bond types = sets of bond coefficients)
% 18 angle types         
% 20 dihedral types      (do not include a bond,angle,dihedral,improper type
% 2 improper types             line if number of bonds,angles,etc is 0)
% 
% -0.5 0.5 xlo xhi       (for periodic systems this is box size,
% -0.5 0.5 ylo yhi        for non-periodic it is min/max extent of atoms)
% -0.5 0.5 zlo zhi       (do not include this line for 2-d simulations)
% 
% Masses
% 
%   1 mass
%   ...
%   N mass                           (N = # of atom types)
% 
% Nonbond Coeffs
% 
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of atom types)
% 
% Bond Coeffs
% 
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of bond types)
% 
% Angle Coeffs
% 
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of angle types)
% 
% Dihedral Coeffs
% 
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of dihedral types)
% 
% Improper Coeffs
% 
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of improper types)
% 
% BondBond Coeffs
% 
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of angle types)
% 
% BondAngle Coeffs
% 
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of angle types)
% 
% MiddleBondTorsion Coeffs
% 
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of dihedral types)
% 
% EndBondTorsion Coeffs
% 
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of dihedral types)
% 
% AngleTorsion Coeffs
% 
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of dihedral types)
% 
% AngleAngleTorsion Coeffs
% 
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of dihedral types)
% 
% BondBond13 Coeffs
% 
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of dihedral types)
% 
% AngleAngle Coeffs
% 
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of improper types)
% 
% Atoms
% 
%   1 molecule-tag atom-type q x y z nx ny nz  (nx,ny,nz are optional -
%   ...                                    see "true flag" input command)
%   ...                
%   N molecule-tag atom-type q x y z nx ny nz  (N = # of atoms)
% 
% Velocities
% 
%   1 vx vy vz
%   ...
%   ...                
%   N vx vy vz                        (N = # of atoms)
% 
% Bonds
% 
%   1 bond-type atom-1 atom-2
%   ...
%   N bond-type atom-1 atom-2         (N = # of bonds)
% 
% Angles
% 
%   1 angle-type atom-1 atom-2 atom-3  (atom-2 is the center atom in angle)
%   ...
%   N angle-type atom-1 atom-2 atom-3  (N = # of angles)
% 
% Dihedrals
% 
%   1 dihedral-type atom-1 atom-2 atom-3 atom-4  (atoms 2-3 form central bond)
%   ...
%   N dihedral-type atom-1 atom-2 atom-3 atom-4  (N = # of dihedrals)
% 
% Impropers
% 
%   1 improper-type atom-1 atom-2 atom-3 atom-4  (atom-2 is central atom)
%   ...
%   N improper-type atom-1 atom-2 atom-3 atom-4  (N = # of impropers)

% Initialize string
Assembled = [];
N_atom_types = 1;   % stable crosslinks, diffusing dynamic, stationary dynamic
N_bond_types = 1;   % Short, dynamics vs long implicit chains

% Header
L = ['Initial Topology',FileTag]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
L = ['\n',num2str(N_atoms),' atoms']; Assembled = [Assembled,L];
L = ['\n',num2str(N_bonds),' bonds']; Assembled = [Assembled,L];
L = ['\n',num2str(0),' angles']; Assembled = [Assembled,L];
L = ['\n',num2str(0),' dihedrals']; Assembled = [Assembled,L];
L = ['\n',num2str(0),' impropers']; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
L = ['\n',num2str(N_atom_types),' atom types']; Assembled = [Assembled,L];
L = ['\n',num2str(N_bond_types),' bond types']; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
L = ['\n',num2str(-Lx/2),' ',num2str(Lx/2),' xlo xhi']; Assembled = [Assembled,L];
L = ['\n',num2str(-Ly/2),' ',num2str(Ly/2),' ylo yhi']; Assembled = [Assembled,L];
L = ['\n',num2str(-Lz/2),' ',num2str(Lz/2),' zlo zhi']; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Append masses
weight = 0.5;
mass_backbone = N_Kuhn*mass_mer + 1/2*Ns*weight*N_Kuhn*mass_mer;
L = '\nMasses';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
for i=1:N_atom_types
    type = i;
    mass = mass_backbone;

    L = ['\n',num2str(type),' ',num2str(mass)];
    Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% Append bond coefficients - for Langevin spring potentials
L = '\nBond Coeffs';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
for i=1:N_bond_types
    type = i;
    coeff1 = 1/2*3*kb*T/b^2;
    coeff2 = b;
    L = ['\n',num2str(type),' ',...
        num2str(coeff1),' ',num2str(coeff2)];
    Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% Append coordinates
L = '\nAtoms';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

for i=1:N_atoms
    tag = Positions(i,dims+1);
    mol_tag = Positions(i,dims+2);
    x = Positions(i,1); y = Positions(i,2);
    if dims==2
        z = 0;
    else
        z = Positions(i,3);
    end
    atom_type = 1; %Stationary stable bonds

    L = ['\n',num2str(tag),' ',num2str(mol_tag),' ',num2str(atom_type),' ',...
        num2str(x),' ',num2str(y),' ',num2str(z)];
    Assembled = [Assembled,L];
end
L = '\n    '; Assembled = [Assembled,L];

% Append initial velocities
L = '\nVelocities';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
for i=1:N_atoms
    tag = Positions(i,dims+1);
    vx = 0; vy = 0; vz = 0;
    L = ['\n',num2str(tag),' ',num2str(vx),' ',num2str(vy),' ',num2str(vz)];
    Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% Append bonds
L = '\n'; Assembled = [Assembled,L];

% Save the topology data for read into LAMMPS
fid = fopen(TopologyFileName,'wt');
fprintf(fid,Assembled);

CheckFig=0;
if CheckFig==1
    PlotNetworkUnwrapped;
    for i=1:size(Positions,1)
        pos = Positions(i,1:dims);
        tag = Positions(i,dims+1);
        scatter3(pos(1),pos(2),pos(3),'g','filled')
        text(pos(1),pos(2),pos(3),['-',num2str(tag)],'FontSize',8)
    end
    xlim([-Inf Inf])
    ylim([-Inf Inf])
    zlim([-Inf Inf])
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveMatlabTopologyData

global FileTag TopologyDir...
    X_SP Y_SP Z_SP Rx_SP Ry_SP Rz_SP...TopologyDirs
    X_SP_Unwrapped Y_SP_Unwrapped Z_SP_Unwrapped ConnectionsSP...

FileName = [TopologyDir,'/X_SP',FileTag,'.txt'];
writematrix(X_SP,FileName);
FileName = [TopologyDir,'/Y_SP',FileTag,'.txt'];
writematrix(Y_SP,FileName);
FileName = [TopologyDir,'/Rx_SP',FileTag,'.txt'];
writematrix(Rx_SP,FileName);
FileName = [TopologyDir,'/Ry_SP',FileTag,'.txt'];
writematrix(Ry_SP,FileName);
FileName = [TopologyDir,'/ConnectionsSP',FileTag,'.txt'];
writematrix(ConnectionsSP,FileName);
FileName = [TopologyDir,'/Unwrapped.X_SP',FileTag,'.txt'];
writematrix(X_SP_Unwrapped,FileName);
FileName = [TopologyDir,'/Unwrapped.Y_SP',FileTag,'.txt'];
writematrix(Y_SP_Unwrapped,FileName);
if ~isempty(Z_SP)
    FileName = [TopologyDir,'/Z_SP',FileTag,'.txt'];
    writematrix(Z_SP,FileName);
    FileName = [TopologyDir,'/Rz_SP',FileTag,'.txt'];
    writematrix(Rz_SP,FileName);
    FileName = [TopologyDir,'/Unwrapped.Z_SP',FileTag,'.txt'];
    writematrix(Z_SP_Unwrapped,FileName);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SizeTheDomain(Np,N_Kuhn,b,Separation)

global Lx Ly Lz PairsPerEdge MinSep PairDistance NominalSeparation

PairsPerEdge = round(Np^(1/3));

MinSep = 2*N_Kuhn*b;
PairDistance = Separation;
NominalSeparation = MinSep+PairDistance;

Lx = NominalSeparation*PairsPerEdge;
Ly = Lx;
Lz = Lx;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitializeNodes

global N Lx Ly Lz Positions Corners dims

InitiateGrid; 

N = size(Positions,1);

% end 

Positions(:,dims+1) = (1:N)';    %Index No.
Positions(:,dims+2) = (1:N)';    %Molecule number

if dims==2
Corners = [-Lx/2 -Ly/2;
            Lx/2 -Ly/2;
            Lx/2  Ly/2;
           -Lx/2  Ly/2];
elseif dims==3
Corners = [-Lx/2 -Ly/2 -Lz/2;
            Lx/2 -Ly/2 -Lz/2;
            Lx/2  Ly/2 -Lz/2;
           -Lx/2  Ly/2 -Lz/2;
           -Lx/2 -Ly/2 Lz/2;
            Lx/2 -Ly/2 Lz/2;
            Lx/2  Ly/2 Lz/2;
           -Lx/2  Ly/2 Lz/2;];
end

CheckCorners = 1;
ShowFigDebug = 1;
if ShowFigDebug==1
    figure(1)
    clf
    hold on
    if dims==2
        scatter(Positions(:,1),Positions(:,2),'k','filled')
        if CheckCorners==1
            scatter(Corners(:,1),Corners(:,2),'r','filled')
        end
    else
        s = scatter3(Positions(:,1),Positions(:,2),Positions(:,3),'k','filled');
        s.SizeData = 6;
        if CheckCorners==1
            s = scatter3(Corners(:,1),Corners(:,2),Corners(:,3),'r','filled');
            s.SizeData = 6;
        end
        zlim([-Lz/2 Lz/2])
        pbaspect([1 1 1])
        view([30,30])
    end
    xlim([-Lx/2 Lx/2])
    ylim([-Ly/2 Ly/2])
    daspect([1 1 1])
    pbaspect([1 1 1])
    title(['N = ',num2str(size(Positions,1))])
    xlabel('x')
    ylabel('y')
    zlabel('z')
end
close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitiateGrid

global Positions Lx Ly Lz PairsPerEdge MinSep PairDistance

XSpacing = zeros(PairsPerEdge*2,1);
for i=1:PairsPerEdge*2
    if i==1
        XSpacing(i) = 0;
    elseif ~mod(i,2)
        XSpacing(i) = XSpacing(i-1)+PairDistance;
    else
        XSpacing(i) = XSpacing(i-1)+MinSep;
    end
end
YSpacing = (0:MinSep:MinSep*(PairsPerEdge-1))';
ZSpacing = (0:MinSep:MinSep*(PairsPerEdge-1))';

[X,Y,Z] = meshgrid(XSpacing(1:2:end),YSpacing,ZSpacing);
X = X(:);%-mean(X(:));
Y = Y(:);%-mean(Y(:));
Z = Z(:);%-mean(Z(:));
Positions_tethers = [X(:) Y(:) Z(:)];

[X,Y,Z] = meshgrid(XSpacing(2:2:end),YSpacing,ZSpacing);
X = X(:);%-mean(X(:));
Y = Y(:);%-mean(Y(:));
Z = Z(:);%-mean(Z(:));
Positions_bindingsites = [X(:) Y(:) Z(:)];

Positions = [Positions_tethers;Positions_bindingsites];
Positions_tethers = Positions_tethers - mean(Positions,1);
Positions_bindingsites = Positions_bindingsites - mean(Positions,1);
Positions = Positions - mean(Positions,1);

Lx = (max(X)+MinSep*5);
Ly = Lx;
Lz = Lx;

CheckFig = 1;
if CheckFig==1
    figure(1)
    clf; hold on
    s = scatter3(Positions_tethers(:,1),...
        Positions_tethers(:,2),...
        Positions_tethers(:,3),'c','filled');
    s.SizeData = 6;
    s = scatter3(Positions_bindingsites(:,1),...
        Positions_bindingsites(:,2),...
        Positions_bindingsites(:,3),'r','filled');
    s.SizeData = 6;
    xlim([-Lx/2 Lx/2])
    ylim([-Ly/2 Ly/2])
    zlim([-Lz/2 Lz/2])
    daspect([1 1 1])
    view(60,30)
end
close all
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotNetwork

global X_SP Y_SP Z_SP Rx_SP Ry_SP Rz_SP Lx Ly Lz ConnectionsSP Positions...
    Corners dims PlotChains N PlotStickers

hold on

%Plot Borders

%Plot Window
C1 = Corners(1,:);
C2 = Corners(2,:);
C3 = Corners(3,:);
C4 = Corners(4,:);
if dims==3
    C5 = Corners(5,:);
    C6 = Corners(6,:);
    C7 = Corners(7,:);
    C8 = Corners(8,:);
end
if dims==2
    plot([C1(1) C2(1)],[C1(2) C2(2)],'k--','LineWidth',1.5)
    plot([C2(1) C3(1)],[C2(2) C3(2)],'k--','LineWidth',1.5)
    plot([C3(1) C4(1)],[C3(2) C4(2)],'k--','LineWidth',1.5)
    plot([C4(1) C1(1)],[C4(2) C1(2)],'k--','LineWidth',1.5)
end
if dims==3
    plot3([C1(1) C2(1)],[C1(2) C2(2)],[C1(3) C2(3)],'k--','LineWidth',1.5)
    plot3([C2(1) C3(1)],[C2(2) C3(2)],[C2(3) C3(3)],'k--','LineWidth',1.5)
    plot3([C3(1) C4(1)],[C3(2) C4(2)],[C3(3) C4(3)],'k--','LineWidth',1.5)
    plot3([C4(1) C1(1)],[C4(2) C1(2)],[C4(3) C1(3)],'k--','LineWidth',1.5)
    plot3([C1(1) C5(1)],[C1(2) C5(2)],[C1(3) C5(3)],'k--','LineWidth',1.5)
    plot3([C2(1) C6(1)],[C2(2) C6(2)],[C2(3) C6(3)],'k--','LineWidth',1.5)
    plot3([C3(1) C7(1)],[C3(2) C7(2)],[C3(3) C7(3)],'k--','LineWidth',1.5)
    plot3([C4(1) C8(1)],[C4(2) C8(2)],[C4(3) C8(3)],'k--','LineWidth',1.5)
    plot3([C5(1) C6(1)],[C5(2) C6(2)],[C5(3) C6(3)],'k--','LineWidth',1.5)
    plot3([C6(1) C7(1)],[C6(2) C7(2)],[C6(3) C7(3)],'k--','LineWidth',1.5)
    plot3([C7(1) C8(1)],[C7(2) C8(2)],[C7(3) C8(3)],'k--','LineWidth',1.5)
    plot3([C8(1) C5(1)],[C8(2) C5(2)],[C8(3) C5(3)],'k--','LineWidth',1.5)
    view([30,45])
end    

if PlotChains==1
    Stickers = Positions(N+1:end,dims+1);
    Pairs = ConnectionsSP(1:N,1:end-1);
    Pairs = Pairs(:);
    X1 = X_SP(1:N,:); X1 = X1(:); X1(ismember(Pairs,Stickers)) = [];
    Rx1 = Rx_SP(1:N,:); Rx1 = Rx1(:); Rx1(ismember(Pairs,Stickers)) = [];
    X2 = X1 + Rx1;
    Y1 = Y_SP(1:N,:); Y1 = Y1(:); Y1(ismember(Pairs,Stickers)) = [];
    Ry1 = Ry_SP(1:N,:); Ry1 = Ry1(:); Ry1(ismember(Pairs,Stickers)) = [];
    Y2 = Y1 + Ry1;
    if dims==3
        Z1 = Z_SP(1:N,:); Z1 = Z1(:); Z1(ismember(Pairs,Stickers)) = [];
        Rz1 = Rz_SP(1:N,:); Rz1 = Rz1(:); Rz1(ismember(Pairs,Stickers)) = [];
        Z2 = Z1 + Rz1;
    end
    Pairs(ismember(Pairs,Stickers)) = [];
    X1(Pairs(:)==0) = [];
    X2(Pairs(:)==0) = [];
    Y1(Pairs(:)==0) = [];
    Y2(Pairs(:)==0) = [];
    if dims==3
        Z1(Pairs(:)==0) = [];
        Z2(Pairs(:)==0) = [];
    end

    Color = 'k';
    LineWidth = 0.5;
    if dims==2
        plot([X1,X2]',[Y1,Y2]',Color,'LineWidth',1)
    else
        plot3([X1,X2]',[Y1,Y2]',[Z1,Z2]',Color,'LineWidth',LineWidth)
    end
end

if PlotStickers==1
    Pairs = ConnectionsSP(N+1:end,1:end-1);
    Pairs = Pairs(:);
    X1 = X_SP(N+1:end,:); X1 = X1(:);  
    Rx1 = Rx_SP(N+1:end,:); Rx1 = Rx1(:); 
    X2 = X1 + Rx1;
    Y1 = Y_SP(N+1:end,:); Y1 = Y1(:); 
    Ry1 = Ry_SP(N+1:end,:); Ry1 = Ry1(:); 
    Y2 = Y1 + Ry1;
    if dims==3
        Z1 = Z_SP(N+1:end,:); Z1 = Z1(:);
        Rz1 = Rz_SP(N+1:end,:); Rz1 = Rz1(:); 
        Z2 = Z1 + Rz1;
    end
    X1(Pairs(:)==0) = [];
    X2(Pairs(:)==0) = [];
    Y1(Pairs(:)==0) = [];
    Y2(Pairs(:)==0) = [];
    if dims==3
        Z1(Pairs(:)==0) = [];
        Z2(Pairs(:)==0) = [];
    end

    Color = 'K';
    LineWidth = 0.5;
    if dims==2
        plot([X1,X2]',[Y1,Y2]','Color',Color,'LineWidth',1);
    else
        plot3([X1,X2]',[Y1,Y2]',[Z1,Z2]','Color',Color,'LineWidth',LineWidth);
    end
end

Fact = 1;
xlim([-Fact*Lx Fact*Lx])
ylim([-Fact*Ly Fact*Ly])
if dims==3
    ylim([-Fact*Lz Fact*Lz])
end

ShowScatter = 1;
if ShowScatter==1
    if dims==2
        s = scatter(Positions(1:N,1),Positions(1:N,2),'k','filled');
    else
        s = scatter3(Positions(1:N,1),Positions(1:N,2),Positions(1:N,3),'k','filled');
    end
    s.SizeData = 2;
end
if PlotStickers==1
    s = scatter3(Positions(N+1:end,1),Positions(N+1:end,2),Positions(N+1:end,3),'r','filled');
    s.SizeData = 2;
end

axis off

set(gcf,'Color','w')
daspect([1 1 1])
pbaspect([1 1 1])
if dims==3
    view([30,30])
end
box on

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotNetworkUnwrapped

global X_SP_Unwrapped Y_SP_Unwrapped Z_SP_Unwrapped Rx_SP Ry_SP Rz_SP...
    Lx Ly Lz ConnectionsSP PosUnwrapped...
    Corners dims PlotChains N PlotStickers

clf
hold on

%Plot Borders

%Plot Window
C1 = Corners(1,:);
C2 = Corners(2,:);
C3 = Corners(3,:);
C4 = Corners(4,:);
if dims==3
    C5 = Corners(5,:);
    C6 = Corners(6,:);
    C7 = Corners(7,:);
    C8 = Corners(8,:);
end
if dims==2
    plot([C1(1) C2(1)],[C1(2) C2(2)],'k--','LineWidth',1.5)
    plot([C2(1) C3(1)],[C2(2) C3(2)],'k--','LineWidth',1.5)
    plot([C3(1) C4(1)],[C3(2) C4(2)],'k--','LineWidth',1.5)
    plot([C4(1) C1(1)],[C4(2) C1(2)],'k--','LineWidth',1.5)
end
if dims==3
    plot3([C1(1) C2(1)],[C1(2) C2(2)],[C1(3) C2(3)],'k--','LineWidth',1.5)
    plot3([C2(1) C3(1)],[C2(2) C3(2)],[C2(3) C3(3)],'k--','LineWidth',1.5)
    plot3([C3(1) C4(1)],[C3(2) C4(2)],[C3(3) C4(3)],'k--','LineWidth',1.5)
    plot3([C4(1) C1(1)],[C4(2) C1(2)],[C4(3) C1(3)],'k--','LineWidth',1.5)
    plot3([C1(1) C5(1)],[C1(2) C5(2)],[C1(3) C5(3)],'k--','LineWidth',1.5)
    plot3([C2(1) C6(1)],[C2(2) C6(2)],[C2(3) C6(3)],'k--','LineWidth',1.5)
    plot3([C3(1) C7(1)],[C3(2) C7(2)],[C3(3) C7(3)],'k--','LineWidth',1.5)
    plot3([C4(1) C8(1)],[C4(2) C8(2)],[C4(3) C8(3)],'k--','LineWidth',1.5)
    plot3([C5(1) C6(1)],[C5(2) C6(2)],[C5(3) C6(3)],'k--','LineWidth',1.5)
    plot3([C6(1) C7(1)],[C6(2) C7(2)],[C6(3) C7(3)],'k--','LineWidth',1.5)
    plot3([C7(1) C8(1)],[C7(2) C8(2)],[C7(3) C8(3)],'k--','LineWidth',1.5)
    plot3([C8(1) C5(1)],[C8(2) C5(2)],[C8(3) C5(3)],'k--','LineWidth',1.5)
    view([30,45])
end    

%Plot molecular chains
% if PlotChains==1
%     XTemp = Positions;
%     XTemp = sortrows(XTemp,size(XTemp,2)-1);  %Sort xl's by which mol they belong to
% 
%     for mol=1:Np
%         
% 
%     end
% 
% end

if PlotChains==1
    Stickers = PosUnwrapped(N+1:end,dims+1);
    Pairs = ConnectionsSP(1:N,1:end-1);
    Pairs = Pairs(:);
    X1 = X_SP_Unwrapped(1:N,:); X1 = X1(:); X1(ismember(Pairs,Stickers)) = [];
    Rx1 = Rx_SP(1:N,:); Rx1 = Rx1(:); Rx1(ismember(Pairs,Stickers)) = [];
    X2 = X1 + Rx1;
    Y1 = Y_SP_Unwrapped(1:N,:); Y1 = Y1(:); Y1(ismember(Pairs,Stickers)) = [];
    Ry1 = Ry_SP(1:N,:); Ry1 = Ry1(:); Ry1(ismember(Pairs,Stickers)) = [];
    Y2 = Y1 + Ry1;
    if dims==3
        Z1 = Z_SP_Unwrapped(1:N,:); Z1 = Z1(:); Z1(ismember(Pairs,Stickers)) = [];
        Rz1 = Rz_SP(1:N,:); Rz1 = Rz1(:); Rz1(ismember(Pairs,Stickers)) = [];
        Z2 = Z1 + Rz1;
    end
    Pairs(ismember(Pairs,Stickers)) = [];
    X1(Pairs(:)==0) = [];
    X2(Pairs(:)==0) = [];
    Y1(Pairs(:)==0) = [];
    Y2(Pairs(:)==0) = [];
    if dims==3
        Z1(Pairs(:)==0) = [];
        Z2(Pairs(:)==0) = [];
    end

    Color = 'k';
    LineWidth = 0.5;
    if dims==2
        plot([X1,X2]',[Y1,Y2]',Color,'LineWidth',1)
    else
        plot3([X1,X2]',[Y1,Y2]',[Z1,Z2]',Color,'LineWidth',LineWidth)
    end
end

if PlotStickers==1
    Pairs = ConnectionsSP(N+1:end,1:end-1);
    Pairs = Pairs(:);
    X1 = X_SP_Unwrapped(N+1:end,:); X1 = X1(:);  
    Rx1 = Rx_SP(N+1:end,:); Rx1 = Rx1(:); 
    X2 = X1 + Rx1;
    Y1 = Y_SP_Unwrapped(N+1:end,:); Y1 = Y1(:); 
    Ry1 = Ry_SP(N+1:end,:); Ry1 = Ry1(:); 
    Y2 = Y1 + Ry1;
    if dims==3
        Z1 = Z_SP_Unwrapped(N+1:end,:); Z1 = Z1(:);
        Rz1 = Rz_SP(N+1:end,:); Rz1 = Rz1(:); 
        Z2 = Z1 + Rz1;
    end
    X1(Pairs(:)==0) = [];
    X2(Pairs(:)==0) = [];
    Y1(Pairs(:)==0) = [];
    Y2(Pairs(:)==0) = [];
    if dims==3
        Z1(Pairs(:)==0) = [];
        Z2(Pairs(:)==0) = [];
    end

    Color = 'r';
    LineWidth = 0.5;
    if dims==2
        plot([X1,X2]',[Y1,Y2]','Color',Color,'LineWidth',1);
    else
        plot3([X1,X2]',[Y1,Y2]',[Z1,Z2]','Color',Color,'LineWidth',LineWidth);
    end
end

Fact = 1;
xlim([-Fact*Lx Fact*Lx])
ylim([-Fact*Ly Fact*Ly])
if dims==3
    ylim([-Fact*Lz Fact*Lz])
end

ShowScatter = 1;
if ShowScatter==1
    if dims==2
        s = scatter(PosUnwrapped(1:N,1),PosUnwrapped(1:N,2),'k','filled');
    else
        s = scatter3(PosUnwrapped(1:N,1),PosUnwrapped(1:N,2),PosUnwrapped(1:N,3),'k','filled');
    end
    s.SizeData = 2;
end
% if PlotStickers==1
%     s = scatter3(Positions(N+1:end,1),Positions(N+1:end,2),Positions(N+1:end,3),'r','filled');
%     s.SizeData = 6;
% end

axis off

set(gcf,'Color','w')
daspect([1 1 1])
pbaspect([1 1 1])
if dims==3
    view([30,30])
end
box on

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function movie2gif(mov, gifFile, varargin)
% ==================
% Matlab movie to GIF Converter.
%
% Syntax: movie2gif(mov, gifFile, prop, value, ...)
% =================================================
% The list of properties is the same like for the command 'imwrite' for the
% file format gif:
%
% 'BackgroundColor' - A scalar integer. This value specifies which index in
%                     the colormap should be treated as the transparent
%                     color for the image and is used for certain disposal
%                     methods in animated GIFs. If X is uint8 or logical,
%                     then indexing starts at 0. If X is double, then
%                     indexing starts at 1.
%
% 'Comment' - A string or cell array of strings containing a comment to be
%             added to the image. For a cell array of strings, a carriage
%             return is added after each row.
%
% 'DelayTime' - A scalar value between 0 and 655 inclusive, that specifies
%               the delay in seconds before displaying the next image.
%
% 'DisposalMethod' - One of the following strings, which sets the disposal
%                    method of an animated GIF: 'leaveInPlace',
%                    'restoreBG', 'restorePrevious', or 'doNotSpecify'.
%
% 'LoopCount' - A finite integer between 0 and 65535 or the value Inf (the
%               default) which specifies the number of times to repeat the
%               animation. By default, the animation loops continuously.
%               For a value of 0, the animation will be played once. For a
%               value of 1, the animation will be played twice, etc.
%
% 'TransparentColor' - A scalar integer. This value specifies which index
%                      in the colormap should be treated as the transparent
%                      color for the image. If X is uint8 or logical, then
%                      indexing starts at 0. If X is double, then indexing
%                      starts at 1
%
% *************************************************************************
% Copyright 2007-2013 by Nicolae Cindea.
if (nargin < 2)
    error('Too few input arguments');
end
if (nargin == 2)
    frameNb = size(mov, 2);
    isFirst = true;
    h = waitbar(0, 'Generate GIF file...');
    for i = 1:frameNb
        waitbar((i-1)/frameNb, h);
        [RGB, ~] = frame2im(mov(i));
        if (exist('rgb2ind', 'file'))
            [IND, map] = rgb2ind(RGB,256);
        else
            [IND, map] = aRGB2IND(RGB);
        end
        if isFirst
            imwrite(IND, map, gifFile, 'gif');
            isFirst=false;
        else
            imwrite(IND, map, gifFile, 'gif', 'WriteMode', 'append');
        end
    end
    close(h);
end
if (nargin > 2)
    h = waitbar(0, 'Generate GIF file...');
    frameNb = size(mov, 2);
    isFirst = true;
    for i = 1:frameNb
        waitbar((i-1)/frameNb, h);
        [RGB, ~] = frame2im(mov(i));
        if (exist('rgb2ind', 'file'))
            [IND, map] = rgb2ind(RGB,256);
        else
            [IND, map] = aRGB2IND(RGB);
        end
        if isFirst
            args = varargin;
            imwrite(IND, map, gifFile, 'gif', args{:});
            isFirst=false;
            
            % supress 'LoopCount' option from the args!!
            args = varargin;
            l = length(args);
            
            posLoopCount = 0;
            for ii = 1:l
                if(ischar(args{ii}))
                    if strcmp(args{ii}, 'LoopCount')
                        posLoopCount = ii;
                    end
                end
            end
            if (posLoopCount)
                args = {args{1:posLoopCount-1}, args{posLoopCount+2:end}};
            end
            
        else
            imwrite(IND, map, gifFile, 'gif', 'WriteMode', 'append', ...
                args{:});
        end
    end
    close(h);
end
end
