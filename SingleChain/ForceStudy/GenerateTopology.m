% Branched newtork topolgoy generator for input to LAMMPS via read_data cmd
% Robert J. Wagner
% 2022-11-29

% This version generates the x-links first and then connects

function GenerateTopology(Package,OverrideTop,OverrideIn,...
    CF,OF,LC,DC,BSOM,DTF,BT)

global dims ShowFiguresDebug...
    OldOrNew PlotStickers MakeMovie...
    TopologyFileName BeadSpringOrMeso...
    InputFileName ConstantsFileName...
    LengthConversion DamperConversion CurrentFolder OutputFolder dtFact...
    Sample Np D N_Kuhn stiffness kbT b dt tau0 damp...
    Lx Ly Lz Corners Corners0 PlotChains BondType

%% Set Basic Toggles
OldOrNew = 1;           %Different looping options for two diff. connectivity algorithms
ShowFiguresDebug = 1;   %1 to plot intermediate figures for debugging
LengthConversion = LC;
DamperConversion = DC;
BeadSpringOrMeso = BSOM;
CurrentFolder = CF;
OutputFolder = OF;
dtFact = DTF;
BondType = BT;

%% Unpack Swept Input Parameters
Sample = Package(1);    %Sample index
Np = Package(2);        %Number of molecules
D = Package(3);         %Diffusion coefficient [m2/s]
N_Kuhn = Package(4);    %Number of Kuhn segments in chain
stiffness = Package(5); %Stiffness of single harmonic bond

kbT = Package(6);       %Thermal energy
b = Package(7);         %Kuhn length
dt = Package(8);        %timestep

tau0 = (b*LengthConversion)^2/D;
damp = kbT/D/DamperConversion;

%% Callout Input Script
% Initializes non-sweeping parameters
dims = 3;   %Dimensionality of system
InputScript(Sample,Np,D,N_Kuhn,stiffness,kbT,b,dt);

%% Make diriectories and set filenames
SetDirAndFileNames;

%% Convert all to correct units
Tab = table(Sample,Np,D,N_Kuhn,stiffness,kbT,b,dt);
writetable(Tab,ConstantsFileName);

%%%%%%%%%%%%%%%%%%%%%%%% CHANGE LOOPING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
if ~isfile(TopologyFileName) || OverrideTop==1

    %% Size the domain
    % Determines correct domain size based on no. xls and nominal density
    SizeTheDomain(Np,N_Kuhn,b);

    %% Initialize the Nodes
    InitializeNodes;

    %% Define the bond(s) and Kuhn segments (for Bead-spring only)
    ConnectTheEnds;

    %% Define window at start of deformation
    if dims==2
        Corners = [-Lx/2 -Ly/2;     %BL
            Lx/2 -Ly/2;     %BR
            Lx/2  Ly/2;     %TR
            -Lx/2  Ly/2];    %TL
    elseif dims==3
        Corners = [-Lx/2 -Ly/2 -Lz/2;   %C1
            +Lx/2 -Ly/2 -Lz/2;  %C2
            +Lx/2 +Ly/2 -Lz/2;  %C3
            -Lx/2 +Ly/2 -Lz/2;  %C4
            -Lx/2 -Ly/2 +Lz/2;  %C5
            +Lx/2 -Ly/2 +Lz/2;  %C6
            +Lx/2 +Ly/2 +Lz/2;  %C7
            -Lx/2 +Ly/2 +Lz/2]; %C8
    end
    Corners0 = Corners;

    %% Equilibrate the bead-spring model
    if BeadSpringOrMeso==0
        % Define Boundaries
        FindBoundaryNodes;  %ID the Boundary Nodes %% UPDATE 3D

        % Replicate bounds with node types
        ReplicateBounds;	%Replicate boundary Nodes accros periodic bounds

        % Generate repulsive connectivity
        ConnectREP;         %Defines Repulsion Connection Matrix

        % Calculate Forces, Step Nodes and Equilibrate Domain to Low Energy State
        EquilibrateTheDomain;   %Defines all forces, steps positions of nodes and
                                %repeats until residual forces diminish to threshold
    end

    ShowFigsDebug = 1;
    if ShowFigsDebug==1
        figure(1); clf; hold on
        PlotStickers=1;
        PlotChains = 1;
        if ShowFiguresDebug==1
            PlotNetwork
        end
    end

    %% Output Topological Data
    SaveTopologyDataLAMMPS(N_Kuhn,b,stiffness);

    MakeMovie=0;
    if MakeMovie==1
        MakeMovies(N_Kuhn);
    end
end

if ~isfile(InputFileName) || OverrideIn==1
    PrintInputFiles(N_Kuhn,b,kbT,dt,damp,CurrentFolder);
end

% toc
% close
clear
clc
clear global

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ConnectTheEnds

global BeadSpringOrMeso Positions N ConnectionsSP...
    X_SP Y_SP Z_SP Rx_SP Ry_SP Rz_SP N_Kuhn b...
    sigmaREP CutoffPartREP CutoffCWAForce SpringCoeff_eq kbT_eq

if BeadSpringOrMeso
    ConnectionsSP = zeros(N,3);
    ConnectionsSP(:,end) = (1:N)';
    ConnectionsSP(1,1) = 2;
    ConnectionsSP(2,1) = 1;
    Rx_SP = zeros(N,size(ConnectionsSP,2)-1);
    Ry_SP = zeros(N,size(ConnectionsSP,2)-1);
    Rz_SP = zeros(N,size(ConnectionsSP,2)-1);
    X_SP = zeros(N,size(ConnectionsSP,2));
    Y_SP = zeros(N,size(ConnectionsSP,2));
    Z_SP = zeros(N,size(ConnectionsSP,2));
    for i=1:size(ConnectionsSP,1)
        for j=1:size(ConnectionsSP,2)-1
            if ConnectionsSP(i,j)~=0
                p1 = i; p2 = ConnectionsSP(i,j);
                Rx_SP(i,j) = Positions(p1,1)-Positions(p2,1);
                Ry_SP(i,j) = Positions(p1,2)-Positions(p2,2);
                Rz_SP(i,j) = Positions(p1,3)-Positions(p2,3);
                X_SP(i,j) = Positions(p2,1);
                Y_SP(i,j) = Positions(p2,2);
                Z_SP(i,j) = Positions(p2,3);
            end
        end
    end
    [X_SP,Y_SP,Z_SP] = ReIndexXY(ConnectionsSP,X_SP,Y_SP,Z_SP,Rx_SP);
else
    add_rng = (N+1:N+N_Kuhn-1)';

    % Interpolate
    add_x = InterpolateChronologically(Positions(1,1),Positions(2,1),N_Kuhn+1);
     
    % break symmetry of beads in yz-plane
    max = b*sqrt(1-0.95);
    min = -max;%0;
    n = N_Kuhn-1;
    add_y = min+rand(1,n)*(max-min);
    add_z = min+rand(1,n)*(max-min);

    Positions(add_rng,1) = add_x(2:end-1)';
    Positions(add_rng,2) = add_y';
    Positions(add_rng,3) = add_z';

    % Re-order the Positions so 1 and end are the end points
    Positions(2,end-1) = add_rng(end);
    Positions(add_rng,end-1) = add_rng-1;
    Positions = sortrows(Positions,size(Positions,2)-1);
    Positions(:,end) = 1;

    N = size(Positions,1);

    % Connect
    ConnectionsSP = zeros(N,3); ConnectionsSP(:,end) = Positions(:,end-1);
    X_SP = zeros(size(ConnectionsSP)); X_SP(:,end) = Positions(:,1);
    Y_SP = zeros(size(ConnectionsSP)); Y_SP(:,end) = Positions(:,2);
    Z_SP = zeros(size(ConnectionsSP)); Z_SP(:,end) = Positions(:,3);
    Rx_SP = zeros(N,size(ConnectionsSP,2)-1);
    Ry_SP = zeros(N,size(ConnectionsSP,2)-1);
    Rz_SP = zeros(N,size(ConnectionsSP,2)-1);
    for i=1:N
        if i==1
            neighbs = 2;
        elseif i==N
            neighbs = N-1;
        else
            neighbs = [i-1,i+1];
        end
        for j=1:length(neighbs)
            p1 = ConnectionsSP(i,end); p2 = neighbs(j);
            pos1 = Positions(p1,1:3); pos2 = Positions(p2,1:3);
            col = find(ConnectionsSP(p1,:)==0,1,'first');
            if ~isempty(col)
                ConnectionsSP(p1,col) = p2;
                X_SP(p1,col) = pos2(1);
                Y_SP(p1,col) = pos2(2);
                Z_SP(p1,col) = pos2(3);
                Rx_SP(p1,col) = pos2(1)-pos1(1);
                Ry_SP(p1,col) = pos2(2)-pos1(2);
                Rz_SP(p1,col) = pos2(3)-pos1(3);
            end
        end
    end
    [X_SP,Y_SP,Z_SP] = ReIndexXY(ConnectionsSP,X_SP,Y_SP,Z_SP,Rx_SP);

    % Adjust force parameters for re-equilibration
    sigmaREP = b; CutoffPartREP = b; CutoffCWAForce = 0.5*b;
    L_eq = b; fact = 0.25;
    SpringCoeff_eq = fact*3*kbT_eq/(L_eq*b);
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = InterpolateChronologically(p1,p2,npts)

if p1>p2
    out = linspace(p2,p1,npts);
    out = fliplr(out);
else
    out = linspace(p1,p2,npts);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PrintInputFiles(N_Kuhn,b,kbT,dt,damp,CurrentFolder)

global TopologyFileName OutputAtom_loc OutputBond_loc InputFileName...
    BeadSpringOrMeso BondType

TurnOnLJ = 0;       % 0 for no volume exclusion (do not need repulsive interactions
if BeadSpringOrMeso==0 % Bead-spring
    damp_eff = damp;
elseif BeadSpringOrMeso==1 % Mesoscale
    damp_eff = damp*N_Kuhn^(2/3);    % Effective damper needed to roughly mimic bead-spring diffusion
elseif BeadSpringOrMeso==2
    damp_eff = damp*N^(2/3);    % Effective damper needed to roughly mimic bead-spring diffusion
end

% Initialize string
Assembled = [];

% Define units
L = '####################################'; Assembled = [Assembled,L];
L = '\n# DEFINE UNITS'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
L = '\nunits        lj'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Define lengthscales
Nb = N_Kuhn*b;           %chain length
Nb2 = Nb*b;             %N*b^2
r0 = sqrt(N_Kuhn)*b;     %mean random walk chain length
T = 1;                  % Temperature - 293 K = 1 in-model unit of temp.
Lmax = 0.4*Nb;          %max LJ cutoff
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# DEFINE LENGTH & FORCE SCALES'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Defining length scales'; Assembled = [Assembled,L];
L = ['\nvariable    b       equal ',num2str(b,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    N       equal ',num2str(N_Kuhn,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    Nb      equal ',num2str(Nb,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    cut_dyn equal ',num2str(0.95*Nb,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    Nb2     equal ',num2str(Nb2,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    r0      equal ',num2str(r0,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    T       equal ',num2str(T,'%.1f')]; Assembled = [Assembled,L];
L = ['\nvariable    Lmax    equal ',num2str(Lmax,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Define Thermodynamic properties
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# DEFINE LENGTH & FORCE SCALES'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
L = ['\nvariable    damp    equal ',num2str(damp_eff,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Define repulsive properties, if any properties
sig_b = b/(2^(1/6));                % backbone LJ length scale
eps_b = kbT;
lambda_b = 0.5;       % for soft core LJ potential
rc = 2*sig_b;              % repulsive cutoff lengthscale
if TurnOnLJ==0
    eps_b = 0;
end
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# DEFINE REPULSION (IF ANY)'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
L = ['\nvariable    sig_b   equal ',num2str(sig_b,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    eps_b   equal ',num2str(eps_b,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    lambda_b equal ',num2str(lambda_b,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    r_c     equal ',num2str(rc,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Define runtime scales and number of steps for each phase
npts = round(Nb/b*1.25);
distances = (linspace(b,0.95*Nb,npts))';
distances = flipud(distances);
incr_fact = 10;
nhold = 8e3*incr_fact;    %number of hold steps at each end-to-end distances
nload = 10;     %number of loading steps
thold = nhold*dt;
tload = nload*dt;
no_holds = length(distances);
no_loads = no_holds-1;
nmax = nhold*no_holds+nload*no_loads;

N_extra_bonds = 2;
N_extra_special = 1000;
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# TIMESCALES & LOADING RATES'; Assembled = [Assembled,L];
L = '\n##########################`##########'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Timestep'; Assembled = [Assembled,L];
L = ['\nvariable    dt      equal ',num2str(dt,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Time for each portion of the simulation'; Assembled = [Assembled,L];
L = ['\nvariable    tload   equal ',num2str(tload,'%.1f')]; Assembled = [Assembled,L];
L = ['\nvariable    thold   equal ',num2str(thold,'%.1f')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Number of timesteps for each sequence'; Assembled = [Assembled,L];
L = ['\nvariable    nhold   equal ',num2str(nhold,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    nload   equal ',num2str(nload,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    nmax    equal ',num2str(nmax,'%.0f')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Running parameters
iout = 20*incr_fact; 
ithermo = 100*iout;
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# SIMULATION SPECIFICS'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Output frequencies'; Assembled = [Assembled,L];
L = ['\nvariable    iout    equal ',num2str(iout,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    ithermo equal ',num2str(ithermo,'%.0f')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# For inter-processor communication'; Assembled = [Assembled,L];
if BeadSpringOrMeso==0 || BeadSpringOrMeso==1
    L = '\nneighbor    ${Lmax}     bin'; Assembled = [Assembled,L];
    L = '\ncomm_modify cutoff      ${Lmax}'; Assembled = [Assembled,L];
elseif BeadSpringOrMeso==2
    L = '\nneighbor    ${Nb}    bin'; Assembled = [Assembled,L];
    L = '\ncomm_modify cutoff   ${Nb}'; Assembled = [Assembled,L];
end
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
if BeadSpringOrMeso==0
	if BondType==0
		L = '\nbond_style   harmonic'; Assembled = [Assembled,L];
	else
		L = '\nbond_style   nonlinear'; Assembled = [Assembled,L];
	end
elseif BeadSpringOrMeso==1
	L = '\nbond_style   pade'; Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

wls_folder = strrep(CurrentFolder,'C:/','/mnt/c/');
% InputFolder1 = ['/mnt/c/Users/rjwag/Documents/Research/NSF/',CurrentFolder,'&'];
InputFolder1 = [wls_folder,'&'];
InputFolder2 = ['/',TopologyFileName];
L = '\n# Import the topology'; Assembled = [Assembled,L];
L = ['\nread_data ',InputFolder1]; Assembled = [Assembled,L];
L = ['\n',InputFolder2,'&']; Assembled = [Assembled,L];
L = ['\n extra/bond/per/atom ',num2str(N_extra_bonds,'%.0f'),...
    ' extra/special/per/atom ',num2str(N_extra_special,'%.0f')]; 
Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Define repulsive potentials between each atom'; Assembled = [Assembled,L];
L = '\npair_style  lj/cut/soft $(2.0) $(1.0) $(v_r_c)'; Assembled = [Assembled,L];
L = '\npair_coeff  1 1	${eps_b} ${sig_b} ${lambda_b}'; Assembled = [Assembled,L];
% if BeadSpringOrMeso==0 
    L = '\npair_coeff  2 2	${eps_b} ${sig_b} ${lambda_b}'; Assembled = [Assembled,L];
% end
L = '\n'; Assembled = [Assembled,L];

L = '\n# Label the groups'; Assembled = [Assembled,L];
rand_seed = randi(9999);
if BeadSpringOrMeso==1 % Mesoscale
    L = '\ngroup	ends    type 1'; Assembled = [Assembled,L];
    L = '\ngroup	ints    type 2'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
elseif BeadSpringOrMeso==0 % Bead-spring
    L = '\ngroup	ends    type 1'; Assembled = [Assembled,L];
    L = '\ngroup	ints    type 2'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
    L = '\nlabelmap	atom	1 bckbn		2 stckr'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    L = '\n# Set Brownian integration'; Assembled = [Assembled,L];
    L = ['\nfix      free    ints    brownian ${T} ',num2str(rand_seed),...
        ' gamma_t ${damp}']; Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% Initial Equilibration
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# REACH STEADY STATE'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];


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

L = '\n# Run the model at various distancs and output force'; Assembled = [Assembled,L];
if BeadSpringOrMeso==0
    end_bead = N_Kuhn+1;
elseif BeadSpringOrMeso==1
    end_bead = 2;
end

for i=1:length(distances)
    new_x_position = distances(i);
    L = ['\nset  atom   ',num2str(end_bead),...     % Set the position of the end atom progressively closer
        '   x   ',num2str(new_x_position,'%.4f'),...
        '   y   ',num2str(0,'%.1f'),...
        '   z   ',num2str(0,'%.1f')]; Assembled = [Assembled,L];
    L = ['\nrun  ',num2str(round(nhold/10))]; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];  
end

if N_Kuhn<=12
    neq = 4*nhold;
elseif N_Kuhn>12 && N_Kuhn<=18
    neq = 8*nhold;
elseif N_Kuhn>18 && N_Kuhn<=36
    neq = 17*nhold;
end

L = ['\nrun  ',num2str(neq)]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

distances = flipud(distances);
for i=1:length(distances)
    new_x_position = distances(i);
    L = ['\nset  atom   ',num2str(end_bead),...     % Set the position of the end atom progressively closer
        '   x   ',num2str(new_x_position,'%.4f'),...
        '   y   ',num2str(0,'%.1f'),...
        '   z   ',num2str(0,'%.1f')]; Assembled = [Assembled,L];
    L = '\nrun  ${nhold}'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];  
end

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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveTopologyDataLAMMPS(N_Kuhn,b,stiffness)

global FileTag Positions dims Lx Ly Lz ConnectionsSP...
    N_atoms N_bonds TopologyFileName BeadSpringOrMeso BondType

mass_mer = 1;  %Arbitrary
kbT_norm = 1;

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
if BeadSpringOrMeso==1 	% Mesoscale
    N_atom_types = 2;  	% Perms and xls
    N_bond_types = 1;   % Short, dynamics vs long implicit chains
elseif BeadSpringOrMeso==0% Bead-spring
    N_atom_types = 2;   % Perms and xls
    N_bond_types = 1;   % Kuhn segments and dynamic bonds
end

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
L = '\nMasses';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
for i=1:N_atom_types
    type = i;
    mass = mass_mer;
    L = ['\n',num2str(type),' ',num2str(mass)];
    Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% Append bond coefficients - for Langevin spring potentials
L = '\nBond Coeffs';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
for i=1:N_bond_types
    type = i;
    if BeadSpringOrMeso==1 % Mesoscale
%        coeff1 = 3/2*kbT_norm/(N_Kuhn*b^2);
%        coeff2 = 0;
        coeff1 = b;
        coeff2 = N_Kuhn;
        L = ['\n',num2str(type),' ',...
            num2str(coeff1),' ',num2str(coeff2)];
    elseif BeadSpringOrMeso==0 % Bead-spring
        if BondType==0
            coeff1 = stiffness*kbT_norm/b^2;
            coeff2 = b;
            L = ['\n',num2str(type),' ',...
                num2str(coeff1),' ',num2str(coeff2)];
        else
            coeff1 = stiffness*kbT_norm;
            coeff2 = b;
            coeff3 = b;
            L = ['\n',num2str(type),' ',...
                num2str(coeff1),' ',num2str(coeff2),' ',num2str(coeff3)];
        end
    end
    Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% Append coordinates
L = '\nAtoms';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
if BeadSpringOrMeso==1 % Mesoscale
    for i=1:N_atoms
        tag = Positions(i,dims+1);
        mol_tag = Positions(i,dims+2);
        x = Positions(i,1); y = Positions(i,2);
        if dims==2
            z = 0;
        else
            z = Positions(i,3);
        end
        atom_type = 1;

        L = ['\n',num2str(tag),' ',num2str(mol_tag),' ',num2str(atom_type),' ',...
            num2str(x),' ',num2str(y),' ',num2str(z)];
        Assembled = [Assembled,L];
    end
elseif BeadSpringOrMeso==0 % Bead-spring
    fixed = [1;N_atoms];
    diffusing = (2:N_atoms-1)';
    for i=1:N_atoms
        tag = Positions(i,dims+1);
        mol_tag = Positions(i,dims+2);
        x = Positions(i,1); y = Positions(i,2);
        if dims==2
            z = 0;
        else
            z = Positions(i,3);
        end

        if ismember(i,fixed)
            atom_type = 1;
        elseif ismember(i,diffusing)
            atom_type = 2;
        else
            warning('Atom does not belong to either type')
        end
        L = ['\n',num2str(tag),' ',num2str(mol_tag),' ',num2str(atom_type),' ',...
            num2str(x),' ',num2str(y),' ',num2str(z)];
        Assembled = [Assembled,L];
    end
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
if BeadSpringOrMeso==0 || BeadSpringOrMeso==1
    L = '\nBonds';  Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
    for i=1:N_bonds
        tag = i;
        atom_1 = UniquePairs(i,1); atom_2 = UniquePairs(i,2);
        bond_type = 1;
        L = ['\n',num2str(tag),' ',num2str(bond_type),' ',...
            num2str(atom_1),' ',num2str(atom_2)];
        Assembled = [Assembled,L];
    end
    L = '\n'; Assembled = [Assembled,L];
end

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
function EquilibrateTheDomain

global ConnectionsSP ConnectionsREP X_SP X_REP Positions Positions0...
    Y_SP Y_REP Z_SP Z_REP Rx_SP Rx_REP Ry_SP Ry_REP Rz_SP Rz_REP...
    ResMax ResAvg NumericalDamper RepulsionRedefFreq dt dims  

F_NET = CalculateNetForces;
% Do not move the first and last nodes - these are fixed
F_NET(1,:) = [0 0 0];
F_NET(end,:) = [0 0 0];

Damper = 10*NumericalDamper;
dX = Damper*F_NET;

Residuals = vecnorm(F_NET,2,2);

%Calculate Residuals
MaxResidual = max(max(Residuals));
AvgResidual = mean(mean(Residuals));
PrevRes = AvgResidual;

k=0;
wb = waitbar(0,'Equilibrating the domain...');
while MaxResidual>ResMax && AvgResidual>ResAvg
    k = k+1;
    
    Positions0 = Positions;
    Positions(:,1:dims) = Positions(:,1:dims) + dX; %Update postions
    
    %Update Positions & Separation Vectors in Connection Matrices
    if ~isempty(X_SP)
        [X_SP,Y_SP,Z_SP,Rx_SP,Ry_SP,Rz_SP] =...
            UpdateNodePositions(ConnectionsSP,X_SP,Y_SP,Z_SP,Rx_SP,Ry_SP,Rz_SP);
    end
    if ~isempty(X_REP)
        [X_REP,Y_REP,Z_REP,Rx_REP,Ry_REP,Rz_REP] =...
            UpdateNodePositions(ConnectionsREP,X_REP,Y_REP,Z_REP,...
            Rx_REP,Ry_REP,Rz_REP);
    end
    %Updates Positions, X_SP, Y_SP, X_REP, Y_REP,...
    %Rx_SP, Ry_SP, Rx_REP and Ry_REP
    
    
    %Implement Eularian condition and edefine repulsive neighbors
    if ~mod(k,RepulsionRedefFreq)
        % Transport Nodes that exit frame to opossite side (only acts if BCs
        % are for Couette flow right now (i.e., parallel plate rheometer)
        ImplementEularianBoundaryConditions
        FindBoundaryNodes
        ReplicateBounds
        ConnectREP
    end
        
    F_NET = CalculateNetForces;
    % Do not move the first and last nodes - these are fixed
    F_NET(1,:) = [0 0 0];
    F_NET(end,:) = [0 0 0];

    dX = Damper*F_NET;
       
    Residuals = vecnorm(F_NET,2,2);
    
    %Calculate Residuals
    MaxResidual = max(max(Residuals));
    AvgResidual = mean(mean(Residuals));
       
    if mod(k,100)==0
        disp(['Residuals = ',num2str(AvgResidual)]);
        waitbar(ResAvg/AvgResidual,wb,'Equilibrating the domain...');

%         PlotStickers = 1;
%         figure(2)
%         clf
%         PlotNetwork
        
        if PrevRes<=AvgResidual %Feedback loop to increase damper as needed (5% at a time)
            Damper = Damper*0.5;
            visc = dt/Damper;
            disp(['Increase visc to ',num2str(visc)]);
        end
        PrevRes = AvgResidual;
    end
    if mod(k,5000)==0
        break
    end
end
close(wb)

end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImplementEularianBoundaryConditions    

%Reinserts nodes across boundaries if they step outside the frame
%ONLY INSTALLED FOR COUETTE FLOW RIGHT NOW

global Positions dims...
    ConnectionsSP X_SP Y_SP Z_SP Rx_SP...
    ConnectionsREP X_REP Y_REP Z_REP Rx_REP Corners 

if dims==2
    CBL = Corners(1,:); CBR = Corners(2,:);
    CTR = Corners(3,:); CTL = Corners(4,:);
    Lx = CBR(1)-CBL(1); Ly = CTR(2)-CBR(2);
    X = Positions(:,1); Y = Positions(:,2);
    dYdX = CTR(2)-CTL(2); dXdY = CTR(1)-CBR(1);

    %Find particles out of left bound and reposition
    XofY_L = FindXofY(Y,CTL,CBL);
    if ~isempty(Positions(Positions(:,1)<XofY_L,:))
        Positions(Positions(:,1)<XofY_L,2) = Positions(Positions(:,1)<XofY_L,2) + dYdX;
        Positions(Positions(:,1)<XofY_L,1) = Positions(Positions(:,1)<XofY_L,1) + Lx;
    end

    %Find particles out of right bound and reposition
    XofY_R = FindXofY(Y,CTR,CBR);
    if ~isempty(Positions(Positions(:,1)>XofY_R,:))
        Positions(Positions(:,1)>XofY_R,2) = Positions(Positions(:,1)>XofY_R,2) - dYdX;
        Positions(Positions(:,1)>XofY_R,1) = Positions(Positions(:,1)>XofY_R,1) - Lx;
    end

    %Find particles out of top bound and reposition
    YofX_B = FindYofX(X,CBL,CBR);
    if ~isempty(Positions(Positions(:,2)<YofX_B,:))
        Positions(Positions(:,2)<YofX_B,1) = Positions(Positions(:,2)<YofX_B,1) + dXdY;
        Positions(Positions(:,2)<YofX_B,2) = Positions(Positions(:,2)<YofX_B,2) + Ly;
    end

    %Find particles out of top bound and reposition
    YofX_T = FindYofX(X,CTL,CTR);
    if ~isempty(Positions(Positions(:,2)>YofX_T,:))
        Positions(Positions(:,2)>YofX_T,1) = Positions(Positions(:,2)>YofX_T,1) - dXdY;
        Positions(Positions(:,2)>YofX_T,2) = Positions(Positions(:,2)>YofX_T,2) - Ly;
    end

elseif dims==3
    C1 = Corners(1,:);
    C2 = Corners(2,:);
    C4 = Corners(4,:);
    C5 = Corners(5,:);
    Lx = norm(C2-C1);
    Ly = norm(C1-C4);
    Lz = norm(C5-C1);
%     X = Positions(:,1); Y = Positions(:,2); Z = Positions(:,3);
    dXdY = C4(1)-C1(1);
    dYdX = C2(2)-C1(2);
    dYdZ = C5(2)-C1(2);
    dZdY = C4(3)-C1(3);
    dZdX = C2(3)-C1(3);
    dXdZ = C5(1)-C1(1);

    %Find particles out of left bound and reposition
%     v1 = C4-C1;
%     v2 = C5-C1;
%     C = C1;
%     Out = FindOutOfBounds(v1,v2,C);

    %Find particles out of left bound and reposition
    v1 = C4-C1;
    v2 = C5-C1;
    C = C1;
    Out = FindOutOfBounds(v1,v2,C);
    Positions(Out,1) = Positions(Out,1) + Lx;
    Positions(Out,2) = Positions(Out,2) + dYdX;
    Positions(Out,3) = Positions(Out,3) + dZdX;

    %Find particles out of right bound and reposition
    v1 = C5-C1;
    v2 = C4-C1;
    C = C2;
    Out = FindOutOfBounds(v1,v2,C);
    Positions(Out,1) = Positions(Out,1) - Lx;
    Positions(Out,2) = Positions(Out,2) - dYdX;
    Positions(Out,3) = Positions(Out,3) - dZdX;

    %Find particles out of front bound and reposition
    v1 = C5-C1;
    v2 = C2-C1;
    C = C1;
    Out = FindOutOfBounds(v1,v2,C);
    Positions(Out,1) = Positions(Out,1) + dXdY;
    Positions(Out,2) = Positions(Out,2) + Ly;
    Positions(Out,3) = Positions(Out,3) + dZdY;

    %Find particles out of front bound and reposition
    v1 = C2-C1;
    v2 = C5-C1;
    C = C4;
    Out = FindOutOfBounds(v1,v2,C);
    Positions(Out,1) = Positions(Out,1) - dXdY;
    Positions(Out,2) = Positions(Out,2) - Ly;
    Positions(Out,3) = Positions(Out,3) - dZdY;

    %Find particles out of front bound and reposition
    v1 = C2-C1;
    v2 = C4-C1;
    C = C1;
    Out = FindOutOfBounds(v1,v2,C);
    Positions(Out,1) = Positions(Out,1) + dXdZ;
    Positions(Out,2) = Positions(Out,2) + dYdZ;
    Positions(Out,3) = Positions(Out,3) + Lz;
    
    %Find particles out of front bound and reposition
    v1 = C4-C1;
    v2 = C2-C1;
    C = C5;
    Out = FindOutOfBounds(v1,v2,C);
    Positions(Out,1) = Positions(Out,1) - dXdZ;
    Positions(Out,2) = Positions(Out,2) - dYdZ;
    Positions(Out,3) = Positions(Out,3) - Lz;
end
if ~isempty(ConnectionsSP)
    [X_SP,Y_SP,Z_SP] = ReIndexXY(ConnectionsSP,X_SP,Y_SP,Z_SP,Rx_SP);
end
if ~isempty(ConnectionsREP)
    [X_REP,Y_REP,Z_REP] = ReIndexXY(ConnectionsREP,X_REP,Y_REP,Z_REP,Rx_REP);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Out = FindOutOfBounds(v1,v2,C)

global Positions dims

n_left = cross(v1,v2)/norm(cross(v1,v2));
n_arr = [n_left(1)*ones(size(Positions,1),1),...
    n_left(2)*ones(size(Positions,1),1),...
    n_left(3)*ones(size(Positions,1),1)];
C1_arr = zeros(size(Positions(:,1:dims)));
C1_arr(:,1) = C(1); C1_arr(:,2) = C(2); C1_arr(:,3) = C(3);
r = Positions(:,1:dims)-C1_arr;
temp = dot(n_arr,r,2);
Out = find(temp<0);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function YofX = FindYofX(X,C2,C1)
    m = (C2(2)-C1(2))/(C2(1)-C1(1)); %Left slope
    YofX = C2(2)-m*(C2(1)-X);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XofY = FindXofY(Y,C2,C1)
    m = (C2(2)-C1(2))/(C2(1)-C1(1)); %Left slope
    XofY = C2(1)-(C2(2)-Y)/m;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,Z,Rx,Ry,Rz] = UpdateNodePositions(Connections,X,Y,Z,Rx,Ry,Rz)

%Updates Positions, X_SP, Y_SP, X_REP, Y_REP,...
%Rx_SP, Ry_SP, Rx_REP and Ry_REP

%No dynamics involved here

global Positions Positions0 dims

dX = Positions(:,1:dims)-Positions0(:,1:dims);

X(Connections(:,1:end-1)==0) = NaN;
Y(Connections(:,1:end-1)==0) = NaN;
Z(Connections(:,1:end-1)==0) = NaN;
Connections(Connections(:,1:end-1)==0) = NaN;

[X,Y,Z] = ReIndexXY(Connections,X,Y,Z,Rx);

%Define displacement of heads of Rx and Ry vectors
dx_Head_temp = dX(:,1);
dy_Head_temp = dX(:,2);
if dims==3
    dz_Head_temp = dX(:,3);
end
% dx_Head_temp = X(:,end)-X_0(:,end);
% dy_Head_temp = Y(:,end)-Y_0(:,end);

%Store dX for heads and tails in matrices the same size as Rx and Ry
dx_Tail = zeros(size(Rx));
dy_Tail = zeros(size(Rx));
dz_Tail = zeros(size(Rx));
dx_Head = zeros(size(Rx));
dy_Head = zeros(size(Rx));
dz_Head = zeros(size(Rx));
for j=1:size(Rx,2)
    TempNodes = Connections(:,j);
    Rows = find(~isnan(TempNodes));
    
    %Allocate dx_Head and dx_Tail seperately
    %Store change in head positions
    dx_Head(Rows,j) = dx_Head_temp(Rows);
    dy_Head(Rows,j) = dy_Head_temp(Rows);
    if dims==3
        dz_Head(Rows,j) = dz_Head_temp(Rows);
    end
    
    %Define tail displacements;
    dx_Tail(Rows,j) = dX(TempNodes(Rows),1);
    dy_Tail(Rows,j) = dX(TempNodes(Rows),2);  
    if dims==3
        dz_Tail(Rows,j) = dX(TempNodes(Rows),3);  
    end
end

Rx = Rx-dx_Tail+dx_Head;
Ry = Ry-dy_Tail+dy_Head;
if dims==3
    Rz = Rz-dz_Tail+dz_Head;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,Z] = ReIndexXY(Connections,X,Y,Z,Rx)

global Positions dims

X(Connections(:,1:end-1)==0) = NaN;
Y(Connections(:,1:end-1)==0) = NaN;
Z(Connections(:,1:end-1)==0) = NaN;
Connections(Connections(:,1:end-1)==0) = NaN;

%Define New Positions in X and Y matrices
X(:,end) = Positions(:,1);
Y(:,end) = Positions(:,2);
if dims==3
    Z(:,end) = Positions(:,3);
end
for j=1:size(Rx,2)
    TempNodes = Connections(:,j);
    Rows = find(~isnan(TempNodes));
    
    %Update X and Y of Tails
    X(Rows,j) = Positions(TempNodes(Rows),1);
    Y(Rows,j) = Positions(TempNodes(Rows),2);  
    if dims==3
        Z(Rows,j) = Positions(TempNodes(Rows),3);
    end
end
        
end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_NET = CalculateNetForces

%Calculates forces based on Rx and Ry
%ADD PRESSURE CONSERVATION HERE - NET FORCE ON EACH NODE JUST NEEDS TO BE
%STORED AS [X,Y] AND ADDED TO F_NET

global N dims

%Spring Forces
F_SP = CalculateSpringForces;
if isempty(F_SP)
    F_SP = zeros(N,dims);
end

%Particle-to-particle forces
F_REP = CalculateRepulsiveForces;

%Wall forces (if there are top and bottom plates
F_WALL = CalculateWallForces;

F_NET = F_SP+F_REP+F_WALL;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_REP = CalculateRepulsiveForces

global Rx_REP Ry_REP Rz_REP F_REP_X F_REP_Y F_REP_Z...
    CutoffCWAForce r0 dims

if dims==2
    Norms = (Rx_REP.^2 + Ry_REP.^2).^0.5;
elseif dims==3
    Norms = (Rx_REP.^2 + Ry_REP.^2 + Rz_REP.^2).^0.5;
end

F_MAG = CalculateInvToAlphaForce(Norms);

rCutoff = CutoffCWAForce*r0;
Force2 = CalculateInvToAlphaForce(rCutoff+0.0001);
Force1 = CalculateInvToAlphaForce(rCutoff);
slope = (Force2-Force1)/(0.0001);
yIntercept = Force2-slope*(rCutoff+0.0001);
F_MAG(Norms<=rCutoff) = Norms(Norms<=rCutoff)*slope+yIntercept;
F_MAG(Norms==0) = NaN;

F_REP_X = -F_MAG.*(Rx_REP./Norms);
F_REP_Y = -F_MAG.*(Ry_REP./Norms);
if dims==3
    F_REP_Z = -F_MAG.*(Rz_REP./Norms);
end

if dims==2
    F_REP = [nansum(F_REP_X,2),nansum(F_REP_Y,2)];
elseif dims==3
    F_REP = [nansum(F_REP_X,2),nansum(F_REP_Y,2),nansum(F_REP_Z,2)];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_WALL = CalculateWallForces

global Positions dims

% Ly = Corners(3,2)-Corners(2,2);

%Tune Wall Force Parameters so that sigma<<r0 and force is not drastic
F_WALL = zeros(size(Positions(:,1:dims)));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_MAG = CalculateInvToAlphaForce(r)

global ForceScaleREP sigmaREP alpha

F_MAG = alpha*ForceScaleREP*(1/sigmaREP - (sigmaREP^alpha)./(r.^(alpha+1)));
F_MAG(r>sigmaREP) = 0;
F_MAG(isinf(F_MAG)) = NaN;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Force = CalculateCWAForce(r)

global ForceScaleREP sigmaREP

eps = ForceScaleREP;
sig = sigmaREP;
ChanlderWeeksLim = 2^(1/6)*sig;

Force = 4*eps*((6*sig^6)./r.^7-(12*sig^12)./r.^13);
Force(r>ChanlderWeeksLim) = 0;
Force(r==0) = NaN;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_SP = CalculateSpringForces

%Intakes X and Y positions of connected springs to calculate force on each
%node due to springs.

global Rx_SP Ry_SP Rz_SP F_SP_X F_SP_Y F_SP_Z dims

if dims==2
    Norms = (Rx_SP.^2+Ry_SP.^2).^0.5;
elseif dims==3
    Norms = (Rx_SP.^2+Ry_SP.^2+Rz_SP.^2).^0.5;
end

F_MAG = CalculateForceVsStretch(Norms);
F_SP_X = F_MAG.*(Rx_SP./Norms);
F_SP_Y = F_MAG.*(Ry_SP./Norms);
if dims==3
    F_SP_Z = F_MAG.*(Rz_SP./Norms);
end

if dims==2
    F_SP = [nansum(F_SP_X,2),nansum(F_SP_Y,2)];
elseif dims==3
    F_SP = [nansum(F_SP_X,2),nansum(F_SP_Y,2),nansum(F_SP_Z,2)];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F =  CalculateForceVsStretch(Norms)

global SpringType r0 SpringCoeff_eq IsometricForce CoeffRatio L CutoffCWAForce...

if SpringType==0 %Linear
    F = -SpringCoeff_eq*Norms;%(r0-Norms);
elseif SpringType==1 %Langevin 
    F = CalculateLangevinForce(Norms);

    %Use Tangent modulus at high stretch for numerical stability
    StretchLimit = 0.90; %Percent
    StretchLimit_r = StretchLimit*L;
    Force1 = CalculateLangevinForce(StretchLimit_r-0.0001);
    Force2 = CalculateLangevinForce(StretchLimit_r+0.0001);
    slope = (Force2-Force1)/0.0002;
    yIntercept = Force2-slope*(StretchLimit_r+0.0001);
    F(Norms>StretchLimit_r) = Norms(Norms>StretchLimit_r)*slope+yIntercept;    
elseif SpringType==2 % Hill Model 
    F = -IsometricForce*exp(-(((Norms-r0)/r0)/SpringCoeff_eq).^2);
    F(Norms>r0) = -IsometricForce*(exp(-(((Norms(Norms>r0)-r0)/r0)/SpringCoeff_eq).^2)+(((Norms(Norms>r0)-r0)/r0)/(CoeffRatio*SpringCoeff_eq)).^2);
    F(Norms==0) = 0;
elseif SpringType==3
    lam = Norms/L;
    F = CalculateWLCForce(lam);
 
	%Use Tangent modulus at high stretch for numerical stability
    StretchLimit = 0.9;
    Force1 = CalculateWLCForce(StretchLimit-0.0001);
    Force2 = CalculateWLCForce(StretchLimit+0.0001);
    slope = (Force2-Force1)/0.0002;
    yIntercept = Force2-slope*(StretchLimit+0.0001);
    F(lam>StretchLimit) = lam(lam>StretchLimit)*slope+yIntercept;  
elseif SpringType==4    %Entropic Langevin Chains with Chanlder-Weeks-Anderson Repulsion
%     lam = Norms/L;
    FL = CalculateLangevinForce(Norms); %Entropic contribution
    
    %Use Tangent modulus at high stretch for numerical stability
%     StretchLimit = 0.90; %Percent
%     StretchLimit_r = StretchLimit*L;
%     Force1 = CalculateLangevinForce(StretchLimit_r-0.0001);
%     Force2 = CalculateLangevinForce(StretchLimit_r+0.0001);
%     slope = (Force2-Force1)/0.0002;
%     yIntercept = Force2-slope*(StretchLimit_r+0.0001);
%     F(Norms>StretchLimit_r) = Norms(Norms>StretchLimit_r)*slope+yIntercept;   
    
    %Use Tangent modulus at low stretch for numerical stability
    FR = -CalculateCWAForce(Norms);  %Volume exclusion contribution
    
    rCutoff = CutoffCWAForce*r0;
    Force2 = -CalculateCWAForce(rCutoff+0.0001);
    slope = -CalculateCWASlope(rCutoff);
    yIntercept = Force2-slope*(rCutoff+0.0001);
    FR(Norms<=rCutoff) = Norms(Norms<=rCutoff)*slope+yIntercept;
    FR(Norms==0) = NaN;

    F = FR+FL;
%     F = F/F0; %Reduce by a factor of F0
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Force = CalculateLangevinForce(r)

global N_Kuhn b kbT

lambda = r/(sqrt(N_Kuhn)*b);
K = kbT/(sqrt(N_Kuhn)*b);

%ENTER EQUATION FOR LANGEVIN FORCE AS FUNCTION OF STRETCH HERE
%REMEMBER THAT STRETCH IS GOING TO BE A VECTOR OR MATRIX SO USE ".*" AND 
%"./" WHERE NEEDED.

Force = -K*lambda.*(lambda.^2 -  3*N_Kuhn)./(lambda.^2 - N_Kuhn);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ConnectREP

%This function defines all repelling neighbor's in one matrix, "ConnectionsREP". 
%Each row of the matrix correspondes to one node, with the final column 
%listing the node number. The columns from 1 to the maximum allowed number 
%of connections lists the connected neighbors' node numbers.

%This function also defines the corresponding position matrices for the
%nodes listed in "ConnectionsREP" as X_REP and Y_REP. 

%It also defines the components of the connection matrices spring lengths
%Rx_REP and Ry_REP, which will be used to update node positions upon force
%displacement

%Finally, it ID's which connections are accross which bounds via binary
%matrices. E.g., RLPairsREP lists 1 where ConnectionsREP corresponds to a
%spring going from the right to left bound. This is useful for plotting
%purposes

global Lx Ly Lz Positions CutoffPartREP dims...
    LeftBnds RightBnds FrontBnds BackBnds BottomBnds TopBnds...
    BottomLeftBnds BottomRightBnds TopLeftBnds TopRightBnds...
    BottomFrontBnds BottomBackBnds TopFrontBnds TopBackBnds...
    FrontLeftBnds FrontRightBnds BackLeftBnds BackRightBnds...
    C1Bnds C2Bnds C3Bnds C4Bnds C5Bnds C6Bnds C7Bnds C8Bnds...
    ConnectionsREP X_REP Y_REP Z_REP Rx_REP Ry_REP Rz_REP NormsREP...
    OldOrNew
    
% For Bulk
NoNodes = size(Positions,1);

if dims==2
    DomainSize = Lx*Ly;
    Density = NoNodes/DomainSize;
    AreaCutoff = pi*CutoffPartREP^2;
    AvgNoNeighbors = Density*AreaCutoff;
elseif dims==3
    DomainSize = Lx*Ly*Lz;
    Density = NoNodes/DomainSize;
    VolumeCutoff = 4/3*pi*CutoffPartREP^3;
    AvgNoNeighbors = Density*VolumeCutoff;
end
AbsMaxConnections = 2*ceil(AvgNoNeighbors);

NeighborsREP = zeros(NoNodes,AbsMaxConnections+1);
Rx_REP = zeros(NoNodes,AbsMaxConnections);
Ry_REP = zeros(NoNodes,AbsMaxConnections);
Rz_REP = zeros(NoNodes,AbsMaxConnections);
DistREP = zeros(NoNodes,AbsMaxConnections);

NeighborsREP(:,end) = Positions(:,end);

% For Right-Left Bnds
NeighborsREP_RightLeft = zeros(NoNodes,AbsMaxConnections+1);
Rx_REPRightLeft = zeros(NoNodes,AbsMaxConnections);
Ry_REPRightLeft = zeros(NoNodes,AbsMaxConnections);
Rz_REPRightLeft = zeros(NoNodes,AbsMaxConnections);
DistREP_RightLeft = zeros(NoNodes,AbsMaxConnections);

NeighborsREP_RightLeft(:,end) = Positions(:,end);

% For Top-Bottom Bnds
NeighborsREP_TopBottom = zeros(NoNodes,AbsMaxConnections+1);
Rx_REPTopBottom = zeros(NoNodes,AbsMaxConnections);
Ry_REPTopBottom = zeros(NoNodes,AbsMaxConnections);
Rz_REPTopBottom = zeros(NoNodes,AbsMaxConnections);
DistREP_TopBottom = zeros(NoNodes,AbsMaxConnections);

NeighborsREP_TopBottom(:,end) = Positions(:,end);

% For Bottom Left-Top Right Bnds
NeighborsREP_BottomLeftTopRight = zeros(NoNodes,AbsMaxConnections+1);
Rx_REPBottomLeftTopRight = zeros(NoNodes,AbsMaxConnections);
Ry_REPBottomLeftTopRight = zeros(NoNodes,AbsMaxConnections);
Rz_REPBottomLeftTopRight = zeros(NoNodes,AbsMaxConnections);
DistREP_BottomLeftTopRight = zeros(NoNodes,AbsMaxConnections);

NeighborsREP_BottomLeftTopRight(:,end) = Positions(:,end);

% For Bottom Right-Top Left Bnds
NeighborsREP_BottomRightTopLeft = zeros(NoNodes,AbsMaxConnections+1);
Rx_REPBottomRightTopLeft = zeros(NoNodes,AbsMaxConnections);
Ry_REPBottomRightTopLeft = zeros(NoNodes,AbsMaxConnections);
Rz_REPBottomRightTopLeft = zeros(NoNodes,AbsMaxConnections);
DistREP_BottomRightTopLeft = zeros(NoNodes,AbsMaxConnections);

NeighborsREP_BottomRightTopLeft(:,end) = Positions(:,end);

if dims==3
    % 1 Additional Face Combo
    % For Back-Front Bnds
    NeighborsREP_BackFront = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPBackFront = zeros(NoNodes,AbsMaxConnections);
    Ry_REPBackFront = zeros(NoNodes,AbsMaxConnections);
    Rz_REPBackFront = zeros(NoNodes,AbsMaxConnections);
    DistREP_BackFront = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_BackFront(:,end) = Positions(:,end);

    % 4 Additional Edge Combos
    % For Front Left-Back Right Bnds
    NeighborsREP_FrontLeftBackRight = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPFrontLeftBackRight = zeros(NoNodes,AbsMaxConnections);
    Ry_REPFrontLeftBackRight = zeros(NoNodes,AbsMaxConnections);
    Rz_REPFrontLeftBackRight = zeros(NoNodes,AbsMaxConnections);
    DistREP_FrontLeftBackRight = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_FrontLeftBackRight(:,end) = Positions(:,end);

    % For Front Right-Back Left Bnds
    NeighborsREP_FrontRightBackLeft = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPFrontRightBackLeft = zeros(NoNodes,AbsMaxConnections);
    Ry_REPFrontRightBackLeft = zeros(NoNodes,AbsMaxConnections);
    Rz_REPFrontRightBackLeft = zeros(NoNodes,AbsMaxConnections);
    DistREP_FrontRightBackLeft = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_FrontRightBackLeft(:,end) = Positions(:,end);

    % For Bottom Front-Top Back Bnds
    NeighborsREP_BottomFrontTopBack = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPBottomFrontTopBack = zeros(NoNodes,AbsMaxConnections);
    Ry_REPBottomFrontTopBack = zeros(NoNodes,AbsMaxConnections);
    Rz_REPBottomFrontTopBack = zeros(NoNodes,AbsMaxConnections);
    DistREP_BottomFrontTopBack = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_BottomFrontTopBack(:,end) = Positions(:,end);

    % For Top Front-Bottom Back Bnds
    NeighborsREP_TopFrontBottomBack = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPTopFrontBottomBack = zeros(NoNodes,AbsMaxConnections);
    Ry_REPTopFrontBottomBack = zeros(NoNodes,AbsMaxConnections);
    Rz_REPTopFrontBottomBack = zeros(NoNodes,AbsMaxConnections);
    DistREP_TopFrontBottomBack = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_TopFrontBottomBack(:,end) = Positions(:,end);

    % 4 Corner Combos
    % For Left Front Bottom-Right Back Top
    NeighborsREP_C1C7 = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPC1C7 = zeros(NoNodes,AbsMaxConnections);
    Ry_REPC1C7 = zeros(NoNodes,AbsMaxConnections);
    Rz_REPC1C7 = zeros(NoNodes,AbsMaxConnections);
    DistREP_C1C7 = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_C1C7(:,end) = Positions(:,end);

    % For Right Front Bottom-Left Back Top
    NeighborsREP_C2C8 = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPC2C8 = zeros(NoNodes,AbsMaxConnections);
    Ry_REPC2C8 = zeros(NoNodes,AbsMaxConnections);
    Rz_REPC2C8 = zeros(NoNodes,AbsMaxConnections);
    DistREP_C2C8 = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_C2C8(:,end) = Positions(:,end);

    % For Right Front Top-Left Back Bottom
    NeighborsREP_C6C4 = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPC6C4 = zeros(NoNodes,AbsMaxConnections);
    Ry_REPC6C4 = zeros(NoNodes,AbsMaxConnections);
    Rz_REPC6C4 = zeros(NoNodes,AbsMaxConnections);
    DistREP_C6C4 = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_C6C4(:,end) = Positions(:,end);

    % For Left Front Top-Right Back Bottom
    NeighborsREP_C5C3 = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPC5C3 = zeros(NoNodes,AbsMaxConnections);
    Ry_REPC5C3 = zeros(NoNodes,AbsMaxConnections);
    Rz_REPC5C3 = zeros(NoNodes,AbsMaxConnections);
    DistREP_C5C3 = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_C5C3(:,end) = Positions(:,end);
end

% tic
if OldOrNew==1
    % Define Connection Matrices in Bulk
    for P1=1:NoNodes
        Pool = (1:NoNodes)'; Pool(Pool==P1) = [];
        ExistingNeighb = (NeighborsREP(P1,1:end-1))';
        ExistingNeighb(ExistingNeighb==0) = [];
        ExistingNeighb = unique(ExistingNeighb);
        Pool(ismember(Pool,ExistingNeighb)) =[];    %Remove already matched pairs
        Pos1 = Positions(P1,1:dims).*ones(size(Pool,1),dims);
        Pos2 = Positions(Pool,1:dims);
        disp = Pos1-Pos2;
        dist = vecnorm(disp,2,2);
        Neighb = Pool(dist<=CutoffPartREP);
        for i=1:length(Neighb)
            P2 = Neighb(i);
            Pos1 = Positions(P1,1:dims);
            Pos2 = Positions(P2,1:dims);
            distance = Pos1-Pos2; %From P2 to P1
            FirstZero_n = find(NeighborsREP(P1,:)==0,1,'first');
            NeighborsREP(P1,FirstZero_n) = P2;
            Rx_REP(P1,FirstZero_n) = distance(1);
            Ry_REP(P1,FirstZero_n) = distance(2);
            if dims==3
                Rz_REP(P1,FirstZero_n) = distance(3);
            end
            DistREP(P1,FirstZero_n) = norm(distance);

            FirstZero_j = find(NeighborsREP(P2,:)==0,1,'first');
            NeighborsREP(P2,FirstZero_j) = P1;
            Rx_REP(P2,FirstZero_j) = -distance(1);
            Ry_REP(P2,FirstZero_j) = -distance(2);
            if dims==3
                Rz_REP(P2,FirstZero_j) = -distance(3);
            end
            DistREP(P2,FirstZero_j) = norm(distance);
        end
    end
else
    for P1=1:NoNodes-1
        for P2=P1+1:NoNodes
            Pos1 = Positions(P1,1:dims);
            Pos2 = Positions(P2,1:dims);
            distance = Pos1-Pos2; %From P2 to P1
            if norm(distance) <= CutoffPartREP %Find Potential Interactions
                FirstZero_n = find(NeighborsREP(P1,:)==0,1,'first');
                NeighborsREP(P1,FirstZero_n) = P2;
                Rx_REP(P1,FirstZero_n) = distance(1);
                Ry_REP(P1,FirstZero_n) = distance(2);
                if dims==3
                    Rz_REP(P1,FirstZero_n) = distance(3);
                end
                DistREP(P1,FirstZero_n) = norm(distance);

                FirstZero_j = find(NeighborsREP(P2,:)==0,1,'first');
                NeighborsREP(P2,FirstZero_j) = P1;
                Rx_REP(P2,FirstZero_j) = -distance(1);
                Ry_REP(P2,FirstZero_j) = -distance(2);
                if dims==3
                    Rz_REP(P2,FirstZero_j) = -distance(3);
                end
                DistREP(P2,FirstZero_j) = norm(distance);
            end
        end
    end
end
% toc
% CutoffPartREP = 10; %% DELETE
%Check accross Right-Left Bounds
Bnds1 = RightBnds;
Bnds2 = LeftBnds;
N1 = size(Bnds1,1);
N2 = size(Bnds2,1);
Shift = zeros(1,dims);
Shift(1) = Lx;
[NeighborsREP_RightLeft,~,...
Rx_REPRightLeft,Ry_REPRightLeft,Rz_REPRightLeft] = ...
    FindRepNeighbsBnds(NeighborsREP_RightLeft,...
    Rx_REPRightLeft,Ry_REPRightLeft,Rz_REPRightLeft,...
    DistREP_RightLeft,N1,N2,Bnds1,Bnds2,Shift);

%Check accross Top-Bottom Bounds
Bnds1 = TopBnds;
Bnds2 = BottomBnds;
N1 = size(Bnds1,1);
N2 = size(Bnds2,1);
Shift = zeros(1,dims);
if dims==2
    Shift(2) = Ly;
elseif dims==3
    Shift(3) = Lz;
end
[NeighborsREP_TopBottom,~,...
Rx_REPTopBottom,Ry_REPTopBottom,Rz_REPTopBottom] = ...
    FindRepNeighbsBnds(NeighborsREP_TopBottom,...
    Rx_REPTopBottom,Ry_REPTopBottom,Rz_REPTopBottom,...
    DistREP_TopBottom,N1,N2,Bnds1,Bnds2,Shift);

%Check accross Bottom Left-Top Right Bounds
Bnds1 = BottomLeftBnds;
Bnds2 = TopRightBnds;
N1 = size(Bnds1,1);
N2 = size(Bnds2,1);
Shift = zeros(1,dims);
if dims==2
    Shift(1) = -Lx;
    Shift(2) = -Ly;
elseif dims==3
    Shift(1) = -Lx;
    Shift(3) = -Lz;
end
[NeighborsREP_BottomLeftTopRight,~,...
Rx_REPBottomLeftTopRight,Ry_REPBottomLeftTopRight,Rz_REPBottomLeftTopRight] = ...
    FindRepNeighbsBnds(NeighborsREP_BottomLeftTopRight,...
    Rx_REPBottomLeftTopRight,Ry_REPBottomLeftTopRight,Rz_REPBottomLeftTopRight,...
    DistREP_BottomLeftTopRight,N1,N2,Bnds1,Bnds2,Shift);

%Check accross Top Left-Bottom Right Bounds
Bnds1 = BottomRightBnds;
Bnds2 = TopLeftBnds;
N1 = size(Bnds1,1);
N2 = size(Bnds2,1);
Shift = zeros(1,dims);
if dims==2
    Shift(1) = +Lx;
    Shift(2) = -Ly;
elseif dims==3
    Shift(1) = +Lx;
    Shift(3) = -Lz;
end
[NeighborsREP_BottomRightTopLeft,~,...
Rx_REPBottomRightTopLeft,Ry_REPBottomRightTopLeft,Rz_REPBottomRightTopLeft] = ...
    FindRepNeighbsBnds(NeighborsREP_BottomRightTopLeft,...
    Rx_REPBottomRightTopLeft,Ry_REPBottomRightTopLeft,Rz_REPBottomRightTopLeft,...
    DistREP_BottomRightTopLeft,N1,N2,Bnds1,Bnds2,Shift);

if dims==3
    % 1 Additional Face Combo
    % Check accross Back-Front bounds
    Bnds1 = FrontBnds;
    Bnds2 = BackBnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(2) = -Ly;
    [NeighborsREP_BackFront,~,...
        Rx_REPBackFront,Ry_REPBackFront,Rz_REPBackFront] = ...
        FindRepNeighbsBnds(NeighborsREP_BackFront,...
        Rx_REPBackFront,Ry_REPBackFront,Rz_REPBackFront,...
        DistREP_BackFront,N1,N2,Bnds1,Bnds2,Shift);

    % 4 Additional Edge Combos
    % Check accross Front Left-Back Right Bnds
    Bnds1 = FrontLeftBnds;
    Bnds2 = BackRightBnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(1) = -Lx;
    Shift(2) = -Ly;
    [NeighborsREP_FrontLeftBackRight,~,...
        Rx_REPFrontLeftBackRight,Ry_REPFrontLeftBackRight,Rz_REPFrontLeftBackRight] = ...
        FindRepNeighbsBnds(NeighborsREP_FrontLeftBackRight,...
        Rx_REPFrontLeftBackRight,Ry_REPFrontLeftBackRight,Rz_REPFrontLeftBackRight,...
        DistREP_FrontLeftBackRight,N1,N2,Bnds1,Bnds2,Shift);

    % Check accross Front Right-Back Left Bnds
    Bnds1 = FrontRightBnds;
    Bnds2 = BackLeftBnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(1) = +Lx;
    Shift(2) = -Ly;
    [NeighborsREP_FrontRightBackLeft,~,...
        Rx_REPFrontRightBackLeft,Ry_REPFrontRightBackLeft,Rz_REPFrontRightBackLeft] = ...
        FindRepNeighbsBnds(NeighborsREP_FrontRightBackLeft,...
        Rx_REPFrontRightBackLeft,Ry_REPFrontRightBackLeft,Rz_REPFrontRightBackLeft,...
        DistREP_FrontRightBackLeft,N1,N2,Bnds1,Bnds2,Shift);

    % Check accross Bottom Front-Top Back Bnds
    Bnds1 = BottomFrontBnds;
    Bnds2 = TopBackBnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(2) = -Ly;
    Shift(3) = -Lz;
    [NeighborsREP_BottomFrontTopBack,~,...
        Rx_REPBottomFrontTopBack,Ry_REPBottomFrontTopBack,Rz_REPBottomFrontTopBack] = ...
        FindRepNeighbsBnds(NeighborsREP_BottomFrontTopBack,...
        Rx_REPBottomFrontTopBack,Ry_REPBottomFrontTopBack,Rz_REPBottomFrontTopBack,...
        DistREP_BottomFrontTopBack,N1,N2,Bnds1,Bnds2,Shift);

    % Check accross Front Top-Back Bottom Bnds
    Bnds1 = TopFrontBnds;
    Bnds2 = BottomBackBnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(2) = -Ly;
    Shift(3) = +Lz;
    [NeighborsREP_TopFrontBottomBack,~,...
        Rx_REPTopFrontBottomBack,Ry_REPTopFrontBottomBack,Rz_REPTopFrontBottomBack] = ...
        FindRepNeighbsBnds(NeighborsREP_TopFrontBottomBack,...
        Rx_REPTopFrontBottomBack,Ry_REPTopFrontBottomBack,Rz_REPTopFrontBottomBack,...
        DistREP_TopFrontBottomBack,N1,N2,Bnds1,Bnds2,Shift);

    % 4 Corners
    % Front Bottom Left-Back Top Right
    Bnds1 = C1Bnds;
    Bnds2 = C7Bnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(1) = -Lx;
    Shift(2) = -Ly;
    Shift(3) = -Lz;
    [NeighborsREP_C1C7,~,...
        Rx_REPC1C7,Ry_REPC1C7,Rz_REPC1C7] = ...
        FindRepNeighbsBnds(NeighborsREP_C1C7,...
        Rx_REPC1C7,Ry_REPC1C7,Rz_REPC1C7,...
        DistREP_C1C7,N1,N2,Bnds1,Bnds2,Shift);

    % Front Bottom Right-Back Top Left
    Bnds1 = C2Bnds;
    Bnds2 = C8Bnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(1) = +Lx;
    Shift(2) = -Ly;
    Shift(3) = -Lz;
    [NeighborsREP_C2C8,~,...
        Rx_REPC2C8,Ry_REPC2C8,Rz_REPC2C8] = ...
        FindRepNeighbsBnds(NeighborsREP_C2C8,...
        Rx_REPC2C8,Ry_REPC2C8,Rz_REPC2C8,...
        DistREP_C2C8,N1,N2,Bnds1,Bnds2,Shift);

    % Front Top Right-Back Bottom Left
    Bnds1 = C6Bnds;
    Bnds2 = C4Bnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(1) = +Lx;
    Shift(2) = -Ly;
    Shift(3) = +Lz;
    [NeighborsREP_C6C4,~,...
        Rx_REPC6C4,Ry_REPC6C4,Rz_REPC6C4] = ...
        FindRepNeighbsBnds(NeighborsREP_C6C4,...
        Rx_REPC6C4,Ry_REPC6C4,Rz_REPC6C4,...
        DistREP_C6C4,N1,N2,Bnds1,Bnds2,Shift);

    % Front Top Left-Back Bottom Right
    Bnds1 = C5Bnds;
    Bnds2 = C3Bnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(1) = -Lx;
    Shift(2) = -Ly;
    Shift(3) = +Lz;
    [NeighborsREP_C5C3,~,...
        Rx_REPC5C3,Ry_REPC5C3,Rz_REPC5C3] = ...
        FindRepNeighbsBnds(NeighborsREP_C5C3,...
        Rx_REPC5C3,Ry_REPC5C3,Rz_REPC5C3,...
        DistREP_C5C3,N1,N2,Bnds1,Bnds2,Shift);
end

% Sort Connected Springs by shortest to longest
if dims==2
    ConnectionsREP = [NeighborsREP(:,1:end-1),...
        NeighborsREP_RightLeft(:,1:end-1),...
        NeighborsREP_TopBottom(:,1:end-1),...
        NeighborsREP_BottomLeftTopRight(:,1:end-1),...
        NeighborsREP_BottomRightTopLeft(:,1:end-1)];
    Rx_REP = [Rx_REP,...
        Rx_REPRightLeft,...
        Rx_REPTopBottom,...
        Rx_REPBottomLeftTopRight,...
        Rx_REPBottomRightTopLeft];
    Ry_REP = [Ry_REP,...
        Ry_REPRightLeft,...
        Ry_REPTopBottom,...
        Ry_REPBottomLeftTopRight,...
        Ry_REPBottomRightTopLeft];
    Rz_REP = [Rz_REP,...
        Rz_REPRightLeft,...
        Rz_REPTopBottom,...
        Rz_REPBottomLeftTopRight,...
        Rz_REPBottomRightTopLeft];
else
    ConnectionsREP = [NeighborsREP(:,1:end-1),...
        NeighborsREP_RightLeft(:,1:end-1),...
        NeighborsREP_TopBottom(:,1:end-1),...
        NeighborsREP_BackFront(:,1:end-1),...
        NeighborsREP_BottomLeftTopRight(:,1:end-1),...
        NeighborsREP_BottomRightTopLeft(:,1:end-1),...
        NeighborsREP_BottomFrontTopBack(:,1:end-1),...
        NeighborsREP_TopFrontBottomBack(:,1:end-1),...
        NeighborsREP_FrontLeftBackRight(:,1:end-1),...
        NeighborsREP_FrontRightBackLeft(:,1:end-1),...
        NeighborsREP_C1C7(:,1:end-1),...
        NeighborsREP_C2C8(:,1:end-1),...
        NeighborsREP_C6C4(:,1:end-1),...
        NeighborsREP_C5C3(:,1:end-1)];
    Rx_REP = [Rx_REP,...
        Rx_REPRightLeft,...
        Rx_REPTopBottom,...
        Rx_REPBackFront,...
        Rx_REPBottomLeftTopRight,...
        Rx_REPBottomRightTopLeft,...
        Rx_REPBottomFrontTopBack,...
        Rx_REPTopFrontBottomBack,...
        Rx_REPFrontLeftBackRight,...
        Rx_REPFrontRightBackLeft,...
        Rx_REPC1C7,...
        Rx_REPC2C8,...
        Rx_REPC6C4,...
        Rx_REPC5C3];
    Ry_REP = [Ry_REP,...
        Ry_REPRightLeft,...
        Ry_REPTopBottom,...
        Ry_REPBackFront,...
        Ry_REPBottomLeftTopRight,...
        Ry_REPBottomRightTopLeft,...
        Ry_REPBottomFrontTopBack,...
        Ry_REPTopFrontBottomBack,...
        Ry_REPFrontLeftBackRight,...
        Ry_REPFrontRightBackLeft,...
        Ry_REPC1C7,...
        Ry_REPC2C8,...
        Ry_REPC6C4,...
        Ry_REPC5C3];
    Rz_REP = [Rz_REP,...
        Rz_REPRightLeft,...
        Rz_REPTopBottom,...
        Rz_REPBackFront,...
        Rz_REPBottomLeftTopRight,...
        Rz_REPBottomRightTopLeft,...
        Rz_REPBottomFrontTopBack,...
        Rz_REPTopFrontBottomBack,...
        Rz_REPFrontLeftBackRight,...
        Rz_REPFrontRightBackLeft,...
        Rz_REPC1C7,...
        Rz_REPC2C8,...
        Rz_REPC6C4,...
        Rz_REPC5C3];
end
if dims==2
    NormsREP = (Rx_REP.^2+Ry_REP.^2).^0.5;
else
    NormsREP = (Rx_REP.^2+Ry_REP.^2+Rz_REP.^2).^0.5;
end   

Nodes = (1:NoNodes)';
ConnectionsREP(:,end) = Nodes;

[X_REP,Y_REP,Z_REP] = DefinePositionMatrix(ConnectionsREP);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,Z] = DefinePositionMatrix(Connections)

global Positions dims

X=zeros(size(Connections));
Y=zeros(size(Connections));
Z=zeros(size(Connections));

X(:,end) = Positions(:,1);
Y(:,end) = Positions(:,2);
if dims==3
    Z(:,end) = Positions(:,3);
end

for j=1:size(X,2)-1
    TempNodes = Connections(:,j);
    NonZeroIndices = find(TempNodes~=0);
    
    X(NonZeroIndices,j) = Positions(TempNodes(NonZeroIndices),1); 
    Y(NonZeroIndices,j) = Positions(TempNodes(NonZeroIndices),2); 
    if dims==3
        Z(NonZeroIndices,j) = Positions(TempNodes(NonZeroIndices),3);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NeighborsREP,DistREP,Rx_REP,Ry_REP,Rz_REP] = ...
    FindRepNeighbsBnds(NeighborsREP,Rx_REP,Ry_REP,Rz_REP,DistREP,...
    N1,N2,Bnds1,Bnds2,Shift)

global dims CutoffPartREP Positions OldOrNew

if OldOrNew==1
    for i=1:N1
        P1 = Bnds1(i);
        Pool = Bnds2; Pool(Pool==P1) = [];
        ExistingNeighb = (NeighborsREP(P1,1:end-1))';
        ExistingNeighb(ExistingNeighb==0) = [];
        ExistingNeighb = unique(ExistingNeighb);
        Pool(ismember(Pool,ExistingNeighb)) =[];    %Remove already matched pairs
        Pos1 = Positions(P1,1:dims).*ones(size(Pool,1),dims);
        Pos2 = Positions(Pool,1:dims) + Shift;
        disp = Pos1-Pos2;
        dist = vecnorm(disp,2,2);
        Neighb = Pool(dist<=CutoffPartREP);
        for j=1:length(Neighb)
            P2 = Neighb(j);
            Pos1 = Positions(P1,1:dims);
            Pos2 = Positions(P2,1:dims) + Shift;
            distance = Pos1-Pos2; %From P2 to P1
            FirstZero_n = find(NeighborsREP(P1,:)==0,1,'first');
            NeighborsREP(P1,FirstZero_n) = P2;
            Rx_REP(P1,FirstZero_n) = distance(1);
            Ry_REP(P1,FirstZero_n) = distance(2);
            if dims==3
                Rz_REP(P1,FirstZero_n) = distance(3);
            end
            DistREP(P1,FirstZero_n) = norm(distance);

            FirstZero_j = find(NeighborsREP(P2,:)==0,1,'first');
            NeighborsREP(P2,FirstZero_j) = P1;
            Rx_REP(P2,FirstZero_j) = -distance(1);
            Ry_REP(P2,FirstZero_j) = -distance(2);
            if dims==3
                Rz_REP(P2,FirstZero_j) = -distance(3);
            end
            DistREP(P2,FirstZero_j) = norm(distance);
        end
    end
else
    for i=1:N1
        P1 = Bnds1(i);
        for j=1:N2
            P2 = Bnds2(j);
            Pos1 = Positions(P1,1:dims);
            Pos2 = Positions(P2,1:dims) + Shift;
            distance = Pos1-Pos2; %From P2 to P1
            if norm(distance) <= CutoffPartREP %Find Potential Interactions
                FirstZero_n = find(NeighborsREP(P1,:)==0,1,'first');
                NeighborsREP(P1,FirstZero_n) = P2;
                Rx_REP(P1,FirstZero_n) = distance(1);
                Ry_REP(P1,FirstZero_n) = distance(2);
                if dims==3
                    Rz_REP(P1,FirstZero_n) = distance(3);
                end
                DistREP(P1,FirstZero_n) = norm(distance);

                FirstZero_j = find(NeighborsREP(P2,:)==0,1,'first');
                NeighborsREP(P2,FirstZero_j) = P1;
                Rx_REP(P2,FirstZero_j) = -distance(1);
                Ry_REP(P2,FirstZero_j) = -distance(2);
                if dims==3
                    Rz_REP(P2,FirstZero_j) = -distance(3);
                end
                DistREP(P2,FirstZero_j) = norm(distance);
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ReplicateBounds

%This function uses the updated boundary node numbers to replicate these
%crosslinks accross the periodic boundaries

global dims Positions Corners...
    LeftBnds RightBnds FrontBnds BackBnds BottomBnds TopBnds...
    BottomLeftBnds BottomRightBnds TopLeftBnds TopRightBnds...
    BottomFrontBnds BottomBackBnds TopFrontBnds TopBackBnds...
    FrontLeftBnds FrontRightBnds BackLeftBnds BackRightBnds...
    C1Bnds C2Bnds C3Bnds C4Bnds C5Bnds C6Bnds C7Bnds C8Bnds...
    PosLeftClns PosRightClns PosFrontClns PosBackClns PosBottomClns PosTopClns...
    PosBottomLeftClns PosBottomRightClns PosTopLeftClns PosTopRightClns...
    PosBottomFrontClns PosBottomBackClns PosTopFrontClns PosTopBackClns...
    PosFrontLeftClns PosFrontRightClns PosBackLeftClns PosBackRightClns...
    PosC1Clns PosC2Clns PosC3Clns PosC4Clns PosC5Clns PosC6Clns PosC7Clns PosC8Clns

%Replicate boundary nodes
if dims==2
    CBL = Corners(1,:);
    CBR = Corners(2,:);
    CTR = Corners(3,:);
    
    X = CBR(1)-CBL(1);
    Y = CTR(2)-CBR(2);
    dXdY = CTR(1)-CBR(1);
    dYdX = CBR(2)-CBL(2);

    % 4 Edges
    PosLeftClns = [Positions(LeftBnds,1) + X,...
                Positions(LeftBnds,2) + dYdX]; 
    PosRightClns = [Positions(RightBnds,1) - X,...
                Positions(RightBnds,2) - dYdX]; 
    PosBottomClns = [Positions(BottomBnds,1) + dXdY,...
                Positions(BottomBnds,2) + Y];
    PosTopClns = [Positions(TopBnds,1) - dXdY,...
                Positions(TopBnds,2) - Y];

    % 4 Corners
    PosBottomLeftClns = [Positions(BottomLeftBnds,1) + X + dXdY,...
                 Positions(BottomLeftBnds,2) + Y + dYdX];
    PosBottomRightClns = [Positions(BottomRightBnds,1) - X + dXdY,...
                 Positions(BottomRightBnds,2) + Y - dYdX];
    PosTopRightClns = [Positions(TopRightBnds,1) - X - dXdY,...
                 Positions(TopRightBnds,2) - Y - dYdX];
    PosTopLeftClns = [Positions(TopLeftBnds,1) + X - dXdY,...
                 Positions(TopLeftBnds,2) - Y + dYdX];
elseif dims==3 % ORTHONORMAL DEFORMATIONS ONLY
    C1 = Corners(1,:);  %Front-Bottom-Left
    C2 = Corners(2,:);  %Front-Bottom-Right
    C4 = Corners(4,:);  %Back-Bottom-Left
    C5 = Corners(5,:);  %Front-Top-Left

    X = C2(1)-C1(1);
    Y = C4(2)-C1(2);
    Z = C5(3)-C1(3);

    % 6 Faces
    PosLeftClns = [Positions(LeftBnds,1) + X,...
        Positions(LeftBnds,2),...
        Positions(LeftBnds,3)];
    PosRightClns = [Positions(RightBnds,1) - X,...
        Positions(RightBnds,2),...
        Positions(RightBnds,3)];
    PosFrontClns = [Positions(FrontBnds,1),...
        Positions(FrontBnds,2) + Y,...
        Positions(FrontBnds,3)];
    PosBackClns = [Positions(BackBnds,1),...
        Positions(BackBnds,2) - Y,...
        Positions(BackBnds,3)];
    PosBottomClns = [Positions(BottomBnds,1),...
        Positions(BottomBnds,2),...
        Positions(BottomBnds,3) + Z];
    PosTopClns = [Positions(TopBnds,1),...
        Positions(TopBnds,2),...
        Positions(TopBnds,3) - Z];

    % 12 Edges
    PosBottomLeftClns = [Positions(BottomLeftBnds,1) + X,...
        Positions(BottomLeftBnds,2),...
        Positions(BottomLeftBnds,3) + Z];
    PosBottomRightClns = [Positions(BottomRightBnds,1) - X,...
        Positions(BottomRightBnds,2),...
        Positions(BottomRightBnds,3) + Z];
    PosTopLeftClns = [Positions(TopLeftBnds,1) + X,...
        Positions(TopLeftBnds,2),...
        Positions(TopLeftBnds,3) - Z];
    PosTopRightClns = [Positions(TopRightBnds,1) - X,...
        Positions(TopRightBnds,2),...
        Positions(TopRightBnds,3) - Z];

    PosBottomFrontClns = [Positions(BottomFrontBnds,1),...
        Positions(BottomFrontBnds,2) + Y,...
        Positions(BottomFrontBnds,3) + Z];
    PosBottomBackClns = [Positions(BottomBackBnds,1),...
        Positions(BottomBackBnds,2) - Y,...
        Positions(BottomBackBnds,3) + Z];
    PosTopFrontClns = [Positions(TopFrontBnds,1),...
        Positions(TopFrontBnds,2) + Y,...
        Positions(TopFrontBnds,3) - Z];
    PosTopBackClns = [Positions(TopBackBnds,1),...
        Positions(TopBackBnds,2) - Y,...
        Positions(TopBackBnds,3) - Z];

    PosFrontLeftClns = [Positions(FrontLeftBnds,1) + X,...
        Positions(FrontLeftBnds,2) + Y,...
        Positions(FrontLeftBnds,3)];
    PosFrontRightClns = [Positions(FrontRightBnds,1) - X,...
        Positions(FrontRightBnds,2) + Y,...
        Positions(FrontRightBnds,3)];
    PosBackLeftClns = [Positions(BackLeftBnds,1) + X,...
        Positions(BackLeftBnds,2) - Y,...
        Positions(BackLeftBnds,3)];
    PosBackRightClns = [Positions(BackRightBnds,1) - X,...
        Positions(BackRightBnds,2) - Y,...
        Positions(BackRightBnds,3)];

    % 8 Corners
    PosC1Clns = [Positions(C1Bnds,1) + X,...
        Positions(C1Bnds,2) + Y,...
        Positions(C1Bnds,3) + Z];
    PosC2Clns = [Positions(C2Bnds,1) - X,...
        Positions(C2Bnds,2) + Y,...
        Positions(C2Bnds,3) + Z];
    PosC3Clns = [Positions(C3Bnds,1) - X,...
        Positions(C3Bnds,2) - Y,...
        Positions(C3Bnds,3) + Z];
    PosC4Clns = [Positions(C4Bnds,1) + X,...
        Positions(C4Bnds,2) - Y,...
        Positions(C4Bnds,3) + Z];
    PosC5Clns = [Positions(C5Bnds,1) + X,...
        Positions(C5Bnds,2) + Y,...
        Positions(C5Bnds,3) - Z];
    PosC6Clns = [Positions(C6Bnds,1) - X,...
        Positions(C6Bnds,2) + Y,...
        Positions(C6Bnds,3) - Z];
    PosC7Clns = [Positions(C7Bnds,1) - X,...
        Positions(C7Bnds,2) - Y,...
        Positions(C7Bnds,3) - Z];
    PosC8Clns = [Positions(C8Bnds,1) + X,...
        Positions(C8Bnds,2) - Y,...
        Positions(C8Bnds,3) - Z];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FindBoundaryNodes

%This funciton identifies the node numbers corresponding to each boundary
%as a function of their positions.

global Lx Ly Lz Positions Corners dims...
    LeftBnds RightBnds FrontBnds BackBnds BottomBnds TopBnds...
    BottomLeftBnds BottomRightBnds TopLeftBnds TopRightBnds...
    BottomFrontBnds BottomBackBnds TopFrontBnds TopBackBnds...
    FrontLeftBnds FrontRightBnds BackLeftBnds BackRightBnds...
    C1Bnds C2Bnds C3Bnds C4Bnds C5Bnds C6Bnds C7Bnds C8Bnds...
    AllBnds

MultFactor = 0.5;

NoNodes = size(Positions,1);
if dims==2
    C1 = Corners(1,:);  %BL
    C2 = Corners(2,:);  %BR
    C3 = Corners(3,:); %TR
    C4 = Corners(4,:); %TL
    
    BottomBnds = PointSlopeProximityFind(C2,C1,MultFactor);
    TopBnds = PointSlopeProximityFind(C4,C3,MultFactor);
    RightBnds = PointSlopeProximityFind(C2,C3,MultFactor);
    LeftBnds = PointSlopeProximityFind(C4,C1,MultFactor);
    
    BottomLeftBnds = BottomBnds(ismember(BottomBnds,LeftBnds));
    BottomRightBnds = BottomBnds(ismember(BottomBnds,RightBnds));
    TopRightBnds = TopBnds(ismember(TopBnds,RightBnds));
    TopLeftBnds = TopBnds(ismember(TopBnds,LeftBnds));
elseif dims==3
    % CONVENTION: 
    % Left-to-Right is (-) to (+) x-axis
    % Front-to-Back is (-) to (+) y-axis
    % Bottom-to-Top is (-) to (+) z-axis
    C1 = Corners(1,:);  %LFB    Left-Front-Bottom
    C7 = Corners(7,:);  %RBT    Right-Back-Top
    
    %For orthotropic boundaries only
    xNorm = [ones(NoNodes,1) zeros(NoNodes,2)];
    yNorm = [zeros(NoNodes,1) ones(NoNodes,1) zeros(NoNodes,1)];
    zNorm = [zeros(NoNodes,2) ones(NoNodes,1)];

    DistToLeft = dot(Positions(:,1:dims)-C1,xNorm,2);
    DistToRight = -dot(Positions(:,1:dims)-C7,xNorm,2);
    DistToFront = dot(Positions(:,1:dims)-C1,yNorm,2);
    DistToBack = -dot(Positions(:,1:dims)-C7,yNorm,2);
    DistToBottom = dot(Positions(:,1:dims)-C1,zNorm,2);
    DistToTop = -dot(Positions(:,1:dims)-C7,zNorm,2);

    % 6 Faces
    LeftBnds = find(DistToLeft<MultFactor*Lx/2);
    RightBnds = find(DistToRight<MultFactor*Lx/2);
    FrontBnds = find(DistToFront<MultFactor*Ly/2);
    BackBnds = find(DistToBack<MultFactor*Ly/2);
    BottomBnds = find(DistToBottom<MultFactor*Lz/2);
    TopBnds = find(DistToTop<MultFactor*Lz/2);

    % 12 Edge Overlaps
    BottomLeftBnds = BottomBnds(ismember(BottomBnds,LeftBnds));
    BottomRightBnds = BottomBnds(ismember(BottomBnds,RightBnds));
    TopLeftBnds = TopBnds(ismember(TopBnds,LeftBnds));
    TopRightBnds = TopBnds(ismember(TopBnds,RightBnds));

    BottomFrontBnds = BottomBnds(ismember(BottomBnds,FrontBnds));
    BottomBackBnds = BottomBnds(ismember(BottomBnds,BackBnds));
    TopFrontBnds = TopBnds(ismember(TopBnds,FrontBnds));
    TopBackBnds = TopBnds(ismember(TopBnds,BackBnds));

    FrontLeftBnds = FrontBnds(ismember(FrontBnds,LeftBnds));
    BackLeftBnds = BackBnds(ismember(BackBnds,LeftBnds));
    FrontRightBnds = FrontBnds(ismember(FrontBnds,RightBnds));
    BackRightBnds = BackBnds(ismember(BackBnds,RightBnds));

    % 8 Corner Overlaps
    C1Bnds = FrontLeftBnds(ismember(FrontLeftBnds,BottomFrontBnds));
    C2Bnds = FrontRightBnds(ismember(FrontRightBnds,BottomRightBnds));
    C3Bnds = BackRightBnds(ismember(BackRightBnds,BottomBackBnds));
    C4Bnds = BackLeftBnds(ismember(BackLeftBnds,BottomBackBnds));
    C5Bnds = FrontLeftBnds(ismember(FrontLeftBnds,TopFrontBnds));
    C6Bnds = FrontRightBnds(ismember(FrontRightBnds,TopFrontBnds));
    C7Bnds = BackRightBnds(ismember(BackRightBnds,TopBackBnds));
    C8Bnds = BackLeftBnds(ismember(BackLeftBnds,TopLeftBnds));

    ShowFigDebug = 0;
    CheckSurfs = 0;
    CheckEdges = 0;
    CheckCorners = 0;
    if ShowFigDebug==1
        if CheckSurfs==1
        hold on
        s = scatter3(Positions(LeftBnds,1),Positions(LeftBnds,2),Positions(LeftBnds,3),'r','filled');
        s.SizeData = 6;
        s = scatter3(Positions(RightBnds,1),Positions(RightBnds,2),Positions(RightBnds,3),'r','filled');
        s.SizeData = 6;

        s = scatter3(Positions(TopBnds,1),Positions(TopBnds,2),Positions(TopBnds,3),'r','filled');
        s.SizeData = 6;
        s = scatter3(Positions(BottomBnds,1),Positions(BottomBnds,2),Positions(BottomBnds,3),'r','filled');
        s.SizeData = 6;

        s = scatter3(Positions(FrontBnds,1),Positions(FrontBnds,2),Positions(FrontBnds,3),'r','filled');
        s.SizeData = 6;
        s = scatter3(Positions(BackBnds,1),Positions(BackBnds,2),Positions(BackBnds,3),'r','filled');
        s.SizeData = 6;
        end

        if CheckEdges==1
            s = scatter3(Positions(FrontLeftBnds,1),Positions(FrontLeftBnds,2),Positions(FrontLeftBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(FrontRightBnds,1),Positions(FrontRightBnds,2),Positions(FrontRightBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(BackLeftBnds,1),Positions(BackLeftBnds,2),Positions(BackLeftBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(BackRightBnds,1),Positions(BackRightBnds,2),Positions(BackRightBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(BottomLeftBnds,1),Positions(BottomLeftBnds,2),Positions(BottomLeftBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(BottomRightBnds,1),Positions(BottomRightBnds,2),Positions(BottomRightBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(TopLeftBnds,1),Positions(TopLeftBnds,2),Positions(TopLeftBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(TopRightBnds,1),Positions(TopRightBnds,2),Positions(TopRightBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(BottomFrontBnds,1),Positions(BottomFrontBnds,2),Positions(BottomFrontBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(BottomBackBnds,1),Positions(BottomBackBnds,2),Positions(BottomBackBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(TopFrontBnds,1),Positions(TopFrontBnds,2),Positions(TopFrontBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(TopBackBnds,1),Positions(TopBackBnds,2),Positions(TopBackBnds,3),'r','filled');
            s.SizeData = 6;
        end

        if CheckCorners==1
            if ~isempty(C1Bnds)
                s = scatter3(Positions(C1Bnds,1),Positions(C1Bnds,2),Positions(C1Bnds,3),'r','filled');
                s.SizeData = 12;
            end
            if ~isempty(C2Bnds)
                s = scatter3(Positions(C2Bnds,1),Positions(C2Bnds,2),Positions(C2Bnds,3),'r','filled');
                s.SizeData = 12;
            end
            if ~isempty(C3Bnds)
                s = scatter3(Positions(C3Bnds,1),Positions(C3Bnds,2),Positions(C3Bnds,3),'r','filled');
                s.SizeData = 12;
            end
            if ~isempty(C4Bnds)
                s = scatter3(Positions(C4Bnds,1),Positions(C4Bnds,2),Positions(C4Bnds,3),'r','filled');
                s.SizeData = 12;
            end
            if ~isempty(C5Bnds)
                s = scatter3(Positions(C5Bnds,1),Positions(C5Bnds,2),Positions(C5Bnds,3),'r','filled');
                s.SizeData = 12;
            end
            if ~isempty(C6Bnds)
                s = scatter3(Positions(C6Bnds,1),Positions(C6Bnds,2),Positions(C6Bnds,3),'r','filled');
                s.SizeData = 12;
            end
            if ~isempty(C7Bnds)
                s = scatter3(Positions(C7Bnds,1),Positions(C7Bnds,2),Positions(C7Bnds,3),'r','filled');
                s.SizeData = 12;
            end
            if ~isempty(C8Bnds)
                s = scatter3(Positions(C8Bnds,1),Positions(C8Bnds,2),Positions(C8Bnds,3),'r','filled');
                s.SizeData = 12;
            end
        end
    end
end

if dims==2
    AllBnds = [LeftBnds;RightBnds;BottomBnds;TopBnds;...
        BottomLeftBnds;BottomRightBnds;TopLeftBnds;TopRightBnds];
elseif dims==3
    AllBnds = [LeftBnds;RightBnds;FrontBnds;BackBnds;BottomBnds;TopBnds;...
        BottomLeftBnds;BottomRightBnds;TopLeftBnds;TopRightBnds;...
        BottomFrontBnds;BottomBackBnds;TopFrontBnds;TopBackBnds;...
        FrontLeftBnds;FrontRightBnds;BackLeftBnds;BackRightBnds;...
        C1Bnds;C2Bnds;C3Bnds;C4Bnds;C5Bnds;C6Bnds;C7Bnds;C8Bnds];
end
AllBnds = unique(AllBnds);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Bnds = PointSlopeProximityFind(C1,C2,MultFactor)

global Positions Corners

CBL = Corners(1,:);
CBR = Corners(2,:);
CTR = Corners(3,:);

%Define min distance between any two corners to determine boundary
%proximity
Lx = CBR(1)-CBL(1);
Ly = CTR(2)-CBR(2);

%Find points on edge closest to Positions
if C1(2)~=C2(2) && C1(1)~=C2(1)
    Slope = (C1(2)-C2(2))/(C1(1)-C2(1));
    
    K = Slope+1/Slope;
    Px = (1/K)*(Positions(:,2)-C1(2)+1/Slope*Positions(:,1)+Slope*C1(1));
    Py = Positions(:,2) + 1/Slope*(Positions(:,1)-Px);
    
    XDistToOrigin = 1/K*(Slope*C1(1)-C1(2));
    YDistToOrigin = -1/Slope*XDistToOrigin;
    DistToOrigin = [XDistToOrigin,YDistToOrigin];
elseif C1(1)==C2(1) %If a vertical line (slope == Infty)
    Px = C1(1)*ones(size(Positions,1),1);
    Py = Positions(:,2);
    DistToOrigin = [Lx/2,0];
elseif C1(2)==C2(2) %If a horizontal line
    Px = Positions(:,1);
    Py = C1(2)*ones(size(Positions,1),1);
    DistToOrigin = [0,Ly/2];
end
Norms = vecnorm([Px,Py]-Positions(:,1:2),2,2);
Bnds = find(Norms<=MultFactor*norm(DistToOrigin));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SizeTheDomain(Np,N_Kuhn,b)

global Lx Ly Lz PairsPerEdge 
PairsPerEdge = round(Np^(1/3));

Lx = 2*N_Kuhn*b;
Ly = Lx;
Lz = Lx;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitializeNodes

global N Lx Ly Lz Positions Corners dims

InitiateGrid; 

N = size(Positions,1);

% end 

Positions(:,dims+1) = (1:N)';   %Index No.
Positions(:,dims+2) = 1;        %Molecule number

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

CheckCorners = 0;
ShowFigDebug = 0;
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
        view([-30,30])
    end
    xlim([-Lx/2 Lx/2])
    ylim([-Ly/2 Ly/2])
    daspect([1 1 1])
    pbaspect([1 1 1])
    title(['N = ',num2str(size(Positions,1))])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
    close
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitiateGrid

global Positions Lx Ly Lz N_Kuhn b

Positions = [0 0 0;
    0.95*N_Kuhn*b 0 0];

CheckFig = 0;
if CheckFig==1
    figure(1)
    clf
    s = scatter3(Positions(:,1),Positions(:,2),Positions(:,3),'k','filled');
    xlim([-Lx/2 Lx/2])
    ylim([-Ly/2 Ly/2])
    zlim([-Lz/2 Lz/2])
    daspect([1 1 1])
    close
end
s.SizeData = 6;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotNetwork

global Lx Ly Lz ConnectionsSP Positions...
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
    froms = ConnectionsSP(1:N,1:end-1);
    froms = froms(:);
    tos = repmat(ConnectionsSP(:,end),size(ConnectionsSP,2)-1,1);
    Pairs = [froms tos]; 
    Pairs = sort(Pairs,2);
    Pairs = unique(Pairs,'rows');
    Pairs(Pairs(:,1)==0,:) = [];

    X1 = Positions(Pairs(:,1),1); 
    X2 = Positions(Pairs(:,2),1); 
    Y1 = Positions(Pairs(:,1),2); 
    Y2 = Positions(Pairs(:,2),2); 
    Z1 = Positions(Pairs(:,1),3); 
    Z2 = Positions(Pairs(:,2),3); 

    Color = 'b';
    LineWidth = 0.5;
    if dims==2
        plot([X1,X2]',[Y1,Y2]',Color,'LineWidth',1)
    else
        plot3([X1,X2]',[Y1,Y2]',[Z1,Z2]',Color,'LineWidth',LineWidth)
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
    s.SizeData = 8;
end
if PlotStickers==1
    s = scatter3(Positions(N+1:end,1),Positions(N+1:end,2),Positions(N+1:end,3),'r','filled');
    s.SizeData = 2;
end

grid on
% axis off

set(gcf,'Color','w')
daspect([1 1 1])
pbaspect([1 1 1])
if dims==3
    view([-30,30])
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
