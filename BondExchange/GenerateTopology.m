% Branched newtork topolgoy generator for input to LAMMPS via read_data cmd
% Robert J. Wagner
% 2022-11-29

% This version generates the x-links first and then connects

function GenerateTopology(Package,OverrideTop,OverrideIn,...
    CF,OF,CDL,TD,UnitType,LC,DC,BSOM,A_tab,LOL)

global dims ShowFiguresDebug Np...
    Ns N OldOrNew PlotStickers MakeMovie...
    TurnOnDynamics TopologyFileName BeadSpringOrMeso...
    SIOrNormalized InputFileName ConstantsFileName...
    LengthConversion DamperConversion CurrentFolder OutputFolder...
    current_dir_lmp dtFact LinearOrLangevin

%% Set Basic Toggles
OldOrNew = 1;           %Different looping options for two diff. connectivity algorithms
ShowFiguresDebug = 1;   %1 to plot intermediate figures for debugging
TurnOnDynamics = TD;
SIOrNormalized = UnitType;
LengthConversion = LC;
DamperConversion = DC;
BeadSpringOrMeso = BSOM;
CurrentFolder = CF;
OutputFolder = OF;
current_dir_lmp = CDL;
LinearOrLangevin = LOL;

N_Kuhn_out = [];

%% Unpack Swept Input Parameters
Sample = Package(1);    %Sample index
Np = Package(2);        %Number of molecules
Ns = Package(3);        %Number of stickeres per tether site
N = Np*Ns;              %Total numbber of permanent crosslinks (i.e., tether sites)
ka = Package(4);        %Activation energy of association
kd = Package(5);        %Activation energy of dissociation
f0 = Package(6);       %Force sensitivity to dissociation
N_Kuhn = Package(7);
b = Package(8);
Separation = Package(9);
phi = Package(10);
kbT = Package(11);

[damp,D,tau0,dtFact] = DefineTimeScale(b,LengthConversion,DamperConversion,...
    BeadSpringOrMeso);

dt = tau0/dtFact;

eaStar = -log(ka*(b*LengthConversion)^2/D);
edStar = -log(kd*(b*LengthConversion)^2/D);

A = A_tab(intersect(find(A_tab(:,1)==N_Kuhn),find(A_tab(:,2)==eaStar)),3);

%% Callout Input Script
% Initializes non-sweeping parameters
dims = 3;   %Dimensionality of system
InputScript(Sample,Np,Ns,N,N_Kuhn,b,phi,eaStar,edStar,D,dt,kbT);

%% Make diriectories and set filenames
SetDirAndFileNames;

%% Convert all to correct units
kb =1;
T = 1;
Tab = table(Np,Ns,N_Kuhn,b,ka,kd,f0,kb,T,dt,damp);
writetable(Tab,ConstantsFileName);

%%%%%%%%%%%%%%%%%%%%%%%% CHANGE LOOPING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
if ~isfile(TopologyFileName) || OverrideTop==1

    %% Size the domain
    % Determines correct domain size based on no. xls and nominal density
    SizeTheDomain(Np,Separation);

    %% Initialize the Nodes
    InitializeNodes(Separation);

    if BeadSpringOrMeso==0
        %% Polymerize Chains
        PolymerizeTheChains;
    elseif BeadSpringOrMeso==1
        %% Add tethered chains
        AddTetheredChains;
    end

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
    PrintInputFiles(Separation,N_Kuhn_out,b_out,ka_out,kd_out,f0_out,T_out,...
        dt_out,damp,tau0,A,CurrentFolder);
end

% toc
close
% clear
% clc
clear global

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PolymerizeTheChains

global dims N Ns Positions ConnectionsSP...
    X_SP Y_SP Z_SP Rx_SP Ry_SP Rz_SP N_Kuhn b...
    StationaryStableNodes DiffusingDynamicNodes DiffusingStableNodes

%% ONLY WORKS IF NS=1. SET UP A LOOP OVER ADD RANGE IF WANT TO MAKE IT GENERAL

N_mers = N*N_Kuhn;
AddRange = (N+1:N+N_mers)';
Positions(AddRange,:) = 0;
Positions(AddRange,dims+1) = AddRange; %Add xl numbers

MaxConn = 2;
ConnectionsSP = zeros(size(Positions,1),MaxConn+1);
ConnectionsSP(:,end) = (1:size(ConnectionsSP,1))';

StationaryStableNodes = (1:N)';
DiffusingDynamicNodes = zeros(Ns*N,1);
DiffusingStableNodes = zeros(length(AddRange)-length(DiffusingDynamicNodes),1);

ct = 0;
stabct = 0; dynct = 0;
for TetherNode=1:N
    for Segment=1:N_Kuhn
        ct = ct+1;
        NewNode = AddRange(ct);
        theta = 2*pi*rand(1);
        varphi = 2*pi*rand(1);
        radius = b;

        rad_proj = radius*cos(varphi);
        dx = rad_proj*cos(theta);
        dy = rad_proj*sin(theta);
        dz = radius*sin(varphi);

        if Segment==1
            OldIndx = find(Positions(:,end-1)==TetherNode,1,'first');
            OldPoint  = Positions(OldIndx,1:dims);
            OldNode = TetherNode;
            stabct = stabct+1;
            DiffusingStableNodes(stabct) = NewNode;
        elseif Segment==N_Kuhn
            OldPoint = Positions(NewNode-1,1:dims);
            OldNode = NewNode-1;
            dynct = dynct+1;
            DiffusingDynamicNodes(dynct) = NewNode;
        else
            OldPoint = Positions(NewNode-1,1:dims);
            OldNode = NewNode-1;
            stabct = stabct+1;
            DiffusingStableNodes(stabct) = NewNode;
        end

        Positions(NewNode,1:dims) = OldPoint + [dx dy dz];
        Positions(NewNode,dims+1) = NewNode;
        Positions(NewNode,dims+2) = TetherNode;

        ConnectionsSP(OldNode,2) = NewNode;
        ConnectionsSP(NewNode,1) = OldNode;
    end
end
% 
% theta = 2*pi*rand(N_mers,1);
% % varphimin = 0;
% varphimax = 2*pi;
% varphi = varphimax*rand(N_mers,1);
% radius = b;%normrnd(0,N_Kuhn*b^2,N_mers,1);
% OldPoints = Positions(1:N,1:dims);
% if dims==2
%     dx = radius.*cos(theta);
%     dy = radius.*sin(theta);
%     NewPoints = OldPoints + [dx dy];
% else
%     radius_proj = radius.*cos(varphi);
%     dx = radius_proj.*cos(theta);
%     dy = radius_proj.*sin(theta);
%     dz = radius.*sin(varphi);
%     NewPoints = OldPoints + [dx dy dz];
% end
% Positions(AddRange,1:dims) = NewPoints;

CheckFig = 1;
if CheckFig==1
    figure(1)
    hold on
    s=scatter3(Positions(1:N,1),Positions(1:N,2),Positions(1:N,3),'k','filled');
    s.SizeData = 6;
    s=scatter3(Positions(DiffusingStableNodes,1),Positions(DiffusingStableNodes,2),Positions(DiffusingStableNodes,3),'b','filled');
    s.SizeData=6;
    s=scatter3(Positions(DiffusingDynamicNodes,1),Positions(DiffusingDynamicNodes,2),Positions(DiffusingDynamicNodes,3),'r','filled');
    s.SizeData=6;
    close
end

% % ONLY WORKS IF NS=1
% ConnectionsSP = zeros(N+Ns*N,Ns+1);
% ConnectionsSP(:,end) = (1:size(ConnectionsSP,1))';
% ConnectionsSP(1:N,1) = ConnectionsSP(N+1:end,2);
% ConnectionsSP(N+1:end,1) = ConnectionsSP(1:N,2);    %Reciprocate

X_SP = zeros(size(ConnectionsSP));
Y_SP = zeros(size(ConnectionsSP));
Z_SP = zeros(size(ConnectionsSP));
for i=1:size(X_SP,1)
    for j=1:size(X_SP,2)
        if ConnectionsSP(i,j)~=0
            X_SP(i,j) = Positions(ConnectionsSP(i,j),1);
            Y_SP(i,j) = Positions(ConnectionsSP(i,j),2);
            Z_SP(i,j) = Positions(ConnectionsSP(i,j),3);
        else
            X_SP(i,j) = NaN;
            Y_SP(i,j) = NaN;
            Z_SP(i,j) = NaN;
        end
    end
end

Rx_SP = zeros(size(Positions,1),MaxConn);
Ry_SP = zeros(size(Positions,1),MaxConn);
Rz_SP = zeros(size(Positions,1),MaxConn);
for i=1:MaxConn
    Rx_SP(:,i) = -X_SP(:,i)+X_SP(:,end);
    Ry_SP(:,i) = -Y_SP(:,i)+Y_SP(:,end);
    Rz_SP(:,i) = -Z_SP(:,i)+Z_SP(:,end);
end
X_SP(:,end) = [];
Y_SP(:,end) = [];
Z_SP(:,end) = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PrintInputFiles(dist,N,b,ka0,kd0,f0,T,dt,damp,tau0,A,CurrentFolder)

global Top_file InputFileName...
    TurnOnDynamics SIOrNormalized BeadSpringOrMeso dtFact LinearOrLangevin...
    OutputDataFolder OutputAtom OutputBond

TurnOnLJ = 0;       % 0 for no volume exclusion (do not need repulsive interactions
if BeadSpringOrMeso==1 % Mesoscale
    damp_eff = damp*N^(2/3);    % Effective damper needed to roughly mimic bead-spring diffusion
    Ldyn = b;               % Attachment cutoff length - smaller for mesoscale
elseif BeadSpringOrMeso==0 % Bead-spring
    damp_eff = damp;
    Ldyn = b;               % Attachment cutoff length - smaller for mesoscale
elseif BeadSpringOrMeso==2 % Scaling Theory
    damp_eff = damp*N^(2/3);    % Effective damper needed to roughly mimic bead-spring diffusion
    Ldyn = N*b;               % Attachment cutoff length - smaller for mesoscale    
end

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
kT = 1;
ea = -kT*log(ka0*tau0);
if BeadSpringOrMeso==2
    N_temp = 2*N;
    ka0 = A/tau0*(3/2/pi/N_temp)^(3/2)*exp(-ea);
    l0 = (2*N_temp*b^2/3)^(0.5);
end

% Define lengthscales
Nb = N*b;           %chain length
Nb2 = Nb*b;         %N*b^2
r0 = sqrt(N)*b;     %mean random walk chain length
Kstiff = 3*kb*T/Nb2;     %chain stiffness (1/2 multiple to capture effective damper)
ftol = Kstiff*b/10000;  %equilibration force tolerance
if BeadSpringOrMeso==0 
    Lmax = 1.1*b;%0.001*36*b;%Nb;%0.4*Nb;       %max LJ cutoff - set very small for no LJs
elseif BeadSpringOrMeso==1 && N<36
    Lmax = 10*b;	%max LJ cutoff - set very small for no LJs
else
    Lmax = 2;
end
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# DEFINE LENGTH & FORCE SCALES'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Defining length scales'; Assembled = [Assembled,L];
L = ['\nvariable    b       equal ',num2str(b,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    N       equal ',num2str(N,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    Nb      equal ',num2str(Nb,'%.2e')]; Assembled = [Assembled,L];
% L = ['\nvariable    cut_dyn equal ',num2str(0.95*Nb,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    cut_dyn equal ',num2str(Nb,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    Nb2     equal ',num2str(Nb2,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    r0      equal ',num2str(r0,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    Kstiff  equal ',num2str(Kstiff,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    ftol    equal ',num2str(ftol,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    Lmax    equal ',num2str(Lmax,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    ka0     equal ',num2str(ka0,'%.2e')]; Assembled = [Assembled,L];
if BeadSpringOrMeso==2
    L = ['\nvariable    l0      equal ',num2str(l0,'%.2e')]; Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% % L = '\n# Defining max stretch'; Assembled = [Assembled,L];
% % L = ['\nvariable    lammax  equal ',num2str(lambda,'%.2f')]; Assembled = [Assembled,L];
% % L = '\n'; Assembled = [Assembled,L];

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
sig_b = 2*b;                % backbone LJ length scale
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
L = ['\nvariable    sig_b   equal ',num2str(sig_b,'%.2e')]; Assembled = [Assembled,L];
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
theal = 2e3*dt*dtFact;
tload = 0;
trelax = 0;
nheal = ceil(theal/dt);
nload = ceil(tload/dt);
nrelax = ceil(trelax/dt);
nmax = nheal+nload+nrelax;
ttot = nmax*dt;

N_extra_bonds = 2;
N_extra_special = 1000;
if BeadSpringOrMeso==2
%     N_extra_neighb = 3e3;
    N_extra_neighb = 10*ceil(4.35*(dist/N/b)^(-3));
end

L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# TIMESCALES & LOADING RATES'; Assembled = [Assembled,L];
L = '\n##########################`##########'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% % L = '\n# Weissenberg number'; Assembled = [Assembled,L];
% % L = ['\nvariable    W       equal ',num2str(W,'%.1f')]; Assembled = [Assembled,L];
% % L = ['\nvariable    lamdot  equal ',num2str(lamdot,'%.2e')]; Assembled = [Assembled,L];
% % L = '\n'; Assembled = [Assembled,L];

L = '\n# Timescales'; Assembled = [Assembled,L];
L = ['\nvariable    taukd   equal ',num2str(taukd,'%.1f')]; Assembled = [Assembled,L];
L = ['\nvariable    tauka   equal ',num2str(tauka,'%.1f')]; Assembled = [Assembled,L];
% L = ['\nvariable    tauW    equal ',num2str(tauW,'%.1f')]; Assembled = [Assembled,L];
L = ['\nvariable    tau_heal equal ',num2str(tau_heal,'%.1f')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Timestep'; Assembled = [Assembled,L];
L = ['\nvariable    dt      equal ',num2str(dt,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Time for each portion of the simulation'; Assembled = [Assembled,L];
L = ['\nvariable    theal   equal ',num2str(theal,'%.1f')]; Assembled = [Assembled,L];
% L = ['\nvariable    tload   equal ',num2str(tload,'%.1f')]; Assembled = [Assembled,L];
% L = ['\nvariable    trelax  equal ',num2str(trelax,'%.1f')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Number of timesteps for each sequence'; Assembled = [Assembled,L];
L = ['\nvariable    nheal   equal ',num2str(nheal,'%.0f')]; Assembled = [Assembled,L];
% L = ['\nvariable    nload   equal ',num2str(nload,'%.0f')]; Assembled = [Assembled,L];
% L = ['\nvariable    nrelax  equal ',num2str(nrelax,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    nmax    equal ',num2str(nmax,'%.0f')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Running parameters
% % % NData = 500; 
% % if BeadSpringOrMeso==0 || BeadSpringOrMeso==1
% %     N_dyn = 1;
% % elseif BeadSpringOrMeso==2
% %     N_dyn = 1;%dtFact;
% % end
% % tdata = dt*dtFact/10*N_dyn;  %Basically measures with frequency matching 1/10x the diffusion timescale
% % ttot = nmax*dt;
% % NData = ceil(ttot/tdata);
% % iout = ceil(nmax/NData); 
% % ithermo = ceil(iout/2);
% Define trial frequency for dynamics

% Sampling timescales
tdyn = tau0/20;
N_dyn = round(tdyn/dt);
disp(['Dynamics will be sampled every ',num2str(N_dyn),' steps'])

% Define output frequency
tdata = tau0/10;
N_data = round(ttot/tdata);       % Number of outputs data points - target ~1k
disp(['The total number of steps will be ',num2str(nmax)])
disp(['The total number of output data points will be ',num2str(N_data)])
iout = round(tdata/dt);

tthermo = tau0/20;
% N_thermo = round(ttot/tthermo);
ithermo = round(tthermo/dt); 
disp(['and the cmd prompt will provide read out every ',num2str(ithermo),' steps.'])


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
    L = '\nneighbor     ${Nb}   bin'; Assembled = [Assembled,L];
    L = '\ncomm_modify  cutoff  ${Lmax}'; Assembled = [Assembled,L];
    L = ['\nneigh_modify   one     ',num2str(N_extra_neighb,'%.0f'),...
        '   page    ',num2str(10*N_extra_neighb,'%.0f')]; Assembled = [Assembled,L];
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

% L = '\n# Define bond type'; Assembled = [Assembled,L];
% L = '\nbond_style    harmonic'; Assembled = [Assembled,L];
% L = '\n'; Assembled = [Assembled,L];

% Define bond styles
if BeadSpringOrMeso==0 || (BeadSpringOrMeso==1 && LinearOrLangevin==0)
    L = '\n# Define bond type'; Assembled = [Assembled,L];
    if BeadSpringOrMeso==1
        L = '\nbond_style    hybrid  harmonic nonlinear'; Assembled = [Assembled,L];
    elseif BeadSpringOrMeso==2
        L = '\nbond_style    harmonic'; Assembled = [Assembled,L];
    elseif BeadSpringOrMeso==0
        L = '\nbond_style    nonlinear'; Assembled = [Assembled,L];
    end
    L = '\n'; Assembled = [Assembled,L];
else
    L = '\nbond_style    hybrid pade nonlinear'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
end

% InputFolder1 = ['/mnt/c/Users/rjwag/Documents/Research/NSF/',CurrentFolder,'&'];
% InputFolder2 = ['/',TopologyFileName];
% L = '\n# Import the topology'; Assembled = [Assembled,L];
% L = ['\nread_data ',InputFolder1]; Assembled = [Assembled,L];
% L = ['\n',InputFolder2,'&']; Assembled = [Assembled,L];

InputFolder1 = ['Input_Data/Initial_Topologies/','&'];
InputFolder2 = Top_file;
L = '\n# Import the topology'; Assembled = [Assembled,L];
L = ['\nread_data ',InputFolder1]; Assembled = [Assembled,L];
L = ['\n',InputFolder2,'&']; Assembled = [Assembled,L];

L = ['\n extra/bond/per/atom ',num2str(N_extra_bonds,'%.0f'),...
    ' extra/special/per/atom ',num2str(N_extra_special,'%.0f')]; 
Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Define repulsive potentials between each atom'; Assembled = [Assembled,L];
if BeadSpringOrMeso==0 || BeadSpringOrMeso==1
L = '\npair_style  lj/cut/soft $(2.0) $(1.0) $(v_r_c)'; Assembled = [Assembled,L];
L = '\npair_coeff  1 1	${eps_r0} ${sig_r0} ${lambda_r0}'; Assembled = [Assembled,L];
L = '\npair_coeff  2 2	${eps_b}  ${sig_b}  ${lambda_b}'; Assembled = [Assembled,L];
elseif BeadSpringOrMeso==2
    L = '\npair_style  lj/cut/soft $(2.0) $(1.0) ${Nb}'; Assembled = [Assembled,L];
    L = '\npair_coeff  1 1	${eps_r0} ${sig_r0} ${Nb}'; Assembled = [Assembled,L];
end
if BeadSpringOrMeso==0 % Bead-spring
    L = '\npair_coeff  3 3	${eps_b}  ${sig_b}  ${lambda_b}'; Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

L = '\n# Label the groups'; Assembled = [Assembled,L];
% Generate random seed
rand_seed = randi(9999);
if BeadSpringOrMeso==1 % Mesoscale
    L = '\ngroup	stationarystable type 1'; Assembled = [Assembled,L];
    L = '\ngroup	diffusingdynamic type 2'; Assembled = [Assembled,L];
    L = '\ngroup	stickers type 2'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
    L = '\nlabelmap	atom	1 bckbn		2 stckr'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    L = '\n# Set Brownian integration'; Assembled = [Assembled,L];
    L = ['\nfix      free diffusingdynamic brownian ${T} ',num2str(rand_seed),' gamma_t ${damp}']; Assembled = [Assembled,L];
elseif BeadSpringOrMeso==0 % Bead-spring
    L = '\ngroup	stationarystable type 1'; Assembled = [Assembled,L];
    L = '\ngroup	diffusingstable type 2'; Assembled = [Assembled,L];
    L = '\ngroup	diffusingdynamic type 3'; Assembled = [Assembled,L];
    L = '\ngroup	stickers type 3'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
    L = '\nlabelmap	atom	1 bckbn		3 stckr'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    L = '\n# Set Brownian integration'; Assembled = [Assembled,L];
    L = ['\nfix      free1 diffusingstable brownian ${T} ',num2str(rand_seed),' gamma_t ${damp}']; Assembled = [Assembled,L];
    L = ['\nfix      free2 diffusingdynamic brownian ${T} ',num2str(rand_seed),' gamma_t ${damp}']; Assembled = [Assembled,L];
elseif BeadSpringOrMeso==2 % Scaling theory
    L = '\ngroup	stationarydynamic type 1'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
    L = '\nlabelmap	atom	1 stckr'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% Initial Equilibration
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# REACH STEADY STATE'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% NO DYNAMICS FOR THIS STUDY
if TurnOnDynamics==1
    if BeadSpringOrMeso==1 % Mesoscale
        L = ['\nfix		dynamic stickers bond/dynamic ',num2str(N_dyn,'%.0f'),...
            ' 2 2 ${ka} ${kd} ${Ldyn} maxbond 1'];% bell ${fsens}'];
        Assembled = [Assembled,L];
    elseif BeadSpringOrMeso==0 % Bead-spring
        L = ['\nfix		dynamic stickers bond/dynamic ',num2str(N_dyn,'%.0f'),...
            ' 3 2 ${ka} ${kd} ${Ldyn} maxbond 1'];% bell ${fsens}'];
        Assembled = [Assembled,L];    
    elseif BeadSpringOrMeso==2 % Scaling Theory - can connect span of chain length
        L = ['\nfix		dynamic stationarydynamic bond/dynamic ',num2str(N_dyn,'%.0f'),...
            ' 1 1 ${ka} ${kd} ${Ldyn} maxbond 1 gauss ${l0} ${ka0}'];
        Assembled = [Assembled,L];   
    end
end
% L = '\nfix		integrate all nve'; Assembled = [Assembled,L];
% L = '\n'; Assembled = [Assembled,L];

L = '\n# Compute bond information'; Assembled = [Assembled,L];
L = '\ncompute 	Pairs all property/local btype batom1 batom2'; Assembled = [Assembled,L];
L = '\ncompute 	Bonds all bond/local dx dy dz dist force'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% OutputAtom_loc = [OutputDataFolder,OutputAtom];
% OutputBond_loc = [OutputDataFolder,OutputBond];
OutputAtom_loc = OutputAtom;
OutputBond_loc = OutputBond;

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
% % L = '\nmin_style cg'; Assembled = [Assembled,L];
% % L = '\nminimize 0.0 ${ftol} 10000 100000'; Assembled = [Assembled,L];
% % L = '\n'; Assembled = [Assembled,L];

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

% L = ['\nfix		printS all print ${iout} "${teq} ${S1} ${S2} ${S3}" file ',...
%     StressEqlFileName,' screen no']; Assembled = [Assembled,L];
% L = '\n'; Assembled = [Assembled,L];

L = '\n# Run initial equilibration'; Assembled = [Assembled,L];
L = '\nrun  ${nheal}'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% NO LOADING OR RELAXATION FOR THIS STUDY
% % % Apply deformation
% % L = '\n####################################'; Assembled = [Assembled,L];
% % L = '\n# APPLY DEFORMATION'; Assembled = [Assembled,L];
% % L = '\n####################################'; Assembled = [Assembled,L];
% % L = '\n'; Assembled = [Assembled,L];
% % L = ['\nfix		load_uniaxial all deform 1 x trate ${lamdot} y trate ',...
% %     '-$(v_lamdot*0.5) z trate -$(v_lamdot*0.5)']; Assembled = [Assembled,L];
% % L = '\n'; Assembled = [Assembled,L];
% % 
% % L = '\n# Output deformation data'; Assembled = [Assembled,L];
% % L = '\n'; Assembled = [Assembled,L];
% % L = '\nvariable	lam equal exp(${lamdot}*${dt}*elapsed)'; Assembled = [Assembled,L];
% % L = '\nvariable	tl equal ${dt}*elapsed'; Assembled = [Assembled,L];
% % L = '\n'; Assembled = [Assembled,L];
% % L = ['\nfix		printS all print ${iout} "${tl} ${lam} ${S1} ${S2} ${S3}" file ',...
% %     StressLdgFileName,' screen no']; Assembled = [Assembled,L];
% % L = '\n'; Assembled = [Assembled,L];
% % % 
% % % L = ['\nfix		printS all print ${iout} "${tl} ${id} ${type} ${x} ${y} ${z}" file ',...
% % %     AtomsLdgFileName,' screen no']; Assembled = [Assembled,L];
% % % L = '\n'; Assembled = [Assembled,L];
% % % 
% % % L = ['\nfix		printS all print ${iout} "${tl} ${index} ${p1} ${p2}" file ',...
% % %     BondsLdgFileName,' screen no']; Assembled = [Assembled,L];
% % % L = '\n'; Assembled = [Assembled,L];
% % 
% % L = '\n# Run applied deformation'; Assembled = [Assembled,L];
% % L = '\nrun  ${nload}'; Assembled = [Assembled,L];
% % L = '\n'; Assembled = [Assembled,L];
% % 
% % % Relax system
% % L = '\n####################################'; Assembled = [Assembled,L];
% % L = '\n# HOLD AT CONSTANT DEFORM. - RELAX'; Assembled = [Assembled,L];
% % L = '\n####################################'; Assembled = [Assembled,L];
% % L = '\n'; Assembled = [Assembled,L];
% % L = '\nunfix	load_uniaxial'; Assembled = [Assembled,L];
% % L = '\n'; Assembled = [Assembled,L];
% % L = '\n# Output relaxation data'; Assembled = [Assembled,L];
% % L = '\nvariable	tr equal ${dt}*elapsed'; Assembled = [Assembled,L];
% % L = '\n'; Assembled = [Assembled,L];
% % L = ['\nfix		printS all print ${iout} "${tr} ${lam} ${S1} ${S2} ${S3}" file '...
% %     StressRlxFileName,' screen no']; Assembled = [Assembled,L];
% % L = '\n'; Assembled = [Assembled,L];
% % % 
% % % L = ['\nfix		printS all print ${iout} "${tr} ${id} ${type} ${x} ${y} ${z}" file ',...
% % %     AtomsRlxFileName,' screen no']; Assembled = [Assembled,L];
% % % L = '\n'; Assembled = [Assembled,L];
% % % 
% % % L = ['\nfix		printS all print ${iout} "${tr} ${index} ${p1} ${p2}" file ',...
% % %     BondsRlxFileName,' screen no']; Assembled = [Assembled,L];
% % % L = '\n'; Assembled = [Assembled,L];
% % 
% % L = '\n# Run relaxation'; Assembled = [Assembled,L];
% % L = '\nrun		${nrelax}'; Assembled = [Assembled,L];

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
    N_atoms N_bonds TopologyFileName Ns BeadSpringOrMeso...
    StationaryStableNodes DiffusingStableNodes...
    DiffusingDynamicNodes LinearOrLangevin

kbT = kb*T;

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
if BeadSpringOrMeso==0 % Bead-spring
    N_atom_types = 3;   % stable crosslinks, diffusing dynamic, stationary dynamic
    N_bond_types = 2;   % Short, dynamics vs long implicit chains
elseif BeadSpringOrMeso==1 % Mesoscale
    N_atom_types = 2;   % Stationary Kuhn bead, diffusing Kuhn bead, diffusing dynamic
    N_bond_types = 2;   % Kuhn segments and dynamic bonds
elseif BeadSpringOrMeso==2 % Scaling theory
    N_atom_types = 1;   % Stationary Kuhn bead, diffusing Kuhn bead, diffusing dynamic
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
weight = 0.5;
mass_sticker = (1-weight)*N_Kuhn*mass_mer;
mass_backbone = N_Kuhn*mass_mer + 1/2*Ns*weight*N_Kuhn*mass_mer;
L = '\nMasses';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
for i=1:N_atom_types
    type = i;
    if BeadSpringOrMeso==1 % Mesoscale
        if type==1
            mass = mass_backbone;
        elseif type==2
            mass = mass_sticker;
        end
    elseif BeadSpringOrMeso==0 % Bead-spring
        if type==1 || type==2
            mass = mass_backbone;
        elseif type==3 
            mass = mass_sticker;
        end
    elseif BeadSpringOrMeso==2 % Scaling theory
        mass = mass_sticker;
    end
    L = ['\n',num2str(type),' ',num2str(mass)];
    Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% Append bond coefficients - for Langevin spring potentials
L = '\nBond Coeffs';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
soften_coeff = 1.02;
for i=1:N_bond_types
    type = i;
    if BeadSpringOrMeso==1 % Mesoscale
        if LinearOrLangevin==0
            if type==1
                coeff1 = 3/2*kbT/(N_Kuhn*b^2);
                coeff2 = 0;
                L = ['\n',num2str(type),' harmonic ',...
                    num2str(coeff1),' ',num2str(coeff2)];
            elseif type==2
                coeff1 = 100*kbT;
                coeff2 = b;
                coeff3 = soften_coeff*b;
                L = ['\n',num2str(type),' nonlinear ',...
                    num2str(coeff1),' ',num2str(coeff2),' ',num2str(coeff3)];
            end
        elseif LinearOrLangevin==1
            if type==1
                coeff1 = b;
                coeff2 = N_Kuhn;
                L = ['\n',num2str(type),' pade ',...
                    num2str(coeff1),' ',num2str(coeff2)];
            elseif type==2
                coeff1 = 100*kbT;
                coeff2 = b;
                coeff3 = soften_coeff*b;
                L = ['\n',num2str(type),' nonlinear ',...
                    num2str(coeff1),' ',num2str(coeff2),' ',num2str(coeff3)];
            end
        end
    elseif BeadSpringOrMeso==0 % Bead-spring
        if type==1
            coeff1 = 800*kbT;
            coeff2 = b;
            coeff3 = b;
        elseif type==2
            coeff1 = 100*kbT;
            coeff2 = b;
            coeff3 = soften_coeff*b;
        end
        L = ['\n',num2str(type),' ',num2str(coeff1),...
            ' ',num2str(coeff2),...
            ' ',num2str(coeff3)];
    elseif BeadSpringOrMeso==2 % Scaling Theory
        coeff1 = 3/2*kbT/(2*N_Kuhn*b^2);
        coeff2 = 0;
        L = ['\n',num2str(type),' ',...
            num2str(coeff1),' ',num2str(coeff2)];
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
        if ismember(i,StationaryStableNodes)
            atom_type = 1; %Stationary stable bonds
        elseif ismember(i,DiffusingDynamicNodes)
            atom_type = 2; %Diffusing dynamic bonds                                 
        end

        L = ['\n',num2str(tag),' ',num2str(mol_tag),' ',num2str(atom_type),' ',...
            num2str(x),' ',num2str(y),' ',num2str(z)];
        Assembled = [Assembled,L];
    end
elseif BeadSpringOrMeso==0 % Bead-spring
    for i=1:N_atoms
        tag = Positions(i,dims+1);
        mol_tag = Positions(i,dims+2);
        x = Positions(i,1); y = Positions(i,2);
        if dims==2
            z = 0;
        else
            z = Positions(i,3);
        end

        if ismember(i,StationaryStableNodes)
            atom_type = 1;
        elseif ismember(i,DiffusingDynamicNodes)
            atom_type = 3;
        elseif ismember(i,DiffusingStableNodes)
            atom_type = 2;
        end

        L = ['\n',num2str(tag),' ',num2str(mol_tag),' ',num2str(atom_type),' ',...
            num2str(x),' ',num2str(y),' ',num2str(z)];
        Assembled = [Assembled,L];
    end
elseif BeadSpringOrMeso==2 % Theory
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
%     Backbones = Positions(1:N,dims+1);
%     Stickers = Positions(N+1:end,dims+1);
%     for i=1:N_bonds
%         tag = i;
%         atom_1 = UniquePairs(i,1); atom_2 = UniquePairs(i,2);
%         if (ismember(atom_1,Backbones) && ismember(atom_2,Stickers)) ||...
%                 (ismember(atom_2,Backbones) && ismember(atom_1,Stickers))
%             bond_type = 1;
%         else
%             bond_type = 2;
%         end
%         L = ['\n',num2str(tag),' ',num2str(bond_type),' ',...
%             num2str(atom_1),' ',num2str(atom_2)];
%         Assembled = [Assembled,L];
%     end
% else
    for i=1:N_bonds
        tag = i;
        atom_1 = UniquePairs(i,1); atom_2 = UniquePairs(i,2);
        bond_type = 1;
        L = ['\n',num2str(tag),' ',num2str(bond_type),' ',...
            num2str(atom_1),' ',num2str(atom_2)];
        Assembled = [Assembled,L];
    end
% end
end
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
function AddTetheredChains

global dims N Positions Ns ConnectionsSP...
    X_SP Y_SP Z_SP Rx_SP Ry_SP Rz_SP N_Kuhn b...
    StationaryStableNodes DiffusingDynamicNodes...
    Lx Ly Lz


%% ONLY WORKS IF NS=1. SET UP A LOOP OVER ADD RANGE IF WANT TO MAKE IT GENERAL

N_stickers = N*Ns;
AddRange = (N+1:N+N_stickers)';
Positions(AddRange,:) = 0;
Positions(AddRange,dims+1) = AddRange; %Add xl numbers

% Rng = (2:2:N)';
% Rng(2:2:end) = [];
% StationaryStableNodes = sort([Rng;Rng-1]);
% StationaryDynamicNodes = (1:N)';
% StationaryDynamicNodes(ismember(StationaryDynamicNodes,StationaryStableNodes)) = [];

StationaryStableNodes = (1:N)';
DiffusingDynamicNodes = AddRange;

theta = 2*pi*rand(N_stickers,1);
% varphimin = 0;
varphimax = 2*pi;
varphi = varphimax*rand(N_stickers,1);
radius = normrnd(0,N_Kuhn*b^2,N_stickers,1);
OldPoints = Positions(StationaryStableNodes,1:dims);
if dims==2
    dx = radius.*cos(theta);
    dy = radius.*sin(theta);
    NewPoints = OldPoints + [dx dy];
else
    radius_proj = radius.*cos(varphi);
    dx = radius_proj.*cos(theta);
    dy = radius_proj.*sin(theta);
    dz = radius.*sin(varphi);
    NewPoints = OldPoints + [dx dy dz];
end
Positions(AddRange,1:dims) = NewPoints;
Positions(AddRange,dims+1) = DiffusingDynamicNodes;
Positions(AddRange,dims+2) = StationaryStableNodes;
Positions(AddRange,end) = 2;

CheckFig = 1;
if CheckFig==1
    figure(1)
    hold on
    s=scatter3(Positions(StationaryStableNodes,1),...
        Positions(StationaryStableNodes,2),...
        Positions(StationaryStableNodes,3),'c','filled');
    s.SizeData = 6;

    s=scatter3(Positions(DiffusingDynamicNodes,1),...
        Positions(DiffusingDynamicNodes,2),...
        Positions(DiffusingDynamicNodes,3),'m','filled');
    s.SizeData = 6;

    view(60,30)
    xlim([-Lx/2 Lx/2])
    ylim([-Ly/2 Ly/2])
    zlim([-Lz/2 Lz/2])
    pbaspect([1 1 1])
    close
end

% ONLY WORKS IF NS=1
ConnectionsSP = zeros(N+Ns*N,Ns+1);
ConnectionsSP(:,end) = (1:size(ConnectionsSP,1))';
ConnectionsSP(StationaryStableNodes,1) = ConnectionsSP(DiffusingDynamicNodes,2);
ConnectionsSP(DiffusingDynamicNodes,1) = ConnectionsSP(StationaryStableNodes,2);    %Reciprocate

X_SP = zeros(size(ConnectionsSP));
Y_SP = zeros(size(ConnectionsSP));
Z_SP = zeros(size(ConnectionsSP));
for i=1:size(X_SP,1)
    for j=1:size(X_SP,2)
        if ConnectionsSP(i,j)~=0
            X_SP(i,j) = Positions(ConnectionsSP(i,j),1);
            Y_SP(i,j) = Positions(ConnectionsSP(i,j),2);
            Z_SP(i,j) = Positions(ConnectionsSP(i,j),3);
        end
    end
end
Rx_SP = -X_SP(:,1)+X_SP(:,2);
Ry_SP = -Y_SP(:,1)+Y_SP(:,2);
Rz_SP = -Z_SP(:,1)+Z_SP(:,2);
X_SP(:,end) = [];
Y_SP(:,end) = [];
Z_SP(:,end) = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SizeTheDomain(Np,Separation)

global Lx Ly Lz PairsPerEdge c0

PairsPerEdge = round(Np^(1/3));

Lx = PairsPerEdge*Separation;
Ly = Lx;
Lz = Lx;

c0 = Np/(Lx*Ly*Lz);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitializeNodes(Separation)

global N Lx Ly Lz Positions Corners dims

% InitiateGrid; 
InitiateGrid(Separation);

N = size(Positions,1);

% end 

Positions(:,end+1) = (1:N)';    %Index No.
Positions(:,end+1) = (1:N)';    %Molecule number
Positions(:,end+1) = 1;     % All ancho points are stable nodes

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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitiateGrid(Separation)

global Positions Lx Ly Lz

XMin = -Lx/2+Separation/2;
XMax = Lx/2-Separation/2;
XSpacing = (XMin:Separation:XMax)';
YSpacing = XSpacing;
ZSpacing = XSpacing;

[X,Y,Z] = meshgrid(XSpacing,YSpacing,ZSpacing);
X = X(:)-mean(X(:));
Y = Y(:)-mean(Y(:));
Z = Z(:)-mean(Z(:));
Positions = [X(:) Y(:) Z(:)];

CheckFig = 1;
if CheckFig==1
    figure(1)
    clf
    s = scatter3(X(:),Y(:),Z(:),'k','filled');
    xlim([-Lx/2 Lx/2])
    ylim([-Ly/2 Ly/2])
    zlim([-Lz/2 Lz/2])
    daspect([1 1 1])
end
s.SizeData = 6;
close

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotNetwork

global X_SP Y_SP Z_SP Rx_SP Ry_SP Rz_SP Lx Ly Lz ConnectionsSP Positions...
    Corners dims PlotChains Np N PlotStickers

hold on
FontSize = 20;
xl_per_mol = N/Np;

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
    Corners dims PlotChains Np N PlotStickers

clf
hold on
FontSize = 20;
xl_per_mol = N/Np;

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
