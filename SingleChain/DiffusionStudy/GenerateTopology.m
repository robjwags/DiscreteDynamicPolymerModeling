% Branched newtork topolgoy generator for input to LAMMPS via read_data cmd
% Robert J. Wagner
% 2022-11-29

% This version generates the x-links first and then connects

function GenerateTopology(Package,OverrideTop,OverrideIn,...
    CF,OF,TD,UnitType,LC,DC,BSOM,dtF)

global dims ShowFiguresDebug Np...
    Ns N OldOrNew PlotStickers MakeMovie...
    TurnOnDynamics TopologyFileName BeadSpringOrMeso...
    SIOrNormalized InputFileName ConstantsFileName...
    LengthConversion DamperConversion CurrentFolder OutputFolder dtFact

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
dtFact = dtF;

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

    %% Add chains
    if BeadSpringOrMeso==0
        % Polymerize Chains
        PolymerizeTheChains;
    else
        % Add tethered chains
        AddTetheredChains;
    end

    figure(1); clf; hold on
    PlotStickers=1;
    if ShowFiguresDebug==1
        PlotNetwork
    end

    %% Output Topological Data
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
function PolymerizeTheChains

global dims N Positions ConnectionsSP...
    X_SP Y_SP Z_SP Rx_SP Ry_SP Rz_SP N_Kuhn b

%% ONLY WORKS IF NS=1. SET UP A LOOP OVER ADD RANGE IF WANT TO MAKE IT GENERAL

N_mers = N*N_Kuhn;
AddRange = (N+1:N+N_mers)';
Positions(AddRange,:) = 0;
Positions(AddRange,dims+1) = AddRange; %Add xl numbers

MaxConn = 2;
ConnectionsSP = zeros(size(Positions,1),MaxConn+1);
ConnectionsSP(:,end) = (1:size(ConnectionsSP,1))';

ct = 0;
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
            OldIndx = find(Positions(:,end)==TetherNode,1,'first');
            OldPoint  = Positions(OldIndx,1:dims);
            OldNode = TetherNode;
        else
            OldPoint = Positions(NewNode-1,1:dims);
            OldNode = NewNode-1;
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
    s=scatter3(Positions(AddRange,1),Positions(AddRange,2),Positions(AddRange,3),'r','filled');
    s.SizeData=6;
    close all
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
function PrintInputFiles(N,b,ka0,kd0,f0,T,dt,damp,CurrentFolder)

global TopologyFileName OutputAtom OutputBond InputFileName...
    TurnOnDynamics SIOrNormalized OutputAtom_loc OutputBond_loc dtFact...
    BeadSpringOrMeso

TurnOnLJ = 0;

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
Lmax = 0.4*Nb;      %max LJ cutoff
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
if BeadSpringOrMeso==0
    L = ['\nvariable    damp    equal ',num2str(damp,'%.2e')]; Assembled = [Assembled,L];
elseif BeadSpringOrMeso==1
    damp_eff = damp*N^(2/3);    %To account for fact that half of chain provides friction
    L = ['\nvariable    damp    equal ',num2str(damp_eff,'%.2e')]; Assembled = [Assembled,L];
end
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
Ldyn = b;                 % Attachment cutoff length
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
tau0 = dt*dtFact;
teq = 60e3*tau0;
neq = ceil(teq/dt);
nmax = neq;

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
L = ['\nvariable    theal   equal ',num2str(teq,'%.1f')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Number of timesteps for each sequence'; Assembled = [Assembled,L];
L = ['\nvariable    neq   equal ',num2str(neq,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    nmax    equal ',num2str(nmax,'%.0f')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Running parameters
tdata = tau0;
ttot = nmax*dt;
n_data = ceil(ttot/tdata);
iout = 10*ceil(nmax/n_data); 
ithermo = 10*iout;
n_dyn = iout;   %check dynamics every N_dyn steps 
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
if BeadSpringOrMeso==0
    L = '\nbond_style    nonlinear'; Assembled = [Assembled,L];
elseif BeadSpringOrMeso==1
    L = '\nbond_style    harmonic'; Assembled = [Assembled,L];
end
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
L = '\npair_style  lj/cut/soft $(2.0) $(1.0) $(v_r_c)'; Assembled = [Assembled,L];
L = '\npair_coeff  1 1	${eps_r0} ${sig_r0} ${lambda_r0}'; Assembled = [Assembled,L];
L = '\npair_coeff  2 2	${eps_b}  ${sig_b}  ${lambda_b}'; Assembled = [Assembled,L];
if BeadSpringOrMeso==0
    L = '\npair_coeff  3 3	${eps_b}  ${sig_b}  ${lambda_b}'; Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

L = '\n# Label the groups'; Assembled = [Assembled,L];
L = '\ngroup	backbones type 1'; Assembled = [Assembled,L];
L = '\ngroup	stickers type 2'; Assembled = [Assembled,L];
if BeadSpringOrMeso==0
    L = '\ngroup	endgroups type 3'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
end
L = '\nlabelmap	atom	1 bckbn		2 stckr'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Set Brownian integration'; Assembled = [Assembled,L];
L = '\nfix      free1 stickers brownian ${T} 1234 gamma_t ${damp}'; Assembled = [Assembled,L];
if BeadSpringOrMeso==0
    L = '\nfix      free2 endgroups brownian ${T} 1234 gamma_t ${damp}'; Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% Initial Equilibration
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# REACH STEADY STATE'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% NO DYNAMICS FOR THIS STUDY
if TurnOnDynamics==1
    if BeadSpringOrMeso==0
        L = ['\nfix		dynamic endgroups bond/dynamic ',num2str(n_dyn,'%.0f'),...
            ' 3 2 ${ka} ${kd} ${Ldyn} maxbond 1'];% bell ${fsens}'];
    elseif BeadSpringOrMeso==1
        L = ['\nfix		dynamic stickers bond/dynamic ',num2str(n_dyn,'%.0f'),...
            ' 2 3 ${ka} ${kd} ${Ldyn} maxbond 1 bell ${fsens}'];
    end
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
L = '\nrun  ${neq}'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

fid = fopen(InputFileName,'wt');
fprintf(fid,Assembled);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b_out,ka_out,kd_out,f0_out,kb_out,T_out,dt_out,damp_out,m_out] =...
    ConvertUnits(b_code,ka_code,kd_code,f0_code,T_code,dt_code,damp_code)

global SIOrNormalized

% Plots surface plot of k_d versus tau_0 and \Delta G/kbT
% Actual units
kb = 1.38e-23;      % J/K
T_out = T_code;         % K
Na = 6.022e23;          % mol^-1
kbT = kb*T_out;     % Thermal energy, J
kd_out = kd_code;       % Hz
ka_out = ka_code;       % Hz
m = 75;                 % g/mol (approximate - EG is ~60, acrylamide is ~76, sytrene is ~100)

CheckTimescales(kbT,Na);

% Conversions
b = 1e-9;           % m per Kuhn segment (ball parked based on typical mer length)
kd_code = 1;        % unit inverse time
kb_code = 1;        % unit energy

MetersPerUnitLength = b/b_code;
KelvinPerDegree = T_out/T_code;
JoulesPerUnitEnergy = kbT/kb_code/T_code;
SecondsPerUnitTime = kd_code/kd_out;
KgPerUnitMass = kb/((MetersPerUnitLength)^2)*(SecondsPerUnitTime^2)* ...
    KelvinPerDegree;
NPerUnitForce = KgPerUnitMass*MetersPerUnitLength/SecondsPerUnitTime^2;
PaPerUnitStress = NPerUnitForce/MetersPerUnitLength^2;
KgSPerUnitDamper = KgPerUnitMass/SecondsPerUnitTime;

Tab = table(MetersPerUnitLength,SecondsPerUnitTime,KelvinPerDegree,...
    KgPerUnitMass,NPerUnitForce,PaPerUnitStress,JoulesPerUnitEnergy,...
    KgSPerUnitDamper);
writetable(Tab,'Unit Conversions.txt')

if SIOrNormalized==0
    % Convert inputs
    kb_out = kb;
    b_out = MetersPerUnitLength*b_code;         % m
    dt_out = dt_code*SecondsPerUnitTime;        % s
    damp_out = damp_code*KgSPerUnitDamper;      % kg/s
    f0_out = f0_code;                           % No conversion as f0 is provided in units of 3kbT/sqrt(N)b
    m_out = m/Na/1000;                          % kg

    TransformLengthScales(MetersPerUnitLength);    %Convert to Angstroms
elseif SIOrNormalized==1
    kb_out = 1;
    b_out = b_code;
    dt_out = dt_code;
    damp_out = damp_code;
    f0_out = f0_code;
    m_out = 0.05;
elseif SIOrNormalized==2
    % Convert inputs
    kb_out = kb*(1e21)*((1e9)^2)/(1e9)^2;
    b_out = MetersPerUnitLength*b_code*1e9;         % nm
    dt_out = dt_code*SecondsPerUnitTime*1e9;        % ns
    damp_out = damp_code*KgSPerUnitDamper*1e21/1e9;     % attogram/nm
    f0_out = f0_code;                                   % No conversion as f0 is provided in units of 3kbT/sqrt(N)b
    m_out = m/Na/1000*1e21;                          % attograms

    kd_out = kd_out/1e9;
    ka_out = ka_out/1e9;

    TransformLengthScales(MetersPerUnitLength);    %Convert to Angstroms
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TransformLengthScales(dxdX)

global dims Positions PosUnwrapped X_SP Y_SP Z_SP Rx_SP Ry_SP Rz_SP...
    Lx Ly Lz

% Transform all pos vecs
X_SP = dxdX*X_SP;
Y_SP = dxdX*Y_SP;
Z_SP = dxdX*Z_SP;
Rx_SP = dxdX*Rx_SP;
Ry_SP = dxdX*Ry_SP;
Rz_SP = dxdX*Rz_SP;

Positions(:,1:dims) = Positions(:,1:dims)*dxdX;
PosUnwrapped(:,1:dims) = PosUnwrapped(:,1:dims)*dxdX;

Lx = Lx*dxdX;
Ly = Ly*dxdX;
Lz = Lz*dxdX;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CheckTimescales(kbT,Na)

Npts = 1000;
FontSize = 20;
LineWidth = 1.5;
alpha = 0;
Jperkcal = 4184;
CheckFig = 0;

% Realistic ranges of segmental relaxation timescale and bond activation
% energies (kcal/mol)

tau0_min = 0;       %ps (https://pubs.acs.org/doi/pdf/10.1021/ma101721d)
tau0_max = 3;       %ps

tau0_range = (logspace(tau0_min,tau0_max,Npts))';   %ps
tau0_range = tau0_range*1e-12;                      %s

Ea_min = 4.76;     	%kcal/mol Hydrogen bonds
Ea_max = 50;        %kcal/mol vitrimeric bonds https://pubs.rsc.org/en/content/articlelanding/2020/ta/d0ta06264b

Ea_range_kcal = (linspace(Ea_min,Ea_max,Npts))';  	%kcal/mol
Ea_range = Ea_range_kcal*Jperkcal;                 	%J/mol
Ea_range = Ea_range/Na;                            	%J

Ratio_range = Ea_range/kbT;     %Activation:thermal energy ratio

[~,Ratio_grid] = meshgrid(tau0_range,Ratio_range);
[tau_grid,Ea_kcal_grid] = meshgrid(tau0_range,Ea_range_kcal);

kd_grid = exp(-Ratio_grid)./tau_grid;

if CheckFig==1
    figure(1); clf; hold on
    s = surf(tau_grid/1e-9,Ea_kcal_grid,kd_grid);
    s.FaceColor = 'interp';
    s.EdgeAlpha = alpha;

    set(gca,'XScale','log')
    set(gca,'ZScale','log')
    set(gca,'ColorScale','log')
    c = colorbar;

    xlim([min(tau0_range/1e-9),max(tau0_range/1e-9)])
    ylim([min(Ea_range_kcal),30])
    caxis([1e-8 3.12e8])

    view(0,90)
    pbaspect([1 1 1])

    set(gca,'FontSize',FontSize/1.5)
    xlabel('$\tau_0$ [ns]','FontSize',FontSize,'Interpreter','latex')
    ylabel('$E_a$ [kcal mol$^{-1}$]','FontSize',FontSize,'Interpreter','latex')
    c.Label.String = '$k_d$ [Hz]';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = FontSize;
    c.Limits = [1e-8 3e8];
    box on
    set(gcf,'Color','w')

    % Add contours
    kds = [0.0001 0.01 1 100 10000];
    for k=1:length(kds)
        kd = kds(k);
        Ea_plot = -kbT*log(tau0_range*kd)*Na/Jperkcal;
        z_plot = kd*2*ones(size(Ea_plot));

        p = plot3(tau0_range/1e-9,Ea_plot,z_plot);
        p.Color = 'k';
        p.LineWidth = LineWidth;
    end
end

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

global FileTag Positions dims Lx Ly Lz Np ConnectionsSP...
    N_atoms N_bonds TopologyFileName N Ns BeadSpringOrMeso

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
if BeadSpringOrMeso==0
    N_atom_types = 3;   % Perms and xls
    N_bond_types = 2;   % Kuhn segments and dynamic bonds
elseif BeadSpringOrMeso==1
    N_atom_types = 2;  % Perms and xls
    N_bond_types = 3;   % Short, dynamics vs long implicit chains
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
    if BeadSpringOrMeso==0
        if type==1
            mass = mass_backbone;
        elseif type==2 || type==3
            mass = mass_sticker;
        end
    elseif BeadSpringOrMeso==1
        if type==1
            mass = mass_backbone;
        elseif type==2
            mass = mass_sticker;
        end
    end
    L = ['\n',num2str(type),' ',num2str(mass)];
    Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% Append bond coefficients - for Langevin spring potentials
L = '\nBond Coeffs';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
for i=1:N_bond_types
    type = i;
    if BeadSpringOrMeso==0
        if type==1
            %         coeff1 = 1.1*3*kb*T/b^2;
            %         coeff2 = 0.85*b;
            coeff1 = 800*kbT;
            coeff2 = b;
            coeff3 = b;
        elseif type==2
            %         coeff1 = 1.1*3*kb*T/b^2;
            %         coeff2 = 0.85*b;
            coeff1 = 800*kbT;
            coeff2 = b;
            coeff3 = b;
        end
        %     L = ['\n',num2str(type),' ',...
        %         num2str(coeff1),' ',num2str(coeff2)];
        L = ['\n',num2str(type),' ',...
            num2str(coeff1),' ',num2str(coeff2),' ',num2str(coeff3)];
    elseif BeadSpringOrMeso==1
        if type==1 || type==2
            coeff1 = 3/2*kb*T/(N_Kuhn*b^2);
            coeff2 = 0;
        elseif type==3
            coeff1 = 3*kb*T/b^2;
            coeff2 = 0;
        end
        L = ['\n',num2str(type),' ',...
            num2str(coeff1),' ',num2str(coeff2)];
    end
    Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% Append coordinates
L = '\nAtoms';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
if BeadSpringOrMeso==0
    N_primary = Np*2;
    TetherGroup  = 1:N_primary;
    EndGroups  = find(ConnectionsSP(:,2)==0);
elseif BeadSpringOrMeso==1
    N_primary = Np*Ns;
end
for i=1:N_atoms
    tag = Positions(i,dims+1);
    mol_tag = Positions(i,dims+2);
    x = Positions(i,1); y = Positions(i,2);
    if dims==2
        z = 0;
    else
        z = Positions(i,3);
    end
    if BeadSpringOrMeso==0
        if ismember(i,EndGroups)
            atom_type = 3;
        elseif ismember(i,TetherGroup)
            atom_type = 1;
        else
            atom_type = 2;
        end
    elseif BeadSpringOrMeso==1
        if i<=N_primary
            atom_type = 1;
        else
            atom_type = 2;
        end
    end
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
L = '\nBonds';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
Backbones = Positions(1:N,dims+1);
% Stickers = Positions(N+1:end,dims+1);
for i=1:N_bonds
    tag = i;
    atom_1 = UniquePairs(i,1); atom_2 = UniquePairs(i,2);
    if BeadSpringOrMeso==0
        bond_type = 1;
    elseif BeadSpringOrMeso
        if ismember(atom_1,Backbones) && ismember(atom_2,Backbones)
            bond_type = 1;
        else
            bond_type = 2;
        end
    end
    L = ['\n',num2str(tag),' ',num2str(bond_type),' ',...
        num2str(atom_1),' ',num2str(atom_2)];
    Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

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
function SaveForceData(N_Kuhn_reset,b_reset,kbT_reset)

global dims Rx_SP Ry_SP Rz_SP kbT N_Kuhn b N_bonds Lx Ly Lz...
    ChainForceFileName

kbT = kbT_reset;
N_Kuhn = N_Kuhn_reset;
b = b_reset;

if dims==2
    V = Lx*Ly;
else
    V = Lx*Ly*Lz;
end
N_mers = N_Kuhn*N_bonds;
V_mer = 4/3*pi*b^3;
V_mers = N_mers*V_mer;

rho = N_mers/V;     %density
phi_p = V_mers/V;   %polymer packing
phi_s = 1-phi_p;      %solvent or free volume

norms = vecnorm([Rx_SP(:),Ry_SP(:),Rz_SP(:)],2,2);
norms(norms==0) = [];
disp(['Max. Bond Length = ',num2str(max(norms))])

% Tabulates force data as follows for LAMMPS:

% # Bond potential for harmonic (one or more comment or blank lines) 
% 
% HAM                             (keyword is the first text on line) 
% N 101 FP 0 0 EQ 0.5             (N, FP, EQ  parameters)                              
%                                 (blank line) 
% 1 0.00 338.0000 1352.0000       (index, bond-length, energy, force) 
% 2 0.01 324.6152 1324.9600 
% ... 
% 101 1.00 338.0000 -1352.0000 

Nb = N_Kuhn*b;
Npts = 100;
r = (linspace(0,2*Nb,Npts))';

f = CalculateLangevinForce(r); 

%Use Tangent modulus at high stretch for numerical stability
StretchLimit = 0.99; %Percent
StretchLimit_r = StretchLimit*Nb;
margin = b/100;
Force1 = CalculateLangevinForce(StretchLimit_r-margin);
Force2 = CalculateLangevinForce(StretchLimit_r+margin);
slope = (Force2-Force1)/(margin*2);
yIntercept = Force2-slope*(StretchLimit_r+margin);
f(r>StretchLimit_r) = r(r>StretchLimit_r)*slope+yIntercept;
f = -f;

U = zeros(size(r));
for i=1:length(r)-1
    U(i+1) = trapz(r(1:i+1),f(1:i+1));
end

FileTag = ['N_',num2str(N_Kuhn),'.b_',num2str(b)];

% Initialize string
Assembled = [];

% Header
L = ['Initial Topology',FileTag]; Assembled = [Assembled,L];
L = '\n    '; Assembled = [Assembled,L];
L = '\nHAM'; Assembled = [Assembled,L];
L = ['\nN ',num2str(Npts)]; Assembled = [Assembled,L];
L = '\n    '; Assembled = [Assembled,L];
for i=1:Npts
    L = ['\n ',num2str(i),' ',num2str(r(i),'%.3f'),' ',...
        num2str(U(i),'%.4f'),' ',num2str(f(i),'%.4f')]; 
    Assembled = [Assembled,L];
end

% Save data
fid = fopen(ChainForceFileName,'wt');
fprintf(fid,Assembled);

CheckFig=1;
if CheckFig==1
    figure(10)
    clf
    hold on

    r0 = sqrt(N_Kuhn)*b;

    subplot(2,1,1)
    hold on
    plot(r,f,'k');
    plot([mean(norms),mean(norms)],[0,max(f)],'k--')
    plot([max(norms),max(norms)],[0,max(f)],'r--')
    plot([r0 r0],[0,max(f)],'b--')
    xlim([0 Nb])
    ylim([0 -Force2*2])
    xlabel('r')
    ylabel('f')

    indx = find(f>-Force2,1,'first');

    subplot(2,1,2)
    hold on
    plot(r,U,'k');
    plot([mean(norms),mean(norms)],[0,max(U)],'k--')
    plot([max(norms),max(norms)],[0,max(U)],'r--')
    plot([r0 r0],[0,max(U)],'b--')
    xlabel('r')
    ylabel('U')
    xlim([0 Nb])
    ylim([0 U(indx)])
    MaxU = trapz(r,f);
    scatter(r(end),MaxU,'kx')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N,b] = SetLengthScales(L,phi)

global Lx Ly Lz dims Positions ConnectionsSP

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

%Define total volume
if dims==2
    V = Lx*Ly;
else
    V = Lx*Ly*Lz;
end
conc = N_atoms/V;
Check = 0;
if Check==1
    disp([num2str(conc),' xls per length cubed'])
end
Vp = V*phi;             %Polymer volume

b = sqrt(6/pi*Vp/N_bonds/L);
N = round(L/b);
b = L/N;                %Conserve for rounding

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveMatlabTopologyData

global FileTag TopologyDir...
    X_SP Y_SP Z_SP Rx_SP Ry_SP Rz_SP...
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
function UnwrapPeriodic

global Positions dims N X_SP Y_SP Z_SP Rx_SP Ry_SP Rz_SP...
    ConnectionsSP PosUnwrapped X_SP_Unwrapped Y_SP_Unwrapped Z_SP_Unwrapped

Np = size(unique(Positions(:,dims+2)),1);
Backbones = Positions(1:N,dims+1);
Stickers = Positions(N+1:end,dims+1);

PosUnwrapped = zeros(size(Positions));
PosUnwrapped(:,dims+1:end) = Positions(:,dims+1:end);

X_SP_Unwrapped = zeros(size(X_SP));
Y_SP_Unwrapped = zeros(size(Y_SP));
Z_SP_Unwrapped = zeros(size(Z_SP));

CheckFigs =0;

xl_ct = 0;
for mol=1:Np
    TempPositions = Positions;
    TempPositions(TempPositions(:,dims+2)~=mol,:) = [];
    TempPositions = sortrows(TempPositions,dims+3);
    for j=1:size(TempPositions,1)
        xl_ct = xl_ct+1;
        xl = TempPositions(j,dims+1);
        if ismember(xl,Backbones)
            if j==1
                Pos_xl = Positions(xl,1:dims);
            else
                indx = find(ConnectionsSP(xl0,:)==xl); 
                dx = Rx_SP(xl0,indx);
                dy = Ry_SP(xl0,indx);
                dX = [dx dy];
                if dims==3
                    dz = Rz_SP(xl0,indx);
                    dX = [dx dy dz];
                end
                Pos_xl = PosUnwrapped(xl0,1:dims) - dX;
            end
            if CheckFigs==1
                if dims==2
                    s = scatter(Pos_xl(1),Pos_xl(2),'c','filled');
                    s.SizeData = 12;
                else
                    s = scatter3(Pos_xl(1),Pos_xl(2),Pos_xl(3),'c','filled');
                    s.SizeData = 12;
                end
            end
            xl0 = xl;
            PosUnwrapped(xl,1:dims) = Pos_xl;
        else
            Bonded_xl = ConnectionsSP(xl,1);
            X_tether = PosUnwrapped(Bonded_xl,1:dims);
            indx = find(ConnectionsSP(xl,:)==Bonded_xl);
            dx = Rx_SP(xl,indx);
            dy = Ry_SP(xl,indx);
            dX = [dx dy];
            if dims==3
                dz = Rz_SP(xl,indx);
                dX = [dx dy dz];
            end
            PosUnwrapped(xl,1:dims) = X_tether + dX;
            if CheckFigs==1
                if dims==2
                    s = scatter(PosUnwrapped(xl,1),PosUnwrapped(xl,2),...
                        'g','filled');
                    s.SizeData = 12;
                else
                    s = scatter3(PosUnwrapped(xl,1),PosUnwrapped(xl,2),...
                        PosUnwrapped(xl,3),'g','filled');
                    s.SizeData = 12;
                end
            end
        end
        
    end
    
    for j=1:size(Stickers,1)
    end
end
[X_SP_Unwrapped,Y_SP_Unwrapped,Z_SP_Unwrapped] =...
    ReIndexXYUnwrapped(ConnectionsSP,...
    X_SP_Unwrapped,Y_SP_Unwrapped,Z_SP_Unwrapped,Rx_SP);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AddTetheredChains

global dims N Positions Ns ConnectionsSP...
    X_SP Y_SP Z_SP Rx_SP Ry_SP Rz_SP N_Kuhn b

%% ONLY WORKS IF NS=1. SET UP A LOOP OVER ADD RANGE IF WANT TO MAKE IT GENERAL

N_stickers = N*Ns;
AddRange = (N+1:N+N_stickers)';
Positions(AddRange,:) = 0;
Positions(AddRange,dims+1) = AddRange; %Add xl numbers

theta = 2*pi*rand(N_stickers,1);
% varphimin = 0;
varphimax = 2*pi;
varphi = varphimax*rand(N_stickers,1);
% radius = normrnd(0,N_Kuhn*b^2,N_stickers,1);

% Pick a radius
radius = zeros(size(varphi));
for i=1:size(varphi,1)
    pos_t = [0 0 0];
    for Segment=1:N_Kuhn
%         ct = ct+1;
        theta = 2*pi*rand(1);
        varphi_t = 2*pi*rand(1);
        norm = b;

        rad_proj = norm*cos(varphi_t);
        dx = rad_proj*cos(theta);
        dy = rad_proj*sin(theta);
        dz = norm*sin(varphi_t);

        pos_t = pos_t+[dx dy dz];
% 
%         if Segment==1
%             OldIndx = find(Positions(:,end)==TetherNode,1,'first');
%             OldPoint  = Positions(OldIndx,1:dims);
%             OldNode = TetherNode;
%         else
%             OldPoint = Positions(NewNode-1,1:dims);
%             OldNode = NewNode-1;
%         end
% 
%         Positions(NewNode,1:dims) = OldPoint + [dx dy dz];
%         Positions(NewNode,dims+1) = NewNode;
%         Positions(NewNode,dims+2) = TetherNode;
% 
%         ConnectionsSP(OldNode,2) = NewNode;
%         ConnectionsSP(NewNode,1) = OldNode;
    end
    radius(i,1) = vecnorm(pos_t,2,2);
end

OldPoints = Positions(1:N,1:dims);
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

CheckFig = 0;
if CheckFig==1
    figure(1)
    hold on
    s=scatter3(Positions(1:N,1),Positions(1:N,2),Positions(1:N,3),'k','filled');
    s.SizeData = 6;
    s=scatter3(Positions(AddRange,1),Positions(AddRange,2),Positions(AddRange,3),'r','filled');
    s.SizeData=6;
    close all
end

% ONLY WORKS IF NS=1
ConnectionsSP = zeros(N+Ns*N,Ns+1);
ConnectionsSP(:,end) = (1:size(ConnectionsSP,1))';
ConnectionsSP(1:N,1) = ConnectionsSP(N+1:end,2);
ConnectionsSP(N+1:end,1) = ConnectionsSP(1:N,2);    %Reciprocate

X_SP = zeros(size(ConnectionsSP));
Y_SP = zeros(size(ConnectionsSP));
Z_SP = zeros(size(ConnectionsSP));
for i=1:size(X_SP,1)
    for j=1:size(X_SP,2)
        X_SP(i,j) = Positions(ConnectionsSP(i,j),1);
        Y_SP(i,j) = Positions(ConnectionsSP(i,j),2);
        Z_SP(i,j) = Positions(ConnectionsSP(i,j),3);
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
function EquilibrateTheDomain

global ConnectionsSP ConnectionsREP X_SP X_REP Positions Positions0...
    Y_SP Y_REP Z_SP Z_REP Rx_SP Rx_REP Ry_SP Ry_REP Rz_SP Rz_REP...
    ResMax ResAvg NumericalDamper RepulsionRedefFreq dt dims  

F_NET = CalculateNetForces;

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
       
    dX = Damper3*F_NET;
       
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

%     %Redefine Corresponding X and Y Matrices (Leave Rx and Ry alone)
%     [X_SP,Y_SP] = ReIndexXY(ConnectionsSP,X_SP,Y_SP,Rx_SP);
%     [X_REP,Y_REP] = ReIndexXY(ConnectionsREP,X_REP,Y_REP,Rx_REP);
elseif dims==3
    C1 = Corners(1,:);
    C2 = Corners(2,:);
    C4 = Corners(4,:);
    C5 = Corners(5,:);
    Lx = norm(C2-C1);
    Ly = norm(C1-C4);
    Lz = norm(C5-C1);
    X = Positions(:,1); Y = Positions(:,2); Z = Positions(:,3);
    dXdY = C4(1)-C1(1);
    dYdX = C2(2)-C1(2);
    dYdZ = C5(2)-C1(2);
    dZdY = C4(3)-C1(3);
    dZdX = C2(3)-C1(3);
    dXdZ = C5(1)-C1(1);

    %Find particles out of left bound and reposition
    v1 = C4-C1;
    v2 = C5-C1;
    C = C1;
    Out = FindOutOfBounds(v1,v2,C);

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
function [X,Y,Z] = ReIndexXYUnwrapped(Connections,X,Y,Z,Rx)

global PosUnwrapped dims

X(Connections(:,1:end-1)==0) = NaN;
Y(Connections(:,1:end-1)==0) = NaN;
Z(Connections(:,1:end-1)==0) = NaN;
Connections(Connections(:,1:end-1)==0) = NaN;

%Define New Positions in X and Y matrices
X(:,end) = PosUnwrapped(:,1);
Y(:,end) = PosUnwrapped(:,2);
if dims==3
    Z(:,end) = PosUnwrapped(:,3);
end
for j=1:size(Rx,2)
    TempNodes = Connections(:,j);
    Rows = find(~isnan(TempNodes));
    
    %Update X and Y of Tails
    X(Rows,j) = PosUnwrapped(TempNodes(Rows),1);
    Y(Rows,j) = PosUnwrapped(TempNodes(Rows),2);  
    if dims==3
        Z(Rows,j) = PosUnwrapped(TempNodes(Rows),3);
    end
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
% if RepulsionType~=2
    F_REP = CalculateRepulsiveForces;
% else
%     F_REP = zeros(size(F_SP));
% end

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
% mu = 0;
% sigma = r0/30;
% ForceScale = 5;

% DistToBounds = zeros(size(Positions(:,1:2)));
% DistToBounds(Positions(:,2)>=0,2) = Ly/2-Positions(Positions(:,2)>=0,2);
% DistToBounds(Positions(:,2)<0,2) = -Ly/2-Positions(Positions(:,2)<0,2);
%     
% F_WALL = CalculateGuassianForce2(DistToBounds,mu,sigma,ForceScale);

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

global SpringType r0 SpringCoeff IsometricForce CoeffRatio L CutoffCWAForce...

if SpringType==0 %Linear
    F = -SpringCoeff*Norms;%(r0-Norms);
elseif SpringType==1 %Langevin 
%     lam = Norms/L; %L is 2X the length of a tether
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
    F = -IsometricForce*exp(-(((Norms-r0)/r0)/SpringCoeff).^2);
    F(Norms>r0) = -IsometricForce*(exp(-(((Norms(Norms>r0)-r0)/r0)/SpringCoeff).^2)+(((Norms(Norms>r0)-r0)/r0)/(CoeffRatio*SpringCoeff)).^2);
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
    lam = Norms/L;
    FL = CalculateLangevinForce(Norms); %Entropic contribution
    
    %Use Tangent modulus at high stretch for numerical stability
    StretchLimit = 0.90; %Percent
    StretchLimit_r = StretchLimit*L;
    Force1 = CalculateLangevinForce(StretchLimit_r-0.0001);
    Force2 = CalculateLangevinForce(StretchLimit_r+0.0001);
    slope = (Force2-Force1)/0.0002;
    yIntercept = Force2-slope*(StretchLimit_r+0.0001);
    F(Norms>StretchLimit_r) = Norms(Norms>StretchLimit_r)*slope+yIntercept;   
    
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
% Force = f(kb,T,b,lambda)

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
[NeighborsREP_RightLeft,DistREP_RightLeft,...
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
[NeighborsREP_TopBottom,DistREP_TopBottom,...
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
[NeighborsREP_BottomLeftTopRight,DistREP_BottomLeftTopRight,...
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
[NeighborsREP_BottomRightTopLeft,DistREP_BottomRightTopLeft,...
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
    [NeighborsREP_BackFront,DistREP_BackFront,...
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
    [NeighborsREP_FrontLeftBackRight,DistREP_FrontLeftBackRight,...
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
    [NeighborsREP_FrontRightBackLeft,DistREP_FrontRightBackLeft,...
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
    [NeighborsREP_BottomFrontTopBack,DistREP_BottomFrontTopBack,...
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
    [NeighborsREP_TopFrontBottomBack,DistREP_TopFrontBottomBack,...
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
    [NeighborsREP_C1C7,DistREP_C1C7,...
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
    [NeighborsREP_C2C8,DistREP_C2C8,...
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
    [NeighborsREP_C6C4,DistREP_C6C4,...
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
    [NeighborsREP_C5C3,DistREP_C5C3,...
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
%     dXdY = C4(1)-C1(1);
%     dYdX = C2(2)-C1(2);
%     dYdZ = C5(2)-C1(2);
%     dZdY = C4(3)-C1(3);
%     dXdZ = C5(1)-C1(1);
%     dZdX = C2(3)-C1(3);

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
%     C2 = Corners(2,:);  %RFB    Right-Front-Bottom
%     C3 = Corners(3,:);  %RBB    Right-Back-Bottom
%     C4 = Corners(4,:);  %LBB    Left-Back-Bottom
%     C5 = Corners(5,:);  %LFT    Left-Front-Top
%     C6 = Corners(6,:);  %RFT    Right-Front-Top
    C7 = Corners(7,:);  %RBT    Right-Back-Top
%     C8 = Corners(8,:);  %LBT    Left-Back-Top
    
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

    ShowFigDebug = 1;
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
% CTL = Corners(1,:);

%Define min distance between any two corners to determine boundary
%proximity
Lx = CBR(1)-CBL(1);
Ly = CTR(2)-CBR(2);
% D1 = norm(CBL-CBR);
% D2 = norm(CBR-CTR);
% D3 = norm(CTR-CBL);
% D4 = norm(CTL-CBR);
% DomainScale = min([D1,D2,D3,D4]);

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
function SizeTheDomainMesoscale

global Lx Ly Lz N EdgeFactor dims AvgParticlesPerEdge

AvgParticlesPerEdge = round(N^(1/dims));
Lx = EdgeFactor*AvgParticlesPerEdge;
Ly = Lx;
if dims==2
    Lz = 0;
else
    Lz = Lx;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SizeTheDomain(Np,N_Kuhn,b,Separation)

global Lx Ly Lz PairsPerEdge MinSep PairDistance NominalSeparation...
    BeadSpringOrMeso AvgParticlesPerEdge dims EdgeFactor N

if BeadSpringOrMeso==0
    PairsPerEdge = round(Np^(1/3));

    MinSep = 2*N_Kuhn*b;
    PairDistance = Separation;
    NominalSeparation = MinSep+PairDistance;

    Lx = NominalSeparation*PairsPerEdge;
    Ly = Lx;
    Lz = Lx;
else
    AvgParticlesPerEdge = round(N^(1/dims));
    Lx = EdgeFactor*AvgParticlesPerEdge;
    Ly = Lx;
    if dims==2
        Lz = 0;
    else
        Lz = Lx;
    end
end

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

global Positions Lx Ly Lz PairsPerEdge MinSep PairDistance BeadSpringOrMeso...
    EdgeFactor AvgParticlesPerEdge dims

if BeadSpringOrMeso==0
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

    [X,Y,Z] = meshgrid(XSpacing,YSpacing,ZSpacing);
    X = X(:)-mean(X(:));
    Y = Y(:)-mean(Y(:));
    Z = Z(:)-mean(Z(:));
    Positions = [X(:) Y(:) Z(:)];

    Lx = 2*(max(X)+MinSep/2);
    Ly = Lx;
    Lz = Lx;

    CheckFig = 0;
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
    close all
elseif BeadSpringOrMeso==1
    Spacing = (0:EdgeFactor:AvgParticlesPerEdge-1)' + EdgeFactor/2 - Lx/2;

    if dims==2
        [X,Y] = meshgrid(Spacing,Spacing);
        Positions = [X(:) Y(:)];
    else
        [X,Y,Z] = meshgrid(Spacing,Spacing,Spacing);
        Positions = [X(:) Y(:) Z(:)];
    end

    CheckFig = 0;
    if CheckFig==1
        figure(1)
        clf
        s = scatter3(X(:),Y(:),Z(:),'k','filled');
    end
    s.SizeData = 6;
    close all
end

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
