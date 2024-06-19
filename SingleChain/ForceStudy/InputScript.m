%2020-01-25 General Network Input Script

function InputScript(Sample,Np,D,N_Kuhn,stiffness,kbT,b,dt)

global r0 tau0 ForceScaleREP muREP sigmaREP...
    CutoffPartREP visc ResMax ResAvg NumericalDamper RepulsionRedefFreq...     
    MakeMovie xSpeed FileTag FileTagCompiled TopologyFileTag FolderName...
    alpha L CutoffCWAForce kbT_eq...
    SpringCoeff SpringType BeadSpringOrMeso...
    OutputFolder OutputDataFolderName LengthConversion BondType

%% Constitutive Spring Properties
r0 = 1;
tau0 = (b*LengthConversion)^2/D;                %diffusion timescale of one Khun segment
kbT_eq = 1;                        %Thermal Energy
L = N_Kuhn*b;                   %Contour Length of Attached Chain or 2x Tether length
SpringCoeff = 3*kbT_eq/(L*b^2);    %Stiffness of chain wrt. stretch
SpringType = 0;                 %0 for linear, 1 for Langevin  

%% Repulsive Force Properties - **ARBITRARY - FOR STABLE HOMOGENIZATION ONLY**
ForceScaleREP = 1.25;       %epsilon from Chandler-Weeks-Anderson (CWA)
muREP = 1;                  %not used in CWA
sigmaREP =  1.25;           %sigma from CWA
CutoffPartREP = sigmaREP;   %Maximum distance at which repulsive contact forces can exists between particles
alpha = 2;                  %exponent for repulsive potential
RepulsionRedefFreq = 100;
CutoffCWAForce = 1/10*L;

%% Computational Inputs
visc = 500;                 %Viscous overdamping parameter
dt0 = 0.075;
ResMax = 0.9;               %Permissible max force on node at eq.
ResAvg = 0.45;              %Permissible avg force on nodes at eq.
NumericalDamper = dt0/visc; %Damper for equilibrating domain at each timestep

%% Make gifs for one sample of each
MakeMovie = 0;	            %1 to make a movie
if Sample==1
    MakeMovie=1;
end
xSpeed = 20;	            %Movie Speed

%% Directory Inputs
if BeadSpringOrMeso==0
    AddOn = 'Bead';
else
    AddOn = 'Meso';
end
if BondType==0
    AddOn = [AddOn,'.Harmonic'];
else
    AddOn = [AddOn,'.Nonlinear'];
end

TopologyFileTag = [AddOn,'.S',num2str(Sample),...
    '.Np',num2str(Np),...
    '.Nk',num2str(N_Kuhn),...
    '.stiff',num2str(stiffness)];
FileTagCompiled = ['.Np',num2str(Np),...
    '.D',num2str(D,'%.1e'),...
    '.Nk',num2str(N_Kuhn),...
    '.Stiff',num2str(stiffness,'%.2f'),...
    '.kbT',num2str(kbT,'%.ef'),...
    '.b',num2str(b,'%.2f'),...
    '.dt',num2str(dt,'%.2e')];

FileTag = [AddOn,'.S',num2str(Sample),FileTagCompiled];
FileTagCompiled = [AddOn,FileTagCompiled];

% Where data will be stored:
FolderName = 'Data/Input_Data';
if ~isfolder(FolderName)
    mkdir(FolderName)
end
OutputDataFolderName = [OutputFolder,'Data/'];
if ~isfolder(OutputDataFolderName)
    mkdir(OutputDataFolderName)
end

end
