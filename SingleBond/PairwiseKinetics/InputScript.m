%2020-01-25 General Network Input Script

function InputScript(EF_temp,Sample,Np,N_temp,...
     ka,kd,f0,dt,damp,N_Kuhn_temp,b_temp,Separation)

global EdgeFactor r0 tau0 ForceScaleREP muREP sigmaREP MaxDist MinDist N...
    CutoffPartREP visc ResMax ResAvg NumericalDamper RepulsionRedefFreq...     
    MakeMovie xSpeed FileTag FileTagCompiled TopologyFileTag FolderName...
    alpha L CutoffCWAForce ConstantsFileTag N_Kuhn b...
    SpringCoeff SpringType kbT BeadSpringOrMeso TurnOnDynamics...
    OutputFolder OutputDataFolder %CurrentFolder OutputFolderName

%% Call Out Inputs and Globalize
EdgeFactor = EF_temp;       %Factor by which domain is decreased linearly
N = N_temp;                 %No. of crosslinks (i.e., nodes) 

%% Constitutive Spring Properties
r0 = 1;
tau0 = 1e-9;                %diffusion timescale of one Khun segment
N_Kuhn = N_Kuhn_temp;
b = b_temp;
kbT = 1;                    %Thermal Energy
L = N_Kuhn*b;               %Contour Length of Attached Chain or 2x Tether length
SpringCoeff = 3*kbT/(L*b^2);    %Stiffness of chain wrt. stretch
SpringType = 0;             %0 for linear, 1 for Langevin  

%% Network Generation Inputs
MinDist = 0.625*(EdgeFactor/0.8)*r0;
MaxDist = 1.8*MinDist;      %Maximum spring connection distance for initial network structure

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
elseif BeadSpringOrMeso==1
    AddOn = 'Meso';
elseif BeadSpringOrMeso==2
    AddOn = 'Thry';
end
if TurnOnDynamics==0
    AddOn = ['.Control.',AddOn];
end

TopologyFileTag = [AddOn,'.S',num2str(Sample),'.Np',num2str(Np),'.N_k',...
    num2str(N_Kuhn),'.d.',num2str(Separation,'%.3f')];
ConstantsFileTag = [AddOn,'.Np',num2str(Np),'.N_k',num2str(N_Kuhn),...
    '.d',num2str(Separation,'%.3f'),'.damp.',num2str(damp,'%.2e')];
FileTagCompiled = ['.Np',num2str(Np),'.N_k',num2str(N_Kuhn),...
        '.d',num2str(Separation,'%.3f'),...
        '.ka',num2str(ka,'%.2e'),'.kd',num2str(kd,'%.2e'),...
        '.f0',num2str(f0,'%.1f'),'.dt',num2str(dt,'%.2e'),...
        '.damp',num2str(damp,'%.2e')];

FileTag = [AddOn,'.S',num2str(Sample),FileTagCompiled];
FileTagCompiled = [AddOn,FileTagCompiled];
FolderName = 'Data/Input_Data';%['M:/Research/NSF/',CurrentFolder,'Data/Input_Data'];
if ~isfolder(FolderName)
    mkdir(FolderName)
end
% OutputFolderName = ['M:/Research/NSF/',CurrentFolder,'Data/'];
% if ~isfolder(OutputFolderName)
%     mkdir(OutputFolderName)
% end
OutputDataFolderName = [OutputFolder,'Data/'];
if ~isfolder(OutputDataFolderName)
    mkdir(OutputDataFolderName)
end

end
