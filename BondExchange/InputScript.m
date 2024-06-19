%2020-01-25 General Network Input Script

% function InputScript(Sample,Np,Nt,Ns,N_temp,N_Kuhn_temp,b_temp,EF_temp,phi,...
%     ka,kd,f0,dt,T,damp,W,lambda)
function InputScript(Sample,Np_temp,Ns_temp,N_temp,N_Kuhn_temp,b_temp,phi_temp,...
    eaStar_temp,edStar_temp,D_temp,dt_temp,kbT_temp)

global EdgeFactor r0 tau0 ForceScaleREP muREP sigmaREP MaxDist MinDist N...
    CutoffPartREP visc ResMax ResAvg NumericalDamper RepulsionRedefFreq...     
    MakeMovie xSpeed FileTag FileTagCompiled TopologyFileTag...
    alpha CutoffCWAForce N_Kuhn b eaStar edStar LengthConversion...
    kbT BeadSpringOrMeso ConstantsFileTag TurnOnDynamics...
    InitialNodesTag phi D dt LinearOrLangevin

%% Call Out Inputs and Globalize
EdgeFactor = 1;       %Factor by which domain is decreased linearly
Np = Np_temp;
Ns = Ns_temp;
N = N_temp;
N_Kuhn = N_Kuhn_temp;
b = b_temp;
phi = phi_temp;
eaStar = eaStar_temp;
edStar = edStar_temp;
D = D_temp;
dt = dt_temp;
kbT = kbT_temp;

%% Constitutive Spring Properties
r0 = 1;
tau0 = (b*LengthConversion)^2/D;                %diffusion timescale of one Khun segment
DefineSpringProperties;
                      
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
CutoffCWAForce = 1/2;%1/10*4;

%% Computational Inputs for Initial Network Generation (only)
visc = 500;                 %Viscous overdamping parameter
dt0 = 0.075;
ResMax = 0.9;               %Permissible max force on node at eq.
ResAvg = 0.45;              %Permissible avg force on nodes at eq.
NumericalDamper = dt0/visc; %Damper for equilibrating domain at each timestep

%% Make gifs for one sample of each
MakeMovie = 0;	            %1 to make a movie
if Sample==1
    MakeMovie=0;
end
xSpeed = 20;	            %Movie Speed


%% Directory Inputs
if BeadSpringOrMeso==0
    AddOn = 'Bead';
elseif BeadSpringOrMeso==1
    AddOn = 'Meso';
	if LinearOrLangevin==1
		AddOn = [AddOn,'.Langevin'];
	end
elseif BeadSpringOrMeso==2
    AddOn = 'Thry';
end
if TurnOnDynamics==0
    AddOn = ['.Control.',AddOn];
end

InitialNodesTag = [AddOn,'.S',num2str(Sample),'.Np',num2str(Np),'.Nk',num2str(N_Kuhn),...
    '.Nx',num2str(N),'.Ns',num2str(Ns),'.phi',num2str(phi,'%.3f')];    
        % This uniquely defines initial tethered node positions
TopologyFileTag = InitialNodesTag;   %This uniquely defines initial topology

FileTagCompiled = ['.Np',num2str(Np),'.Nk',num2str(N_Kuhn),...
    '.ea',num2str(eaStar,'%.2f'),'.ed',num2str(edStar,'%.2e'),...
    '.dt',num2str(dt,'%.2e'),'.phi',num2str(phi,'%.3f')];

ConstantsFileTag = FileTagCompiled;
FileTag = [AddOn,'.S',num2str(Sample),FileTagCompiled];
FileTagCompiled = [AddOn,FileTagCompiled];


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DefineSpringProperties

global N_Kuhn_eq b_eq L_eq SpringCoeff_eq SpringType kbT_eq 

kbT_eq = 5e-2;                 %Thermal Energy
L_eq = 3.5;                    %Contour Length of Attached Chain or 2x Tether length
b_eq = 0.005;                  %Length of Kuhn Segment
N_Kuhn_eq = L_eq/b_eq;               %Number of Kuhn Segments
SpringCoeff_eq = 3*kbT_eq/(L_eq*b_eq);	%Stiffness of chain
SpringType = 0;             %0 for linear, 1 for Langevin, 2 for Hill Model,
                            %3 for semiflexible filaments
                            %4 for Langevin with CWA-like repulsion
   
end
