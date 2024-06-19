%2020-01-25 General Network Input Script

function DefineMatlabNetworkInitiationInputs

global EdgeFactor sample r0 ForceScaleREP muREP sigmaREP MaxDist MinDist...
    CutoffPartREP visc ResMax ResAvg NumericalDamper RepulsionRedefFreq...     
    MakeMovie xSpeed alpha CutoffCWAForce

%% Call Out Inputs and Globalize
EdgeFactor = 1;       %Factor by which domain is decreased linearly


%% Constitutive Spring Properties
r0 = 1;
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
if sample==1
    MakeMovie=0;
end
xSpeed = 20;	            %Movie Speed


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
