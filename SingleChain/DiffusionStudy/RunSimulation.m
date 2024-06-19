function RunSimulation(Package,Override,ToggleDynamics,LC,DC,BSOM,CF,OF,...
    LAMMPSExecutablePath)

global InputFileName OutputAtom_loc OutputBond_loc ...
    OutputAtom OutputBond TurnOnDynamics...
    LengthConversion DamperConversion BeadSpringOrMeso CurrentFolder...
    OutputDir OutputFolder

TurnOnDynamics = ToggleDynamics;
LengthConversion = LC;
DamperConversion = DC;
BeadSpringOrMeso = BSOM;
CurrentFolder = CF;
OutputFolder = OF;

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
InputScript(EdgeFactor,Sample,Np,N,ka,kd,f0,dt,damp,N_Kuhn,b,Separation);

%% Make diriectories and set filenames
SetDirAndFileNames;

if ~isfile(OutputAtom) || ~isfile(OutputBond) || Override==1
%     command = ['wsl /mnt/c/Users/rjwag/Documents/lammps2/lammps/build/lmp -in ',...
%         InputFileName];
    command = ['wsl ',LAMMPSExecutablePath,' -in ',InputFileName];
    system(command); 

    if isfile(OutputAtom_loc)
        movefile(OutputAtom_loc,OutputDir)
        movefile(OutputBond_loc,OutputDir)
    end
end

clear global

end


