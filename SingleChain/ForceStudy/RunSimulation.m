function RunSimulation(Package,Override,ToggleDynamics,LC,DC,BSOM,CF,OF,BT,...
    LAMMPSExecutablePath)

global InputFileName OutputAtom_loc OutputBond_loc ...
    OutputAtom OutputBond TurnOnDynamics...
    LengthConversion DamperConversion BeadSpringOrMeso CurrentFolder OutputFolder...
    OutputDir BondType

TurnOnDynamics = ToggleDynamics;
LengthConversion = LC;
DamperConversion = DC;
BeadSpringOrMeso = BSOM;
CurrentFolder = CF;
OutputFolder = OF;
BondType = BT;

%% Unpack Swept Input Parameters
% Swept parameters
Sample = Package(1);    %Sample index
Np = Package(2);        %Number of molecules
D = Package(3);         %Diffusion coefficient [m2/s]
N_Kuhn = Package(4);    %Number of Kuhn segments in chain
stiffness = Package(5); %Stiffness of single harmonic bond

kbT = Package(6);       %Thermal energy
b = Package(7);         %Kuhn length
dt = Package(8);        %timestep

%% Callout Input Script
% Initializes non-sweeping parameters
InputScript(Sample,Np,D,N_Kuhn,stiffness,kbT,b,dt);

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


