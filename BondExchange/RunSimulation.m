function RunSimulation(Package,Override,ToggleDynamics,LC,DC,BSOM,CF,OF,...
    CDL,lmp_path,LOL)

global InputFileName OutputAtom_loc OutputBond_loc ...
    OutputAtom OutputBond TurnOnDynamics LinearOrLangevin...
    LengthConversion DamperConversion BeadSpringOrMeso CurrentFolder...
    OutputFolder current_dir_lmp OutputDataFolder RawDataFileName dims


TurnOnDynamics = ToggleDynamics;
LengthConversion = LC;
DamperConversion = DC;
BeadSpringOrMeso = BSOM;
CurrentFolder = CF;
OutputFolder = OF;
current_dir_lmp = CDL;
LinearOrLangevin = LOL;

%% Unpack Swept Input Parameters
Sample = Package(1);    %Sample index
Np = Package(2);        %Number of molecules
Ns = Package(3);        %Number of stickeres per tether site
N = Np*Ns;              %Total numbber of permanent crosslinks (i.e., tether sites)
ka = Package(4);        %Activation energy of association
kd = Package(5);        %Activation energy of dissociation
N_Kuhn = Package(7);
b = Package(8);
phi = Package(10);
kbT = Package(11);

[~,D,tau0,dtFact] = DefineTimeScale(b,LengthConversion,DamperConversion,...
    BeadSpringOrMeso);

dt = tau0/dtFact;

eaStar = -log(ka*(b*LengthConversion)^2/D);
edStar = -log(kd*(b*LengthConversion)^2/D);

%% Callout Input Script
% Initializes non-sweeping parameters
dims = 3;   %Dimensionality of system
InputScript(Sample,Np,Ns,N,N_Kuhn,b,phi,eaStar,edStar,D,dt,kbT);

%% Make diriectories and set filenames
SetDirAndFileNames;
DefineCompiledFileNames;

OutputAtom_loc = [OutputDataFolder,OutputAtom];
OutputBond_loc = [OutputDataFolder,OutputBond];
if ((~isfile(OutputAtom_loc) || ~isfile(OutputBond_loc)) &&...
        ~isfile(RawDataFileName)) || Override==1
    setenv('OMP_NUM_THREADS', '2');
    if BeadSpringOrMeso==0
        command = ['wsl mpirun -np 6 ',lmp_path,' -in ',...
            InputFileName];
    else
        command = ['wsl mpirun -np 6 ',lmp_path,' -in ',...
            InputFileName];
    end
    system(command); 

    if isfile(OutputAtom)
        movefile(OutputAtom,OutputDataFolder)
        movefile(OutputBond,OutputDataFolder)
    end
end

clear global

end


