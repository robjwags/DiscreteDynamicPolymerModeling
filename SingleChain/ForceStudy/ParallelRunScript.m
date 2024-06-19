function ParallelRunScript

clear
clear global
clc
fclose all;
close all
close all hidden

%% Unit conversions
LengthConversion = 3.74e-9;     %meters per unit code length
DamperConversion = 2.89e-4;     %N*s/m per unit damper in code
ForceConversion = 1.08e-12;     %N per unit force in code

%% Work flow controls
% Toggles - runs these processes as called out - will skip for cases where files already exists
GenerateTopologies = 1;     % 1 to generate input files 
RunLAMMPS = 1;              % 1 to run LAMMPS simulations 
CompileEnsembleData = 1;    % 1 to parse raw data and compute outputs
PostProcess = 1;            % 1 to post-process the parsed/compiled data and generate output plots
CompareModels = 1;          % 1 to compare bead-spring model to implicit Langeving model

ToggleDynamics = 0;         %1 for dynamics 0 for control systems (no dynamics) DO NOT TOUCH

% Overrides - reruns these processes even if files already exists
OverrideTopologies = 0;     %1 to override initiated topologies
OverrideInputScripts = 1;   %1 to rewrite the input scripts
OverrideRun = 0;            %1 to rerun simulations
OverrideCompile = 0;        %1 to override parsing of LAMMPS outputs
OverrideCompute = 0;        %1 to override computation of stress, end-to-end vecs, MSDs, etc.

% Directory and processor control
NoProcessors = 1;
CurrentFolder = DefineCurrentFolder;    % Ends with /
OutputDrive = 'M:';
OutputFolder = strrep(CurrentFolder, 'C:', OutputDrive);
LAMMPSExecutablePath = '/mnt/c/Users/rjwag/Documents/lammps/build/lmp';

% Set which outputs to calculate
CalculateForces = 1;        % 1 to compute averge force in bonds
CalculateEndtoEnd = 1;      % 1 to compile end-to-end distribution data and 
CalculateAlignmentt = 1;    % 1 to compute metric tensors (r \otimes r) of chains that elucidates alignment

for BeadSpringOrMeso=0:1    % Loop over bead-spring (0) and mesoscale (1) models
    %% Input Parameters
    Np = 1;
    b = 0.1667;     % Kuhn length
    kbT = 293*1.38e-23;     % Thermal energy in Joules
    
    %% Sweeping Parameters
    b_SI = b*LengthConversion;  % m
    N_Kuhn = [12 18 36];    % Number of Kuhn lengths per chain (swept)
    D_nom = 1e-10;          % nominal diffusion coefficient of a monomer in m^2/s
    D = D_nom*10.^[-2 0] ;  % m^2/s
    tau0 = (b_SI^2)./D;     % s
    damps = kbT./D;                     % kg/s or N s/m
    damps = damps/DamperConversion;     % single momonmer damping coefficient (effective in mesoscle depends on N)
    dtFact = 320;                       % Factor by which tau0 is divided to attain timestep size
    dt = 1/dtFact*min(tau0);            % NOTE THAT Diffusion Coeff, damper, and timestep are all linked
    D = D(1); damps = damps(1);
    Stiffnesses_SI = [100 200 400 800]*kbT/b_SI^2;  % Force scale that defines Kuhn bond stiffness (in SI units)
    Stiffnesses = Stiffnesses_SI/kbT*b_SI^2;        % In normalized units of the code
    BondType = 1;                       % 0 for Harmonic, 1 for nonlinear (keep as 1)
   
    if BeadSpringOrMeso==0
        Samples = 1:15;     % Sample set to be swept over
    else
        Samples = 1;
    end
    NoParam = 8;
    
    %% Compile Input Parameters
    Perms = length(N_Kuhn)*length(damps)*length(Stiffnesses)*length(Samples);
    Package = zeros(Perms,NoParam); ct = 0;
    for i=1:length(Samples)
        for j=1:length(Np)
            for k=1:length(D)
                for l=1:length(N_Kuhn)
                    for m=1:length(Stiffnesses)
                        ct = ct+1;

                        % Swept parameters
                        Package(ct,1) = Samples(i);
                        Package(ct,2) = Np(j);
                        Package(ct,3) = D(k);
                        Package(ct,4) = N_Kuhn(l);
                        Package(ct,5) = Stiffnesses(m);

                        % Constants
                        Package(ct,6) = kbT;
                        Package(ct,7) = b;
                        Package(ct,8) = dt;
                   end
                end
            end
        end
    end
   
    %% Run parfor to generate topologies
    if GenerateTopologies==1
        if NoProcessors>1
            delete(gcp('nocreate'))
            parpool(NoProcessors)
            parfor n=1:size(Package,1)
                GenerateTopology(Package(n,:),OverrideTopologies,...
                    OverrideInputScripts,CurrentFolder,OutputFolder,...
                    LengthConversion,DamperConversion,...
                    BeadSpringOrMeso,dtFact,BondType);
            end
        else
            wb2 = waitbar(0,'Generating all input scripts...');
            for m1=1:size(Package,1)
                GenerateTopology(Package(m1,:),OverrideTopologies,...
                    OverrideInputScripts,CurrentFolder,OutputFolder,...
                    LengthConversion,DamperConversion,...
                    BeadSpringOrMeso,dtFact,BondType);
                waitbar(m1/size(Package,1),wb2,'Generating all input scripts...')
            end
            close(wb2)
        end
    end

    %% Loop over LAMMPS sims
    if RunLAMMPS==1
        wb = waitbar(0,'Running all jobs...');
        for m1=1:size(Package)
            RunSimulation(Package(m1,:),OverrideRun,ToggleDynamics,...
                LengthConversion,DamperConversion,BeadSpringOrMeso,...
                CurrentFolder,OutputFolder,...
                BondType,LAMMPSExecutablePath);
            waitbar(m1/size(Package,1),wb,'Running all jobs...')
        end
        close(wb)
    end

    %% Post Process Data
    if CompileEnsembleData==1
        CompileData(Package,OverrideCompile,OverrideCompute,...
            LengthConversion,DamperConversion,BeadSpringOrMeso,...
            CalculateForces,CalculateEndtoEnd,...
            CalculateAlignmentt,CurrentFolder,OutputFolder,...
            NoProcessors,BondType);
    end
end

%% Post Process Data
if PostProcess==1
    PostProcessData(Package,ToggleDynamics,...
        LengthConversion,DamperConversion,ForceConversion,...
        BeadSpringOrMeso,CompareModels,...
        CurrentFolder,OutputFolder,BondType);
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CurrentFolder = DefineCurrentFolder

CurrentFolder = pwd;    
CurrentFolder = strrep(CurrentFolder, '\', '/');

end