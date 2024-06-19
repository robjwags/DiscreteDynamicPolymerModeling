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

%% Work flow controls
% Toggles - runs these processes as called out - will skip for cases where files already exists
GenerateTopologies = 1; 
RunLAMMPS = 1;
CompileEnsembleData = 1;
PostProcess = 1;
CompareModels = 1;

ToggleDynamics = 0;             % 1 for dynamics 0 for control systems (no dynamics) - DO NOT CHANGE
MakeHistos = 1;

% Overrides - reruns these processes even if files already exists
OverrideTopologies = 0;         % 1 to initiate topologies
OverrideInputScripts = 0;       % 1 to rewrite the input scripts
OverrideRun = 0;                % 1 to rerun simulations
OverrideCompile = 0;            % 1 to override parsing of LAMMPS outputs
OverrideCompute = 0;            % 1 to override computation of stress, end-to-end vecs, MSDs, etc.
OverridePostProcess = 0;        % 1 to override generation of stress plots
OverrideHistogramMovies = 0;    % 1 to overwrite the histogram movies
OverrideMSDFitData = 0;         % 1 to override MSD plots
OverrideCompareData = 1;        % 1 to override Comparison plots

% Set which outputs to calculate
CalculateStress = 0;            % 1 to compute virial stress-stretch/time data
CalculateEndtoEnd = 1;          % 1 to compile end-to-end distribution data and 
CalculateMSD = 1;               % 1 to compute the mean-square displacement (i.e. diffusion) data
CalculateBondKinetics = 0;      % 1 to compute the bond kinetics and exchange rates
CalculateAlignmentt = 0;        % 1 to compute metric tensors (r \otimes r) of chains that elucidates alignment
MakeMoviesHistos = 0;           % 1 to generate .MOV's of end-to-end distributions
WRTStretchOrLength = 0;         % 0 to make histos wrt. chain stretch (lambda), and 1 for length
HistoUnits = 0;                 % 0 for normalized, 1 for SI (nm)

NoProcessors = 1;
CurrentFolder = DefineCurrentFolder;
OutputDrive = 'M:';
OutputFolder = strrep(CurrentFolder, 'C:', OutputDrive);
LAMMPSExecutablePath = '/mnt/c/Users/rjwag/Documents/lammps/build/lmp';

for BeadSpringOrMeso=0:1
    %% Input Parameters
    Np = 5^3;       % Number of polymer pairs (make a cubic number for periodic RVEs)
    if BeadSpringOrMeso==1
        Np = 11^3;
    end
    Ns = 1;         % Number of free stickers on each tether site
    T = 1;          % Temperatures
    L = 1;          % Nominal spacing (only matters with dynamics)
    ka = 1;         % Normalized bond activation energy for association in 
                    % units of [kbT]
    kd = 1;         % Normalized bond activation energy for dissociation in 
                    % units of [kbT]
    f0 = 10000;     % Force sensitivity for dissociation in Bell's model in 
                    % units of [3kbT/(\sqrt(N) b)] Set very high to turn off
                    % bond dissociation
    b = 0.1667;     % Kuhn length
    
    %% Sweeping Parameters
    N_Kuhn = [12 18 36]; % Complete set initial investigated was [12 18 24 30 36] with consistent trends
    Distances = 1;                      % Separation between tether sites in units of Nb
    

    [dts,damps,D,dtFact] = DefineTimeStep(b,LengthConversion,...
        DamperConversion,BeadSpringOrMeso);
   
    Samples = 1;    % goal is to have enoug sampling from Np
    NoParam = 10;
    UnitType = 1;   %0 for SI units, 1 for normalized (LJ in LAMMPS), 2 for nano
    
    %% Compile Input Parameters
    Perms = length(Np)*length(Ns)*length(T)*length(L)*length(ka)*length(kd)*...
        length(f0)*length(damps)*length(N_Kuhn)*length(b)*length(Distances);
    Package = zeros(Perms,NoParam); ct = 0;
    for i=1:length(Samples)
     for j=1:length(Np)
      for k=1:length(Ns)
       for l=1:length(T)
        for m=1:length(L)
         for n=1:length(ka)
          for o=1:length(kd)
           for p=1:length(f0)
            for q=1:length(damps)
             for r=1:length(N_Kuhn)
              for s=1:length(b)
               for t=1:length(Distances)
                   ct = ct+1;
                   Package(ct,1) = Samples(i);
                   Package(ct,2) = Np(j);
                   Package(ct,3) = Ns(k);
                   Package(ct,4) = T(l);
                   Package(ct,5) = L(m);
                   Package(ct,6) = ka(n);
                   Package(ct,7) = kd(o);
                   Package(ct,8) = f0(p);
                   % All parameters attached to damper, timestep, diffusion
                   % coefficient
                   Package(ct,9) = dts(q);
                   Package(ct,10) = damps(q);
                   Package(ct,11) = D(q);
                   % All parameters attached to N_Kuhns
                   Package(ct,12) = N_Kuhn(r);
                   Package(ct,13) = b(s);
                   % Distances in correct units
                   Package(ct,14) = N_Kuhn(r)*b(s)*Distances(t);
               end
              end
             end
            end
           end
          end
         end
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
                    ToggleDynamics,UnitType,LengthConversion,DamperConversion,...
                    BeadSpringOrMeso,dtFact);
            end
        else
            for m1=1:size(Package,1)
                GenerateTopology(Package(m1,:),OverrideTopologies,...
                    OverrideInputScripts,CurrentFolder,OutputFolder,...
                    ToggleDynamics,UnitType,LengthConversion,DamperConversion,...
                    BeadSpringOrMeso,dtFact);
            end
        end
    end

    %% Loop over LAMMPS sims
    if RunLAMMPS==1
        for m1=1:size(Package)
            RunSimulation(Package(m1,:),OverrideRun,ToggleDynamics,...
                LengthConversion,DamperConversion,BeadSpringOrMeso,...
                CurrentFolder,OutputFolder,LAMMPSExecutablePath);
        end
    end

    %% Post Process Data
    if CompileEnsembleData==1
        CompileData(Package,OverrideCompile,OverrideCompute,                                                            ...
            ToggleDynamics,...
            LengthConversion,DamperConversion,BeadSpringOrMeso,...
            CalculateStress,CalculateEndtoEnd,CalculateMSD,...
            CalculateBondKinetics,CalculateAlignmentt,...
            CurrentFolder,OutputFolder,NoProcessors);
    end
end

%% Post Process Data
if PostProcess==1
    PostProcessData(Package,OverridePostProcess,ToggleDynamics,...
        MakeHistos,MakeMoviesHistos,OverrideHistogramMovies,...
        OverrideMSDFitData,WRTStretchOrLength,HistoUnits,...
        LengthConversion,DamperConversion,BeadSpringOrMeso,CompareModels,...
        OverrideCompareData,CurrentFolder,OutputFolder);
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CurrentFolder = DefineCurrentFolder

CurrentFolder = pwd;    
CurrentFolder = strrep(CurrentFolder, '\', '/');

end