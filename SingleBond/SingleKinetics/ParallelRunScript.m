function ParallelRunScript

clear
clear global
clc
fclose all;
close all
close all hidden

%% Unit conversions
LengthConversion = 3.74e-9; % meters per unit code length
DamperConversion = 2.89e-4; % N*s/m per unit damper in code

%% Work flow controls
% Toggles - runs these processes as called out - will skip for cases where files already exists
GenerateTopologies = 1; 
RunLAMMPS = 1;
CompileEnsembleData = 1;
PostProcess = 1;

% Overrides - reruns these processes even if files already exists
OverrideTopologies = 0;     % 1 to initiate topologies
OverrideInputScripts = 0;   % 1 to rewrite the input scripts
OverrideRun = 0;            % 1 to rerun simulations
OverrideCompile = 0;        % 1 to override compilation of sim data
OverrideCompute = 0;        % 1 to override computation of stress, end-to-end vecs, MSDs, etc.
OverridePostProcess = 0;    % 1 to override generation of stress plots

% Set which outputs to calculate
CalculateStress = 0;            % 1 to compute virial stress-stretch/time data
CalculateEndtoEnd = 0;          % 1 to compile end-to-end distribution data and 
CalculateMSD = 0;               % 1 to compute the mean-square displacement (i.e. diffusion) data
CalculateBondKinetics = 1;      % 1 to compute the bond kinetics and exchange rates
CalculateAlignmentt = 0;        % 1 to compute metric tensors (r \otimes r) of chains that elucidates alignment
FitPowerOrExponent = 0;         % 0 for fitting inverse power law, 1 for exponential decay to the ka_eff relation
AppendScalingData = 1;          % 1 to append data from discrete implementation of scaling theory to plots
TimeRun = 0;                    % 1 to time how long it takes to run each of models

NoProcessors = 1;

% Set directories
CurrentFolder = DefineCurrentFolder;
OutputDrive = 'M:';                 % Set to same as current drive, if no separate drive. Recommended 2-4 TB storage
OutputFolder = strrep(CurrentFolder, 'C:', OutputDrive);
LAMMPSExecutablePath = '/mnt/c/Users/rjwag/Documents/lammps/build/lmp';


%% Input Parameters
ToggleDynamics = 1;         % 1 for dynamics 0 for control systems (no dynamics) - SET TO 1
Np = 6^3;       % Number of polymer pairs(make a cubic number) %DO NOT CHANGE
Ns = 1;         % Number of free stickers on each tether site
T = 1;          % Temperatures
L = 1;          % Nominal spacing (only matters with dynamics)
kbT = 293*1.38e-23;     % Thermal energy in Joules
ea = logspace(-2,0,3)*kbT;  % Normalized bond activation energy for association in 
                    %units of [kbT]
Ratios = 1e3;   % ratio of ka:kd
f0 = 10000;     % Force sensitivity for dissociation in Bell's model in 
                % units of [3kbT/(\sqrt(N) b)] Set very high to turn off
                % bond dissociation
b = 0.1667;     % Kuhn length
D_nom = 1e-10;  % nominal diffusion coefficient of a monomer in m^2/s
Fact = 0.5;            

%% Sweeping Parameters
NormLength = b*LengthConversion;            % m
N_Kuhn = [12 18 36];
D = D_nom*Fact.^0;                          % m^2/s
DiffusionTimeScales = (NormLength^2)./D;    % s
damps = kbT./D;                             % kg/s or N s/m
damps = damps/DamperConversion;             % this is for a single monomerq (effective damper will depend on N)
dtFact = 20;        % Sets dt as 1/dtFact*DiffusionTimeScale
dts = 1/dtFact*DiffusionTimeScales;             % NOTE THAT Diffusion Coeff, damper, and timestep are all linked
Distances = linspace(0.0125,0.5,10);        % Separation between tether sites in units of Nb
ka = 1/DiffusionTimeScales*exp(-ea/kbT);
kd = 1/DiffusionTimeScales*exp(-0.1)/Ratios;

Samples = 1;    % goal is to have enoug sampling from Np
NoParam = 10;
UnitType = 1;   % 0 for SI units, 1 for normalized (LJ in LAMMPS), 2 for nano

if TimeRun==1
    % Remove all but fastest and slowest sim cases
    ka = ka([1,end]);
    N_Kuhn = N_Kuhn(1:2:end);
    Distances = Distances([1,end]);
    N_Kuhn = N_Kuhn([1,end]);
end

%% Index pre-factor, A, for scaling theory
A_tab = [12     0.01    0.5324;
         12     0.1     0.5611;
         12     1       0.8859;
         18     0.01    0.5300;
         18     0.1     0.5501;
         18     1       0.8633;
         36     0.01    0.7075;
         36     0.1     0.7872;
         36     1       1.3493];

%% Compile Input Parameters
Perms = length(Np)*length(Ns)*length(T)*length(L)*length(ka)*...
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

if TimeRun==1
    Package = Package([3,6],:);
    T_run = zeros(3,size(Package,1));
    T_compile = zeros(3,size(Package,1));
end

% Loop over the two model types
for BeadSpringOrMeso=0:2
    % BeadSpringOrMeso:
    %   Set 0 for bead-spring model
    %   Set 1 for mesoscale model
    %   Set 2 for scaling theory model
    %% Run parfor to generate topologies
    if GenerateTopologies==1
        if NoProcessors>1
            delete(gcp('nocreate'))
            parpool(NoProcessors)
            parfor n=1:size(Package,1)
                GenerateTopology(Package(n,:),OverrideTopologies,...
                    OverrideInputScripts,CurrentFolder,OutputFolder,...
                    ToggleDynamics,UnitType,LengthConversion,DamperConversion,...
                    BeadSpringOrMeso,dtFact,A_tab);
            end
            delete(gcp('nocreate'))
        else
            for n=1:size(Package,1)
                GenerateTopology(Package(n,:),OverrideTopologies,...
                    OverrideInputScripts,CurrentFolder,OutputFolder,...
                    ToggleDynamics,UnitType,LengthConversion,DamperConversion,...
                    BeadSpringOrMeso,dtFact,A_tab);
            end
        end
    end

    %% Loop over LAMMPS sims
    if RunLAMMPS==1
        for m1=1:size(Package)
            tic
            RunSimulation(Package(m1,:),OverrideRun,ToggleDynamics,...
                LengthConversion,DamperConversion,BeadSpringOrMeso,...
                CurrentFolder,OutputFolder,LAMMPSExecutablePath);

            if TimeRun==1
                T_run(BeadSpringOrMeso+1,m1) = toc;
            else
                toc
            end
        end
    end

    %% Post Process Data
    if CompileEnsembleData==1
        if TimeRun==1
            NoProcessors = 1;
            T_compile(BeadSpringOrMeso+1,:) = ...
                CompileData(Package,OverrideCompile,OverrideCompute,ToggleDynamics,...
                LengthConversion,DamperConversion,BeadSpringOrMeso,...
                CalculateStress,CalculateEndtoEnd,CalculateMSD,...
                CalculateBondKinetics,CalculateAlignmentt,...
                CurrentFolder,OutputFolder,...
                NoProcessors);
        else
            [~] = ...
                CompileData(Package,OverrideCompile,OverrideCompute,ToggleDynamics,...
                LengthConversion,DamperConversion,BeadSpringOrMeso,...
                CalculateStress,CalculateEndtoEnd,CalculateMSD,...
                CalculateBondKinetics,CalculateAlignmentt,...
                CurrentFolder,OutputFolder,...
                NoProcessors);
        end
    end
end

if TimeRun==1
    tab = table(T_run,T_compile);
    writetable(tab,'Performance Metrics.txt')
end

%% Post Process Data
if PostProcess==1
    PostProcessData(Package,OverridePostProcess,ToggleDynamics,...
        LengthConversion,DamperConversion,BeadSpringOrMeso,...
        FitPowerOrExponent,AppendScalingData,CurrentFolder,OutputFolder);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CurrentFolder = DefineCurrentFolder

CurrentFolder = pwd;    
CurrentFolder = strrep(CurrentFolder, '\', '/');

end