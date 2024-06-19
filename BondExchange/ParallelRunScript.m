function ParallelRunScript

clear
clear global
clc
fclose all;
close all
close all hidden

%% Unit conversions
LengthConversion = 3.74e-9;     % meters per unit code length
DamperConversion = 2.89e-4;     % N*s/m per unit damper in code

%% Work flow controls
% Toggles - runs these processes as called out - will skip for cases where files already exists
GenerateTopologies = 0; 
RunLAMMPS = 0;
CompileEnsembleData = 0;        % Set to 0 if compiling using stampede2
PostProcess = 1;

% Overrides - reruns these processes even if files already exists
OverrideTopologies = 0;         % 1 to initiate topologies
OverrideInputScripts = 0;       % 1 to rewrite the input scripts
OverrideRun = 0;                % 1 to rerun simulations
OverrideCompile = 0;            % 1 to override compilation of sim data
OverrideCompute = 0;            % 1 to override computation of stress, end-to-end vecs, MSDs, etc.
OverridePostProcess = 0;        % 1 to override generation of stress plots

% Set which outputs to calculate
CalculateStress = 0;            % 1 to compute virial stress-stretch/time data
CalculateEndtoEnd = 0;          % 1 to compile end-to-end distribution data and 
CalculateMSD = 0;               % 1 to compute the mean-square displacement (i.e. diffusion) data
CalculateBondKinetics = 1;      % 1 to compute the bond kinetics and exchange rates
CalculateAlignment = 0;       	% 1 to compute metric tensors (r \otimes r) of chains that elucidates alignment
TraceABond = 0;                 % 1 to trace a bond through time - particularly those who experience exchange
AppendScalingData = 0;          % 1 to append data from discrete implementation of scaling theory to plots
ToggleDynamics = 1;             % 1 for dynamics 0 for control systems (no dynamics) - LEAVE SET TO 1
LinearOrLangevin = 1;           % 1 to make bonds of mesoscale Langevin and 0 for Gaussian - LEAVE SET TO 1

% Input current director, same directory configured for WSL (for LAMMPS), 
% output director (same as current if want same hard drive and path), 
% and the path to the lmp executable
current_dir = 'C:/Users/16102/Documents/Research/NSF/2024_BondExchange/2024_01_15/';   % End with '/'
current_dir_lmp = '/mnt/c/Users/16102/Documents/Research/NSF/2024_BondExchange/2024_01_15/';
output_dir = current_dir;    % End with '/'
lmp_path = '/mnt/c/Users/16102/Documents/lammps/lammps/build/lmp';  % End with '/'

NoProcessors = 1;


%% Input Parameters
N_edge = 7;     % Must set to an even number so particles at bounds can mate accross
Np = N_edge^3;  % Number of polymer bonding sets (must make a cubic number)
                % this is now the number of polymers with adjacent binding
                % sites with which to connect
Ns = 1;         % Number of free stickers on each tether site
kbT = 293*1.38e-23;  % Thermal energy in Joules
ea = 0.01*kbT;	% Normalized bond activation energy for association in 
                % units of [kbT]
Ratios = 1e3;   % ratio of ka:kd
f0 = 10000;     % Force sensitivity for dissociation in Bell's model in 
                % units of [3kbT/(\sqrt(N) b)] Set very high to turn off
                % bond dissociation
b = 0.1667;     % Kuhn length        

%% Sweeping Parameters
N_Kuhn = [12 18 36];                        % Only checking extremes for this study

% Bond kinetics and loading rate parameters
[~,~,tau0,~] = DefineTimeScale(b,LengthConversion,DamperConversion,0); % Define timescales
ka = 1/tau0*exp(-ea/kbT);
kd = 1/tau0*exp(-0.1)./Ratios;

phi = logspace(-2,0,15);                   % packing volume of polymer assuming v_chain ~ Nb3
phi(end) = [];
phi(1:2) = [];

Samples = 1;    % goal is to have enoug sampling from Np
NoParam = 10;
UnitType = 1;   % 0 for SI units, 1 for normalized (LJ in LAMMPS), 2 for nano

%% Index pre-factor, A, for scaling theory
A_tab = [12     0.01    2.5678;
         18     0.01    2.2498;
         36     0.01    1.9322];

%% Compile Input Parameters
Perms = length(Np)*length(Ns)*length(ka)*length(kd)*...
    length(f0)*length(N_Kuhn)*length(b)*length(phi);
Package = zeros(Perms,NoParam); ct = 0;
for i=1:length(Samples)
 for j=1:length(Np)
  for k=1:length(Ns)
   for l=1:length(ka)
    for m=1:length(kd)
     for n=1:length(f0)
      for o=1:length(N_Kuhn)
       for p=1:length(b)
        for q=1:length(phi)
         	ct = ct+1;
           	v_poly = Np*N_Kuhn(o)*b(p)^3;       % polymer volume
        	v_dom = v_poly/phi(q);           	% domain volume
            Lx = (v_dom)^(1/3);
            Separation = Lx/(Np^(1/3));

            Package(ct,1) = Samples(i);
            Package(ct,2) = Np(j);
            Package(ct,3) = Ns(k);
            Package(ct,4) = ka(l);
            Package(ct,5) = kd(m);
            Package(ct,6) = f0(n);

            % All parameters attached to N_Kuhns
            Package(ct,7) = N_Kuhn(o);
            Package(ct,8) = b(p);
               
            % Distances in correct units
            Package(ct,9) = Separation;
            Package(ct,10) = phi(q);
            
            Package(ct,11) = kbT;
        end
       end
      end
     end
    end
   end
  end
 end
end

% Loop over the two model types
for BeadSpringOrMeso=0:1
    
    %% Run parfor to generate topologies
    if GenerateTopologies==1
        if NoProcessors>1
            delete(gcp('nocreate'))
            parpool(NoProcessors)
            parfor i=1:size(Package,1)
                GenerateTopology(Package(i,:),OverrideTopologies,...
                    OverrideInputScripts,current_dir,output_dir,current_dir_lmp,...
                    ToggleDynamics,UnitType,LengthConversion,DamperConversion,...
                    BeadSpringOrMeso,A_tab,LinearOrLangevin);
            end
        else
            wb = waitbar(0,'Generating all topologies and input files');
            for i=1:size(Package,1)
                waitbar(i/size(Package,1),wb,'Generating all topologies and input files');
                GenerateTopology(Package(i,:),OverrideTopologies,...
                    OverrideInputScripts,current_dir,output_dir,current_dir_lmp,...
                    ToggleDynamics,UnitType,LengthConversion,DamperConversion,...
                    BeadSpringOrMeso,A_tab,LinearOrLangevin);
            end
            close(wb)
        end
    end

    %% Loop over LAMMPS sims
    if RunLAMMPS==1
        wb = waitbar(0,'Running all packages');
        for m1=1:size(Package,1)
            waitbar(m1/size(Package,1),wb,'Running all packages');
            RunSimulation(Package(m1,:),OverrideRun,ToggleDynamics,...
                LengthConversion,DamperConversion,BeadSpringOrMeso,...
                current_dir,output_dir,current_dir_lmp,lmp_path,...
                LinearOrLangevin);
        end
        close(wb)
    end

    %% Post Process Data
    if CompileEnsembleData==1
        CompileData(Package,OverrideCompile,OverrideCompute,ToggleDynamics,...
            LengthConversion,DamperConversion,BeadSpringOrMeso,...
            CalculateStress,CalculateEndtoEnd,CalculateMSD,...
            CalculateBondKinetics,CalculateAlignment,current_dir,output_dir,...
            NoProcessors,LinearOrLangevin);
    end
end

%% Post Process Data
if PostProcess==1
    PostProcessData(Package,OverridePostProcess,ToggleDynamics,...
        LengthConversion,DamperConversion,BeadSpringOrMeso,AppendScalingData,...
        current_dir,LinearOrLangevin,TraceABond);
end

end

