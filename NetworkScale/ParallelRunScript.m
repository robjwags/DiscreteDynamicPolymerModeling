 function InitiateRunCompile(beadspring_folder,storage_drive,package,...
    BeadSpringOrMeso,GenerateTopologies,RunLAMMPS,CompileEnsembleData,...
    LinearOrLangevin,ToggleDynamics,UnitType,...
    length_conversion,damper_conversion,Overrides)

n_processors = 1;

% Set which outputs to calculate during computation
CalculateStress = 1;            % 1 to compute virial stress-stretch/time data
CalculateEndtoEnd = 1;          % 1 to compile end-to-end distribution data and 
CalculateMSD = 1;               % 1 to compute the mean-square displacement (i.e. diffusion) data
CalculateBondKinetics = 1;      % 1 to compute the bond kinetics and exchange rates
CalculateAlignment = 1;         % 1 to compute metric tensors (r \otimes r) of chains that elucidates alignment
ComputeNetworkScale = 1;        % 1 to compute virial stress of bead-spring at network-scale and 0 to compute at pairwise bond-scale

% Define the overrides
OverrideTopologies = Overrides.OverrideTopologies;   
OverrideInputScripts = Overrides.OverrideInputScripts;
OverrideRun = Overrides.OverrideRun;
OverrideCompile = Overrides.OverrideCompile;
OverrideCompute = Overrides.OverrideCompute;


%% Run parfor to generate topologies
if GenerateTopologies==1
    if n_processors>1
        delete(gcp('nocreate'))
        parpool(n_processors)
        parfor n=1:size(package,1)
            GenerateTopology(package(n,:),OverrideTopologies,...
                OverrideInputScripts,CurrentFolder,...
                ToggleDynamics,UnitType,BeadSpringOrMeso,...
                length_conversion,damper_conversion,LinearOrLangevin);
        end
    else
        for n=1:size(package,1)
            GenerateTopology(package(n,:),OverrideTopologies,...
                OverrideInputScripts,CurrentFolder,...
                ToggleDynamics,UnitType,BeadSpringOrMeso,...
                length_conversion,damper_conversion,LinearOrLangevin);
        end
    end
end

% Generate list of inputs
string = '';
for n=1:size(package,1)
    output = GenerateListOfInputs(package(n,:),...
        CurrentFolder,...
        ToggleDynamics,UnitType,BeadSpringOrMeso,...
        length_conversion,damper_conversion,LinearOrLangevin);
    string = [string,['"/projects/bcfy/wagner2/NetworkScale/2023_12_12/',...
        output,'" ']];
end
fid = fopen('input_list.txt','wt');
fprintf(fid, string);
fclose(fid);

%% Loop over LAMMPS sims
if RunLAMMPS==1
    wb = waitbar(0,'Running all simulations');
    for n=1:size(package)
        tic
        start_time = cputime;
        waitbar(n/size(package,1),wb,'Running all simulations')
        RunSimulation(package(n,:),OverrideRun,ToggleDynamics,...
            length_conversion,damper_conversion,BeadSpringOrMeso,...
            CurrentFolder,LinearOrLangevin);
        if TimeRun==1
            T_run(BeadSpringOrMeso+1,n) = toc;
        else
            toc
        end
        run_time = cputime-start_time;
        disp(['run time = ',num2str(run_time),' s'])
    end
    close(wb)
end

%% Post Process Data
if CompileEnsembleData==1
    %         CompileData(Package,OverrideCompile,ToggleDynamics);
    T_compile(BeadSpringOrMeso+1,:) = ...
        CompileData(package,OverrideCompile,OverrideCompute,ToggleDynamics,...
        length_conversion,damper_conversion,BeadSpringOrMeso,...
        CalculateStress,CalculateEndtoEnd,CalculateMSD,...
        CalculateBondKinetics,CalculateAlignment,CurrentFolder,...
        n_processors,ComputeNetworkScale,LinearOrLangevin);
end

if TimeRun==1
    tab = table(T_run,T_compile);
    writetable(tab,'Performance Metrics.txt')
end

end