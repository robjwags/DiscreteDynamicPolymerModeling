function RunSimulation(package,Parameters,Controls,Directories,BSOM,...
    OverrideRun)

global input_filename TurnOnDynamics BeadSpringOrMeso LinearOrLangevin...
    OutputAtom OutputBond...
    OutputAtom_eq OutputBond_eq...
    OutputAtom_ld OutputBond_ld...
    OutputAtom_rlx OutputBond_rlx...
    current_folder output_folder raw_data_filename ...
    Np Nt N b ka dt eaStar edStar...
    length_conversion damper_conversion...
    samples N_Kuhns phis kds eq_time_factors Weissenbergs...
    sample N_Kuhn phi kd eq_time_factor Weissenberg...
    stress_eq_filename stress_ldg_filename stress_rlx_filename

TurnOnDynamics = Controls.ToggleDynamics;
LinearOrLangevin = Controls.LinearOrLangevin;
length_conversion = Parameters.length_conversion;
damper_conversion = Parameters.damper_conversion;
BeadSpringOrMeso = BSOM;
current_folder = Directories.current_folder;

UnpackParameters(package,Parameters);

% Since RunSimulation.m is swept at a functional level above this script,
% individual values of swept parameters are the same as the ranges defined
% in UnpackParameters;
sample = samples;
eq_time_factor = eq_time_factors;
N_Kuhn = N_Kuhns;
phi = phis;
kd = kds;
Weissenberg = Weissenbergs;

[~,D,tau0,dtFact] = DefineTimeScale(b,length_conversion,damper_conversion,...
    BeadSpringOrMeso,phi,N_Kuhn);
dt = tau0/dtFact;

N = Np*Nt;      % is the number of total tethers for equilibration

eaStar = -log(ka*(b*length_conversion)^2/D);
edStar = -log(kd*(b*length_conversion)^2/D);

Directories = DefineFolders(Controls);
DefineFileNames(Controls);

% files to check - if these exists, skip simulation step unless
% OverrideRun==1
if BeadSpringOrMeso==0
    check_atoms_file = [output_folder,'/',OutputAtom];
    check_bonds_file = [output_folder,'/',OutputAtom];
elseif BeadSpringOrMeso==1
    check_atoms_file = [output_folder,'/',OutputAtom_rlx];
    check_bonds_file = [output_folder,'/',OutputAtom_rlx];
end

% check run conditions and run simulations
if ((~isfile(check_atoms_file) || ~isfile(check_bonds_file)) &&...
        ~isfile(raw_data_filename)) || OverrideRun==1
    setenv('OMP_NUM_THREADS', '2');
    if BeadSpringOrMeso==0
        command = ['wsl mpirun -np 8 ',Directories.lammps_exe,' -in ',...
            input_filename];
    else
        command = ['wsl mpirun -np 8 ',Directories.lammps_exe,' -in ',...
            input_filename];
    end
    tic; disp('Started run')
    system(command); 
    toc; disp ('Ended run')

    % Move files to appropriate ouput directory
    if isfile(OutputAtom)
        movefile(OutputAtom,output_folder)
        movefile(OutputBond,output_folder)
    end
    if isfile(OutputAtom_eq)
        movefile(OutputAtom_eq,output_folder)
        movefile(OutputBond_eq,output_folder)
    end
    if isfile(OutputAtom_ld)
        movefile(OutputAtom_ld,output_folder)
        movefile(OutputBond_ld,output_folder)
    end
    if isfile(OutputAtom_rlx)
        movefile(OutputAtom_rlx,output_folder)
        movefile(OutputBond_rlx,output_folder)
    end
    if isfile(stress_eq_filename)
        movefile(stress_eq_filename,output_folder)
    end
    if isfile(stress_ldg_filename)
        movefile(stress_ldg_filename,output_folder)
    end
    if isfile(stress_rlx_filename)
        movefile(stress_rlx_filename,output_folder)
    end
end

clear global

end


