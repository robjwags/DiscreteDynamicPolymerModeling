function InitiateRunCompile(BeadSpringOrMeso,Parameters)

Controls = DefineControlSwitches(BeadSpringOrMeso);
Directories = DefineFolders(Controls);
package = Parameters.package;

%% Run parfor to generate topologies
if Controls.GenerateTopologies==1 || Controls.RerunFailed==1
    %     Controls.RerunFailed = 0; %% DELETE
    if Controls.RerunFailed==1
        rerun_package = IdentifyFailedRuns(BeadSpringOrMeso,Controls);
        if ~isempty(rerun_package)
            Controls.OverrideInputScripts = 1;
            Controls.OverrideTopologies = 1;
            package = rerun_package;
            regenerate_combos = package(:,[1,3,4]);
            [~,indx] = unique(regenerate_combos,'rows');
            package = rerun_package(indx,:); % only regenerate for unique combos of sample, N_Kuhn, and phi
        end
    end
    if Controls.n_processors>1
        delete(gcp('nocreate'))
        parpool(n_processors)
        OverrideTopologies = Controls.OverrideTopologies;
        OverrideInputScripts = Controls.OverrideInputScripts;
        parfor n=1:size(package,1)
            GenerateTopology(package,Parameters,...
                BeadSpringOrMeso,OverrideTopologies,...
                OverrideInputScripts,Directories,n,package);
        end
    else
        wb = waitbar(0,'Generating all input files...');
        for n=1:size(package,1)
            waitbar(n/size(package,1),wb,'Generating all input files...')
            GenerateTopology(package(n,:),Parameters,...
                BeadSpringOrMeso,Controls.OverrideTopologies,...
                Controls.OverrideInputScripts,Directories,n,package);
        end
        close(wb)
    end
end

if Controls.PrepareForStampede==1
    warning(['Run procedure for this model was set up for stampede3. '...
        'Make sure pertinent data is available and properly stored on ',...
        'hard drive for successful compilation.'])
else
    %% Loop over LAMMPS sims
    if Controls.RunSimulations==1 || Controls.RerunFailed==1
        if Controls.RerunFailed==1 && ~isempty(rerun_package)
            Controls.OverrideRun = 1;
        end
        if Controls.RunTimingCase==1
            tic
        end
        wb = waitbar(0,'Running all simulations...');
        for n=1:size(package)
            waitbar(n/size(package,1),wb,'Running all simulations...')
            start_time = cputime;

            RunSimulation(package(n,:),Parameters,Controls,Directories,...
                BeadSpringOrMeso,Controls.OverrideRun);
            end_time = cputime;
        end
        % if intention is to time runs
        if Controls.RunTimingCase==1
            run_time = toc;
            cpu_run_time = end_time-start_time;
            prefix = 'run';
            SaveRunTimes(run_time,cpu_run_time,BeadSpringOrMeso,prefix);
        end
        close(wb)
    end
end
package = Parameters.package;

%% Compile Ensemble Data
if (Controls.CompileData==1 || Controls.RerunFailed==1) && Controls.RunTimingCase~=1
    if Controls.RecompileFailed==1
        rerun_package = IdentifyFailedCompiles(BeadSpringOrMeso,Controls);
        if ~isempty(rerun_package)
            Controls.OverrideCompile = 1;
            Controls.OverrideCompute = 1;
            n_reruns = size(rerun_package,1);
            n_samps = 3;
            rerun_package = repmat(rerun_package,n_samps,1); % repeat for all 3 samples
            for i=1:n_samps
                rerun_package((i-1)*n_reruns+1:i*n_reruns,1) = i;
            end
            package = unique(rerun_package,'rows');
        end
    end
    if Controls.RunTimingCase==1
        tic
        start_time = cputime;
    end
    [~] = CompileData(package,Parameters,Controls,Directories,...
        BeadSpringOrMeso,Controls.OverrideCompute);
    if Controls.RunTimingCase==1
        end_time = cputime;
        compile_time = toc;
        cpu_compile_time = end_time-start_time;
        prefix = 'compile';
        SaveRunTimes(compile_time,cpu_compile_time,BeadSpringOrMeso,prefix);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveRunTimes(run_times,cpu_run_time,BeadSpringOrMeso,prefix)

if BeadSpringOrMeso==0
    mod_type = 'bead_spring';
elseif BeadSpringOrMeso==1
    mod_type = 'mesoscale';
end
if isnan(run_times) % if no run times stored
    % Do not save
    disp(['No ',prefix,' times to record'])
else
    if sum(isnan(run_times))==0 % if all run times have been recorded succesfully
        tag = [prefix,'_times.'];
    else % gather a sampling from the successful runs but don't override
        % pre-existing data
        s = CountPreexistingTempRunTimeFiles(mod_type);
        tag = [prefix,'_times_samp_',num2str(s),'.'];
    end
    writetable(table(run_times,cpu_run_time),[tag,mod_type,'.txt'])
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = CountPreexistingTempRunTimeFiles(mod_type)

directory = pwd;                    % path of directory to be searched 
files_folders = dir(directory);     % files and folders in directory
files_in_dir = files_folders(~([files_folders.isdir]));  % Returns only the files in the directory                    
string_of_interest = 'Run_Times_samp_';
n_files = length(files_in_dir);
i=1;
check1 = zeros(n_files,1);
check2 = zeros(n_files,1);
while(i<=n_files)
    file_name = files_in_dir(i).name; % name of current file
    if contains(file_name,string_of_interest)
        check1(i) = 1;
    end
    if contains(file_name,mod_type)
        check2(i) = 1;
    end
    i = i+1;
end
check = check1+check2;
check(check~=2) = 0;
n_containing = sum(check)/2;
s = n_containing+1;

end
