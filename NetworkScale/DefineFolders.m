function Directories = DefineFolders(Controls)

global output_foldername_prefix

% Returns structure "Directories" containing all the directory names

% Manages the directories 

Directories.current_folder = DefineCurrentFolder;   % Defines current folder
Directories.beadspring_folder = '2024_03_30_BS/';   % End with '/'
    % sub-directory in which updated bead-spring input data and src code 
    % is stored

Directories.mesoscale_folder = '2024_03_30_MS/';    % End with '/'
    % sub-directory in which updated mesoscale input data and src code is 
    % stored

if Controls.RunLoadingRateSweep==1
    Directories.original_beadspring_folder = Directories.beadspring_folder;
    Directories.original_mesoscale_folder = Directories.mesoscale_folder;
    Directories.beadspring_folder = ...
        ['2024_04_21_LdRtSweep/',Directories.beadspring_folder];
    Directories.mesoscale_folder = ...
        ['2024_04_21_LdRtSweep/',Directories.mesoscale_folder];   
    output_foldername_prefix = 'Loading Rate Sweep/';
end
if Controls.RunChainLengthSweep==1
    Directories.original_beadspring_folder = Directories.beadspring_folder;
    Directories.original_mesoscale_folder = Directories.mesoscale_folder;
    Directories.beadspring_folder = ...
        ['2024_03_23_LengthSweep/',Directories.beadspring_folder];
    Directories.mesoscale_folder = ...
        ['2024_03_23_LengthSweep/',Directories.mesoscale_folder];   
    output_foldername_prefix = 'Length Sweep/';
end
if Controls.RunDetachmentRateSweep==1
    Directories.original_beadspring_folder = Directories.beadspring_folder;
    Directories.original_mesoscale_folder = Directories.mesoscale_folder;
    Directories.beadspring_folder = ...
        ['2024_04_21_DtchRtSweep/',Directories.beadspring_folder];
    Directories.mesoscale_folder = ...
        ['2024_04_21_DtchRtSweep/',Directories.mesoscale_folder];   
    output_foldername_prefix = 'Detachment Rate Sweep/';
end
if Controls.RunTheLJCase==1
    Directories.beadspring_folder = ...
        ['2024_03_16_LJ_Study/',Directories.beadspring_folder];
    Directories.mesoscale_folder = ...
        ['2024_03_16_LJ_Study/',Directories.mesoscale_folder];   
    output_foldername_prefix = 'LJ Study/Without LJ/';
end
if Controls.RunOscillatory==1
    Directories.original_beadspring_folder = Directories.beadspring_folder;
    Directories.original_mesoscale_folder = Directories.mesoscale_folder;
    Directories.beadspring_folder = ...
        ['2024_05_08_Oscillatory/',Directories.beadspring_folder];
    Directories.mesoscale_folder = ...
        ['2024_05_08_Oscillatory/',Directories.mesoscale_folder];   
    output_foldername_prefix = 'Frequency Sweep/';
end
if Controls.RunLargeDeformation==1
    Directories.original_beadspring_folder = Directories.beadspring_folder;
    Directories.original_mesoscale_folder = Directories.mesoscale_folder;
    Directories.beadspring_folder = ...
        ['2024_05_08_LargeDeform/',Directories.beadspring_folder];
    Directories.mesoscale_folder = ...
        ['2024_05_08_LargeDeform/',Directories.mesoscale_folder];   
    output_foldername_prefix = 'Large Deformation/';
end
if Controls.RunTimingCase==1
    Directories.original_beadspring_folder = Directories.beadspring_folder;
    Directories.original_mesoscale_folder = Directories.mesoscale_folder;
    Directories.beadspring_folder = ...
        ['2024_05_11_RunTimes/',Directories.beadspring_folder];
    Directories.mesoscale_folder = ...
        ['2024_05_11_RunTimes/',Directories.mesoscale_folder];   
    output_foldername_prefix = 'Compare Run Times/';
end


Directories.storage_drive = 'M:/';                  % End with '/'
    % if desire is to store large data sets on a seperate internal drive,
    % set this path. Data will be stored in sub-directory at the path
    % <[storage_drive,current_folder]>. If desire is to store on the 
    % current drive, set this path to empty (i.e., <StorageDrive = [];).

Directories.lammps_exe = '/mnt/c/Users/rjwag/Documents/lammps/build/lmp';
    % path to lammps executable used to run simulations

Directories.WSL_path = '/mnt/c/Users/rjwag/Documents/';
    % path to where current directory is stored in WSL-compatible syntax

if Controls.PrepareForStampede==1
    if Controls.RunChainLengthSweep==1
        Directories.WSL_path = '/work2/99999/robwags/stampede3/Length_sweep/Input_scripts/';
    elseif Controls.RunLoadingRateSweep==1
        Directories.WSL_path = '/work2/99999/robwags/stampede3/Loading_Rate_sweep/Input_scripts/';
    elseif Controls.RunDetachmentRateSweep==1
        Directories.WSL_path = '/work2/99999/robwags/stampede3/Detachment_Rate_sweep/Input_scripts/';
    elseif Controls.RunOscillatory
        Directories.WSL_path = '/work2/99999/robwags/stampede3/Frequency_sweep/Input_scripts/';
    elseif Controls.RunLargeDeformation
        Directories.WSL_path = '/work2/99999/robwags/stampede3/Large_deform/Input_scripts/';
    elseif Controls.RunTimingCase
        Directories.WSL_path = '/work2/99999/robwags/stampede3/Compare_runtimes/Input_scripts/';
    else
        Directories.WSL_path = '/work2/99999/robwags/stampede3/Full_sweep/Input_scripts/';
    end
end

% Check drive names
CheckDrives(Directories);

% Define all important sub diretories
Directories = DefineSubDirectories(Directories,Controls);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Directories = DefineSubDirectories(Directories,Controls)

current_folder = Directories.current_folder;
storage_drive = Directories.storage_drive;
for BeadSpringOrMeso=0:1
    if BeadSpringOrMeso==0
        temp_folder = Directories.beadspring_folder;
    elseif BeadSpringOrMeso==1
        temp_folder = Directories.mesoscale_folder;
    end

    % where .in files files will be stored
    input_folder = [temp_folder,'Input_Data'];
    % will be stored
    if ~isfolder(input_folder)
        mkdir(input_folder)
    end

    % were the raw matlab topologies will be stored
    topology_folder = [temp_folder,'Input_Data/Initial_Topologies'];
    if ~isfolder(topology_folder)
        mkdir(topology_folder)
    end

    % where raw .dump files will be stored
    output_folder = [storage_drive,current_folder,temp_folder,...
        'Data/Outputs'];
    if ~isfolder(output_folder)
        mkdir(output_folder)
    end


    % where compiled and post-processed data gotten from .dump files will be
    % stored
    compiled_data_folder = [storage_drive,current_folder,temp_folder,...
        'Data/Compiled Outputs'];
    if ~isfolder(compiled_data_folder)
        mkdir(compiled_data_folder)
    end

    if BeadSpringOrMeso==0
        Directories.input_folder_bs = input_folder;
        Directories.topology_folder_bs = topology_folder;
        Directories.output_folder_bs = output_folder;
        Directories.compiled_data_folder_bs = compiled_data_folder;
    elseif BeadSpringOrMeso==1
        Directories.input_folder_ms = input_folder;
        Directories.topology_folder_ms = topology_folder;
        Directories.output_folder_ms = output_folder;
        Directories.compiled_data_folder_ms = compiled_data_folder;
    end
    
    if Controls.RunChainLengthSweep==1 || Controls.RunLoadingRateSweep==1 ...
            || Controls.RunDetachmentRateSweep==1% then also define the original folders
        if BeadSpringOrMeso==0
            temp_folder = Directories.original_beadspring_folder;
        elseif BeadSpringOrMeso==1
            temp_folder = Directories.original_mesoscale_folder;
        end

        input_folder_original = [temp_folder,'Input_Data'];
        topology_folder_original = [temp_folder,'Input_Data/Initial_Topologies'];
        output_folder_original = [storage_drive,current_folder,temp_folder,...
            'Data/Outputs'];
        compiled_data_folder_original = [storage_drive,current_folder,temp_folder,...
            'Data/Compiled Outputs'];

        if BeadSpringOrMeso==0
            Directories.input_folder_bs_original = input_folder_original;
            Directories.topology_folder_bs_original = topology_folder_original;
            Directories.output_folder_bs_original = output_folder_original;
            Directories.compiled_data_folder_bs_original = compiled_data_folder_original;
        elseif BeadSpringOrMeso==1
            Directories.input_folder_ms_original = input_folder_original;
            Directories.topology_folder_ms_original = topology_folder_original;
            Directories.output_folder_ms_original = output_folder_original;
            Directories.compiled_data_folder_ms_original = compiled_data_folder_original;
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CheckDrives(Directories)

current_folder = Directories.current_folder;
beadspring_folder = Directories.beadspring_folder;
mesoscale_folder = Directories.mesoscale_folder;
storage_drive = Directories.storage_drive;


if ~isempty(storage_drive)
    if ~isfolder([storage_drive,current_folder])
        warning(['Indicated folder on storage drive does not exists. '...
            'Either data has not been created or wrong directory has ' ...
            'been indicated.']);
        mkdir([storage_drive,current_folder])
    end
    if ~isfolder([storage_drive,current_folder,beadspring_folder])
        warning(['Folder containing bead-spring raw data not found. ' ...
            'Either data has not been created or wrong directory had ' ...
            'been indicated.']);
    end
    if ~isfolder([storage_drive,current_folder,mesoscale_folder])
        warning(['Folder containing mesoscale raw data not found. ' ...
            'Either data has not been created or wrong directory had ' ...
            'been indicated.']);
    end
    if ~isfolder(beadspring_folder)
        warning(['Folder containing bead-spring src code not found. ' ...
            'Check name of directory']);
    end
    if ~isfolder(mesoscale_folder)
        warning(['Folder containing mesoscale src code not found. ' ...
            'Check name of directory']);
    end
else
    if ~isfolder(beadspring_folder)
        warning(['Folder containing bead-spring src code and raw data ' ...
            'not found. Either data has not been created or wrong ' ...
            'directory had been indicated.']);
    end
    if ~isfolder(mesoscale_folder)
        warning(['Folder containing mesoscale raw data not found. ' ...
            'Either data has not been created or wrong directory had ' ...
            'been indicated.']);
    end
end

lmp_exe_temp = Directories.lammps_exe;
rem = extractBefore(lmp_exe_temp,'/Users');
lmp_exe_temp = ['C:/',strrep(lmp_exe_temp,rem,'')];
if ~isfile(lmp_exe_temp)
    warning(['lmp executable for running lammps simulations not found. ' ...
        'Check the path set on line 18 is correct.'])
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [current_subdir] = DefineCurrentFolder

% Defines subdirectory (below level of 'Document/<first subdirectory>) for
% data storage purposes
current_folder = pwd;
current_subdir = extractAfter(current_folder,'Documents\');
current_subdir = [strrep(current_subdir,'\','/'),'/'];

end
