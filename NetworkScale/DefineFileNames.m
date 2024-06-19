function FileNames = DefineFileNames(Controls)

global sample Np N_Kuhn N Nt Ns phi eaStar edStar dt eq_time_factor...
    Weissenberg lambda...
    file_tag topology_file_tag initial_nodes_tag compiled_file_tag...
    current_folder output_folder compiled_data_folder topology_folder...
    input_folder...
    BeadSpringOrMeso LinearOrLangevin TurnOnDynamics...
    langeving_force_filename FENE_force_filename...
    full_topology_filename partial_topology_filename...
    initial_node_pos_filename initial_matlab_filename...
    input_filename...
    OutputAtom OutputAtom_eq OutputAtom_ld OutputAtom_rlx...
    OutputBond OutputBond_eq OutputBond_ld OutputBond_rlx...
    stress_ldg_filename stress_rlx_filename stress_eq_filename...
    raw_data_filename raw_atoms_filename raw_bonds_filename...
    stress_data_filename msd_data_filename endtoend_data_filename...
    storeage_loss_filename...
    alignment_data_filename bond_kinetics_filename...
    time_stretch_data_filename clustering_data_filename...
    postprocess_movie_file_name omega


% Define file tags - unique identifiers distinguishing data from each
% parameter combo. and sample
if BeadSpringOrMeso==0
    add_on = 'Bead';
elseif BeadSpringOrMeso==1
    add_on = 'Meso';
	if LinearOrLangevin==1
		add_on = [add_on,'.Langevin'];
	end
elseif BeadSpringOrMeso==2
    add_on = 'Thry';
end
if TurnOnDynamics==0
    add_on = ['.Control.',add_on];
end

initial_nodes_tag = [add_on,'.S',num2str(sample),'.Np',num2str(Np),'.Nk',num2str(N_Kuhn),...
    '.Nx',num2str(N),'.Nt',num2str(Nt),'.Ns',num2str(Ns),'.phi',num2str(phi,'%.3f')];    
        % This uniquely defines initial tethered node positions
topology_file_tag = initial_nodes_tag;   %This uniquely defines initial topology

compiled_file_tag = ['.Np',num2str(Np),'.Nt',num2str(Nt),'.Nk',num2str(N_Kuhn),...
    '.ea',num2str(eaStar,'%.2f'),'.ed',num2str(edStar,'%.2e'),...
    '.dt',num2str(dt,'%.2e'),'.teq',num2str(eq_time_factor),...
    '.W',num2str(Weissenberg,'%.4f'),...
    '.lam',num2str(lambda),'.phi',num2str(phi,'%.3f')];

% constants_file_tag = compiled_file_tag;
file_tag = [add_on,'.S',num2str(sample),compiled_file_tag];
compiled_file_tag = [add_on,compiled_file_tag];

if Controls.RunOscillatory==1
    file_tag = [file_tag,'.omega.',num2str(omega,'%.2e')];
    compiled_file_tag = [compiled_file_tag,'.omega.',num2str(omega,'%.2e')];
end

Directories = DefineFolders(Controls);

current_folder = Directories.current_folder;
if BeadSpringOrMeso==0
    input_folder = Directories.input_folder_bs;
    output_folder = Directories.output_folder_bs;
    compiled_data_folder = Directories.compiled_data_folder_bs;
    topology_folder = Directories.topology_folder_bs;
elseif BeadSpringOrMeso==1
    input_folder = Directories.input_folder_ms;
    output_folder = Directories.output_folder_ms;
    compiled_data_folder = Directories.compiled_data_folder_ms;
    topology_folder = Directories.topology_folder_ms;
end

langeving_force_filename = [input_folder,'/Langevinforce',num2str(N_Kuhn),'.txt'];
FENE_force_filename = [input_folder,'/FENEforce.txt']; 

% Set topology filenames
initial_node_pos_filename = [topology_folder,'/Positions.',initial_nodes_tag];
partial_topology_filename = ['Topology.',topology_file_tag,'.txt']; % excludes the folder name
full_topology_filename = [topology_folder,'/',partial_topology_filename];
if BeadSpringOrMeso==0
    initial_tab = erase(initial_nodes_tag,'Bead.');
elseif BeadSpringOrMeso==1
    initial_tab = erase(initial_nodes_tag,'Meso.');
    if LinearOrLangevin==1
        initial_tab = erase(initial_tab,'Langevin.');
    end
end
initial_matlab_filename = [topology_folder,'/','Initial.',initial_tab,'.m'];

if TurnOnDynamics==1
    add_on = '';
else
    add_on = '.Control';
end


% .in filenames
input_filename = [input_folder,'/','Input',file_tag,'.in'];

% .dump filenames
OutputAtom = ['atoms',file_tag,'.dump'];
OutputBond = ['bonds',file_tag,'.dump'];

OutputAtom_eq = ['atoms',file_tag,'_eq.dump'];
OutputBond_eq = ['bonds',file_tag,'_eq.dump'];

OutputAtom_ld = ['atoms',file_tag,'_load.dump'];
OutputBond_ld = ['bonds',file_tag,'_load.dump'];

OutputAtom_rlx = ['atoms',file_tag,'_rlx.dump'];
OutputBond_rlx = ['bonds',file_tag,'_rlx.dump'];

% .txt stress output filenames
% stress_eq_filename = [input_folder,'/stress_eql',add_on,file_tag,'.txt'];
% stress_ldg_filename = [input_folder,'/stress_ldg',add_on,file_tag,'.txt'];
% stress_rlx_filename = [input_folder,'/stress_rlx',add_on,file_tag,'.txt'];
stress_eq_filename = ['stress_eql',add_on,file_tag,'.txt'];
stress_ldg_filename = ['stress_ldg',add_on,file_tag,'.txt'];
stress_rlx_filename = ['stress_rlx',add_on,file_tag,'.txt'];

% .mat structure for compiled outputs
raw_data_filename = [compiled_data_folder,'/Raw',add_on,file_tag,'.m'];
    % Includes atom, bond, and domain boundaries information
raw_atoms_filename = [compiled_data_folder,'/Raw.Atoms',add_on,file_tag,'.m'];
    % Includes atom, bond, and domain boundaries information
raw_bonds_filename = [compiled_data_folder,'/Raw.Bonds',add_on,file_tag,'.m'];
    % Includes atom, bond, and domain boundaries information
stress_data_filename = [compiled_data_folder,'/Stress',add_on,compiled_file_tag,'.m'];
    % Includes virial stress information
msd_data_filename = [compiled_data_folder,'/MSD',add_on,compiled_file_tag,'.m'];
    % Includes time, as well as diffusion information by atom type (e.g. sticker,
    % backbone, etc.) and spherical direction (i.e., tangential vs radial)
endtoend_data_filename = [compiled_data_folder,'/EndtoEnd',add_on,compiled_file_tag,'.m'];
    % Includes compiled end-to-end vectors and corresponding orientational
    % metric tensors for network
alignment_data_filename = [compiled_data_folder,'/MetricTensor',add_on,compiled_file_tag,'.m'];
    % Includes r \otimes r components for chain alignment
bond_kinetics_filename = [compiled_data_folder,'/Kinetics',add_on,compiled_file_tag,'.m'];
    % Includes compiled end-to-end vectors and corresponding orientational
    % metric tensors for network
time_stretch_data_filename = [compiled_data_folder,'/TimeAndStretch',add_on,compiled_file_tag,'.m'];
clustering_data_filename = [compiled_data_folder,'/Clustering',add_on,compiled_file_tag,'.m'];

if Controls.RunOscillatory==1
    storeage_loss_filename = [compiled_data_folder,'/DMA',add_on,compiled_file_tag,'.m'];
end

% Network-scale movie name during deformation
postprocess_movie_file_name = file_tag;

% Store everything in data structure
% Filenames to save force-extension relations
FileNames.langeving_force_filename = langeving_force_filename;
FileNames.FENE_force_filename = FENE_force_filename;

% Filenames to save topological initiation data from matlab
FileNames.initial_node_pos_filename = initial_node_pos_filename;
FileNames.full_topology_filename = full_topology_filename;
FileNames.initial_matlab_filename = initial_matlab_filename;

% Save .in filename for lammps
FileNames.input_filename = input_filename;

% Save .dump filenames from lammps
FileNames.OutputAtom = OutputAtom;
FileNames.OutputBond = OutputBond;
FileNames.OutputAtom_eq = OutputAtom_eq;
FileNames.OutputBond_eq = OutputBond_eq;
FileNames.OutputAtom_ld = OutputAtom_ld;
FileNames.OutputBond_ld = OutputBond_ld;
FileNames.OutputAtom_rlx = OutputAtom_rlx;
FileNames.OutputBond_rlx = OutputBond_rlx;

% Save stress output filenames 
FileNames.stress_eq_filename = stress_eq_filename;
FileNames.stress_ldg_filename = stress_ldg_filename;
FileNames.stress_rlx_filename = stress_rlx_filename;

% Save compiled output filenames
FileNames.raw_data_filename = raw_data_filename;
FileNames.raw_atoms_filename = raw_atoms_filename;
FileNames.raw_bonds_filename = raw_bonds_filename;
FileNames.stress_data_filename = stress_data_filename;
FileNames.msd_data_filename = msd_data_filename;
FileNames.endtoend_data_filename = endtoend_data_filename;
FileNames.alignment_data_filename = alignment_data_filename;
FileNames.bond_kinetics_filename = bond_kinetics_filename;
FileNames.time_stretch_data_filename = time_stretch_data_filename;

FileNames.postprocess_movie_file_name = postprocess_movie_file_name;

if Controls.RunChainLengthSweep==1
    if BeadSpringOrMeso==0
        % move already-generated files from original input folder to current
        MoveFilesWithSubstring(Directories.input_folder_bs_original,...
            Directories.input_folder_bs,file_tag)
        MoveFilesWithSubstring(Directories.topology_folder_bs_original,...
            Directories.topology_folder_bs,initial_nodes_tag)
        MoveFilesWithSubstring(Directories.topology_folder_bs_original,...
            Directories.topology_folder_bs,topology_file_tag)
        MoveFilesWithSubstring(Directories.topology_folder_bs_original,...
            Directories.topology_folder_bs,initial_tab)

        % move already-generated files from original output folder to current
        MoveFilesWithSubstring(Directories.output_folder_bs_original,...
            Directories.output_folder_bs,file_tag)

        % move already-generated files from original compiled folder to current
        MoveFilesWithSubstring(Directories.compiled_data_folder_bs_original,...
            Directories.compiled_data_folder_bs,file_tag)

        MoveFilesWithSubstring(Directories.compiled_data_folder_bs_original,...
            Directories.compiled_data_folder_bs,compiled_file_tag)

    elseif BeadSpringOrMeso==1
        % move already-generated files from original input folder to current
        MoveFilesWithSubstring(Directories.input_folder_ms_original,...
            Directories.input_folder_ms,file_tag)
        MoveFilesWithSubstring(Directories.topology_folder_ms_original,...
            Directories.topology_folder_ms,initial_nodes_tag)
        MoveFilesWithSubstring(Directories.topology_folder_ms_original,...
            Directories.topology_folder_ms,topology_file_tag)
        MoveFilesWithSubstring(Directories.topology_folder_ms_original,...
            Directories.topology_folder_ms,initial_tab)

        % move already-generated files from original output folder to current
        MoveFilesWithSubstring(Directories.output_folder_ms_original,...
            Directories.output_folder_ms,file_tag)

        % move already-generated files from original compiled folder to current
        MoveFilesWithSubstring(Directories.compiled_data_folder_ms_original,...
            Directories.compiled_data_folder_ms,file_tag)

        MoveFilesWithSubstring(Directories.compiled_data_folder_ms_original,...
            Directories.compiled_data_folder_ms,compiled_file_tag)

    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MoveFilesWithSubstring(source_folder,target_folder,check_string)

% Find all files containing the specified string
files = dir(fullfile(source_folder, ['*' check_string '*']));

% Move each file to the target folder
for i = 1:length(files)
    source_file = fullfile(source_folder, files(i).name);
    target_file = fullfile(target_folder, files(i).name);
    if isfile(source_file) && ~isfile(target_file)
        copyfile(source_file, target_file);
    end
end

end
