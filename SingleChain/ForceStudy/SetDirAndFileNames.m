function SetDirAndFileNames

global FolderName FileTag TopologyFileTag...
    TopologyDir TopologyFileName ConstantsFileName...
    InputFileName OutputAtom OutputBond...
    OutputDataFolderName OutputDir...
    OutputAtom_loc OutputBond_loc...
    RawDataFileName ForceDataFileName...
    EndToEndDataFileName AlignmentDataFileName...
    FileTagCompiled TimeStretchDataFileName...
    
InputDir = FolderName;
OutputDir = [OutputDataFolderName,'Outputs'];
TopologyDir = [OutputDataFolderName,'Matlab Topological Data'];
if ~isfolder(InputDir)
    mkdir(InputDir)
end
if ~isfolder(OutputDir)
    mkdir(OutputDir)
end
if ~isfolder(TopologyDir)
    mkdir(TopologyDir)
end

% Set topology filenames
TopologyFileName = [FolderName,'/Topology.',TopologyFileTag,'.txt'];

% .in filename
InputFileName = [FolderName,'/','Input.',FileTag,'.in'];

% Constantss filename
ConstantsFileName = [FolderName,'/','Constants.',FileTag,'.txt'];

% Output filenames
OutputAtom_loc = ['atoms.',FileTag,'.dump'];
OutputBond_loc = ['bonds.',FileTag,'.dump'];
OutputAtom = [OutputDir,'/',OutputAtom_loc];
OutputBond = [OutputDir,'/',OutputBond_loc];


CompiledFolderName = [OutputDataFolderName,'Compiled Outputs'];
if ~isfolder(CompiledFolderName)
    mkdir(CompiledFolderName)
end

AddOn = [];

RawDataFileName = [CompiledFolderName,'/Raw.',AddOn,FileTag,'.m'];
    % Includes atom, bond, and domain boundaries information
EndToEndDataFileName = [CompiledFolderName,'/EndtoEnd.',AddOn,FileTagCompiled,'.m'];
    % Includes end-to-end vectors of all bonds
ForceDataFileName = [CompiledFolderName,'/Force.',AddOn,FileTagCompiled,'.m'];
    % Includes time, stretch, and virial stress information
AlignmentDataFileName = [CompiledFolderName,'/Alignment.',AddOn,FileTagCompiled,'.m'];
    % Includes r \otimes r components for chain alignment
TimeStretchDataFileName = [CompiledFolderName,'/TimeAndStretch.',AddOn,FileTagCompiled,'.m'];


end