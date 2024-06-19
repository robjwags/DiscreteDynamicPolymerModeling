%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DefineCompiledFileNames

global FileTag RawDataFileName ForceDataFileName AlignmentDataFileName...
    FileTagCompiled OutputFolderName TimeStretchDataFileName

FolderName = [OutputFolderName,'Compiled Outputs'];
if ~isfolder(FolderName)
    mkdir(FolderName)
end

AddOn = [];

RawDataFileName = [FolderName,'/Raw',AddOn,FileTag,'.m'];
    % Includes atom, bond, and domain boundaries information
ForceDataFileName = [FolderName,'/Force',AddOn,FileTagCompiled,'.m'];
    % Includes time, stretch, and virial stress information
AlignmentDataFileName = [FolderName,'/MericTensor',AddOn,FileTagCompiled,'.m'];
    % Includes r \otimes r components for chain alignment
TimeStretchDataFileName = [FolderName,'/TimeAndStretch',AddOn,FileTagCompiled,'.m'];

end

