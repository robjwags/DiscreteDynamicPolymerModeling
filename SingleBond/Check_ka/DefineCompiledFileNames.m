%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DefineCompiledFileNames

global FileTag RawDataFileName BondKineticsDataFileName...
    FileTagCompiled OutputFolderName 

FolderName = [OutputFolderName,'Compiled Outputs'];
if ~isfolder(FolderName)
    mkdir(FolderName)
end

AddOn = [];

RawDataFileName = [FolderName,'/Raw',AddOn,FileTag,'.m'];
    % Includes atom, bond, and domain boundaries information
    
BondKineticsDataFileName = [FolderName,'/Kinetics',AddOn,FileTagCompiled,'.m'];
    % Includes compiled end-to-end vectors and corresponding orientational
    % metric tensors for network
    
end

