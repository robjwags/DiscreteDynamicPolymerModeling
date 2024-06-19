%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DefineCompiledFileNames

global FileTag RawDataFileName StressDataFileName MSDDataFileName...
    EndToEndDataFileName BondKineticsDataFileName AlignmentDataFileName...
    FileTagCompiled OutputFolderName TimeStretchDataFileName...
    RawAtomsDataFileName RawBondsDataFileName

FolderName = [OutputFolderName,'Compiled Outputs'];
if ~isfolder(FolderName)
    mkdir(FolderName)
end

% if TurnOnDynamics==0
%     AddOn = '.Control';
% else
    AddOn = [];
% end
RawDataFileName = [FolderName,'/Raw',AddOn,FileTag,'.m'];
    % Includes atom, bond, and domain boundaries information
RawAtomsDataFileName = [FolderName,'/Raw.Atoms',AddOn,FileTag,'.m'];
    % Includes atom, bond, and domain boundaries information
RawBondsDataFileName = [FolderName,'/Raw.Bonds',AddOn,FileTag,'.m'];
    % Includes atom, bond, and domain boundaries information
StressDataFileName = [FolderName,'/Stress',AddOn,FileTagCompiled,'.m'];
    % Includes time, stretch, and virial stress information
MSDDataFileName = [FolderName,'/MSD',AddOn,FileTagCompiled,'.m'];
    % Includes time, as well as diffusion information by atom type (e.g. sticker,
    % backbone, etc.) and spherical direction (i.e., tangential vs radial)
EndToEndDataFileName = [FolderName,'/EndtoEnd',AddOn,FileTagCompiled,'.m'];
    % Includes compiled end-to-end vectors and corresponding orientational
    % metric tensors for network
AlignmentDataFileName = [FolderName,'/MericTensor',AddOn,FileTagCompiled,'.m'];
    % Includes r \otimes r components for chain alignment
BondKineticsDataFileName = [FolderName,'/Kinetics',AddOn,FileTagCompiled,'.m'];
    % Includes compiled end-to-end vectors and corresponding orientational
    % metric tensors for network
TimeStretchDataFileName = [FolderName,'/TimeAndStretch',AddOn,FileTagCompiled,'.m'];

end

