function SetDirAndFileNames

global FileTag TopologyFileTag BeadSpringOrMeso LinearOrLangevin...
    TopologyDir TopologyFileName Top_file ChainForceFileName ChainFileTag...
    InputFileName OutputAtom OutputBond...
    StressLdgFileName StressRlxFileName StressEqlFileName...
    AtomsEqlFileName AtomsLdgFileName AtomsRlxFileName...
    BondsEqlFileName BondsLdgFileName BondsRlxFileName...
    TurnOnDynamics ConstantsFileName FileTagCompiled...
    CurrentFolder OutputFolder OutputDataFolder CompiledOutputsFolderName...
    InputFolderName InitialNodesTag InitialNodePosFileName...
    N_Kuhn LangevinForceFileName FENEForceFileName InitialMatlabFileName...
    

InputFolderName = [CurrentFolder,'Input_Data/'];
if ~isfolder(InputFolderName)
    mkdir(InputFolderName)
end

LangevinForceFileName = [InputFolderName,'Langevinforce',num2str(N_Kuhn),'.txt'];
FENEForceFileName = [InputFolderName,'FENEforce.txt']; 

OutputDataFolder = [OutputFolder,'Data/Outputs/'];
if ~isfolder(OutputDataFolder)
    mkdir(OutputDataFolder)
end

CompiledOutputsFolderName = [OutputFolder,'Data/Compiled Outputs'];
if ~isfolder(CompiledOutputsFolderName)
    mkdir(CompiledOutputsFolderName)
end

% TopologyDir = ['M:/Research/NSF/',CurrentFolder,'Data/Initial_Topologies'];
TopologyDir = [InputFolderName,'Initial_Topologies'];
if ~isfolder(TopologyDir)
    mkdir(TopologyDir)
end

% Set constants filenames
if ~isfolder([TopologyDir,'/Matlab Topological Data'])
    mkdir([TopologyDir,'/Matlab Topological Data'])
end
ConstantsFileName = ...
    [TopologyDir,'/Matlab Topological Data/Constants',FileTagCompiled,'.txt'];

% Set topology filenames
InitialNodePosFileName = [TopologyDir,'/Positions.',InitialNodesTag];
TopologyFileName = [TopologyDir,'/Topology.',TopologyFileTag,'.txt'];
Top_file = ['Topology.',TopologyFileTag,'.txt'];
if BeadSpringOrMeso==0
%     initial_tab = erase(FileTag,'Bead.');
    initial_tab = erase(InitialNodesTag,'Bead.');
elseif BeadSpringOrMeso==1
%     initial_tab = erase(FileTag,'Meso.');
    initial_tab = erase(InitialNodesTag,'Meso.');
    if LinearOrLangevin==1
        initial_tab = erase(initial_tab,'Langevin.');
    end
end
InitialMatlabFileName = [TopologyDir,'/','Initial.',initial_tab,'.m'];

% Set force table filenames
ChainForceFileName = [InputFolderName,'Force',ChainFileTag,'.table'];

if TurnOnDynamics==1
    AddOn = '';
else
    AddOn = '.Control';
end

% .in filenames
InputFileName = ['Input_Data/Input',FileTag,'.in'];

OutputAtom = ['atoms',FileTag,'.dump'];
OutputBond = ['bonds',FileTag,'.dump'];

StressEqlFileName = ['Input_Data/','stress_eql',AddOn,FileTag,'.txt'];
StressLdgFileName = ['Input_Data/','stress_ldg',AddOn,FileTag,'.txt'];
StressRlxFileName = ['Input_Data/','stress_rlx',AddOn,FileTag,'.txt'];

AtomsEqlFileName = ['Input_Data/','atoms_eql',AddOn,FileTag,'.txt'];
AtomsLdgFileName = ['Input_Data/','atoms_ldg',AddOn,FileTag,'.txt'];
AtomsRlxFileName = ['Input_Data/','atoms_rlx',AddOn,FileTag,'.txt'];

BondsEqlFileName = ['Input_Data/','bonds_eql',AddOn,FileTag,'.txt'];
BondsLdgFileName = ['Input_Data/','bonds_ldg',AddOn,FileTag,'.txt'];
BondsRlxFileName = ['Input_Data/','bonds_rlx',AddOn,FileTag,'.txt'];

end
