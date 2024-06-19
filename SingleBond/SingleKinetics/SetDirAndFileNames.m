function SetDirAndFileNames

global FolderName FileTag TopologyFileTag...
    TopologyDir TopologyFileName ChainForceFileName...
    InputFileName OutputAtom OutputBond...
    ConstantsFileName...%ConstantsFileTag...
    FileTagCompiled OutputDataFolderName OutputDir...
    OutputAtom_loc OutputBond_loc

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

% Set constants filenames
ConstantsFileName = ...
    [OutputDir,'/Constants',FileTagCompiled,'.txt'];

% Set topology filenames
TopologyFileName = [InputDir,'/Topology',TopologyFileTag,'.txt'];

% Set force table filenames
ChainForceFileName = [InputDir,'/Force',TopologyFileTag,'.txt'];

% .in filenames
InputFileName= [FolderName,'/','Input',FileTag,'.in'];

OutputAtom_loc = ['atoms',FileTag,'.dump'];
OutputBond_loc = ['bonds',FileTag,'.dump'];
OutputAtom = [OutputDir,'/',OutputAtom_loc];
OutputBond = [OutputDir,'/',OutputBond_loc];

% StressEqlFileName = [OutputDir,'/stress_eql',FileTag,'.txt'];
% StressLdgFileName = [OutputDir,'/stress_ldg',FileTag,'.txt'];
% StressRlxFileName = [OutputDir,'/stress_rlx',FileTag,'.txt'];
% 
% AtomsEqlFileName = [OutputDir,'/atoms_eql',FileTag,'.txt'];
% AtomsLdgFileName = [OutputDir,'/atoms_ldg',FileTag,'.txt'];
% AtomsRlxFileName = [OutputDir,'/atoms_rlx',FileTag,'.txt'];
% 
% BondsEqlFileName = [OutputDir,'/bonds_eql',FileTag,'.txt'];
% BondsLdgFileName = [OutputDir,'/bonds_ldg',FileTag,'.txt'];
% BondsRlxFileName = [OutputDir,'/bonds_rlx',FileTag,'.txt'];

end