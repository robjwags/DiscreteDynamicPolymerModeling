function Main

clear
clear global
clc
fclose all;
close all
close all hidden


%% Work flow controls
% Define which processes of the workflow will be conducted - will skip for 
% cases where files already exists
Controls = DefineControlSwitches(1);% Introduced here to alert user to it. 
                                    % This function is called out as needed 
                                    % thorugh the remainder of the src code 


%% Directory control
% use this function to set the path to where the bead-spring and mesoscale 
% src code and output data are stored, as well as the lammps executable 
% used to run the simulations
Directories = DefineFolders(Controls);


%% Define the sweeping parameters and package into one array, 'package'
Parameters = DefineparameterPackage(Controls);



%% Initiate topologies, run simulations, and compile data
if Controls.RunOscillatory==1 || Controls.RunLargeDeformation==1
    ModelTypes = 1;
else
    ModelTypes = 0:1;
end
% for both model types
for ModelType=ModelTypes
    % loop over model type (0 for bead-spring, 1 for mesoscale)
    % This flag is called BeadSpringOrMeso in all other src code 
    InitiateRunCompile(ModelType,Parameters);
end


%% Post process the data
if Controls.RunTimingCase~=1
    PostProcessData(Parameters,Controls,Directories,ModelTypes)
else
    PlotTimingCases(Directories,ModelTypes);
end

end
