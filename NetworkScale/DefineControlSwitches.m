function Controls = DefineControlSwitches(BSOM)

% This is a master script in which the toggles for all procedures are
% flipped on or off (e.g., whether to run the simulations by setting
% "Controls.RunSimulations = 1;")

% Returns structure "Controls" containing all flags that, if set to 1,
% control whether or not a process or feature is toggled on.

Controls.n_processors = 1;          % number of processors to use for 
                                    % generating initial topologies in 
                                    % parallel     



%% Define which portions of the workflow will be run
Controls.GenerateTopologies = 1;    % 1 to generate raw data for inputs to 
                                    % LAMMPS 
Controls.RunSimulations = 1;        % 1 to run the simulations
Controls.CompileData = 1;           % 1 to compile all of the data from 
                                    % the runs
Controls.PostProcess = 1;           % 1 to post process the data
Controls.RerunFailed = 0;           % Roughly 1-3% of mesoscale cases don't 
                                    % converge due to steep gradients in 
                                    % the dynamic bonds. Set this to 1 to 
                                    % automatically identify which 
                                    % mesoscale cases didn't converge, 
                                    % reinitiate the topologies for these,
                                    % and then rerun these instances. Set 
                                    % to 0 if raw output data is no longer 
                                    % stored but compiled data is and don't 
                                    % want to rerun everything.
Controls.RecompileFailed = 0;       % Same logic as rerun failed, but for 
                                    % compilation of data



%% Togle simulation condition options
Controls.LinearOrLangevin = 1;      % 1 to make bonds of mesoscale Langevin 
                                    % and 0 for Gaussian
                                    % Note: Linear springs could not 
                                    % reproduce the force-extension relation
                                    % of the bead-spring model
Controls.ToggleDynamics = 1;        % 1 for dynamics 0 for control systems 
                                    % (no dynamics)
Controls.UnitType = 1;              % do not change

%% Run the Chain length sweep
Controls.RunChainLengthSweep = 1;   % Set 1 to run a sweep over:
                                    % N_Kuhn \in {12,18,24,30,36} to ID
                                    % where (at which chain lengths) the \
                                    % discrepancy in stress responses 
                                    % occurs between models. Will
                                    % automatically set:
                                    % phi = 0.204;
                                    % W = 0.1;
                                    % kd = 0;

%% Run loading rate sweep
Controls.RunLoadingRateSweep = 0;   % Set 1 to run a sweep over:
                                    % edot \in {0.01,0.018,0.032,0.056,0.1}/tao0 to ID
                                    % where (at which loading rates) the \
                                    % discrepancy in stress responses 
                                    % occurs between models. Will
                                    % automatically set:
                                    % phi = 0.204;
                                    % N = 12;
                                    % kd = 0;

%% Run detachment rate sweep
Controls.RunDetachmentRateSweep = 0;% Set 1 to run a sweep over:
                                    % kd \in {0 0.001,0.0032,0.01,0.032,0.1}/tao0 to ID
                                    % where (at which loading rates) the \
                                    % discrepancy in stress responses 
                                    % occurs between models. Will
                                    % automatically set:
                                    % phi = 0.204;
                                    % N = 12;
                                    % kd = 0;


%% Run test case
Controls.RunTimingCase = 0;         % Set 1 to run the simple case of a 
                                    % network with N=12, phi=0.204,
                                    % kd=0.01/tau0, and epsdot = 0.01/tau0, 
                                    % using both models. Find the
                                    % wall-clock time of each and compare.
                                    % Set Np=60, Samples=1, and teq=1

%% Run Large-scale mesocale simulations
Controls.RunOscillatory = 0;        % Set 1 to run a larg scale mesoscale
                                    % model undergoing oscillatory strain
                                    % for 17 distinct frequencies

Controls.RunLargeDeformation = 0; 	% Set 1 to run a larg scale mesoscale
                                    % model undergoing 20x elongation
									


%% Prepare for a sweep in Stampede3
if BSOM==0 || Controls.RunTimingCase==1 % Prepare all bead-spring sims on Stampede3
    Controls.PrepareForStampede = 1;% Set 1 to write all .in files so 
                                    % that they call out the correct input
                                    % folders when pulling in initial
                                    % topologies using read_data
else
    Controls.PrepareForStampede = 0;
    if Controls.RunOscillatory==1 ...
		|| Controls.RunLargeDeformation==1
        Controls.PrepareForStampede = 1;
    end
end

%% Run the Bead-spring model with LJ potentials between beads
Controls.RunTheLJCase = 0;          % 1 to run a seperate MS & BS cases 
                                    % without LJ potential. By default phi 
                                    % equals 0.204 and 0.5 with LJ exclusion
                                    % (lj/cut/soft). Since the BS model was
                                    % run with LJ potential by default,
                                    % this turns off the LJ and runs 9
                                    % samples of each model type for
                                    % comparison. It then makes the
                                    % comparison in PostProcess.m

%% Make sure don't have two specialty sweeps running at once
if Controls.RunChainLengthSweep==1
    Controls.RunTheLJCase = 0;      
    Controls.RunLoadingRateSweep = 0;
    Controls.RunDetachmentRateSweep = 0;
elseif Controls.RunTheLJCase==1
    Controls.RunChainLengthSweep = 0;      
    Controls.RunLoadingRateSweep = 0; 
    Controls.RunDetachmentRateSweep = 0;
elseif Controls.RunLoadingRateSweep==1
    Controls.RunChainLengthSweep = 0;      
    Controls.RunTheLJCase = 0; 
    Controls.RunDetachmentRateSweep = 0;
elseif Controls.RunDetachmentRateSweep==1
    Controls.RunChainLengthSweep = 0;
    Controls.RunTheLJCase = 0;
    Controls.RunLoadingRateSweep = 0;
end

if Controls.RunChainLengthSweep + Controls.RunTheLJCase + ...
        Controls.RunLoadingRateSweep + Controls.RunDetachmentRateSweep + ...
        Controls.RunOscillatory + Controls.RunLargeDeformation + ...
        Controls.RunTimingCase > 1 
    error(['Trying to run more than one specialty run case at once.',...
        ' Check Lines 46 to 98 of DefineConrolSwitches.m and make sure',...
        'that no more than one toggle option is set to 1.'])
end

%% Set which outputs to calculate during compilation and computational
% analysis
Controls.CalculateStress = 1;       % 1 to compute virial stress-stretch/
                                    % time data
Controls.CalculateEndtoEnd = 1;     % 1 to compile end-to-end distribution 
                                    % data and 
Controls.CalculateMSD = 1;          % 1 to compute the mean-square 
                                    % displacement (i.e. diffusion) data
Controls.CalculateBondKinetics = 1; % 1 to compute the bond kinetics and 
                                    % exchange rates
Controls.CalculateAlignment = 1;    % 1 to compute metric tensors 
                                    % (r \otimes r) of chains that 
                                    % elucidates alignment
Controls.ComputeNetworkScale = 1;   % 1 to compute virial stress of bead-
                                    % spring at network-scale and 0 to 
                                    % compute at pairwise bond-scale
Controls.CalculateClusteringMetrics = 1;    % 1 to quantify clustering via 
                                    % measures of self-connected molecules,
                                    % neighbor-connected molecules, and
                                    % distribution of neighboring
Controls.MakeMovieDuringCompute = 0;    %1 to plot 3D deformation of sample
                                    % one network and store in a .AVI                                     
                                    % connections

%% For the large deformation cases turn off most of the calculations
if Controls.RunOscillatory==1 || Controls.RunLargeDeformation==1
    Controls.CalculateStress = 1;
    Controls.CalculateEndtoEnd = 1;
    Controls.CalculateMSD = 0;
    Controls.CalculateBondKinetics = 1;
    Controls.CalculateAlignment = 1;
    Controls.ComputeNetworkScale = 0;
    Controls.CalculateClusteringMetrics = 0;
    Controls.MakeMovieDuringCompute = 0;
end

%% Overrides - reruns these processes even if files already exists
% By default, code will skip computationally intensive steps that have 
% already been run. Set the swtiches below to 1 to override this skipping
% procedure
Controls.OverrideTopologies = 0;    % 1 to rewrite topology files for lammps
Controls.OverrideInputScripts = 0;  % 1 to rewrite the input scripts
Controls.OverrideRun = 0;           % 1 to rerun simulations
Controls.OverrideCompile = 0;       % 1 to override compilation of sim data
Controls.OverrideCompute = 0;       % 1 to override computation of stress, 
                                    % end-to-end vecs, MSDs, etc.
Controls.OverrideStorageLoss = 0;   % 1 to override computation of stress, 
                                    % end-to-end vecs, MSDs, etc.
Controls.OverrideConsolidate = 0;   % 1 to override consolidation of 
                                    % plotting data between different
                                    % equilibration times
Controls.OverridePostProcess = 1;   % 1 to override generation of stress 
                                    % plots
Controls.OverrideMovie = 0;         % 1 to rewrite video of 3D deformation

if Controls.RunTheLJCase==1
    Controls.MakeMovieDuringCompute = 0;
end


end