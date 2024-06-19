function compile_times = CompileData(package,Parameters,Controls,Directories,...
        BeadSpringOrMeso,OverrideCompute)

% Collects data from all LAMMPS runs and compiles it for ensemble averaging
% or collective histogram outputs

% Compiles stress-time responses
% Stress-strain responses
% End-to-end data for histogram binning
% MSDs of distal sticers and permanent crosslinks (tether sites)
% Bond dynamics data


if Controls.n_processors==1
    compile_times = LoopCompileData(package,Parameters,Controls,Directories,...
        BeadSpringOrMeso,OverrideCompute);
else % ONLY WORKS IF Samples=1;
    package_types = DefineUniquePackageTypes(package);  
%     Controls.TimeRun = 0;
    delete(gcp('nocreate'))
    parpool(Controls.n_processors)
    parfor n=1:size(package_types,1)
        [~] = LoopCompileData(package,Parameters,Controls,Directories,...
            BeadSpringOrMeso,OverrideCompute);
    end
    delete(gcp('nocreate'))
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compile_times = LoopCompileData(package,Parameters,Controls, ...
    Directories,BSOM,OverrideCompute)

% Collects data from all LAMMPS runs and compiles it for ensemble averaging
% or collective histogram outputs

% Compiles stress-time responses
% Stress-strain responses

% End-to-end data for histogram binning
% Bond dynamics data

global BeadSpringOrMeso TurnOnDynamics LinearOrLangevin...
    output_folder raw_data_filename storeage_loss_filename...
    CalculateStress CalculateEndtoEnd...
    CalculateMSD CalculateBondKinetics CalculateAlignment...
    CalculateClusteringMetrics...
    OutputAtom OutputBond...
    OutputAtom_eq OutputBond_eq...
    length_conversion damper_conversion...
    Np Nt b lambda ka...
    samples eq_time_factors N_Kuhns phis kds Weissenbergs omegas...
    sample eq_time_factor N_Kuhn phi kd Weissenberg omega...
    D tau0 dtFact dt N eaStar edStar...
    E_prime E_dbl_prime freq delta sig0...
    E_prime_mean E_dbl_prime_mean freq_mean delta_mean sig0_mean...
    E_prime_se E_dbl_prime_se freq_se delta_se sig0_se

BeadSpringOrMeso = BSOM;

% Simulation control switches
TurnOnDynamics = Controls.ToggleDynamics;
LinearOrLangevin = Controls.LinearOrLangevin;

% Control switches for which things to calculate
CalculateStress = Controls.CalculateStress;
CalculateEndtoEnd = Controls.CalculateEndtoEnd;
CalculateMSD = Controls.CalculateMSD;
CalculateBondKinetics = Controls.CalculateBondKinetics;
CalculateAlignment = Controls.CalculateAlignment;
CalculateClusteringMetrics = Controls.CalculateClusteringMetrics;

% Define directories
DefineFileNames(Controls);

% Define the parameters
UnpackParameters(package,Parameters);

simulation_ct = 0;
if Controls.n_processors==1
    total_sims = size(package,1);
    wb = waitbar(0,'Compiling all data & computing outputs...');
    compile_times = NaN*ones(size(package,1),1);
else
    sim_type_no = package(end-1);
    n_package_types = package(end);
    disp(['Parameter combo. ',num2str(sim_type_no),' of ',num2str(n_package_types),' running']);
    compile_times = [];
end

%% Unpack Swept Input Parameters
if Controls.RunOscillatory~=1
    for i=1:length(eq_time_factors)
     for j=1:length(N_Kuhns)
      for k=1:length(phis)
       for l=1:length(kds)
        for m=1:length(Weissenbergs)
            eq_time_factor = eq_time_factors(i);
            N_Kuhn = N_Kuhns(j);
            phi = phis(k);
            kd = kds(l);
            Weissenberg = Weissenbergs(m);
    
            % Define timescales and binding energy scales
            [~,D,tau0,dtFact] = DefineTimeScale(b,length_conversion,damper_conversion,...
                BeadSpringOrMeso,phi,N_Kuhn);
            dt = tau0/dtFact;
    
            N = Np*Nt;      % is the number of total tethers for equilibration
    
            eaStar = -log(ka*(b*length_conversion)^2/D);
            edStar = -log(kd*(b*length_conversion)^2/D);
    
            package_temp = package(ismember(package(:,2),eq_time_factor),:);
            package_temp = package_temp(ismember(package_temp(:,3),N_Kuhn),:);
            package_temp = package_temp(ismember(package_temp(:,4),phi),:);
            package_temp = package_temp(ismember(package_temp(:,5),kd),:);
            package_temp = package_temp(ismember(package_temp(:,6),Weissenberg),:);
    
            if ~isempty(package_temp) && size(package_temp,1)~=length(samples)
                warning('Not enough package outputs for number of samples')
            end
                    
            % Compile the Raw Data
            ct = 0;
            for pkg=1:size(package_temp,1)
                ct = ct+1;
                simulation_ct = simulation_ct+1;
                if Controls.n_processors==1
                    pcnt_complete = simulation_ct/total_sims;
                    waitbar(pcnt_complete,wb,'Compiling all data & computing outputs...')
                    start_time = cputime;
                end
    
                %% Unpack Swept Input Parameters
                sample = package_temp(pkg,1);    %Sample index
    
                %% Make diriectories and set filenames
                DefineFolders(Controls);
                DefineFileNames(Controls);
    
                %% Size the domain to know the volume
                SizeTheDomain;
    
                if ~isfile(raw_data_filename) || Controls.OverrideCompile
                    %% Extract data from atoms.dump and bonds.dump files
                    tic
                    if BeadSpringOrMeso==0 || Controls.RunOscillatory==1 ||...
                            Controls.RunLargeDeformation==1
                        parse_atoms_file = [output_folder,'/',OutputAtom];
                        parse_bonds_file = [output_folder,'/',OutputBond];
                    elseif BeadSpringOrMeso==1
                        parse_atoms_file = [output_folder,'/',OutputAtom_eq];
                        parse_bonds_file = [output_folder,'/',OutputBond_eq];
                    end
    
                    ParseDumpFiles(parse_atoms_file,parse_bonds_file,dt,tau0,...
                        lambda,Weissenberg,kd,ka,Controls);
                    toc
                end
    %             if Controls.n_processors==1 && Controls.TimeRun==1
    %                 end_time = cputime;
    %                 if end_time-start_time>5  % if compilation was conducted
    %                     compile_times(1,ct) = end_time-start_time;
    %                 end
    %             end
            end
    
            if ~isempty(package_temp)
                %% Check if need to compute outputs
                compute = CheckIfNeedToComputeData(OverrideCompute,Controls);
                if compute
                    sample = 1;
    
                    %% set directories and filenames
                    DefineFolders(Controls);
                    DefineFileNames(Controls);
    
                    dat = load(raw_data_filename,'-mat','timesteps','Corners');
                    timesteps = dat.timesteps;
                    Corners = dat.Corners;
                    n_steps = size(timesteps,1);
                    if n_steps~=length(Corners)
                        n_steps = length(Corners);
                    end
                    n_samples = size(package_temp,1);
                    %                 end
    
                    if Controls.RunOscillatory==1
                        CountSteps(package_temp,Controls);
                    end
                    PreallocateOutputs(n_steps,n_samples,Np);
    
                    %% Compute outputs
    %                 samples = 1:size(package_temp,1);
                    for samp=1:size(package_temp,1)
                        %% Unpack Swept Input Parameters
                        sample = package_temp(samp,1);    %Sample index
    
    
    
                        %% set directories and filenames
                        DefineFolders(Controls);
                        DefineFileNames(Controls);
    
                        %% Calculate outputs
                        tic
                        ComputeOutputs(OverrideCompute,samp,...
                            Directories,Controls);
                        toc
                    end
                    SaveComputedOutputs(OverrideCompute,Controls);
                end
                if Controls.RunOscillatory==1
                    if ~isfile(storeage_loss_filename) || Controls.OverrideStorageLoss==1
                        FitStorageAndLoss(Controls,Directories,kd);
                    end
                end
            end
        end
       end
      end
     end
    end
else
    E_prime_mean = zeros(length(omegas),1);
    E_dbl_prime_mean = zeros(length(omegas),1);
    freq_mean = zeros(length(omegas),1);
    delta_mean = zeros(length(omegas),1);
    sig0_mean = zeros(length(omegas),1);

    E_prime_se = zeros(length(omegas),1);
    E_dbl_prime_se = zeros(length(omegas),1);
    freq_se = zeros(length(omegas),1);
    delta_se = zeros(length(omegas),1);
    sig0_se = zeros(length(omegas),1);

    for i=1:length(eq_time_factors)
     for j=1:length(N_Kuhns)
      for k=1:length(phis)
       for l=1:length(kds)
        for m=1:length(Weissenbergs)
         for n=1:length(omegas)
            eq_time_factor = eq_time_factors(i);
            N_Kuhn = N_Kuhns(j);
            phi = phis(k);
            kd = kds(l);
            Weissenberg = Weissenbergs(m);
            omega = omegas(n);
    
            % Define timescales and binding energy scales
            [~,D,tau0,dtFact] = DefineTimeScale(b,length_conversion,damper_conversion,...
                BeadSpringOrMeso,phi,N_Kuhn);
            dt = tau0/dtFact;
    
            N = Np*Nt;      % is the number of total tethers for equilibration
    
            eaStar = -log(ka*(b*length_conversion)^2/D);
            edStar = -log(kd*(b*length_conversion)^2/D);
    
            package_temp = package(ismember(package(:,2),eq_time_factor),:);
            package_temp = package_temp(ismember(package_temp(:,3),N_Kuhn),:);
            package_temp = package_temp(ismember(package_temp(:,4),phi),:);
            package_temp = package_temp(ismember(package_temp(:,5),kd),:);
            package_temp = package_temp(ismember(package_temp(:,6),Weissenberg),:);
            if Controls.RunOscillatory==1
                package_temp = package_temp(ismember(package_temp(:,9),omega),:);
            end
    
            if ~isempty(package_temp) && size(package_temp,1)~=length(samples)
                warning('Not enough package outputs for number of samples')
            end
                    
            % Compile the Raw Data
            ct = 0;
            for pkg=1:size(package_temp,1)
                ct = ct+1;
                simulation_ct = simulation_ct+1;
                if Controls.n_processors==1
                    pcnt_complete = simulation_ct/total_sims;
                    waitbar(pcnt_complete,wb,'Compiling all data & computing outputs...')
                    start_time = cputime;
                end
    
                %% Unpack Swept Input Parameters
                sample = package_temp(pkg,1);    %Sample index
    
                %% Make diriectories and set filenames
                DefineFolders(Controls);
                DefineFileNames(Controls);
    
                %% Size the domain to know the volume
                SizeTheDomain;
    
                if ~isfile(raw_data_filename) || Controls.OverrideCompile
                    %% Extract data from atoms.dump and bonds.dump files
                    tic
                    if BeadSpringOrMeso==0 || Controls.RunOscillatory==1 ||...
                            Controls.RunLargeDeformation==1
                        parse_atoms_file = [output_folder,'/',OutputAtom];
                        parse_bonds_file = [output_folder,'/',OutputBond];
                    elseif BeadSpringOrMeso==1
                        parse_atoms_file = [output_folder,'/',OutputAtom_eq];
                        parse_bonds_file = [output_folder,'/',OutputBond_eq];
                    end
    
                    ParseDumpFiles(parse_atoms_file,parse_bonds_file,dt,tau0,...
                        lambda,Weissenberg,kd,ka,Controls);
                    toc
                end
            end
    
            if ~isempty(package_temp)
                %% Check if need to compute outputs
                compute = CheckIfNeedToComputeData(OverrideCompute,Controls);
                if compute
                    sample = 1;
    
                    %% set directories and filenames
                    DefineFolders(Controls);
                    DefineFileNames(Controls);
    
                    dat = load(raw_data_filename,'-mat','timesteps','Corners');
                    timesteps = dat.timesteps;
                    Corners = dat.Corners;
                    n_steps = size(timesteps,1);
                    if n_steps~=length(Corners)
                        n_steps = length(Corners);
                    end
                    n_samples = size(package_temp,1);
                    %                 end
    
                    if Controls.RunOscillatory==1
                        CountSteps(package_temp,Controls);
                    end
                    PreallocateOutputs(n_steps,n_samples,Np);
    
                    %% Compute outputs
    %                 samples = 1:size(package_temp,1);
                    for samp=1:size(package_temp,1)
                        %% Unpack Swept Input Parameters
                        sample = package_temp(samp,1);    %Sample index
    
    
    
                        %% set directories and filenames
                        DefineFolders(Controls);
                        DefineFileNames(Controls);
    
                        %% Calculate outputs
                        tic
                        ComputeOutputs(OverrideCompute,samp,...
                            Directories,Controls);
                        toc
                    end
                    SaveComputedOutputs(OverrideCompute,Controls);
                end
                if Controls.RunOscillatory==1
                    if ~isfile(storeage_loss_filename) || Controls.OverrideStorageLoss==1
                            if n==1
                                E_prime_all = cell(length(omegas),1);
                                E_dbl_prime_all = cell(length(omegas),1);
                                freq = zeros(length(omegas),1);
                                freq_all = cell(length(omegas),1);
                                delta_all = cell(length(omegas),1);
                                sig0_all = cell(length(omegas),1);
                            end

                            FitStorageAndLoss(Controls,Directories,kd);

                            E_prime_all{n} = E_prime;
                            E_dbl_prime_all{n} = E_dbl_prime;
                            freq_all{n} = freq;
                            delta_all{n} = delta;
                            sig0_all{n} = sig0;
                            
                            indx = intersect(find(E_prime>0),find(E_dbl_prime>0));
                            E_prime_mean(n) = nanmean(E_prime(indx));
                            E_dbl_prime_mean(n) = nanmean(E_dbl_prime(indx));
                            freq_mean(n) = nanmean(freq(:));
                            delta_mean(n) = nanmean(delta(:));
                            sig0_mean(n) = nanmean(sig0(:));

                            E_prime_se(n) = nanstd(E_prime(indx),1)/sqrt(length(E_prime(indx)));
                            E_dbl_prime_se(n) = nanstd(E_dbl_prime(indx),1)/sqrt(length(E_dbl_prime(indx)));
                            freq_se(n) = nanstd(freq(:),1)/sqrt(length(freq(:)));
                            delta_se(n) = nanstd(delta(:),1)/sqrt(length(delta(:)));
                            sig0_se(n) = nanstd(sig0(:),1)/sqrt(length(sig0(:)));
                    end
                end
            end
         end
         if ~isfile(storeage_loss_filename) || Controls.OverrideStorageLoss==1
            SaveStorageAndLossData(E_prime_all,E_dbl_prime_all,freq_all,delta_all,sig0_all);
         end
        end
       end
      end
     end
    end
end

if Controls.n_processors==1
    close(wb)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveStorageAndLossData(E_prime_all,E_dbl_prime_all,freq_all,delta_all,sig0_all)

global storeage_loss_filename...
    E_prime_mean E_dbl_prime_mean freq_mean delta_mean sig0_mean...
    E_prime_se E_dbl_prime_se freq_se delta_se sig0_se

StorageAndLoss.E_prime_all = E_prime_all;
StorageAndLoss.E_dbl_prime_all = E_dbl_prime_all;
StorageAndLoss.freq_all = freq_all;
StorageAndLoss.delta_all = delta_all;
StorageAndLoss.sig0_all = sig0_all;

StorageAndLoss.E_prime_mean = E_prime_mean;
StorageAndLoss.E_dbl_prime_mean = E_dbl_prime_mean;
StorageAndLoss.freq_mean = freq_mean;
StorageAndLoss.delta_mean = delta_mean;
StorageAndLoss.sig0_mean = sig0_mean;

StorageAndLoss.E_prime_se = E_prime_se;
StorageAndLoss.E_dbl_prime_se = E_dbl_prime_se;
StorageAndLoss.freq_se = freq_se;
StorageAndLoss.delta_se = delta_se;
StorageAndLoss.sig0_se = sig0_se;

save(storeage_loss_filename,'-struct','StorageAndLoss')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FitStorageAndLoss(Controls,Directories,kd0)

global time stretch sig11 sig22 sig12...
    time_stretch_data_filename stress_data_filename storeage_loss_filename...
    E_prime E_dbl_prime freq output_foldername_prefix...
    N_Kuhn Np phi lambda force_conversion length_conversion

Stress = load(stress_data_filename,'-mat');
sig11 = Stress.sig11;
sig12 = Stress.sig12;
sig22 = Stress.sig22;

samples = size(sig11,2);

TimeStretch = load(time_stretch_data_filename,'-mat');
time = TimeStretch.time;
stretch = TimeStretch.stretch;

% Compute max in-plane shear stress
tau12 = (sig11-sig22)/2*sin(2*pi/4);% + sig12*cos(2*pi/4);
eps11 = stretch-1;
eps22 = -eps11;
gam12 = (eps11-eps22)*sin(2*pi/4);%*cos(2*pi/4);

for i=1:samples
    sample = i;
    ComputeStorageAndLossModuli(tau12(:,sample),gam12,time,samples,sample);
%     ComputeStorageAndLossModuli(sig11(:,sample),eps11,time,samples,sample);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CountSteps(package_temp,Controls)

global raw_data_filename sample n_success_all

n_success = zeros(size(package_temp,1),1);

for samp=1:size(package_temp,1)
    sample = package_temp(samp,1);    %Sample index

    %% set directories and filenames
    DefineFolders(Controls);
    DefineFileNames(Controls);

    dat = load(raw_data_filename,'-mat','timesteps','Corners');
    timesteps = dat.timesteps;
%     if BeadSpringOrMeso==1
%         dat = load(raw_data_filename,'-mat','time','Corners');
%         time{samp} = dat.time;
%     end
    n_success(samp) = length(timesteps);
end

n_success_all = min(n_success);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SizeTheDomain

global Lx Ly Lz MaxDist MinDist...
    phi b Np Nt  N_Kuhn sigmaREP CutoffPartREP

v_bead = pi/6*b^3;
N_backbone = Nt+(Nt-1)*(N_Kuhn-1);  %number of backbone beads per chain
N_sidegroups = Nt*N_Kuhn;           %number of sidegroup beads per chain
N_beads = Np*(N_backbone+N_sidegroups);
v_poly = N_beads*v_bead;
v_domain = v_poly/phi;

% Redefine nominal spacing for repulsion forces
tether_conc = Nt*Np/v_domain;
sigmaREP = 1.2*(tether_conc)^(-1/3);  %This will set the repulsion lengths scale for homogenization
CutoffPartREP = sigmaREP;

% Redefine nominal spacing for initial positioning
MaxDist = (tether_conc)^(-1/3);
MinDist = 0.7*(tether_conc)^(-1/3);

Lx = v_domain^(1/3);
Ly = v_domain^(1/3);
Lz = v_domain^(1/3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function package_types = DefineUniquePackageTypes(package)

% this algorithm identies all of the unique combinations of actual sweeping
% parameters (as opposed to sample numbers) so that they may be looped over
% in parallel and compiled into structures containing the data for each
% sample in one place

package_types = zeros(1e3,size(package,2));

eq_time_factors = unique(package(:,2));        
N_Kuhns = unique(package(:,3));
phis = unique(package(:,4));
kds = unique(package(:,5));        
Weissenbergs = unique(package(:,6));        

ct = 0;
for i=1:length(eq_time_factors)
 for j=1:length(N_Kuhns)
  for k=1:length(phis)
   for l=1:length(kds)
    for m=1:length(Weissenbergs)
        package_temp = ...
            package(ismember(package(:,2),eq_time_factors(i)),:);
        package_temp = ...
            package_temp(ismember(package_temp(:,3),N_Kuhns(j)),:);
        package_temp = ...
            package_temp(ismember(package_temp(:,4),phis(k)),:);
        package_temp = ...
            package_temp(ismember(package_temp(:,5),kds(l)),:);
        package_temp = ...
            package_temp(ismember(package_temp(:,6),Weissenbergs(m)),:);

        if ~isempty(package_temp)
            ct = ct+1;
            package_types(ct,:) = package_temp;
        end
    end
   end
  end
 end
end

package_types(package_types(:,1)==0,:) = [];
package_types(:,end+1) = (1:size(package_types,1))';
package_types(:,end+1) = size(package_types,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveComputedOutputs(OverrideCompute,Controls)

global CalculateStress CalculateEndtoEnd...
    CalculateMSD CalculateBondKinetics CalculateAlignment...
    stress_data_filename msd_data_filename...
    endtoend_data_filename bond_kinetics_filename alignment_data_filename...
    clustering_data_filename...
    time_stretch_data_filename...
    time stretch V sig11 sig22 sig33 sig12 sig23 sig31...
    sig11_vir sig22_vir sig33_vir sig12_vir sig23_vir sig31_vir...
    sig11_lmp sig22_lmp sig33_lmp...
    bondconc...
    msd_st msd_st_err...
    msd_st_r0 msd_st_r0_err...
    msd_st_r00 msd_st_r00_err...
    msd_st_t0 msd_st_t0_err...
    msd_st_t00 msd_st_t00_err...
    msd_tether msd_err_tether...
    ka_out ka1_out kd_out Na Nd...
    no_attachments no_detachments...
    rx ry rz bondtypes...
    n_bonds n_to_self n_to_other c_alpha...
    rxr11_st rxr22_st rxr33_st rxr12_st rxr23_st rxr31_st...
    rxr11_st_err rxr22_st_err rxr33_st_err rxr12_st_err rxr23_st_err rxr31_st_err...
    attached_bond_lifetimes detached_bond_lifetimes renormalized_bond_lifetimes...
    open_repeat_lifetimes open_exchange_lifetimes...
    total_exchange_events...
    total_repeat_events...
    chain_concentration...
    omega_min omega_max E_prime E_dbl_prime freq...
    E_prime_all E_dbl_prime_all freq_raw...
    storeage_loss_filename % mean_exchange_rate se_exchange_rate

if (~isfile(time_stretch_data_filename) || OverrideCompute)
    TimeAndStretch.time = time;
    TimeAndStretch.stretch = stretch;
    TimeAndStretch.V = V;
    save(time_stretch_data_filename,'-struct','TimeAndStretch')
end

if (~isfile(clustering_data_filename) || OverrideCompute)
    ClusteringData.n_bonds = n_bonds;
    ClusteringData.n_to_self = n_to_self;
    ClusteringData.n_to_other = n_to_other;
    ClusteringData.c_alpha = c_alpha;
    save(clustering_data_filename,'-struct','ClusteringData')
end

if CalculateStress && (~isfile(stress_data_filename) || OverrideCompute)
    Stress.sig11 = sig11;
    Stress.sig22 = sig22;
    Stress.sig33 = sig33;
    Stress.sig12 = sig12;
    Stress.sig23 = sig23;
    Stress.sig31 = sig31;
    Stress.bondconc = bondconc;

    Stress.sig11_vir = sig11_vir;
    Stress.sig22_vir = sig22_vir;
    Stress.sig33_vir = sig33_vir;
    Stress.sig12_vir = sig12_vir;
    Stress.sig23_vir = sig23_vir;
    Stress.sig31_vir = sig31_vir;

    if ~isempty(sig11_lmp)
        Stress.sig11_lmp = sig11_lmp;
        Stress.sig22_lmp = sig22_lmp;
        Stress.sig33_lmp = sig33_lmp;
%         Stress.sig12_lmp = sig12_lmp;
%         Stress.sig23_lmp = sig23_lmp;
%         Stress.sig31_lmp = sig31_lmp;
    end

    save(stress_data_filename,'-struct','Stress')
end

if CalculateEndtoEnd && (~isfile(endtoend_data_filename) || OverrideCompute)
    EndtoEnd.rx = rx;
    EndtoEnd.ry = ry;
    EndtoEnd.rz = rz;
    EndtoEnd.bondtypes = bondtypes;

    save(endtoend_data_filename,'-struct','EndtoEnd')
end

if CalculateMSD && (~isfile(msd_data_filename) || OverrideCompute)
    MSD.msd = msd_st;
    MSD.msd_err = msd_st_err;
    MSD.msd_st_r0 = msd_st_r0;
    MSD.msd_st_r0_err = msd_st_r0_err;
    MSD.msd_st_r00 = msd_st_r00;
    MSD.msd_st_r00_err = msd_st_r00_err;
    MSD.msd_st_t0 = msd_st_t0;
    MSD.msd_st_t0_err = msd_st_t0_err;
    MSD.msd_st_t00 = msd_st_t00;
    MSD.msd_st_t00_err = msd_st_t00_err;
    MSD.msd_tether = msd_tether;
    MSD.msd_err_tether = msd_err_tether;

    save(msd_data_filename,'-struct','MSD')
end

if CalculateBondKinetics && (~isfile(bond_kinetics_filename) || OverrideCompute)
    BondKinetics.ka = ka_out;
    BondKinetics.ka1 = ka1_out;
    BondKinetics.kd = kd_out;
    BondKinetics.Na = Na;
    BondKinetics.Nd = Nd;
    BondKinetics.no_attachments = no_attachments;
    BondKinetics.no_detachments = no_detachments;
    BondKinetics.AttachedLifetimes = attached_bond_lifetimes;
    BondKinetics.DetachedLifetimes = detached_bond_lifetimes;
    BondKinetics.RenormalizedLifetimes = renormalized_bond_lifetimes;
    BondKinetics.OpenRepeatLifetimes = open_repeat_lifetimes;
    BondKinetics.OpenExchangeLifetimes = open_exchange_lifetimes;
    BondKinetics.ExchangeEvents = total_exchange_events;
    BondKinetics.RepeatEvents = total_repeat_events;
    BondKinetics.Concentration = chain_concentration;
    
    save(bond_kinetics_filename,'-struct','BondKinetics')
end

if CalculateAlignment && (~isfile(alignment_data_filename) || OverrideCompute)
    RxR.rxr11 = rxr11_st;
    RxR.rxr22 = rxr22_st;
    RxR.rxr33 = rxr33_st;
    RxR.rxr12 = rxr12_st;
    RxR.rxr23 = rxr23_st;
    RxR.rxr31 = rxr31_st;
    RxR.rxr11_err = rxr11_st_err;
    RxR.rxr22_err = rxr22_st_err;
    RxR.rxr33_err = rxr33_st_err;
    RxR.rxr12_err = rxr12_st_err;
    RxR.rxr23_err = rxr23_st_err;
    RxR.rxr31_err = rxr31_st_err;

    save(alignment_data_filename,'-struct','RxR')
end

if Controls.RunOscillatory==1
    StorageLoss.min_freq = omega_min;
    StorageLoss.max_freq = omega_max;
    StorageLoss.E_prime_all = E_prime_all;
    StorageLoss.E_dbl_prime_all = E_dbl_prime_all;
    StorageLoss.freq_raw = freq_raw;

    save(storeage_loss_filename,'-struct','StorageLoss')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ComputeOutputs(OverrideCompute,samp,Directories,Controls)

global Np N_Kuhn b kbT dt phi Weissenberg eq_time_factor edStar...
    raw_data_filename CalculateStress CalculateEndtoEnd...
    CalculateMSD CalculateBondKinetics CalculateAlignment...
    CalculateClusteringMetrics...
    stress_data_filename msd_data_filename...
    endtoend_data_filename bond_kinetics_filename alignment_data_filename...
    clustering_data_filename...
    sig11 sig22 sig33 sig12 sig23 sig31...
    sig11_vir sig22_vir sig33_vir sig12_vir sig23_vir sig31_vir...
    sig11_lmp sig22_lmp sig33_lmp...
    bondconc time stretch V...
    msd_st msd_st_err...
    msd_st_r0 msd_st_r0_err...
    msd_st_r00 msd_st_r00_err...
    msd_st_t0 msd_st_t0_err...
    msd_st_t00 msd_st_t00_err...
    msd_tether msd_err_tether...
    ka_out ka1_out kd_out Na Nd...
    rx ry rz meso_pairs...
    n_bonds n_to_self n_to_other c_alpha... 
    rxr11_st rxr22_st rxr33_st rxr12_st rxr23_st rxr31_st...
    rxr11_st_err rxr22_st_err rxr33_st_err rxr12_st_err rxr23_st_err rxr31_st_err...
    attached_bond_lifetimes detached_bond_lifetimes renormalized_bond_lifetimes...
    open_repeat_lifetimes open_exchange_lifetimes...
    total_exchange_events mean_exchange_time se_exchange_time...
    total_repeat_events mean_repeat_time se_repeat_time...
    exchange_rate mean_exchange_rate se_exchange_rate...
    Lx Ly Lz length_conversion chain_concentration...
    BeadSpringOrMeso no_attachments no_detachments...
    postprocess_movie_file_name...
    output_foldername_prefix...
    omega_min omega_max tf tau0 samples n_success_all...
    stress_ldg_filename...
    E_prime E_dbl_prime freq...
    E_prime_all E_dbl_prime_all freq_raw

% Domain data
dat = load(raw_data_filename,'-mat','timesteps','Corners');
timesteps = dat.timesteps;
if BeadSpringOrMeso==1
    dat = load(raw_data_filename,'-mat','time','Corners');
    time = dat.time;
end
Corners = dat.Corners;
meso_pairs =[];

StartMSDMeasurePct = 0.001;
MSDStartIndx = ceil(StartMSDMeasurePct/100*length(timesteps));
MSDStartStep = timesteps(MSDStartIndx);

%% Initialize dynamics variables
timestep0 = 0;
% dynamic_pairs0 = [];
all_dynamic_pairs = [];
dynamic_pairs0 = [];
N_bound0 = 0;
N_free0 = Np;
N_virgin = Np;

%% Load data
x = [];
if (CalculateEndtoEnd && (~isfile(endtoend_data_filename) || OverrideCompute))...
        || CalculateAlignment && (~isfile(alignment_data_filename) || OverrideCompute)...
        || CalculateStress && (~isfile(stress_data_filename) || OverrideCompute)...
        || CalculateBondKinetics && (~isfile(bond_kinetics_filename) || OverrideCompute)

    dat = load(raw_data_filename,'-mat','type','id','mol','x','y','z',...
        'btype','batom1','batom2','rx','ry','rz');
    types = dat.type;
    mol = dat.mol;
    batom1 = dat.batom1;
    batom2 = dat.batom2;
    id = dat.id;
    x = dat.x; y = dat.y; z = dat.z;
    rx_all = dat.rx; ry_all = dat.ry; rz_all = dat.rz;
    bond_types = dat.btype;
end
if CalculateMSD && (~isfile(msd_data_filename) || OverrideCompute)
    dat = load(raw_data_filename,'-mat','mol');
    mol = dat.mol;
    if isempty(x)
        dat = load(raw_data_filename,'-mat','x','y','z','type','id');
        x = dat.x;
        y = dat.y;
        z = dat.z;
        types = dat.type;
        id = dat.id;
    end
end
if CalculateBondKinetics && (~isfile(bond_kinetics_filename) || OverrideCompute)
    dat = load(raw_data_filename,'-mat','bid','batom1','batom2');
    p1 = dat.batom1;
    p2 = dat.batom2;
    bond_id = dat.bid;
end
clear dat

if CalculateStress && (~isfile(stress_data_filename) || OverrideCompute)
    ld_stress_filename = [Directories.output_folder_ms,'/',stress_ldg_filename];
    if isfile(ld_stress_filename)
        lmp_stress_dat = dlmread(ld_stress_filename,'',1,0);
        lmp_stress_dat = lmp_stress_dat(1:n_success_all,:);
        sig11_lmp(:,samp) = lmp_stress_dat(:,3);
        sig22_lmp(:,samp) = lmp_stress_dat(:,4);
        sig33_lmp(:,samp) = lmp_stress_dat(:,5);
    end
end

ColorCoding = 1;    % 1 for absolute chain length, 2 for chain length
% in x-direction, 3 for chain alignment with
% x-direction
if ColorCoding==1
    add_on = 'stretch.';
elseif ColorCoding==2
    add_on = 'dx.';
elseif ColorCoding==3
    add_on = 'x_align.';
end
making_movie = 0;
movie_dir = [output_foldername_prefix,'Deformation Videos'];
movie_name = [movie_dir,'/',add_on,postprocess_movie_file_name];
if Controls.MakeMovieDuringCompute==1 && samp==1 && eq_time_factor==1 &&...
        (~isfile([movie_name,'.avi']) || Controls.OverrideMovie==1)
    if ~isfolder(movie_dir)
        mkdir(movie_dir)
    end

    v = VideoWriter([movie_dir,'/',add_on,postprocess_movie_file_name],...
        'Uncompressed AVI');
    v.FrameRate = 20;
    i_mov = 3;
    open(v);
    making_movie = 1;
end

%% Initialize dynamics variables
wb2 = waitbar(0,'Computing stress, end-to-end, MSD, & dynamics data...');

if length(time)>length(timesteps)
    time(length(timesteps)+1:end) = [];
end

if Controls.RunOscillatory==1
    timesteps = timesteps(1:n_success_all);
    time = time(1:n_success_all);
end

for i=1:length(timesteps)
    timestep = timesteps(i);
    if i==1
        Lx0 = abs(Corners(i,3)-Corners(i,2));
    end
    if BeadSpringOrMeso==0
        time(i) = timestep*dt;
    end
    stretch(i) = abs(Corners(i,3)-Corners(i,2))/Lx0;

    Length = length_conversion*abs(Corners(i,3)-Corners(i,2));
    Width = length_conversion*abs(Corners(i,5)-Corners(i,4));
    Height = length_conversion*abs(Corners(i,7)-Corners(i,6));
    V(i) = Length*Width*Height;

    waitbar(i/length(timesteps),wb2,...
        'Computing stress, end-to-end, MSD, & dynamics data...')
    rs = []; %AtomDataTemp = []; BondDataTemp = [];


    %% Store end-to-end components & Bond types
    if CalculateEndtoEnd && (~isfile(endtoend_data_filename) || OverrideCompute)
        rs = ExtractEndToEnd(types{i},id{i},mol{i},x{i},y{i},z{i},...
            bond_types{i},batom1{i},batom2{i},...
            rx_all{i},ry_all{i},rz_all{i},Np,i,samp,Corners(i,:));
        rx(i,:,samp) = rs(:,1); 
        ry(i,:,samp) = rs(:,2);
        rz(i,:,samp) = rs(:,3);
    end

    if making_movie==1 && ~mod(i,i_mov)
        tit = ['$t$ = ',num2str(time(i)/1e-9,'%.1f'),' ns',...
            ', $\lambda$ = ',num2str(stretch(i),'%.2f'),...
            ', $N$ = ',num2str(N_Kuhn),...
            ', $W$ = ',num2str(Weissenberg,'%.2f'),...
            ', $\phi$ = ',num2str(phi,'%.1f'),...
            ', $\epsilon_d^*$ = ',num2str(edStar,'%.2f')];

        if i==i_mov
            Length0 = Length;
        end
        PlotTheNetwork(rs,rx,ry,rz,x,y,z,id,types,i,mol,bond_types,batom1,batom2,...
            rx_all,ry_all,rz_all,Corners,samp,Length,Width,Height,Length0,tit,ColorCoding);
        frame = getframe(gcf);
        writeVideo(v,frame)

    end
    
    %% Store end-to-end components & Bond types
    if CalculateClusteringMetrics && ...
            (~isfile(clustering_data_filename) || OverrideCompute)
        [n_bonds(i,:,samp),n_to_self(i,:,samp),...
            n_to_other(i,:,samp),c_alpha(i,:,samp)] =...
            QuantifyClustering(id{i},mol{i},bond_types{i},...
            batom1{i},batom2{i});
    end
    
    %% Compute metric tensors by bond type
    if CalculateAlignment && (~isfile(alignment_data_filename) || OverrideCompute)
        if isempty(rs)
            rs = ExtractEndToEnd(types{i},id{i},mol{i},x{i},y{i},z{i},...
            bond_types{i},batom1{i},batom2{i},...
            rx_all{i},ry_all{i},rz_all{i},...
            Np,i,samp,Corners(i,:));
            rx(i,:,samp) = rs(:,1);
            ry(i,:,samp) = rs(:,2);
            rz(i,:,samp) = rs(:,3);
        end
        [rxr11_st(i,samp),rxr22_st(i,samp),rxr33_st(i,samp),...
            rxr12_st(i,samp),rxr23_st(i,samp),rxr31_st(i,samp),...
            rxr11_st_err(i,samp),rxr22_st_err(i,samp),rxr33_st_err(i,samp),...
            rxr12_st_err(i,samp),rxr23_st_err(i,samp),rxr31_st_err(i,samp)] = ...
            ComputeMetricTensor(rx(i,:,samp),ry(i,:,samp),rz(i,:,samp));
    end

    %% Compute stress
    if CalculateStress && (~isfile(stress_data_filename) || OverrideCompute)
        if isempty(rs)          
            rs = ExtractEndToEnd(types{i},id{i},mol{i},x{i},y{i},z{i},...
                bond_types{i},batom1{i},batom2{i},...
                rx_all{i},ry_all{i},rz_all{i},Np,i,samp,Corners(i,:));
            rx(i,:,samp) = rs(:,1);
            ry(i,:,samp) = rs(:,2);
            rz(i,:,samp) = rs(:,3);
        end
        forces = ComputeForces(rs,N_Kuhn,b,kbT);

        [sigma,sigma_vir,bondconc(i,samp)] = ComputeStressAndBondConc(forces,rs,Corners(i,:),...
            rx_all{i},ry_all{i},rz_all{i},bond_types{i});
        sig11(i,samp) = sigma(1);
        sig22(i,samp) = sigma(2);
        sig33(i,samp) = sigma(3);
        sig12(i,samp) = sigma(4);
        sig23(i,samp) = sigma(5);
        sig31(i,samp) = sigma(6);
        
        if BeadSpringOrMeso==0
            sig11_vir(i,samp) = sigma_vir(1);
            sig22_vir(i,samp) = sigma_vir(2);
            sig33_vir(i,samp) = sigma_vir(3);
            sig12_vir(i,samp) = sigma_vir(4);
            sig23_vir(i,samp) = sigma_vir(5);
            sig31_vir(i,samp) = sigma_vir(6);
        end
    end

    %% Compute MSDs (by atom type)
    if CalculateMSD && (~isfile(msd_data_filename) || OverrideCompute)
        if timestep<=MSDStartStep
            msd_st(i) = 0;
            msd_st_err(i) = 0;
        else
            stickertype = 2;   % End groups (MSD should follow Rouse diffusion)
            tethertype = 1;
            prev_indx = i-1;
            [msd_st(i,samp),msd_st_err(i,samp),...
                msd_st_r0(i,samp),msd_st_r0_err(i,samp),...
                msd_st_t0(i,samp),msd_st_t0_err(i,samp),...
                msd_st_r00(i,samp),msd_st_r00_err(i,samp),...
                msd_st_t00(i,samp),msd_st_t00_err(i,samp),...
                msd_tether(i,samp),msd_err_tether(i,samp)] = ...
                ComputeMSDs(x{i},y{i},z{i},...
                x{prev_indx},y{prev_indx},z{prev_indx},...
                x{MSDStartIndx},y{MSDStartIndx},z{MSDStartIndx},...
                Corners(i,:),...
                types{i},types{prev_indx},types{MSDStartIndx},...
                id{i},id{prev_indx},id{MSDStartIndx},...
                mol{i},mol{prev_indx},mol{MSDStartIndx},...
                stickertype,tethertype);
        end
    end

    %% Compute bond dynamics rates
    if CalculateBondKinetics && (~isfile(bond_kinetics_filename) || OverrideCompute)
        tiptyp = 2;
        dynbondtype = 2;
        if i==1
            dt_temp = dt;
            Dt = 0;
        else
            dt_temp = (time(i)-time(i-1))/(timesteps(i)-timesteps(i-1));
            Dt = time(i)-time(i-1);
        end
        [kd_out(i,samp),ka_out(i,samp),ka1_out(i,samp),...
            dynamic_pairs0,N_bound0,N_free0,timestep0,all_dynamic_pairs,...
            N_virgin,no_attachments(i,samp),no_detachments(i,samp)] = ...
            ComputeBondDynamics(timestep,i,dt_temp,Dt,...
            p1,p2,bond_types,bond_id,...
            types{i},dynamic_pairs0,...
            N_bound0,N_free0,...
            timestep0,tiptyp,dynbondtype,all_dynamic_pairs,N_virgin);
        Na(i,samp) = N_bound0;
        Nd(i,samp) = N_free0;
% %         %% DELETE
% %         if Weissenberg==0.1 && N_Kuhn==12 && phi==0.5 && edStar~=Inf
% % %             figure(3); 
% % %             if i==1
% % %                 clf;
% % %             end
% % %             hold on
% % %             plot(1:i,movmean(ka_out(1:i,samp),100))
% % 
% %             figure(4); 
% %             if i==1
% %                 clf;
% %             end
% %             hold on 
% %             if ~mod(i,10)
% %                 scatter(i,dt_temp)
% %             end
% % 
% %             figure(5); clf; hold on
% %             plot(1:i,movmean(no_attachments(1:i,samp),100))
% %             plot(1:i,movmean(no_detachments(1:i,samp),100))
% %         end
    end
end
close(wb2)

if making_movie==1
    close(v);
end

clear x y z mol types id rx_all ry_all rz_all p1 p2 bond_types bond_id

if phi==0.5 && N_Kuhn==12 && BeadSpringOrMeso==1 && Weissenberg==0.1
    disp('check ka')
end

% Refer to Stampede2/CompileDataStamped_v1.m for old method
if CalculateBondKinetics && (~isfile(bond_kinetics_filename) || OverrideCompute)
    attached_bond_lifetimes = [];
    detached_bond_lifetimes = [];
    open_exchange_lifetimes = [];
    open_repeat_lifetimes = [];
    renormalized_bond_lifetimes = [];

    chain_concentration = Np/(Lx*Ly*Lz)/length_conversion^3;
    
    ever_bonded = all_dynamic_pairs(:,1:2);  % all nodes that are bonded at any point in time
    ever_bonded = unique(ever_bonded(:));
    N_bonded = length(ever_bonded);
    exchange_events = zeros(length(timesteps),N_bonded);
    all_repeat_times = zeros(length(timesteps),N_bonded);
    for i=1:N_bonded  %for 1 to the number of nodes that experience bonds at any point in time
        temp_node = ever_bonded(i);
        indx1 = find(all_dynamic_pairs(:,1)==temp_node);    %find all instances in which node is first member of pair
        indx2 = find(all_dynamic_pairs(:,2)==temp_node);    %find all instances in which node is second member of pair
        all = all_dynamic_pairs([indx1;indx2],:);   %combine these instances
        all = sortrows(all,3);      %sort by time
        temp_pairs = all(:,1:2);    %isolate the pairs
        for t=1:size(temp_pairs,1)
            indx = find(temp_pairs==temp_node);
            if indx==2
                temp_pairs(t,2) = temp_pairs(t,1);
                temp_pairs(t,1) = temp_node;
            end
        end
        all(:,1:2) = temp_pairs;

        % Define attached bond lifetimes
        delta_t = diff(all(:,3));   %define difference in time between detachment event and subsequent attachment event
        delta_neighb = diff(all(:,2));
        if ~isempty(delta_t)
            attached_times = delta_t(1:2:end);  %Defines the attached bond lifetimes
            if length(delta_t)>1
                detached_times = delta_t(2:2:end);  %Defines the detached bond lifetimes
                repeat_times = detached_times(delta_neighb(2:2:end)==0);

                repeat_indices = find(delta_neighb(2:2:end)==0)*2+1;
                repeat_timesteps = all(repeat_indices,end);
                all_repeat_times(repeat_timesteps,i) = repeat_times;

                exchange_times = detached_times(delta_neighb(2:2:end)~=0);
                exchange_indices = find(delta_neighb(2:2:end)~=0)*2+1;
                exchange_timesteps = all(exchange_indices,end);
                if ~isempty(exchange_timesteps)
                    exchange_events(exchange_timesteps,i) = 1;
                end
            else
                detached_times = [];
                repeat_times = [];
                exchange_times = [];
            end




            attached_bond_lifetimes = cat(1,attached_bond_lifetimes,...
                attached_times);
            detached_bond_lifetimes = cat(1,detached_bond_lifetimes,...
                detached_times);
            open_exchange_lifetimes = cat(1,open_exchange_lifetimes,...
                exchange_times);
            open_repeat_lifetimes = cat(1,open_repeat_lifetimes,...
                repeat_times);
            
            % Define renormalized bond lifetimes
            start_t = all(1,3);
            for t=1:length(delta_neighb)
                temp = delta_neighb(t);
                if temp~=0 && ~mod(t,2)
                    renormalized_time = all(t,3)-start_t;
                    renormalized_bond_lifetimes =...
                        cat(1,renormalized_bond_lifetimes,...
                        renormalized_time);
                    next = t+1;
                    if next<=size(all,1)
                        start_t = all(next,3);
                    end
                end
            end
        end
    end
           
    total_exchange_events = sum(exchange_events,2); % exchange events wrt. time
    mean_exchange_time = mean(open_exchange_lifetimes);
    se_exchange_time = std(open_exchange_lifetimes)/sqrt(length(open_exchange_lifetimes));

    repeat_events = all_repeat_times;
    repeat_events(repeat_events~=0) = 1;
    total_repeat_events = sum(repeat_events,2);
    mean_repeat_time = mean(open_repeat_lifetimes);
    se_repeat_time = std(open_repeat_lifetimes)/sqrt(length(open_repeat_lifetimes));

    if length(total_exchange_events)<size(Na,1)
        rng = (length(total_exchange_events)+1:size(Na,1));
        total_exchange_events(rng) = zeros(length(rng),1);
    elseif size(Na,1)<length(total_exchange_events)
        rng = (size(Na,1)+1:length(total_exchange_events));
        total_exchange_events(rng) = [];
    end
    exchange_rate = total_exchange_events./Na/dt;
    mean_exchange_rate = nanmean(exchange_rate(:));
    se_exchange_rate = nanstd(exchange_rate(:))/sqrt(length(~isnan(exchange_rate(:))));
end

% Plot the raw stresses
figure(2);
if length(samples)>6
    if samp==1
        clf; hold on
    end
    colors_all = colormap('jet');
    indices_clr = round((linspace(1,length(colors_all),length(samples))))';
    color = colors_all(indices_clr(samp),:);
else
    if samp==1
        clf; hold on
        color = [0 1 0];
    elseif samp==2
        color = [0 1 0.5];
    elseif samp==3
        color = [0 1 1];
    elseif samp==4
        color = [0 0 1];
    elseif samp==5
        color = [0.5 0 1];
    elseif samp==6
        color = [1 0 0];
    end
end
p = plot(time*1e9,(sig11(:,samp)-sig11(1,samp))/1e3);
p.Color = color;
p = plot(time*1e9,(sig22(:,samp)-sig22(1,samp))/1e3);
p.Color = color*0.5;
p = plot(time*1e9,(sig12(:,samp)-sig12(1,samp))/1e3);
p.Color = color*0.25;

if BeadSpringOrMeso==0
    subdir = Directories.beadspring_folder;
elseif BeadSpringOrMeso==1
    subdir = Directories.mesoscale_folder;
end
foldername = [output_foldername_prefix,subdir,'Raw Stress'];
if ~isfolder(foldername)
    mkdir(foldername)
end

xlabel('time (ns)')
ylabel('stress (kPa)')
filename = [foldername,'/N_',num2str(N_Kuhn),...
    '.ed_',num2str(edStar,'%.2e'),...
    '.teq_',num2str(eq_time_factor,'%.2f'),...
    '.W_',num2str(Weissenberg,'%.3f'),...
    '.phi_',num2str(phi,'%.3f')];
saveas(gcf,[filename,'.png'])
saveas(gcf,[filename,'.fig'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ComputeStorageAndLossModuli(stress,strain,time,samples,sample)

global E_prime E_dbl_prime freq omega time_freq delta sig0

% Interpolate all of the data
n_pts = 1e5;
time_interp = (linspace(min(time),max(time),n_pts))';
strain = interp1(time,strain,time_interp);
stress = interp1(time,stress,time_interp);
% omega = interp1(time,omega,time_interp);
time = time_interp;

% Find all sampling starting points (everwhere strain transitions from
% negative to positive)
strain_sign = strain./abs(strain);
samp_indices = [find(strain_sign==1,1,'first');find(diff(strain_sign)==2)];%;length(strain_sign)];
if length(samp_indices)==1
    samp_indices = [samp_indices;length(strain_sign)];
elseif samp_indices(end)/length(strain_sign)<0.95
    samp_indices = [samp_indices;length(strain_sign)];
end
% samp_indices = [1;find(diff(strain_sign)==2)];

n_samps = length(samp_indices);

% Calculate the time-dependent frequency
% freq = omega/2/pi;

if sample==1
    E_prime = zeros(n_samps-1,samples);
    E_dbl_prime = zeros(n_samps-1,samples);
    freq = zeros(n_samps-1,samples);
    time_freq = zeros(n_samps-1,samples);
    delta = zeros(n_samps-1,samples);
    sig0 = zeros(n_samps-1,samples);
end
% if sample==1
%     E_prime = zeros(n_samps,samples);
%     E_dbl_prime = zeros(n_samps,samples);
%     freq = zeros(n_samps,samples);
%     time_freq = zeros(n_samps,samples);
%     delta = zeros(n_samps,samples);
%     sig0 = zeros(n_samps,samples);
% end

mean_stress = mean(stress,2);

eps0 = max(strain)-min(strain);
sig0_guess = max(mean_stress)-min(mean_stress);
delta_guess = pi/4;

% for i=2:n_samps-1
for i=1:n_samps-1
% for i=1:n_samps
    samp_indx = samp_indices(i);
    time_t = time(samp_indx);
    strain_t = strain(samp_indx);
    stress_t = mean_stress(samp_indx);
    freq_t = omega;
    time_freq(i,sample) = time_t;

    indx_low = samp_indices(i);
    indx_up = samp_indices(i+1);

    rng = (indx_low:indx_up)';
    
    figure(10); clf; 
    subplot(2,1,1); hold on
    plot(time,strain,'k')
    plot(time(rng),strain(rng),'r')
    scatter(time_t,strain_t)

    subplot(2,1,2); hold on
    plot(time,mean_stress,'k')
    plot(time(rng),mean_stress(rng),'r')
    scatter(time_t,stress_t)

    % Fit the shifted sinuosoid to the stress data
    [delta(i,sample),sig0(i,sample),freq(i,sample)]...
        = FitSinusoidFindDelta(freq_t,time(rng),strain(rng),...
        mean_stress(rng),delta_guess,sig0_guess,eps0);

    E_prime(i,sample) = sig0(i,sample)/eps0*cos(delta(i,sample));
    E_dbl_prime(i,sample) = sig0(i,sample)/eps0*sin(delta(i,sample));

    if abs(delta(i,sample))>pi/2
        disp(delta(i,sample))
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delta,sig0,freq] = FitSinusoidFindDelta(freq_guess,time,strain,...
        stress,delta_guess,sig0_guess,eps0)

global tau0

time = time-time(1);
stress = stress-mean(stress);

window_size = length(strain);
Y = fft(strain,window_size);
P2 = abs(Y/window_size);
P1 = P2(1:window_size/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = (0:(window_size/2))/window_size * (1/mean(diff(time)));

% Identify the main frequency
[~, idx_max] = max(P1);
main_frequency = f(idx_max);

sig0 = max(stress)-min(stress);

freq_nom = freq_guess;
no_fit_attempts = 100;
n_pts = 51;
R2 = 0;
ct = 0;
wb = waitbar(0,'Fitting strain...');
while R2<0.995
    ct = ct+1;
    waitbar(ct/no_fit_attempts,wb,'Fitting strain...')

    if ct==1
        sig_freq = 1/2*freq_guess;
    else
        sig_freq = sig_freq*0.95;
    end
    freq_rng = [(linspace(freq_nom-sig_freq,freq_nom+sig_freq,n_pts))';freq_nom];
    freq_rng = unique(freq_rng); 
    freq_rng(freq_rng<0) = [];
%     freq_rng(freq_rng<1/2*freq_guess) = [];

    Perms = zeros(length(freq_rng),1);
    i = 0;
    for mi=1:length(freq_rng)
        i = i+1;
        Perms(i,1) = freq_rng(mi);
    end
    R2_all = zeros(size(Perms,1),1);
    freq_all = zeros(size(Perms,1),1);

    for i=1:size(Perms,1)
        freq_all(i) = Perms(i,1);

        ft = eps0/2*sin(2*pi*freq_all(i)/tau0*time);

        RSS = sum((strain-ft).^2);
        TSS = sum((strain-mean(strain)).^2);
        R2_all(i) = 1-RSS/TSS;
    end
    difference = abs(R2_all-1);
    indx = find(difference==min(difference),1,'first');

    freq_nom = freq_all(indx);
    R2 = R2_all(indx);

    if ~mod(ct,100)
        figure(100); clf; hold on
        scatter(time,strain,'k','filled')
        plot(time,eps0/2*sin(2*pi*freq_nom/tau0*time),'r--')
        set(gcf,'color','w')
    end
    if ct>no_fit_attempts
        break;
    end
end
close(wb)

R2_strain = R2;

% fit the stress

no_fit_attempts = 10;
sig_mean = mean(stress);
delta_nom = -delta_guess;
n_pts = 21;%101;
R2 = 0;
ct = 0;
sig0_nom = sig0_guess;
wb = waitbar(0,'Fitting Exponential...');
while R2<0.995
    ct = ct+1;
    waitbar(ct/no_fit_attempts,wb,'Fitting Exponential...')

    if ct==1
        sig_delta = 1/8*delta_nom;
        delta_rng = (linspace(-pi,pi,n_pts))';
        sig_sig0 = 1/8*sig0_nom;
    else
        sig_delta = sig_delta*0.95;
        sig_sig0 = sig_sig0*0.95;
    end
    delta_rng = unique(delta_rng); 
    delta_rng(abs(delta_rng)>pi) = []; % delta for stress cannot precede delta for strain
    sig0_rng = [(linspace(sig0_nom-sig_sig0,sig0_nom+sig_sig0,n_pts))';sig0_nom];
    sig0_rng = unique(sig0_rng); 
    sig0_rng(sig0_rng<=0) = [];

    Perms = zeros(length(delta_rng),1);
    i = 0;
    for mi=1:length(delta_rng)
        for ni=1:length(sig0_rng)
            i = i+1;
            Perms(i,1) = delta_rng(mi);
            Perms(i,2) = sig0_rng(ni);
        end
    end
    R2_all = zeros(size(Perms,1),1);
    delta_all = zeros(size(Perms,1),1);
    sig0_all = zeros(size(Perms,1),1);

    for i=1:size(Perms,1)
        delta_all(i) = Perms(i,1);
        sig0_all(i) = Perms(i,2);

        ft = sig0_all(i)/2*sin(2*pi*freq_nom/tau0*time+delta_all(i))+sig_mean;

        RSS = sum((stress-ft).^2);
        TSS = sum((stress-mean(stress)).^2);
        R2_all(i) = 1-RSS/TSS;
    end
    difference = abs(R2_all-1);
    indx = find(difference==min(difference),1,'first');

    delta_nom = delta_all(indx);
    sig0_nom = sig0_all(indx);
    R2 = R2_all(indx);

    if ~mod(ct,100)
        figure(100); clf; hold on
        scatter(time,stress,'k','filled')
        plot(time,sig0_nom/2*sin(2*pi*freq_nom/tau0*time+delta_nom)+sig_mean,'r--')
        set(gcf,'color','w')
    end
    if ct>no_fit_attempts
        break;
    end
end
close(wb)
delta_stress = delta_nom;
sig0 = sig0_nom;
freq = freq_nom;
R2_stress = R2;

delta = delta_stress;

t_lag = delta/(2*pi*freq/tau0);

check_fig = 1;
if check_fig==1
    figure(100); clf
    subplot(2,1,1); hold on
    scatter(time,strain,'k','filled')
    plot(time,eps0/2*sin(2*pi*freq_nom/tau0*time),'r--')
    xlim([min(time) max(time)*1.1])

    strain_ft = eps0/2*sin(2*pi*freq_nom/tau0*time);
    tmp_rng = (1:round(length(strain_ft)));
    max_indx_strain = find(strain_ft(tmp_rng)==max(strain_ft(tmp_rng)),1,'first');

    plot([time(max_indx_strain) time(max_indx_strain)],[-eps0/2 eps0/2],'k--')

    if t_lag>0
        plot([time(max_indx_strain)-t_lag,time(max_indx_strain)],[-eps0/4,-eps0/4],'k--')
    else
        plot([time(max_indx_strain),time(max_indx_strain)-t_lag],[-eps0/4,-eps0/4],'k--')
    end

    subplot(2,1,2); hold on
    scatter(time,stress,'k','filled')
    plot(time,sig0/2*sin(2*pi*freq_nom/tau0*time+delta_stress)+sig_mean,'r--')
    xlim([min(time) max(time)*1.1])

    g1 = sig0/2*cos(delta_stress)*sin(2*pi*freq/tau0*time);
    plot(time,g1,'m--')

    g2 = sig0/2*sin(delta_stress)*cos(2*pi*freq/tau0*time);
    plot(time,g2,'b-.')

    scatter(time,g1+g2,'rx')

%     E_prime = sig0/eps0*cos(delta);
%     E_dbl_prime = sig0/eps0*sin(delta);
%     delta_calc = atan(E_dbl_prime/E_prime);
%     disp(delta)
%     disp(delta_calc)

    stress_ft = sig0/2*sin(2*pi*freq_nom/tau0*time+delta_stress)+sig_mean;
    tmp_rng = (1:round(length(stress_ft)));
    max_indx_stress = find(stress_ft(tmp_rng)==max(stress_ft(tmp_rng)),1,'first');

    plot([time(max_indx_stress) time(max_indx_stress)],[min(stress) max(stress)],'k--')
    set(gcf,'color','w')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,k,R2,delta_A,delta_k] = FitExp(A_nom,k_nom,x,y)

global NoFitAttempts

NoPts = 21;
R2 = 0;
ct = 0;
wb5 = waitbar(0,'Fitting Exponential...');
figure(100)
while R2<0.995
    ct = ct+1;
    waitbar(ct/NoFitAttempts,wb5,'Fitting Exponential...')

    if ct==1
        sig_A = 1/8*A_nom;
        sig_k = 1/8*k_nom;
    else
        sig_A = sig_A*0.95;
        sig_k = sig_k*0.95;
    end
    A_rng = [(linspace(A_nom-sig_A,A_nom+sig_A,NoPts))';A_nom];
    A_rng = unique(A_rng); 
    k_rng = [(linspace(k_nom-sig_k,k_nom+sig_k,NoPts))';k_nom];
    k_rng = unique(k_rng);

    Perms = zeros(length(A_rng)*length(k_rng),1);
    i = 0;
    for mi=1:length(A_rng)
        for bi=1:length(k_rng)
            i = i+1;
            Perms(i,1) = A_rng(mi);
            Perms(i,2) = k_rng(bi);
        end
    end
    R2_all = zeros(size(Perms,1),1);
    A_all = zeros(size(Perms,1),1);
    k_all = zeros(size(Perms,1),1);

    for i=1:size(Perms,1)
        A_all(i) = Perms(i,1);
        k_all(i) = Perms(i,2);

        ft = A_all(i)*exp(x*k_all(i));

        RSS = sum((y-ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_all(i) = 1-RSS/TSS;
    end
    diff = abs(R2_all-1);
    indx = find(diff==min(diff),1,'first');

    A_nom = A_all(indx);
    k_nom = k_all(indx);
    R2 = R2_all(indx);

    if ~mod(ct,20)
        figure(100); clf; hold on
        scatter(x,y,'k','filled')
        plot(x,A_nom*exp(x*k_nom),'k--')
    end
    if ct>NoFitAttempts
        break;
    end
end
A = A_nom;
k = k_nom;
figure(100); close
close(wb5)

find_conf_int = 1;
if find_conf_int==1 %Perturb each parameter to find range in which chi2<=1
    R2_temp = R2;
    A_temp = A;
    dA = 0.01*A;
    while R2_temp>0.05
        A_temp = A_temp+dA;
        y_ft = A_temp*exp(x*k);
        RSS = sum((y-y_ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_temp = 1-RSS/TSS;
    end
    delta_A_pos = abs(A_temp-A);
    R2_temp = R2;
    k_temp = k;
    dtau = 0.01*k;
    while R2_temp>0.05
        k_temp = k_temp+dtau;
        y_ft = A*exp(x*k_temp);
        RSS = sum((y-y_ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_temp = 1-RSS/TSS;
    end
    delta_k_pos = abs(k_temp-k);

    R2_temp = R2;
    A_temp = A;
    dA = 0.01*A;
    while R2_temp>0.05
        A_temp = A_temp-dA;
        y_ft = A_temp*exp(x/k);
        RSS = sum((y-y_ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_temp = 1-RSS/TSS;
    end
    delta_A_neg = abs(A_temp-A);
    R2_temp = R2;
    k_temp = k;
    dtau = 0.01*k;
    while R2_temp>0.05
        k_temp = k_temp-dtau;
        y_ft = A*exp(x/k_temp);
        RSS = sum((y-y_ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_temp = 1-RSS/TSS;
    end
    delta_k_neg = abs(k_temp-k);

    delta_A = (delta_A_neg + delta_A_pos)/2;
    delta_k = (delta_k_neg + delta_k_pos)/2;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function indices = FindClosestIndices(in,samps)

% Validate input
if isempty(in) || isempty(samps)
    error('Input vectors must not be empty.');
end

% Initialize the output indices array
indices = zeros(length(samps), 1);

% Loop through each sample point and find the closest index in omega
for i = 1:length(samps)
    [~, idx] = min(abs(in - samps(i)));
    indices(i) = idx;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function omega_samples = GenerateEvenOmegaSamps(omega_min,omega_max,n_samp)

% Validate input
if omega_min <= 0 || omega_max <= 0
    error('omega_min and omega_max must be positive.');
end
if omega_min >= omega_max
    error('omega_min must be less than omega_max.');
end
if n_samp < 1
    error('n_samp must be at least 1.');
end

% Generate n_samp+2 evenly spaced points in log space for omega
log_omega_samples = linspace(log(omega_min), log(omega_max), n_samp + 2);

% Convert log spaced points back to the original omega space
omega_samples = exp(log_omega_samples);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function samples = FindSamplingPointsExp(min_val,max_val,n_samp)
    % Validate input
    if min_val <= 0 || max_val <= 0
        error('omega_min and omega_max must be positive.');
    end
    if min_val >= max_val
        error('omega_min must be less than omega_max.');
    end
    if n_samp < 1
        error('n_samp must be at least 1.');
    end
    
    % Generate log-spaced samples
    samples = linspace(log(min_val),log(max_val),n_samp+2);
    
    % Exponentiate to get the samples in the original space
    samples = exp(samples);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function samples = FindSamplingPointsLog(min_val,max_val,n)

% Validate input
if min_val <= 0 || max_val <= 0
    error('omega_min and omega_max must be positive.');
end
if min_val >= max_val
    error('omega_min must be less than omega_max.');
end
if n < 1
    error('n_samp must be at least 1.');
end

% Generate log-spaced samples
samples = logspace(log10(min_val), log10(max_val), n + 2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotTheNetwork(rs,rx,ry,rz,x,y,z,id,types,i,mol,bond_types,batom1,batom2,...
    rx_all,ry_all,rz_all,Corners,samp,Length,Width,Height,Length0,tit,ColorCoding)

global N_Kuhn meso_pairs length_conversion BeadSpringOrMeso

% Define color bar inputs
if ColorCoding==1   % absolute chain length
    color_lab = '$\boldmath{r} /(\sqrt N b)$';
    cbar_ticks = [0 1/sqrt(N_Kuhn) 1];
    cbar_tick_labels = {'0','1','$\sqrt N$'};
elseif ColorCoding==2   % length in x-direction 
    color_lab = '$(\boldmath{r} \cdot \hat{\boldmath{x}})/(\sqrt N b)$';
    cbar_ticks = [0 1/sqrt(N_Kuhn) 1];
    cbar_tick_labels = {'0','1','$\sqrt N$'};
elseif ColorCoding==3   % alignment with x-direction abs(x \cdot r)
    color_lab = '$(\boldmath{r} \cdot \hat{\boldmath{x}})/ | \boldmath{r} |$';
    cbar_ticks = [0 0.5 1];
    cbar_tick_labels = {'0.0','0.5','1.0'};
end

% Extract end-to-end vecs if not already available
if isempty(rs)
    rs = ExtractEndToEnd(types{i},id{i},mol{i},x{i},y{i},z{i},...
        bond_types{i},batom1{i},batom2{i},...
        rx_all{i},ry_all{i},rz_all{i},...
        Np,i,samp,Corners(i,:));
    rx(i,:,samp) = rs(:,1);
    ry(i,:,samp) = rs(:,2);
    rz(i,:,samp) = rs(:,3);
end
x_temp = x{i};
y_temp = y{i};
z_temp = z{i};
id_temp = id{i};
type_temp = types{i};
rx_temp = rx(i,:,samp);
ry_temp = ry(i,:,samp);
rz_temp = rz(i,:,samp);

% if BeadSpringOrMeso==1 % then sort the meso_pairs by bond type so all 
%     % tether-sticker pairs have tether in the first columne
%     meso_pair_node_types = zeros(size(meso_pairs));
%     for i=1:size(meso_pairs,1)
%         for j=1:size(meso_pairs,2)
% %             meso_pair_node_types(:,j) = ...
% %                 type_temp(ismember(meso_pairs(:,j),id_temp));
%             id_ij = id_temp(meso_pairs(i,j));
%             meso_pair_node_types(i,j) = type_temp(id_temp==id_ij);
%         end
%     end
%     [~,indx] = sortrows(meso_pair_node_types);
% end

figure(100); clf; hold on

% Plot edges of bounding box
C1 = [-Length/2 -Width/2 -Height/2]/length_conversion;
C2 = [Length/2 -Width/2 -Height/2]/length_conversion;
C3 = [Length/2 Width/2 -Height/2]/length_conversion;
C4 = [-Length/2 Width/2 -Height/2]/length_conversion;
C5 = [-Length/2 -Width/2 Height/2]/length_conversion;
C6 = [Length/2 -Width/2 Height/2]/length_conversion;
C7 = [Length/2 Width/2 Height/2]/length_conversion;
C8 = [-Length/2 Width/2 Height/2]/length_conversion;
brdr_width = 0.5;
brdr_clr = [0.5 0.5 0.5];
plot3([C1(1) C2(1)],[C1(2) C2(2)],[C1(3) C2(3)],'Color',brdr_clr,'LineStyle','-','LineWidth',brdr_width)
plot3([C2(1) C3(1)],[C2(2) C3(2)],[C2(3) C3(3)],'Color',brdr_clr,'LineStyle','-','LineWidth',brdr_width)
plot3([C3(1) C4(1)],[C3(2) C4(2)],[C3(3) C4(3)],'Color',brdr_clr,'LineStyle','-','LineWidth',brdr_width)
plot3([C4(1) C1(1)],[C4(2) C1(2)],[C4(3) C1(3)],'Color',brdr_clr,'LineStyle','-','LineWidth',brdr_width)
plot3([C1(1) C5(1)],[C1(2) C5(2)],[C1(3) C5(3)],'Color',brdr_clr,'LineStyle','-','LineWidth',brdr_width)
plot3([C2(1) C6(1)],[C2(2) C6(2)],[C2(3) C6(3)],'Color',brdr_clr,'LineStyle','-','LineWidth',brdr_width)
plot3([C3(1) C7(1)],[C3(2) C7(2)],[C3(3) C7(3)],'Color',brdr_clr,'LineStyle','-','LineWidth',brdr_width)
plot3([C4(1) C8(1)],[C4(2) C8(2)],[C4(3) C8(3)],'Color',brdr_clr,'LineStyle','-','LineWidth',brdr_width)
plot3([C5(1) C6(1)],[C5(2) C6(2)],[C5(3) C6(3)],'Color',brdr_clr,'LineStyle','-','LineWidth',brdr_width)
plot3([C6(1) C7(1)],[C6(2) C7(2)],[C6(3) C7(3)],'Color',brdr_clr,'LineStyle','-','LineWidth',brdr_width)
plot3([C7(1) C8(1)],[C7(2) C8(2)],[C7(3) C8(3)],'Color',brdr_clr,'LineStyle','-','LineWidth',brdr_width)
plot3([C8(1) C5(1)],[C8(2) C5(2)],[C8(3) C5(3)],'Color',brdr_clr,'LineStyle','-','LineWidth',brdr_width)

% Plot the springs
line_width = 4;

% Define the starts and ends of each spring
temp_indx = zeros(size(meso_pairs,1),1);
temp_indx2 = zeros(size(meso_pairs,1),1);
for prs=1:size(meso_pairs,1)
    temp_indx(prs) = find(id_temp==meso_pairs(prs,1));
    temp_indx2(prs) = find(id_temp==meso_pairs(prs,2));
end
if BeadSpringOrMeso==0
    mult = 1;
elseif BeadSpringOrMeso==1
    mult = -1;
end
x_froms = x_temp(temp_indx);
x_tos = x_froms + mult*rx_temp';
y_froms = y_temp(temp_indx);
y_tos = y_froms + mult*ry_temp';
z_froms = z_temp(temp_indx);
z_tos = z_froms + mult*rz_temp';

% [p2_replace,rplc_indx2,p_int2] = ...
p_int2 = PlotTheSprings(x_froms,x_tos,y_froms,y_tos,z_froms,z_tos,...
    rx_temp,ry_temp,rz_temp,...
    Length,Width,Height,ColorCoding,line_width);

% Show the crosslink positions
temp_indx = zeros(size(meso_pairs,1),1);
for prs=1:size(meso_pairs,1)
    temp_indx(prs) = find(id_temp==meso_pairs(prs,2));
end
if BeadSpringOrMeso==0
    mult = -1;
elseif BeadSpringOrMeso==1
    mult = 1;
end
x_froms = x_temp(temp_indx);
x_tos = x_froms + mult*rx_temp';
y_froms = y_temp(temp_indx);
y_tos = y_froms + mult*ry_temp';
z_froms = z_temp(temp_indx);
z_tos = z_froms + mult*rz_temp';
% x_froms = x_temp(temp_indx);
% x_tos = x_froms - rx_temp';
% y_froms = y_temp(temp_indx);
% y_tos = y_froms - ry_temp';
% z_froms = z_temp(temp_indx);
% z_tos = z_froms - rz_temp';

% [p1_replace,rplc_indx1,p_int1] = ...
p_int1 = PlotTheSprings(x_froms,x_tos,y_froms,y_tos,z_froms,z_tos,...
    rx_temp,ry_temp,rz_temp,...
    Length,Width,Height,ColorCoding,line_width);

% Show where periodic nodes cross boundary
% Plot the nodes as points
ShowStickers = 1;
th_type = 1;
st_type = 2;
tether_clr = 'c';
sticker_clr = 'c';
node_size = 20;
bnd_size = 15;
if ShowStickers==1
    all = unique(meso_pairs(:));
    x_all = x_temp(ismember(id_temp,all));
    y_all = y_temp(ismember(id_temp,all));
    z_all = z_temp(ismember(id_temp,all));
    type_all = type_temp(ismember(id_temp,all));

    % plot the tethers
    X1 = [x_all(type_all==th_type),y_all(type_all==th_type),z_all(type_all==th_type)];
    s = scatter3(X1(:,1),X1(:,2),X1(:,3),tether_clr,'filled');
    s.SizeData = node_size;
    s.MarkerEdgeColor = 'k';

    % plot the stickers
    X1 = [x_all(type_all==st_type),y_all(type_all==st_type),z_all(type_all==st_type)];
    s = scatter3(X1(:,1),X1(:,2),X1(:,3),sticker_clr,'filled');
    s.SizeData = node_size;
    s.MarkerEdgeColor = 'k';
    s.LineWidth = 0.25;

    % plot where bonds cross the boundaries
    s = scatter3(p_int1(:,1),p_int1(:,2),p_int1(:,3),'filled');
    s.MarkerFaceColor = [0.85 0.85 0.85];
    s.SizeData = bnd_size;
    s.LineWidth = 0.1;
    s.MarkerEdgeColor = 'none';
    s.MarkerEdgeColor = [0.5 0.5 0.5];

    s = scatter3(p_int2(:,1),p_int2(:,2),p_int2(:,3),'filled');
    s.MarkerFaceColor = [0.85 0.85 0.85];
    s.SizeData = bnd_size;
    s.LineWidth = 0.1;
    s.MarkerEdgeColor = [0.5 0.5 0.5];
end


view(45,30)
daspect([1 1 1])
xlim([-3/2*Length0 3/2*Length0]/length_conversion)
ylim([-1/2*Length0 1/2*Length0]/length_conversion)
zlim([-1.01*Length0 1/2*Length0]/length_conversion)

set(gcf,'Position',[100 100 1200 920])
set(gcf,'color','k')
set(gca,'Color','k')
box off
axis off
xticks([])
yticks([])
zticks([])

% Compute order parameter and add to title
r_mag = vecnorm([rx_temp' ry_temp' rz_temp'],2,2);
rx_hat = rx_temp'./r_mag;
ry_hat = ry_temp'./r_mag;
rz_hat = rz_temp'./r_mag;
order_param = vecnorm(sum([rx_hat ry_hat rz_hat],1),2,2)/size(rx_hat,1);

tit = [tit,', $\varphi$ = ',num2str(order_param,'%.2f')];

title(tit,'FontSize',16,'Interpreter','latex','Color','w')

cbar = colorbar;
cbar.Color = 'w';
cbar.Label.String = color_lab;
cbar.Label.Interpreter = 'latex';
cbar.Label.Color = 'w';
cbar.Label.FontSize = 20;
cbar.Ticks = cbar_ticks;
cbar.TickLabels = cbar_tick_labels;
cbar.TickLabelInterpreter = 'latex';
set(gca,'FontSize',20);

% Project quivers onto xy-plane
% r_mag = vecnorm([rx_temp',ry_temp',rz_temp'],2,2);
r_mag2d = vecnorm([rx_temp',ry_temp'],2,2);
qu = rx_temp'./r_mag2d;
qv = ry_temp'./r_mag2d;
qw = zeros(size(r_mag2d));
qx = x_froms;
qy = y_froms;
qz = -Length0/length_conversion*ones(size(qy));

% color_frac = vecnorm([qu qv qw],2,2)
% c_rng = (linspace(0,1,100))';
% c = [c_rng flipud(c_rng) zeros(size(c_rng))];
% c = hot;

% c_rng = (linspace(0,1,100))';
% c_rng = (logspace(-2,0,100))';
% c_rng(1:10) = [];
% c = [c_rng c_rng c_rng];

% c_rng = (linspace(0,1,100))';
c_rng = (logspace(-2,0,10))';
% c = [ones(size(c_rng)) c_rng c_rng];
c_rng(1:2) = [];
c = [c_rng zeros(size(c_rng)) zeros(size(c_rng))];

colors = DefineColors(c,abs(rx_temp'./r_mag));

for i=1:size(colors,1)
    q = quiver3(qx(i),qy(i),qz(i),qu(i),qv(i),qw(i),'Color',colors(i,:));
    q.LineWidth = 1.5;
    q.AutoScaleFactor = 0.2;
    q.MaxHeadSize = 1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p_int = ...
    PlotTheSprings(x_froms,x_tos,y_froms,y_tos,z_froms,z_tos,...
    rx_temp,ry_temp,rz_temp,...
    Length,Width,Height,ColorCoding,line_width)

global length_conversion N_Kuhn b

% Crop any springs going outside of the box
p1 = [x_froms y_froms z_froms]; % starting points for each spring
p2 = [x_tos y_tos z_tos];       % ending points for each spring

% wherever p2 is out of bounds and p1 is in, redefine p2
[p2_replace,p_int] = FindOutOfBoundPoints(p1,p2,...
    Length/length_conversion,Width/length_conversion,Height/length_conversion);

% Define the color scale
r_mag = vecnorm([rx_temp',ry_temp',rz_temp'],2,2);
r_max = 0.75*N_Kuhn*b;
if ColorCoding==1   % absolute chain length
    color_frac = r_mag/r_max;
elseif ColorCoding==2   % length in x-direction 
    color_frac = abs(rx_temp'/r_max);
elseif ColorCoding==3   % alignment with x-direction abs(x \cdot r)
    color_frac = abs(rx_temp'./r_mag);
end
color_frac(color_frac>=1) = 1;
c = parula;
colors = DefineColors(c,color_frac);

% Plot the bonds
for i=1:size(colors,1)
    color = colors(i,:);
    plot3([p1(i,1) p2_replace(i,1)],...
        [p1(i,2) p2_replace(i,2)],...
        [p1(i,3) p2_replace(i,3)],...
        'LineWidth',line_width,'Color',color);
%     plot3([p1(i,1) p2(i,1)],...
%         [p1(i,2) p2(i,2)],...
%         [p1(i,3) p2(i,3)],...
%         'LineWidth',line_width,'Color',color);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p2_replace,p_int] = ...
    FindOutOfBoundPoints(p1,p2,Length,Width,Height)

% Find replacemnt points for p2s outside of bounds
p_in = p1(~IsInside(p2,Length,Width,Height),:); % p1s whose p2s are out of bounds
p_out = p2(~IsInside(p2,Length,Width,Height),:);% p2s that are out of bounds
% replacement_indices_p2 = find(~IsInside(p2,Length,Width,Height)); % rows of p2s that are out of bounds

% Find points where the lines through p_in and p_out intersect with bounds
p_int = zeros(size(p_out,1),3); % p2_new will be the replacement values for p2
for i=1:size(p_out,1)
    p_in_temp = p_in(i,:);
    p_out_temp = p_out(i,:);
    p_int(i,:) = ...
        FindBoundaryIntersection(p_in_temp,p_out_temp,Length,Width,Height);
end

p2_replace = p2;
p2_replace(~IsInside(p2,Length,Width,Height),:) = p_int;

% p2_replace = p2;
% p2_replace(replacement_indices_p2,:) = p_int;
% 
% % Run this algorithm one more time
% p_in = p1(~IsInside(p2,Length,Width,Height),:); % p1s whose p2s are out of bounds
% p_out = p2_replace(~IsInside(p2_replace,Length,Width,Height),:);% p2s that are out of bounds
% replacement_indices_p2 = find(~IsInside(p2_replace,Length,Width,Height)); % rows of p2s that are out of bounds
% 
% % Find points where the lines through p_in and p_out intersect with bounds
% p_int2 = zeros(size(p_out,1),3); % p2_new will be the replacement values for p2
% for i=1:size(p_out,1)
%     p_in_temp = p_in(i,:);
%     p_out_temp = p_out(i,:);
%     p_int2(i,:) = ...
%         FindBoundaryIntersection(p_in_temp,p_out_temp,Length,Width,Height);
% end
% 
% p2_replace = p2;
% p2_replace(replacement_indices_p2,:) = p_int2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pt_int = ...
    FindBoundaryIntersection(p_in,p_out,Length,Width,Height)


% Extract line parameters
x_in = p_in(1); y_in = p_in(2); z_in = p_in(3);
pt_int = p_out;

ct = 0;
while abs(pt_int(1))>Length/2 || abs(pt_int(2))>Width/2 || abs(pt_int(3))>Height/2
    ct = ct+1;
    x_out = pt_int(1); y_out = pt_int(2); z_out = pt_int(3);

    % Define vector of line
    u_line = (p_out-p_in)/vecnorm(p_out-p_in,2,2);  % director of line
    pt_line = p_out;                                % pt defining line

    % Define the plane with point and norm
    [u_plane,pt_plane] = DefineThePlane(x_out,y_out,z_out,Length,Width,Height);

    % Based on point and plane, find the intersection
    [pt_int,~] = LinePlaneIntersection(u_line,pt_line,u_plane,pt_plane);

    if ct>3 % Should not have to repeat more than number of dimensions
        break;
    end
end

check_figs=0;
if check_figs ==1
    view(45,30); daspect([1 1 1]); xlabel('x'); ylabel('y'); zlabel('z');
    scatter3(x_in,y_in,z_in,'g','filled')
    scatter3(x_out,y_out,z_out,'r','filled')
    quiver3(pt_line(1),pt_line(2),pt_line(3),u_line(1),u_line(2),u_line(3),'b')
    quiver3(pt_plane(1),pt_plane(2),pt_plane(3),u_plane(1),u_plane(2),u_plane(3),'b')
    scatter3(pt_int(1),pt_int(2),pt_int(3),...
        'k','filled');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u_plane,pt_plane] = DefineThePlane(x_out,y_out,z_out,...
    Length,Width,Height)

% Define which plane the bond passes through
if x_out<-Length/2 || x_out>Length/2
    max_dim = 1;
elseif y_out<-Width/2 || y_out>Width/2
    max_dim = 2;
elseif z_out<-Height/2 || z_out>Height/2
    max_dim = 3;
end
% max_dim = find(abs(p_out)==max(abs(p_out)));
switch max_dim
    case 1
        if x_out>0 % if particle is beyond x-bound
            pt_plane = [Length Width Height]/2;
        elseif x_out<0 % if particle is below x-bound
            pt_plane = -[Length Width Height]/2;
        end
        u_plane = [1 0 0];
    case 2
        if y_out>0 % if particle is beyond y-bound
            pt_plane = [Length Width Height]/2;
        elseif y_out<0 % if particle is below y-bound
            pt_plane = -[Length Width Height]/2;
        end
        u_plane = [0 1 0];
    case 3
        if z_out>0 % if particle is beyond z-bound
            pt_plane = [Length Width Height]/2;
        elseif z_out<0 % if particle is below z-bound
            pt_plane = -[Length Width Height]/2;
        end
        u_plane = [0 0 1];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I,rc] = LinePlaneIntersection(u, N, n, M, verbose)
%% line_plane_intersection : function to compute the intersection point
% between the (N,u) line and the (M,n) plane of the 3D space.
%
% Author & support : nicolas.douillet (at) free.fr, 2019-2020.
%
% Syntax
%
% [I,rc] = line_plane_intersection(u, N, n, M);
% [I,rc] = line_plane_intersection(u, N, n, M, verbose);
%
%
% Description
%
% [I,rc] = line_plane_intersection(u, N, n, M) computes the coordinates of I,
% the intersection point between the line (u,N) and the plane (n,M).
% In the most generic case, I is a point in the 3D space, but when
% the line is stricly parallel to the plane, I is the empty set, and when
% the line is included in the plane, I is a function handle corresponding
% to the system of parametric equations of the line.
%
% [I,rc] = line_plane_intersection(u, N, n, M, verbose) displays a message in
% console when verbose is set either to logical true or real numeric 1, and
% doesn't when it is set to logical false or real numeric 0.
%
%
% Principle
%
% Based on solving Descartes plane equation :
%
% ax + by + cz + d = 0, where n = [a, b, c] is a vector normal to the plane,
%
% combined with the parametric equations system of a 3D line :
%
% x(t) = x0 + at 
% y(t) = y0 + bt
% z(t) = z0 + ct
%
% where N0 = [x0, y0, z0] is a point belonging to the line, and u = [a, b, c], a vector directing this line.
%
%
% Input arguments
%
% - u : real row or column vector double. numel(u) = 3. One director vector of the parametric line.
%
% - N : real row or column vector double. numel(N) = 3. One point belonging to the line.
%
% - n : real row or column vector double. numel(n) = 3. One normal vector to the plane.
%
% - M : real row or column vector double. numel(M) = 3. One point belonging to the plane.
%
% - verbose : logical *true (1)/false(0), to enable/disable the verbose mode.
%
%
% Output arguments
%
% - I = [xI yI zI], real row or column vector double, the intersection point.
%
% - rc : return code, integer scalar doublein the set {1,2,3}.
%        0 : void / [] intersection
%        1 : point intersection (unique).
%        2 : line intersection
%
%        rc return code is necessary to distinguish between cases where
%        (N,u) line and the (M,n) plane intersection is a single point
%        and where it is the line itself.

%% Input parsing
assert(nargin > 3,'Not enough input arguments.');
assert(nargin < 6,'Too many input arguments.');
if nargin < 5    
    verbose = true;    
else    
    assert(islogical(verbose) || isreal(verbose),'verbose must be of type either logical or real numeric.');    
end
assert(isequal(size(u),size(N),size(n),size(M)),'Inputs u, M, n, and M must have the same size.');
assert(isequal(numel(u),numel(N),numel(n),numel(M),3),'Inputs u, M, n, and M must have the same number of elements (3).');
assert(isequal(ndims(u),ndims(N),ndims(n),ndims(M)),'Inputs u, M, n, and M must have the same number of dimensions.');
%% Body
% Plane offset parameter
d = -dot(n,M);
% Specific cases treatment
if ~dot(n,u) % n & u perpendicular vectors
    if dot(n,N) + d == 0 % N in P => line belongs to the plane
        if verbose
            disp('(N,u) line belongs to the (M,n) plane. Their intersection is the whole (N,u) line.');
        end
        I = M;
        rc = 2;
    else % line // to the plane
        if verbose
            disp('(N,u) line is parallel to the (M,n) plane. Their intersection is the empty set.');
        end
        I = [];
        rc = 0;
    end
else
    
    % Parametric line parameter t
    t = - (d + dot(n,N)) / dot(n,u);
    
    % Intersection coordinates
    I = N + u*t;
    
    rc = 1;
    
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inside = IsInside(point, Length, Width, Height)

x = point(:,1);
y = point(:,2);
z = point(:,3);

inside = (x >= -Length/2 & x <= Length/2) & ...
    (y >= -Width/2 & y <= Width/2) & ...
    (z >= -Height/2 & z <= Height/2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function colors = DefineColors(map,val)

x_interp = (linspace(0,1,length(map)))';
colors = [interp1(x_interp,map(:,1),val),...
    interp1(x_interp,map(:,2),val),...
    interp1(x_interp,map(:,3),val)];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kd,ka,ka1,dynamic_pairs0,N_bound0,N_free0,timestep0,...
    all_dynamic_pairs,N_virgin_bonds,no_attachments,no_detachments] = ...
    ComputeBondDynamics(timestep,i,dt,Dt,...
    p1,p2,bond_types,bond_id,...
    types,dynamic_pairs0,...
    N_bound0,N_free0,timestep0,tiptyp,dynbondtype,...
    all_dynamic_pairs,N_virgin_bonds)

prev_indx = i-1;

% Define timestep
% Dt = (timestep-timestep0)*dt;
Time = timestep*dt;

% Define current dynamic pairs
if ~isempty(p1{i})
    pairs = [p1{i} p2{i} bond_types{i} bond_id{i}];
    dynamic_pairs = pairs(pairs(:,3)==dynbondtype,:);
    pair_col = sort(dynamic_pairs(:,1:2),2);
    dynamic_pairs(:,1:2) = pair_col;
    dynamic_pairs = sortrows(dynamic_pairs,1);
    dynamic_pairs(1:2:end,:) = [];

    % compute number of bound and free stickers
    N_stickers = length(types(types==tiptyp));
    N_bound = size(dynamic_pairs,1);
    N_free = N_stickers-N_bound;

    % Define previous dynamic pairs
%     if i==1
%         dynamic_pairs0 = [];
%     else
    if ~isempty(dynamic_pairs0)
        pairs0 = [p1{prev_indx} p2{prev_indx} bond_types{prev_indx} bond_id{prev_indx}];
        dynamic_pairs0 = pairs0(pairs0(:,3)==dynbondtype,:);
        pair_col = sort(dynamic_pairs0(:,1:2),2);
        dynamic_pairs0(:,1:2) = pair_col;
        dynamic_pairs0 = sortrows(dynamic_pairs0,1);
        dynamic_pairs0(1:2:end,:) = [];
    end

    % ID previously existing and new dynamic bonds
    current_pairs = dynamic_pairs(:,1:2);
    if ~isempty(dynamic_pairs0)
        prev_pairs = dynamic_pairs0(:,1:2);
        new_pairs = current_pairs(~ismember(current_pairs,prev_pairs,'rows'),:);
    else
        prev_pairs = [];
        new_pairs = current_pairs;
    end

    % Define number of attachments this step & corresponding attachment rate
    if i==1  || isempty(new_pairs)
        ka = 0;
        no_attachments = 0;
    else
        no_attachments = size(new_pairs,1);   
        ka = no_attachments/N_free0/Dt;
    end

    % Find first-time new pairs and define ka1
    if ~isempty(new_pairs) && isempty(all_dynamic_pairs)
        N_first_time_attachments = size(new_pairs,1);
        ka1 = N_first_time_attachments/N_virgin_bonds/Dt;
    elseif ~isempty(new_pairs) && ~isempty(all_dynamic_pairs)
        np = sort(new_pairs,2);
        ap = sort(all_dynamic_pairs(:,1:2),2);

        % delete pairs that already existed at some point
        np(ismember(np,ap,'rows'),:) = [];

        N_first_time_attachments = size(np,1);
        ka1 = N_first_time_attachments/N_virgin_bonds/Dt;
    else
        N_first_time_attachments = 0;
        ka1 = 0;
    end
    N_virgin_bonds = N_virgin_bonds-N_first_time_attachments;

    % Store all pairs history
    if ~isempty(new_pairs) % Add the attachment events
        all_dynamic_pairs = [all_dynamic_pairs;[new_pairs,...
            Time*ones(size(new_pairs,1),1),...
            ones(size(new_pairs,1),1),...
            i*ones(size(new_pairs,1),1)]];  % 1 in fourth position indicates onset of a new pair
    end

    % ID lost bonds (i.e., detachments)
    if isempty(dynamic_pairs0)
        deleted_pairs = [];
    else
        surviving_nodes1 = ismember(prev_pairs(:,1),current_pairs(:,1));
        surviving_nodes2 = ismember(prev_pairs(:,2),current_pairs(:,2));
        indx = surviving_nodes1 + surviving_nodes2;
        deleted_pairs = prev_pairs(indx~=2,:);
    end

    no_detachments =  size(deleted_pairs,1);

    if i==1 || N_bound0==0
        kd = 0;
    else
        kd = no_detachments/N_bound0/Dt;
    end

    % disp(['ka = ',num2str(ka),', kd = ',num2str(kd)])
    if ~isempty(deleted_pairs) % Add the deletion events
        all_dynamic_pairs = [all_dynamic_pairs;[deleted_pairs,...
            Time*ones(size(deleted_pairs,1),1),...
            zeros(size(deleted_pairs,1),1),...
            i*ones(size(deleted_pairs,1),1)]]; % 1 in fourth position indicates deletion of a pair
    end

    % Set old values
    dynamic_pairs0 = dynamic_pairs;
    N_bound0 = N_bound;
    N_free0 = N_free;
    timestep0 = timestep;
else
    ka1 = 0;
    ka = 0;
    kd = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [msd,msd_err,msd0_r,msd0_r_err,msd0_t,msd0_t_err,...
    msd00_r,msd00_r_err,msd00_t,msd00_t_err,...
    msd_tether,msd_err_tether]...
    = ComputeMSDs(x,y,z,x0,y0,z0,...
    x00,y00,z00,Corners,...
    types,types0,types00,id,id0,id00,mol,mol0,mol00,...
    stickertype,tethertype)

Lx = abs(Corners(5)-Corners(2));
Ly = abs(Corners(6)-Corners(3));
Lz = abs(Corners(7)-Corners(4));

% if BeadSpringOrMeso==1
    % Check that these are aligned
    AtomData00 = [x00 y00 z00 types00 id00 mol00];
    AtomData0 = [x0 y0 z0 types0 id0 mol0];
    AtomDataC = [x y z types id mol];

    % Sort AtomDatas by particle ID
    AtomData00 = sortrows(AtomData00,5);
    AtomData0 = sortrows(AtomData0,5);
    AtomDataC = sortrows(AtomDataC,5);

    Pos00 = AtomData00(:,1:3);
    Pos = AtomDataC(:,1:3);

    % Find tether particle of corresponding molecule
    AtomData_st = AtomDataC(AtomDataC(:,4)==tethertype,:);
    AtomData_th = AtomDataC(AtomDataC(:,4)==stickertype,:);
    AtomData_st = sortrows(AtomData_st,6);
    AtomData_th = sortrows(AtomData_th,6);
    R = AtomData_st(:,1:3)-AtomData_th(:,1:3);

    % if sum(ismember(AtomData_st(:,end),AtomData_th(:,end),'rows'))~=size(AtomData_st,1)
    %     [AtomData_st(:,end),AtomData_th(:,end)];
    % end

    AtomData0_st = AtomData0(AtomData0(:,4)==tethertype,:);
    AtomData0_th = AtomData0(AtomData0(:,4)==stickertype,:);
    AtomData0_st = sortrows(AtomData0_st,6);
    AtomData0_th = sortrows(AtomData0_th,6);
    R0 = AtomData0_st(:,1:3)-AtomData0_th(:,1:3);

    AtomData00_st = AtomData00(AtomData00(:,4)==tethertype,:);
    AtomData00_th = AtomData00(AtomData00(:,4)==stickertype,:);
    AtomData00_st = sortrows(AtomData00_st,6);
    AtomData00_th = sortrows(AtomData00_th,6);
    R00 = AtomData00_st(:,1:3)-AtomData00_th(:,1:3);

%     if sum(ismember(AtomData00_st(:,end),AtomData00_th(:,end),'rows'))~=size(AtomData00_st,1)
%         [AtomData00_st(:,end),AtomData00_th(:,end)];
%     end

    %Correct for periodic bounds by moving R0 back, as needed
    dr0_temp = R-R0;
    R0(dr0_temp(:,1)>Lx/2,1) = R0(dr0_temp(:,1)>Lx/2,1)+Lx;
    R0(dr0_temp(:,1)<-Lx/2,1) = R0(dr0_temp(:,1)<-Lx/2,1)-Lx;
    R0(dr0_temp(:,2)>Lx/2,2) = R0(dr0_temp(:,2)>Lx/2,2)+Lx;
    R0(dr0_temp(:,2)<-Lx/2,2) = R0(dr0_temp(:,2)<-Lx/2,2)-Lx;
    R0(dr0_temp(:,3)>Lx/2,3) = R0(dr0_temp(:,3)>Lx/2,3)+Lx;
    R0(dr0_temp(:,3)<-Lx/2,3) = R0(dr0_temp(:,3)<-Lx/2,3)-Lx;

    theta0 = acos(dot(R0,R,2)./(vecnorm(R0,2,2).*vecnorm(R,2,2)));
    dR0_r = vecnorm(R,2,2)-vecnorm(R0,2,2);
    dR0_t = vecnorm(R0,2,2).*theta0;
    dR0_r(isnan(dR0_r)) = [];
    dR0_t(isnan(dR0_t)) = [];
    msd0_r = mean(dR0_r.^2);
    msd0_r_err = std(dR0_r.^2)/sqrt(size(dR0_r,1));
    msd0_t = mean(dR0_t.^2);
    msd0_t_err = std(dR0_t.^2)/sqrt(size(dR0_t,1));

    %Correct for periodic bounds by moving R00 back, as needed
    dr00_temp = R-R00;
    R00(dr00_temp(:,1)>Lx/2,1) = R00(dr00_temp(:,1)>Lx/2,1)+Lx;
    R00(dr00_temp(:,1)<-Lx/2,1) = R00(dr00_temp(:,1)<-Lx/2,1)-Lx;
    R00(dr00_temp(:,2)>Lx/2,2) = R00(dr00_temp(:,2)>Lx/2,2)+Lx;
    R00(dr00_temp(:,2)<-Lx/2,2) = R00(dr00_temp(:,2)<-Lx/2,2)-Lx;
    R00(dr00_temp(:,3)>Lx/2,3) = R00(dr00_temp(:,3)>Lx/2,3)+Lx;
    R00(dr00_temp(:,3)<-Lx/2,3) = R00(dr00_temp(:,3)<-Lx/2,3)-Lx;

    theta00 = acos(dot(R00,R,2)./(vecnorm(R00,2,2).*vecnorm(R,2,2)));
    dR00_r = vecnorm(R,2,2)-vecnorm(R00,2,2);
    dR00_t = vecnorm(R00,2,2).*theta00;
    dR00_r(isnan(dR00_r)) = [];
    dR00_t(isnan(dR00_t)) = [];
    msd00_r = mean(dR00_r.^2);
    msd00_r_err = std(dR00_r.^2)/sqrt(size(dR00_r,1));
    msd00_t = mean(dR00_t.^2);
    msd00_t_err = std(dR00_t.^2)/sqrt(size(dR00_t,1));

    % Eliminate partilces of wrong type
    Pos00(AtomData00(:,4)~=stickertype,:) = [];
    Pos(AtomDataC(:,4)~=stickertype,:) = [];

    % N = size(Pos,1);
    dr = Pos-Pos00;

    % Adjust for periodic bounds
    if ~isempty(dr(abs(dr(:,1))>0.5*Lx,1)) || ~isempty(dr(abs(dr(:,2))>0.5*Ly,2)) ||...
            ~isempty(dr(abs(dr(:,3))>0.5*Lz,3))
        dr(abs(dr(:,1))>0.5*Lx,1);
    end
    dr(abs(dr(:,1))>0.5*Lx,1) = dr(abs(dr(:,1))>0.5*Lx,1)-Lx*sign(dr(abs(dr(:,1))>0.5*Lx,1));
    dr(abs(dr(:,2))>0.5*Ly,2) = dr(abs(dr(:,2))>0.5*Ly,2)-Ly*sign(dr(abs(dr(:,2))>0.5*Ly,2));
    dr(abs(dr(:,3))>0.5*Lz,3) = dr(abs(dr(:,3))>0.5*Lz,3)-Lz*sign(dr(abs(dr(:,3))>0.5*Lz,3));

    MSD_temp = mean((vecnorm(dr,2,2)).^2);
    MSD_err = std((vecnorm(dr,2,2)).^2)/sqrt(size(dr,1));

    msd = MSD_temp;
    msd_err = MSD_err;

    % Calculate MSD for tethers
    Pos00 = AtomData00(:,1:3);
    Pos = AtomDataC(:,1:3);

    % Eliminate partilces of wrong type
    Pos00(AtomData00(:,4)~=tethertype,:) = [];
    Pos(AtomDataC(:,4)~=tethertype,:) = [];

    % N = size(Pos,1);
    dr = Pos-Pos00;

    % Adjust for periodic bounds
    if ~isempty(dr(abs(dr(:,1))>0.5*Lx,1)) || ~isempty(dr(abs(dr(:,2))>0.5*Ly,2)) ||...
            ~isempty(dr(abs(dr(:,3))>0.5*Lz,3))
        dr(abs(dr(:,1))>0.5*Lx,1);
    end
    dr(abs(dr(:,1))>0.5*Lx,1) = dr(abs(dr(:,1))>0.5*Lx,1)-Lx*sign(dr(abs(dr(:,1))>0.5*Lx,1));
    dr(abs(dr(:,2))>0.5*Ly,2) = dr(abs(dr(:,2))>0.5*Ly,2)-Ly*sign(dr(abs(dr(:,2))>0.5*Ly,2));
    dr(abs(dr(:,3))>0.5*Lz,3) = dr(abs(dr(:,3))>0.5*Lz,3)-Lz*sign(dr(abs(dr(:,3))>0.5*Lz,3));

    MSD_temp = mean((vecnorm(dr,2,2)).^2);
    MSD_err = std((vecnorm(dr,2,2)).^2)/sqrt(size(dr,1));

    msd_tether = MSD_temp;
    msd_err_tether = MSD_err;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigma,sigma_vir,conc] = ComputeStressAndBondConc(fc,rc,Corners,...
    rx,ry,rz,bond_types)

global length_conversion BeadSpringOrMeso kbT b...

% Compute volume
L = length_conversion*abs(Corners(3)-Corners(2));
W = length_conversion*abs(Corners(5)-Corners(4));
H = length_conversion*abs(Corners(7)-Corners(6));
V = L*W*H;

rc = rc*length_conversion;

% compute the network-scale stress on a per-chain basis
% Dyad f \otimes r
dyad = [fc(:,1).*rc(:,1),...
    fc(:,2).*rc(:,2),...
    fc(:,3).*rc(:,3),...
    fc(:,1).*rc(:,2),...
    fc(:,2).*rc(:,3),...
    fc(:,3).*rc(:,1)];

n_chains = size(fc,1);
conc = n_chains/V;

sigma = 1/(2*V)*sum(dyad,1);

% compute the virial stress on a per Kuhn segment basis
if BeadSpringOrMeso==0
    rx_bckbn = rx(bond_types~=2);
    ry_bckbn = ry(bond_types~=2);
    rz_bckbn = rz(bond_types~=2);
    rvir = [rx_bckbn ry_bckbn rz_bckbn]*length_conversion;
    norms = vecnorm(rvir,2,2);

    % Based on LAMMPS nonlinear potential
    stiffness = 800;
    E = stiffness*kbT;  % Converted back to SI units
    l0 = b*length_conversion;
    r0 = b*length_conversion;
    r_samp = linspace(0,2*b*length_conversion,50);
    psi = E*((r_samp-r0).^2)./(l0^2-(r_samp-r0).^2);
    f_samp = gradient(psi)./gradient(r_samp);
    force_mags = interp1(r_samp,f_samp,norms);
    r_hat = rvir./norms;
    fvir = force_mags.*r_hat;
    f_hat = fvir./vecnorm(fvir,2,2);

    show_debug_fig = 0;
    if show_debug_fig==1
        figure(1e3); clf; hold on
        rng = (1:5)';

        % plot force and r vectors greater than b (tensile => should cause positive stress)
        r_hat_temp = r_hat;
        r_hat_temp = r_hat_temp(norms>b*length_conversion,:);
        f_hat_temp = f_hat;
        f_hat_temp = f_hat_temp(norms>b*length_conversion,:);
        quiver3(zeros(size(rng)),...
            zeros(size(rng)),...
            zeros(size(rng)),...
            r_hat_temp(rng,1),r_hat_temp(rng,2),r_hat_temp(rng,3),'g');
        quiver3(zeros(size(rng)),...
            zeros(size(rng)),...
            zeros(size(rng)),...
            f_hat_temp(rng,1),f_hat_temp(rng,2),f_hat_temp(rng,3),'c');

        % plot force and r vectors shorter than b (compressive => should cause negative stress)
        r_hat_temp = r_hat;
        r_hat_temp = r_hat_temp(norms<b*length_conversion,:);
        f_hat_temp = f_hat;
        f_hat_temp = f_hat_temp(norms<b*length_conversion,:);
        quiver3(zeros(size(rng)),...
            zeros(size(rng)),...
            zeros(size(rng)),...
            r_hat_temp(rng,1),r_hat_temp(rng,2),r_hat_temp(rng,3),'r');
        quiver3(zeros(size(rng)),...
            zeros(size(rng)),...
            zeros(size(rng)),...
            f_hat_temp(rng,1),f_hat_temp(rng,2),f_hat_temp(rng,3),'m');

        view(135,30)
        box on
    end

    % Dyad f \otimes r
    dyad_vir = [fvir(:,1).*rvir(:,1),...
        fvir(:,2).*rvir(:,2),...
        fvir(:,3).*rvir(:,3),...
        fvir(:,1).*rvir(:,2),...
        fvir(:,2).*rvir(:,3),...
        fvir(:,3).*rvir(:,1)];

    sigma_vir = 1/(2*V)*sum(dyad_vir,1);
else
    sigma_vir = [];     % sigma_vir is empty for mesoscale
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function forces = ComputeForces(r,N_Kuhn,b,kbT)

global length_conversion BeadSpringOrMeso

% types = types';
b = b*length_conversion;
r = r*length_conversion;
norms = vecnorm(r,2,2);
r_hat = [r(:,1)./norms r(:,2)./norms r(:,3)./norms];

if BeadSpringOrMeso==1
    lam = norms/N_Kuhn/b;
    numer = lam.*(3 - lam.^2);
    denom = 1 - lam.^2;
    force_mags = kbT*numer./denom/b;

elseif BeadSpringOrMeso==0
    lambda = norms/(sqrt(N_Kuhn)*b);
    force_mags = kbT*lambda/(sqrt(N_Kuhn)*b).*...
        ((lambda.^2-3*N_Kuhn)./(lambda.^2-N_Kuhn));
end

forces = force_mags.*r_hat;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g11,g22,g33,g12,g23,g31,...
    g11_err,g22_err,g33_err,...
    g12_err,g23_err,g31_err] = ...
    ComputeMetricTensor(rx,ry,rz)

norms = (rx.^2 +ry.^2 +rz.^2).^(1/2);
rx_hat = rx./norms; ry_hat = ry./norms; rz_hat = rz./norms;

g11_temp = rx_hat.*rx_hat;
g22_temp = ry_hat.*ry_hat;
g33_temp = rz_hat.*rz_hat;
g12_temp = rx_hat.*ry_hat;
g23_temp = ry_hat.*rz_hat;
g31_temp = rz_hat.*rx_hat;

g11 = mean(g11_temp);
g22 = mean(g22_temp);
g33 = mean(g33_temp);
g12 = mean(g12_temp);
g23 = mean(g23_temp);
g31 = mean(g31_temp);

g11_err = std(g11_temp)/sqrt(length(g11_temp));
g22_err = std(g22_temp)/sqrt(length(g22_temp));
g33_err = std(g33_temp)/sqrt(length(g33_temp));
g12_err = std(g12_temp)/sqrt(length(g12_temp));
g23_err = std(g23_temp)/sqrt(length(g23_temp));
g31_err = std(g31_temp)/sqrt(length(g31_temp));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n_bonds,n_to_self,n_to_other,c_alpha]...
    = QuantifyClustering(id,mol,bond_types,batom1,batom2)

dyn_type = 2;
atoms1 = batom1(bond_types==dyn_type,:);
atoms2 = batom2(bond_types==dyn_type,:);

n_bonds = size(atoms1,1);

pair_data = [];
for i=1:n_bonds
    atom1 = atoms1(i);
    atom2 = atoms2(i);
    
    mol1 = mol(id==atom1);
    mol2 = mol(id==atom2);

    pair_data_temp = [atom1 atom2 mol1 mol2];
    pair_data = cat(1,pair_data,pair_data_temp);
end

% for each molecule:
n_mol = length(unique(mol));
n_bonds = zeros(n_mol,1);   % number of sticker bonds for each molecule
n_to_self = zeros(n_mol,1); % number of sticker bonds to self for each 
                            % molecule
n_to_other = zeros(n_mol,1);% number of sticker bonds to other molecules
                            % for each molecule
n_other = zeros(n_mol,1);   % number of unique, non-self neighbors
c_alpha = zeros(n_mol,1);   % clustering coefficient for molecule

if ~isempty(pair_data)
    % Sort pair_data and eliminate redundants
    pair_data(:,1:2) = sort(pair_data(:,1:2),2);
    pair_data(:,3:4) = sort(pair_data(:,3:4),2);
    pair_data = unique(pair_data,'rows');

    for mol_temp=1:n_mol
        % Find all bonded pairs with mol_temp as a member
        temp_pairs1 = pair_data(pair_data(:,3)==mol_temp,:);
        temp_pairs2 = pair_data(pair_data(:,4)==mol_temp,:);

        temp_pairs = [temp_pairs1;temp_pairs2];

        % Calculate total number bonded neighbors
        n_bonds(mol_temp) = size(temp_pairs,1);

        % Calculate number of bonds to self
        i1 = find(temp_pairs(:,3)==mol_temp);
        i2 = find(temp_pairs(:,4)==mol_temp);
        self_bonded_pairs = temp_pairs(intersect(i1,i2),:);
        n_to_self(mol_temp) = size(self_bonded_pairs,1);

        % Calculate number of bonds to others
        nonself_bonded_pairs = ...
            temp_pairs(sum(temp_pairs(:,3:4),2)~=2*mol_temp,:);
        n_to_other(mol_temp) = size(nonself_bonded_pairs,1);

        % Calculate number of unique, non-self bonded neighbors
        neighbs = nonself_bonded_pairs(:,3:4);
        neighbs(neighbs==mol_temp) = [];
        others = unique(neighbs);
        n_other(mol_temp) = length(others);

        % Calculate clustering coefficient at scale of molecules
        k_alpha = n_other(mol_temp);    % number of unique bonds to other 
                                        % molecukes
        p1 = ismember(pair_data(:,3),neighbs);
        p2 = ismember(pair_data(:,4),neighbs);
        temp_indx = p1+p2;
        temp_indx(temp_indx~=2) = 0;
        temp_indx(temp_indx==2) = 1;
        attached_neighbor_pairs = pair_data(temp_indx==1,3:4);
        attached_neighbor_pairs = unique(attached_neighbor_pairs,'rows');

%         if ~isempty(attached_neighbor_pairs)
%             attached_neighbor_pairs;
%         end

        T_alpha = size(attached_neighbor_pairs,1);  % numberr of mutual
                                        % attachments between attached
                                        % neighbors of mol_temp
        c_alpha(mol_temp) = 2*T_alpha/(k_alpha*(k_alpha-1));
    end
end

n_bonds = n_bonds';
n_to_self = n_to_self';
n_to_other = n_to_other';
c_alpha = c_alpha';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rs = ExtractEndToEnd(types,id,mol,x,y,z,bond_types,batom1,batom2,...
    rx,ry,rz,Np,i,samp,Corners)

global BeadSpringOrMeso bondtypes Nt Ns meso_pairs N_Kuhn b...
%      Lx Ly Lz meso_pairs0
 
CheckFig = 0;

L = abs(Corners(3)-Corners(2));
W = abs(Corners(5)-Corners(4));
H = abs(Corners(7)-Corners(6));

if BeadSpringOrMeso==1
    rx = rx(bond_types~=2);
    ry = ry(bond_types~=2);
    rz = rz(bond_types~=2);
    rs = [rx ry rz];
    
    meso_pairs = [batom1(bond_types==1), batom2(bond_types==1)];
else
    if isempty(meso_pairs) 
        % Initialize all vectors
        n_tethers = Nt-1;
        n_stickers = Nt*Ns; % per molecule
        Nsegments = Np*(n_tethers+n_stickers);
        rx_out = zeros(Nsegments,1);
        ry_out = zeros(Nsegments,1);
        rz_out = zeros(Nsegments,1);

        batom1_out = zeros(Nsegments,1);
        batom2_out = zeros(Nsegments,1);

        % Define segments of side-branches
        ct = 0;
        all_pairs = [batom1 batom2];
        all_pairs(bond_types==2,:) = [];
        mols = unique(mol);
        for i=1:length(mols)
            molecule = mols(i);     % define the molecule
            temp_types = types(mol==molecule);      % define the particle types belonging to the molecule
            temp_id = id(mol==molecule);            % find particles belonging to the molecule
            [members,sort_indx] = sort(temp_id);    % sort particles belonging to the molecule
            temp_types = temp_types(sort_indx);     % sort the particle types

            % Isolate pairs belonging to the molecule
            indx1 = find(ismember(all_pairs(:,1),members));     % find all bound pairs belonging to the molecule
            indx2 = find(ismember(all_pairs(:,2),members));
            indx_mps = intersect(indx1,indx2);
            mol_pairs = all_pairs(indx_mps,:);      % bonds belonging to the molecule (i.e., Kuhn segments)
            mol_rx = rx(indx_mps,:);
            mol_ry = ry(indx_mps,:);
            mol_rz = rz(indx_mps,:);

            tethers = members(1:Nt);
            sticker_rng = Nt+(1:n_stickers);
            stickers = members(sticker_rng);
            others = members;
            others(ismember(others,tethers)) = [];
            others(ismember(others,stickers)) = [];
            if isempty(temp_types(sticker_rng)==1)
                warning('Incorrect stickers selected')
            end

            % Debugging fig
            if CheckFig==1
                figure(100); clf; hold on
                tether_pos = [x(ismember(id,tethers)),...
                    y(ismember(id,tethers)),...
                    z(ismember(id,tethers))];
                sticker_pos = [x(ismember(id,stickers)),...
                    y(ismember(id,stickers)),...
                    z(ismember(id,stickers))];
                other_pos = [x(ismember(id,others)),...
                    y(ismember(id,others)),...
                    z(ismember(id,others))];
                s = scatter3(tether_pos(:,1),tether_pos(:,2),tether_pos(:,3),'filled');
                s.MarkerFaceColor = 'c'; s.SizeData = 10;
                s = scatter3(sticker_pos(:,1),sticker_pos(:,2),sticker_pos(:,3),'filled');
                s.MarkerFaceColor = 'r'; s.SizeData = 10;
                s = scatter3(other_pos(:,1),other_pos(:,2),other_pos(:,3),'filled');
                s.MarkerFaceColor = [0.5 0.5 0.5]; s.SizeData = 10;
                view([30,45])
            end
            % Define summated vector of bead-springs between adjacent tethers
            mol_pairs = sort(mol_pairs,2);  % sort the Kuhn segments' end particles from left-to-right
            [mol_pairs,unq_indx] = unique(mol_pairs,'rows');  % This should be unnecessary
            mol_rx = mol_rx(unq_indx);
            mol_ry = mol_ry(unq_indx);
            mol_rz = mol_rz(unq_indx);
            flag = 1;
            while ~isempty(mol_pairs)  %while there remains molecule pairs
                % ID first tether in line
                if flag  %if its the first particle of the chain segment
                    indx = find(ismember(mol_pairs(:,1),tethers),1,'first'); % find the first instance of a tether in the Kuhn segments list
                    from = mol_pairs(indx,1);   % Going from...
                    to = mol_pairs(indx,2);     % ... to
                    from0 = from;

                    [rx_temp,ry_temp,rz_temp] =  ...
                        Calcdr(from,to,x,y,z,id,mol_rx,mol_ry,mol_rz,indx,L,W,H);

                    flag = 0; add_vecs = 0;
                    rxs = []; rys = []; rzs = [];   % Reset the rxs, rys and rzs for this chain segment
                    ct = ct+1;
                else
                    indx3 = find(ismember(mol_pairs(:,1),from),1);
                    indx4 = find(ismember(mol_pairs(:,2),from),1);
                    indx = [indx3;indx4];
                    to  = mol_pairs(indx,mol_pairs(indx,:)~=from);

                    [rx_temp,ry_temp,rz_temp] =  ...
                        Calcdr(from,to,x,y,z,id,mol_rx,mol_ry,mol_rz,indx,L,W,H);

                    if ismember(to,tethers) || ismember(to,stickers)
                        flag = 1; add_vecs = 1;
                    end
                end

                % Append vector components to stored
                rxs = cat(1,rxs,rx_temp);
                rys = cat(1,rys,ry_temp);
                rzs = cat(1,rzs,rz_temp);

                if add_vecs
                    rx_out(ct) = sum(rxs);
                    ry_out(ct) = sum(rys);
                    rz_out(ct) = sum(rzs);

                    rx_check = x(id==to)-x(id==from0);
                    ry_check = y(id==to)-y(id==from0);
                    rz_check = z(id==to)-z(id==from0);

                    bound_fact = 2; % Bigger bound factor means more 
                                    % sensitively going to treat springs as
                                    % if they span the periodic boundaries

                    if rx_check<-L/bound_fact
                        rx_check = rx_check+L;
                    elseif rx_check>L/bound_fact
                        rx_check = rx_check-L;
                    end
%                     rx_check(rx_check<-L/bound_fact) = rx_check(rx_check<-L/bound_fact)+L;
%                     rx_check(rx_check>L/bound_fact) = rx_check(rx_check>L/bound_fact)-L;

                    if ry_check<-W/bound_fact
                        ry_check = ry_check+W;
                    elseif ry_check>W/bound_fact
                        ry_check = ry_check-W;
                    end
%                     ry_check(ry_check<-W/bound_fact) = ry_check(ry_check<-W/bound_fact)+W;
%                     ry_check(ry_check>W/bound_fact) = ry_check(ry_check>W/bound_fact)-W;

                    if rz_check<-H/bound_fact
                        rz_check = rz_check+H;
                    elseif rz_check>H/bound_fact
                        rz_check = rz_check-H;
                    end
%                      rz_check(rz_check<-H/bound_fact) = rz_check(rz_check<-H/bound_fact)+H;
%                     rz_check(rz_check>H/bound_fact) = rz_check(rz_check>H/bound_fact)-H;

                    if round(rx_check,2)~=round(rx_out(ct),2) || ...
                            round(ry_check,2)~=round(ry_out(ct),2) || ...
                            round(rz_check,2)~=round(rz_out(ct),2)
                        warning('Mismatch in values for two methods - check round off')
                    end

                    batom1_out(ct) = from0;
                    batom2_out(ct) = to;
                end

                if CheckFig==1
                    x1 = x(id==from);
                    y1 = y(id==from);
                    z1 = z(id==from);
                    x2 = x1+rxs(end);
                    y2 = y1+rys(end);
                    z2 = z1+rzs(end);
                    plot3([x1 x2],[y1 y2],[z1 z2],'k-')
                    if add_vecs
                        x4 = x(id==to);
                        y4 = y(id==to);
                        z4 = z(id==to);
                        x3 = x4-rx_out(ct);
                        y3 = y4-ry_out(ct);
                        z3 = z4-rz_out(ct);
                        plot3([x4 x3],[y4 y3],[z4 z3],'r')
                    end
                end

                % Initialize reference bead for next iteration
                from = to;      % the new particle 'going from' is the old 'going to'
                mol_pairs(indx,:) = [];
                mol_rx(indx,:) = [];
                mol_ry(indx,:) = [];
                mol_rz(indx,:) = [];
            end
        end
        rs = [rx_out ry_out rz_out];
        meso_pairs = [batom1_out batom2_out];
    else
        [~,indx] = sort(id,1);
        x_sort = x(indx); 
        y_sort = y(indx); 
        z_sort = z(indx); 
        x_from = x_sort(meso_pairs(:,1));
        y_from = y_sort(meso_pairs(:,1));
        z_from = z_sort(meso_pairs(:,1));
        x_to = x_sort(meso_pairs(:,2));
        y_to = y_sort(meso_pairs(:,2));
        z_to = z_sort(meso_pairs(:,2));
        rx_out = x_to-x_from;
        ry_out = y_to-y_from;
        rz_out = z_to-z_from;

        % adjust for periodic boundaries
        indx_less = find(rx_out<-L/2);
        indx_more = find(rx_out>L/2);
        if ~isempty(indx_less)
            rx_out(indx_less) = rx_out(indx_less)+L;
        end
        if ~isempty(indx_more)
            rx_out(indx_more) = rx_out(indx_more)-L;
        end

        indx_less = find(ry_out<-W/2);
        indx_more = find(ry_out>W/2);
        if ~isempty(indx_less)
            ry_out(indx_less) = ry_out(indx_less)+W;
        end
        if ~isempty(indx_more)
            ry_out(indx_more) = ry_out(indx_more)-W;
        end

        indx_less = find(rz_out<-H/2);
        indx_more = find(rz_out>H/2);
        if ~isempty(indx_less)
            rz_out(indx_less) = rz_out(indx_less)+H;
        end
        if ~isempty(indx_more)
            rz_out(indx_more) = rz_out(indx_more)-H;
        end

        if ~isempty(rx_out(abs(rx_out)>N_Kuhn*b))
            warning('spring too long in x-direction');
        end
        if ~isempty(ry_out(abs(ry_out)>N_Kuhn*b))
            warning('spring too long in y-direction');
        end
        if ~isempty(rz_out(abs(rz_out)>N_Kuhn*b))
            warning('spring too long in z-direction');
        end
% 
%         rx_out(rx_out<-L/2) = rx_out(rx_out<-L/2)+L;
%         rx_out(rx_out>L/2) = rx_out(rx_out>L/2)-L;
% 
%         ry_out(ry_out<-W/2) = ry_out(ry_out<-W/2)+W;
%         ry_out(ry_out>W/2) = ry_out(ry_out>W/2)-W;    
% 
%         rz_out(rz_out<-H/2) = rz_out(rz_out<-H/2)+H;
%         rz_out(rz_out>H/2) = rz_out(rz_out>H/2)-H;
% 
%         rx_check = rx_out(abs(rx_out)>1.1*N_Kuhn*b);
%         ry_check = ry_out(abs(ry_out)>1.1*N_Kuhn*b);
%         
%         
%         rz_check = rz_out(abs(rz_out)>1.1*N_Kuhn*b);
%         if ~isempty(rx_check) || ~isempty(ry_check) || ~isempty(rz_check) 
%             warning('One of the springs is too long')
%         end
        rs = [rx_out ry_out rz_out];
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rx_temp,ry_temp,rz_temp] =  ...
    Calcdr(from,to,x,y,z,id,mol_rx,mol_ry,mol_rz,indx,L,W,H)

NewMethod = 1;
if NewMethod==1
    % New method
    p1 = [x(id==from),y(id==from),z(id==from)];
    p2 = [x(id==to),y(id==to),z(id==to)];
    r_temp = p2-p1;

    % Adjust for periodic bounds
    if r_temp(1)>L/2
        r_temp(1) = r_temp(1)-L;
    end
    if r_temp(1)<-L/2
        r_temp(1) = r_temp(1)+L;
    end
    if r_temp(2)>W/2
        r_temp(2) = r_temp(2)-W;
    end
    if r_temp(2)<-W/2
        r_temp(2) = r_temp(2)+W;
    end
    if r_temp(3)>H/2
        r_temp(3) = r_temp(3)-H;
    end
    if r_temp(3)<-H/2
        r_temp(3) = r_temp(3)+H;
    end
    rx_temp = r_temp(1);
    ry_temp = r_temp(2);
    rz_temp = r_temp(3);
else
    % Old method
    if from>to                  % if going from bigger to smaller, sign is positive, else negative
        sign = 1;
    else
        sign = -1;
    end
    % Alternate sign check - if there exists a particle in the
    % position of +/- or +/- accross any of the periodic
    % bounds, then it is the correct position/sign. Otherwise,
    % swap it and check again
    rx_temp = sign*mol_rx(indx);
    ry_temp = sign*mol_ry(indx);
    rz_temp = sign*mol_rz(indx);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PreallocateOutputs(n_steps,n_samples,Np)

global Nt Ns...
    sig11 sig22 sig33 sig12 sig23 sig31...
    sig11_err sig22_err sig33_err sig12_err sig23_err sig31_err...
    bondconc time stretch V...
    msd_st msd_st_err...
    msd_st_r0 msd_st_r0_err...
    msd_st_r00 msd_st_r00_err...
    msd_st_t0 msd_st_t0_err...
    msd_st_t00 msd_st_t00_err...
    msd_tether msd_err_tether...
    ka_out ka1_out kd_out Na Nd...
    rx ry rz bondtypes...
    n_bonds n_to_self n_to_other c_alpha...
    rxr11_st rxr22_st rxr33_st rxr12_st rxr23_st rxr31_st...
    rxr11_st_err rxr22_st_err rxr33_st_err rxr12_st_err rxr23_st_err rxr31_st_err...
    no_attachments no_detachments n_success_all

% Stress
sig11 = zeros(n_success_all,n_samples);
sig22 = zeros(n_success_all,n_samples);
sig33 = zeros(n_success_all,n_samples);
sig12 = zeros(n_success_all,n_samples);
sig23 = zeros(n_success_all,n_samples);
sig31 = zeros(n_success_all,n_samples);
sig11_err = zeros(n_success_all,n_samples);
sig22_err = zeros(n_success_all,n_samples);
sig33_err = zeros(n_success_all,n_samples);
sig12_err = zeros(n_success_all,n_samples);
sig23_err = zeros(n_success_all,n_samples);
sig31_err = zeros(n_success_all,n_samples);
time = zeros(n_success_all,1);
stretch = zeros(n_success_all,1);
V = zeros(n_success_all,1);
bondconc = zeros(n_success_all,n_samples);

% MSD
msd_st = zeros(n_success_all,n_samples);
msd_st_err = zeros(n_success_all,n_samples);
msd_st_r0 = zeros(n_success_all,n_samples);
msd_st_r0_err = zeros(n_success_all,n_samples);
msd_st_r00 = zeros(n_success_all,n_samples);
msd_st_r00_err = zeros(n_success_all,n_samples);
msd_st_t0 = zeros(n_success_all,n_samples);
msd_st_t0_err = zeros(n_success_all,n_samples);
msd_st_t00 = zeros(n_success_all,n_samples);
msd_st_t00_err = zeros(n_success_all,n_samples);
msd_tether = zeros(n_success_all,n_samples);
msd_err_tether = zeros(n_success_all,n_samples);

% Kinetics
ka_out = zeros(n_success_all,n_samples);
ka1_out = zeros(n_success_all,n_samples);
kd_out = zeros(n_success_all,n_samples);
Na = zeros(n_success_all,n_samples);
Nd = zeros(n_success_all,n_samples);
no_attachments = zeros(n_success_all,n_samples);
no_detachments = zeros(n_success_all,n_samples);

n_stickers = Nt*Ns;
n_tethers = Nt-1;
n_chains = Np*(n_tethers+n_stickers);        %Np is the number of pairs 

n_bonds = zeros(n_success_all,Np,n_samples);
n_to_self = zeros(n_success_all,Np,n_samples);
n_to_other = zeros(n_success_all,Np,n_samples);
c_alpha = zeros(n_success_all,Np,n_samples);

% Note - [rx,ry,rz] is resereved for the end-to-end vectors of a full chain
% occuring betweetween two crosslink sites, not individual Kuhn segments
rx = zeros(n_success_all,n_chains,n_samples);
ry = zeros(n_success_all,n_chains,n_samples);
rz = zeros(n_success_all,n_chains,n_samples);
bondtypes = zeros(n_success_all,n_chains,n_samples);

rxr11_st = zeros(n_success_all,n_samples);
rxr22_st = zeros(n_success_all,n_samples);
rxr33_st = zeros(n_success_all,n_samples);
rxr12_st = zeros(n_success_all,n_samples);
rxr23_st = zeros(n_success_all,n_samples);
rxr31_st = zeros(n_success_all,n_samples);
rxr11_st_err = zeros(n_success_all,n_samples);
rxr22_st_err = zeros(n_success_all,n_samples);
rxr33_st_err = zeros(n_success_all,n_samples);
rxr12_st_err = zeros(n_success_all,n_samples);
rxr23_st_err = zeros(n_success_all,n_samples);
rxr31_st_err = zeros(n_success_all,n_samples);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compute = CheckIfNeedToComputeData(OverrideCompute,Controls)

global stress_data_filename msd_data_filename endtoend_data_filename...
    bond_kinetics_filename alignment_data_filename...
    time_stretch_data_filename storeage_loss_filename...
    CalculateStress CalculateEndtoEnd...
    CalculateMSD CalculateBondKinetics CalculateAlignment

compute = 0;
if (~isfile(stress_data_filename) && CalculateStress) ||...
        (~isfile(endtoend_data_filename) && CalculateEndtoEnd) ||...
        (~isfile(msd_data_filename) && CalculateMSD)  ||...
        (~isfile(bond_kinetics_filename) && CalculateBondKinetics)  ||...
        (~isfile(alignment_data_filename) && CalculateAlignment)  ||...
        ~isfile(time_stretch_data_filename)  ||...
        OverrideCompute
    compute = 1;
end
% if Controls.RunOscillatory==1 && (~isfile(storeage_loss_filename) || ...
%         Controls.OverrideStorageLoss==1)
%     compute = 1;
% end

end
