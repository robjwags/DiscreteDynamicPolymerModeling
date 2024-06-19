function CompileData(Package,OverrideCompile,OverrideCompute,ToggleDynamics,...
    LC,DC,CS,CE,CM,CB,CA,CF,OF,NoProcessors)

% Collects data from all LAMMPS runs and compiles it for ensemble averaging
% or collective histogram outputs

% Compiles stress-time responses
% Stress-strain responses

% End-to-end data for histogram binning
% Bond dynamics data

PackageTypes = DefineUniquePackageTypes(Package);

if NoProcessors==1
    LoopCompileData(Package,OverrideCompile,OverrideCompute,ToggleDynamics,...
        LC,DC,CS,CE,CM,CB,CA,CF,OF,NoProcessors)
else
    delete(gcp('nocreate'))
    parpool(NoProcessors)
    parfor n=1:size(PackageTypes,1)
        LoopCompileData(PackageTypes(n,:),OverrideCompile,OverrideCompute,...
            ToggleDynamics,LC,DC,CS,CE,CM,CB,CA,CF,OF,NoProcessors)
    end
    delete(gcp('nocreate'))
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LoopCompileData(Package,OverrideCompile,OverrideCompute,ToggleDynamics,...
    LC,DC,CS,CE,CM,CB,CA,CF,OF,NoProcessors)

% Collects data from all LAMMPS runs and compiles it for ensemble averaging
% or collective histogram outputs

% Compiles stress-time responses
% Stress-strain responses

% End-to-end data for histogram binning
% Bond dynamics data

global LineWidth TurnOnDynamics FontSize...
    LengthConversion DamperConversion RawDataFileName...
    CalculateStress CalculateEndtoEnd...
    CalculateMSD CalculateBondKinetics CalculateAlignment...
    CurrentFolder OutputFolder

LineWidth = 1.5;
FontSize = 20;
TurnOnDynamics = ToggleDynamics;
LengthConversion = LC;
DamperConversion = DC;
% BeadSpringOrMeso = BSOM;
CalculateStress = CS;
CalculateEndtoEnd = CE;
CalculateMSD = CM;
CalculateBondKinetics = CB;
CalculateAlignment = CA;
CurrentFolder = CF;
OutputFolder = OF;

SimsCt = 0;
if NoProcessors==1
    TotalSims = size(Package,1);
    wb1 = waitbar(0,'Compiling all Data...');
else
    SimTypeNo = Package(end-1);
    NoPackageTypes = Package(end);
    disp(['Parameter combo. ',num2str(SimTypeNo),' of ',num2str(NoPackageTypes),' running']);
end

% Samples = unique(Package(:,1));
Nps = unique(Package(:,2));        %Number of molecules
Nss = unique(Package(:,3));        %Number of stickeres per tether site
Ts = unique(Package(:,4));         %Temperature
Nbs = unique(Package(:,5));        %Contour length of chains
kas = unique(Package(:,6));        %Activation energy of association
kds = unique(Package(:,7));        %Activation energy of dissociation
f0s = unique(Package(:,8));       %Force sensitivity to dissociation
dts = unique(Package(:,9));       %timestep size
damps = unique(Package(:,10));    %damping coefficient in units [mass/time]
N_Kuhns = unique(Package(:,12));
bs = unique(Package(:,13));
Separations = unique(Package(:,14));

for i=1:length(dts)
 for j=1:length(Nps)
  for k=1:length(Nss)
   for l=1:length(Ts)
    for m=1:length(Nbs)
     for n=1:length(kas)
      for o=1:length(kds)
       for p=1:length(f0s)
        for q=1:length(damps)
         for r=1:length(N_Kuhns)
          for s=1:length(bs)
           for t=1:length(Separations)

            dt = dts(i);
            Np = Nps(j);
            Ns = Nss(k);
            N = Np*Ns;             
            T = Ts(l);
            Nb = Nbs(m);
            ka = kas(n);
            kd = kds(o);
            f0 = f0s(p);
            damp = damps(q);
            EdgeFactor = Nb;
%             kbT = 293*1.38e-23; % Joules

            N_Kuhn = N_Kuhns(r);
            b = bs(s);
            Separation = Separations(t);

            PackageTemp = Package(ismember(Package(:,2),Np),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,3),Ns),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,4),T),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,5),Nb),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,6),ka),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,7),kd),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,8),f0),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,9),dt),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,10),damp),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,12),N_Kuhn),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,13),b),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,14),Separation),:);

            %% Compile the Raw Data
            for pkg=1:size(PackageTemp,1)
                SimsCt = SimsCt+1;
                if NoProcessors==1
                    PcntComplete = SimsCt/TotalSims;
                    waitbar(PcntComplete,wb1,'Compiling all Data...')
                end

                %% Unpack Swept Input Parameters
                Sample = PackageTemp(pkg,1);    %Sample index

                %% Callout Input Script
                % Initializes non-sweeping parameters
                InputScript(EdgeFactor,Sample,Np,N,ka,kd,f0,dt,damp,...
                    N_Kuhn,b,Separation);

                %% Make diriectories and set filenames
                SetDirAndFileNames;
                DefineCompiledFileNames;

                if ~isfile(RawDataFileName) || OverrideCompile
                    %% Extract data from atoms.dump file
                    timesteps = UnpackTimeSteps;

                    %% Extract data from bonds.dump file
                    Bonds = UnpackBondsDump(timesteps);

                    %% Save File for Sample
                    SaveSampleFiles(Bonds,timesteps);
                end
            end
            
            if ~isempty(PackageTemp)
                %% Check if need to compute outputs
                compute = CheckIfNeedToComputeData(OverrideCompute);
                if compute
                    for samp=1:size(PackageTemp,1)
                        %% Unpack Swept Input Parameters
                        Sample = PackageTemp(pkg,1);    %Sample index

                        %% Callout Input Script
                        % Initializes non-sweeping parameters
                        InputScript(EdgeFactor,Sample,Np,N,ka,kd,f0,dt,damp,...
                            N_Kuhn,b,Separation);

                        %% Make diriectories and set filenames
                        SetDirAndFileNames;
                        DefineCompiledFileNames;

                        %% Import Bonds Data
                        RawData = load(RawDataFileName,'-mat');
                        BondData = RawData.Bonds;
                        Timesteps = RawData.Timesteps;

                        %% Allocate sizes of outputs
                        if samp==1
                            N_steps = size(Timesteps,1);
                            N_samples = size(PackageTemp,1);
                            PreallocateOutputs(N_steps,N_samples);
                        end

                        %% Calculate outputs
                        ComputeOutputs(BondData,Np,dt,samp,OverrideCompute);
                    end
                    SaveComputedOutputs(OverrideCompute);
                end
            end

           end
          end
         end
        end
       end
      end
     end
    end
   end
  end
 end
end
if NoProcessors==1
    close(wb1)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PackageTypes = DefineUniquePackageTypes(Package)

PackageTypes = zeros(1e3,size(Package,2));

Nps = unique(Package(:,2));        %Number of molecules
Nss = unique(Package(:,3));        %Number of stickeres per tether site
Ts = unique(Package(:,4));         %Temperature
Nbs = unique(Package(:,5));        %Contour length of chains
kas = unique(Package(:,6));        %Activation energy of association
kds = unique(Package(:,7));        %Activation energy of dissociation
f0s = unique(Package(:,8));       %Force sensitivity to dissociation
dts = unique(Package(:,9));       %timestep size
damps = unique(Package(:,10));    %damping coefficient in units [mass/time]
N_Kuhns = unique(Package(:,12));
bs = unique(Package(:,13));
Separations = unique(Package(:,14));

ct = 0;
for i=1:length(dts)
 for j=1:length(Nps)
  for k=1:length(Nss)
   for l=1:length(Ts)
    for m=1:length(Nbs)
     for n=1:length(kas)
      for o=1:length(kds)
       for p=1:length(f0s)
        for q=1:length(damps)
         for r=1:length(N_Kuhns)
          for s=1:length(bs)
           for t=1:length(Separations)

            dt = dts(i);
            Np = Nps(j);
            Ns = Nss(k);
            T = Ts(l);
            Nb = Nbs(m);
            ka = kas(n);
            kd = kds(o);
            f0 = f0s(p);
            damp = damps(q);

            N_Kuhn = N_Kuhns(r);
            b = bs(s);
            Separation = Separations(t);

            PackageTemp = Package(ismember(Package(:,2),Np),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,3),Ns),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,4),T),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,5),Nb),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,6),ka),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,7),kd),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,8),f0),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,9),dt),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,10),damp),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,12),N_Kuhn),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,13),b),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,14),Separation),:);

            if ~isempty(PackageTemp)
                ct = ct+1;
                PackageTypes(ct,:) = PackageTemp;
            end

           end
          end
         end
        end
       end
      end
     end
    end
   end
  end
 end
end

PackageTypes(PackageTypes(:,1)==0,:) = [];
PackageTypes(:,end+1) = (1:size(PackageTypes,1))';
PackageTypes(:,end+1) = size(PackageTypes,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveComputedOutputs(OverrideCompute)

global CalculateBondKinetics BondKineticsDataFileName...
    time ka_out ka1_out kd_out Na Nd...
    attached_bond_lifetimes detached_bond_lifetimes tau_a1

if CalculateBondKinetics && (~isfile(BondKineticsDataFileName) || OverrideCompute)
    BondKinetics.time = time;
    BondKinetics.ka = ka_out;
    BondKinetics.ka1 = ka1_out;
    BondKinetics.kd = kd_out;
    BondKinetics.Na = Na;
    BondKinetics.Nd = Nd;
    BondKinetics.FirstAttachment = tau_a1;
    BondKinetics.AttachedLifetimes = attached_bond_lifetimes;
    BondKinetics.DetachedLifetimes = detached_bond_lifetimes;

    save(BondKineticsDataFileName,'-struct','BondKinetics')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ComputeOutputs(BondData,Np,dt,samp,OverrideCompute)

global CalculateBondKinetics BondKineticsDataFileName time...
    ka_out ka1_out kd_out Na Nd...
    attached_bond_lifetimes detached_bond_lifetimes tau_a1

% timestep id bond_type node1 node2 rx ry rz norm force
timesteps = unique(BondData(:,1));

%% Initialize dynamics variables
timestep0 = 0;
dynamic_pairs0 = [];
all_dynamic_pairs = [];
N_bound0 = 0;
N_open0 = 0;
N_virgin = Np;

%% Initialize dynamics variables
wb2 = waitbar(0,'Computing dynamics data...');
for i=1:length(timesteps)
    timestep = timesteps(i);
    time(i) = timestep*dt;
    waitbar(i/length(timesteps),wb2,...
        'Computing dynamics data...')
    BondDataTemp = BondData(BondData(:,1)==timestep,:);

    %% Compute bond dynamics rates
    if CalculateBondKinetics && (~isfile(BondKineticsDataFileName) || OverrideCompute)
        dynbondtype = 1;
        [kd_out(i,samp),ka_out(i,samp),ka1_out(i,samp),...
            dynamic_pairs0,N_bound0,N_open0,timestep0,all_dynamic_pairs,...
            N_virgin] = ...
            ComputeBondDynamics(timestep,i,dt,BondDataTemp,...
            dynamic_pairs0,N_bound0,N_open0,...
            timestep0,Np,dynbondtype,all_dynamic_pairs,N_virgin);

        Na(i,samp) = N_bound0;
        Nd(i,samp) = N_open0;
    end
end
close(wb2)

% Calculate the time to first encounter, as well as the subsequent
% lifetimes

% Sort the pairs and reassemble
if ~isempty(all_dynamic_pairs)
    just_pairs = all_dynamic_pairs(:,1:2);
    just_times = all_dynamic_pairs(:,3:4);
    just_pairs = sort(just_pairs,2);
    all_dynamic_pairs = [just_pairs,just_times]; % append times back onto pairs

    % Define unique pairs
    unique_pairs = unique(all_dynamic_pairs(:,1:2),'rows');

    % Initialize time to first attachment
    tau_a1 = zeros(size(unique_pairs,1),1);
    attached_bond_lifetimes = [];
    detached_bond_lifetimes = [];
    for pair_indx = 1:size(unique_pairs,1)
        pair = unique_pairs(pair_indx,:);

        % Identify rows belonging to unique pair
        IndxL = find(ismember(all_dynamic_pairs(:,1),pair(1)));
        IndxR = find(ismember(all_dynamic_pairs(:,2),pair(2)));
        Indx = intersect(IndxL,IndxR);

        % Isolate the instances of this pair
        pairs_temp = all_dynamic_pairs(Indx,:);

        % Sort the rows by event time (note that first event and all odd rows
        % must be attachment events while the evens are detachments)
        pairs_temp = sortrows(pairs_temp,3);

        % Isolate attachment and detachment events
        tau_a1(pair_indx,1) = pairs_temp(1,3); %time to first attachment

        all_event_times = pairs_temp(:,3:end);
        all_event_durations = diff(all_event_times,1,1);

        attached_bond_lifetimes = [attached_bond_lifetimes;all_event_durations(1:2:end,1)];
        detached_bond_lifetimes = [detached_bond_lifetimes;all_event_durations(2:2:end,1)];
    end
else
    attached_bond_lifetimes =[];
    detached_bond_lifetimes = [];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kd,ka,ka1,dynamic_pairs0,N_bound0,N_open0,timestep0,...
    all_dynamic_pairs,N_virgin_bonds] = ...
    ComputeBondDynamics(timestep,i,dt,BondsTemp,...
    dynamic_pairs0,N_bound0,N_open0,timestep0,Np,dynbondtype,...
    all_dynamic_pairs,N_virgin_bonds)

Dt = (timestep-timestep0)*dt;
Time = timestep*dt;

if ~isempty(BondsTemp)
    % Isolate dynamic bonds
    pairs = BondsTemp(:,4:5);
    bondtypes = BondsTemp(:,3);
    dynamic_pairs = pairs(bondtypes==dynbondtype,:);

    % Delete redundant dyanmic pairs
    dynamic_pairs = sort(dynamic_pairs,2);
    dynamic_pairs = sortrows(dynamic_pairs,1);
    dynamic_pairs(2:2:end,:) = [];

    % compute number of bound and free stickers
    N_pairs = Np;
    N_bound = size(dynamic_pairs,1);

    N_open = N_pairs-N_bound;

    % ID previously existing and new dynamic bonds
    if isempty(dynamic_pairs0)
        %     old_pairs = [];
        new_pairs = dynamic_pairs;
    else
        old_nodes1 = ismember(dynamic_pairs(:,1),dynamic_pairs0(:,1));
        old_nodes2 = ismember(dynamic_pairs(:,2),dynamic_pairs0(:,2));
        indx = old_nodes1 + old_nodes2; %new if this sums to 2
        new_pairs = dynamic_pairs(indx~=2,:);
    end


    % Find first-time new pairs and define ka1
    if ~isempty(new_pairs) && isempty(all_dynamic_pairs)
        no_first_time_attachments = size(new_pairs,1);
        ka1 = no_first_time_attachments/N_virgin_bonds/Dt;
    elseif ~isempty(new_pairs) && ~isempty(all_dynamic_pairs)
        np = sort(new_pairs,2);
        ap = sort(all_dynamic_pairs(:,1:2),2);

        indxL = find(ismember(np(:,1),ap(:,1)));
        indxR = find(ismember(np(:,2),ap(:,2)));
        indx = intersect(indxL,indxR);
        np(indx,:) = [];

        no_first_time_attachments = size(np,1);
        ka1 = no_first_time_attachments/N_virgin_bonds/Dt;
    else
        no_first_time_attachments = 0;
        ka1 = 0;
    end
    N_virgin_bonds = N_virgin_bonds-no_first_time_attachments;

    % Store all pairs history
    if ~isempty(new_pairs) % Add the attachment events
        all_dynamic_pairs = [all_dynamic_pairs;[new_pairs,...
            Time*ones(size(new_pairs,1),1),...
            ones(size(new_pairs,1),1)]];
    end

    no_attachments = size(new_pairs,1);   

    if i==1
        ka = 0;
    else
        ka = no_attachments/N_open0/Dt;
    end

    if ka<0 || N_open0<0 || Dt<0
        warning('ka<0, should always be >= 0')
    end

    % ID lost bonds (i.e., detachments)
    if isempty(dynamic_pairs0)
        %     surviving_pairs = [];
        deleted_pairs = [];
    else
        surviving_nodes1 = ismember(dynamic_pairs0(:,1),dynamic_pairs(:,1));
        surviving_nodes2 = ismember(dynamic_pairs0(:,2),dynamic_pairs(:,2));
        indx = surviving_nodes1 + surviving_nodes2;
        deleted_pairs = dynamic_pairs0(indx~=2,:);
    end

    no_detachments =  size(deleted_pairs,1);

    if i==1 || N_bound0==0
        kd = 0;
    else
        kd = no_detachments/N_bound0/Dt;
    end

    if ~isempty(deleted_pairs) % Add the deletion events
        all_dynamic_pairs = [all_dynamic_pairs;[deleted_pairs,...
            Time*ones(size(deleted_pairs,1),1),...
            zeros(size(deleted_pairs,1),1)]];
    end

    % Set old values
    dynamic_pairs0 = dynamic_pairs;
    N_bound0 = N_bound;
    N_open0 = N_open;
    timestep0 = timestep;
else
    ka = 0;
    kd = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PreallocateOutputs(N_steps,N_samples)

global bondconc time ka_out ka1_out kd_out Na Nd...

time = zeros(N_steps,1);
bondconc = zeros(N_steps,N_samples);

% Kinetics
ka_out = zeros(N_steps,N_samples);
ka1_out = zeros(N_steps,N_samples);
kd_out = zeros(N_steps,N_samples);
Na = zeros(N_steps,N_samples);
Nd = zeros(N_steps,N_samples);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compute = CheckIfNeedToComputeData(OverrideCompute)

global BondKineticsDataFileName CalculateBondKinetics

compute = 0;
if (~isfile(BondKineticsDataFileName) && CalculateBondKinetics) ||...
        OverrideCompute
    compute = 1;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveSampleFiles(Bonds,timesteps)

global RawDataFileName 

Raw.Bonds = Bonds;
Raw.Timesteps = timesteps;
save(RawDataFileName,'-struct','Raw','-v7.3')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Bonds = UnpackBondsDump(timestep)

global OutputBond

N_output = length(timestep);

% Extract periodic domain boundaries
fmt1=repmat('%f',1,6);             % first data section format
fid=fopen(OutputBond,'r');     % open file
l=fgetl(fid);                      % read first line
ct = 0;
wb = waitbar(0,'Extracting bonds data...');
while ~feof(fid)                   % loop through file by record
    if strfind(l,'ITEM: ENTRIES index c_Pairs[1] c_Pairs[2] c_Pairs[3] c_Bonds[1] c_Bonds[2] c_Bonds[3] c_Bonds[4] c_Bonds[5]')
        ct = ct+1;
        Bonds_array{:,:,ct} = cell2mat(textscan(fid,'%f'));
        if ~mod(ct,50)
            PcntComplete = ct/N_output;
            waitbar(PcntComplete,wb,'Extracting bonds data...')
        end
    end                                 % break when find first section
    l=fgetl(fid);                       % next record
end
close(wb)

% Unpack bond information
N_var = 9;      % id bond_type node1 node2 rx ry rz norm force
Bonds = zeros(1,N_var+1);
ct = 0;
wb = waitbar(0,'Unpacking bond data...');
StartIndx = 0;
for i=1:N_output-1
    if ~mod(i,50)
        PcntComplete = i/N_output;
        waitbar(PcntComplete,wb,'Unpacking bond data...')
    end
    timestep_temp = timestep(i);
    Bonds_temp = Bonds_array{:,:,i};
    N_bonds = length(Bonds_temp)/N_var;
    Range = StartIndx+(1:N_bonds)';

    ID_indx = (1:N_var:size(Bonds_temp,1))';
    type_indx = ID_indx+1;
    node1_indx = type_indx+1;
    node2_indx = node1_indx+1;
    rx_indx = node2_indx+1;
    ry_indx = rx_indx+1;
    rz_indx = ry_indx+1;
    norm_indx = rz_indx+1;
    force_indx = norm_indx+1;

    Bonds(Range,1) = timestep_temp*ones(size(Range));
    Bonds(Range,2) = Bonds_temp(ID_indx);
    Bonds(Range,3) = Bonds_temp(type_indx);
    Bonds(Range,4) = Bonds_temp(node1_indx);
    Bonds(Range,5) = Bonds_temp(node2_indx);
    Bonds(Range,6) = Bonds_temp(rx_indx);
    Bonds(Range,7) = Bonds_temp(ry_indx);
    Bonds(Range,8) = Bonds_temp(rz_indx);
    Bonds(Range,9) = Bonds_temp(norm_indx);
    Bonds(Range,10) = Bonds_temp(force_indx);

    StartIndx = size(Bonds,1);
end
close(wb)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function timestep = UnpackTimeSteps

global OutputAtom

Corners = []; Atoms = [];

% Extract timesteps
fmt1=repmat('%f',1,6);             % first data section format
fid=fopen(OutputAtom,'r');          % open file
l=fgetl(fid);                      % read first line
ct = 0;
while ~feof(fid)                   % loop through file by record
    if strfind(l,'ITEM: TIMESTEP')
        ct = ct+1;
        timestep(ct,1) = cell2mat(textscan(fid,'%f'));
    end                                 % break when find first section
    l=fgetl(fid);                       % next record
end
N_output = length(timestep);

end
