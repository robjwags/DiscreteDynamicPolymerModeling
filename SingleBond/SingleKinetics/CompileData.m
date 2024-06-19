function T_comp = CompileData(Package,OverrideCompile,OverrideCompute,ToggleDynamics,...
    LC,DC,BSOM,CS,CE,CM,CB,CA,CF,OF,NoProcessors)

% Collects data from all LAMMPS runs and compiles it for ensemble averaging
% or collective histogram outputs

% Compiles stress-time responses
% Stress-strain responses

% End-to-end data for histogram binning
% Bond dynamics data

PackageTypes = DefineUniquePackageTypes(Package);

if NoProcessors==1
    T_comp = LoopCompileData(Package,OverrideCompile,OverrideCompute,ToggleDynamics,...
        LC,DC,BSOM,CS,CE,CM,CB,CA,CF,OF,NoProcessors);
else
    tic
    delete(gcp('nocreate'))
    parpool(NoProcessors)
    parfor n=1:size(PackageTypes,1)
        [~] = LoopCompileData(PackageTypes(n,:),OverrideCompile,OverrideCompute,...
            ToggleDynamics,LC,DC,BSOM,CS,CE,CM,CB,CA,CF,O,NoProcessors);
    end
    delete(gcp('nocreate'))
    T_comp = toc;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T_comp = LoopCompileData(Package,OverrideCompile,OverrideCompute,ToggleDynamics,...
    LC,DC,BSOM,CS,CE,CM,CB,CA,CF,OF,NoProcessors)

% Collects data from all LAMMPS runs and compiles it for ensemble averaging
% or collective histogram outputs

% Compiles stress-time responses
% Stress-strain responses

% End-to-end data for histogram binning
% Bond dynamics data

global LineWidth TurnOnDynamics FontSize...
    LengthConversion DamperConversion...
    BeadSpringOrMeso RawDataFileName...
    CalculateStress CalculateEndtoEnd...
    CalculateMSD CalculateBondKinetics CalculateAlignment...
    CurrentFolder OutputFolder OutputAtom OutputBond

LineWidth = 1.5;
FontSize = 20;
TurnOnDynamics = ToggleDynamics;
LengthConversion = LC;
DamperConversion = DC;
BeadSpringOrMeso = BSOM;
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

T = zeros(1,size(Package,1));

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
            N = Np*Ns;             
            T = Ts(l);
            Nb = Nbs(m);
            ka = kas(n);
            kd = kds(o);
            f0 = f0s(p);
            damp = damps(q);
            EdgeFactor = Nb;
            kbT = 293*1.38e-23; % Joules

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
                tic
                ct = ct+1;
                D = unique(PackageTemp(11));
                SimsCt = SimsCt+1;
                if NoProcessors==1
                    PcntComplete = SimsCt/TotalSims;
                    waitbar(PcntComplete,wb1,'Compiling all Data...')
                end
%                 else
%                     if ~mod(SimsCt)
%                     disp(['Parameter combo. ',num2str(SimTypeNo),' of ',...
%                         num2str(NoPackageTypes),' is ',...
%                         num2str(PcntComplete*100,'%.2f'),'% complete']);
%                 end

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
                    [timesteps,Corners,Atoms] = UnpackAtomsDump;

                    %% Extract data from bonds.dump file
                    Bonds = UnpackBondsDump(timesteps);

                    %% Save File for Sample
                    SaveSampleFiles(Corners,Atoms,Bonds);
                end
                
                % Delete raw data once replaced by compressed MATLAB format
                if isfile(OutputAtom) && isfile(RawDataFileName)
                    fclose all;
                    delete(OutputAtom)
                end
                if isfile(OutputBond) && isfile(RawDataFileName)
                    fclose all;
                    delete(OutputBond)
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
                        Corners = RawData.Corners;
                        BondData = RawData.Bonds;
                        AtomData = RawData.Atoms;


                        %% Allocate sizes of outputs
                        if samp==1
                            N_steps = size(Corners,1);
                            N_samples = size(PackageTemp,1);
%                             N_seg = length(unique(BondData(:,2)));
                            PreallocateOutputs(N_steps,N_samples,Np,N_Kuhn);
                        end

                        %% Calculate outputs
                        ComputeOutputs(Corners,BondData,AtomData,...
                            Np,N_Kuhn,b,kbT,ka,D,Separation,dt,samp,...
                            OverrideCompute);
                    end
                    SaveComputedOutputs(OverrideCompute);
                end
                T_comp(1,ct) = toc;
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
            N = Np*Ns;             
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

global CalculateStress CalculateEndtoEnd...
    CalculateMSD CalculateBondKinetics CalculateAlignment...
    StressDataFileName MSDDataFileName...
    EndToEndDataFileName BondKineticsDataFileName AlignmentDataFileName...
    TimeStretchDataFileName...
    time stretch sig11 sig22 sig33 sig12 sig23 sig31...
    bondconc...
    msd_st msd_st_err...
    msd_st_r0 msd_st_r0_err...
    msd_st_r00 msd_st_r00_err...
    msd_st_t0 msd_st_t0_err...
    msd_st_t00 msd_st_t00_err...
    ka_out ka1_out kd_out k_rpt Na Nd...
    rx ry rz bondtypes...
    rxr11_st rxr22_st rxr33_st rxr12_st rxr23_st rxr31_st...
    rxr11_st_err rxr22_st_err rxr33_st_err rxr12_st_err rxr23_st_err rxr31_st_err...
    attached_bond_lifetimes detached_bond_lifetimes tau_a1


if (~isfile(TimeStretchDataFileName) || OverrideCompute)
    TimeAndStretch.time = time;
    TimeAndStretch.stretch = stretch;
    save(TimeStretchDataFileName,'-struct','TimeAndStretch')
end

if CalculateStress && (~isfile(StressDataFileName) || OverrideCompute)
    Stress.sig11 = sig11;
    Stress.sig22 = sig22;
    Stress.sig33 = sig33;
    Stress.sig12 = sig12;
    Stress.sig23 = sig23;
    Stress.sig31 = sig31;
    Stress.bondconc = bondconc;

    save(StressDataFileName,'-struct','Stress')
end

if CalculateEndtoEnd && (~isfile(EndToEndDataFileName) || OverrideCompute)
    EndtoEnd.rx = rx;
    EndtoEnd.ry = ry;
    EndtoEnd.rz = rz;
    EndtoEnd.bondtypes = bondtypes;

    save(EndToEndDataFileName,'-struct','EndtoEnd')
end

if CalculateMSD && (~isfile(MSDDataFileName) || OverrideCompute)
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

    save(MSDDataFileName,'-struct','MSD')
end

if CalculateBondKinetics && (~isfile(BondKineticsDataFileName) || OverrideCompute)
    BondKinetics.ka = ka_out;
    BondKinetics.ka1 = ka1_out;
    BondKinetics.k_rpt = k_rpt;
    BondKinetics.kd = kd_out;
    BondKinetics.Na = Na;
    BondKinetics.Nd = Nd;
    BondKinetics.FirstAttachment = tau_a1;
    BondKinetics.AttachedLifetimes = attached_bond_lifetimes;
    BondKinetics.DetachedLifetimes = detached_bond_lifetimes;

    save(BondKineticsDataFileName,'-struct','BondKinetics')
end

if CalculateAlignment && (~isfile(AlignmentDataFileName) || OverrideCompute)
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

    save(AlignmentDataFileName,'-struct','RxR')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ComputeOutputs(Corners,BondData,AtomData,Np,N_Kuhn,b,kbT,ka,D,...
    dist,dt,samp,...
    OverrideCompute)

global CalculateStress CalculateEndtoEnd...
    CalculateMSD CalculateBondKinetics CalculateAlignment...
    StressDataFileName MSDDataFileName...
    EndToEndDataFileName BondKineticsDataFileName AlignmentDataFileName...
    BeadSpringOrMeso sig11 sig22 sig33 sig12 sig23 sig31...
    bondconc time stretch...
    msd_st msd_st_err...
    msd_st_r0 msd_st_r0_err...
    msd_st_r00 msd_st_r00_err...
    msd_st_t0 msd_st_t0_err...
    msd_st_t00 msd_st_t00_err...
    ka_out ka1_out kd_out k_rpt Na Nd...
    rx ry rz bondtypes...
    rxr11_st rxr22_st rxr33_st rxr12_st rxr23_st rxr31_st...
    rxr11_st_err rxr22_st_err rxr33_st_err rxr12_st_err rxr23_st_err rxr31_st_err...
    attached_bond_lifetimes detached_bond_lifetimes tau_a1...
    LengthConversion RawDataFileName

tau0 = (b*LengthConversion)^2/D;
eaStar = -log(ka*tau0);

% timestep id bond_type node1 node2 rx ry rz norm force
timesteps = unique(BondData(:,1));
if size(Corners,1)~=length(timesteps)
    Diff = size(Corners,1)-length(timesteps);
    if Diff<0
        timesteps(end+Diff+1:end) = [];
    end
end

% if CalculateBondKinetics && (~isfile(BondKineticsDataFileName) || OverrideCompute)
%     dat = load(RawDataFileName,'-mat','bid','batom1','batom2');
%     p1 = dat.batom1;
%     p2 = dat.batom2;
%     bond_id = dat.bid;
% end

StartMSDMeasurePct = 0.001;
MSDStartIndx = ceil(StartMSDMeasurePct/100*length(timesteps));
MSDStartStep = timesteps(MSDStartIndx);
%% Initialize dynamics variables
timestep0 = 0;
dynamic_pairs0 = [];
all_dynamic_pairs = [];
N_bound0 = 0;
N_free0 = 0;
N_virgin = Np;

%% Initialize dynamics variables
wb2 = waitbar(0,'Computing stress, end-to-end, MSD, & dynamics data...');
for i=1:length(timesteps)
    timestep = timesteps(i);
    if i==1
        Lx0 = abs(Corners(2)-Corners(5));
    end
    time(i) = timestep*dt;
    stretch(i) = abs(Corners(2)-Corners(5))/Lx0;
    waitbar(i/length(timesteps),wb2,...
        'Computing stress, end-to-end, MSD, & dynamics data...')
    rs = []; AtomDataTemp = []; BondDataTemp = [];

    %% Store end-to-end components & Bond types
    if CalculateEndtoEnd && (~isfile(EndToEndDataFileName) || OverrideCompute)
        AtomDataTemp = AtomData(AtomData(:,1)==timestep,:);
        BondDataTemp = BondData(BondData(:,1)==timestep,:);
        rs = ExtractEndToEnd(AtomDataTemp,BondDataTemp,Np,N_Kuhn,i,samp);
        rx(i,:,samp) = rs(:,1); 
        ry(i,:,samp) = rs(:,2); 
        rz(i,:,samp) = rs(:,3);
    end
    
    %% Compute metric tensors by bond type
    if CalculateAlignment && (~isfile(AlignmentDataFileName) || OverrideCompute)
        if isempty(AtomDataTemp) || isempty(BondDataTemp)
            AtomDataTemp = AtomData(AtomData(:,1)==timestep,:);
            BondDataTemp = BondData(BondData(:,1)==timestep,:);
        end
        if isempty(rs)
            rs = ExtractEndToEnd(AtomDataTemp,BondDataTemp,Np,N_Kuhn,i,samp);
        end
        [rxr11_st(i,samp),rxr22_st(i,samp),rxr33_st(i,samp),...
            rxr12_st(i,samp),rxr23_st(i,samp),rxr31_st(i,samp),...
            rxr11_st_err(i,samp),rxr22_st_err(i,samp),rxr33_st_err(i,samp),...
            rxr12_st_err(i,samp),rxr23_st_err(i,samp),rxr31_st_err(i,samp)] = ...
            ComputeMetricTensor(rx(i,:,samp),ry(i,:,samp),rz(i,:,samp));
    end

    %% Compute stress
    if CalculateStress && (~isfile(StressDataFileName) || OverrideCompute)
        if isempty(AtomDataTemp) || isempty(BondDataTemp)
            AtomDataTemp = AtomData(AtomData(:,1)==timestep,:);
            BondDataTemp = BondData(BondData(:,1)==timestep,:);
        end
        forces = ComputeForces(bondtypes(i,:,samp),rs,N_Kuhn,b,kbT);
        [sigma,bondconc(i,samp)] = ComputeStressAndBondConc(forces,rs,Corners(i,:));
        sig11(i,samp) = sigma(1);
        sig22(i,samp) = sigma(2);
        sig33(i,samp) = sigma(3);
        sig12(i,samp) = sigma(4);
        sig23(i,samp) = sigma(5);
        sig31(i,samp) = sigma(6);
    end

    %% Compute MSDs (by atom type)
    if CalculateMSD && (~isfile(MSDDataFileName) || OverrideCompute) &&...
            BeadSpringOrMeso~=2
        if isempty(AtomDataTemp) || isempty(BondDataTemp)
            AtomDataTemp = AtomData(AtomData(:,1)==timestep,:);
            BondDataTemp = BondData(BondData(:,1)==timestep,:);
        end
        if timestep<=MSDStartStep
            msd_st(i) = 0;
            msd_st_err(i) = 0;
            PrevStep = MSDStartStep;
        else
            if BeadSpringOrMeso==0
                stickertype = 3;   %End groups (MSD should follow Rouse diffusion)
                tethertype = 1;
            elseif BeadSpringOrMeso==1
                stickertype = 2;
                tethertype = 1;
            end
            [msd_st(i,samp),msd_st_err(i,samp),...
                msd_st_r0(i,samp),msd_st_r0_err(i,samp),...
                msd_st_t0(i,samp),msd_st_t0_err(i,samp),...
                msd_st_r00(i,samp),msd_st_r00_err(i,samp),...
                msd_st_t00(i,samp),msd_st_t00_err(i,samp),PrevStep] = ...
                ComputeMSDs(timestep,AtomData,AtomDataTemp,Corners(i,:),...
                MSDStartStep,PrevStep,stickertype,tethertype);
        end
    end

    %% Compute bond dynamics rates
    if CalculateBondKinetics && (~isfile(BondKineticsDataFileName) || OverrideCompute)
        if isempty(AtomDataTemp) || isempty(BondDataTemp)
            AtomDataTemp = AtomData(AtomData(:,1)==timestep,:);
            BondDataTemp = BondData(BondData(:,1)==timestep,:);
        end
        if BeadSpringOrMeso==0
            tiptyp = 3;     % type of the stickers
            dynbondtype = 2;
        elseif BeadSpringOrMeso==1
            tiptyp = 3;     % type of the stickers
            dynbondtype = 3;
        elseif BeadSpringOrMeso==2
            tiptyp = 1;     % type of the stickers
            dynbondtype = 1;
        end
        if isempty(AtomDataTemp) || isempty(BondDataTemp)
            AtomDataTemp = AtomData(AtomData(:,1)==timestep,:);
            BondDataTemp = BondData(BondData(:,1)==timestep,:);
        end
        [kd_out(i,samp),ka_out(i,samp),ka1_out(i,samp),k_rpt(i,samp),...
            dynamic_pairs0,N_bound0,N_free0,timestep0,all_dynamic_pairs,...
            N_virgin] = ...
            ComputeBondDynamics(timestep,i,dt,BondDataTemp,...
            AtomDataTemp,dynamic_pairs0,N_bound0,N_free0,...
            timestep0,tiptyp,dynbondtype,all_dynamic_pairs,N_virgin,BeadSpringOrMeso);
        Na(i,samp) = N_bound0;
        Nd(i,samp) = N_free0;
    end


    % 
    % if CalculateBondKinetics && (~isfile(BondKineticsDataFileName) || OverrideCompute)
    %     tiptyp = 2;
    %     dynbondtype = 2;
    %     if i==1
    %         dt_temp = dt;
    %         Dt = 0;
    %     else
    %         dt_temp = (time(i)-time(i-1))/(timesteps(i)-timesteps(i-1));
    %         Dt = time(i)-time(i-1);
    %     end
    %     [kd_out(i,samp),ka_out(i,samp),ka1_out(i,samp),...
    %         dynamic_pairs0,N_bound0,N_free0,timestep0,all_dynamic_pairs,...
    %         N_virgin,no_attachments(i,samp),no_detachments(i,samp)] = ...
    %         ComputeBondDynamics(timestep,i,dt_temp,Dt,...
    %         p1,p2,bond_types,bond_id,...
    %         types{i},dynamic_pairs0,...
    %         N_bound0,N_free0,...
    %         timestep0,tiptyp,dynbondtype,all_dynamic_pairs,N_virgin);
    %     Na(i,samp) = N_bound0;
    %     Nd(i,samp) = N_free0;
    % end
end

if ~isempty(all_dynamic_pairs)
    all_dynamic_pairs(:,1:2) = sort(all_dynamic_pairs(:,1:2),2);

    ka_check = zeros(size(timesteps));
    kd_check = zeros(size(timesteps));
    ka1_check = zeros(size(timesteps));
    krpt_check = zeros(size(timesteps));
    N_attached = zeros(size(timesteps));
    N_detached = zeros(size(timesteps));
    N_d_virgin = zeros(size(timesteps));
    N_d_priors = zeros(size(timesteps));
    N_attached(1) = Na(1);
    N_detached(1) = Nd(1);
    N_d_virgin(1) = Nd(1);
    for i=2:length(timesteps)
        Dt = time(i)-time(i-1);
        temp_time = time(i);
        % Nd_temp = Nd(i);
        % Na_temp = Na(i);

        % pull current and prior pairs. Identify unique pairs
        temp_data = all_dynamic_pairs(all_dynamic_pairs(:,3)==temp_time,:); % data at step i
        prior_data = all_dynamic_pairs(all_dynamic_pairs(:,3)<temp_time,:); % all data before step i
        prior_pairs = prior_data(:,1:2); % pairs that existed prior (might still exist)
        [~,unique_prior_pairs_indx] = unique(prior_pairs(:,1:2),'rows');
        unique_prior_pairs = prior_pairs(unique_prior_pairs_indx,:); % unique sets of prior_pairs

        uptonow_data = [prior_data;temp_data]; % all data up to and including step i
        uptonow_pairs = uptonow_data(:,1:2);
        [~,unique_uptonow_pairs_indx] = unique(uptonow_pairs(:,1:2),'rows');
        unique_uptonow_pairs = uptonow_pairs(unique_uptonow_pairs_indx,:); % unique sets of prior_pairs

        % Find wherever repeat events occured
        repeats = temp_data(ismember(temp_data(:,1:2),unique_prior_pairs(:,1:2),'rows'),:); %
        % if ~isempty(repeats)
        % disp(repeats);
        % end
        % repeats = temp_data(repeat_indices);
        % [~,repeat_indices] = ismember(unique_prior_pairs(:,1:2),temp_data(:,1:2),'rows'); %
        % repeats = temp_data(repeat_indices);

        % Find how many bonds are detached that were attached prior
        N_d_repeats = 0;
        for j=1:size(unique_prior_pairs)
            check_dat = uptonow_data(ismember(uptonow_data(:,1:2),unique_prior_pairs(j,1:2),'rows'),:);
            % check_dat = uptonow_data(indx_membs,:);
            if ~mod(size(check_dat,1),2) % if size of check_dat is even, then ended on a detachment event and N_d_repeats = N_d_repeats+1;
                N_d_repeats = N_d_repeats+1;
            end
        end

        attachments = temp_data(temp_data(:,end)==1,:);
        detachments = temp_data(temp_data(:,end)==0,:);

        n_attachments = size(attachments,1); % number of attachment events occuring at step i
        n_detachments = size(detachments,1); % number of detachment events occuring at step i
        n_repeats = size(repeats(repeats(:,end)==1,:),1); % number of repeat attachment events occuring at step i
        n_first = n_attachments-n_repeats; % number of first-time attachment events occuring at step i

        N_attached(i) = N_attached(i-1)+n_attachments-n_detachments; % number of bonds attached at step i
        N_detached(i) = N_detached(1)-N_attached(i); % number of bonds detached at step i
        % N_detached(i) = N_detached(i-1)+n_detachments-n_attachments;
        N_d_virgin(i) = N_detached(1)-size(unique_uptonow_pairs,1); % number of detached bonds that have never been attached
        N_d_priors(i) = N_d_repeats; % number of detached bonds that have been attached prior

        ka_check(i) = n_attachments/N_detached(i)/Dt;
        krpt_check(i) = n_repeats/N_d_priors(i)/Dt;
        ka1_check(i) = n_first/N_d_virgin(i)/Dt;
        kd_check(i) = n_detachments/N_attached(i)/Dt;
    end
    % k_rpt = krpt_check;
    k_rpt(isnan(k_rpt)) = 0;
end

check_fig = 0;
if check_fig==1
    figure(1); clf; hold on
    k0 = (3/2/pi/N_Kuhn)^(3/2)*exp(-eaStar);
    l0 = sqrt(2/3*N_Kuhn*b^2);
    smooth_ka = smooth(ka_out,100);
    ka_in = k0/tau0*exp(-dist^2/l0^2);
    scatter(time,ka_out)
    scatter(time,smooth_ka)
    plot(time,ka_in*ones(size(time)),'k--')
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
    %     attachment_times = pairs_temp(pairs_temp(:,end)==1,3); %should have all 1's at end
    %     detachment_times = pairs_temp(pairs_temp(:,end)==0,3); %should have all 0's at end

    all_event_times = pairs_temp(:,3:end);
    all_event_durations = diff(all_event_times,1,1);

    attached_bond_lifetimes = [attached_bond_lifetimes;all_event_durations(1:2:end,1)];
    detached_bond_lifetimes = [detached_bond_lifetimes;all_event_durations(2:2:end,1)];
end
% ka1 = 1./tau_a1;
else
    attached_bond_lifetimes =[];
    detached_bond_lifetimes = [];
end


% Crop data, as needed
rem_indx = find(time(2:end)==0);
CropOutputs(rem_indx,samp);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CropOutputs(rem_indx,samp)

global sig11 sig22 sig33 sig12 sig23 sig31...
    bondconc time stretch...
    msd_st msd_st_err...
    msd_st_r0 msd_st_r0_err...
    msd_st_r00 msd_st_r00_err...
    msd_st_t0 msd_st_t0_err...
    msd_st_t00 msd_st_t00_err...
    ka_out ka1_out kd_out k_rpt Na Nd...
    rx ry rz bondtypes...
    rxr11_st rxr22_st rxr33_st rxr12_st rxr23_st rxr31_st...
    rxr11_st_err rxr22_st_err rxr33_st_err rxr12_st_err rxr23_st_err rxr31_st_err...

sig11(rem_indx,samp) = NaN;
sig22(rem_indx,samp) = NaN;
sig33(rem_indx,samp) = NaN; 
sig12(rem_indx,samp) = NaN; 
sig23(rem_indx,samp) = NaN; 
sig31(rem_indx,samp) = NaN;
bondconc(rem_indx,samp) = NaN; 
time(rem_indx,samp) = NaN; 
stretch(rem_indx,samp) = NaN;
msd_st(rem_indx,samp) = NaN; 
msd_st_err(rem_indx,samp) = NaN;
msd_st_r0(rem_indx,samp) = NaN; 
msd_st_r0_err(rem_indx,samp) = NaN;
msd_st_r00(rem_indx,samp) = NaN; 
msd_st_r00_err(rem_indx,samp) = NaN;
msd_st_t0(rem_indx,samp) = NaN; 
msd_st_t0_err(rem_indx,samp) = NaN;
msd_st_t00(rem_indx,samp) = NaN; 
msd_st_t00_err(rem_indx,samp) = NaN;
ka_out(rem_indx,samp) = NaN; 
ka1_out(rem_indx,samp) = NaN; 
k_rpt(rem_indx,samp) = NaN; 
kd_out(rem_indx,samp) = NaN; 
Na(rem_indx,samp) = NaN; 
Nd(rem_indx,samp) = NaN;
rx(rem_indx,samp) = NaN; 
ry(rem_indx,samp) = NaN; 
rz(rem_indx,samp) = NaN; 
bondtypes(rem_indx,samp) = NaN;
rxr11_st(rem_indx,samp) = NaN; 
rxr22_st(rem_indx,samp) = NaN; 
rxr33_st(rem_indx,samp) = NaN; 
rxr12_st(rem_indx,samp) = NaN; 
rxr23_st(rem_indx,samp) = NaN; 
rxr31_st(rem_indx,samp) = NaN;
rxr11_st_err(rem_indx,samp) = NaN; 
rxr22_st_err(rem_indx,samp) = NaN; 
rxr33_st_err(rem_indx,samp) = NaN; 
rxr12_st_err(rem_indx,samp) = NaN; 
rxr23_st_err(rem_indx,samp) = NaN; 
rxr31_st_err(rem_indx,samp) = NaN;

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [kd,ka,ka1,dynamic_pairs0,N_bound0,N_free0,timestep0,...
%     all_dynamic_pairs,N_virgin_bonds,no_attachments,no_detachments] = ...
%     ComputeBondDynamics(timestep,i,dt,Dt,...
%     p1,p2,bond_types,bond_id,...
%     types,dynamic_pairs0,...
%     N_bound0,N_free0,timestep0,tiptyp,dynbondtype,...
%     all_dynamic_pairs,N_virgin_bonds)
% 
% prev_indx = i-1;
% 
% % Define timestep
% % Dt = (timestep-timestep0)*dt;
% Time = timestep*dt;
% 
% % Define current dynamic pairs
% if ~isempty(p1{i})
%     pairs = [p1{i} p2{i} bond_types{i} bond_id{i}];
%     dynamic_pairs = pairs(pairs(:,3)==dynbondtype,:);
%     pair_col = sort(dynamic_pairs(:,1:2),2);
%     dynamic_pairs(:,1:2) = pair_col;
%     dynamic_pairs = sortrows(dynamic_pairs,1);
%     dynamic_pairs(1:2:end,:) = [];
% 
%     % compute number of bound and free stickers
%     N_stickers = length(types(types==tiptyp));
%     N_bound = size(dynamic_pairs,1);
%     N_free = N_stickers-N_bound;
% 
%     % Define previous dynamic pairs
% %     if i==1
% %         dynamic_pairs0 = [];
% %     else
%     if ~isempty(dynamic_pairs0)
%         pairs0 = [p1{prev_indx} p2{prev_indx} bond_types{prev_indx} bond_id{prev_indx}];
%         dynamic_pairs0 = pairs0(pairs0(:,3)==dynbondtype,:);
%         pair_col = sort(dynamic_pairs0(:,1:2),2);
%         dynamic_pairs0(:,1:2) = pair_col;
%         dynamic_pairs0 = sortrows(dynamic_pairs0,1);
%         dynamic_pairs0(1:2:end,:) = [];
%     end
% 
%     % ID previously existing and new dynamic bonds
%     current_pairs = dynamic_pairs(:,1:2);
%     if ~isempty(dynamic_pairs0)
%         prev_pairs = dynamic_pairs0(:,1:2);
%         new_pairs = current_pairs(~ismember(current_pairs,prev_pairs,'rows'),:);
%     else
%         prev_pairs = [];
%         new_pairs = current_pairs;
%     end
% 
%     % Define number of attachments this step & corresponding attachment rate
%     if i==1  || isempty(new_pairs)
%         ka = 0;
%         no_attachments = 0;
%     else
%         no_attachments = size(new_pairs,1);   
%         ka = no_attachments/N_free0/Dt;
%     end
% 
%     % Find first-time new pairs and define ka1
%     if ~isempty(new_pairs) && isempty(all_dynamic_pairs)
%         N_first_time_attachments = size(new_pairs,1);
%         ka1 = N_first_time_attachments/N_virgin_bonds/Dt;
%     elseif ~isempty(new_pairs) && ~isempty(all_dynamic_pairs)
%         np = sort(new_pairs,2);
%         ap = sort(all_dynamic_pairs(:,1:2),2);
% 
%         % delete pairs that already existed at some point
%         np(ismember(np,ap,'rows'),:) = [];
% 
%         N_first_time_attachments = size(np,1);
%         ka1 = N_first_time_attachments/N_virgin_bonds/Dt;
%     else
%         N_first_time_attachments = 0;
%         ka1 = 0;
%     end
%     N_virgin_bonds = N_virgin_bonds-N_first_time_attachments;
% 
%     % Store all pairs history
%     if ~isempty(new_pairs) % Add the attachment events
%         all_dynamic_pairs = [all_dynamic_pairs;[new_pairs,...
%             Time*ones(size(new_pairs,1),1),...
%             ones(size(new_pairs,1),1),...
%             i*ones(size(new_pairs,1),1)]];  % 1 in fourth position indicates onset of a new pair
%     end
% 
%     % ID lost bonds (i.e., detachments)
%     if isempty(dynamic_pairs0)
%         deleted_pairs = [];
%     else
%         surviving_nodes1 = ismember(prev_pairs(:,1),current_pairs(:,1));
%         surviving_nodes2 = ismember(prev_pairs(:,2),current_pairs(:,2));
%         indx = surviving_nodes1 + surviving_nodes2;
%         deleted_pairs = prev_pairs(indx~=2,:);
%     end
% 
%     no_detachments =  size(deleted_pairs,1);
% 
%     if i==1 || N_bound0==0
%         kd = 0;
%     else
%         kd = no_detachments/N_bound0/Dt;
%     end
% 
%     % disp(['ka = ',num2str(ka),', kd = ',num2str(kd)])
%     if ~isempty(deleted_pairs) % Add the deletion events
%         all_dynamic_pairs = [all_dynamic_pairs;[deleted_pairs,...
%             Time*ones(size(deleted_pairs,1),1),...
%             zeros(size(deleted_pairs,1),1),...
%             i*ones(size(deleted_pairs,1),1)]]; % 1 in fourth position indicates deletion of a pair
%     end
% 
%     % Set old values
%     dynamic_pairs0 = dynamic_pairs;
%     N_bound0 = N_bound;
%     N_free0 = N_free;
%     timestep0 = timestep;
% else
%     ka1 = 0;
%     ka = 0;
%     kd = 0;
% end
% 
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kd,ka,ka1,k_rpt,dynamic_pairs0,N_bound0,N_free0,timestep0,...
    all_dynamic_pairs,N_virgin_bonds] = ...
    ComputeBondDynamics(timestep,i,dt,BondsTemp,AtomsTemp,...
    dynamic_pairs0,N_bound0,N_free0,timestep0,tiptyp,dynbondtype,...
    all_dynamic_pairs,N_virgin_bonds,BeadSpringOrMeso)

Dt = (timestep-timestep0)*dt;
Time = timestep*dt;

if ~isempty(BondsTemp)
    % Isolate current data
%     BondsTemp = BondData(BondData(:,1)==timestep,:);
%     AtomsTemp = AtomData(AtomData(:,1)==timestep,:);

    % Isolate dynamic bonds
    pairs = BondsTemp(:,4:5);
    bondtypes = BondsTemp(:,3);
    dynamic_pairs = pairs(bondtypes==dynbondtype,:);

    % Delete redundant dynamic pairs
    dynamic_pairs = sort(dynamic_pairs,2);
    dynamic_pairs = sortrows(dynamic_pairs,1);
    dynamic_pairs(2:2:end,:) = [];

    % compute number of bound and free stickers
    N_stickers = length(AtomsTemp(AtomsTemp(:,end-1)==tiptyp));
    N_bound = size(dynamic_pairs,1);

    if BeadSpringOrMeso==2
        N_free = N_stickers-2*N_bound; % No free stickers
    else
        N_free = N_stickers-N_bound;
    end

    % ID previously existing and new dynamic bonds
    if isempty(dynamic_pairs0)
        %     old_pairs = [];
        new_pairs = dynamic_pairs;
    else
        old_nodes1 = ismember(dynamic_pairs(:,1),dynamic_pairs0(:,1));
        old_nodes2 = ismember(dynamic_pairs(:,2),dynamic_pairs0(:,2));
        indx = old_nodes1 + old_nodes2; %new if this sums to 2
        %     old_pairs = dynamic_pairs(indx==2,:);
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
    no_rpt_attachments = size(new_pairs,1)-no_first_time_attachments;
    k_rpt = no_rpt_attachments/N_free0/Dt;
    N_virgin_bonds = N_virgin_bonds-no_first_time_attachments;

    % Store all pairs history
    if ~isempty(new_pairs) % Add the attachment events
        % columns 1:2 are stickers 1 and 2, column 3 is the time of the
        % event and column 4 is 1 for attachment, 0 for detachment
        all_dynamic_pairs = [all_dynamic_pairs;[new_pairs,...
            Time*ones(size(new_pairs,1),1),...
            ones(size(new_pairs,1),1)]];
    end

    no_attachments = size(new_pairs,1);   

    if i==1
        ka = 0;
    else
        if BeadSpringOrMeso==2
            ka = no_attachments/(N_free0/2)/Dt; % Normalize by number of free potential bonds 
                                                % (i.e. half the number of stickers)
        else
            ka = no_attachments/N_free0/Dt;
        end
    end

    if ka<0 || N_free0<0 || Dt<0
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
        %     surviving_pairs = dynamic_pairs0(indx==2);
        deleted_pairs = dynamic_pairs0(indx~=2,:);
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
            zeros(size(deleted_pairs,1),1)]];
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
    msd00_r,msd00_r_err,msd00_t,msd00_t_err,PrevStep]...
    = ComputeMSDs(timestep,AtomData,AtomDataC,Corners,MSDStartStep,PrevStep,...
    stickertype,tethertype)

% if timestep==30 || timestep==35
%     timestep;
% end

Lx = Corners(5)-Corners(2);
Ly = Corners(6)-Corners(3);
Lz = Corners(7)-Corners(4);

% % Sort AtomData by ID number
% AtomData = sortrows(AtomData,1);

% Check that these are aligned
AtomData00 = AtomData(AtomData(:,1)==MSDStartStep,:);   %First step
AtomData0 = AtomData(AtomData(:,1)==PrevStep,:);   %Previous step
% AtomDataC = AtomData(AtomData(:,1)==timestep,:);

% Sort AtomDatas by particle ID
AtomData00 = sortrows(AtomData00,2);
AtomData0 = sortrows(AtomData0,2);
AtomDataC = sortrows(AtomDataC,2);

Pos00 = AtomData00(:,3:5);
Pos = AtomDataC(:,3:5);

% Find tether particle of corresponding molecule
AtomData_st = AtomDataC(AtomDataC(:,end-1)==tethertype,:);
AtomData_th = AtomDataC(AtomDataC(:,end-1)==stickertype,:);
AtomData_st = sortrows(AtomData_st,10);
AtomData_th = sortrows(AtomData_th,10);
R = AtomData_st(:,3:5)-AtomData_th(:,3:5);

if sum(ismember(AtomData_st(:,end),AtomData_th(:,end),'rows'))~=size(AtomData_st,1)
    [AtomData_st(:,end),AtomData_th(:,end)];
end

AtomData0_st = AtomData0(AtomData0(:,end-1)==tethertype,:);
AtomData0_th = AtomData0(AtomData0(:,end-1)==stickertype,:);
AtomData0_st = sortrows(AtomData0_st,10);
AtomData0_th = sortrows(AtomData0_th,10);
R0 = AtomData0_st(:,3:5)-AtomData0_th(:,3:5);

% if sum(ismember(AtomData0_st(:,end),AtomData0_th(:,end),'rows'))~=size(AtomData0_st,1)
%     [AtomData0_st(:,end),AtomData0_th(:,end)];
% end

AtomData00_st = AtomData00(AtomData00(:,end-1)==tethertype,:);
AtomData00_th = AtomData00(AtomData00(:,end-1)==stickertype,:);
AtomData00_st = sortrows(AtomData00_st,10);
AtomData00_th = sortrows(AtomData00_th,10);
R00 = AtomData00_st(:,3:5)-AtomData00_th(:,3:5);

if sum(ismember(AtomData00_st(:,end),AtomData00_th(:,end),'rows'))~=size(AtomData00_st,1)
    [AtomData00_st(:,end),AtomData00_th(:,end)];
end

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

% if max(theta0)>pi
%     theta0;
% end

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

% if max(theta00)>pi
%     theta00;
% end

% Eliminate partilces of wrong type
Pos00(AtomData00(:,end-1)~=stickertype,:) = [];
Pos(AtomDataC(:,end-1)~=stickertype,:) = [];

% N = size(Pos,1);
dr = Pos-Pos00;

% Maxdr = max(vecnorm(dr,2,2));
% if Maxdr>8
%     disp(Maxdr);
% end

% Adjust for periodic bounds
if ~isempty(dr(abs(dr(:,1))>0.5*Lx,1)) || ~isempty(dr(abs(dr(:,2))>0.5*Ly,2)) ||...
        ~isempty(dr(abs(dr(:,3))>0.5*Lz,3))
    dr(abs(dr(:,1))>0.5*Lx,1);
end
dr(abs(dr(:,1))>0.5*Lx,1) = dr(abs(dr(:,1))>0.5*Lx,1)-Lx*sign(dr(abs(dr(:,1))>0.5*Lx,1));
dr(abs(dr(:,2))>0.5*Ly,2) = dr(abs(dr(:,2))>0.5*Ly,2)-Ly*sign(dr(abs(dr(:,2))>0.5*Ly,2));
dr(abs(dr(:,3))>0.5*Lz,3) = dr(abs(dr(:,3))>0.5*Lz,3)-Lz*sign(dr(abs(dr(:,3))>0.5*Lz,3));

% for type=1:2
%     dr_temp = dr(AtomData0(:,end)==type,:);
MSD_temp = mean((vecnorm(dr,2,2)).^2);
MSD_err = std((vecnorm(dr,2,2)).^2)/sqrt(size(dr,1));

msd = MSD_temp;
msd_err = MSD_err;
% end

PrevStep = timestep;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigma,conc] = ComputeStressAndBondConc(f,r,Corners)

global LengthConversion

% Compute volume
LWH = LengthConversion*abs(Corners(2:4)-Corners(5:7));
r = r*LengthConversion;
V = prod(LWH);

% Dyad f \otimes r
dyad = [f(:,1).*r(:,1),...
    f(:,2).*r(:,2),...
    f(:,3).*r(:,3),...
    f(:,1).*r(:,2),...
    f(:,2).*r(:,3),...
    f(:,3).*r(:,1)];

N_bonds = size(f,1);
conc = N_bonds/V;

sigma = 1/(2*V)*sum(dyad,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function forces = ComputeForces(types,r,N_Kuhn,b,kbT)

global LengthConversion BeadSpringOrMeso

types = types';
b = b*LengthConversion;
r = r*LengthConversion;
norms = vecnorm(r,2,2);
r_hat = [r(:,1)./norms r(:,2)./norms r(:,3)./norms];

force_mags = zeros(size(types));
Stiffness1 = 3*kbT/(N_Kuhn*b^2);
Stiffness3 = 3*kbT/(b^2)/2; % Divide by 2 because dynamic bonds are double counted in TNT bond package

if BeadSpringOrMeso==0
    force_mags(types==4) = Stiffness1*norms(types==4);
    force_mags(types==2) = Stiffness1*norms(types==2);
else
    force_mags(types==1) = Stiffness1*norms(types==1);
    force_mags(types==2) = Stiffness1*norms(types==2);
    force_mags(types==3) = Stiffness3*norms(types==3);
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
function rs = ExtractEndToEnd(AtomDataTemp,BondDataTemp,Np,N_Kuhn,i,samp)

global BeadSpringOrMeso bondtypes StationaryStableNodes...
    StationaryDynamicNodes DiffusingDynamicNodes DiffusingStableNodes 

atomtypes = AtomDataTemp(:,end-1);
atomnumbers = AtomDataTemp(:,2);
% atompositions = AtomDataTemp(:,3:5);
% molnumber = AtomDataTemp(:,end);

if BeadSpringOrMeso==1
    StationaryStableNodes = atomnumbers(atomtypes==1);
    StationaryDynamicNodes = atomnumbers(atomtypes==2);
    DiffusingDynamicNodes = atomnumbers(atomtypes==3);
    
    BondDataTemp(BondDataTemp(:,3)~=2,:) = []; %Only concerned with end-to-end of chains
    rs = BondDataTemp(:,6:8);
    if size(BondDataTemp,1)~=Np
        rs = BondDataTemp(:,6:8);
    end
    bondtypes(i,:,samp) = BondDataTemp(:,3);

    
    %             rs = BondData(BondData(:,1)==timestep,6:8);
    %             bondtypes(i,:,samp) = BondData(BondData(:,1)==timestep,3);
else
    %             atomtypes = AtomData(AtomData(:,1)==timestep,end-1);
    %             atomnumbers = AtomData(AtomData(:,1)==timestep,2);
    %             atompositions = AtomData(AtomData(:,1)==timestep,3:5);

    StationaryStableNodes = atomnumbers(atomtypes==1);
    StationaryDynamicNodes = atomnumbers(atomtypes==2);
    DiffusingDynamicNodes = atomnumbers(atomtypes==3);
    DiffusingStableNodes = atomnumbers(atomtypes==4);

%     NChains = length(StationaryStableNodes);
%     rs = zeros(NChains,3);

    % Isolate 
    StationaryStableAtomData = ...
        AtomDataTemp(ismember(atomnumbers,StationaryStableNodes),:);
    DiffusingDynamicAtomData = ...
        AtomDataTemp(ismember(atomnumbers,DiffusingDynamicNodes),:);

    % If these matrices are different sizes, call error
    if size(StationaryStableAtomData,1)~=size(DiffusingDynamicAtomData,1)
        error('Different number of stable and dynamic chain ends')
    end

    % Sort by molecule number
    tethers = sortrows(StationaryStableAtomData,10);
    stickers = sortrows(DiffusingDynamicAtomData,10);

    % Define end-to-end vectors
    rs = stickers(:,3:5) - tethers(:,3:5);
    
%     for chainnumber=1:NChains
%     for i=1:length(StationaryStableNodes)
%         chainnumber = StationaryStableNodes(i);
%         tether = chainnumber;
%         sticker = NChains*2 + N_Kuhn*chainnumber; %Goes past chain number 
% 
%         type_t = atomtypes(atomnumbers==tether,1);
%         type_s = atomtypes(atomnumbers==sticker,1);
%         if type_t~=1 || type_s~=4
%             error('Check atom types')
%         end
% 
%         pos_t = atompositions(atomnumbers==tether,:);
%         pos_s = atompositions(atomnumbers==sticker,:);
%         rs(chainnumber,:) = pos_s-pos_t;
%     end
    bondtypes(i,:,samp) = 4;%BondDataTemp(:,3);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PreallocateOutputs(N_steps,N_samples,Np,N_Kuhn)

global BeadSpringOrMeso...
    sig11 sig22 sig33 sig12 sig23 sig31...
    sig11_err sig22_err sig33_err sig12_err sig23_err sig31_err...
    bondconc time stretch...
    msd_st msd_st_err...
    msd_st_r0 msd_st_r0_err...
    msd_st_r00 msd_st_r00_err...
    msd_st_t0 msd_st_t0_err...
    msd_st_t00 msd_st_t00_err...
    ka_out ka1_out kd_out k_rpt Na Nd...
    rx ry rz bondtypes...
    rxr11_st rxr22_st rxr33_st rxr12_st rxr23_st rxr31_st...
    rxr11_st_err rxr22_st_err rxr33_st_err rxr12_st_err rxr23_st_err rxr31_st_err

% Stress
sig11 = zeros(N_steps,N_samples);
sig22 = zeros(N_steps,N_samples);
sig33 = zeros(N_steps,N_samples);
sig12 = zeros(N_steps,N_samples);
sig23 = zeros(N_steps,N_samples);
sig31 = zeros(N_steps,N_samples);
sig11_err = zeros(N_steps,N_samples);
sig22_err = zeros(N_steps,N_samples);
sig33_err = zeros(N_steps,N_samples);
sig12_err = zeros(N_steps,N_samples);
sig23_err = zeros(N_steps,N_samples);
sig31_err = zeros(N_steps,N_samples);
time = zeros(N_steps,1);
stretch = zeros(N_steps,1);
bondconc = zeros(N_steps,N_samples);

% MSD
msd_st = zeros(N_steps,N_samples);
msd_st_err = zeros(N_steps,N_samples);
msd_st_r0 = zeros(N_steps,N_samples);
msd_st_r0_err = zeros(N_steps,N_samples);
msd_st_r00 = zeros(N_steps,N_samples);
msd_st_r00_err = zeros(N_steps,N_samples);
msd_st_t0 = zeros(N_steps,N_samples);
msd_st_t0_err = zeros(N_steps,N_samples);
msd_st_t00 = zeros(N_steps,N_samples);
msd_st_t00_err = zeros(N_steps,N_samples);

% Kinetics
ka_out = zeros(N_steps,N_samples);
ka1_out = zeros(N_steps,N_samples);
k_rpt = zeros(N_steps,N_samples);
kd_out = zeros(N_steps,N_samples);
Na = zeros(N_steps,N_samples);
Nd = zeros(N_steps,N_samples);

% End-to-end vector data
% if BeadSpringOrMeso==0
%     N_chains = N_seg/N_Kuhn;    %rx, ry, rz, for all N_seg segments
%                             %corresponding standdar erros of the mean
% else
%     N_chains = N_seg;
% end

N_chains = Np;        %Np is the number of pairs 

% Note - [rx,ry,rz] is resereved for the end-to-end vectors of a full chain
% occuring betweetween two crosslink sites, not individual Kuhn segments
rx = zeros(N_steps,N_chains,N_samples);
ry = zeros(N_steps,N_chains,N_samples);
rz = zeros(N_steps,N_chains,N_samples);
bondtypes = zeros(N_steps,N_chains,N_samples);

rxr11_st = zeros(N_steps,N_samples);
rxr22_st = zeros(N_steps,N_samples);
rxr33_st = zeros(N_steps,N_samples);
rxr12_st = zeros(N_steps,N_samples);
rxr23_st = zeros(N_steps,N_samples);
rxr31_st = zeros(N_steps,N_samples);
rxr11_st_err = zeros(N_steps,N_samples);
rxr22_st_err = zeros(N_steps,N_samples);
rxr33_st_err = zeros(N_steps,N_samples);
rxr12_st_err = zeros(N_steps,N_samples);
rxr23_st_err = zeros(N_steps,N_samples);
rxr31_st_err = zeros(N_steps,N_samples);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compute = CheckIfNeedToComputeData(OverrideCompute)

global StressDataFileName MSDDataFileName EndToEndDataFileName...
    BondKineticsDataFileName AlignmentDataFileName...
    TimeStretchDataFileName...
    CalculateStress CalculateEndtoEnd...
    CalculateMSD CalculateBondKinetics CalculateAlignment

compute = 0;
if (~isfile(StressDataFileName) && CalculateStress) ||...
        (~isfile(MSDDataFileName) && CalculateEndtoEnd) ||...
        (~isfile(EndToEndDataFileName) && CalculateMSD)  ||...
        (~isfile(BondKineticsDataFileName) && CalculateBondKinetics)  ||...
        (~isfile(AlignmentDataFileName) && CalculateAlignment)  ||...
        ~isfile(TimeStretchDataFileName)  ||...
        OverrideCompute
    compute = 1;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotIntStressStrain(steps,time,Np,...
    sig11_ensemble, sig22_ensemble,sig33_ensemble)

global FileTag FontSize TurnOnDynamics

[sig11,sig11_err] = CalcMeanAndSE(sig11_ensemble,2);
[sig22,sig22_err] = CalcMeanAndSE(sig22_ensemble,2);
[sig33,sig33_err] = CalcMeanAndSE(sig33_ensemble,2);

figure(2); clf; hold on
TempStress = sig11; TempErr = sig11_err;
Color = 'k'; Style = '-'; TempStep = 1;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);
Color = 'b'; Style = '-'; TempStep = 2;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);
Color = 'c'; Style = '-'; TempStep = 3;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);

TempStress = sig22; TempErr = sig22_err;
Color = 'k'; Style = '--'; TempStep = 1;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);
Color = 'b'; TempStep = 2;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);
Color = 'c'; TempStep = 3;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);

TempStress = sig33; TempErr = sig33_err;
Color = 'k'; Style = ':'; TempStep = 1;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);
Color = 'b'; TempStep = 2;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);
Color = 'c'; TempStep = 3;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);

set(gca,'FontSize',FontSize/1.5)
xlabel('$t$ [s]','FontSize',FontSize,'Interpreter','latex')
ylabel('$\sigma$ [kPa]','FontSize',FontSize,'Interpreter','latex')
pbaspect([1 1 1])
set(gcf,'Color','w')

if ~isfolder('Output Plots')
    mkdir('Output Plots')
end
if TurnOnDynamics==0
    FileName = 'Control.Stress Check';
    title(['Dynamics OFF, $N_p$ = ',num2str(Np)],'FontSize',FontSize/1.5,'Interpreter','latex')
else
    FileName = 'Stress Check';
    title(['Dynamics ON, $N_p$ = ',num2str(Np)],'FontSize',FontSize/1.5,'Interpreter','latex')
end
CheckStressFileName = ['Output Plots/',FileName,FileTag,'.png'];
saveas(gcf,CheckStressFileName)
CheckStressFileName = ['Output Plots/',FileName,FileTag,'.fig'];
saveas(gcf,CheckStressFileName)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveCompiledFiles(sig11_ensemble,sig22_ensemble,sig33_ensemble,...
    time,stretch,steps)

global CompiledStress11FileName CompiledStress22FileName...
    CompiledStress33FileName

writematrix([sig11_ensemble,time,stretch,steps],CompiledStress11FileName)
writematrix([sig22_ensemble,time,stretch,steps],CompiledStress22FileName)
writematrix([sig33_ensemble,time,stretch,steps],CompiledStress33FileName)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveSampleFiles(Corners,Atoms,Bonds)

global RawDataFileName %AtomsFileName BondsFileName CornersFileName 

% structArray = cell2struct(cellArray, fields, dim)

Raw.Corners = Corners;
Raw.Atoms = Atoms;
Raw.Bonds = Bonds;
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
function [timestep,Corners,Atoms] = UnpackAtomsDump

global OutputAtom

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

% Extract periodic domain boundaries
fmt1=repmat('%f',1,6);             % first data section format
fid=fopen(OutputAtom,'r');     % open file
l=fgetl(fid);                      % read first line
ct = 0;
wb = waitbar(0,'Extracting domain boundaries...');
while ~feof(fid)                   % loop through file by record
    if strfind(l,'ITEM: BOX BOUNDS pp pp pp')
        ct = ct+1;
        Corners_array{:,:,ct} = cell2mat(textscan(fid,'%f'));
        if ~mod(ct,50)
            PcntComplete = ct/N_output;
            waitbar(PcntComplete,wb,'Extracting domain boundaries...')
        end
    end                                 % break when find first section
    l=fgetl(fid);                       % next record
end
close(wb)

% Unpack corners data
N_var = 6;      % x1 x2 y1 y2 z1 z2
Corners = zeros(N_output,N_var+1);  %+1 for timesteps
Corners(:,1) = timestep;
wb = waitbar(0,'Unpacking boundaries...');
for i=1:N_output
    if ~mod(i,50)
        PcntComplete = i/N_output;
        waitbar(PcntComplete,wb,'Unpacking boundaries...')
    end
    Corners_temp = Corners_array{:,:,i};
    Corners(i,2) = Corners_temp(1);
    Corners(i,3) = Corners_temp(3);
    Corners(i,4) = Corners_temp(5);
    Corners(i,5) = Corners_temp(2);
    Corners(i,6) = Corners_temp(4);
    Corners(i,7) = Corners_temp(6);
end
close(wb)

% Extract atom index, position, velocity, and type
% information
fmt1=repmat('%f',1,6);             % first data section format
fid=fopen(OutputAtom,'r');     % open file
l=fgetl(fid);                      % read first line
ct = 0;
wb = waitbar(0,'Extracting atom data...');
while ~feof(fid)                   % loop through file by record
    if strfind(l,'ITEM: ATOMS id x y z vx vy vz type')
        ct = ct+1;
        if ~mod(ct,50)
            PcntComplete = ct/N_output;
            waitbar(PcntComplete,wb,'Extracting atom data...')
        end
        Atoms_array{:,:,ct} = cell2mat(textscan(fid,'%f'));
    end                                 % break when find first section
    l=fgetl(fid);                       % next record
end
close(wb)

% Unpack atom information
N_var = 9;      % id x y z vx vy vz type mol
Atoms = zeros(1,N_var+1);
ct = 0;
wb = waitbar(0,'Unpacking atom data...');
StartIndx = 0;
for i=1:N_output
    if ~mod(i,50)
        PcntComplete = i/N_output;
        waitbar(PcntComplete,wb,'Unpacking atom data...')
    end
    timestep_temp = timestep(i);
    Atoms_temp = Atoms_array{:,:,i};
    N_atoms = length(Atoms_temp)/N_var;
    Range = StartIndx+(1:N_atoms)';
    
    ID_indx = (1:N_var:size(Atoms_temp,1))';
    x_indx = ID_indx+1;
    y_indx = x_indx+1;
    z_indx = y_indx+1;
    vx_indx = z_indx+1;
    vy_indx = vx_indx+1;
    vz_indx = vy_indx+1;
    type_indx = vz_indx+1;
    mol_indx = type_indx+1;
   
    Atoms(Range,1) = timestep_temp*ones(size(Range));
    Atoms(Range,2) = Atoms_temp(ID_indx);
    Atoms(Range,3) = Atoms_temp(x_indx);
    Atoms(Range,4) = Atoms_temp(y_indx);
    Atoms(Range,5) = Atoms_temp(z_indx);
    Atoms(Range,6) = Atoms_temp(vx_indx);
    Atoms(Range,7) = Atoms_temp(vy_indx);
    Atoms(Range,8) = Atoms_temp(vz_indx);
    Atoms(Range,9) = Atoms_temp(type_indx);
    Atoms(Range,10) = Atoms_temp(mol_indx);

    StartIndx = size(Atoms,1);
end
close(wb)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function time = ExtractTime(EqSize)

global StressEqlFileName

% Import data
% During equilibration
EquilibrationData = dlmread(['Data/Outputs/',StressEqlFileName],'',1,0);
if size(EquilibrationData,1)~=EqSize
    Diff = size(EquilibrationData,1)-EqSize;
    if Diff>0
        EquilibrationData(end-Diff+1:end,:) = [];
    else
        EquilibrationData(end:end+Diff,:) = NaN;
    end
end
time = EquilibrationData(:,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EqSize = DeterminAllocationSize(PackageTemp,...
    EdgeFactor,Np,N,ka,kd,f0,dt,damp,N_Kuhn,b,Separation)

global StressEqlFileName  

% Initializes non-sweeping parameters
InputScript(EdgeFactor,1,Np,N,ka,kd,f0,dt,damp,...
    N_Kuhn,b,Separation);
SetDirAndFileNames;

%Compute shortest data sizes
EqSizes = zeros(size(PackageTemp,1),1);
% LdSizes = zeros(size(PackageTemp,1),1);
% RlSizes = zeros(size(PackageTemp,1),1);

for n=1:size(PackageTemp,1)

    %% Unpack Swept Input Parameters
    Sample = PackageTemp(n,1);    %Sample index

    %% Callout Input Script
    % Initializes non-sweeping parameters
    InputScript(EdgeFactor,Sample,Np,N,ka,kd,f0,dt,damp,...
        N_Kuhn,b,Separation);

    %% Make diriectories and set filenames
    SetDirAndFileNames;

    %% Import data

    % During equilibration
    EquilibrationData = dlmread(['Data/Outputs/',StressEqlFileName],'',1,0);
    EqSizes(n) = size(EquilibrationData,1);

%     % During deformation
%     LoadingData = dlmread(['Outputs/',StressLdgFileName],'',1,0);
%     LdSizes(n) = size(LoadingData,1);
% 
%     % During relaxation
%     RelaxationData = dlmread(['Outputs/',StressRlxFileName],'',1,0);
%     RlSizes(n) =  size(RelaxationData,1);
end

EqSize = min(EqSizes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,f] = PlotCurve(X,Y,Err,Color,Style)

global LineWidth

ErrX = [X;flipud(X)];
ErrY = [Y+Err;flipud(Y-Err)];
f = fill(ErrX,ErrY,Color);
f.FaceAlpha = 0.25;
f.LineStyle = 'none';

p = plot(X,Y);
p.Color = Color;
p.LineWidth = LineWidth;
p.LineStyle = Style;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mean,SE] = CalcMeanAndSE(Data,Dim)

Mean = nanmean(Data,Dim);
SE = nanstd(Data,0,Dim)/sqrt(size(Data,Dim));

end
