function CompileData(Package,OverrideCompile,OverrideCompute,ToggleDynamics,...
    LC,DC,BSOM,CS,CE,CM,CB,CA,CD,OD,NoProcessors,LoL)

% Collects data from all LAMMPS runs and compiles it for ensemble averaging
% or collective histogram outputs

% Compiles stress-time responses
% Stress-strain responses

% End-to-end data for histogram binning
% Bond dynamics data

PackageTypes = DefineUniquePackageTypes(Package);

if NoProcessors==1
    LoopCompileData(Package,OverrideCompile,OverrideCompute,ToggleDynamics,...
        LC,DC,BSOM,CS,CE,CM,CB,CA,CD,OD,NoProcessors,LoL)
else
    delete(gcp('nocreate'))
    parpool(NoProcessors)
    parfor n=1:size(PackageTypes,1)
        LoopCompileData(PackageTypes(n,:),OverrideCompile,OverrideCompute,...
            ToggleDynamics,LC,DC,BSOM,CS,CE,CM,CB,CA,CD,OD,NoProcessors,LoL)
    end
    delete(gcp('nocreate'))
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LoopCompileData(Package,OverrideCompile,OverrideCompute,ToggleDynamics,...
    LC,DC,BSOM,CS,CE,CM,CB,CA,CF,OF,NoProcessors,LoL)

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
    CurrentFolder OutputFolder dtFact...
    OutputAtom OutputBond OutputDataFolder LinearOrLangevin

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
LinearOrLangevin = LoL;

if BeadSpringOrMeso==2
    CalculateStress = 0;            % 1 to compute virial stress-stretch/time data
    CalculateEndtoEnd = 0;          % 1 to compile end-to-end distribution data and
    CalculateMSD = 0;               % 1 to compute the mean-square displacement (i.e. diffusion) data
    CalculateAlignment = 0;         % 1 to compute metric tensors (r \otimes r) of chains that elucidates alignment
end

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
kas = unique(Package(:,4));         %Temperature
kds = unique(Package(:,5));        %Contour length of chains
f0s = unique(Package(:,6));        %Activation energy of association
N_Kuhns = unique(Package(:,7));        %Activation energy of dissociation
bs = unique(Package(:,8));
Separations = unique(Package(:,9));
phis = unique(Package(:,10));
kbTs = unique(Package(:,11));

for j=1:length(Nps)
 for k=1:length(Nss)
  for l=1:length(kas)
   for m=1:length(kds)
    for n=1:length(f0s)
     for o=1:length(N_Kuhns)
      for p=1:length(bs)
       for q=1:length(Separations)
        for r=1:length(phis)
         for s=1:length(kbTs)
            Np = Nps(j);
            Ns = Nss(k);
            N = Np*Ns;             
            ka = kas(l);
            kd = kds(m);
            f0 = f0s(n);
            N_Kuhn = N_Kuhns(o);
            b = bs(p);
            Separation = Separations(q);
            phi = phis(r);
            kbT = kbTs(s);
                       
            [~,D,tau0,dtFact] = DefineTimeScale(b,LengthConversion,DamperConversion,...
                 BeadSpringOrMeso);
             
            dt = tau0/dtFact;
             
            eaStar = -log(ka*(b*LengthConversion)^2/D);
            edStar = -log(kd*(b*LengthConversion)^2/D);


            PackageTemp = Package(ismember(Package(:,2),Np),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,3),Ns),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,4),ka),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,5),kd),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,6),f0),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,7),N_Kuhn),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,8),b),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,9),Separation),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,10),phi),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,11),kbT),:);

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
                InputScript(Sample,Np,Ns,N,N_Kuhn,b,phi,eaStar,edStar,D,dt,kbT);

                %% Make diriectories and set filenames
                SetDirAndFileNames;
                DefineCompiledFileNames;
                SizeTheDomain(Np,Separation);


                if ~isfile(RawDataFileName) || OverrideCompile
                    %% Extract data from atoms.dump and bonds.dump files
                    tic
                    OutputAtom_loc = [OutputDataFolder,OutputAtom];
                    OutputBond_loc = [OutputDataFolder,OutputBond];
                    ParseDumpFiles(OutputAtom_loc,OutputBond_loc,dt,dtFact);

                    %% Save File for Sample
%                     save(RawDataFileName,'-struct','RawData','-v7.3')
                    toc
%                     delete(OutputAtom); delete(OutputBond);
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
                        InputScript(Sample,Np,Ns,N,N_Kuhn,b,phi,eaStar,edStar,D,dt,kbT);

                        %% Make diriectories and set filenames
                        SetDirAndFileNames;
                        DefineCompiledFileNames;

                        %% Allocate sizes of outputs
                        if samp==1
                            dat = load(RawDataFileName,'-mat','Corners');
                            Corners = dat.Corners;
                            N_steps = size(Corners,1);
                            N_samples = size(PackageTemp,1);
                            PreallocateOutputs(N_steps,N_samples,Np);
                        end

                        %% Calculate outputs
                        tic
                        ComputeOutputs(Np,N_Kuhn,b,kbT,dt,samp,OverrideCompute);
                        toc
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
if NoProcessors==1
    close(wb1)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SizeTheDomain(Np,Separation)

global Lx Ly Lz PairsPerEdge c0

PairsPerEdge = round(Np^(1/3));

Lx = PairsPerEdge*Separation;
Ly = Lx;
Lz = Lx;

c0 = Np/(Lx*Ly*Lz);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PackageTypes = DefineUniquePackageTypes(Package)

PackageTypes = zeros(1e3,size(Package,2));

Nps = unique(Package(:,2));        %Number of molecules
Nss = unique(Package(:,3));        %Number of stickeres per tether site
kas = unique(Package(:,4));        %Activation energy of association
kds = unique(Package(:,5));        %Activation energy of dissociation
f0s = unique(Package(:,6));       %Force sensitivity to dissociation
N_Kuhns = unique(Package(:,7));
bs = unique(Package(:,8));
Separations = unique(Package(:,9));
phis = unique(Package(:,10));
kbTs = unique(Package(:,11));

ct = 0;
for j=1:length(Nps)
 for k=1:length(Nss)
  for l=1:length(kas)
   for m=1:length(kds)
    for n=1:length(f0s)
     for o=1:length(N_Kuhns)
      for p=1:length(bs)
       for q=1:length(Separations)
        for r=1:length(phis)
         for s=1:length(kbTs)

            Np = Nps(j);
            Ns = Nss(k);
            ka = kas(l);
            kd = kds(m);
            f0 = f0s(n);
            N_Kuhn = N_Kuhns(o);
            b = bs(p);
            Separation = Separations(q);
            phi = phis(r);
            kbT = kbTs(s);

            PackageTemp = Package(ismember(Package(:,2),Np),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,3),Ns),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,4),ka),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,5),kd),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,6),f0),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,7),N_Kuhn),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,8),b),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,9),Separation),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,10),phi),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,11),kbT),:);

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
    ka_out ka1_out kd_out Na Nd...
    rx ry rz bondtypes...
    rxr11_st rxr22_st rxr33_st rxr12_st rxr23_st rxr31_st...
    rxr11_st_err rxr22_st_err rxr33_st_err rxr12_st_err rxr23_st_err rxr31_st_err...
    attached_bond_lifetimes detached_bond_lifetimes renormalized_bond_lifetimes...
    open_repeat_lifetimes open_exchange_lifetimes...
    total_exchange_events...
    total_repeat_events...
    chain_concentration adjusted_bond_lifetimes % mean_exchange_rate se_exchange_rate

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
    BondKinetics.kd = kd_out;
    BondKinetics.Na = Na;
    BondKinetics.Nd = Nd;
    BondKinetics.AttachedLifetimes = attached_bond_lifetimes;
    BondKinetics.DetachedLifetimes = detached_bond_lifetimes;
    BondKinetics.RenormalizedLifetimes = renormalized_bond_lifetimes;
    BondKinetics.AdjustedBondLifetimes = adjusted_bond_lifetimes;
    BondKinetics.OpenRepeatLifetimes = open_repeat_lifetimes;
    BondKinetics.OpenExchangeLifetimes = open_exchange_lifetimes;
    BondKinetics.ExchangeEvents = total_exchange_events;
    BondKinetics.RepeatEvents = total_repeat_events;
    BondKinetics.Concentration = chain_concentration;
    
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
function ComputeOutputs(Np,N_Kuhn,b,kbT,dt,samp,...
    OverrideCompute)

global RawDataFileName CalculateStress CalculateEndtoEnd...
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
    ka_out ka1_out kd_out Na Nd...
    rx ry rz bondtypes...
    rxr11_st rxr22_st rxr33_st rxr12_st rxr23_st rxr31_st...
    rxr11_st_err rxr22_st_err rxr33_st_err rxr12_st_err rxr23_st_err rxr31_st_err...
    attached_bond_lifetimes detached_bond_lifetimes...
    open_repeat_lifetimes open_exchange_lifetimes...
    detached_lifetimes_to_first...
    renormalized_bond_lifetimes adjusted_bond_lifetimes...
    total_exchange_events mean_exchange_time se_exchange_time...
    total_repeat_events mean_repeat_time se_repeat_time...
    exchange_rate mean_exchange_rate se_exchange_rate...
    Lx Ly Lz LengthConversion chain_concentration

% Domain data
dat = load(RawDataFileName,'-mat','timesteps','Corners');
timesteps = dat.timesteps;
Corners = dat.Corners;

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
tic
x = [];
if (CalculateEndtoEnd && (~isfile(EndToEndDataFileName) || OverrideCompute))...
        || CalculateAlignment && (~isfile(AlignmentDataFileName) || OverrideCompute)...
        || CalculateStress && (~isfile(StressDataFileName) || OverrideCompute)
    dat = load(RawDataFileName,'-mat','type','id','x','y','z',...
        'btype','rx','ry','rz');
    types = dat.type;
    id = dat.id;
    x = dat.x; y = dat.y; z = dat.z;
    rx_all = dat.rx; ry_all = dat.ry; rz_all = dat.rz;
    bond_types = dat.btype;
end
if CalculateMSD && (~isfile(MSDDataFileName) || OverrideCompute)
    dat = load(RawDataFileName,'-mat','mol');
    mol = dat.mol;
    if isempty(x)
        dat = load(RawDataFileName,'-mat','x','y','z','type','id');
        x = dat.x;
        y = dat.y;
        z = dat.z;
        types = dat.type;
        id = dat.id;
    end
end
if CalculateBondKinetics && (~isfile(BondKineticsDataFileName) || OverrideCompute)
    dat = load(RawDataFileName,'-mat','bid','batom1','batom2','btype','type','id');
    p1 = dat.batom1;
    p2 = dat.batom2;
    bond_id = dat.bid;
    bond_types = dat.btype;
    id = dat.id;
    types = dat.type;
    attachment_hist = zeros(Np,length(timesteps)+1);
    if BeadSpringOrMeso==0
        tiptyp = 3;
        dynbondtype = 2;
    elseif BeadSpringOrMeso==1
        tiptyp = 2;
        dynbondtype = 2;
    elseif BeadSpringOrMeso==2
        tiptyp = 1;
        dynbondtype = 1;
    end
    all_p = [id{1},types{1}];
    dyn_p = all_p(all_p(:,2)==tiptyp,1);
    attachment_hist(:,1) = sort(dyn_p);
end
toc
tic
clear dat
toc

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
    rs = []; %AtomDataTemp = []; BondDataTemp = [];


    %% Store end-to-end components & Bond types
    if CalculateEndtoEnd && (~isfile(EndToEndDataFileName) || OverrideCompute)
        rs = ExtractEndToEnd(types{i},id{i},x{i},y{i},z{i},...
            bond_types{i},rx_all{i},ry_all{i},rz_all{i},...
            Np,N_Kuhn,i,samp);
        rx(i,:,samp) = rs(:,1); 
        ry(i,:,samp) = rs(:,2); 
        rz(i,:,samp) = rs(:,3);
    end
    
    %% Compute metric tensors by bond type
    if CalculateAlignment && (~isfile(AlignmentDataFileName) || OverrideCompute)
        if isempty(rs)
            rs = ExtractEndToEnd(types{i},id{i},x{i},y{i},z{i},...
                bond_types{i},rx_all{i},ry_all{i},rz_all{i},...
                Np,N_Kuhn,i,samp);
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
    if CalculateStress && (~isfile(StressDataFileName) || OverrideCompute)
        if isempty(rs)          
            rs = ExtractEndToEnd(types{i},id{i},x{i},y{i},z{i},...
                bond_types{i},rx_all{i},ry_all{i},rz_all{i},...
                Np,N_Kuhn,i,samp);
            rx(i,:,samp) = rs(:,1);
            ry(i,:,samp) = rs(:,2);
            rz(i,:,samp) = rs(:,3);
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
    if CalculateMSD && (~isfile(MSDDataFileName) || OverrideCompute)
        if timestep<=MSDStartStep
            msd_st(i) = 0;
            msd_st_err(i) = 0;
%             PrevStep = MSDStartStep;
        else
            if BeadSpringOrMeso==0
                stickertype = 3;   %End groups (MSD should follow Rouse diffusion)
                tethertype = 1;
            else
                stickertype = 2;
                tethertype = 1;
            end
            prev_indx = i-1;
            [msd_st(i,samp),msd_st_err(i,samp),...
                msd_st_r0(i,samp),msd_st_r0_err(i,samp),...
                msd_st_t0(i,samp),msd_st_t0_err(i,samp),...
                msd_st_r00(i,samp),msd_st_r00_err(i,samp),...
                msd_st_t00(i,samp),msd_st_t00_err(i,samp)] = ...
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
    if CalculateBondKinetics && (~isfile(BondKineticsDataFileName) || OverrideCompute)
%         [kd_out(i,samp),ka_out(i,samp),ka1_out(i,samp),...
%             dynamic_pairs0,N_bound0,N_free0,timestep0,all_dynamic_pairs,...
%             N_virgin] = ...
%             ComputeBondDynamics(timestep,i,dt,...
%             p1,p2,bond_types,bond_id,...
%             types{i},dynamic_pairs0,...
%             N_bound0,N_free0,...
%             timestep0,tiptyp,dynbondtype,all_dynamic_pairs,N_virgin);
%         Na(i,samp) = N_bound0;
%         Nd(i,samp) = N_free0;
        
        % Alt method for computing bond lifetime outputs
        all_pairs = [p1{i}, p2{i}, bond_types{i}];
        dyn_pairs = all_pairs(all_pairs(:,end)==dynbondtype,1:2);
        for j=1:size(dyn_pairs,1)
            attachment_hist(attachment_hist(:,1)==dyn_pairs(j,1),i+1) = dyn_pairs(j,2);
            attachment_hist(attachment_hist(:,1)==dyn_pairs(j,2),i+1) = dyn_pairs(j,1);
        end
    end
end
close(wb2)

clear x y z mol types id rx_all ry_all rz_all p1 p2 bond_types bond_id

Delta_t = time(2)-time(1);
% Refer to Stampede2/CompileDataStamped_v1.m for old method
if CalculateBondKinetics && (~isfile(BondKineticsDataFileName) || OverrideCompute)
    attachment_hist = attachment_hist';
    att_mat = attachment_hist(2:end,:); att_mat(att_mat~=0) = 1;
    Na = sum(att_mat,2);
    Nd = Np-Na;
    
    % Initialize counts of events
    first_attachments = zeros(size(att_mat));
    attachments = zeros(size(att_mat));
    detachments = zeros(size(att_mat));
    repeats = zeros(size(att_mat));
    exchanges = zeros(size(att_mat));
    
    % Initialize sample sets of different lifetime type populations
    detached_lifetimes_to_first = [];
    detached_bond_lifetimes = [];
    attached_bond_lifetimes = [];
    open_exchange_lifetimes = [];
    open_repeat_lifetimes = [];
    for j=1:Np
        history = attachment_hist(2:end,j);
        changes = diff(history);        % demarks when sticker changed state.
        % First instance always marks when
        % the sticker first bonded.
        transitions = find(changes~=0);
        if ~isempty(transitions)
            first_attachments(transitions(1),j) = 1;	% marks the first time
            % the molecule became
            % attached
            time_between = diff([0;transitions])*Delta_t;    % time between adjacent
            % events
            detached_lifetimes_to_first = cat(1,detached_lifetimes_to_first,...
                time_between(1));
            for k=2:length(transitions)
                temp_time = time_between(k);        % time from prior event to now
                state1 = history(transitions(k));	% particle to which is
                % bonded or 0 if open
                state2 = history(transitions(k)+1);	% particle to which is
                % bonded or 0 if open
                
                if state1==0 % then it is an attachment event
                    attachments(transitions(k),j) = 1;  % marks new attachment event
                    % If it is an attachment event, must also check who the
                    % prior neighbor was to see if repeat or exchange
                    prior_neighbor = history(transitions(k-1)-1);
                    
                    % add a detached lifetime to the sample set
                    detached_bond_lifetimes = cat(1,detached_bond_lifetimes,temp_time);
                    if prior_neighbor==state2
                        repeats(transitions(k),j) = 1;
                        
                        % add a detached lifetime prior to repeat attachment to sample set
                        open_repeat_lifetimes = cat(1,open_repeat_lifetimes,...
                            temp_time);
                    else
                        exchanges(transitions(k),j) = 1;
                        
                        % add a detached lifetime prior to exchange attachment to sample set
                        open_exchange_lifetimes = cat(1,open_exchange_lifetimes,...
                            temp_time);
                    end
                elseif state2==0 % then it is a detachment event
                    detachments(transitions(k),j) = 1;  % marks new detachment event
                    
                    % add an attached lifetime to the sample set
                    attached_bond_lifetimes = cat(1,attached_bond_lifetimes,temp_time);
                    
                else % if both are zero, must have been an exchange event
                    % => both attachment and detachment occured
                    attachments(transitions(k),j) = 1;  % marks new attachment event
                    detachments(transitions(k)-1,j) = 1;  % marks new detachment event
                    exchanges(transitions(k),j) = 1;	% marks exchange event
                    
                    % add a detached lifetime to sample set
                    detached_bond_lifetimes = cat(1,detached_bond_lifetimes,temp_time);
                    
                    % add a detached lifetime prior to exchange attachment to sample set
                    open_exchange_lifetimes = cat(1,open_exchange_lifetimes,temp_time);
                end
            end
        end
    end
    
    kd_out = sum(detachments,2)./Na/(dt*diff(timesteps(1:2)));
    ka_out = sum(attachments,2)./Nd/(dt*diff(timesteps(1:2)));
    ka1_out = sum(first_attachments,2)./Nd/(dt*diff(timesteps(1:2)));
    
    % Find the adjusted and renormalized bond lifetimes
    adjusted_bond_lifetimes = [];
    renormalized_bond_lifetimes = [];
    comp_att = zeros(size(att_mat));    % This will be a composite matrix 
                                        % that holds the history for every
                                        % single bonds. It will contain
                                        % NaNs where no events occur, 1's
                                        % where first time attachments
                                        % occur, 0's where detachments
                                        % occur, 2's where repeat
                                        % attachments occur, and 3's where
                                        % exchange attachments occur
    comp_att(:,:) = NaN;
    comp_att(first_attachments==1) = 1;
    comp_att(detachments==1) = 0;
    comp_att(repeats==1) = 2;
    comp_att(exchanges==1) = 3;
    for j=1:Np
        history = comp_att(:,j);	% contains all of the history for a single
                                    % molecule as decribed above
        new_bonds = find(history==1 | history==3);  % IDs anywhere a new
                                                    % bond pair was formed
        if ~isempty(new_bonds)
            rnm_times_temp = diff(new_bonds)*Delta_t;    % These are renormalized 
                                                    % bond lifetimes
            renormalized_bond_lifetimes = cat(1,renormalized_bond_lifetimes,...
                rnm_times_temp);
            
            % To find adjusted bond lifetimes, go from 1->3 or 3->3 and back to
            % most recent break event prior to the latter
            breaks = find(history==0);   % IDs anywhere bond breakage occured
            for k=2:length(new_bonds)
                new_bond_temp = new_bonds(k);
                breaks_temp = breaks(breaks<new_bond_temp);
                adj_time_temp = (breaks_temp(end) - new_bonds(k-1))*Delta_t;
                adjusted_bond_lifetimes = cat(1,adjusted_bond_lifetimes,...
                    adj_time_temp);
            end
        end
    end
%     attached_bond_lifetimes = [];
%     detached_bond_lifetimes = [];
%     open_exchange_lifetimes = [];
%     open_repeat_lifetimes = [];
%     renormalized_bond_lifetimes = [];
% 
    chain_concentration = Np/(Lx*Ly*Lz)/LengthConversion^3;
%     
%     ever_bonded = all_dynamic_pairs(:,1:2);  % all nodes that are bonded at any point in time
%     ever_bonded = unique(ever_bonded(:));
%     N_bonded = length(ever_bonded);
%     exchange_events = zeros(length(timesteps),N_bonded);
%     all_repeat_times = zeros(length(timesteps),N_bonded);
%     for i=1:N_bonded  %for 1 to the number of nodes that experience bonds at any point in time
%         temp_node = ever_bonded(i);
%         indx1 = find(all_dynamic_pairs(:,1)==temp_node);    %find all instances in which node is first member of pair
%         indx2 = find(all_dynamic_pairs(:,2)==temp_node);    %find all instances in which node is second member of pair
%         all = all_dynamic_pairs([indx1;indx2],:);   %combine these instances
%         all = sortrows(all,3);      %sort by time
%         temp_pairs = all(:,1:2);    %isolate the pairs
%         for t=1:size(temp_pairs,1)
%             indx = find(temp_pairs==temp_node);
%             if indx==2
%                 temp_pairs(t,2) = temp_pairs(t,1);
%                 temp_pairs(t,1) = temp_node;
%             end
%         end
%         all(:,1:2) = temp_pairs;
% 
%         % Define attached bond lifetimes
%         delta_t = diff(all(:,3));   %define difference in time between detachment event and subsequent attachment event
%         delta_neighb = diff(all(:,2));
%         attached_bond_lifetimes = [];
%         detached_bond_lifetimes = [];
%         open_exchange_lifetimes = [];
%         open_repeat_lifetimes = [];
%         renormalized_bond_lifetimes = [];
%         if ~isempty(delta_t)
%             attached_times = delta_t(1:2:end);  %Defines the attached bond lifetimes
%             detached_times = delta_t(2:2:end);  %Defines the detached bond lifetimes
%             
%             repeat_times = detached_times(delta_neighb(2:2:end)==0);
%             repeat_indices = find(delta_neighb(2:2:end)==0)*2+1; 
%             repeat_timesteps = all(repeat_indices,end);
%             all_repeat_times(repeat_timesteps,i) = repeat_times;
% 
%             exchange_times = detached_times(delta_neighb(2:2:end)~=0);
%             exchange_indices = find(delta_neighb(2:2:end)~=0)*2+1;
%             exchange_timesteps = all(exchange_indices,end);
%             if ~isempty(exchange_timesteps)
%                 exchange_events(exchange_timesteps,i) = 1;
%             end
% 
%             attached_bond_lifetimes = cat(1,attached_bond_lifetimes,attached_times);
%             detached_bond_lifetimes = cat(1,detached_bond_lifetimes,detached_times);
%             open_exchange_lifetimes = cat(1,open_exchange_lifetimes,detached_times);
%             open_repeat_lifetimes = cat(1,open_repeat_lifetimes,repeat_times);
%             
%             % Define renormalized bond lifetimes
%             start_t = all(1,3);
%             for t=1:length(delta_neighb)
%                 temp = delta_neighb(t);
%                 if temp~=0 && ~mod(t,2)
%                     renormalized_time = all(t,3)-start_t;
%                     renormalized_bond_lifetimes = ...
%                         cat(1,renormalized_bond_lifetimes,renormalized_time);
%                     next = t+1;
%                     if next<=size(all,1)
%                         start_t = all(next,3);
%                     end
%                 end
%             end
%         end
%     end
           
    total_exchange_events = sum(exchanges,2); % exchange events wrt. time
    mean_exchange_time = mean(open_exchange_lifetimes);
    se_exchange_time = std(open_exchange_lifetimes)/sqrt(length(open_exchange_lifetimes));

%     repeats = all_repeat_times;
%     repeats(repeats~=0) = 1;
    total_repeat_events = sum(repeats,2);
    mean_repeat_time = mean(open_repeat_lifetimes);
    se_repeat_time = std(open_repeat_lifetimes)/sqrt(length(open_repeat_lifetimes));

    exchange_rate = total_exchange_events./Na/dt;
    mean_exchange_rate = nanmean(exchange_rate(:));
    se_exchange_rate = nanstd(exchange_rate(:))/sqrt(length(~isnan(exchange_rate(:))));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kd,ka,ka1,dynamic_pairs0,N_bound0,N_free0,timestep0,...
    all_dynamic_pairs,N_virgin_bonds] = ...
    ComputeBondDynamics(timestep,i,dt,...
    p1,p2,bond_types,bond_id,...
    types,dynamic_pairs0,...
    N_bound0,N_free0,timestep0,tiptyp,dynbondtype,...
    all_dynamic_pairs,N_virgin_bonds)

% Importantly, returns "all_dynamic_pairs", which is
%       [node1 node2 bond_type 
prev_indx = i-1;

% Define timestep
Dt = (timestep-timestep0)*dt;
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
    else
        no_attachments = size(new_pairs,1);   %divide by 2 because Sam's code doesn't delete redundant bonds
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
    msd00_r,msd00_r_err,msd00_t,msd00_t_err]...
    = ComputeMSDs(x,y,z,x0,y0,z0,...
    x00,y00,z00,Corners,...
    types,types0,types00,id,id0,id00,mol,mol0,mol00,...
    stickertype,tethertype)

Lx = abs(Corners(5)-Corners(2));
Ly = abs(Corners(6)-Corners(3));
Lz = abs(Corners(7)-Corners(4));

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
function rs = ExtractEndToEnd(types,id,x,y,z,bond_types,rx,ry,rz,...
    Np,N_Kuhn,i,samp)

global BeadSpringOrMeso bondtypes % Lx Ly Lz

if BeadSpringOrMeso==1
    rx = rx(bond_types~=2);
    ry = ry(bond_types~=2);
    rz = rz(bond_types~=2);
    rs = [rx ry rz];
    
    bondtypes(i,:,samp) = bond_types(bond_types~=2);
else
    atomtypes = types;
    atomnumbers = id;
    atompositions = [x y z];

    NChains = Np;
    rs = zeros(NChains,3);
    for chainnumber=1:NChains
        tether = chainnumber;
        sticker = NChains + N_Kuhn*chainnumber;

        type_t = atomtypes(atomnumbers==tether,1);
        type_s = atomtypes(atomnumbers==sticker,1);
        if type_t~=1 || type_s~=3
            error('Check atom types')
        end

        pos_t = atompositions(atomnumbers==tether,:);
        pos_s = atompositions(atomnumbers==sticker,:);
        rs(chainnumber,:) = pos_s-pos_t;
    end
    bondtypes(i,:,samp) = 4;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PreallocateOutputs(N_steps,N_samples,Np)

global sig11 sig22 sig33 sig12 sig23 sig31...
    sig11_err sig22_err sig33_err sig12_err sig23_err sig31_err...
    bondconc time stretch...
    msd_st msd_st_err...
    msd_st_r0 msd_st_r0_err...
    msd_st_r00 msd_st_r00_err...
    msd_st_t0 msd_st_t0_err...
    msd_st_t00 msd_st_t00_err...
    ka_out ka1_out kd_out Na Nd...
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
kd_out = zeros(N_steps,N_samples);
Na = zeros(N_steps,N_samples);
Nd = zeros(N_steps,N_samples);
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
        (~isfile(EndToEndDataFileName) && CalculateEndtoEnd) ||...
        (~isfile(MSDDataFileName) && CalculateMSD)  ||...
        (~isfile(BondKineticsDataFileName) && CalculateBondKinetics)  ||...
        (~isfile(AlignmentDataFileName) && CalculateAlignment)  ||...
        ~isfile(TimeStretchDataFileName)  ||...
        OverrideCompute
    compute = 1;
end
end
