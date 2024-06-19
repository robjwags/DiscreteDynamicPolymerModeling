function CompileData(Package,OverrideCompile,OverrideCompute,...
    LC,DC,BSOM,CFo,CE,CA,CF,OF,NoProcessors,BT)

% Collects data from all LAMMPS runs and compiles it for ensemble averaging
% or collective histogram outputs

% Compiles stress-time responses
% Stress-strain responses

% End-to-end data for histogram binning
% Bond dynamics data

PackageTypes = DefineUniquePackageTypes(Package);

if NoProcessors==1
    LoopCompileData(Package,OverrideCompile,OverrideCompute,...
        LC,DC,BSOM,CFo,CE,CA,CF,OF,NoProcessors,BT,[])
else
    delete(gcp('nocreate'))
    parpool(NoProcessors)
    parfor n=1:size(PackageTypes,1)
        LoopCompileData(Package,OverrideCompile,OverrideCompute,...
            LC,DC,BSOM,CFo,CE,CA,CF,OF,NoProcessors,BT,PackageTypes(n,:))
    end
    delete(gcp('nocreate'))
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PackageTypes = DefineUniquePackageTypes(Package)

PackageTypes = zeros(1e3,size(Package,2));

Samples = unique(Package(:,1));
Nps = unique(Package(:,2));       
Ds = unique(Package(:,3));       
N_Kuhns = unique(Package(:,4));      
Stiffnesses = unique(Package(:,5));       
kbTs = unique(Package(:,6));       
bs = unique(Package(:,7));      
dts = unique(Package(:,8));      

ct = 0;
for i=1:length(Nps)
 for j=1:length(Ds)
  for k=1:length(N_Kuhns)
   for l=1:length(Stiffnesses)
    for m=1:length(kbTs)
     for n=1:length(bs)
      for o=1:length(dts)

        Np = Nps(i);
        D = Ds(j);
        N_Kuhn = N_Kuhns(k);
        Stiffness = Stiffnesses(l);
        kbT = kbTs(m);
        b = bs(n);
        dt = dts(o);

        PackageTemp = Package(ismember(Package(:,1),Samples(1)),:);
        PackageTemp = Package(ismember(PackageTemp(:,2),Np),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,3),D),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,4),N_Kuhn),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,5),Stiffness),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,6),kbT),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,7),b),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,8),dt),:);

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

PackageTypes(PackageTypes(:,1)==0,:) = [];
PackageTypes(:,end+1) = (1:size(PackageTypes,1))';
PackageTypes(:,end+1) = size(PackageTypes,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LoopCompileData(Package,OverrideCompile,OverrideCompute,...
    LC,DC,BSOM,CFo,CE,CA,CF,OF,NoProcessors,BT,PackageTypes)

global LineWidth FontSize...
    LengthConversion DamperConversion...
    BeadSpringOrMeso RawDataFileName...
    CalculateForce CalculateEndtoEnd...
    CalculateAlignment...
    CurrentFolder OutputFolder BondType

LineWidth = 1.5;
FontSize = 20;
LengthConversion = LC;
DamperConversion = DC;
BeadSpringOrMeso = BSOM;
CalculateForce = CFo;
CalculateEndtoEnd = CE;
CalculateAlignment = CA;
CurrentFolder = CF;
OutputFolder = OF;
BondType = BT;

SimsCt = 0;
if NoProcessors==1
    TotalSims = size(Package,1);
    wb1 = waitbar(0,'Compiling all Data...');
else
    SimTypeNo = PackageTypes(end-1);
    NoPackageTypes = PackageTypes(end);
    disp(['Parameter combo. ',num2str(SimTypeNo),' of ',num2str(NoPackageTypes),' running']);
end

%% Unpack Swept Input Parameters
if NoProcessors==1
    PT = Package;
else
    PT = PackageTypes;
end
Nps = unique(PT(:,2));        %Number of molecules
Ds = unique(PT(:,3));         %Diffusion coefficient [m2/s]
N_Kuhns = unique(PT(:,4));    %Number of Kuhn segments in chain
stiffnesses = unique(PT(:,5)); %Stiffness of single harmonic bond

kbTs = unique(PT(:,6));       %Thermal energy
bs = unique(PT(:,7));         %Kuhn length
dts = unique(PT(:,8));        %timestep


for i=1:length(Nps)
 for j=1:length(Ds)
  for k=1:length(N_Kuhns)
   for l=1:length(stiffnesses)
    for m=1:length(kbTs)
     for n=1:length(bs)
      for o=1:length(dts)

        Np = Nps(i);
        D = Ds(j);
        N_Kuhn = N_Kuhns(k);
        stiffness = stiffnesses(l);
        kbT = kbTs(m);
        b = bs(n);
        dt = dts(o);

        PackageTemp = Package(ismember(Package(:,2),Np),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,3),D),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,4),N_Kuhn),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,5),stiffness),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,6),kbT),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,7),b),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,8),dt),:);

        %% Compile the Raw Data
        for pkg=1:size(PackageTemp,1)
            SimsCt = SimsCt+1;
            if NoProcessors==1
                PcntComplete = SimsCt/TotalSims;
                waitbar(PcntComplete,wb1,'Compiling all Data...')
            end

            Sample = PackageTemp(pkg,1);    %Sample index

            %% Callout Input Script
            % Initializes non-sweeping parameters
            InputScript(Sample,Np,D,N_Kuhn,stiffness,kbT,b,dt);

            %% Make diriectories and set filenames
            SetDirAndFileNames;

            if ~isfile(RawDataFileName) || OverrideCompile==1
                %% Extract data from atoms.dump file
                [timesteps,Corners,Atoms] = UnpackAtomsDump;

                %% Extract data from bonds.dump file
                Bonds = UnpackBondsDump(timesteps);

                %% Save File for Sample
                SaveSampleFiles(Corners,Atoms,Bonds);
            end
        end

        if ~isempty(PackageTemp)
            compute = CheckIfNeedToComputeData(OverrideCompute);
            if compute
                for samp=1:size(PackageTemp,1)
                    %% Unpack Swept Input Parameters
                    Sample = PackageTemp(samp,1);    %Sample index

                    %% Callout Input Script
                    % Initializes non-sweeping parameters
                    InputScript(Sample,Np,D,N_Kuhn,stiffness,kbT,b,dt);

                    %% Make diriectories and set filenames
                    SetDirAndFileNames;

                    %% Import Bonds Data
                    RawData = load(RawDataFileName,'-mat');
                    Corners = RawData.Corners;
                    BondData = RawData.Bonds;

                    %% Allocate sizes of outputs
                    if samp==1
                        N_steps = size(Corners,1);
                        N_samples = size(PackageTemp,1);
                        N_seg = length(unique(BondData(:,2)));
                        PreallocateOutputs(N_steps,N_samples,N_seg,N_Kuhn);
                    end

                    %% Calculate outputs
                    ComputeOutputs(Corners,BondData,...
                        stiffness,N_Kuhn,b,kbT,dt,samp,OverrideCompute);
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

if NoProcessors==1
    close(wb1)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveComputedOutputs(OverrideCompute)

global CalculateForce CalculateEndtoEnd CalculateAlignment...
    ForceDataFileName EndToEndDataFileName AlignmentDataFileName...
    TimeStretchDataFileName...
    time stretch...
    fx fy fz...
    fx_end fy_end fz_end...
    rx ry rz...
    rxr11_st rxr22_st rxr33_st rxr12_st rxr23_st rxr31_st...
    rxr11_st_err rxr22_st_err rxr33_st_err rxr12_st_err rxr23_st_err rxr31_st_err

if (~isfile(TimeStretchDataFileName) || OverrideCompute)
    TimeAndStretch.time = time;
    TimeAndStretch.stretch = stretch;
    save(TimeStretchDataFileName,'-struct','TimeAndStretch')
end

if CalculateForce && (~isfile(ForceDataFileName) || OverrideCompute)
    Forces.fx = fx;
    Forces.fy = fy;
    Forces.fz = fz;

    Forces.fx_end = fx_end;
    Forces.fy_end = fy_end;
    Forces.fz_end = fz_end;

    save(ForceDataFileName,'-struct','Forces')
end

if CalculateEndtoEnd && (~isfile(EndToEndDataFileName) || OverrideCompute)
    EndtoEnd.rx = rx;
    EndtoEnd.ry = ry;
    EndtoEnd.rz = rz;

    save(EndToEndDataFileName,'-struct','EndtoEnd')
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
function PreallocateOutputs(N_steps,N_samples,N_seg,N_Kuhn)

global BeadSpringOrMeso time stretch...
    fx fy fz...
    fx_end fy_end fz_end...
    rx ry rz...
    rxr11_st rxr22_st rxr33_st rxr12_st rxr23_st rxr31_st...
    rxr11_st_err rxr22_st_err rxr33_st_err rxr12_st_err rxr23_st_err rxr31_st_err

% End-to-end vector data
% if BeadSpringOrMeso==0
    N_chains = N_seg;    %rx, ry, rz, for all N_seg segments
% else
%     N_chains = N_seg/N_Kuhn;
% end

% Force vectors
fx_end = zeros(N_steps,N_samples);
fy_end = zeros(N_steps,N_samples);
fz_end = zeros(N_steps,N_samples);

fx = zeros(N_steps,N_chains,N_samples);
fy = zeros(N_steps,N_chains,N_samples);
fz = zeros(N_steps,N_chains,N_samples);

% End-to-end vectors
rx = zeros(N_steps,N_chains,N_samples);
ry = zeros(N_steps,N_chains,N_samples);
rz = zeros(N_steps,N_chains,N_samples);

% Time and stretch
time = zeros(N_steps,1);
stretch = zeros(N_steps,1);

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

global ForceDataFileName EndToEndDataFileName...
    AlignmentDataFileName...
    TimeStretchDataFileName...
    CalculateForce CalculateEndtoEnd...
    CalculateAlignment...

compute = 0;
if (~isfile(ForceDataFileName) && CalculateForce)  ||...
        (~isfile(EndToEndDataFileName) && CalculateEndtoEnd)  ||...
        (~isfile(AlignmentDataFileName) && CalculateAlignment)  ||...
        ~isfile(TimeStretchDataFileName)  ||...
        OverrideCompute
    compute = 1;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ComputeOutputs(Corners,BondData,stiffness,N_Kuhn,b,kbT,dt,samp,...
    OverrideCompute)

global CalculateForce CalculateEndtoEnd CalculateAlignment...
    ForceDataFileName EndToEndDataFileName AlignmentDataFileName...
    time stretch...
    fx fy fz...
    fx_end fy_end fz_end...
    rx ry rz...
    rxr11_st rxr22_st rxr33_st rxr12_st rxr23_st rxr31_st...
    rxr11_st_err rxr22_st_err rxr33_st_err rxr12_st_err rxr23_st_err rxr31_st_err

% timestep id bond_type node1 node2 rx ry rz norm force
timesteps = unique(BondData(:,1));
if size(Corners,1)~=length(timesteps)
    Diff = size(Corners,1)-length(timesteps);
    if Diff<0
        timesteps(end+Diff+1:end) = [];
    end
end

%% Initialize dynamics variables
wb2 = waitbar(0,'Computing stress, end-to-end, MSD, & dynamics data...');
for i=1:length(timesteps)
    timestep = timesteps(i);
    time(i) = timestep*dt;
    waitbar(i/length(timesteps),wb2,...
        'Computing stress, end-to-end, MSD, & dynamics data...')
    rs = []; 
    BondDataTemp = [];

    %% Store end-to-end components & Bond types
    if CalculateEndtoEnd && (~isfile(EndToEndDataFileName) || OverrideCompute)
        BondDataTemp = BondData(BondData(:,1)==timestep,:);
        [rs,stretch(i)] = ExtractEndToEnd(BondDataTemp,N_Kuhn,b);
        rx(i,:,samp) = rs(:,1); 
        ry(i,:,samp) = rs(:,2); 
        rz(i,:,samp) = rs(:,3);
    end
    
    %% Compute metric tensors by bond type
    if CalculateAlignment && (~isfile(AlignmentDataFileName) || OverrideCompute)
        if isempty(BondDataTemp)
            BondDataTemp = BondData(BondData(:,1)==timestep,:);
        end
        if isempty(rs)
            [rs,stretch(i)] = ExtractEndToEnd(BondDataTemp,N_Kuhn,b);
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
    if CalculateForce && (~isfile(ForceDataFileName) || OverrideCompute)
        if isempty(BondDataTemp)
            BondDataTemp = BondData(BondData(:,1)==timestep,:);
        end
        if isempty(rs)
            [rs,stretch(i)] = ExtractEndToEnd(BondDataTemp,N_Kuhn,b);
            rx(i,:,samp) = rs(:,1);
            ry(i,:,samp) = rs(:,2);
            rz(i,:,samp) = rs(:,3);
        end
        forces = ComputeForces(rs,N_Kuhn,b,kbT,stiffness);
        fx(i,:,samp) = forces(:,1);
        fy(i,:,samp) = forces(:,2);
        fz(i,:,samp) = forces(:,3);

        f_end = ComputeForces(rs(end,:),N_Kuhn,b,kbT,stiffness);
        fx_end(i,:,samp) = f_end(1);
        fy_end(i,:,samp) = f_end(2);
        fz_end(i,:,samp) = f_end(3);
    end
end
close(wb2)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rs,stretch] = ExtractEndToEnd(BondDataTemp,N_Kuhn,b)

rs = BondDataTemp(:,6:8);
rnet = sum(rs,1);
stretch = vecnorm(rnet,2,2)/(sqrt(N_Kuhn)*b);

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
function forces = ComputeForces(r,N_Kuhn,b,kbT,stiffness)

global LengthConversion BeadSpringOrMeso BondType

b = b*LengthConversion;
r = r*LengthConversion;
norms = vecnorm(r,2,2);
r_hat = [r(:,1)./norms r(:,2)./norms r(:,3)./norms];

if BeadSpringOrMeso==0
    if BondType==0
        K = stiffness*kbT/b^2;  % Converted back to SI units
        force_mags = K*(norms-b);
    else
        E = stiffness*kbT;  % Converted back to SI units
        L = b;
        r0 = b;
        r_samp = linspace(0,2*b,50);
        Psi = E*((r_samp-r0).^2)./(L^2-(r_samp-r0).^2);
        f_samp = gradient(Psi)./gradient(r_samp);
        force_mags = interp1(r_samp,f_samp,norms);
    end
else
%     K = 3*kbT/(N_Kuhn*b^2);
%     force_mags = K*norms;

    lam = norms/N_Kuhn/b;
    numer = lam*(3 - lam.^2);
    denom = 1 - lam.^2;
    force_mags = -kbT*numer/denom/b;
end
forces = force_mags.*r_hat;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveSampleFiles(Corners,Atoms,Bonds)

global RawDataFileName %AtomsFileName BondsFileName CornersFileName 

% structArray = cell2struct(cellArray, fields, dim)

Raw.Corners = Corners;
Raw.Atoms = Atoms;
Raw.Bonds = Bonds;
save(RawDataFileName,'-struct','Raw')

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
for i=1:N_output
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
