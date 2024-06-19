% Branched newtork topolgoy generator for input to LAMMPS via read_data cmd
% Robert J. Wagner
% 2022-11-29

% This version generates the x-links first and then connects

function GenerateTopology(package,Parameters,BSOM,...
    OverrideTopologies,OverrideInputScripts,Directories,n,package_tot)

global samples eq_time_factors N_Kuhns phis kds Weissenbergs...
    sample eq_time_factor N_Kuhn phi kd Weissenberg...
    omega omegas...
    dt Np kbT kb T b ka...
    length_conversion damper_conversion...
    damp D tau0 dtFact...
    dims ShowFiguresDebug...
    OldOrNew PlotStickers MakeMovie...
    TurnOnDynamics full_topology_filename...
    UnitType input_filename BeadSpringOrMeso...
    current_folder WSL_path...
    check_Kuhn_lengths LinearOrLangevin TurnOffLJ initial_matlab_filename...
    Lx Ly Lz...
    

Controls = DefineControlSwitches(BSOM);
if Controls.RunTheLJCase == 1
    TurnOffLJ = 0;
else
    TurnOffLJ = 1;
end

% Unpack package and parameters into globals
BeadSpringOrMeso = BSOM;
UnpackParameters(package,Parameters);

% Since GenerateTopology is swept at a functional level above this script,
% individual values of swept parameters are the same as the ranges defined
% in UnpackParameters;
sample = samples;
eq_time_factor = eq_time_factors;
N_Kuhn = N_Kuhns;
phi = phis;
kd = kds;
Weissenberg = Weissenbergs;
omega = omegas;

[damp,D,tau0,dtFact] = DefineTimeScale(b,length_conversion,damper_conversion,...
    BeadSpringOrMeso,phi,N_Kuhn);
dt = tau0/dtFact;

DefineNumbersOfTethersAndParticles;
DefineNormalizedBindingEnergies;

% Set Toggles for MATLAB topology generation
OldOrNew = 1;           % Different looping options for two diff.
% connectivity algorithms (ignore)
ShowFiguresDebug = 0;   % 1 to plot intermediate figures for debugging
TurnOnDynamics = Controls.ToggleDynamics;
UnitType = Controls.UnitType;
current_folder = Directories.current_folder;
WSL_path = Directories.WSL_path;
LinearOrLangevin = Controls.LinearOrLangevin;
check_Kuhn_lengths = 0;



%% Callout Input Script
% Initializes non-sweeping parameters
dims = 3;           %Dimensionality of system
DefineMatlabNetworkInitiationInputs;

%% Make diriectories and set filenames
DefineFileNames(Controls);    % Can either rely on globals or output structure
% with all filenames

if ~isfile(full_topology_filename) || OverrideTopologies==1

    if BeadSpringOrMeso==0 || ~isfile(initial_matlab_filename) ...
            || OverrideTopologies==1

        %         if ~isfile(initial_matlab_filename) || OverrideTopologies==1
        %% Initialize the positions of all tethers
        SetTheInitialTetherPositions(OverrideTopologies,Controls);

        %% Polymerize the backbone and equilibrate
        EquilibrateTheBackBone(OverrideTopologies);

        %% Polymerize and equilibrate the branches
        PolymerzieAndEquilibrateBranches(OverrideTopologies);

        if Controls.RunLargeDeformation~=1
            InterpolateBeadPositions(OverrideTopologies);
        end

        %% Unwrap nodes periodicity
        UnwrapPeriodic(OverrideTopologies);

        %% Output the Langevin chain table
        if LinearOrLangevin==1
            DefineLangevinForces(N_Kuhn,b);
            DefineFENEForces(b);
        end
        %         end
        SaveMatlabData;         % For seeding the mesoscale model

    end
    if BeadSpringOrMeso==1      % Take bead-spring initial config and delete
        % intermediate beads
        %% Import the initial position daata
        InitialData = load(initial_matlab_filename,'-mat');
        InitialPos = InitialData.Positions;
        id = InitialPos(:,1);
        mol = InitialPos(:,2);
        types = InitialPos(:,3);
        x = InitialPos(:,4);
        y = InitialPos(:,5);
        z = InitialPos(:,6);
        PosUnwrapped = InitialData.PosUnwrapped;

        InitialSize = InitialData.Boundaries;
        Lx = InitialSize(1); Ly = InitialSize(2); Lz = InitialSize(3);

        InitialBonds = InitialData.Bonds;
        bond_types = InitialBonds(:,2);
        batom1 = InitialBonds(:,3);
        batom2 = InitialBonds(:,4);

        %% Define end-to-end connections and positions for the mesoscale
        % using the network configuration and conformation from the
        % bead-spring inititation procedure.
        if Controls.RunLargeDeformation~=1
            ReestablishEndtoEnds(types,id,mol,x,y,z,bond_types,batom1,batom2,...
                Np,Lx,Ly,Lz,PosUnwrapped)
        end
    end

    %% Output Topological Data
    PlotStickers=1;
    SaveTopologyDataLAMMPS(N_Kuhn,kb,T,b,Lx,Ly,Lz);

    if MakeMovie==1
        MakeMovies(N_Kuhn);
    end

end

% %% Output force-distance table for Langevin & Units
% SaveForceData;

%% CHECK WHAT KBT SHOULD BE REVISIT
if ~isfile(input_filename) || OverrideInputScripts==1 ...
        || OverrideTopologies==1
    PrintInputFiles(ka,kd,kbT,T,Controls,n,package_tot);
end

% toc
close all
clear
clc
clear global

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rs,batom1_out,batom2_out] = ReestablishEndtoEnds(types,id,mol,...
    x,y,z,bond_types,batom1,batom2,Np,L,W,H,PosUnw)

global Nt Ns PosUnwrapped ConnectionsSP N_atoms N_bonds...

CheckFig = 0;

% Define rx,ry,rz
rx = zeros(size(batom1));
ry = zeros(size(batom1));
rz = zeros(size(batom1));
for i=1:length(batom1)
    p1 = batom1(i); p2 = batom2(i);
    X1 = [x(id==p1),y(id==p1),z(id==p1)];
    X2 = [x(id==p2),y(id==p2),z(id==p2)];
    r = X2-X1;
    rx(i) = r(1); ry(i) = r(2); rz(i) = r(3);
end

% Initialize all vectors
Ntethers = Nt-1;
Nstickers = Nt*Ns; % per molecule
Nsegments = Np*(Ntethers+Nstickers);
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
    sticker_rng = Nt+(1:Nstickers);
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
            %% STORE BATOMS BTYPES

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
                daspect([1 1 1])
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
x_out = x(types~=3);
% y_out = y(types~=3);
% z_out = z(types~=3);

% Re-assemble into PosUnwrapped and ConnectionsSP
N_atoms = size(x_out,1);
N_bonds = size(batom1_out,1);

% Update unwrapped positions
PosUnwrapped = PosUnw;
rem_ids = find(types==3);
PosUnwrapped(ismember(PosUnwrapped(:,4),rem_ids),:) = [];

% Update connections
ConnectionsSP = zeros(size(PosUnwrapped,1),5);
ConnectionsSP(:,end) = (1:N_atoms)';

for i=1:N_bonds
    from = batom1_out(i);
    to = batom2_out(i);

    ConnectionsSP(from,find(ConnectionsSP(from,:)==0,1,'first')) = to;
    ConnectionsSP(to,find(ConnectionsSP(to,:)==0,1,'first')) = from;
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
function DefineFENEForces(b)

global FENE_force_filename

Npts = 100;

eps = 100;
r0 = b;
lam = b;

r = (linspace(0.01*b,1.99*b,Npts))';
E = eps*(r-r0).^2./(lam^2-(r-r0).^2);
f = gradient(E)./gradient(r);

check_fig = 1;
if check_fig==1
    figure(1e4); clf; hold on
    plot(r,f,'k')
    plot(r,E,'r');
    ylim([-1e3 1e3])
    close(gcf)
end

CreateBondTableFile(FENE_force_filename,r,E,f)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DefineLangevinForces(N_Kuhn,b)

global langeving_force_filename

Npts = 100;
pcnt_tangent = 97.5;
r = (linspace(0,pcnt_tangent/100*N_Kuhn*b,Npts))';
lambda = r/(sqrt(N_Kuhn)*b);
f_lang = lambda/(sqrt(N_Kuhn)*b).*...
    ((lambda.^2-3*N_Kuhn)./(lambda.^2-N_Kuhn));
f_lin = 3*r/(N_Kuhn*b^2);

% add linear slope
p1 = [r(end-1) f_lang(end-1)];
p2 = [r(end) f_lang(end)];
m = (p2(2)-p1(2))/(p2(1)-p1(1));
y_int = p1(2)-m*p1(1);

r = (linspace(0,1.1*N_Kuhn*b,Npts))';
lambda = r/(sqrt(N_Kuhn)*b);
f_lang = lambda/(sqrt(N_Kuhn)*b).*...
    ((lambda.^2-3*N_Kuhn)./(lambda.^2-N_Kuhn));

f_lang(r>=pcnt_tangent/100*N_Kuhn*b) = m*(r(r>=pcnt_tangent/100*N_Kuhn*b))+y_int;

E = zeros(size(r));
for i=2:length(r)
    E(i) = trapz(r(1:i),f_lang(1:i));
end
f_check = gradient(E)./gradient(r);

check_fig = 1;
if check_fig==1
    figure(1e4); clf; hold on
    plot(r,f_lin,'k--')
    plot(r,f_lang,'k');
    plot(r,f_check,'r')
    close(gcf)
end

CreateBondTableFile(langeving_force_filename,r,E,f_lang)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CreateBondTableFile(filename, bondLength, energy, force)

% Open the file for writing
fileID = fopen(filename, 'w');

% Write the header lines
fprintf(fileID, '# Bond potential for harmonic (one or more comment or blank lines)\n');
fprintf(fileID, 'HAM\n');
fprintf(fileID, '\n');
fprintf(fileID, ['N ',num2str(size(bondLength,1)),'\n']);
fprintf(fileID, '\n');

% Write the data
nDataPoints = length(bondLength);
for i = 1:nDataPoints
    fprintf(fileID, '%d %.2f %.4f %.4f\n', i, bondLength(i), energy(i), force(i));
    %     fprintf(fileID, '%.2f %.4f\n', bondLength(i), energy(i));
end

% Close the file
fclose(fileID);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InterpolateBeadPositions(Override)

global N_bead dims b CutoffCWAForce Positions X_SP Y_SP Z_SP...
    ConnectionsSP sigmaREP CutoffPartREP Rx_SP Ry_SP Rz_SP...
    b_eq L_eq SpringCoeff_eq kbT_eq ResAvg ResMax check_Kuhn_lengths



if size(Positions,1)~=N_bead || Override

    ImplementEularianBoundaryConditions;
    FindBoundaryNodes;
    ReplicateBounds;

    all_pairs = zeros(size(ConnectionsSP,1)*(size(ConnectionsSP,2)-1),2);
    for i=1:size(ConnectionsSP,2)-1
        rng = size(ConnectionsSP,1)*(i-1)+1:size(ConnectionsSP,1)*i;
        froms = ConnectionsSP(:,end);
        tos = ConnectionsSP(:,i);
        pairs = [tos,froms];
        all_pairs(rng,:) = pairs;

    end
    all_pairs(all_pairs(:,1)==0,:) = [];
    all_pairs = sort(all_pairs,2);
    all_pairs = unique(all_pairs,'rows');

    N_bonds = size(all_pairs,1);
    for i=1:N_bonds
        % Interpolate N_beads between each pair of connected nodes and
        % append these to the Positions, PositionsUnwrapped

        %% ACCOUNT FOR PERIODIC BOUNDARIES IF CROSSES A BOUNDARY
        temp_pair = all_pairs(i,:);
        p1 = temp_pair(1); p2 = temp_pair(2);
        [add_x_all,add_y_all,add_z_all] = InterpolateBeads(p1,p2);
        add_x = add_x_all(2:end-1);
        add_y = add_y_all(2:end-1);
        add_z = add_z_all(2:end-1);

        strt = size(Positions,1)+1;
        add_rng = (strt:strt+length(add_x)-1)';
        Positions(add_rng,1) = add_x';
        Positions(add_rng,2) = add_y';
        Positions(add_rng,3) = add_z';

        % ID molecule number
        id_prev = find(Positions(:,4)==0,1,'first');
        add_ids = (id_prev:id_prev+length(add_x)-1)';
        Positions(add_rng,4) = add_ids;

        % Add mer number
        mol_no = Positions(p1,5);
        mer_rows = find(Positions(:,5)==mol_no);
        mer_ct = length(mer_rows);
        add_mers = (mer_ct+1:mer_ct+length(add_x))';
        Positions(add_rng,end) = add_mers;

        % Molecule numbers
        Positions(add_rng,5) = mol_no;

        % Adjust ConnectionsSP
        % 1) Remove old pairs from ConnectionsSP
        ConnectionsSP(p1,ConnectionsSP(p1,:)==p2) = 0;
        ConnectionsSP(p2,ConnectionsSP(p2,:)==p1) = 0;
        ConnectionsSP(add_rng,:) = 0;
        ConnectionsSP(add_rng,end) = add_ids;

        % 2) Define list of beads with new neighbors in order of
        % connectivity and populate ConnectionsSP
        beads_with_new = [p1;add_ids;p2];
        for j=1:length(beads_with_new)-1
            Rx = diff(add_x_all);
            Ry = diff(add_y_all);
            Rz = diff(add_z_all);

            new_p1 = beads_with_new(j);
            new_p2 = beads_with_new(j+1);
            col = find(ConnectionsSP(new_p1,:)==0,1,'first');
            ConnectionsSP(new_p1,col) = new_p2;
            pos1 = Positions(new_p1,1:dims);
            pos2 = Positions(new_p2,1:dims);
            Rx_SP(new_p1,col) = -Rx(j);
            Ry_SP(new_p1,col) = -Ry(j);
            Rz_SP(new_p1,col) = -Rz(j);
            X_SP(new_p1,col) = pos1(1);
            Y_SP(new_p1,col) = pos1(2);
            Z_SP(new_p1,col) = pos1(3);

            col = find(ConnectionsSP(new_p2,:)==0,1,'first');
            ConnectionsSP(new_p2,col) = new_p1;
            Rx_SP(new_p2,col) = Rx(j);
            Ry_SP(new_p2,col) = Ry(j);
            Rz_SP(new_p2,col) = Rz(j);
            X_SP(new_p2,col) = pos2(1);
            Y_SP(new_p2,col) = pos2(2);
            Z_SP(new_p2,col) = pos2(3);
        end
    end

    % Adjust force parameters for re-equilibration
    sigmaREP = b; CutoffPartREP = b; CutoffCWAForce = 0.5*b;
    L_eq = b; fact = 0.25;
    SpringCoeff_eq = fact*3*kbT_eq/(L_eq*b_eq);
    ResMax = ResMax/2;
    ResAvg = ResAvg/2;
    check_Kuhn_lengths = 1;

    ConnectREP
    ImplementEularianBoundaryConditions
    FindBoundaryNodes
    ReplicateBounds
    ConnectREP
    EquilibrateTheDomain;

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [add_x,add_y,add_z] = InterpolateBeads(p1,p2)

global LeftBnds RightBnds FrontBnds BackBnds BottomBnds TopBnds...
    Lx Ly Lz N_Kuhn Corners Positions

pos1 = Positions(p1,1:3);
pos2 = Positions(p2,1:3);
pos1_0 = pos1;
pos2_0 = pos2;

if ~isempty(intersect(LeftBnds,RightBnds)) || ...
        ~isempty(intersect(TopBnds,BottomBnds)) || ...
        ~isempty(intersect(FrontBnds,BackBnds))
    warning('Overlap in boundary nodes accross opposite bounds')
end

flag1 = 0;
flag2 = 0;
flag3 = 0;
flag4 = 0;
flag5 = 0;
flag6 = 0;
if ismember(p1,LeftBnds) && ismember(p2,RightBnds)
    pos1 = pos1 + [Lx 0 0];
    flag1 = 1;
elseif ismember(p1,RightBnds) && ismember(p2,LeftBnds)
    pos1 = pos1 - [Lx 0 0];
    flag2 = 1;
end

if ismember(p1,BackBnds) && ismember(p2,FrontBnds)
    pos1 = pos1 - [0 Ly 0];
    flag3 = 1;
elseif ismember(p1,FrontBnds) && ismember(p2,BackBnds)
    pos1 = pos1 + [0 Ly 0];
    flag4 = 1;
end

if ismember(p1,BottomBnds) && ismember(p2,TopBnds)
    pos1 = pos1 + [0 0 Lz];
    flag5 = 1;
elseif ismember(p1,TopBnds) && ismember(p2,BottomBnds)
    pos1 = pos1 - [0 0 Lz];
    flag6 = 1;
end

% add_x = InterpolateChronologically(pos1(1),pos2(1),N_Kuhn-1);
% add_y = InterpolateChronologically(pos1(2),pos2(2),N_Kuhn-1);
% add_z = InterpolateChronologically(pos1(3),pos2(3),N_Kuhn-1);
add_x = InterpolateChronologically(pos1(1),pos2(1),N_Kuhn);
add_y = InterpolateChronologically(pos1(2),pos2(2),N_Kuhn);
add_z = InterpolateChronologically(pos1(3),pos2(3),N_Kuhn);

CheckFigs = 0;
if CheckFigs && ...
        (flag1 || flag2 || flag3 || flag4 || flag5 || flag6)
    figure(100); clf; hold on
    %Plot Window
    C1 = Corners(1,:);
    C2 = Corners(2,:);
    C3 = Corners(3,:);
    C4 = Corners(4,:);
    C5 = Corners(5,:);
    C6 = Corners(6,:);
    C7 = Corners(7,:);
    C8 = Corners(8,:);

    plot3([C1(1) C2(1)],[C1(2) C2(2)],[C1(3) C2(3)],'k--','LineWidth',1.5)
    plot3([C2(1) C3(1)],[C2(2) C3(2)],[C2(3) C3(3)],'k--','LineWidth',1.5)
    plot3([C3(1) C4(1)],[C3(2) C4(2)],[C3(3) C4(3)],'k--','LineWidth',1.5)
    plot3([C4(1) C1(1)],[C4(2) C1(2)],[C4(3) C1(3)],'k--','LineWidth',1.5)
    plot3([C1(1) C5(1)],[C1(2) C5(2)],[C1(3) C5(3)],'k--','LineWidth',1.5)
    plot3([C2(1) C6(1)],[C2(2) C6(2)],[C2(3) C6(3)],'k--','LineWidth',1.5)
    plot3([C3(1) C7(1)],[C3(2) C7(2)],[C3(3) C7(3)],'k--','LineWidth',1.5)
    plot3([C4(1) C8(1)],[C4(2) C8(2)],[C4(3) C8(3)],'k--','LineWidth',1.5)
    plot3([C5(1) C6(1)],[C5(2) C6(2)],[C5(3) C6(3)],'k--','LineWidth',1.5)
    plot3([C6(1) C7(1)],[C6(2) C7(2)],[C6(3) C7(3)],'k--','LineWidth',1.5)
    plot3([C7(1) C8(1)],[C7(2) C8(2)],[C7(3) C8(3)],'k--','LineWidth',1.5)
    plot3([C8(1) C5(1)],[C8(2) C5(2)],[C8(3) C5(3)],'k--','LineWidth',1.5)
    view([30,45])

    scatter3(pos1_0(1),pos1_0(2),pos1_0(3),'r')
    scatter3(pos1(1),pos1(2),pos1(3),'r','filled')
    scatter3(pos2_0(1),pos2_0(2),pos2_0(3),'g','filled')
    scatter3(add_x,add_y,add_z,'c','filled')

    daspect([1 1 1])
    set(gcf,'color','w')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = InterpolateChronologically(p1,p2,npts)

if p1>p2
    out = linspace(p2,p1,npts);
    out = fliplr(out);
else
    out = linspace(p1,p2,npts);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  PolymerzieAndEquilibrateBranches(Override)

global initial_node_pos_filename...
    Positions Positions0...
    X_SP Y_SP Z_SP...
    X_REP Y_REP Z_REP...
    ConnectionsSP ConnectionsREP sigmaREP CutoffPartREP...
    Rx_SP Ry_SP Rz_SP...
    Rx_REP Ry_REP Rz_REP...
    Corners Corners0...
    ResMax ResAvg N_meso...
    PlotStickers ShowFiguresDebug

% Add branches
dat = load(initial_node_pos_filename,'-mat');
Positions = dat.Positions;
if size(Positions,1)~=N_meso || Override
    PolymerizeBranches;

    %Equilibrate to reposition nodes based on newly introduced stickers
    ResMax = 1.5;
    ResAvg = 0.9;
    EquilibrateTheDomain;

    InitialTetherPos.Positions = Positions;
    InitialTetherPos.Positions0 = Positions0;
    InitialTetherPos.Corners0 = Corners0;
    InitialTetherPos.Corners = Corners;
    InitialTetherPos.X_SP = X_SP;
    InitialTetherPos.Y_SP = Y_SP;
    InitialTetherPos.Z_SP = Z_SP;
    InitialTetherPos.X_REP = X_REP;
    InitialTetherPos.Y_REP = Y_REP;
    InitialTetherPos.Z_REP = Z_REP;
    InitialTetherPos.ConnectionsSP = ConnectionsSP;
    InitialTetherPos.ConnectionsREP = ConnectionsREP;
    InitialTetherPos.Rx_SP = Rx_SP;
    InitialTetherPos.Ry_SP = Ry_SP;
    InitialTetherPos.Rz_SP = Rz_SP;
    InitialTetherPos.Rx_REP = Rx_REP;
    InitialTetherPos.Ry_REP = Ry_REP;
    InitialTetherPos.Rz_REP = Rz_REP;
    InitialTetherPos.sigmaREP = sigmaREP;
    InitialTetherPos.CutoffPartREP = CutoffPartREP;

    save(initial_node_pos_filename,'-struct','InitialTetherPos','-v7.3')
else
    Positions = dat.Positions;
    Positions0 = dat.Positions0;
    X_SP = dat.X_SP;
    Y_SP = dat.Y_SP;
    Z_SP = dat.Z_SP;
    X_REP = dat.X_REP;
    Y_REP = dat.Y_REP;
    Z_REP = dat.Z_REP;
    ConnectionsSP = dat.ConnectionsSP;
    ConnectionsREP = dat.ConnectionsREP;
    Rx_SP = dat.Rx_SP;
    Ry_SP = dat.Ry_SP;
    Rz_SP = dat.Rz_SP;
    Rx_REP = dat.Rx_REP;
    Ry_REP = dat.Ry_REP;
    Rz_REP = dat.Rz_REP;
    sigmaREP = dat.sigmaREP;
    CutoffPartREP = dat.CutoffPartREP;
end

PlotStickers = 1;
if ShowFiguresDebug==1
    figure(4)
    clf
    PlotNetwork
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EquilibrateTheBackBone(Override)

global initial_node_pos_filename...
    Positions Positions0...
    X_SP Y_SP Z_SP...
    X_REP Y_REP Z_REP...
    ConnectionsSP ConnectionsREP sigmaREP CutoffPartREP...
    Rx_SP Ry_SP Rz_SP...
    Rx_REP Ry_REP Rz_REP...
    Corners Corners0...
    PlotChains ShowFiguresDebug

dat = load(initial_node_pos_filename,'-mat');
if isfield(dat,'X_SP') && ~Override
    Positions = dat.Positions;
    Positions0 = dat.Positions0;
    X_SP = dat.X_SP;
    Y_SP = dat.Y_SP;
    Z_SP = dat.Z_SP;
    X_REP = dat.X_REP;
    Y_REP = dat.Y_REP;
    Z_REP = dat.Z_REP;
    ConnectionsSP = dat.ConnectionsSP;
    ConnectionsREP = dat.ConnectionsREP;
    Rx_SP = dat.Rx_SP;
    Ry_SP = dat.Ry_SP;
    Rz_SP = dat.Rz_SP;
    Rx_REP = dat.Rx_REP;
    Ry_REP = dat.Ry_REP;
    Rz_REP = dat.Rz_REP;
    sigmaREP = dat.sigmaREP;
    CutoffPartREP = dat.CutoffPartREP;
else
    % Polymerize the backbone
    PolymerizeChains;

    % Equilibrate to alleviate excess force in any long springs
    EquilibrateTheDomain;

    % Append the the data
    InitialTetherPos.Positions = Positions;
    InitialTetherPos.Positions0 = Positions0;
    InitialTetherPos.Corners0 = Corners0;
    InitialTetherPos.Corners = Corners;
    InitialTetherPos.X_SP = X_SP;
    InitialTetherPos.Y_SP = Y_SP;
    InitialTetherPos.Z_SP = Z_SP;
    InitialTetherPos.X_REP = X_REP;
    InitialTetherPos.Y_REP = Y_REP;
    InitialTetherPos.Z_REP = Z_REP;
    InitialTetherPos.ConnectionsSP = ConnectionsSP;
    InitialTetherPos.ConnectionsREP = ConnectionsREP;
    InitialTetherPos.Rx_SP = Rx_SP;
    InitialTetherPos.Ry_SP = Ry_SP;
    InitialTetherPos.Rz_SP = Rz_SP;
    InitialTetherPos.Rx_REP = Rx_REP;
    InitialTetherPos.Ry_REP = Ry_REP;
    InitialTetherPos.Rz_REP = Rz_REP;
    InitialTetherPos.sigmaREP = sigmaREP;
    InitialTetherPos.CutoffPartREP = CutoffPartREP;

    save(initial_node_pos_filename,'-struct','InitialTetherPos','-v7.3')
end


PlotChains = 1;
if ShowFiguresDebug==1
    figure(3)
    clf
    PlotNetwork
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetTheInitialTetherPositions(Override,Controls)

global dims Corners Lx Ly Lz Corners0 PlotChains PlotStickers...
    ShowFiguresDebug initial_node_pos_filename Positions...
    sigmaREP CutoffPartREP Rx_REP Ry_REP Rz_REP

if ~isfile(initial_node_pos_filename) || Override
    % Size the domain
    % Determines correct domain size based on no. xls and nominal density
    SizeTheDomain;

    % Initialize the Nodes
    InitializeNodes(Controls);

    % Define Boundaries
    FindBoundaryNodes;  %ID the Boundary Nodes %% UPDATE 3D

    % Replicate bounds with node types
    ReplicateBounds;	%Replicate boundary Nodes accros periodic bounds

    % Generate repulsive connectivity
    ConnectREP;         %Defines Repulsion Connection Matrix

    % Define window at start of deformation
    if dims==2
        Corners = [-Lx/2 -Ly/2;     %BL
            Lx/2 -Ly/2;     %BR
            Lx/2  Ly/2;     %TR
            -Lx/2  Ly/2];    %TL
    elseif dims==3
        Corners = [-Lx/2 -Ly/2 -Lz/2;   %C1
            +Lx/2 -Ly/2 -Lz/2;  %C2
            +Lx/2 +Ly/2 -Lz/2;  %C3
            -Lx/2 +Ly/2 -Lz/2;  %C4
            -Lx/2 -Ly/2 +Lz/2;  %C5
            +Lx/2 -Ly/2 +Lz/2;  %C6
            +Lx/2 +Ly/2 +Lz/2;  %C7
            -Lx/2 +Ly/2 +Lz/2]; %C8
    end
    Corners0 = Corners;

    PlotChains = 0;
    PlotStickers = 0;
    if ShowFiguresDebug==1
        figure(1)
        clf
        PlotNetwork
    end

    % Calculate Forces, Step Nodes and Equilibrate Domain to Low Energy State
    EquilibrateTheDomain;   %Defines all forces, steps positions of nodes and
    %repeats until residual forces diminish to threshold

    if ShowFiguresDebug==1
        figure(2)
        clf
        PlotNetwork
    end

    InitialTetherPos.Corners0 = Corners0;
    InitialTetherPos.Corners = Corners;
    InitialTetherPos.Positions = Positions;
    InitialTetherPos.sigmaREP = sigmaREP;
    InitialTetherPos.CutoffPartREP = CutoffPartREP;
    InitialTetherPos.Rx_REP = Rx_REP;
    InitialTetherPos.Ry_REP = Ry_REP;
    InitialTetherPos.Rz_REP = Rz_REP;

    save(initial_node_pos_filename,'-struct','InitialTetherPos','-v7.3')
else
    dat = load(initial_node_pos_filename,'-mat');
    Corners = dat.Corners;
    Corners0 = dat.Corners0;
    Positions = dat.Positions;
    sigmaREP = dat.sigmaREP;
    CutoffPartREP = dat.CutoffPartREP;
    Rx_REP = dat.Rx_REP;
    Ry_REP = dat.Ry_REP;
    Rz_REP = dat.Rz_REP;
    C1 = Corners(1,:);
    C2 = Corners(2,:);
    C4 = Corners(4,:);
    C5 = Corners(5,:);
    Lx = norm(C2-C1);
    Ly = norm(C1-C4);
    Lz = norm(C5-C1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PrintInputFiles(ka0,kd0,kbT,T,Controls,n,package_tot)

global full_topology_filename partial_topology_filename...
    OutputAtom OutputBond...
    OutputAtom_eq OutputBond_eq...
    OutputAtom_ld OutputBond_ld...
    OutputAtom_rlx OutputBond_rlx...
    input_filename...
    stress_ldg_filename stress_rlx_filename stress_eq_filename TurnOnDynamics...
    UnitType energy_conversion...
    BeadSpringOrMeso LinearOrLangevin TurnOffLJ...
    N_Kuhn b tau0 damp dt lambda Weissenberg eq_time_factor phi...
    current_folder WSL_path Nt...
    omega_min omega_max tf Lx Ly Lz Np omega nmax_tot

SizeTheDomain;

if BeadSpringOrMeso==1 % Mesoscale
    damp_st = damp*N_Kuhn/2;        % Effective damper needed to roughly
                                    % mimic bead-spring diffusion of
                                    % stickers
    damp_th = damp*N_Kuhn/2*Nt;     % Effective damper assumed needed to
                                    % roughly mimic bead-spring diffusion
                                    % of tethers. Nt is the functionality
                                    % of the crosslink sites at the
                                    % tethers, so this states that each
                                    % tether lumps half of its attached
                                    % chains.
    Ldyn = b;                       % Bond attachment cutoff lengthscale
elseif BeadSpringOrMeso==0 % Bead-spring
    damp_th = damp;
    damp_st = damp;
    Ldyn = b;
end
kbT = kbT*energy_conversion;

% Initialize string
Assembled = [];

% Define units
L = '####################################'; Assembled = [Assembled,L];
L = '\n# DEFINE UNITS'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
if UnitType==0
    L = '\nunits        si'; Assembled = [Assembled,L];
elseif UnitType==1
    L = '\nunits        lj'; Assembled = [Assembled,L];
else
    L = '\nunits        nano'; Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% Define lengthscales
Nb = N_Kuhn*b;              %chain length
Nb2 = Nb*b;                 %N*b^2
r0 = sqrt(N_Kuhn)*b;        %mean random walk chain length
Kstiff = 3*kbT/Nb2;         %chain stiffness

if BeadSpringOrMeso==0
    Lmax = 1.1*b;           %max LJ cutoff - set very small for no LJs
else
    Lmax = Nb;              %max LJ cutoff - set very small for no LJs
    % This will cause buggy periodic boundaries if
    % set too low relative to the maximum
    % interaction potential length scale
end

ftol = Kstiff*b/10000;  %equilibration force tolerance

L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# DEFINE LENGTH & FORCE SCALES'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Defining length scales'; Assembled = [Assembled,L];
L = ['\nvariable    Lx       equal ',num2str(Lx,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    Ly       equal ',num2str(Ly,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    Lz       equal ',num2str(Lz,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    b       equal ',num2str(b,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    N       equal ',num2str(N_Kuhn,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    Nb      equal ',num2str(Nb,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    Nb2     equal ',num2str(Nb2,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    r0      equal ',num2str(r0,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    Kstiff  equal ',num2str(Kstiff,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    ftol    equal ',num2str(ftol,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    Lmax    equal ',num2str(Lmax,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Defining max stretch'; Assembled = [Assembled,L];
L = ['\nvariable    lammax  equal ',num2str(lambda,'%.2f')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Define Thermodynamic properties
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# DEFINE LENGTH & FORCE SCALES'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
L = ['\nvariable    T       equal ',num2str(T,'%.0f')]; Assembled = [Assembled,L];
% L = ['\nvariable    damp    equal ',num2str(damp_eff,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    damp_th equal ',num2str(damp_th,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    damp_st equal ',num2str(damp_st,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Define repulsive properties, if any properties
lambda_r0 = 0.5;            % for soft core LJ potential
sig_r0 = b/(2^(1/6));       % backbone LJ length scale
eps_r0 = kbT;               % backbone LJ energy parameter
sig_b = b/(2^(1/6));        % backbone LJ length scale
eps_b = kbT;                % backbone LJ energy parameter
lambda_b = lambda_r0;       % for soft core LJ potential
rc = 2*sig_r0;              % repulsive cutoff lengthscale

% Turn LJ off for the comparison study only
if TurnOffLJ==0
    eps_r0 = 0;
    eps_b = 0;
    rc = b;
end

L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# DEFINE REPULSION (IF ANY)'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
% L = ['\nvariable    sig_r0  equal ',num2str(sig_r0,'%.2e')]; Assembled = [Assembled,L];
% L = ['\nvariable    eps_r0  equal ',num2str(eps_r0,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    sig_b  equal ',num2str(sig_r0,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    eps_b  equal ',num2str(eps_r0,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    lambda_r0 equal ',num2str(lambda_r0,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    sig_b   equal ',num2str(sig_b,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    eps_b   equal ',num2str(eps_b,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    lambda_b equal ',num2str(lambda_b,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    r_c     equal ',num2str(rc,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Define dynamic bond properties
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# DYNAMIC BOND PARAMETERS'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Kinetic rates'; Assembled = [Assembled,L];
L = ['\nvariable    ka  equal ',num2str(ka0,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    kd  equal ',num2str(kd0,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Attachment cutoff distance'; Assembled = [Assembled,L];
L = ['\nvariable    Ldyn    equal ',num2str(Ldyn,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];


% Define all pertinent timescales, output intervals, etc.
[dt_eq,dt_load,dt_relax,...
    neq,nload,nrelax,nmax,ttot,~,...
    i_dyn_eq,i_dyn_load,i_dyn_relax,...
    N_data,...
    iout_eq,iout_load,iout_relax,...
    ithermo_eq,~,~,...
    lamdot,taukd,tauka,tauW,teq,tload,trelax] =...
    DefineInputTimescales(tau0,Weissenberg,kd0,ka0,lambda,dt,...
    BeadSpringOrMeso,Controls);

% % if omega==1e-4
% %     nload = nload*4; % Test a case of running a much longer simulation
% %     nmax = nmax*4;
% % end
% disp(nmax*dt)

if N_Kuhn==36 && phi==0.5
    disp(N_Kuhn)
    disp(phi)
    disp(Weissenberg)
    disp(eq_time_factor)
    disp(neq*dt)
    disp(nload*dt)
    disp(nrelax*dt)
end

N_extra_bonds = 2;
N_extra_special = 1000;
if BeadSpringOrMeso==0
    N_extra_neighbors = 2000;
else
    if phi>0.4
        N_extra_neighbors = 10000;
    else
        N_extra_neighbors = 2500;
    end
end
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# TIMESCALES & LOADING RATES'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Weissenberg number'; Assembled = [Assembled,L];
L = ['\nvariable    Weiss   equal ',num2str(Weissenberg,'%.1f')]; Assembled = [Assembled,L];
L = ['\nvariable    lamdot  equal ',num2str(lamdot,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Timescales'; Assembled = [Assembled,L];
L = ['\nvariable    taukd   equal ',num2str(taukd,'%.1f')]; Assembled = [Assembled,L];
L = ['\nvariable    tauka   equal ',num2str(tauka,'%.1f')]; Assembled = [Assembled,L];
L = ['\nvariable    tauW    equal ',num2str(tauW,'%.1f')]; Assembled = [Assembled,L];
L = ['\nvariable    taueq equal ',num2str(teq,'%.1f')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Timestep'; Assembled = [Assembled,L];
L = ['\nvariable    dt        equal ',num2str(dt,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    dt_eq     equal ',num2str(dt_eq,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    dt_load   equal ',num2str(dt_load,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    dt_relax  equal ',num2str(dt_relax,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Time for each portion of the simulation'; Assembled = [Assembled,L];
L = ['\nvariable    teq   equal ',num2str(teq,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    tload   equal ',num2str(tload,'%.2e')]; Assembled = [Assembled,L];
L = ['\nvariable    trelax  equal ',num2str(trelax,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Number of timesteps for each sequence'; Assembled = [Assembled,L];
L = ['\nvariable    neq   equal ',num2str(neq,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    nload   equal ',num2str(nload,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    nrelax  equal ',num2str(nrelax,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    nmax    equal ',num2str(nmax,'%.0f')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Output frequencies'; Assembled = [Assembled,L];
L = ['\nvariable    iout_eq     equal ',num2str(iout_eq,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    iout_load     equal ',num2str(iout_load,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    iout_relax    equal ',num2str(iout_relax,'%.0f')]; Assembled = [Assembled,L];
L = ['\nvariable    ithermo     equal ',num2str(ithermo_eq,'%.0f')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['For W = ',num2str(Weissenberg),' with an equilibration time factor of ',num2str(eq_time_factor),','])
disp(['the total simulated time will be ',num2str(ttot,'%.2e'),' s,']);
nmaxout = addComma(nmax);
disp(['and the total number of steps will be ',nmaxout,'.']);

% Define output frequency
disp(['The total number of output data points will be ',num2str(N_data)])
disp(['and the cmd prompt will provide read out every ',num2str(ithermo_eq),' steps.'])
diary('Run time variables.txt')


% Running parameters
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# SIMULATION SPECIFICS'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# For inter-processor communication'; Assembled = [Assembled,L];
L = ['\nneighbor    ',num2str(Lmax,'%.2e'),' bin']; Assembled = [Assembled,L];
L = ['\ncomm_modify cutoff  ',num2str(Lmax,'%.2e')]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Initialize system
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n# INITIALIZE THE SYSTEM'; Assembled = [Assembled,L];
L = '\n####################################'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

L = '\n# Define atom style'; Assembled = [Assembled,L];
L = '\natom_style   bond'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Define bond styles
if BeadSpringOrMeso==0 || (BeadSpringOrMeso==1 && LinearOrLangevin==0)
    L = '\n# Define bond type'; Assembled = [Assembled,L];
    if BeadSpringOrMeso==1
        L = '\nbond_style    hybrid  harmonic nonlinear'; Assembled = [Assembled,L];
    elseif BeadSpringOrMeso==2
        L = '\nbond_style    harmonic'; Assembled = [Assembled,L];
    elseif BeadSpringOrMeso==0
        L = '\nbond_style    nonlinear'; Assembled = [Assembled,L];
    end
    L = '\n'; Assembled = [Assembled,L];
else
    L = '\nbond_style    hybrid pade nonlinear'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
end

if Controls.PrepareForStampede==1
    InputFolder1 = [WSL_path,'&'];
    InputFolder2 = partial_topology_filename;
    L = '\n# Import the topology'; Assembled = [Assembled,L];
    L = ['\nread_data ',InputFolder1]; Assembled = [Assembled,L];
    L = ['\n',InputFolder2,'&']; Assembled = [Assembled,L];
    L = ['\n extra/bond/per/atom ',num2str(N_extra_bonds,'%.0f'),...
        ' extra/special/per/atom ',num2str(N_extra_special,'%.0f')];
    Assembled = [Assembled,L];
else
    InputFolder1 = [WSL_path,current_folder,'&'];
    InputFolder2 = full_topology_filename;
    L = '\n# Import the topology'; Assembled = [Assembled,L];
    L = ['\nread_data ',InputFolder1]; Assembled = [Assembled,L];
    L = ['\n',InputFolder2,'&']; Assembled = [Assembled,L];
    L = ['\n extra/bond/per/atom ',num2str(N_extra_bonds,'%.0f'),...
        ' extra/special/per/atom ',num2str(N_extra_special,'%.0f')];
    Assembled = [Assembled,L];
end
L = ['\nneigh_modify    one ',num2str(N_extra_neighbors)];
Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Define bond styles
L = '\n# Define repulsive potentials between each atom'; Assembled = [Assembled,L];
L = '\npair_style  lj/cut/soft $(2.0) $(1.0) $(v_r_c)'; Assembled = [Assembled,L];
L = '\npair_coeff  1 1	${eps_b}  ${sig_b}  ${lambda_b}'; Assembled = [Assembled,L];
L = '\npair_coeff  2 2	${eps_b}  ${sig_b}  ${lambda_b}'; Assembled = [Assembled,L];
if BeadSpringOrMeso==0
    L = '\npair_coeff  3 3 ${eps_b}  ${sig_b}  ${lambda_b}'; Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

L = '\n# Label the groups'; Assembled = [Assembled,L];
L = '\ngroup	backbones type 1'; Assembled = [Assembled,L];
L = '\ngroup	stickers type 2'; Assembled = [Assembled,L];
if BeadSpringOrMeso==0
    L = '\ngroup	chains type 3'; Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];
L = '\nlabelmap	atom	1 bckbn		2 stckr'; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Generate random seed
rand_seed = randi(9999);
L = '\n# Set Brownian integration'; Assembled = [Assembled,L];
L = ['\nfix      part_th    backbones   brownian ${T} ',num2str(rand_seed),...
    ' gamma_t ${damp_th}']; Assembled = [Assembled,L];
L = ['\nfix      part_st    stickers    brownian ${T} ',num2str(rand_seed),...
    ' gamma_t ${damp_st}']; Assembled = [Assembled,L];
if BeadSpringOrMeso==0
    L = ['\nfix      part_ch    chains    brownian ${T} ',num2str(rand_seed),...
        ' gamma_t ${damp_st}']; Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% Define output files and run conditions

if Controls.RunLargeDeformation==1
    % Initial Equilibration
    L = '\n####################################'; Assembled = [Assembled,L];
    L = '\n# REACH STEADY STATE'; Assembled = [Assembled,L];
    L = '\n####################################'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    L = '\ntimestep ${dt_eq}'; Assembled = [Assembled,L]; % Set timestep to dt_eq
    if TurnOnDynamics==1
        L = ['\nfix		dynamic stickers bond/dynamic ',num2str(i_dyn_eq,'%.0f'),...
            ' 2 2 ${ka} ${kd} ${Ldyn} maxbond 1'];
        Assembled = [Assembled,L];
    end
    L = '\n'; Assembled = [Assembled,L];

    L = '\n# Compute bond information'; Assembled = [Assembled,L];
    L = '\ncompute 	Pairs all property/local btype batom1 batom2'; Assembled = [Assembled,L];
    L = '\ncompute 	Bonds all bond/local dx dy dz dist force'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    L = '\n# Generate .dump files'; Assembled = [Assembled,L];
    L = ['\ndump atomsdump_eq all custom ${iout_eq} ',OutputAtom,' id x y z vx vy vz type mol'];
    Assembled = [Assembled,L];
    L = ['\ndump bondsdump_eq all local ${iout_eq} ',OutputBond,' index c_Pairs[*] c_Bonds[*]'];
    Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    L = '\n# Outputs in command window'; Assembled = [Assembled,L];
    L = '\nthermo_style    custom step time dt bonds atoms etotal'; Assembled = [Assembled,L];
    L = '\nthermo          ${ithermo}'; Assembled = [Assembled,L];
    L = '\nthermo_modify	flush yes'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    L = '\nvariable	lam equal 1'; Assembled = [Assembled,L];

    L = '\n# Compute virial stress on each atom and sum throughout system'; Assembled = [Assembled,L];
    L = '\ncompute  sigma all pressure NULL bond'; Assembled = [Assembled,L];
    L = '\ncompute  sigma_pr all pressure NULL pair'; Assembled = [Assembled,L];
    L = '\nvariable	S1 equal c_sigma[1]'; Assembled = [Assembled,L];
    L = '\nvariable	S2 equal c_sigma[2]'; Assembled = [Assembled,L];
    L = '\nvariable	S3 equal c_sigma[3]'; Assembled = [Assembled,L];
    L = '\nvariable	S12 equal c_sigma[4]'; Assembled = [Assembled,L];
    L = '\nvariable	S23 equal c_sigma[5]'; Assembled = [Assembled,L];
    L = '\nvariable	S32 equal c_sigma[6]'; Assembled = [Assembled,L];
    L = '\nvariable	Sb11 equal c_sigma_pr[1]'; Assembled = [Assembled,L];
    L = '\nvariable	Sb22 equal c_sigma_pr[2]'; Assembled = [Assembled,L];
    L = '\nvariable	Sb33 equal c_sigma_pr[3]'; Assembled = [Assembled,L];
    L = '\nvariable	Sb12 equal c_sigma_pr[4]'; Assembled = [Assembled,L];
    L = '\nvariable	Sb23 equal c_sigma_pr[5]'; Assembled = [Assembled,L];
    L = '\nvariable	Sb32 equal c_sigma_pr[6]'; Assembled = [Assembled,L];
    L = '\nvariable	teq equal ${dt_eq}*elapsed'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    L = ['\nfix		printS all print ${iout_eq} "${teq} ${lam} ${S1} ${S2} ${S3}" file ',...
        stress_eq_filename,' screen no']; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    % Run initial equilibration
    L = '\n####################################'; Assembled = [Assembled,L];
    L = '\n# RUN INITIAL EQUILIBRATION'; Assembled = [Assembled,L];
    L = '\n####################################'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
    L = '\nrun  ${neq}'; Assembled = [Assembled,L]; % Run neq timesteps
    L = '\n'; Assembled = [Assembled,L];


    % Apply deformation
    L = '\n'; Assembled = [Assembled,L];
    L = '\n####################################'; Assembled = [Assembled,L];
    L = '\n# APPLY DEFORMATION'; Assembled = [Assembled,L];
    L = '\n####################################'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
    L = '\ntimestep ${dt_load}'; Assembled = [Assembled,L]; % Set timestep to dt_load
    if TurnOnDynamics==1
        L = ['\nfix		dynamic stickers bond/dynamic ',num2str(i_dyn_load,'%.0f'),...
            ' 2 2 ${ka} ${kd} ${Ldyn} maxbond 1'];
        Assembled = [Assembled,L];
    end
    L = ['\nfix		load_uniaxial all deform 1 x trate ${lamdot} y trate ',...
        '-$(v_lamdot*0.5) z trate -$(v_lamdot*0.5)']; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    L = '\n# Output deformation stress data'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
    L = '\nvariable	lam equal exp(${lamdot}*${dt_load}*elapsed)'; Assembled = [Assembled,L];
    L = '\nvariable	tl equal ${dt_load}*elapsed'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
    L = ['\nfix		printS all print ${iout_load} "${tl} ${lam} ${S1} ${S2} ${S3}" file ',...
        stress_ldg_filename,' screen no']; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    % Run loading phase
    L = '\n'; Assembled = [Assembled,L];
    L = '\n####################################'; Assembled = [Assembled,L];
    L = '\n# RUN LOADING'; Assembled = [Assembled,L];
    L = '\n####################################'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
    L = '\nrun  ${nload}'; Assembled = [Assembled,L]; % Run nload timesteps
    L = '\n'; Assembled = [Assembled,L];

elseif Controls.RunOscillatory==1
%     % Initial Equilibration
%     L = '\n####################################'; Assembled = [Assembled,L];
%     L = '\n# REACH STEADY STATE'; Assembled = [Assembled,L];
%     L = '\n####################################'; Assembled = [Assembled,L];
%     L = '\n'; Assembled = [Assembled,L];
% 
%     L = '\ntimestep ${dt_eq}'; Assembled = [Assembled,L]; % Set timestep to dt_eq
%     if TurnOnDynamics==1
%         L = ['\nfix		dynamic stickers bond/dynamic ',num2str(i_dyn_eq,'%.0f'),...
%             ' 2 2 ${ka} ${kd} ${Ldyn} maxbond 1'];
%         Assembled = [Assembled,L];
%     end
%     L = '\n'; Assembled = [Assembled,L];
% 
%     L = '\n# Compute bond information'; Assembled = [Assembled,L];
%     L = '\ncompute 	Pairs all property/local btype batom1 batom2'; Assembled = [Assembled,L];
%     L = '\ncompute 	Bonds all bond/local dx dy dz dist force'; Assembled = [Assembled,L];
%     L = '\n'; Assembled = [Assembled,L];
% 
%     L = '\n# Generate .dump files'; Assembled = [Assembled,L];
%     L = ['\ndump atomsdump_eq all custom ${iout_load} ',OutputAtom,' id x y z vx vy vz type mol'];
%     Assembled = [Assembled,L];
%     L = ['\ndump bondsdump_eq all local ${iout_load} ',OutputBond,' index c_Pairs[*] c_Bonds[*]'];
%     Assembled = [Assembled,L];
%     L = '\n'; Assembled = [Assembled,L];
% 
%     L = '\n# Outputs in command window'; Assembled = [Assembled,L];
%     L = '\nthermo_style    custom step time dt bonds atoms etotal'; Assembled = [Assembled,L];
%     L = '\nthermo          ${ithermo}'; Assembled = [Assembled,L];
%     L = '\nthermo_modify	flush yes'; Assembled = [Assembled,L];
%     L = '\n'; Assembled = [Assembled,L];
% 
%     L = '\nvariable	lam equal 1'; Assembled = [Assembled,L];
% 
%     L = '\n# Compute virial stress on each atom and sum throughout system'; Assembled = [Assembled,L];
%     L = '\ncompute  sigma all pressure NULL bond'; Assembled = [Assembled,L];
%     L = '\ncompute  sigma_pr all pressure NULL pair'; Assembled = [Assembled,L];
%     L = '\nvariable	S1 equal c_sigma[1]'; Assembled = [Assembled,L];
%     L = '\nvariable	S2 equal c_sigma[2]'; Assembled = [Assembled,L];
%     L = '\nvariable	S3 equal c_sigma[3]'; Assembled = [Assembled,L];
%     L = '\nvariable	S12 equal c_sigma[4]'; Assembled = [Assembled,L];
%     L = '\nvariable	S23 equal c_sigma[5]'; Assembled = [Assembled,L];
%     L = '\nvariable	S32 equal c_sigma[6]'; Assembled = [Assembled,L];
%     L = '\nvariable	Sb11 equal c_sigma_pr[1]'; Assembled = [Assembled,L];
%     L = '\nvariable	Sb22 equal c_sigma_pr[2]'; Assembled = [Assembled,L];
%     L = '\nvariable	Sb33 equal c_sigma_pr[3]'; Assembled = [Assembled,L];
%     L = '\nvariable	Sb12 equal c_sigma_pr[4]'; Assembled = [Assembled,L];
%     L = '\nvariable	Sb23 equal c_sigma_pr[5]'; Assembled = [Assembled,L];
%     L = '\nvariable	Sb32 equal c_sigma_pr[6]'; Assembled = [Assembled,L];
%     L = '\nvariable	teq equal ${dt_eq}*elapsed'; Assembled = [Assembled,L];
%     L = '\n'; Assembled = [Assembled,L];

    % Initial Equilibration
    L = '\n####################################'; Assembled = [Assembled,L];
    L = '\n# REACH STEADY STATE'; Assembled = [Assembled,L];
    L = '\n####################################'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    if TurnOnDynamics==1
        L = ['\nfix		dynamic stickers bond/dynamic ',num2str(i_dyn_eq,'%.0f'),...
            ' 2 2 ${ka} ${kd} ${Ldyn} maxbond 1'];
        Assembled = [Assembled,L];
    end

    L = '\n# Compute bond information'; Assembled = [Assembled,L];
    L = '\ncompute 	Pairs all property/local btype batom1 batom2'; Assembled = [Assembled,L];
    L = '\ncompute 	Bonds all bond/local dx dy dz dist force'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    L = '\n# Outputs files'; Assembled = [Assembled,L];
    L = ['\ndump mydump all custom ${iout_eq} ',OutputAtom,' id x y z vx vy vz type mol'];
    Assembled = [Assembled,L];
    L = ['\ndump bondsdump all local ${iout_eq} ',OutputBond,' index c_Pairs[*] c_Bonds[*]'];
    Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    L = '\n# Outputs in command window'; Assembled = [Assembled,L];
    L = '\nthermo_style    custom step time dt bonds atoms etotal'; Assembled = [Assembled,L];
    L = '\nthermo          ${ithermo}'; Assembled = [Assembled,L];
    L = '\nthermo_modify	flush yes'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    L = '\n# Set timestep and minimization conditions'; Assembled = [Assembled,L];
    L = '\ntimestep ${dt}'; Assembled = [Assembled,L];

    L = '\n# Compute virial stress on each atom and sum throughout system';
    Assembled = [Assembled,L];
    L = '\nvariable	lam equal 1'; Assembled = [Assembled,L];


%     L = '\n# Output equilibration data'; Assembled = [Assembled,L];
%     L = '\ncompute  sigma all pressure NULL bond'; Assembled = [Assembled,L];
%     L = '\ncompute  sigma_pr all pressure NULL pair'; Assembled = [Assembled,L];
%     L = '\nvariable	S1 equal c_sigma[1]'; Assembled = [Assembled,L];
%     L = '\nvariable	S2 equal c_sigma[2]'; Assembled = [Assembled,L];
%     L = '\nvariable	S3 equal c_sigma[3]'; Assembled = [Assembled,L];
%     L = '\nvariable	S12 equal c_sigma[4]'; Assembled = [Assembled,L];
%     L = '\nvariable	S23 equal c_sigma[5]'; Assembled = [Assembled,L];
%     L = '\nvariable	S32 equal c_sigma[6]'; Assembled = [Assembled,L];
%     L = '\nvariable	Sb11 equal c_sigma_pr[1]'; Assembled = [Assembled,L];
%     L = '\nvariable	Sb22 equal c_sigma_pr[2]'; Assembled = [Assembled,L];
%     L = '\nvariable	Sb33 equal c_sigma_pr[3]'; Assembled = [Assembled,L];
%     L = '\nvariable	Sb12 equal c_sigma_pr[4]'; Assembled = [Assembled,L];
%     L = '\nvariable	Sb23 equal c_sigma_pr[5]'; Assembled = [Assembled,L];
%     L = '\nvariable	Sb32 equal c_sigma_pr[6]'; Assembled = [Assembled,L];
%     L = '\nvariable	teq equal ${dt}*elapsed'; Assembled = [Assembled,L];
%     L = '\n'; Assembled = [Assembled,L];
% 
%     L = ['\nfix		printS all print ${iout_eq} "${teq} ${lam} ${S1} ${S2} ${S3} ${S12} ${S23} ${S32}" file ',...
%         stress_eq_filename,' screen no']; Assembled = [Assembled,L];
%     L = '\n'; Assembled = [Assembled,L];

    L = '\n# Run initial equilibration'; Assembled = [Assembled,L];
    L = '\nrun  ${neq}'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    
    % Apply deformation
    L = '\n'; Assembled = [Assembled,L];
    L = '\n####################################'; Assembled = [Assembled,L];
    L = '\n# APPLY DEFORMATION'; Assembled = [Assembled,L];
    L = '\n####################################'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
    L = '\ntimestep ${dt_load}'; Assembled = [Assembled,L]; % Set timestep to dt_load
    if TurnOnDynamics==1
        L = ['\nfix		dynamic stickers bond/dynamic ',num2str(i_dyn_load,'%.0f'),...
            ' 2 2 ${ka} ${kd} ${Ldyn} maxbond 1'];
        Assembled = [Assembled,L];
    end

    L = '\n# Expoonentially increasing frequency'; Assembled = [Assembled,L];
    L = ['\nvariable    eps_0       equal ',num2str(lambda-1,'%.3f')]; Assembled = [Assembled,L];
    L = ['\nvariable    omega_min   equal ',num2str(omega_min/tau0,'%.2e')]; Assembled = [Assembled,L];
    L = ['\nvariable    omega_max   equal ',num2str(omega_max/tau0,'%.2e')]; Assembled = [Assembled,L];
    L = ['\nvariable    omega       equal ',num2str(omega/tau0,'%.2e')]; Assembled = [Assembled,L];
%     L = ['\nvariable    tf          equal ',num2str(tf*tau0,'%.2e')]; Assembled = [Assembled,L];
    L = ['\nvariable    tf          equal ',num2str(tf*tau0/omega,'%.2e')]; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

    L = '\n# Define important intermediate constants and variables'; Assembled = [Assembled,L];
    L = '\nvariable    K1          equal 2*PI*${omega_min}'; Assembled = [Assembled,L];
    L = '\nvariable    K2          equal log(${omega_max}/${omega_min})'; Assembled = [Assembled,L];
    L = '\nvariable    K3          equal exp(${K2}*time/${tf})'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];

%     L = '\n# Define frequency, sinusoidal stretches/strains, and sinusoidal strain rates'; Assembled = [Assembled,L];
%     L = '\nvariable     omega       equal ${omega_min}*${K3}'; Assembled = [Assembled,L];
%     L = '\nvariable     eps11       equal ${eps_0}*sin(${K1}*${K3}*time)'; Assembled = [Assembled,L];
%     L = '\nvariable     deps11dt    equal ${eps_0}/${tf}*${K1}*(${K2}*time+${tf})*${K3}*cos(${K1}*${K3}*time)'; Assembled = [Assembled,L];
%     L = '\nvariable     deps22dt    equal ${deps11dt}*(-1/2)'; Assembled = [Assembled,L];
%     L = '\nvariable     lam1        equal ${eps11}+1'; Assembled = [Assembled,L];
%     L = '\nvariable     lam2        equal 1/sqrt(${lam1})'; Assembled = [Assembled,L];  
%     L = '\nvariable     eps22       equal ${lam2}-1'; Assembled = [Assembled,L];  
%     L = '\nvariable     dxt         equal ${eps11}*${Lx}'; Assembled = [Assembled,L];  
%     L = '\nvariable     dyt         equal ${eps22}*${Ly}'; Assembled = [Assembled,L];  
%     L = '\nvariable     dzt         equal ${eps22}*${Lz}'; Assembled = [Assembled,L]; 
%     L = '\n'; Assembled = [Assembled,L];
%     
%     L = '\n# Check values'; Assembled = [Assembled,L];
%     L = '\nprint "Current value of omega: ${omega}"'; Assembled = [Assembled,L];    
%     L = '\nprint "Current value of lam1: ${lam1}"'; Assembled = [Assembled,L];    
%     L = '\nprint "Current value of lam2: ${lam2}"'; Assembled = [Assembled,L];    
%     L = '\nprint "Current value of eps11: ${eps11}"'; Assembled = [Assembled,L];    
%     L = '\nprint "Current value of eps22: ${eps22}"'; Assembled = [Assembled,L];   
%     L = '\nprint "Current value of deps11t: ${deps11dt}"'; Assembled = [Assembled,L];   
%     L = '\nprint "Current value of deps22t: ${deps22dt}"'; Assembled = [Assembled,L];   
%     L = '\nprint "Current value of dxt: ${dxt}"'; Assembled = [Assembled,L];    
%     L = '\nprint "Current value of dyt: ${dyt}"'; Assembled = [Assembled,L];    
%     L = '\nprint "Current value of dzt: ${dzt}"'; Assembled = [Assembled,L];    
%     L = '\nprint "Current value of time: ${dt_load}*step"'; Assembled = [Assembled,L];    
%     L = '\n'; Assembled = [Assembled,L];

    L = '\n# Apply sinusoidal pure shear displacement with variable frequency'; Assembled = [Assembled,L];
    L = '\nvariable displace equal "v_Lx/2 * v_eps_0 * sin(2*PI * v_omega * (step*dt - v_teq))"'; Assembled = [Assembled,L];
    L = '\nvariable rate equal "v_Lx/2 * v_eps_0 * 2*PI*v_omega * cos(2*PI * v_omega * (step*dt - v_teq))"'; Assembled = [Assembled,L];
%     L = '\nvariable displace2 equal "-v_Lx/2 * v_eps_0 * sin(2*PI * v_omega * step*dt)"'; Assembled = [Assembled,L];
%     L = '\nvariable rate2 equal "-v_Lx/2 * v_eps_0 * 2*PI*v_omega * cos(2*PI * v_omega * step*dt)"'; Assembled = [Assembled,L];
    L = '\nvariable displace2 equal "-v_displace"'; Assembled = [Assembled,L];
    L = '\nvariable rate2 equal "-v_rate"'; Assembled = [Assembled,L];
    
    L = '\nfix      def_x all deform 1 x variable v_displace v_rate remap x'; Assembled = [Assembled,L];
    L = '\nfix      def_y all deform 1 y variable v_displace2 v_rate2 remap x'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];


    L = '\n# Output deformation stress data'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
    L = '\nvariable	lam equal exp(${lamdot}*${dt_load}*elapsed)'; Assembled = [Assembled,L];
    L = '\nvariable	tl equal ${dt_load}*elapsed'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];


    % Run loading phase
    L = '\n'; Assembled = [Assembled,L];
    L = '\n####################################'; Assembled = [Assembled,L];
    L = '\n# RUN LOADING'; Assembled = [Assembled,L];
    L = '\n####################################'; Assembled = [Assembled,L];
    L = '\n'; Assembled = [Assembled,L];
    L = '\nrun  ${nload}'; Assembled = [Assembled,L]; % Run nload timesteps
    L = '\n'; Assembled = [Assembled,L];

else % if regular stress relaxation benchmarking study
    if BeadSpringOrMeso==0 % Only one .dump file
        % Initial Equilibration
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n# REACH STEADY STATE'; Assembled = [Assembled,L];
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        if TurnOnDynamics==1
            L = ['\nfix		dynamic stickers bond/dynamic ',num2str(i_dyn_eq,'%.0f'),...
                ' 2 2 ${ka} ${kd} ${Ldyn} maxbond 1'];
            Assembled = [Assembled,L];
        end

        L = '\n# Compute bond information'; Assembled = [Assembled,L];
        L = '\ncompute 	Pairs all property/local btype batom1 batom2'; Assembled = [Assembled,L];
        L = '\ncompute 	Bonds all bond/local dx dy dz dist force'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = '\n# Outputs files'; Assembled = [Assembled,L];
        L = ['\ndump mydump all custom ${iout_eq} ',OutputAtom,' id x y z vx vy vz type mol'];
        Assembled = [Assembled,L];
        L = ['\ndump bondsdump all local ${iout_eq} ',OutputBond,' index c_Pairs[*] c_Bonds[*]'];
        Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = '\n# Outputs in command window'; Assembled = [Assembled,L];
        L = '\nthermo_style    custom step time dt bonds atoms etotal'; Assembled = [Assembled,L];
        L = '\nthermo          ${ithermo}'; Assembled = [Assembled,L];
        L = '\nthermo_modify	flush yes'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = '\n# Set timestep and minimization conditions'; Assembled = [Assembled,L];
        L = '\ntimestep ${dt}'; Assembled = [Assembled,L];

        L = '\n# Compute virial stress on each atom and sum throughout system';
        Assembled = [Assembled,L];
        L = '\nvariable	lam equal 1'; Assembled = [Assembled,L];


        L = '\n# Output equilibration data'; Assembled = [Assembled,L];
        L = '\ncompute  sigma all pressure NULL bond'; Assembled = [Assembled,L];
        L = '\ncompute  sigma_pr all pressure NULL pair'; Assembled = [Assembled,L];
        L = '\nvariable	S1 equal c_sigma[1]'; Assembled = [Assembled,L];
        L = '\nvariable	S2 equal c_sigma[2]'; Assembled = [Assembled,L];
        L = '\nvariable	S3 equal c_sigma[3]'; Assembled = [Assembled,L];
        L = '\nvariable	S12 equal c_sigma[4]'; Assembled = [Assembled,L];
        L = '\nvariable	S23 equal c_sigma[5]'; Assembled = [Assembled,L];
        L = '\nvariable	S32 equal c_sigma[6]'; Assembled = [Assembled,L];
        L = '\nvariable	Sb11 equal c_sigma_pr[1]'; Assembled = [Assembled,L];
        L = '\nvariable	Sb22 equal c_sigma_pr[2]'; Assembled = [Assembled,L];
        L = '\nvariable	Sb33 equal c_sigma_pr[3]'; Assembled = [Assembled,L];
        L = '\nvariable	Sb12 equal c_sigma_pr[4]'; Assembled = [Assembled,L];
        L = '\nvariable	Sb23 equal c_sigma_pr[5]'; Assembled = [Assembled,L];
        L = '\nvariable	Sb32 equal c_sigma_pr[6]'; Assembled = [Assembled,L];
        L = '\nvariable	teq equal ${dt}*elapsed'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = ['\nfix		printS all print ${iout_eq} "${teq} ${lam} ${S1} ${S2} ${S3} ${S12} ${S23} ${S32}" file ',...
            stress_eq_filename,' screen no']; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = '\n# Run initial equilibration'; Assembled = [Assembled,L];
        L = '\nrun  ${neq}'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        % Apply deformation
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n# APPLY DEFORMATION'; Assembled = [Assembled,L];
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];
        L = ['\nfix		load_uniaxial all deform 1 x trate ${lamdot} y trate ',...
            '-$(v_lamdot*0.5) z trate -$(v_lamdot*0.5)']; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = '\n# Output deformation data'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];
        L = '\nvariable	lam equal exp(${lamdot}*${dt}*elapsed)'; Assembled = [Assembled,L];
        L = '\nvariable	tl equal ${dt}*elapsed'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = ['\nfix		printS all print ${iout_eq} "${tl} ${lam} ${S1} ${S2} ${S3}" file ',...
            stress_ldg_filename,' screen no']; Assembled = [Assembled,L];

        L = '\n'; Assembled = [Assembled,L];

        L = '\n# Run applied deformation'; Assembled = [Assembled,L];
        L = '\nrun  ${nload}'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        % Relax system
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n# HOLD AT CONSTANT DEFORM. - RELAX'; Assembled = [Assembled,L];
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];
        L = '\nunfix	load_uniaxial'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];
        L = '\n# Output relaxation data'; Assembled = [Assembled,L];
        L = '\nvariable	tr equal ${dt}*elapsed'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = ['\nfix		printS all print ${iout_eq} "${tr} ${lam} ${S1} ${S2} ${S3}" file '...
            stress_rlx_filename,' screen no']; Assembled = [Assembled,L];

        L = '\n'; Assembled = [Assembled,L];

        L = '\n# Run relaxation'; Assembled = [Assembled,L];
        L = '\nrun		${nrelax}'; Assembled = [Assembled,L];

    elseif BeadSpringOrMeso==1
        % Initial Equilibration
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n# REACH STEADY STATE'; Assembled = [Assembled,L];
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = '\ntimestep ${dt_eq}'; Assembled = [Assembled,L]; % Set timestep to dt_eq
        if TurnOnDynamics==1
            L = ['\nfix		dynamic stickers bond/dynamic ',num2str(i_dyn_eq,'%.0f'),...
                ' 2 2 ${ka} ${kd} ${Ldyn} maxbond 1'];
            Assembled = [Assembled,L];
        end
        L = '\n'; Assembled = [Assembled,L];

        L = '\n# Compute bond information'; Assembled = [Assembled,L];
        L = '\ncompute 	Pairs all property/local btype batom1 batom2'; Assembled = [Assembled,L];
        L = '\ncompute 	Bonds all bond/local dx dy dz dist force'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = '\n# Eq. outputs files'; Assembled = [Assembled,L];
        L = ['\ndump atomsdump_eq all custom ${iout_eq} ',OutputAtom_eq,' id x y z vx vy vz type mol'];
        Assembled = [Assembled,L];
        L = ['\ndump bondsdump_eq all local ${iout_eq} ',OutputBond_eq,' index c_Pairs[*] c_Bonds[*]'];
        Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = '\n# Outputs in command window'; Assembled = [Assembled,L];
        L = '\nthermo_style    custom step time dt bonds atoms etotal'; Assembled = [Assembled,L];
        L = '\nthermo          ${ithermo}'; Assembled = [Assembled,L];
        L = '\nthermo_modify	flush yes'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = '\nvariable	lam equal 1'; Assembled = [Assembled,L];

        L = '\n# Compute virial stress on each atom and sum throughout system'; Assembled = [Assembled,L];
        L = '\ncompute  sigma all pressure NULL bond'; Assembled = [Assembled,L];
        L = '\ncompute  sigma_pr all pressure NULL pair'; Assembled = [Assembled,L];
        L = '\nvariable	S1 equal c_sigma[1]'; Assembled = [Assembled,L];
        L = '\nvariable	S2 equal c_sigma[2]'; Assembled = [Assembled,L];
        L = '\nvariable	S3 equal c_sigma[3]'; Assembled = [Assembled,L];
        L = '\nvariable	S12 equal c_sigma[4]'; Assembled = [Assembled,L];
        L = '\nvariable	S23 equal c_sigma[5]'; Assembled = [Assembled,L];
        L = '\nvariable	S32 equal c_sigma[6]'; Assembled = [Assembled,L];
        L = '\nvariable	Sb11 equal c_sigma_pr[1]'; Assembled = [Assembled,L];
        L = '\nvariable	Sb22 equal c_sigma_pr[2]'; Assembled = [Assembled,L];
        L = '\nvariable	Sb33 equal c_sigma_pr[3]'; Assembled = [Assembled,L];
        L = '\nvariable	Sb12 equal c_sigma_pr[4]'; Assembled = [Assembled,L];
        L = '\nvariable	Sb23 equal c_sigma_pr[5]'; Assembled = [Assembled,L];
        L = '\nvariable	Sb32 equal c_sigma_pr[6]'; Assembled = [Assembled,L];
        L = '\nvariable	teq equal ${dt_eq}*elapsed'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = ['\nfix		printS all print ${iout_eq} "${teq} ${lam} ${S1} ${S2} ${S3}" file ',...
            stress_eq_filename,' screen no']; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        % Run initial equilibration
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n# RUN INITIAL EQUILIBRATION'; Assembled = [Assembled,L];
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];
        L = '\nrun  ${neq}'; Assembled = [Assembled,L]; % Run neq timesteps
        L = '\n'; Assembled = [Assembled,L];


        % Apply deformation
        L = '\n'; Assembled = [Assembled,L];
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n# APPLY DEFORMATION'; Assembled = [Assembled,L];
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];
        L = '\ntimestep ${dt_load}'; Assembled = [Assembled,L]; % Set timestep to dt_load
        if TurnOnDynamics==1
            L = ['\nfix		dynamic stickers bond/dynamic ',num2str(i_dyn_load,'%.0f'),...
                ' 2 2 ${ka} ${kd} ${Ldyn} maxbond 1'];
            Assembled = [Assembled,L];
        end
        L = ['\nfix		load_uniaxial all deform 1 x trate ${lamdot} y trate ',...
            '-$(v_lamdot*0.5) z trate -$(v_lamdot*0.5)']; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = '\n# Deformation outputs files'; Assembled = [Assembled,L];
        L = ['\ndump atomsdump_ld all custom ${iout_load} ',OutputAtom_ld,' id x y z vx vy vz type mol'];
        Assembled = [Assembled,L];
        L = ['\ndump bondsdump_ld all local ${iout_load} ',OutputBond_ld,' index c_Pairs[*] c_Bonds[*]'];
        Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = '\n# Output deformation data'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];
        L = '\nvariable	lam equal exp(${lamdot}*${dt_load}*elapsed)'; Assembled = [Assembled,L];
        L = '\nvariable	tl equal ${dt_load}*elapsed'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];
        L = ['\nfix		printS all print ${iout_load} "${tl} ${lam} ${S1} ${S2} ${S3}" file ',...
            stress_ldg_filename,' screen no']; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        % Run loading phase
        L = '\n'; Assembled = [Assembled,L];
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n# RUN LOADING'; Assembled = [Assembled,L];
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];
        L = '\nrun  ${nload}'; Assembled = [Assembled,L]; % Run nload timesteps
        L = '\n'; Assembled = [Assembled,L];

        % Relax system
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n# HOLD AT CONSTANT DEFORM. - RELAX'; Assembled = [Assembled,L];
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];
        L = '\ntimestep ${dt_relax}'; Assembled = [Assembled,L]; % Set timestep to dt_eq
        if TurnOnDynamics==1
            L = ['\nfix		dynamic stickers bond/dynamic ',num2str(i_dyn_relax,'%.0f'),...
                ' 2 2 ${ka} ${kd} ${Ldyn} maxbond 1'];
            Assembled = [Assembled,L];
        end
        L = '\n'; Assembled = [Assembled,L];
        L = '\nunfix	load_uniaxial'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = '\n# Relaxation outputs files'; Assembled = [Assembled,L];
        L = ['\ndump atomsdump_rlx all custom ${iout_relax} ',OutputAtom_rlx,' id x y z vx vy vz type mol'];
        Assembled = [Assembled,L];
        L = ['\ndump bondssdump_rlx all local ${iout_relax} ',OutputBond_rlx,' index c_Pairs[*] c_Bonds[*]'];
        Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        L = '\n# Output relaxation data'; Assembled = [Assembled,L];
        L = '\nvariable	tr equal ${dt_relax}*elapsed'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];
        L = ['\nfix		printS all print ${iout_relax} "${tr} ${lam} ${S1} ${S2} ${S3}" file '...
            stress_rlx_filename,' screen no']; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];

        % Run relaxation phase
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n# RUN RELAXATION'; Assembled = [Assembled,L];
        L = '\n####################################'; Assembled = [Assembled,L];
        L = '\n'; Assembled = [Assembled,L];
        L = '\nrun  ${nrelax}'; Assembled = [Assembled,L]; % Run nrelax timesteps
        L = '\n'; Assembled = [Assembled,L];
    end
end

% Save the .in file
fid = fopen(input_filename,'wt');
fprintf(fid,Assembled);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function numOut = addComma(numIn)
jf=java.text.DecimalFormat; % comma for thousands, three decimal places
numOut= char(jf.format(numIn)); % omit "char" if you want a string out
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeMovies(N)

global dims Np Nt Ns TopologyFileTag

% Plot periodic and unwrapped domains
figure(4)
clf
PlotNetwork

figure(5)
clf
PlotNetworkUnwrapped;

if dims==2
    if ~isfolder('Charts')
        mkdir('Charts')
    end
    figure(4)
    PlotName = ['Charts/Np_',num2str(Np),'.Nt_',num2str(Nt),...
        '.Nb_',num2str(N),'.z_',num2str(Ns),'.png'];
    saveas(gcf,PlotName);

    figure(5)
    PlotName = ['Charts/Unwrapped.Np_',num2str(Np),'.Nt_',num2str(Nt),...
        '.Nb_',num2str(N),'.z_',num2str(Ns),'.png'];
    saveas(gcf,PlotName);
end

if dims==3
    figure(4)
    Frames = 120;
    ang0 = 45;
    anglerange = (linspace(ang0,ang0+360,Frames))';
    for i=1:Frames
        ang = anglerange(i);
        set(gca,'CameraViewAngleMode','Manual')
        axis equal
        view(ang,30)
        F1 = getframe(gcf);
        %         writeVideo(v,F1)
        mov(i) = F1;
    end
    if ~isfolder('Gifs')
        mkdir('Gifs')
    end
    GifName = ['Gifs/',TopologyFileTag,'.gif'];
    movie2gif(mov, GifName, 'DelayTime', 0,'LoopCount',Inf)
end

figure(5)
%         if ~isfolder('Movies')
%             mkdir('Movies')
%         end
%     MovieName = ['Movies/Unwrapped.Np_',num2str(Np),'.Nt_',num2str(N),...
%        '.Nb_',num2str(Nt),'.z_',num2str(Ns)];
%     if MakeMovie==1
%         v = VideoWriter(MovieName);
%         v.FrameRate = 10;
%         open(v)
%     end
ang0 = 45;
anglerange = (linspace(ang0,ang0+360,Frames))';
for i=1:Frames
    ang = anglerange(i);
    set(gca,'CameraViewAngleMode','Manual')
    axis equal
    view(ang,30)
    F1 = getframe(gcf);
    %         writeVideo(v,F1)
    mov(i) = F1;
end
%     close(v)
if ~isfolder('Gifs')
    mkdir('Gifs')
end
GifName = ['Gifs/Unwrapped_',TopologyFileTag,'.gif'];
movie2gif(mov, GifName, 'DelayTime', 0,'LoopCount',Inf)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveMatlabData

global Positions PosUnwrapped dims Lx Ly Lz Np Nt ConnectionsSP...
    N_atoms N_bonds Ns initial_matlab_filename

% if ~isfile(initial_matlab_filename)
    % Define N_atoms and N_bonds
    N_atoms = size(Positions,1);
    Froms = zeros(N_atoms*(size(ConnectionsSP,2)-1),1);
    Tos = Froms;
    for i=1:size(ConnectionsSP,2)-1
        rng = ((N_atoms*i-N_atoms+1):(N_atoms*i))';
        Froms(rng) = ConnectionsSP(:,i);
        Tos(rng) = ConnectionsSP(:,end);
    end
    Pairs = [Froms Tos];
    UniquePairs = sort(Pairs,2);
    UniquePairs(UniquePairs(:,1)==0,:) = [];
    UniquePairs = unique(UniquePairs,'rows');
    N_bonds = size(UniquePairs,1);

    % Append coordinates
    N_primary = Np*Nt;      % No tether sites
    N_stickers = Np*Nt*Ns;  % No sticker sites
    tether_rng = 1:N_primary;       % indices of the tether nodes
    sticker_rng = N_primary+1:N_primary+N_stickers; % indices of the sticker nodes
    rest_rng = N_primary+N_stickers+1:size(Positions,1); % indices of the in-between backbone beads

    tags = PosUnwrapped(:,dims+1);
    mol_tags = PosUnwrapped(:,dims+2);
    x = PosUnwrapped(:,1); y = PosUnwrapped(:,2);
    if dims==2
        z = zeros(size(x));
    else
        z = PosUnwrapped(:,3);
    end

    atom_type = zeros(size(x));
    for i=1:N_atoms
        if ismember(i,tether_rng) %% || ismember(i,rest_rng)
            atom_type(i) = 1;
        elseif ismember(i,sticker_rng)
            atom_type(i) = 2;
        elseif ismember(i,rest_rng)
            atom_type(i) = 3;
        else
            warning('Node is not a member of any population of beads')
        end
    end
    InitialPositions = [tags mol_tags atom_type x y z];

    % Append bonds
    tags = (1:N_bonds)';
    atom_1 = UniquePairs(:,1);
    atom_2 = UniquePairs(:,2);
    bond_types = ones(size(tags));
    InitialBonds = [tags bond_types atom_1 atom_2];

    % Define structure
    MatlabData.Positions = InitialPositions;
    MatlabData.Bonds = InitialBonds;
    MatlabData.Boundaries = [Lx Ly Lz];
    MatlabData.PosUnwrapped = PosUnwrapped;

    save(initial_matlab_filename,'-struct','MatlabData','-v7.3')
% else
%     dat = load(initial_matlab_filename,'-mat');
%     InitialPositions = dat.Positions;
%     InitialBonds = dat.Bonds;
%     bounds = dat.Boundaries;
%     Lx = bounds(1);
%     Ly = bounds(2);
%     Lz = bounds(3);
%     PosUnwrapped = dat.PosUnwrapped ;
% end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveTopologyDataLAMMPS(N_Kuhn,kb,T,b,Lx,Ly,Lz)

global FileTag PosUnwrapped dims Np Nt ConnectionsSP... %Positions
    N_atoms N_bonds full_topology_filename Ns BeadSpringOrMeso LinearOrLangevin

% Define N_atoms and N_bonds
% kbT = energy_conversion*kbT;
kbT = kb*T;
mass_mer = 1;

N_atoms = size(PosUnwrapped,1);
Froms = zeros(N_atoms*(size(ConnectionsSP,2)-1),1);
Tos = Froms;
for i=1:size(ConnectionsSP,2)-1
    rng = ((N_atoms*i-N_atoms+1):(N_atoms*i))';
    Froms(rng) = ConnectionsSP(:,i);
    Tos(rng) = ConnectionsSP(:,end);
end
Pairs = [Froms Tos];
UniquePairs = sort(Pairs,2);
UniquePairs(UniquePairs(:,1)==0,:) = [];
UniquePairs = unique(UniquePairs,'rows');
N_bonds = size(UniquePairs,1);

% Formats output data as follows for LAMMPS:

% LAMMPS Description           (1st line of file)
%
% 100 atoms         (this must be the 3rd line, 1st 2 lines are ignored)
% 95 bonds                (# of bonds to be simulated)
% 50 angles               (include these lines even if number = 0)
% 30 dihedrals
% 20 impropers
%
% 5 atom types           (# of nonbond atom types)
% 10 bond types          (# of bond types = sets of bond coefficients)
% 18 angle types
% 20 dihedral types      (do not include a bond,angle,dihedral,improper type
% 2 improper types             line if number of bonds,angles,etc is 0)
%
% -0.5 0.5 xlo xhi       (for periodic systems this is box size,
% -0.5 0.5 ylo yhi        for non-periodic it is min/max extent of atoms)
% -0.5 0.5 zlo zhi       (do not include this line for 2-d simulations)
%
% Masses
%
%   1 mass
%   ...
%   N mass                           (N = # of atom types)
%
% Nonbond Coeffs
%
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of atom types)
%
% Bond Coeffs
%
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of bond types)
%
% Angle Coeffs
%
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of angle types)
%
% Dihedral Coeffs
%
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of dihedral types)
%
% Improper Coeffs
%
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of improper types)
%
% BondBond Coeffs
%
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of angle types)
%
% BondAngle Coeffs
%
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of angle types)
%
% MiddleBondTorsion Coeffs
%
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of dihedral types)
%
% EndBondTorsion Coeffs
%
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of dihedral types)
%
% AngleTorsion Coeffs
%
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of dihedral types)
%
% AngleAngleTorsion Coeffs
%
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of dihedral types)
%
% BondBond13 Coeffs
%
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of dihedral types)
%
% AngleAngle Coeffs
%
%   1 coeff1 coeff2 ...
%   ...
%   N coeff1 coeff2 ...              (N = # of improper types)
%
% Atoms
%
%   1 molecule-tag atom-type q x y z nx ny nz  (nx,ny,nz are optional -
%   ...                                    see "true flag" input command)
%   ...
%   N molecule-tag atom-type q x y z nx ny nz  (N = # of atoms)
%
% Velocities
%
%   1 vx vy vz
%   ...
%   ...
%   N vx vy vz                        (N = # of atoms)
%
% Bonds
%
%   1 bond-type atom-1 atom-2
%   ...
%   N bond-type atom-1 atom-2         (N = # of bonds)
%
% Angles
%
%   1 angle-type atom-1 atom-2 atom-3  (atom-2 is the center atom in angle)
%   ...
%   N angle-type atom-1 atom-2 atom-3  (N = # of angles)
%
% Dihedrals
%
%   1 dihedral-type atom-1 atom-2 atom-3 atom-4  (atoms 2-3 form central bond)
%   ...
%   N dihedral-type atom-1 atom-2 atom-3 atom-4  (N = # of dihedrals)
%
% Impropers
%
%   1 improper-type atom-1 atom-2 atom-3 atom-4  (atom-2 is central atom)
%   ...
%   N improper-type atom-1 atom-2 atom-3 atom-4  (N = # of impropers)

% Initialize string
Assembled = [];
if BeadSpringOrMeso==0
    N_atom_types = 3;   % Stickers, tethers, and backbone
else
    N_atom_types = 2;   % Stickers and backbone
end
N_bond_types = 2;       % Dynamic vs stable

% Header
L = ['Initial Topology',FileTag]; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
L = ['\n',num2str(N_atoms),' atoms']; Assembled = [Assembled,L];
L = ['\n',num2str(N_bonds),' bonds']; Assembled = [Assembled,L];
L = ['\n',num2str(0),' angles']; Assembled = [Assembled,L];
L = ['\n',num2str(0),' dihedrals']; Assembled = [Assembled,L];
L = ['\n',num2str(0),' impropers']; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
L = ['\n',num2str(N_atom_types),' atom types']; Assembled = [Assembled,L];
L = ['\n',num2str(N_bond_types),' bond types']; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
L = ['\n',num2str(-Lx/2),' ',num2str(Lx/2),' xlo xhi']; Assembled = [Assembled,L];
L = ['\n',num2str(-Ly/2),' ',num2str(Ly/2),' ylo yhi']; Assembled = [Assembled,L];
L = ['\n',num2str(-Lz/2),' ',num2str(Lz/2),' zlo zhi']; Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];

% Append masses
mass_backbone = mass_mer;
mass_sticker = mass_mer;
L = '\nMasses';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
for i=1:N_atom_types
    type = i;
    if type==1 || type==3
        mass = mass_backbone;
    elseif type==2
        mass = mass_sticker;
    end
    L = ['\n',num2str(type),' ',num2str(mass)];
    Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% Append bond coefficients - for Langevin spring potentials
L = '\nBond Coeffs';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
soften_coeff = 1.05;    % This slight softening proves necessary for
% numerical stability of the mesocale model
for i=1:N_bond_types
    type = i;
    if BeadSpringOrMeso==1 % Mesoscale
        if LinearOrLangevin==0
            if type==1
                coeff1 = 3/2*kbT/(N_Kuhn*b^2);
                coeff2 = 0;
                L = ['\n',num2str(type),' harmonic ',...
                    num2str(coeff1),' ',num2str(coeff2)];
            elseif type==2
                coeff1 = 100*kbT;
                coeff2 = b;
                coeff3 = soften_coeff*b;
                L = ['\n',num2str(type),' nonlinear ',...
                    num2str(coeff1),' ',num2str(coeff2),' ',num2str(coeff3)];
            end
        elseif LinearOrLangevin==1
            if type==1
                coeff1 = b;
                coeff2 = N_Kuhn;
                L = ['\n',num2str(type),' pade ',...
                    num2str(coeff1),' ',num2str(coeff2)];
            elseif type==2
                coeff1 = 100*kbT;
                coeff2 = b;
                coeff3 = soften_coeff*b;
                L = ['\n',num2str(type),' nonlinear ',...
                    num2str(coeff1),' ',num2str(coeff2),' ',num2str(coeff3)];
            end
        end
    elseif BeadSpringOrMeso==0 % Bead-spring
        if type==1
            coeff1 = 800*kbT;
            coeff2 = b;
            coeff3 = b;
        elseif type==2
            coeff1 = 100*kbT;
            coeff2 = b;
            coeff3 = b;
        end
        L = ['\n',num2str(type),' ',num2str(coeff1),...
            ' ',num2str(coeff2),...
            ' ',num2str(coeff3)];
    elseif BeadSpringOrMeso==2 % Scaling Theory
        coeff1 = 3/2*kbT/(2*N_Kuhn*b^2);
        coeff2 = 0;
        L = ['\n',num2str(type),' ',...
            num2str(coeff1),' ',num2str(coeff2)];
    end
    Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% Append coordinates
L = '\nAtoms';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
N_primary = Np*Nt;      % No tether sites
N_stickers = Np*Nt*Ns;  % No sticker sites
tether_rng = 1:N_primary;       % indices of the tether nodes
sticker_rng = N_primary+1:N_primary+N_stickers; % indices of the sticker nodes
if BeadSpringOrMeso==0
    rest_rng = N_primary+N_stickers+1:size(PosUnwrapped,1); % indices of the in-between backbone beads
else
    rest_rng = [];
end
for i=1:N_atoms
    tag = PosUnwrapped(i,dims+1);
    mol_tag = PosUnwrapped(i,dims+2);
    x = PosUnwrapped(i,1); y = PosUnwrapped(i,2);
    if dims==2
        z = 0;
    else
        z = PosUnwrapped(i,3);
    end
    if ismember(i,tether_rng) %% || ismember(i,rest_rng)
        atom_type = 1;
    elseif ismember(i,sticker_rng)
        atom_type = 2;
    elseif ismember(i,rest_rng)
        atom_type = 3;
    else
        warning('Node is not a member of any population of beads')
    end
    L = ['\n',num2str(tag),' ',num2str(mol_tag),' ',num2str(atom_type),' ',...
        num2str(x),' ',num2str(y),' ',num2str(z)];
    Assembled = [Assembled,L];
end
L = '\n    '; Assembled = [Assembled,L];

% Append initial velocities
L = '\nVelocities';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
for i=1:N_atoms
    tag = PosUnwrapped(i,dims+1);
    vx = 0; vy = 0; vz = 0;
    L = ['\n',num2str(tag),' ',num2str(vx),' ',num2str(vy),' ',num2str(vz)];
    Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

% Append bonds
L = '\nBonds';  Assembled = [Assembled,L];
L = '\n'; Assembled = [Assembled,L];
for i=1:N_bonds
    tag = i;
    atom_1 = UniquePairs(i,1); atom_2 = UniquePairs(i,2);
    bond_type = 1;
    L = ['\n',num2str(tag),' ',num2str(bond_type),' ',...
        num2str(atom_1),' ',num2str(atom_2)];
    Assembled = [Assembled,L];
end
L = '\n'; Assembled = [Assembled,L];

fid = fopen(full_topology_filename,'wt');
fprintf(fid,Assembled);

CheckFig=0;
if CheckFig==1
    PlotNetworkUnwrapped;
    for i=1:size(PosUnwrapped,1)
        pos = PosUnwrapped(i,1:dims);
        tag = PosUnwrapped(i,dims+1);
        scatter3(pos(1),pos(2),pos(3),'g','filled')
        text(pos(1),pos(2),pos(3),['-',num2str(tag)],'FontSize',8)
    end
    xlim([-Inf Inf])
    ylim([-Inf Inf])
    zlim([-Inf Inf])
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UnwrapPeriodic(Override)

global PosUnwrapped...
    ShowFiguresDebug...
    Positions Positions0 Corners0 Corners...
    X_SP_Unwrapped Y_SP_Unwrapped Z_SP_Unwrapped X_SP Y_SP Z_SP...
    X_REP Y_REP Z_REP ConnectionsSP ConnectionsREP...
    Rx_SP Ry_SP Rz_SP Rx_REP Ry_REP Rz_REP sigmaREP CutoffPartREP...
    initial_node_pos_filename

dat = load(initial_node_pos_filename,'-mat');
if ~isfield(dat,'PosUnwrapped') || Override

    UnwrapBeadSpring;

    % Save data
    InitialTetherPos.PosUnwrapped = PosUnwrapped;
    InitialTetherPos.Positions = Positions;
    InitialTetherPos.Positions0 = Positions0;
    InitialTetherPos.Corners0 = Corners0;
    InitialTetherPos.Corners = Corners;
    InitialTetherPos.X_SP_Unwrapped = X_SP_Unwrapped;
    InitialTetherPos.Y_SP_Unwrapped = Y_SP_Unwrapped;
    InitialTetherPos.Z_SP_Unwrapped = Z_SP_Unwrapped;
    InitialTetherPos.X_SP = X_SP;
    InitialTetherPos.Y_SP = Y_SP;
    InitialTetherPos.Z_SP = Z_SP;
    InitialTetherPos.X_REP = X_REP;
    InitialTetherPos.Y_REP = Y_REP;
    InitialTetherPos.Z_REP = Z_REP;
    InitialTetherPos.ConnectionsSP = ConnectionsSP;
    InitialTetherPos.ConnectionsREP = ConnectionsREP;
    InitialTetherPos.Rx_SP = Rx_SP;
    InitialTetherPos.Ry_SP = Ry_SP;
    InitialTetherPos.Rz_SP = Rz_SP;
    InitialTetherPos.Rx_REP = Rx_REP;
    InitialTetherPos.Ry_REP = Ry_REP;
    InitialTetherPos.Rz_REP = Rz_REP;
    InitialTetherPos.sigmaREP = sigmaREP;
    InitialTetherPos.CutoffPartREP = CutoffPartREP;

    save(initial_node_pos_filename,'-struct','InitialTetherPos','-v7.3')
else
    Positions = dat.Positions;
    Positions0 = dat.Positions0;
    PosUnwrapped = dat.PosUnwrapped;
    X_SP_Unwrapped = dat.X_SP_Unwrapped;
    Y_SP_Unwrapped = dat.Y_SP_Unwrapped;
    Z_SP_Unwrapped = dat.Z_SP_Unwrapped;
    X_SP = dat.X_SP;
    Y_SP = dat.Y_SP;
    Z_SP = dat.Z_SP;
    X_REP = dat.X_REP;
    Y_REP = dat.Y_REP;
    Z_REP = dat.Z_REP;
    ConnectionsSP = dat.ConnectionsSP;
    ConnectionsREP = dat.ConnectionsREP;
    Rx_SP = dat.Rx_SP;
    Ry_SP = dat.Ry_SP;
    Rz_SP = dat.Rz_SP;
    Rx_REP = dat.Rx_REP;
    Ry_REP = dat.Ry_REP;
    Rz_REP = dat.Rz_REP;
    sigmaREP = dat.sigmaREP;
    CutoffPartREP = dat.CutoffPartREP;
end


if ShowFiguresDebug==1
    figure(5)
    clf
    PlotNetworkUnwrapped
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UnwrapBeadSpring

global dims Np...
    PosUnwrapped...
    Positions  X_SP_Unwrapped Y_SP_Unwrapped Z_SP_Unwrapped...
    X_SP Y_SP Z_SP...
    ConnectionsSP Rx_SP Ry_SP Rz_SP


PosUnwrapped = zeros(size(Positions));
PosUnwrapped(:,dims+1:end) = Positions(:,dims+1:end);

X_SP_Unwrapped = zeros(size(X_SP));
Y_SP_Unwrapped = zeros(size(Y_SP));
Z_SP_Unwrapped = zeros(size(Z_SP));

CheckFigs = 0;
if CheckFigs==1
    figure(100); clf; hold on
    PlotNetwork
end

color_rng = (linspace(0,1,Np))';
colors = [color_rng flipud(color_rng) flipud(color_rng)];
for mol=1:Np
    temp_positions = Positions;
    temp_positions(temp_positions(:,dims+2)~=mol,:) = [];
    members = temp_positions(:,dims+1);
    members = sort(members);
    temp_connections = ConnectionsSP(members,:);

    froms = temp_connections(:,1:end-1); froms = froms(:);
    tos = repmat(temp_connections(:,end),size(temp_connections,2)-1,1);
    pairs = [froms tos]; pairs(pairs(:,1)==0,:) = [];
    pairs = sort(pairs,2); pairs = unique(pairs,'rows');

    already_positioned = pairs(1,1);        % Establish a positiioned particle
    p0 = already_positioned;
    PosUnwrapped(p0,:) = Positions(p0,:);
    i = 0;
    while ~isempty(pairs)
        % ID an already positioned member of the pairs list that still has
        % neighbors
        refs = intersect(pairs(:),already_positioned);
        ref = refs(1);

        % Grab first pair in pairs list containing the ref particle
        pair_indx = min([find(pairs(:,1)==ref,1,'first');...
            find(pairs(:,2)==ref,1,'first')]);
        pair = pairs(pair_indx,:);
        new = pair(pair~=ref);

        % Establish reference position
        pos_ref = PosUnwrapped(ref,1:dims); % reference position

        % Position new node's location based on connectivity vectors
        col = find(ConnectionsSP(new,:)==ref);
        dx = Rx_SP(new,col);
        dy = Ry_SP(new,col);
        dX = [dx dy];
        if dims==3
            dz = Rz_SP(new,col);
            dX = [dx dy dz];
        end
        PosUnwrapped(new,1:dims) = pos_ref + dX;

        if CheckFigs==1
            if i==1
                s = scatter3(PosUnwrapped(ref,1),...
                    PosUnwrapped(ref,2),...
                    PosUnwrapped(ref,3),'filled');
                s.SizeData = 8;
                s.MarkerFaceColor = colors(mol,:);
            end
            s = scatter3(PosUnwrapped(new,1),...
                PosUnwrapped(new,2),...
                PosUnwrapped(new,3),'filled');
            s.SizeData = 8;
            s.MarkerFaceColor = colors(mol,:);
        end

        % Remove pair that has been addressed
        already_positioned = cat(1,already_positioned,new);
        pairs(pair_indx,:) = [];
    end
end
[X_SP_Unwrapped,Y_SP_Unwrapped,Z_SP_Unwrapped] =...
    ReIndexXYUnwrapped(ConnectionsSP,...
    X_SP_Unwrapped,Y_SP_Unwrapped,Z_SP_Unwrapped,Rx_SP);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PolymerizeBranches

global dims N Lx Ly Lz Positions Ns ConnectionsSP...
    X_SP Y_SP Z_SP Rx_SP Ry_SP Rz_SP

if dims==2
    rho = N/(Lx*Ly);     %Should match that of FCC based on initiation
else
    rho = N/(Lx*Ly*Lz);     %Should match that of FCC based on initiation
end
xi = rho^(-1/dims);
radius = 0.5*xi;    %Set less than mesh size for stability

% Add [N x Ns] rows to Positions, ConnectionsSP, X-Z, Rx-Rz
% Recall: Positions = [X, Y, Z, xl number, molecule number, mer number]
N_stickers = N*Ns;
AddRange = (N+1:N+N_stickers)';
Positions(AddRange,:) = 0;
Positions(AddRange,dims+1) = AddRange; %Add xl numbers
ConnectionsSP(AddRange,:) = 0;
ConnectionsSP(AddRange,end) = AddRange;
X_SP(AddRange,:) = 0; Y_SP(AddRange,:) = 0; Z_SP(AddRange,:) = 0;
Rx_SP(AddRange,:) = 0; Ry_SP(AddRange,:) = 0; Rz_SP(AddRange,:) = 0;

%% ADD TO RX_SP RY_SP RZ_SP CONNECTIONSSP CONNECTIONSREP

% Generate random stickers at 1/2 distance of existing mesh size
indx = N;
PosTemp = Positions;
PosTemp = sortrows(PosTemp,dims+2); %Sort by molecule number
PosTemp(PosTemp(:,dims+2)==0,:) = [];
wb = waitbar(0,'Polymerizing the branches');
for i=1:N
    waitbar(i/N,wb)
    mol = PosTemp(i,dims+2);
    for j=1:Ns
        indx = indx+1;
        Positions(indx,dims+2) = mol;
        point = PosTemp(i,1:dims);
        theta = 2*pi*rand();
        %         varphi = pi*rand();

        if dims==2
            dx = radius*cos(theta);
            dy = radius*sin(theta);
            newpoint = [point(1)+dx,point(2)+dy];
        else
            varphimin = -pi/2;
            varphimax = pi/2;
            varphi = varphimin+rand()*(varphimax-varphimin);
            radius_projected = radius*cos(varphi);
            dx = radius_projected*cos(theta);
            dy = radius_projected*sin(theta);
            dz = radius*sin(varphi);
            newpoint = [point(1)+dx,...
                point(2)+dy,...
                point(3)+dz];
        end

        %Store positions data, ConnectionsSP, X_SP, Y_SP, Z_SP, Rx_SP, Ry_SP, Rz_SP

        % Updat Positions
        Positions(indx,1:dims) = newpoint;
        sisters = Positions((Positions(:,dims+2)==mol),dims+1:end);
        xl_nos = sisters(:,end);
        Positions(indx,end) = max(xl_nos)+1;

        %Update  Connections SP
        P1 = Positions(indx,dims+1);
        P2 = PosTemp(i,dims+1);
        disp = Positions(P2,1:dims)-Positions(P1,1:dims);
        %Index the neighbors
        FirstZero_n = find(ConnectionsSP(P1,:)==0,1,'first');
        ConnectionsSP(P1,FirstZero_n) = P2;
        Rx_SP(P1,FirstZero_n) = -disp(1);
        Ry_SP(P1,FirstZero_n) = -disp(2);
        if dims==3
            Rz_SP(P1,FirstZero_n) = -disp(3);
        end
        FirstZero_j = find(ConnectionsSP(P2,:)==0,1,'first');
        ConnectionsSP(P2,FirstZero_j) = P1;
        Rx_SP(P2,FirstZero_j) = +disp(1);
        Ry_SP(P2,FirstZero_j) = +disp(2);
        if dims==3
            Rz_SP(P2,FirstZero_j) = +disp(3);
        end
        [X_SP,Y_SP,Z_SP] = ReIndexXY(ConnectionsSP,X_SP,Y_SP,Z_SP,Rx_SP);

        %Store
        Check = 0;
        if Check==1
            figure(100)
            hold on
            PlotNetwork;
            if dims==2
                s = scatter(point(1),point(2),'k','filled');
                s.SizeData = 6;
                s = scatter(newpoint(1),newpoint(2),'r','filled');
                s.SizeData = 6;
            elseif dims==3
                s = scatter3(point(1),point(2),point(3),'k','filled');
                s.SizeData = 6;
                s = scatter3(newpoint(1),newpoint(2),newpoint(3),'r','filled');
                s.SizeData = 6;
            end
        end
    end

end
close(wb)

%Update ConnectionsREP
FindBoundaryNodes;
ReplicateBounds;
ConnectREP;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PolymerizeChains

global dims N Nt Positions Rx_SP Ry_SP Rz_SP ConnectionsSP...
    MaxZ X_SP Y_SP Z_SP...

MaxZ = 4;
ConnectionsSP = zeros(size(Positions,1),MaxZ+1);
ConnectionsSP(:,end) = (1:N)';
Rx_SP = zeros(size(Positions,1),MaxZ);
Ry_SP = zeros(size(Positions,1),MaxZ);
Rz_SP = zeros(size(Positions,1),MaxZ);
X_SP = zeros(size(Rx_SP));
Y_SP = zeros(size(Rx_SP));
Z_SP = zeros(size(Rx_SP));
%Loop over all xl and initiate a new molecule every Nt
k = 0;
mol = 0;
Pool = (1:N)';
Positions(:,end+1:end+2) = zeros(size(Positions,1),2);
FindBoundaryNodes;      %Update boundary nodes since equilibration

wb = waitbar(0,'Polymerizing the crosslinks...');
while k<N
    k = k+1;
    PcntComplete = k/N;
    if ~mod(k,10)
        waitbar(PcntComplete,wb,'Polymerizing the crosslinks...')
    end

    if ~mod(k-1,Nt)   %If first xl in molecule - pick at random
        mol = mol+1;
        seq = 1;                %Position in molecule
        mn = 1;
        mx = size(Pool,1);
        indx = round((mx-mn).*rand() + mn);
        xl = Pool(indx);
    else %Find closest and link
        PoolArray = FindPeriodicPool(xl,Pool);

        % Remove already-bonded nodes from PoolArray
        All_xl = (1:N)';
        Removed = All_xl(~ismember(All_xl,Pool));
        PoolArray(ismember(PoolArray(:,1),Removed),:) = [];

        X1 = Positions(xl,1:dims).*ones(size(PoolArray,1),dims);
        X2 = PoolArray(:,2:dims+1);%Positions(PoolArray(:,1),2:dims+1);
        disp = X2-X1;
        dist = vecnorm(disp,2,2);
        SortMat = [dist PoolArray(:,1) disp];
        SortMat = sortrows(SortMat,1);  %sort by distance
        P1 = xl;                %First xl in pairing
        xl = SortMat(1,2);
        P2 = xl;                %Second xl in pairing
        seq = seq+1;
        disp_closest = SortMat(1,3:end);

        Check = 0;
        if Check==1
            if dims==2
                scatter(Positions(P1,1),Positions(P1,2),'r','filled')
                scatter(Positions(P2,1),Positions(P2,2),'c','filled')
            elseif dims==3
                scatter3(Positions(P1,1),Positions(P1,2),Positions(P1,3),'r','filled')
                scatter3(Positions(P2,1),Positions(P2,2),Positions(P2,3),'c','filled')
            end
        end
        %Index the neighbors
        FirstZero_n = find(ConnectionsSP(P1,:)==0,1,'first');
        ConnectionsSP(P1,FirstZero_n) = P2;
        Rx_SP(P1,FirstZero_n) = -disp_closest(1);
        Ry_SP(P1,FirstZero_n) = -disp_closest(2);
        if dims==3
            Rz_SP(P1,FirstZero_n) = -disp_closest(3);
        end
        FirstZero_j = find(ConnectionsSP(P2,:)==0,1,'first');
        ConnectionsSP(P2,FirstZero_j) = P1;
        Rx_SP(P2,FirstZero_j) = +disp_closest(1);
        Ry_SP(P2,FirstZero_j) = +disp_closest(2);
        if dims==3
            Rz_SP(P2,FirstZero_j) = +disp_closest(3);
        end
    end
    Positions(xl,end-1) = mol;  %Demarks which molecule this xl belongs to
    Positions(xl,end) = seq;    %Demarks which position in molecule xl occupies
    Pool(Pool==xl) = [];

    for j=1:MaxZ
        X_SP(ConnectionsSP(:,j)~=0,j) = Positions(ConnectionsSP(:,j)~=0,1);
        Y_SP(ConnectionsSP(:,j)~=0,j) = Positions(ConnectionsSP(:,j)~=0,2);
        if dims==3
            Z_SP(ConnectionsSP(:,j)~=0,j) = Positions(ConnectionsSP(:,j)~=0,3);
        end
    end
end

close(wb)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PoolArray = FindPeriodicPool(xl,Pool)

global dims AllBnds Lx Ly Lz Color Positions...
    LeftBnds RightBnds FrontBnds BackBnds BottomBnds TopBnds...
    BottomLeftBnds BottomRightBnds TopLeftBnds TopRightBnds...
    BottomFrontBnds BottomBackBnds TopFrontBnds TopBackBnds...
    FrontLeftBnds FrontRightBnds BackLeftBnds BackRightBnds...
    C1Bnds C2Bnds C3Bnds C4Bnds C5Bnds C6Bnds C7Bnds C8Bnds...

% ID if ismember of boundary nodes;
PoolArray = zeros(length(Pool),dims+1); %Allocates nodes into sortable array for distance check
PoolArray(:,1) = Pool;
PoolArray(:,2:end) = Positions(Pool,1:dims);
if ismember(xl,AllBnds) %Append neighboring ghost nodes as needed
    %Add ghost nodes to pool
    %Check which boundary(ies) its on and append ghost nodes and their positions
    if ismember(xl,LeftBnds)
        % Add opposite face
        OppBnd = RightBnds; Color = 'r';
        if dims==2
            Shift = [-Lx 0];
        elseif dims==3
            Shift = [-Lx 0 0];
        end
        PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        % If ismember of one of four edges on left, add opposite edge
        if ismember(xl,BottomLeftBnds)
            OppBnd = TopRightBnds; Color = 'g';
            if dims==2
                Shift = [-Lx -Ly];
            else
                Shift = [-Lx 0 -Lz];
            end
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
        if ismember(xl,TopLeftBnds)
            OppBnd = BottomRightBnds; Color = 'g';
            if dims==2
                Shift = [-Lx +Ly];
            else
                Shift = [-Lx 0 +Lz];
            end
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
        if ismember(xl,FrontLeftBnds)
            OppBnd = BackRightBnds; Color = 'c';
            Shift = [-Lx -Ly 0];
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
        if ismember(xl,BackLeftBnds)
            OppBnd = FrontRightBnds; Color = 'c';
            Shift = [-Lx +Ly 0];
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
        % If ismember of one of four corners on left, add opposite corner
        if ismember(xl,C1Bnds)
            OppBnd = C7Bnds; Color = 'm';
            Shift = [-Lx -Ly -Lz];
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
        if ismember(xl,C4Bnds)
            OppBnd = C6Bnds; Color = 'm';
            Shift = [-Lx +Ly -Lz];
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
        if ismember(xl,C8Bnds)
            OppBnd = C2Bnds; Color = 'm';
            Shift = [-Lx +Ly +Lz];
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
        if ismember(xl,C5Bnds)
            OppBnd = C3Bnds; Color = 'm';
            Shift = [-Lx -Ly +Lz];
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
    end

    if ismember(xl,RightBnds)
        % Add opposite face
        OppBnd = LeftBnds; Color = 'r';
        if dims==2
            Shift = [+Lx 0];
        elseif dims==3
            Shift = [+Lx 0 0];
        end
        PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        % If ismember of one of four edges on left, add opposite edge
        if ismember(xl,BottomRightBnds)
            OppBnd = TopLeftBnds; Color = 'g';
            if dims==2
                Shift = [+Lx -Ly];
            else
                Shift = [+Lx 0 -Lz];
            end
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
        if ismember(xl,TopRightBnds)
            OppBnd = BottomLeftBnds; Color = 'g';
            if dims==2
                Shift = [+Lx +Ly];
            else
                Shift = [+Lx 0 +Lz];
            end
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
        if ismember(xl,FrontRightBnds)
            OppBnd = BackLeftBnds; Color = 'c';
            Shift = [+Lx -Ly 0];
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
        if ismember(xl,BackRightBnds)
            OppBnd = FrontLeftBnds; Color = 'c';
            Shift = [+Lx +Ly 0];
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
        % If ismember of one of four corners on left, add opposite corner
        if ismember(xl,C2Bnds)
            OppBnd = C8Bnds; Color = 'm';
            Shift = [+Lx -Ly -Lz];
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
        if ismember(xl,C3Bnds)
            OppBnd = C5Bnds; Color = 'm';
            Shift = [+Lx +Ly -Lz];
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
        if ismember(xl,C7Bnds)
            OppBnd = C1Bnds; Color = 'm';
            Shift = [+Lx +Ly +Lz];
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
        if ismember(xl,C6Bnds)
            OppBnd = C4Bnds; Color = 'm';
            Shift = [+Lx -Ly +Lz];
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
    end

    if ismember(xl,BottomBnds)
        % Add opposite face
        OppBnd = TopBnds; Color = 'r';
        if dims==2
            Shift = [0 -Ly];
        elseif dims==3
            Shift = [0 0 -Lz];
        end
        PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        % If ismember of one of two remaining edges on bottom, add opposite edge
        if ismember(xl,BottomFrontBnds)
            OppBnd = TopBackBnds; Color = 'g';
            Shift = [0 -Ly -Lz];
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
        if ismember(xl,BottomBackBnds)
            OppBnd = TopFrontBnds; Color = 'g';
            Shift = [0 +Ly -Lz];
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
    end

    if ismember(xl,TopBnds)
        % Add opposite face
        OppBnd = BottomBnds; Color = 'r';
        if dims==2
            Shift = [0 +Ly];
        elseif dims==3
            Shift = [0 0 +Lz];
        end
        PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        % If ismember of one of two remaining edges on bottom, add opposite edge
        if ismember(xl,TopFrontBnds)
            OppBnd = BottomBackBnds; Color = 'g';
            Shift = [0 -Ly +Lz];
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
        if ismember(xl,TopBackBnds)
            OppBnd = BottomFrontBnds; Color = 'g';
            Shift = [0 +Ly +Lz];
            PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
        end
    end

    if ismember(xl,FrontBnds)
        % Add opposite face
        OppBnd = BackBnds; Color = 'r';
        Shift = [0 -Ly 0];
        PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
    end

    if ismember(xl,BackBnds)
        % Add opposite face
        OppBnd = FrontBnds; Color = 'r';
        Shift = [0 +Ly 0];
        PoolArray = AddOppBnds(PoolArray,OppBnd,Shift);
    end
end

if ismember(0,PoolArray(:,1))
    PoolArray(PoolArray(:,1)==0,:) = [];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PoolArray = AddOppBnds(PoolArray,OppBnd,Shift)

global dims Positions Color

Range = (length(PoolArray)+1:length(PoolArray)+length(OppBnd))';
PoolArray(Range,1) = OppBnd;
PoolArray(Range,2:dims+1) = Positions(OppBnd,1:dims) + Shift;

ShowFigDebug = 0;
if ShowFigDebug==1
    if dims==2
        s = scatter(Positions(OppBnd,1)+Shift(1),Positions(OppBnd,2)+Shift(2),...
            Color,'filled');
    elseif dims==3
        s = scatter3(Positions(OppBnd,1)+Shift(1),Positions(OppBnd,2)+Shift(2),...
            Positions(OppBnd,3)+Shift(3),Color,'filled');
    end
    s.SizeData = 6;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EquilibrateTheDomain

global ConnectionsSP ConnectionsREP X_SP X_REP Positions Positions0...
    Y_SP Y_REP Z_SP Z_REP Rx_SP Rx_REP Ry_SP Ry_REP Rz_SP Rz_REP...
    ResMax ResAvg NumericalDamper RepulsionRedefFreq dt dims...
    check_Kuhn_lengths b


F_NET = CalculateNetForces;

Damper = 10*NumericalDamper;
dX = Damper*F_NET;

Residuals = vecnorm(F_NET,2,2);

%Calculate Residuals
MaxResidual = max(max(Residuals));
AvgResidual = mean(mean(Residuals));
PrevRes = AvgResidual;

allowable_pcnt = 15;
if check_Kuhn_lengths
    norms_Kuhn = DefineKuhnLengths;
    min_b = min(norms_Kuhn);
    max_b = max(norms_Kuhn);

    % Allowable deviation from b
    min_allowable_b = b-allowable_pcnt/100*b;
    max_allowable_b = b+allowable_pcnt/100*b;

    if min_b<min_allowable_b || max_b>max_allowable_b
        check_b_condition = 1;
    else
        check_b_condition = 0;
    end

else
    check_b_condition = 0;
end

k=0;
wb = waitbar(0,'Equilibrating the domain...');
while (MaxResidual>ResMax && AvgResidual>ResAvg) || check_b_condition
    k = k+1;

    Positions0 = Positions;
    Positions(:,1:dims) = Positions(:,1:dims) + dX; %Update postions

    %Update Positions & Separation Vectors in Connection Matrices
    if ~isempty(X_SP)
        [X_SP,Y_SP,Z_SP,Rx_SP,Ry_SP,Rz_SP] =...
            UpdateNodePositions(ConnectionsSP,X_SP,Y_SP,Z_SP,Rx_SP,Ry_SP,Rz_SP);
    end
    if ~isempty(X_REP)
        [X_REP,Y_REP,Z_REP,Rx_REP,Ry_REP,Rz_REP] =...
            UpdateNodePositions(ConnectionsREP,X_REP,Y_REP,Z_REP,...
            Rx_REP,Ry_REP,Rz_REP);
    end
    %Updates Positions, X_SP, Y_SP, X_REP, Y_REP,...
    %Rx_SP, Ry_SP, Rx_REP and Ry_REP


    %Implement Eularian condition and edefine repulsive neighbors
    if ~mod(k,RepulsionRedefFreq)
        % Transport Nodes that exit frame to opossite side (only acts if BCs
        % are for Couette flow right now (i.e., parallel plate rheometer)
        ImplementEularianBoundaryConditions
        FindBoundaryNodes
        ReplicateBounds
        ConnectREP
    end

    F_NET = CalculateNetForces;

    dX = Damper*F_NET;

    Residuals = vecnorm(F_NET,2,2);

    %Calculate Residuals
    MaxResidual = max(max(Residuals));
    AvgResidual = mean(mean(Residuals));

    if mod(k,100)==0
        disp(['Residuals = ',num2str(AvgResidual)]);
        waitbar(ResAvg/AvgResidual,wb,'Equilibrating the domain...');

        if PrevRes<=AvgResidual %Feedback loop to increase damper as needed (5% at a time)
            Damper = Damper*0.5;
            visc = dt/Damper;
            disp(['Increase visc to ',num2str(visc)]);
        end
        PrevRes = AvgResidual;
    end

    if check_Kuhn_lengths
        norms_Kuhn = DefineKuhnLengths;
        min_b = min(norms_Kuhn);
        max_b = max(norms_Kuhn);

        % Allowable deviation from b
        min_allowable_b = b-allowable_pcnt/100*b;
        max_allowable_b = b+allowable_pcnt/100*b;

        if min_b<min_allowable_b || max_b>max_allowable_b
            check_b_condition = 1;
        else
            break  % If the b length criteria are met, stop equilibrating
        end
    else
        check_b_condition = 0;
    end


    if mod(k,5000)==0
        break
    end
end
close(wb)

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImplementEularianBoundaryConditions

%Reinserts nodes across boundaries if they step outside the frame
global Positions dims...
    ConnectionsSP X_SP Y_SP Z_SP Rx_SP...
    ConnectionsREP X_REP Y_REP Z_REP Rx_REP Corners

if dims==2
    CBL = Corners(1,:); CBR = Corners(2,:);
    CTR = Corners(3,:); CTL = Corners(4,:);
    Lx = CBR(1)-CBL(1); Ly = CTR(2)-CBR(2);
    X = Positions(:,1); Y = Positions(:,2);
    dYdX = CTR(2)-CTL(2); dXdY = CTR(1)-CBR(1);

    %Find particles out of left bound and reposition
    XofY_L = FindXofY(Y,CTL,CBL);
    if ~isempty(Positions(Positions(:,1)<XofY_L,:))
        Positions(Positions(:,1)<XofY_L,2) = Positions(Positions(:,1)<XofY_L,2) + dYdX;
        Positions(Positions(:,1)<XofY_L,1) = Positions(Positions(:,1)<XofY_L,1) + Lx;
    end

    %Find particles out of right bound and reposition
    XofY_R = FindXofY(Y,CTR,CBR);
    if ~isempty(Positions(Positions(:,1)>XofY_R,:))
        Positions(Positions(:,1)>XofY_R,2) = Positions(Positions(:,1)>XofY_R,2) - dYdX;
        Positions(Positions(:,1)>XofY_R,1) = Positions(Positions(:,1)>XofY_R,1) - Lx;
    end

    %Find particles out of top bound and reposition
    YofX_B = FindYofX(X,CBL,CBR);
    if ~isempty(Positions(Positions(:,2)<YofX_B,:))
        Positions(Positions(:,2)<YofX_B,1) = Positions(Positions(:,2)<YofX_B,1) + dXdY;
        Positions(Positions(:,2)<YofX_B,2) = Positions(Positions(:,2)<YofX_B,2) + Ly;
    end

    %Find particles out of top bound and reposition
    YofX_T = FindYofX(X,CTL,CTR);
    if ~isempty(Positions(Positions(:,2)>YofX_T,:))
        Positions(Positions(:,2)>YofX_T,1) = Positions(Positions(:,2)>YofX_T,1) - dXdY;
        Positions(Positions(:,2)>YofX_T,2) = Positions(Positions(:,2)>YofX_T,2) - Ly;
    end
elseif dims==3
    C1 = Corners(1,:);
    C2 = Corners(2,:);
    C4 = Corners(4,:);
    C5 = Corners(5,:);
    Lx = norm(C2-C1);
    Ly = norm(C1-C4);
    Lz = norm(C5-C1);
    X = Positions(:,1); Y = Positions(:,2); Z = Positions(:,3);
    dXdY = C4(1)-C1(1);
    dYdX = C2(2)-C1(2);
    dYdZ = C5(2)-C1(2);
    dZdY = C4(3)-C1(3);
    dZdX = C2(3)-C1(3);
    dXdZ = C5(1)-C1(1);

    %Find particles out of left bound and reposition
    v1 = C4-C1;
    v2 = C5-C1;
    C = C1;
    Out = FindOutOfBounds(v1,v2,C);

    %Find particles out of left bound and reposition
    v1 = C4-C1;
    v2 = C5-C1;
    C = C1;
    Out = FindOutOfBounds(v1,v2,C);
    Positions(Out,1) = Positions(Out,1) + Lx;
    Positions(Out,2) = Positions(Out,2) + dYdX;
    Positions(Out,3) = Positions(Out,3) + dZdX;

    %Find particles out of right bound and reposition
    v1 = C5-C1;
    v2 = C4-C1;
    C = C2;
    Out = FindOutOfBounds(v1,v2,C);
    Positions(Out,1) = Positions(Out,1) - Lx;
    Positions(Out,2) = Positions(Out,2) - dYdX;
    Positions(Out,3) = Positions(Out,3) - dZdX;

    %Find particles out of front bound and reposition
    v1 = C5-C1;
    v2 = C2-C1;
    C = C1;
    Out = FindOutOfBounds(v1,v2,C);
    Positions(Out,1) = Positions(Out,1) + dXdY;
    Positions(Out,2) = Positions(Out,2) + Ly;
    Positions(Out,3) = Positions(Out,3) + dZdY;

    %Find particles out of front bound and reposition
    v1 = C2-C1;
    v2 = C5-C1;
    C = C4;
    Out = FindOutOfBounds(v1,v2,C);
    Positions(Out,1) = Positions(Out,1) - dXdY;
    Positions(Out,2) = Positions(Out,2) - Ly;
    Positions(Out,3) = Positions(Out,3) - dZdY;

    %Find particles out of front bound and reposition
    v1 = C2-C1;
    v2 = C4-C1;
    C = C1;
    Out = FindOutOfBounds(v1,v2,C);
    Positions(Out,1) = Positions(Out,1) + dXdZ;
    Positions(Out,2) = Positions(Out,2) + dYdZ;
    Positions(Out,3) = Positions(Out,3) + Lz;

    %Find particles out of front bound and reposition
    v1 = C4-C1;
    v2 = C2-C1;
    C = C5;
    Out = FindOutOfBounds(v1,v2,C);
    Positions(Out,1) = Positions(Out,1) - dXdZ;
    Positions(Out,2) = Positions(Out,2) - dYdZ;
    Positions(Out,3) = Positions(Out,3) - Lz;
end
if ~isempty(ConnectionsSP)
    [X_SP,Y_SP,Z_SP] = ReIndexXY(ConnectionsSP,X_SP,Y_SP,Z_SP,Rx_SP);
end
if ~isempty(ConnectionsREP)
    [X_REP,Y_REP,Z_REP] = ReIndexXY(ConnectionsREP,X_REP,Y_REP,Z_REP,Rx_REP);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Out = FindOutOfBounds(v1,v2,C)

global Positions dims

n_left = cross(v1,v2)/norm(cross(v1,v2));
n_arr = [n_left(1)*ones(size(Positions,1),1),...
    n_left(2)*ones(size(Positions,1),1),...
    n_left(3)*ones(size(Positions,1),1)];
C1_arr = zeros(size(Positions(:,1:dims)));
C1_arr(:,1) = C(1); C1_arr(:,2) = C(2); C1_arr(:,3) = C(3);
r = Positions(:,1:dims)-C1_arr;
temp = dot(n_arr,r,2);
Out = find(temp<0);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function YofX = FindYofX(X,C2,C1)
m = (C2(2)-C1(2))/(C2(1)-C1(1)); %Left slope
YofX = C2(2)-m*(C2(1)-X);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XofY = FindXofY(Y,C2,C1)
m = (C2(2)-C1(2))/(C2(1)-C1(1)); %Left slope
XofY = C2(1)-(C2(2)-Y)/m;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,Z,Rx,Ry,Rz] = UpdateNodePositions(Connections,X,Y,Z,Rx,Ry,Rz)

%Updates Positions, X_SP, Y_SP, X_REP, Y_REP,...
%Rx_SP, Ry_SP, Rx_REP and Ry_REP

%No dynamics involved here

global Positions Positions0 dims

dX = Positions(:,1:dims)-Positions0(:,1:dims);

X(Connections(:,1:end-1)==0) = NaN;
Y(Connections(:,1:end-1)==0) = NaN;
Z(Connections(:,1:end-1)==0) = NaN;
Connections(Connections(:,1:end-1)==0) = NaN;

[X,Y,Z] = ReIndexXY(Connections,X,Y,Z,Rx);

%Define displacement of heads of Rx and Ry vectors
dx_Head_temp = dX(:,1);
dy_Head_temp = dX(:,2);
if dims==3
    dz_Head_temp = dX(:,3);
end
% dx_Head_temp = X(:,end)-X_0(:,end);
% dy_Head_temp = Y(:,end)-Y_0(:,end);

%Store dX for heads and tails in matrices the same size as Rx and Ry
dx_Tail = zeros(size(Rx));
dy_Tail = zeros(size(Rx));
dz_Tail = zeros(size(Rx));
dx_Head = zeros(size(Rx));
dy_Head = zeros(size(Rx));
dz_Head = zeros(size(Rx));
for j=1:size(Rx,2)
    TempNodes = Connections(:,j);
    Rows = find(~isnan(TempNodes));

    %Allocate dx_Head and dx_Tail seperately
    %Store change in head positions
    dx_Head(Rows,j) = dx_Head_temp(Rows);
    dy_Head(Rows,j) = dy_Head_temp(Rows);
    if dims==3
        dz_Head(Rows,j) = dz_Head_temp(Rows);
    end

    %Define tail displacements;
    dx_Tail(Rows,j) = dX(TempNodes(Rows),1);
    dy_Tail(Rows,j) = dX(TempNodes(Rows),2);
    if dims==3
        dz_Tail(Rows,j) = dX(TempNodes(Rows),3);
    end
end

Rx = Rx-dx_Tail+dx_Head;
Ry = Ry-dy_Tail+dy_Head;
if dims==3
    Rz = Rz-dz_Tail+dz_Head;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,Z] = ReIndexXYUnwrapped(Connections,X,Y,Z,Rx)

global PosUnwrapped dims

X(Connections(:,1:end-1)==0) = NaN;
Y(Connections(:,1:end-1)==0) = NaN;
Z(Connections(:,1:end-1)==0) = NaN;
Connections(Connections(:,1:end-1)==0) = NaN;

%Define New Positions in X and Y matrices
X(:,end) = PosUnwrapped(:,1);
Y(:,end) = PosUnwrapped(:,2);
if dims==3
    Z(:,end) = PosUnwrapped(:,3);
end
for j=1:size(Rx,2)
    TempNodes = Connections(:,j);
    Rows = find(~isnan(TempNodes));

    %Update X and Y of Tails
    X(Rows,j) = PosUnwrapped(TempNodes(Rows),1);
    Y(Rows,j) = PosUnwrapped(TempNodes(Rows),2);
    if dims==3
        Z(Rows,j) = PosUnwrapped(TempNodes(Rows),3);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,Z] = ReIndexXY(Connections,X,Y,Z,Rx)

global Positions dims

X(Connections(:,1:end-1)==0) = NaN;
Y(Connections(:,1:end-1)==0) = NaN;
Z(Connections(:,1:end-1)==0) = NaN;
Connections(Connections(:,1:end-1)==0) = NaN;

%Define New Positions in X and Y matrices
X(:,end) = Positions(:,1);
Y(:,end) = Positions(:,2);
if dims==3
    Z(:,end) = Positions(:,3);
end
for j=1:size(Rx,2)
    TempNodes = Connections(:,j);
    Rows = find(~isnan(TempNodes));

    %Update X and Y of Tails
    X(Rows,j) = Positions(TempNodes(Rows),1);
    Y(Rows,j) = Positions(TempNodes(Rows),2);
    if dims==3
        Z(Rows,j) = Positions(TempNodes(Rows),3);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_NET = CalculateNetForces
%Calculates forces based on Rx and Ry

global N dims

%Spring Forces
F_SP = CalculateSpringForces;
if isempty(F_SP)
    F_SP = zeros(N,dims);
end

%Particle-to-particle forces
F_REP = CalculateRepulsiveForces;

%Wall forces (if there are top and bottom plates
F_WALL = CalculateWallForces;

F_NET = F_SP+F_REP+F_WALL;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_REP = CalculateRepulsiveForces

global Rx_REP Ry_REP Rz_REP F_REP_X F_REP_Y F_REP_Z...
    CutoffCWAForce r0 dims

if dims==2
    Norms = (Rx_REP.^2 + Ry_REP.^2).^0.5;
elseif dims==3
    Norms = (Rx_REP.^2 + Ry_REP.^2 + Rz_REP.^2).^0.5;
end

F_MAG = CalculateInvToAlphaForce(Norms);

rCutoff = CutoffCWAForce*r0;
Force2 = CalculateInvToAlphaForce(rCutoff+0.0001);
Force1 = CalculateInvToAlphaForce(rCutoff);
slope = (Force2-Force1)/(0.0001);
yIntercept = Force2-slope*(rCutoff+0.0001);
F_MAG(Norms<=rCutoff) = Norms(Norms<=rCutoff)*slope+yIntercept;
F_MAG(Norms==0) = NaN;

F_REP_X = -F_MAG.*(Rx_REP./Norms);
F_REP_Y = -F_MAG.*(Ry_REP./Norms);
if dims==3
    F_REP_Z = -F_MAG.*(Rz_REP./Norms);
end

if dims==2
    F_REP = [nansum(F_REP_X,2),nansum(F_REP_Y,2)];
elseif dims==3
    F_REP = [nansum(F_REP_X,2),nansum(F_REP_Y,2),nansum(F_REP_Z,2)];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_WALL = CalculateWallForces

global Positions dims

F_WALL = zeros(size(Positions(:,1:dims)));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_MAG = CalculateInvToAlphaForce(r)

global ForceScaleREP sigmaREP alpha

F_MAG = alpha*ForceScaleREP*(1/sigmaREP - (sigmaREP^alpha)./(r.^(alpha+1)));
F_MAG(r>sigmaREP) = 0;
F_MAG(isinf(F_MAG)) = NaN;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Force = CalculateCWAForce(r)

global ForceScaleREP sigmaREP

eps = ForceScaleREP;
sig = sigmaREP;
ChanlderWeeksLim = 2^(1/6)*sig;

Force = 4*eps*((6*sig^6)./r.^7-(12*sig^12)./r.^13);
Force(r>ChanlderWeeksLim) = 0;
Force(r==0) = NaN;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_SP = CalculateSpringForces

%Intakes X and Y positions of connected springs to calculate force on each
%node due to springs.

global Rx_SP Ry_SP Rz_SP F_SP_X F_SP_Y F_SP_Z dims

if dims==2
    Norms = (Rx_SP.^2+Ry_SP.^2).^0.5;
elseif dims==3
    Norms = (Rx_SP.^2+Ry_SP.^2+Rz_SP.^2).^0.5;
end

F_MAG = CalculateForceVsStretch(Norms);
F_SP_X = F_MAG.*(Rx_SP./Norms);
F_SP_Y = F_MAG.*(Ry_SP./Norms);
if dims==3
    F_SP_Z = F_MAG.*(Rz_SP./Norms);
end

if dims==2
    F_SP = [nansum(F_SP_X,2),nansum(F_SP_Y,2)];
elseif dims==3
    F_SP = [nansum(F_SP_X,2),nansum(F_SP_Y,2),nansum(F_SP_Z,2)];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F =  CalculateForceVsStretch(Norms)

global SpringType r0 SpringCoeff_eq IsometricForce CoeffRatio L CutoffCWAForce...

if SpringType==0 %Linear
    F = -SpringCoeff_eq*Norms;%(r0-Norms);
elseif SpringType==1 %Langevin
    %     lam = Norms/L; %L is 2X the length of a tether
    F = CalculateLangevinForce(Norms);

    %Use Tangent modulus at high stretch for numerical stability
    StretchLimit = 0.90; %Percent
    StretchLimit_r = StretchLimit*L;
    Force1 = CalculateLangevinForce(StretchLimit_r-0.0001);
    Force2 = CalculateLangevinForce(StretchLimit_r+0.0001);
    slope = (Force2-Force1)/0.0002;
    yIntercept = Force2-slope*(StretchLimit_r+0.0001);
    F(Norms>StretchLimit_r) = Norms(Norms>StretchLimit_r)*slope+yIntercept;
elseif SpringType==2 % Hill Model
    F = -IsometricForce*exp(-(((Norms-r0)/r0)/SpringCoeff_eq).^2);
    F(Norms>r0) = -IsometricForce*(exp(-(((Norms(Norms>r0)-r0)/r0)/SpringCoeff_eq).^2)+(((Norms(Norms>r0)-r0)/r0)/(CoeffRatio*SpringCoeff_eq)).^2);
    F(Norms==0) = 0;
elseif SpringType==3
    lam = Norms/L;
    F = CalculateWLCForce(lam);

    %Use Tangent modulus at high stretch for numerical stability
    StretchLimit = 0.9;
    Force1 = CalculateWLCForce(StretchLimit-0.0001);
    Force2 = CalculateWLCForce(StretchLimit+0.0001);
    slope = (Force2-Force1)/0.0002;
    yIntercept = Force2-slope*(StretchLimit+0.0001);
    F(lam>StretchLimit) = lam(lam>StretchLimit)*slope+yIntercept;
elseif SpringType==4    %Entropic Langevin Chains with Chanlder-Weeks-Anderson Repulsion
    FL = CalculateLangevinForce(Norms); %Entropic contribution

    %Use Tangent modulus at high stretch for numerical stability
    StretchLimit = 0.90; %Percent
    StretchLimit_r = StretchLimit*L;
    Force1 = CalculateLangevinForce(StretchLimit_r-0.0001);
    Force2 = CalculateLangevinForce(StretchLimit_r+0.0001);
    slope = (Force2-Force1)/0.0002;
    yIntercept = Force2-slope*(StretchLimit_r+0.0001);
    F(Norms>StretchLimit_r) = Norms(Norms>StretchLimit_r)*slope+yIntercept;

    %Use Tangent modulus at low stretch for numerical stability
    FR = -CalculateCWAForce(Norms);  %Volume exclusion contribution

    rCutoff = CutoffCWAForce*r0;
    Force2 = -CalculateCWAForce(rCutoff+0.0001);
    slope = -CalculateCWASlope(rCutoff);
    yIntercept = Force2-slope*(rCutoff+0.0001);
    FR(Norms<=rCutoff) = Norms(Norms<=rCutoff)*slope+yIntercept;
    FR(Norms==0) = NaN;

    F = FR+FL;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Force = CalculateLangevinForce(r)

global N_Kuhn b kbT

lambda = r/(sqrt(N_Kuhn)*b);
K = kbT/(sqrt(N_Kuhn)*b);

Force = -K*lambda.*(lambda.^2 -  3*N_Kuhn)./(lambda.^2 - N_Kuhn);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ConnectREP

%This function defines all repelling neighbor's in one matrix, "ConnectionsREP".
%Each row of the matrix correspondes to one node, with the final column
%listing the node number. The columns from 1 to the maximum allowed number
%of connections lists the connected neighbors' node numbers.

%This function also defines the corresponding position matrices for the
%nodes listed in "ConnectionsREP" as X_REP and Y_REP.

%It also defines the components of the connection matrices spring lengths
%Rx_REP and Ry_REP, which will be used to update node positions upon force
%displacement

%Finally, it ID's which connections are accross which bounds via binary
%matrices. E.g., RLPairsREP lists 1 where ConnectionsREP corresponds to a
%spring going from the right to left bound. This is useful for plotting
%purposes

global Lx Ly Lz Positions CutoffPartREP dims...
    LeftBnds RightBnds FrontBnds BackBnds BottomBnds TopBnds...
    BottomLeftBnds BottomRightBnds TopLeftBnds TopRightBnds...
    BottomFrontBnds BottomBackBnds TopFrontBnds TopBackBnds...
    FrontLeftBnds FrontRightBnds BackLeftBnds BackRightBnds...
    C1Bnds C2Bnds C3Bnds C4Bnds C5Bnds C6Bnds C7Bnds C8Bnds...
    ConnectionsREP X_REP Y_REP Z_REP Rx_REP Ry_REP Rz_REP NormsREP...
    OldOrNew

% For Bulk
NoNodes = size(Positions,1);

if dims==2
    DomainSize = Lx*Ly;
    Density = NoNodes/DomainSize;
    AreaCutoff = pi*CutoffPartREP^2;
    AvgNoNeighbors = Density*AreaCutoff;
elseif dims==3
    DomainSize = Lx*Ly*Lz;
    Density = NoNodes/DomainSize;
    VolumeCutoff = 4/3*pi*CutoffPartREP^3;
    AvgNoNeighbors = Density*VolumeCutoff;
end
AbsMaxConnections = 2*ceil(AvgNoNeighbors);

NeighborsREP = zeros(NoNodes,AbsMaxConnections+1);
Rx_REP = zeros(NoNodes,AbsMaxConnections);
Ry_REP = zeros(NoNodes,AbsMaxConnections);
Rz_REP = zeros(NoNodes,AbsMaxConnections);
DistREP = zeros(NoNodes,AbsMaxConnections);

NeighborsREP(:,end) = Positions(:,end);

% For Right-Left Bnds
NeighborsREP_RightLeft = zeros(NoNodes,AbsMaxConnections+1);
Rx_REPRightLeft = zeros(NoNodes,AbsMaxConnections);
Ry_REPRightLeft = zeros(NoNodes,AbsMaxConnections);
Rz_REPRightLeft = zeros(NoNodes,AbsMaxConnections);
DistREP_RightLeft = zeros(NoNodes,AbsMaxConnections);

NeighborsREP_RightLeft(:,end) = Positions(:,end);

% For Top-Bottom Bnds
NeighborsREP_TopBottom = zeros(NoNodes,AbsMaxConnections+1);
Rx_REPTopBottom = zeros(NoNodes,AbsMaxConnections);
Ry_REPTopBottom = zeros(NoNodes,AbsMaxConnections);
Rz_REPTopBottom = zeros(NoNodes,AbsMaxConnections);
DistREP_TopBottom = zeros(NoNodes,AbsMaxConnections);

NeighborsREP_TopBottom(:,end) = Positions(:,end);

% For Bottom Left-Top Right Bnds
NeighborsREP_BottomLeftTopRight = zeros(NoNodes,AbsMaxConnections+1);
Rx_REPBottomLeftTopRight = zeros(NoNodes,AbsMaxConnections);
Ry_REPBottomLeftTopRight = zeros(NoNodes,AbsMaxConnections);
Rz_REPBottomLeftTopRight = zeros(NoNodes,AbsMaxConnections);
DistREP_BottomLeftTopRight = zeros(NoNodes,AbsMaxConnections);

NeighborsREP_BottomLeftTopRight(:,end) = Positions(:,end);

% For Bottom Right-Top Left Bnds
NeighborsREP_BottomRightTopLeft = zeros(NoNodes,AbsMaxConnections+1);
Rx_REPBottomRightTopLeft = zeros(NoNodes,AbsMaxConnections);
Ry_REPBottomRightTopLeft = zeros(NoNodes,AbsMaxConnections);
Rz_REPBottomRightTopLeft = zeros(NoNodes,AbsMaxConnections);
DistREP_BottomRightTopLeft = zeros(NoNodes,AbsMaxConnections);

NeighborsREP_BottomRightTopLeft(:,end) = Positions(:,end);

if dims==3
    % 1 Additional Face Combo
    % For Back-Front Bnds
    NeighborsREP_BackFront = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPBackFront = zeros(NoNodes,AbsMaxConnections);
    Ry_REPBackFront = zeros(NoNodes,AbsMaxConnections);
    Rz_REPBackFront = zeros(NoNodes,AbsMaxConnections);
    DistREP_BackFront = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_BackFront(:,end) = Positions(:,end);

    % 4 Additional Edge Combos
    % For Front Left-Back Right Bnds
    NeighborsREP_FrontLeftBackRight = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPFrontLeftBackRight = zeros(NoNodes,AbsMaxConnections);
    Ry_REPFrontLeftBackRight = zeros(NoNodes,AbsMaxConnections);
    Rz_REPFrontLeftBackRight = zeros(NoNodes,AbsMaxConnections);
    DistREP_FrontLeftBackRight = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_FrontLeftBackRight(:,end) = Positions(:,end);

    % For Front Right-Back Left Bnds
    NeighborsREP_FrontRightBackLeft = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPFrontRightBackLeft = zeros(NoNodes,AbsMaxConnections);
    Ry_REPFrontRightBackLeft = zeros(NoNodes,AbsMaxConnections);
    Rz_REPFrontRightBackLeft = zeros(NoNodes,AbsMaxConnections);
    DistREP_FrontRightBackLeft = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_FrontRightBackLeft(:,end) = Positions(:,end);

    % For Bottom Front-Top Back Bnds
    NeighborsREP_BottomFrontTopBack = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPBottomFrontTopBack = zeros(NoNodes,AbsMaxConnections);
    Ry_REPBottomFrontTopBack = zeros(NoNodes,AbsMaxConnections);
    Rz_REPBottomFrontTopBack = zeros(NoNodes,AbsMaxConnections);
    DistREP_BottomFrontTopBack = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_BottomFrontTopBack(:,end) = Positions(:,end);

    % For Top Front-Bottom Back Bnds
    NeighborsREP_TopFrontBottomBack = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPTopFrontBottomBack = zeros(NoNodes,AbsMaxConnections);
    Ry_REPTopFrontBottomBack = zeros(NoNodes,AbsMaxConnections);
    Rz_REPTopFrontBottomBack = zeros(NoNodes,AbsMaxConnections);
    DistREP_TopFrontBottomBack = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_TopFrontBottomBack(:,end) = Positions(:,end);

    % 4 Corner Combos
    % For Left Front Bottom-Right Back Top
    NeighborsREP_C1C7 = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPC1C7 = zeros(NoNodes,AbsMaxConnections);
    Ry_REPC1C7 = zeros(NoNodes,AbsMaxConnections);
    Rz_REPC1C7 = zeros(NoNodes,AbsMaxConnections);
    DistREP_C1C7 = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_C1C7(:,end) = Positions(:,end);

    % For Right Front Bottom-Left Back Top
    NeighborsREP_C2C8 = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPC2C8 = zeros(NoNodes,AbsMaxConnections);
    Ry_REPC2C8 = zeros(NoNodes,AbsMaxConnections);
    Rz_REPC2C8 = zeros(NoNodes,AbsMaxConnections);
    DistREP_C2C8 = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_C2C8(:,end) = Positions(:,end);

    % For Right Front Top-Left Back Bottom
    NeighborsREP_C6C4 = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPC6C4 = zeros(NoNodes,AbsMaxConnections);
    Ry_REPC6C4 = zeros(NoNodes,AbsMaxConnections);
    Rz_REPC6C4 = zeros(NoNodes,AbsMaxConnections);
    DistREP_C6C4 = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_C6C4(:,end) = Positions(:,end);

    % For Left Front Top-Right Back Bottom
    NeighborsREP_C5C3 = zeros(NoNodes,AbsMaxConnections+1);
    Rx_REPC5C3 = zeros(NoNodes,AbsMaxConnections);
    Ry_REPC5C3 = zeros(NoNodes,AbsMaxConnections);
    Rz_REPC5C3 = zeros(NoNodes,AbsMaxConnections);
    DistREP_C5C3 = zeros(NoNodes,AbsMaxConnections);

    NeighborsREP_C5C3(:,end) = Positions(:,end);
end

% tic
if OldOrNew==1
    % Define Connection Matrices in Bulk
    for P1=1:NoNodes
        Pool = (1:NoNodes)'; Pool(Pool==P1) = [];
        ExistingNeighb = (NeighborsREP(P1,1:end-1))';
        ExistingNeighb(ExistingNeighb==0) = [];
        ExistingNeighb = unique(ExistingNeighb);
        Pool(ismember(Pool,ExistingNeighb)) =[];    %Remove already matched pairs
        Pos1 = Positions(P1,1:dims).*ones(size(Pool,1),dims);
        Pos2 = Positions(Pool,1:dims);
        disp = Pos1-Pos2;
        dist = vecnorm(disp,2,2);
        Neighb = Pool(dist<=CutoffPartREP);
        for i=1:length(Neighb)
            P2 = Neighb(i);
            Pos1 = Positions(P1,1:dims);
            Pos2 = Positions(P2,1:dims);
            distance = Pos1-Pos2; %From P2 to P1
            FirstZero_n = find(NeighborsREP(P1,:)==0,1,'first');
            NeighborsREP(P1,FirstZero_n) = P2;
            Rx_REP(P1,FirstZero_n) = distance(1);
            Ry_REP(P1,FirstZero_n) = distance(2);
            if dims==3
                Rz_REP(P1,FirstZero_n) = distance(3);
            end
            DistREP(P1,FirstZero_n) = norm(distance);

            FirstZero_j = find(NeighborsREP(P2,:)==0,1,'first');
            NeighborsREP(P2,FirstZero_j) = P1;
            Rx_REP(P2,FirstZero_j) = -distance(1);
            Ry_REP(P2,FirstZero_j) = -distance(2);
            if dims==3
                Rz_REP(P2,FirstZero_j) = -distance(3);
            end
            DistREP(P2,FirstZero_j) = norm(distance);
        end
    end
else
    for P1=1:NoNodes-1
        for P2=P1+1:NoNodes
            Pos1 = Positions(P1,1:dims);
            Pos2 = Positions(P2,1:dims);
            distance = Pos1-Pos2; %From P2 to P1
            if norm(distance) <= CutoffPartREP %Find Potential Interactions
                FirstZero_n = find(NeighborsREP(P1,:)==0,1,'first');
                NeighborsREP(P1,FirstZero_n) = P2;
                Rx_REP(P1,FirstZero_n) = distance(1);
                Ry_REP(P1,FirstZero_n) = distance(2);
                if dims==3
                    Rz_REP(P1,FirstZero_n) = distance(3);
                end
                DistREP(P1,FirstZero_n) = norm(distance);

                FirstZero_j = find(NeighborsREP(P2,:)==0,1,'first');
                NeighborsREP(P2,FirstZero_j) = P1;
                Rx_REP(P2,FirstZero_j) = -distance(1);
                Ry_REP(P2,FirstZero_j) = -distance(2);
                if dims==3
                    Rz_REP(P2,FirstZero_j) = -distance(3);
                end
                DistREP(P2,FirstZero_j) = norm(distance);
            end
        end
    end
end
% toc
% CutoffPartREP = 10; %% DELETE
%Check accross Right-Left Bounds
Bnds1 = RightBnds;
Bnds2 = LeftBnds;
N1 = size(Bnds1,1);
N2 = size(Bnds2,1);
Shift = zeros(1,dims);
Shift(1) = Lx;
[NeighborsREP_RightLeft,DistREP_RightLeft,...
    Rx_REPRightLeft,Ry_REPRightLeft,Rz_REPRightLeft] = ...
    FindRepNeighbsBnds(NeighborsREP_RightLeft,...
    Rx_REPRightLeft,Ry_REPRightLeft,Rz_REPRightLeft,...
    DistREP_RightLeft,N1,N2,Bnds1,Bnds2,Shift);

%Check accross Top-Bottom Bounds
Bnds1 = TopBnds;
Bnds2 = BottomBnds;
N1 = size(Bnds1,1);
N2 = size(Bnds2,1);
Shift = zeros(1,dims);
if dims==2
    Shift(2) = Ly;
elseif dims==3
    Shift(3) = Lz;
end
[NeighborsREP_TopBottom,DistREP_TopBottom,...
    Rx_REPTopBottom,Ry_REPTopBottom,Rz_REPTopBottom] = ...
    FindRepNeighbsBnds(NeighborsREP_TopBottom,...
    Rx_REPTopBottom,Ry_REPTopBottom,Rz_REPTopBottom,...
    DistREP_TopBottom,N1,N2,Bnds1,Bnds2,Shift);

%Check accross Bottom Left-Top Right Bounds
Bnds1 = BottomLeftBnds;
Bnds2 = TopRightBnds;
N1 = size(Bnds1,1);
N2 = size(Bnds2,1);
Shift = zeros(1,dims);
if dims==2
    Shift(1) = -Lx;
    Shift(2) = -Ly;
elseif dims==3
    Shift(1) = -Lx;
    Shift(3) = -Lz;
end
[NeighborsREP_BottomLeftTopRight,DistREP_BottomLeftTopRight,...
    Rx_REPBottomLeftTopRight,Ry_REPBottomLeftTopRight,Rz_REPBottomLeftTopRight] = ...
    FindRepNeighbsBnds(NeighborsREP_BottomLeftTopRight,...
    Rx_REPBottomLeftTopRight,Ry_REPBottomLeftTopRight,Rz_REPBottomLeftTopRight,...
    DistREP_BottomLeftTopRight,N1,N2,Bnds1,Bnds2,Shift);

%Check accross Top Left-Bottom Right Bounds
Bnds1 = BottomRightBnds;
Bnds2 = TopLeftBnds;
N1 = size(Bnds1,1);
N2 = size(Bnds2,1);
Shift = zeros(1,dims);
if dims==2
    Shift(1) = +Lx;
    Shift(2) = -Ly;
elseif dims==3
    Shift(1) = +Lx;
    Shift(3) = -Lz;
end
[NeighborsREP_BottomRightTopLeft,DistREP_BottomRightTopLeft,...
    Rx_REPBottomRightTopLeft,Ry_REPBottomRightTopLeft,Rz_REPBottomRightTopLeft] = ...
    FindRepNeighbsBnds(NeighborsREP_BottomRightTopLeft,...
    Rx_REPBottomRightTopLeft,Ry_REPBottomRightTopLeft,Rz_REPBottomRightTopLeft,...
    DistREP_BottomRightTopLeft,N1,N2,Bnds1,Bnds2,Shift);

if dims==3
    % 1 Additional Face Combo
    % Check accross Back-Front bounds
    Bnds1 = FrontBnds;
    Bnds2 = BackBnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(2) = -Ly;
    [NeighborsREP_BackFront,DistREP_BackFront,...
        Rx_REPBackFront,Ry_REPBackFront,Rz_REPBackFront] = ...
        FindRepNeighbsBnds(NeighborsREP_BackFront,...
        Rx_REPBackFront,Ry_REPBackFront,Rz_REPBackFront,...
        DistREP_BackFront,N1,N2,Bnds1,Bnds2,Shift);

    % 4 Additional Edge Combos
    % Check accross Front Left-Back Right Bnds
    Bnds1 = FrontLeftBnds;
    Bnds2 = BackRightBnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(1) = -Lx;
    Shift(2) = -Ly;
    [NeighborsREP_FrontLeftBackRight,DistREP_FrontLeftBackRight,...
        Rx_REPFrontLeftBackRight,Ry_REPFrontLeftBackRight,Rz_REPFrontLeftBackRight] = ...
        FindRepNeighbsBnds(NeighborsREP_FrontLeftBackRight,...
        Rx_REPFrontLeftBackRight,Ry_REPFrontLeftBackRight,Rz_REPFrontLeftBackRight,...
        DistREP_FrontLeftBackRight,N1,N2,Bnds1,Bnds2,Shift);

    % Check accross Front Right-Back Left Bnds
    Bnds1 = FrontRightBnds;
    Bnds2 = BackLeftBnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(1) = +Lx;
    Shift(2) = -Ly;
    [NeighborsREP_FrontRightBackLeft,DistREP_FrontRightBackLeft,...
        Rx_REPFrontRightBackLeft,Ry_REPFrontRightBackLeft,Rz_REPFrontRightBackLeft] = ...
        FindRepNeighbsBnds(NeighborsREP_FrontRightBackLeft,...
        Rx_REPFrontRightBackLeft,Ry_REPFrontRightBackLeft,Rz_REPFrontRightBackLeft,...
        DistREP_FrontRightBackLeft,N1,N2,Bnds1,Bnds2,Shift);

    % Check accross Bottom Front-Top Back Bnds
    Bnds1 = BottomFrontBnds;
    Bnds2 = TopBackBnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(2) = -Ly;
    Shift(3) = -Lz;
    [NeighborsREP_BottomFrontTopBack,DistREP_BottomFrontTopBack,...
        Rx_REPBottomFrontTopBack,Ry_REPBottomFrontTopBack,Rz_REPBottomFrontTopBack] = ...
        FindRepNeighbsBnds(NeighborsREP_BottomFrontTopBack,...
        Rx_REPBottomFrontTopBack,Ry_REPBottomFrontTopBack,Rz_REPBottomFrontTopBack,...
        DistREP_BottomFrontTopBack,N1,N2,Bnds1,Bnds2,Shift);

    % Check accross Front Top-Back Bottom Bnds
    Bnds1 = TopFrontBnds;
    Bnds2 = BottomBackBnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(2) = -Ly;
    Shift(3) = +Lz;
    [NeighborsREP_TopFrontBottomBack,DistREP_TopFrontBottomBack,...
        Rx_REPTopFrontBottomBack,Ry_REPTopFrontBottomBack,Rz_REPTopFrontBottomBack] = ...
        FindRepNeighbsBnds(NeighborsREP_TopFrontBottomBack,...
        Rx_REPTopFrontBottomBack,Ry_REPTopFrontBottomBack,Rz_REPTopFrontBottomBack,...
        DistREP_TopFrontBottomBack,N1,N2,Bnds1,Bnds2,Shift);

    % 4 Corners
    % Front Bottom Left-Back Top Right
    Bnds1 = C1Bnds;
    Bnds2 = C7Bnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(1) = -Lx;
    Shift(2) = -Ly;
    Shift(3) = -Lz;
    [NeighborsREP_C1C7,DistREP_C1C7,...
        Rx_REPC1C7,Ry_REPC1C7,Rz_REPC1C7] = ...
        FindRepNeighbsBnds(NeighborsREP_C1C7,...
        Rx_REPC1C7,Ry_REPC1C7,Rz_REPC1C7,...
        DistREP_C1C7,N1,N2,Bnds1,Bnds2,Shift);

    % Front Bottom Right-Back Top Left
    Bnds1 = C2Bnds;
    Bnds2 = C8Bnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(1) = +Lx;
    Shift(2) = -Ly;
    Shift(3) = -Lz;
    [NeighborsREP_C2C8,DistREP_C2C8,...
        Rx_REPC2C8,Ry_REPC2C8,Rz_REPC2C8] = ...
        FindRepNeighbsBnds(NeighborsREP_C2C8,...
        Rx_REPC2C8,Ry_REPC2C8,Rz_REPC2C8,...
        DistREP_C2C8,N1,N2,Bnds1,Bnds2,Shift);

    % Front Top Right-Back Bottom Left
    Bnds1 = C6Bnds;
    Bnds2 = C4Bnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(1) = +Lx;
    Shift(2) = -Ly;
    Shift(3) = +Lz;
    [NeighborsREP_C6C4,DistREP_C6C4,...
        Rx_REPC6C4,Ry_REPC6C4,Rz_REPC6C4] = ...
        FindRepNeighbsBnds(NeighborsREP_C6C4,...
        Rx_REPC6C4,Ry_REPC6C4,Rz_REPC6C4,...
        DistREP_C6C4,N1,N2,Bnds1,Bnds2,Shift);

    % Front Top Left-Back Bottom Right
    Bnds1 = C5Bnds;
    Bnds2 = C3Bnds;
    N1 = size(Bnds1,1);
    N2 = size(Bnds2,1);
    Shift = zeros(1,dims);
    Shift(1) = -Lx;
    Shift(2) = -Ly;
    Shift(3) = +Lz;
    [NeighborsREP_C5C3,DistREP_C5C3,...
        Rx_REPC5C3,Ry_REPC5C3,Rz_REPC5C3] = ...
        FindRepNeighbsBnds(NeighborsREP_C5C3,...
        Rx_REPC5C3,Ry_REPC5C3,Rz_REPC5C3,...
        DistREP_C5C3,N1,N2,Bnds1,Bnds2,Shift);
end

% Sort Connected Springs by shortest to longest
if dims==2
    ConnectionsREP = [NeighborsREP(:,1:end-1),...
        NeighborsREP_RightLeft(:,1:end-1),...
        NeighborsREP_TopBottom(:,1:end-1),...
        NeighborsREP_BottomLeftTopRight(:,1:end-1),...
        NeighborsREP_BottomRightTopLeft(:,1:end-1)];
    Rx_REP = [Rx_REP,...
        Rx_REPRightLeft,...
        Rx_REPTopBottom,...
        Rx_REPBottomLeftTopRight,...
        Rx_REPBottomRightTopLeft];
    Ry_REP = [Ry_REP,...
        Ry_REPRightLeft,...
        Ry_REPTopBottom,...
        Ry_REPBottomLeftTopRight,...
        Ry_REPBottomRightTopLeft];
    Rz_REP = [Rz_REP,...
        Rz_REPRightLeft,...
        Rz_REPTopBottom,...
        Rz_REPBottomLeftTopRight,...
        Rz_REPBottomRightTopLeft];
else
    ConnectionsREP = [NeighborsREP(:,1:end-1),...
        NeighborsREP_RightLeft(:,1:end-1),...
        NeighborsREP_TopBottom(:,1:end-1),...
        NeighborsREP_BackFront(:,1:end-1),...
        NeighborsREP_BottomLeftTopRight(:,1:end-1),...
        NeighborsREP_BottomRightTopLeft(:,1:end-1),...
        NeighborsREP_BottomFrontTopBack(:,1:end-1),...
        NeighborsREP_TopFrontBottomBack(:,1:end-1),...
        NeighborsREP_FrontLeftBackRight(:,1:end-1),...
        NeighborsREP_FrontRightBackLeft(:,1:end-1),...
        NeighborsREP_C1C7(:,1:end-1),...
        NeighborsREP_C2C8(:,1:end-1),...
        NeighborsREP_C6C4(:,1:end-1),...
        NeighborsREP_C5C3(:,1:end-1)];
    Rx_REP = [Rx_REP,...
        Rx_REPRightLeft,...
        Rx_REPTopBottom,...
        Rx_REPBackFront,...
        Rx_REPBottomLeftTopRight,...
        Rx_REPBottomRightTopLeft,...
        Rx_REPBottomFrontTopBack,...
        Rx_REPTopFrontBottomBack,...
        Rx_REPFrontLeftBackRight,...
        Rx_REPFrontRightBackLeft,...
        Rx_REPC1C7,...
        Rx_REPC2C8,...
        Rx_REPC6C4,...
        Rx_REPC5C3];
    Ry_REP = [Ry_REP,...
        Ry_REPRightLeft,...
        Ry_REPTopBottom,...
        Ry_REPBackFront,...
        Ry_REPBottomLeftTopRight,...
        Ry_REPBottomRightTopLeft,...
        Ry_REPBottomFrontTopBack,...
        Ry_REPTopFrontBottomBack,...
        Ry_REPFrontLeftBackRight,...
        Ry_REPFrontRightBackLeft,...
        Ry_REPC1C7,...
        Ry_REPC2C8,...
        Ry_REPC6C4,...
        Ry_REPC5C3];
    Rz_REP = [Rz_REP,...
        Rz_REPRightLeft,...
        Rz_REPTopBottom,...
        Rz_REPBackFront,...
        Rz_REPBottomLeftTopRight,...
        Rz_REPBottomRightTopLeft,...
        Rz_REPBottomFrontTopBack,...
        Rz_REPTopFrontBottomBack,...
        Rz_REPFrontLeftBackRight,...
        Rz_REPFrontRightBackLeft,...
        Rz_REPC1C7,...
        Rz_REPC2C8,...
        Rz_REPC6C4,...
        Rz_REPC5C3];
end
if dims==2
    NormsREP = (Rx_REP.^2+Ry_REP.^2).^0.5;
else
    NormsREP = (Rx_REP.^2+Ry_REP.^2+Rz_REP.^2).^0.5;
end

Nodes = (1:NoNodes)';
ConnectionsREP(:,end) = Nodes;

[X_REP,Y_REP,Z_REP] = DefinePositionMatrix(ConnectionsREP);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,Z] = DefinePositionMatrix(Connections)

global Positions dims

X=zeros(size(Connections));
Y=zeros(size(Connections));
Z=zeros(size(Connections));

X(:,end) = Positions(:,1);
Y(:,end) = Positions(:,2);
if dims==3
    Z(:,end) = Positions(:,3);
end

for j=1:size(X,2)-1
    TempNodes = Connections(:,j);
    NonZeroIndices = find(TempNodes~=0);

    X(NonZeroIndices,j) = Positions(TempNodes(NonZeroIndices),1);
    Y(NonZeroIndices,j) = Positions(TempNodes(NonZeroIndices),2);
    if dims==3
        Z(NonZeroIndices,j) = Positions(TempNodes(NonZeroIndices),3);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NeighborsREP,DistREP,Rx_REP,Ry_REP,Rz_REP] = ...
    FindRepNeighbsBnds(NeighborsREP,Rx_REP,Ry_REP,Rz_REP,DistREP,...
    N1,N2,Bnds1,Bnds2,Shift)

global dims CutoffPartREP Positions OldOrNew

if OldOrNew==1
    for i=1:N1
        P1 = Bnds1(i);
        Pool = Bnds2; Pool(Pool==P1) = [];
        ExistingNeighb = (NeighborsREP(P1,1:end-1))';
        ExistingNeighb(ExistingNeighb==0) = [];
        ExistingNeighb = unique(ExistingNeighb);
        Pool(ismember(Pool,ExistingNeighb)) =[];    %Remove already matched pairs
        Pos1 = Positions(P1,1:dims).*ones(size(Pool,1),dims);
        Pos2 = Positions(Pool,1:dims) + Shift;
        disp = Pos1-Pos2;
        dist = vecnorm(disp,2,2);
        Neighb = Pool(dist<=CutoffPartREP);
        for j=1:length(Neighb)
            P2 = Neighb(j);
            Pos1 = Positions(P1,1:dims);
            Pos2 = Positions(P2,1:dims) + Shift;
            distance = Pos1-Pos2; %From P2 to P1
            FirstZero_n = find(NeighborsREP(P1,:)==0,1,'first');
            NeighborsREP(P1,FirstZero_n) = P2;
            Rx_REP(P1,FirstZero_n) = distance(1);
            Ry_REP(P1,FirstZero_n) = distance(2);
            if dims==3
                Rz_REP(P1,FirstZero_n) = distance(3);
            end
            DistREP(P1,FirstZero_n) = norm(distance);

            FirstZero_j = find(NeighborsREP(P2,:)==0,1,'first');
            NeighborsREP(P2,FirstZero_j) = P1;
            Rx_REP(P2,FirstZero_j) = -distance(1);
            Ry_REP(P2,FirstZero_j) = -distance(2);
            if dims==3
                Rz_REP(P2,FirstZero_j) = -distance(3);
            end
            DistREP(P2,FirstZero_j) = norm(distance);
        end
    end
else
    for i=1:N1
        P1 = Bnds1(i);
        for j=1:N2
            P2 = Bnds2(j);
            Pos1 = Positions(P1,1:dims);
            Pos2 = Positions(P2,1:dims) + Shift;
            distance = Pos1-Pos2; %From P2 to P1
            if norm(distance) <= CutoffPartREP %Find Potential Interactions
                FirstZero_n = find(NeighborsREP(P1,:)==0,1,'first');
                NeighborsREP(P1,FirstZero_n) = P2;
                Rx_REP(P1,FirstZero_n) = distance(1);
                Ry_REP(P1,FirstZero_n) = distance(2);
                if dims==3
                    Rz_REP(P1,FirstZero_n) = distance(3);
                end
                DistREP(P1,FirstZero_n) = norm(distance);

                FirstZero_j = find(NeighborsREP(P2,:)==0,1,'first');
                NeighborsREP(P2,FirstZero_j) = P1;
                Rx_REP(P2,FirstZero_j) = -distance(1);
                Ry_REP(P2,FirstZero_j) = -distance(2);
                if dims==3
                    Rz_REP(P2,FirstZero_j) = -distance(3);
                end
                DistREP(P2,FirstZero_j) = norm(distance);
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ReplicateBounds

%This function uses the updated boundary node numbers to replicate these
%crosslinks accross the periodic boundaries

global dims Positions Corners...
    LeftBnds RightBnds FrontBnds BackBnds BottomBnds TopBnds...
    BottomLeftBnds BottomRightBnds TopLeftBnds TopRightBnds...
    BottomFrontBnds BottomBackBnds TopFrontBnds TopBackBnds...
    FrontLeftBnds FrontRightBnds BackLeftBnds BackRightBnds...
    C1Bnds C2Bnds C3Bnds C4Bnds C5Bnds C6Bnds C7Bnds C8Bnds...
    PosLeftClns PosRightClns PosFrontClns PosBackClns PosBottomClns PosTopClns...
    PosBottomLeftClns PosBottomRightClns PosTopLeftClns PosTopRightClns...
    PosBottomFrontClns PosBottomBackClns PosTopFrontClns PosTopBackClns...
    PosFrontLeftClns PosFrontRightClns PosBackLeftClns PosBackRightClns...
    PosC1Clns PosC2Clns PosC3Clns PosC4Clns PosC5Clns PosC6Clns PosC7Clns PosC8Clns

%Replicate boundary nodes
if dims==2
    CBL = Corners(1,:);
    CBR = Corners(2,:);
    CTR = Corners(3,:);

    X = CBR(1)-CBL(1);
    Y = CTR(2)-CBR(2);
    dXdY = CTR(1)-CBR(1);
    dYdX = CBR(2)-CBL(2);

    % 4 Edges
    PosLeftClns = [Positions(LeftBnds,1) + X,...
        Positions(LeftBnds,2) + dYdX];
    PosRightClns = [Positions(RightBnds,1) - X,...
        Positions(RightBnds,2) - dYdX];
    PosBottomClns = [Positions(BottomBnds,1) + dXdY,...
        Positions(BottomBnds,2) + Y];
    PosTopClns = [Positions(TopBnds,1) - dXdY,...
        Positions(TopBnds,2) - Y];

    % 4 Corners
    PosBottomLeftClns = [Positions(BottomLeftBnds,1) + X + dXdY,...
        Positions(BottomLeftBnds,2) + Y + dYdX];
    PosBottomRightClns = [Positions(BottomRightBnds,1) - X + dXdY,...
        Positions(BottomRightBnds,2) + Y - dYdX];
    PosTopRightClns = [Positions(TopRightBnds,1) - X - dXdY,...
        Positions(TopRightBnds,2) - Y - dYdX];
    PosTopLeftClns = [Positions(TopLeftBnds,1) + X - dXdY,...
        Positions(TopLeftBnds,2) - Y + dYdX];
elseif dims==3 % ORTHONORMAL DEFORMATIONS ONLY
    C1 = Corners(1,:);  %Front-Bottom-Left
    C2 = Corners(2,:);  %Front-Bottom-Right
    C4 = Corners(4,:);  %Back-Bottom-Left
    C5 = Corners(5,:);  %Front-Top-Left

    X = C2(1)-C1(1);
    Y = C4(2)-C1(2);
    Z = C5(3)-C1(3);

    % 6 Faces
    PosLeftClns = [Positions(LeftBnds,1) + X,...
        Positions(LeftBnds,2),...
        Positions(LeftBnds,3)];
    PosRightClns = [Positions(RightBnds,1) - X,...
        Positions(RightBnds,2),...
        Positions(RightBnds,3)];
    PosFrontClns = [Positions(FrontBnds,1),...
        Positions(FrontBnds,2) + Y,...
        Positions(FrontBnds,3)];
    PosBackClns = [Positions(BackBnds,1),...
        Positions(BackBnds,2) - Y,...
        Positions(BackBnds,3)];
    PosBottomClns = [Positions(BottomBnds,1),...
        Positions(BottomBnds,2),...
        Positions(BottomBnds,3) + Z];
    PosTopClns = [Positions(TopBnds,1),...
        Positions(TopBnds,2),...
        Positions(TopBnds,3) - Z];

    % 12 Edges
    PosBottomLeftClns = [Positions(BottomLeftBnds,1) + X,...
        Positions(BottomLeftBnds,2),...
        Positions(BottomLeftBnds,3) + Z];
    PosBottomRightClns = [Positions(BottomRightBnds,1) - X,...
        Positions(BottomRightBnds,2),...
        Positions(BottomRightBnds,3) + Z];
    PosTopLeftClns = [Positions(TopLeftBnds,1) + X,...
        Positions(TopLeftBnds,2),...
        Positions(TopLeftBnds,3) - Z];
    PosTopRightClns = [Positions(TopRightBnds,1) - X,...
        Positions(TopRightBnds,2),...
        Positions(TopRightBnds,3) - Z];

    PosBottomFrontClns = [Positions(BottomFrontBnds,1),...
        Positions(BottomFrontBnds,2) + Y,...
        Positions(BottomFrontBnds,3) + Z];
    PosBottomBackClns = [Positions(BottomBackBnds,1),...
        Positions(BottomBackBnds,2) - Y,...
        Positions(BottomBackBnds,3) + Z];
    PosTopFrontClns = [Positions(TopFrontBnds,1),...
        Positions(TopFrontBnds,2) + Y,...
        Positions(TopFrontBnds,3) - Z];
    PosTopBackClns = [Positions(TopBackBnds,1),...
        Positions(TopBackBnds,2) - Y,...
        Positions(TopBackBnds,3) - Z];

    PosFrontLeftClns = [Positions(FrontLeftBnds,1) + X,...
        Positions(FrontLeftBnds,2) + Y,...
        Positions(FrontLeftBnds,3)];
    PosFrontRightClns = [Positions(FrontRightBnds,1) - X,...
        Positions(FrontRightBnds,2) + Y,...
        Positions(FrontRightBnds,3)];
    PosBackLeftClns = [Positions(BackLeftBnds,1) + X,...
        Positions(BackLeftBnds,2) - Y,...
        Positions(BackLeftBnds,3)];
    PosBackRightClns = [Positions(BackRightBnds,1) - X,...
        Positions(BackRightBnds,2) - Y,...
        Positions(BackRightBnds,3)];

    % 8 Corners
    PosC1Clns = [Positions(C1Bnds,1) + X,...
        Positions(C1Bnds,2) + Y,...
        Positions(C1Bnds,3) + Z];
    PosC2Clns = [Positions(C2Bnds,1) - X,...
        Positions(C2Bnds,2) + Y,...
        Positions(C2Bnds,3) + Z];
    PosC3Clns = [Positions(C3Bnds,1) - X,...
        Positions(C3Bnds,2) - Y,...
        Positions(C3Bnds,3) + Z];
    PosC4Clns = [Positions(C4Bnds,1) + X,...
        Positions(C4Bnds,2) - Y,...
        Positions(C4Bnds,3) + Z];
    PosC5Clns = [Positions(C5Bnds,1) + X,...
        Positions(C5Bnds,2) + Y,...
        Positions(C5Bnds,3) - Z];
    PosC6Clns = [Positions(C6Bnds,1) - X,...
        Positions(C6Bnds,2) + Y,...
        Positions(C6Bnds,3) - Z];
    PosC7Clns = [Positions(C7Bnds,1) - X,...
        Positions(C7Bnds,2) - Y,...
        Positions(C7Bnds,3) - Z];
    PosC8Clns = [Positions(C8Bnds,1) + X,...
        Positions(C8Bnds,2) - Y,...
        Positions(C8Bnds,3) - Z];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FindBoundaryNodes

%This funciton identifies the node numbers corresponding to each boundary
%as a function of their positions.

global Lx Ly Lz Positions Corners dims...
    LeftBnds RightBnds FrontBnds BackBnds BottomBnds TopBnds...
    BottomLeftBnds BottomRightBnds TopLeftBnds TopRightBnds...
    BottomFrontBnds BottomBackBnds TopFrontBnds TopBackBnds...
    FrontLeftBnds FrontRightBnds BackLeftBnds BackRightBnds...
    C1Bnds C2Bnds C3Bnds C4Bnds C5Bnds C6Bnds C7Bnds C8Bnds...
    AllBnds

MultFactor = 0.5;

NoNodes = size(Positions,1);
if dims==2
    C1 = Corners(1,:);  %BL
    C2 = Corners(2,:);  %BR
    C3 = Corners(3,:); %TR
    C4 = Corners(4,:); %TL

    BottomBnds = PointSlopeProximityFind(C2,C1,MultFactor);
    TopBnds = PointSlopeProximityFind(C4,C3,MultFactor);
    RightBnds = PointSlopeProximityFind(C2,C3,MultFactor);
    LeftBnds = PointSlopeProximityFind(C4,C1,MultFactor);

    BottomLeftBnds = BottomBnds(ismember(BottomBnds,LeftBnds));
    BottomRightBnds = BottomBnds(ismember(BottomBnds,RightBnds));
    TopRightBnds = TopBnds(ismember(TopBnds,RightBnds));
    TopLeftBnds = TopBnds(ismember(TopBnds,LeftBnds));
elseif dims==3
    % CONVENTION:
    % Left-to-Right is (-) to (+) x-axis
    % Front-to-Back is (-) to (+) y-axis
    % Bottom-to-Top is (-) to (+) z-axis
    C1 = Corners(1,:);  %LFB    Left-Front-Bottom
    %     C2 = Corners(2,:);  %RFB    Right-Front-Bottom
    %     C3 = Corners(3,:);  %RBB    Right-Back-Bottom
    %     C4 = Corners(4,:);  %LBB    Left-Back-Bottom
    %     C5 = Corners(5,:);  %LFT    Left-Front-Top
    %     C6 = Corners(6,:);  %RFT    Right-Front-Top
    C7 = Corners(7,:);  %RBT    Right-Back-Top
    %     C8 = Corners(8,:);  %LBT    Left-Back-Top

    %For orthotropic boundaries only
    xNorm = [ones(NoNodes,1) zeros(NoNodes,2)];
    yNorm = [zeros(NoNodes,1) ones(NoNodes,1) zeros(NoNodes,1)];
    zNorm = [zeros(NoNodes,2) ones(NoNodes,1)];

    DistToLeft = dot(Positions(:,1:dims)-C1,xNorm,2);
    DistToRight = -dot(Positions(:,1:dims)-C7,xNorm,2);
    DistToFront = dot(Positions(:,1:dims)-C1,yNorm,2);
    DistToBack = -dot(Positions(:,1:dims)-C7,yNorm,2);
    DistToBottom = dot(Positions(:,1:dims)-C1,zNorm,2);
    DistToTop = -dot(Positions(:,1:dims)-C7,zNorm,2);

    % 6 Faces
    LeftBnds = find(DistToLeft<MultFactor*Lx/2);
    RightBnds = find(DistToRight<MultFactor*Lx/2);
    FrontBnds = find(DistToFront<MultFactor*Ly/2);
    BackBnds = find(DistToBack<MultFactor*Ly/2);
    BottomBnds = find(DistToBottom<MultFactor*Lz/2);
    TopBnds = find(DistToTop<MultFactor*Lz/2);

    % 12 Edge Overlaps
    BottomLeftBnds = BottomBnds(ismember(BottomBnds,LeftBnds));
    BottomRightBnds = BottomBnds(ismember(BottomBnds,RightBnds));
    TopLeftBnds = TopBnds(ismember(TopBnds,LeftBnds));
    TopRightBnds = TopBnds(ismember(TopBnds,RightBnds));

    BottomFrontBnds = BottomBnds(ismember(BottomBnds,FrontBnds));
    BottomBackBnds = BottomBnds(ismember(BottomBnds,BackBnds));
    TopFrontBnds = TopBnds(ismember(TopBnds,FrontBnds));
    TopBackBnds = TopBnds(ismember(TopBnds,BackBnds));

    FrontLeftBnds = FrontBnds(ismember(FrontBnds,LeftBnds));
    BackLeftBnds = BackBnds(ismember(BackBnds,LeftBnds));
    FrontRightBnds = FrontBnds(ismember(FrontBnds,RightBnds));
    BackRightBnds = BackBnds(ismember(BackBnds,RightBnds));

    % 8 Corner Overlaps
    C1Bnds = FrontLeftBnds(ismember(FrontLeftBnds,BottomFrontBnds));
    C2Bnds = FrontRightBnds(ismember(FrontRightBnds,BottomRightBnds));
    C3Bnds = BackRightBnds(ismember(BackRightBnds,BottomBackBnds));
    C4Bnds = BackLeftBnds(ismember(BackLeftBnds,BottomBackBnds));
    C5Bnds = FrontLeftBnds(ismember(FrontLeftBnds,TopFrontBnds));
    C6Bnds = FrontRightBnds(ismember(FrontRightBnds,TopFrontBnds));
    C7Bnds = BackRightBnds(ismember(BackRightBnds,TopBackBnds));
    C8Bnds = BackLeftBnds(ismember(BackLeftBnds,TopLeftBnds));

    ShowFigDebug = 1;
    CheckSurfs = 0;
    CheckEdges = 0;
    CheckCorners = 0;
    if ShowFigDebug==1
        if CheckSurfs==1
            hold on
            s = scatter3(Positions(LeftBnds,1),Positions(LeftBnds,2),Positions(LeftBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(RightBnds,1),Positions(RightBnds,2),Positions(RightBnds,3),'r','filled');
            s.SizeData = 6;

            s = scatter3(Positions(TopBnds,1),Positions(TopBnds,2),Positions(TopBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(BottomBnds,1),Positions(BottomBnds,2),Positions(BottomBnds,3),'r','filled');
            s.SizeData = 6;

            s = scatter3(Positions(FrontBnds,1),Positions(FrontBnds,2),Positions(FrontBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(BackBnds,1),Positions(BackBnds,2),Positions(BackBnds,3),'r','filled');
            s.SizeData = 6;
        end

        if CheckEdges==1
            s = scatter3(Positions(FrontLeftBnds,1),Positions(FrontLeftBnds,2),Positions(FrontLeftBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(FrontRightBnds,1),Positions(FrontRightBnds,2),Positions(FrontRightBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(BackLeftBnds,1),Positions(BackLeftBnds,2),Positions(BackLeftBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(BackRightBnds,1),Positions(BackRightBnds,2),Positions(BackRightBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(BottomLeftBnds,1),Positions(BottomLeftBnds,2),Positions(BottomLeftBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(BottomRightBnds,1),Positions(BottomRightBnds,2),Positions(BottomRightBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(TopLeftBnds,1),Positions(TopLeftBnds,2),Positions(TopLeftBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(TopRightBnds,1),Positions(TopRightBnds,2),Positions(TopRightBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(BottomFrontBnds,1),Positions(BottomFrontBnds,2),Positions(BottomFrontBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(BottomBackBnds,1),Positions(BottomBackBnds,2),Positions(BottomBackBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(TopFrontBnds,1),Positions(TopFrontBnds,2),Positions(TopFrontBnds,3),'r','filled');
            s.SizeData = 6;
            s = scatter3(Positions(TopBackBnds,1),Positions(TopBackBnds,2),Positions(TopBackBnds,3),'r','filled');
            s.SizeData = 6;
        end

        if CheckCorners==1
            if ~isempty(C1Bnds)
                s = scatter3(Positions(C1Bnds,1),Positions(C1Bnds,2),Positions(C1Bnds,3),'r','filled');
                s.SizeData = 12;
            end
            if ~isempty(C2Bnds)
                s = scatter3(Positions(C2Bnds,1),Positions(C2Bnds,2),Positions(C2Bnds,3),'r','filled');
                s.SizeData = 12;
            end
            if ~isempty(C3Bnds)
                s = scatter3(Positions(C3Bnds,1),Positions(C3Bnds,2),Positions(C3Bnds,3),'r','filled');
                s.SizeData = 12;
            end
            if ~isempty(C4Bnds)
                s = scatter3(Positions(C4Bnds,1),Positions(C4Bnds,2),Positions(C4Bnds,3),'r','filled');
                s.SizeData = 12;
            end
            if ~isempty(C5Bnds)
                s = scatter3(Positions(C5Bnds,1),Positions(C5Bnds,2),Positions(C5Bnds,3),'r','filled');
                s.SizeData = 12;
            end
            if ~isempty(C6Bnds)
                s = scatter3(Positions(C6Bnds,1),Positions(C6Bnds,2),Positions(C6Bnds,3),'r','filled');
                s.SizeData = 12;
            end
            if ~isempty(C7Bnds)
                s = scatter3(Positions(C7Bnds,1),Positions(C7Bnds,2),Positions(C7Bnds,3),'r','filled');
                s.SizeData = 12;
            end
            if ~isempty(C8Bnds)
                s = scatter3(Positions(C8Bnds,1),Positions(C8Bnds,2),Positions(C8Bnds,3),'r','filled');
                s.SizeData = 12;
            end
        end
    end
end

if dims==2
    AllBnds = [LeftBnds;RightBnds;BottomBnds;TopBnds;...
        BottomLeftBnds;BottomRightBnds;TopLeftBnds;TopRightBnds];
elseif dims==3
    AllBnds = [LeftBnds;RightBnds;FrontBnds;BackBnds;BottomBnds;TopBnds;...
        BottomLeftBnds;BottomRightBnds;TopLeftBnds;TopRightBnds;...
        BottomFrontBnds;BottomBackBnds;TopFrontBnds;TopBackBnds;...
        FrontLeftBnds;FrontRightBnds;BackLeftBnds;BackRightBnds;...
        C1Bnds;C2Bnds;C3Bnds;C4Bnds;C5Bnds;C6Bnds;C7Bnds;C8Bnds];
end
AllBnds = unique(AllBnds);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Bnds = PointSlopeProximityFind(C1,C2,MultFactor)

global Positions Corners

CBL = Corners(1,:);
CBR = Corners(2,:);
CTR = Corners(3,:);
% CTL = Corners(1,:);

%Define min distance between any two corners to determine boundary
%proximity
Lx = CBR(1)-CBL(1);
Ly = CTR(2)-CBR(2);
% D1 = norm(CBL-CBR);
% D2 = norm(CBR-CTR);
% D3 = norm(CTR-CBL);
% D4 = norm(CTL-CBR);
% DomainScale = min([D1,D2,D3,D4]);

%Find points on edge closest to Positions
if C1(2)~=C2(2) && C1(1)~=C2(1)
    Slope = (C1(2)-C2(2))/(C1(1)-C2(1));

    K = Slope+1/Slope;
    Px = (1/K)*(Positions(:,2)-C1(2)+1/Slope*Positions(:,1)+Slope*C1(1));
    Py = Positions(:,2) + 1/Slope*(Positions(:,1)-Px);

    XDistToOrigin = 1/K*(Slope*C1(1)-C1(2));
    YDistToOrigin = -1/Slope*XDistToOrigin;
    DistToOrigin = [XDistToOrigin,YDistToOrigin];
elseif C1(1)==C2(1) %If a vertical line (slope == Infty)
    Px = C1(1)*ones(size(Positions,1),1);
    Py = Positions(:,2);
    DistToOrigin = [Lx/2,0];
elseif C1(2)==C2(2) %If a horizontal line
    Px = Positions(:,1);
    Py = C1(2)*ones(size(Positions,1),1);
    DistToOrigin = [0,Ly/2];
end
Norms = vecnorm([Px,Py]-Positions(:,1:2),2,2);
Bnds = find(Norms<=MultFactor*norm(DistToOrigin));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SizeTheDomain

global Lx Ly Lz MaxDist MinDist...
    dims phi b Np Nt  N_Kuhn sigmaREP CutoffPartREP

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

Lx = v_domain^(1/dims);
Ly = v_domain^(1/dims);
if dims==2
    Lz = 0;
else
    Lz = Lx;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitializeNodes(Controls)

global Lx Ly Lz Positions Corners dims Np %N

ShowFigDebug = 0;

if dims==2 %Generate Poissons's Distribution
    InitiateNetworkPoissons2D; %Function that generates Poisson's dist. network
else
    if Controls.RunLargeDeformation==1 || Np>1e3
        InitiateNetworkGrid3D;
    else
        InitiateNetworkPoissons3D; %Function that generates Poisson's dist. network
    end
end

N = size(Positions,1);

% end

Positions(:,end+1) = (1:N)';    %Index No.

if dims==2
    Corners = [-Lx/2 -Ly/2;
        Lx/2 -Ly/2;
        Lx/2  Ly/2;
        -Lx/2  Ly/2];
elseif dims==3
    Corners = [-Lx/2 -Ly/2 -Lz/2;
        Lx/2 -Ly/2 -Lz/2;
        Lx/2  Ly/2 -Lz/2;
        -Lx/2  Ly/2 -Lz/2;
        -Lx/2 -Ly/2 Lz/2;
        Lx/2 -Ly/2 Lz/2;
        Lx/2  Ly/2 Lz/2;
        -Lx/2  Ly/2 Lz/2;];
end

CheckCorners = 0;
if ShowFigDebug==1
    figure(1)
    clf
    hold on
    if dims==2
        scatter(Positions(:,1),Positions(:,2),'k','filled')
        if CheckCorners==1
            scatter(Corners(:,1),Corners(:,2),'r','filled')
        end
    else
        s = scatter3(Positions(:,1),Positions(:,2),Positions(:,3),'k','filled');
        s.SizeData = 6;
        if CheckCorners==1
            s = scatter3(Corners(:,1),Corners(:,2),Corners(:,3),'r','filled');
            s.SizeData = 6;
        end
        zlim([-Lz/2 Lz/2])
        pbaspect([1 1 1])
        view([30,30])
    end
    xlim([-Lx/2 Lx/2])
    ylim([-Ly/2 Ly/2])
    daspect([1 1 1])
    pbaspect([1 1 1])
    title(['N = ',num2str(size(Positions,1))])
    xlabel('x')
    ylabel('y')
    zlabel('z')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitiateNetworkPoissons2D

global Positions Lx Ly N MaxDist MinDist

%start with a single random point
dmax = MaxDist;
dmin = MinDist;


% minDistance = spacing + 2*r ;
mRandomness = 0.1;
% numCC = [] ;

%initial conditions
point= [ 0 0 ] ;

%we want minDistance to be small when density is high and big when
%density is small.

%packing is how many points around each point we'll test for a new
%location. A high randomness means we'll have a low number (more
%random) and a low randomness means a high number (more uniform).
Q = 100;
packing = Q + floor( (1.0 - mRandomness)*90 ) ;

numC = 1 ;
ck = 1 ;

wb = waitbar(0,'Initializing the nodes...');
while( ck > 0 )

    if numC>length(point)
        break
    end

    pointCen = point(numC,:) ;
    pointR = [] ;
    for i=1: packing
        radius = dmin +(dmax-dmin)*rand() ;
        angle  = 2*pi*rand() ;
        newpoint = [ pointCen(1)+radius*cos(angle) pointCen(2)+radius*sin(angle) ] ;
        nx = newpoint(1) ; ny = newpoint(2) ;
        % consider boundary domain
        if(nx<-Lx/2 || nx >Lx/2)
        elseif(ny<-Ly/2 || ny >Ly/2)
        else
            pointR = [ pointR  ; newpoint ] ;
        end
    end

    point ;

    % Eliminate the circle points if the distance between preoccupy points
    % and the circle points are smaller than minimun distance.
    count = 0 ;
    pointR1= [] ;
    for i=1: length(pointR)
        x1 = pointR(i,1) ; y1 = pointR(i,2) ;
        d1 = sqrt( ( x1 - point(:,1) ).^2 + ( y1 - point(:,2) ).^2 ) ;
        Re1 = [] ; Re1 = find(d1<dmin ) ;
        Lre1 = length(Re1) ;
        if Lre1 == 0
            count = count + 1 ;
            pointR1(count,:) = pointR(i,:)  ;
        end
    end

    nprt = size(pointR1,1) ;

    if(nprt == 0 )

    elseif(nprt == 1 )
        pointR = zeros(nprt,2) ;
        pointR = pointR1 ;
        point = [ point; pointR ] ;

    else
        pointR = zeros(nprt,2) ;
        pointR = pointR1 ;

        IDstack = [] ; pointR2 = [] ; IDs = [];
        for i=1: length(pointR)

            if (i==1)
                x1 = pointR(i,1) ; y1 = pointR(i,2) ;
                d =  sqrt( ( x1 - pointR(:,1) ).^2 + ( y1 - pointR(:,2) ).^2 ) ;

                Remain = find(d<dmin) ;
                IDstack = [ IDstack Remain']   ;
                IDs = [IDs 1 ] ;
            else
                sum1 = sum(ismember(IDstack,i)) ;
                if sum1 == 0
                    x1 = pointR(i,1) ; y1 = pointR(i,2) ;
                    d =  sqrt( ( x1 - pointR(:,1) ).^2 + ( y1 - pointR(:,2) ).^2 ) ;

                    Remain = find(d<dmin) ;
                    IDstack = [ IDstack Remain'] ;
                    IDs = [IDs i ] ;
                end
            end
        end
        IDs ;
        length(IDs) ;
        count1= 0 ;
        for i=1: length(pointR)

            if i==1
                count1=count1+1 ;
                pointR2(count1,:) = pointR(i,:) ;
            else
                sum1 = sum(ismember(IDs,i)) ;
                if sum1 == 1
                    count1=count1+1;
                    pointR2(count1,:) = pointR(i,:) ;
                end
            end

        end
        nprt1 = length(pointR2) ;
        pointR = zeros(nprt1,2) ;
        pointR = pointR2;

        point = [ point; pointR ] ;

    end
    % disp('==========================================')
    numC = numC+1 ;
    if ~mod(numC,100)
        waitbar(numC/length(point),wb,'Initializing the nodes...')
    end
end
close(wb)

Positions  = point  ;

% NumPoints = [ numCC ; size(Positions,1) ];

yRelevantPoints = find(abs(Positions(:,1)) <= abs(Positions(:,2)));
xRelevantPoints = find(abs(Positions(:,1)) > abs(Positions(:,2)));
ydistFromOrigin = abs(Positions(yRelevantPoints,2));
xdistFromOrigin = abs(Positions(xRelevantPoints,1));
ydistFromOrigin = [yRelevantPoints,ydistFromOrigin];
xdistFromOrigin = [xRelevantPoints,xdistFromOrigin];
NormDistFromOrigin = [ydistFromOrigin;xdistFromOrigin];
NormDistFromOrigin = sortrows(NormDistFromOrigin,2);
RemoveIndices = NormDistFromOrigin(N+1:end,1);
Positions(RemoveIndices,:) = [];
% NumPoints = size(Positions,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitiateNetworkPoissons3D

global Positions Lx Ly Lz MaxDist MinDist N

%start with a single random point
dmax = MaxDist;
dmin = MinDist;

% minDistance = spacing + 2*r ;
mRandomness = 0.1;
% numCC = [] ;

%initial conditions
point= [ 0 0 0 ] ;

%we want minDistance to be small when density is high and big when
%density is small.

%packing is how many points around each point we'll test for a new
%location. A high randomness means we'll have a low number (more
%random) and a low randomness means a high number (more uniform).
Q = 100;
packing = Q + floor( (1.0 - mRandomness)*90 ) ;

numC = 1 ;
ck = 1 ;

wb = waitbar(0,'Initializing the nodes...');
while( ck > 0 )

    if numC>length(point)
        break
    end

    pointCen = point(numC,:) ;
    pointR = [] ;
    for i=1: packing
        radius = dmin +(dmax-dmin)*rand();
        theta = 2*pi*rand();
        varphimin = -pi/2;
        varphimax = pi/2;
        varphi = varphimin+rand()*(varphimax-varphimin);

        %         varphi = pi*rand();
        radius_projected = radius*cos(varphi);
        dx = radius_projected*cos(theta);
        dy = radius_projected*sin(theta);
        dz = radius*sin(varphi);
        newpoint = [pointCen(1)+dx,...
            pointCen(2)+dy,...
            pointCen(3)+dz] ;
        nx = newpoint(1); ny = newpoint(2); nz = newpoint(3);
        % consider boundary domain
        if(nx<-Lx/2 || nx>Lx/2)
        elseif(ny<-Ly/2 || ny>Ly/2)
        elseif(nz<-Lz/2 || nz>Lz/2)
        else
            pointR = [ pointR  ; newpoint ] ;
        end
    end

    point ;

    % Eliminate the circle points if the distance between preoccupy points
    % and the circle points are smaller than minimun distance.
    count = 0 ;
    pointR1= [] ;
    for i=1:size(pointR,1)
        x1 = pointR(i,1); y1 = pointR(i,2); z1 = pointR(i,3);
        d1 = sqrt((x1 - point(:,1)).^2 + (y1 - point(:,2)).^2 + ...
            (z1 - point(:,3)).^2);
        Re1 = [] ; Re1 = find(d1<dmin );
        Lre1 = length(Re1) ;
        if Lre1 == 0
            count = count + 1 ;
            pointR1(count,:) = pointR(i,:) ;
        end
    end

    nprt = size(pointR1,1) ;

    if(nprt == 0 )

    elseif(nprt == 1 )
        pointR = zeros(nprt,2) ;
        pointR = pointR1 ;
        point = [ point; pointR ] ;

    else
        pointR = zeros(nprt,2) ;
        pointR = pointR1 ;

        IDstack = [] ; pointR2 = [] ; IDs = [];
        for i=1:size(pointR,1)

            if (i==1)
                x1 = pointR(i,1); y1 = pointR(i,2); z1 = pointR(i,3);
                d =  sqrt((x1-pointR(:,1)).^2 + (y1 - pointR(:,2)).^2 + ...
                    (z1 - pointR(:,3)).^2) ;

                Remain = find(d<dmin) ;
                IDstack = [ IDstack Remain']   ;
                IDs = [IDs 1 ] ;
            else
                sum1 = sum(ismember(IDstack,i)) ;
                if sum1 == 0
                    x1 = pointR(i,1); y1 = pointR(i,2); z1 = pointR(i,3);
                    d =  sqrt((x1 - pointR(:,1)).^2 + ...
                        (y1 - pointR(:,2)).^2 + ...
                        (z1 - pointR(:,3)).^2) ;

                    Remain = find(d<dmin) ;
                    IDstack = [ IDstack Remain'] ;
                    IDs = [IDs i ] ;
                end
            end
        end
        IDs ;
        length(IDs) ;
        count1= 0 ;
        for i=1: length(pointR)

            if i==1
                count1=count1+1 ;
                pointR2(count1,:) = pointR(i,:) ;
            else
                sum1 = sum(ismember(IDs,i)) ;
                if sum1 == 1
                    count1=count1+1;
                    pointR2(count1,:) = pointR(i,:) ;
                end
            end

        end
        nprt1 = length(pointR2) ;
        pointR = zeros(nprt1,2) ;
        pointR = pointR2;

        point = [ point; pointR ] ;

    end
    % disp('==========================================')
    numC = numC+1 ;

    if ~mod(numC,100)
        waitbar(numC/length(point),wb,'Initializing the nodes...')
    end
end
close(wb)

Positions  = point  ;

% NumPoints = [ numCC ; size(Positions,1) ];

xRelevantPoints = []; yRelevantPoints = []; zRelevantPoints = [];
for i=1:size(point,1)
    MaxComp = max(abs(Positions(i,:)));
    if MaxComp==abs(Positions(i,1))
        xRelevantPoints = [xRelevantPoints;i];
    elseif MaxComp==abs(Positions(i,2))
        yRelevantPoints = [yRelevantPoints;i];
    else
        zRelevantPoints = [zRelevantPoints;i];
    end
end
xdistFromOrigin = abs(Positions(xRelevantPoints,1));
ydistFromOrigin = abs(Positions(yRelevantPoints,2));
zdistFromOrigin = abs(Positions(zRelevantPoints,3));

xdistFromOrigin = [xRelevantPoints,xdistFromOrigin];
ydistFromOrigin = [yRelevantPoints,ydistFromOrigin];
zdistFromOrigin = [zRelevantPoints,zdistFromOrigin];

NormDistFromOrigin = [ydistFromOrigin;xdistFromOrigin;zdistFromOrigin];
NormDistFromOrigin = sortrows(NormDistFromOrigin,2);
RemoveIndices = NormDistFromOrigin(N+1:end,1);
Positions(RemoveIndices,:) = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitiateNetworkGrid3D

global Positions Lx Ly Lz N Nt Ns Np

pts_per_edge = ceil(N^(1/3));

dist = Lx/(pts_per_edge);
x_in = ((-Lx+dist)/2:dist:(Lx-dist)/2)';

[x,y,z] = meshgrid(x_in,x_in,x_in);

Positions_tmp = [x(:),y(:),z(:)];
Positions = Positions_tmp;

% Delete extraneous points at random
if size(Positions_tmp,1)>N
    n_to_del = size(Positions_tmp,1)-N;
    rem_indx = randperm(size(Positions_tmp,1),n_to_del);
    Positions(rem_indx,:) = [];
end
   
% Perturb the points slightly
% Generate random noise in the range (-d, d) for each element
max_dist = 0.2*dist;
noise = 2 * max_dist * rand(size(Positions)) - max_dist;

% Add noise to each element of the Positions matrix
Positions = Positions + noise;

CheckFig = 1;
if CheckFig==1
    figure(1); clf; hold on
    s = scatter3(Positions(:,1),Positions(:,2),Positions(:,3),'k','filled');
    s.SizeData = 2;
    view(60,30);
    daspect([1 1 1])
    plot3([Lx/2 Lx/2],[Ly/2 Ly/2],[-Lz/2 Lz/2],'r'  )
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotNetwork

global X_SP Y_SP Z_SP Rx_SP Ry_SP Rz_SP Lx Ly Lz ConnectionsSP Positions...
    Corners dims PlotChains Np N PlotStickers Nt Ns

hold on
% FontSize = 20;
% xl_per_mol = N/Np;
line_width_bd = 0.5;
line_width_bnd = 1.5;

%Plot Borders

%Plot Window
C1 = Corners(1,:);
C2 = Corners(2,:);
C3 = Corners(3,:);
C4 = Corners(4,:);
if dims==3
    C5 = Corners(5,:);
    C6 = Corners(6,:);
    C7 = Corners(7,:);
    C8 = Corners(8,:);
end
if dims==2
    plot([C1(1) C2(1)],[C1(2) C2(2)],'k-','LineWidth',line_width_bd)
    plot([C2(1) C3(1)],[C2(2) C3(2)],'k-','LineWidth',line_width_bd)
    plot([C3(1) C4(1)],[C3(2) C4(2)],'k-','LineWidth',line_width_bd)
    plot([C4(1) C1(1)],[C4(2) C1(2)],'k-','LineWidth',line_width_bd)
end
if dims==3
    plot3([C1(1) C2(1)],[C1(2) C2(2)],[C1(3) C2(3)],'k-','LineWidth',line_width_bd)
    plot3([C2(1) C3(1)],[C2(2) C3(2)],[C2(3) C3(3)],'k-','LineWidth',line_width_bd)
    plot3([C3(1) C4(1)],[C3(2) C4(2)],[C3(3) C4(3)],'k-','LineWidth',line_width_bd)
    plot3([C4(1) C1(1)],[C4(2) C1(2)],[C4(3) C1(3)],'k-','LineWidth',line_width_bd)
    plot3([C1(1) C5(1)],[C1(2) C5(2)],[C1(3) C5(3)],'k-','LineWidth',line_width_bd)
    plot3([C2(1) C6(1)],[C2(2) C6(2)],[C2(3) C6(3)],'k-','LineWidth',line_width_bd)
    plot3([C3(1) C7(1)],[C3(2) C7(2)],[C3(3) C7(3)],'k-','LineWidth',line_width_bd)
    plot3([C4(1) C8(1)],[C4(2) C8(2)],[C4(3) C8(3)],'k-','LineWidth',line_width_bd)
    plot3([C5(1) C6(1)],[C5(2) C6(2)],[C5(3) C6(3)],'k-','LineWidth',line_width_bd)
    plot3([C6(1) C7(1)],[C6(2) C7(2)],[C6(3) C7(3)],'k-','LineWidth',line_width_bd)
    plot3([C7(1) C8(1)],[C7(2) C8(2)],[C7(3) C8(3)],'k-','LineWidth',line_width_bd)
    plot3([C8(1) C5(1)],[C8(2) C5(2)],[C8(3) C5(3)],'k-','LineWidth',line_width_bd)
    view([30,45])
end

if PlotChains==1
    NStickers = Nt*Ns*Np;
    NPrimary = Nt*Np;
    if size(Positions,1)<NStickers+NPrimary
        Stickers = [];
    else
        Stickers = Positions(N+1:N+NStickers,dims+1);
    end
    Froms = ConnectionsSP(1:end,1:end-1);
    Froms = Froms(:);
    Tos = ConnectionsSP(:,end); Tos = repmat(Tos,size(ConnectionsSP,2)-1,1);
    Pairs = [Froms Tos];
    Pairs(Pairs(:,1)==0,:) = [];
    Pairs(Pairs(:,2)==0,:) = [];
    Pairs(ismember(Pairs(:,1),Stickers),:) = [];
    Pairs(ismember(Pairs(:,2),Stickers),:) = [];

    select_rng = (1:size(X_SP,1));
    select_rng(Stickers) = [];
    X1 = X_SP(select_rng,:); X1 = X1(:);
    Rx1 = Rx_SP(select_rng,:); Rx1 = Rx1(:);
    X2 = X1 + Rx1;
    Y1 = Y_SP(select_rng,:); Y1 = Y1(:);
    Ry1 = Ry_SP(select_rng,:); Ry1 = Ry1(:);
    Y2 = Y1 + Ry1;
    if dims==3
        Z1 = Z_SP(select_rng,:); Z1 = Z1(:);
        Rz1 = Rz_SP(select_rng,:); Rz1 = Rz1(:);
        Z2 = Z1 + Rz1;
    end
    Pairs(ismember(Pairs,Stickers)) = [];
    X1(Pairs(:)==0) = [];
    X2(Pairs(:)==0) = [];
    Y1(Pairs(:)==0) = [];
    Y2(Pairs(:)==0) = [];
    if dims==3
        Z1(Pairs(:)==0) = [];
        Z2(Pairs(:)==0) = [];
    end

    Color = [0.5 0.5 0.5];
    LineWidth = line_width_bnd;
    if dims==2
        plot([X1,X2]',[Y1,Y2]','Color',Color,'LineWidth',1);
    else
        plot3([X1,X2]',[Y1,Y2]',[Z1,Z2]','Color',Color,'LineWidth',LineWidth);
    end
    Color = 'k';
    LineWidth = line_width_bnd/2;
    if dims==2
        plot([X1,X2]',[Y1,Y2]','Color',Color,'LineWidth',1);
    else
        plot3([X1,X2]',[Y1,Y2]',[Z1,Z2]','Color',Color,'LineWidth',LineWidth);
    end
end

if PlotStickers==1
    NStickers = Nt*Ns*Np;
    NPrimary = Nt*Np;
    if size(Positions,1)<NStickers+NPrimary
        Stickers = [];
    else
        Stickers = Positions(N+1:N+NStickers,dims+1);
    end
    Pairs = ConnectionsSP(Stickers,1:end-1);
    Pairs = Pairs(:);
    X1 = X_SP(Stickers,:); X1 = X1(:);
    Rx1 = Rx_SP(Stickers,:); Rx1 = Rx1(:);
    X2 = X1 + Rx1;
    Y1 = Y_SP(Stickers,:); Y1 = Y1(:);
    Ry1 = Ry_SP(Stickers,:); Ry1 = Ry1(:);
    Y2 = Y1 + Ry1;
    if dims==3
        Z1 = Z_SP(Stickers,:); Z1 = Z1(:);
        Rz1 = Rz_SP(Stickers,:); Rz1 = Rz1(:);
        Z2 = Z1 + Rz1;
    end
    X1(Pairs(:)==0) = [];
    X2(Pairs(:)==0) = [];
    Y1(Pairs(:)==0) = [];
    Y2(Pairs(:)==0) = [];
    if dims==3
        Z1(Pairs(:)==0) = [];
        Z2(Pairs(:)==0) = [];
    end

    Color = 'r';
    LineWidth = line_width_bnd;
    if dims==2
        plot([X1,X2]',[Y1,Y2]',Color,'LineWidth',1);
    else
        plot3([X1,X2]',[Y1,Y2]',[Z1,Z2]',Color,'LineWidth',LineWidth);
    end
end

Fact = 1;
xlim([-Fact*Lx Fact*Lx])
ylim([-Fact*Ly Fact*Ly])
if dims==3
    ylim([-Fact*Lz Fact*Lz])
end

ShowScatter = 1;
if ShowScatter==1
    if dims==2
        s = scatter(Positions(1:N,1),Positions(1:N,2),'k','filled');
    else
        s = scatter3(Positions(1:N,1),Positions(1:N,2),Positions(1:N,3),'k','filled');
    end
    s.SizeData = 2;
end
% if PlotStickers==1
%     s = scatter3(Positions(N+1:end,1),Positions(N+1:end,2),Positions(N+1:end,3),'r','filled');
%     s.SizeData = 6;
% end



axis off

set(gcf,'Color','w')
daspect([1 1 1])
pbaspect([1 1 1])
if dims==3
    view([40,30])
    set(gca,'Projection','perspective')
end
box on

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotNetworkUnwrapped

global X_SP_Unwrapped Y_SP_Unwrapped Z_SP_Unwrapped Rx_SP Ry_SP Rz_SP...
    Lx Ly Lz ConnectionsSP PosUnwrapped...
    Corners dims PlotChains Np Nt Ns N PlotStickers

clf
hold on

%Plot Borders
line_width_bd = 1;
%Plot Window
C1 = Corners(1,:);
C2 = Corners(2,:);
C3 = Corners(3,:);
C4 = Corners(4,:);
if dims==3
    C5 = Corners(5,:);
    C6 = Corners(6,:);
    C7 = Corners(7,:);
    C8 = Corners(8,:);
end
if dims==2
    plot([C1(1) C2(1)],[C1(2) C2(2)],'k-','LineWidth',line_width_bd)
    plot([C2(1) C3(1)],[C2(2) C3(2)],'k-','LineWidth',line_width_bd)
    plot([C3(1) C4(1)],[C3(2) C4(2)],'k-','LineWidth',line_width_bd)
    plot([C4(1) C1(1)],[C4(2) C1(2)],'k-','LineWidth',line_width_bd)
end
if dims==3
    plot3([C1(1) C2(1)],[C1(2) C2(2)],[C1(3) C2(3)],'k-','LineWidth',line_width_bd)
    plot3([C2(1) C3(1)],[C2(2) C3(2)],[C2(3) C3(3)],'k-','LineWidth',line_width_bd)
    plot3([C3(1) C4(1)],[C3(2) C4(2)],[C3(3) C4(3)],'k-','LineWidth',line_width_bd)
    plot3([C4(1) C1(1)],[C4(2) C1(2)],[C4(3) C1(3)],'k-','LineWidth',line_width_bd)
    plot3([C1(1) C5(1)],[C1(2) C5(2)],[C1(3) C5(3)],'k-','LineWidth',line_width_bd)
    plot3([C2(1) C6(1)],[C2(2) C6(2)],[C2(3) C6(3)],'k-','LineWidth',line_width_bd)
    plot3([C3(1) C7(1)],[C3(2) C7(2)],[C3(3) C7(3)],'k-','LineWidth',line_width_bd)
    plot3([C4(1) C8(1)],[C4(2) C8(2)],[C4(3) C8(3)],'k-','LineWidth',line_width_bd)
    plot3([C5(1) C6(1)],[C5(2) C6(2)],[C5(3) C6(3)],'k-','LineWidth',line_width_bd)
    plot3([C6(1) C7(1)],[C6(2) C7(2)],[C6(3) C7(3)],'k-','LineWidth',line_width_bd)
    plot3([C7(1) C8(1)],[C7(2) C8(2)],[C7(3) C8(3)],'k-','LineWidth',line_width_bd)
    plot3([C8(1) C5(1)],[C8(2) C5(2)],[C8(3) C5(3)],'k-','LineWidth',line_width_bd)
    view([30,45])
end

if PlotChains==1
    NStickers = Nt*Ns*Np;
    NPrimary = Nt*Np;
    if size(PosUnwrapped,1)<NStickers+NPrimary
        Stickers = [];
    else
        Stickers = PosUnwrapped(N+1:N+NStickers,dims+1);
    end
    Froms = ConnectionsSP(1:end,1:end-1);
    Froms = Froms(:);
    Tos = ConnectionsSP(:,end); Tos = repmat(Tos,size(ConnectionsSP,2)-1,1);
    Pairs = [Froms Tos];
    Pairs(Pairs(:,1)==0,:) = [];
    Pairs(Pairs(:,2)==0,:) = [];
    Pairs(ismember(Pairs(:,1),Stickers),:) = [];
    Pairs(ismember(Pairs(:,2),Stickers),:) = [];

    select_rng = (1:size(X_SP_Unwrapped,1));
    select_rng(Stickers) = [];
    X1 = X_SP_Unwrapped(select_rng,:); X1 = X1(:);
    Rx1 = Rx_SP(select_rng,:); Rx1 = Rx1(:);
    X2 = X1 + Rx1;
    Y1 = Y_SP_Unwrapped(select_rng,:); Y1 = Y1(:);
    Ry1 = Ry_SP(select_rng,:); Ry1 = Ry1(:);
    Y2 = Y1 + Ry1;
    if dims==3
        Z1 = Z_SP_Unwrapped(select_rng,:); Z1 = Z1(:);
        Rz1 = Rz_SP(select_rng,:); Rz1 = Rz1(:);
        Z2 = Z1 + Rz1;
    end
    Pairs(ismember(Pairs,Stickers)) = [];
    X1(Pairs(:)==0) = [];
    X2(Pairs(:)==0) = [];
    Y1(Pairs(:)==0) = [];
    Y2(Pairs(:)==0) = [];
    if dims==3
        Z1(Pairs(:)==0) = [];
        Z2(Pairs(:)==0) = [];
    end

    Color = 'k';
    LineWidth = 0.5;
    if dims==2
        plot([X1,X2]',[Y1,Y2]',Color,'LineWidth',1)
    else
        plot3([X1,X2]',[Y1,Y2]',[Z1,Z2]',Color,'LineWidth',LineWidth)
    end
end

if PlotStickers==1
    NStickers = Nt*Ns*Np;
    NPrimary = Nt*Np;
    if size(PosUnwrapped,1)<NStickers+NPrimary
        Stickers = [];
    else
        Stickers = PosUnwrapped(N+1:N+NStickers,dims+1);
    end
    Pairs = ConnectionsSP(Stickers,1:end-1);
    Pairs = Pairs(:);
    X1 = X_SP_Unwrapped(Stickers,:); X1 = X1(:);
    Rx1 = Rx_SP(Stickers,:); Rx1 = Rx1(:);
    X2 = X1 + Rx1;
    Y1 = Y_SP_Unwrapped(Stickers,:); Y1 = Y1(:);
    Ry1 = Ry_SP(Stickers,:); Ry1 = Ry1(:);
    Y2 = Y1 + Ry1;
    if dims==3
        Z1 = Z_SP_Unwrapped(Stickers,:); Z1 = Z1(:);
        Rz1 = Rz_SP(Stickers,:); Rz1 = Rz1(:);
        Z2 = Z1 + Rz1;
    end
    X1(Pairs(:)==0) = [];
    X2(Pairs(:)==0) = [];
    Y1(Pairs(:)==0) = [];
    Y2(Pairs(:)==0) = [];
    if dims==3
        Z1(Pairs(:)==0) = [];
        Z2(Pairs(:)==0) = [];
    end

    Color = 'r';
    LineWidth = 0.5;
    if dims==2
        plot([X1,X2]',[Y1,Y2]','Color',Color,'LineWidth',1);
    else
        plot3([X1,X2]',[Y1,Y2]',[Z1,Z2]','Color',Color,'LineWidth',LineWidth);
    end
end

Fact = 1;
xlim([-Fact*Lx Fact*Lx])
ylim([-Fact*Ly Fact*Ly])
if dims==3
    ylim([-Fact*Lz Fact*Lz])
end

ShowScatter = 1;
if ShowScatter==1
    if dims==2
        s = scatter(PosUnwrapped(1:N,1),PosUnwrapped(1:N,2),'k','filled');
    else
        s = scatter3(PosUnwrapped(1:N,1),PosUnwrapped(1:N,2),PosUnwrapped(1:N,3),'k','filled');
    end
    s.SizeData = 2;
end

axis off

set(gcf,'Color','w')
daspect([1 1 1])
pbaspect([1 1 1])
if dims==3
    view([30,30])
end
box on

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function norms_Kuhn = DefineKuhnLengths

global X_SP Y_SP Z_SP Rx_SP Ry_SP Rz_SP...
    ConnectionsSP Positions dims Np Nt Ns N

NStickers = Nt*Ns*Np;
NPrimary = Nt*Np;
if size(Positions,1)<NStickers+NPrimary
    Stickers = [];
else
    Stickers = Positions(N+1:N+NStickers,dims+1);
end
Pairs = ConnectionsSP(Stickers,1:end-1);
Pairs = Pairs(:);
X1 = X_SP(Stickers,:); X1 = X1(:);
Rx1 = Rx_SP(Stickers,:); Rx1 = Rx1(:);
X2 = X1 + Rx1;
Y1 = Y_SP(Stickers,:); Y1 = Y1(:);
Ry1 = Ry_SP(Stickers,:); Ry1 = Ry1(:);
Y2 = Y1 + Ry1;
if dims==3
    Z1 = Z_SP(Stickers,:); Z1 = Z1(:);
    Rz1 = Rz_SP(Stickers,:); Rz1 = Rz1(:);
    Z2 = Z1 + Rz1;
end
X1(Pairs(:)==0) = [];
X2(Pairs(:)==0) = [];
Y1(Pairs(:)==0) = [];
Y2(Pairs(:)==0) = [];
Z1(Pairs(:)==0) = [];
Z2(Pairs(:)==0) = [];

vec = [X1-X2,Y1-Y2,Z1-Z2];
norms_Kuhn = vecnorm(vec,2,2);
norms_Kuhn(norms_Kuhn==0) = [];
norms_Kuhn(isnan(norms_Kuhn)) = [];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function movie2gif(mov, gifFile, varargin)
% ==================
% Matlab movie to GIF Converter.
%
% Syntax: movie2gif(mov, gifFile, prop, value, ...)
% =================================================
% The list of properties is the same like for the command 'imwrite' for the
% file format gif:
%
% 'BackgroundColor' - A scalar integer. This value specifies which index in
%                     the colormap should be treated as the transparent
%                     color for the image and is used for certain disposal
%                     methods in animated GIFs. If X is uint8 or logical,
%                     then indexing starts at 0. If X is double, then
%                     indexing starts at 1.
%
% 'Comment' - A string or cell array of strings containing a comment to be
%             added to the image. For a cell array of strings, a carriage
%             return is added after each row.
%
% 'DelayTime' - A scalar value between 0 and 655 inclusive, that specifies
%               the delay in seconds before displaying the next image.
%
% 'DisposalMethod' - One of the following strings, which sets the disposal
%                    method of an animated GIF: 'leaveInPlace',
%                    'restoreBG', 'restorePrevious', or 'doNotSpecify'.
%
% 'LoopCount' - A finite integer between 0 and 65535 or the value Inf (the
%               default) which specifies the number of times to repeat the
%               animation. By default, the animation loops continuously.
%               For a value of 0, the animation will be played once. For a
%               value of 1, the animation will be played twice, etc.
%
% 'TransparentColor' - A scalar integer. This value specifies which index
%                      in the colormap should be treated as the transparent
%                      color for the image. If X is uint8 or logical, then
%                      indexing starts at 0. If X is double, then indexing
%                      starts at 1
%
% *************************************************************************
% Copyright 2007-2013 by Nicolae Cindea.
if (nargin < 2)
    error('Too few input arguments');
end
if (nargin == 2)
    frameNb = size(mov, 2);
    isFirst = true;
    h = waitbar(0, 'Generate GIF file...');
    for i = 1:frameNb
        waitbar((i-1)/frameNb, h);
        [RGB, ~] = frame2im(mov(i));
        if (exist('rgb2ind', 'file'))
            [IND, map] = rgb2ind(RGB,256);
        else
            [IND, map] = aRGB2IND(RGB);
        end
        if isFirst
            imwrite(IND, map, gifFile, 'gif');
            isFirst=false;
        else
            imwrite(IND, map, gifFile, 'gif', 'WriteMode', 'append');
        end
    end
    close(h);
end
if (nargin > 2)
    h = waitbar(0, 'Generate GIF file...');
    frameNb = size(mov, 2);
    isFirst = true;
    for i = 1:frameNb
        waitbar((i-1)/frameNb, h);
        [RGB, ~] = frame2im(mov(i));
        if (exist('rgb2ind', 'file'))
            [IND, map] = rgb2ind(RGB,256);
        else
            [IND, map] = aRGB2IND(RGB);
        end
        if isFirst
            args = varargin;
            imwrite(IND, map, gifFile, 'gif', args{:});
            isFirst=false;

            % supress 'LoopCount' option from the args!!
            args = varargin;
            l = length(args);

            posLoopCount = 0;
            for ii = 1:l
                if(ischar(args{ii}))
                    if strcmp(args{ii}, 'LoopCount')
                        posLoopCount = ii;
                    end
                end
            end
            if (posLoopCount)
                args = {args{1:posLoopCount-1}, args{posLoopCount+2:end}};
            end

        else
            imwrite(IND, map, gifFile, 'gif', 'WriteMode', 'append', ...
                args{:});
        end
    end
    close(h);
end
end