function ParseDumpFiles(afile,bfile,dt,dtFact)

%% Parsing atoms

global RawDataFileName BeadSpringOrMeso

tau0 = dtFact*dt;
if BeadSpringOrMeso==0 || BeadSpringOrMeso==1
    N_dyn = 1;
elseif BeadSpringOrMeso==2
    N_dyn = dtFact;
end
theal = 2e3*dt*dtFact;
tload = 0;
trelax = 0;
nheal = ceil(theal/dt);
nload = ceil(tload/dt);
nrelax = ceil(trelax/dt);
nmax = nheal+nload+nrelax;
ttot = nmax*dt;

% tdata = dt*dtFact/10*N_dyn;  %Basically measures with frequency matching 1/10x the diffusion timescale
% NData = ceil(ttot/tdata);
tdata = tau0/10;
NData = round(ttot/tdata);       % Number of outputs data points - target ~1k

fid = fopen(afile);

ii = 0; Atoms = [];
Corners = zeros(NData,7);
timesteps = zeros(NData,1);
wb = waitbar(0,'Parsing atom data...');
while 1 == 1

    ii = ii + 1;

    waitbar(ii/NData,wb,'Parsing atom data...');

    % Skip first line
    if fid==-1
       fid; 
    end
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end

    % Timestep
%     tstep(ii) = str2double(fgetl(fid));
    tstep_temp = str2double(fgetl(fid));

    % Skip
    tline = fgetl(fid);

    % Number of atoms in the system
    natoms(ii) = str2double(fgetl(fid));

    % Skip header
    tline = fgetl(fid);

    % Simulation dimensions
%     xlims{ii} = str2double(strsplit(fgetl(fid)));
%     ylims{ii} = str2double(strsplit(fgetl(fid)));
%     zlims{ii} = str2double(strsplit(fgetl(fid)));
    xlims_temp = str2double(strsplit(fgetl(fid)));
    ylims_temp = str2double(strsplit(fgetl(fid)));
    zlims_temp = str2double(strsplit(fgetl(fid)));

    Corners(ii,1) = tstep_temp;
    Corners(ii,2:3) = xlims_temp;
    Corners(ii,4:5) = ylims_temp;
    Corners(ii,6:7) = zlims_temp;

    % Column legend
    columns = fgetl(fid);

    % Read the rest of the file 
    C = textscan(fid,'%f %f %f %f %f %f %f %f %f');

    id{ii}   = C{1};
    x{ii}    = C{2};
    y{ii}    = C{3};
    z{ii}    = C{4};
    vx{ii}   = C{5};
    vy{ii}   = C{6};
    vz{ii}   = C{7};
    type{ii} = C{8};
    mol{ii}  = C{9};
    timesteps(ii,1) = tstep_temp;

%     id_temp   = C{1};
%     x_temp    = C{2};
%     y_temp    = C{3};
%     z_temp    = C{4};
%     vx_temp   = C{5};
%     vy_temp   = C{6};
%     vz_temp   = C{7};
%     type_temp = C{8};
%     mol_temp  = C{9};
%     timestep_temp = tstep_temp*ones(size(id_temp));
% 
%     Atoms = [Atoms;[timestep_temp,id_temp,...
%         x_temp,y_temp,z_temp,...
%         vx_temp,vy_temp,vz_temp,...
%         type_temp,mol_temp]];
end
Corners(timesteps(2:end)==0,:) = [];
timesteps(timesteps(2:end)==0) = [];

close(wb)
fclose(fid);

tic
RawAtomsData.timesteps = timesteps;
RawAtomsData.Corners = Corners;
save(RawDataFileName,'-struct','RawAtomsData','-v7.3')
toc

save (RawDataFileName,'id','-append');
save (RawDataFileName,'x','-append');
save (RawDataFileName,'y','-append');
save (RawDataFileName,'z','-append');
save (RawDataFileName,'type','-append');
save (RawDataFileName,'mol','-append');

clear timesteps Corners id x y z vx vy vz type mol

%% Parsing bonds
fid = fopen(bfile);

wb = waitbar(0,'Parsing bond data...');
nbonds = zeros(NData,1); tbstep = zeros(NData,1);
ii = 0; Bonds = [];
while 1 == 1

    ii = ii + 1;
    waitbar(ii/NData,wb,'Parsing bond data...');

    % Skip first line
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end

    % Timestep
    tbstep(ii,1) = str2double(fgetl(fid));

    % Skip
    tline = fgetl(fid);

    % Number of atoms in the system
    nbonds(ii,1) = str2double(fgetl(fid));
    if nbonds(ii) == 0
        % Skip
        tline = fgetl(fid);
        tline = fgetl(fid);
        tline = fgetl(fid);
        tline = fgetl(fid);
        tline = fgetl(fid);
        continue
    end

    % Skip header
    tline = fgetl(fid);

    % Simulation dimensions
    xlimsb{ii} = str2double(strsplit(fgetl(fid)));
    ylimsb{ii} = str2double(strsplit(fgetl(fid)));
    zlimsb{ii} = str2double(strsplit(fgetl(fid)));

    % Column legend
    columns = fgetl(fid);

    % Read the rest of the file 
    C = textscan(fid,'%f %f %f %f %f %f %f %f %f');

    bid{ii}   = C{1};
    btype{ii}   = C{2};
    batom1{ii}  = C{3};
    batom2{ii}  = C{4};
    rx{ii}      = C{5};
    ry{ii}      = C{6};
    rz{ii}      = C{7};
    bLength{ii} = C{8};
    bForce{ii}  = C{9};
end
close(wb)
fclose(fid);

save (RawDataFileName,'btype','-append');
save (RawDataFileName,'bid','-append');
save (RawDataFileName,'batom1','-append');
save (RawDataFileName,'batom2','-append');
save (RawDataFileName,'rx','-append');
save (RawDataFileName,'ry','-append');
save (RawDataFileName,'rz','-append');
save (RawDataFileName,'bForce','-append');

clear btype bid batom1  batom2 rx ry rz bForce

end
