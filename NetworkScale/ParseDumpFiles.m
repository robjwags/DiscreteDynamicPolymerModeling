function ParseDumpFiles(afile,bfile,dt,tau0,lambda,Weissenberg,kd0,ka0,Controls)

%% Parsing atoms

global raw_data_filename BeadSpringOrMeso 

% Define all pertinent timescales, output intervals, etc.
[dt_eq,dt_load,dt_relax,neq,nload,~,nmax,~,~,~,~,~,~,iout_eq,~,~,~,~,~,...
    ~,~,~,~,~,~,~] =...
    DefineInputTimescales(tau0,Weissenberg,...
    kd0,ka0,lambda,dt,BeadSpringOrMeso,Controls);

n_data = floor(nmax/iout_eq);
fid = fopen(afile);

% ID cutoffs for each regime and each dt
pcnt_eq = neq/nmax;
pcnt_ld = (neq+nload)/nmax;

ii = 0; 

Corners = zeros(n_data,7);
timesteps = zeros(n_data,1);
time = zeros(n_data,1);

wb = waitbar(0,'Parsing atom data...');
while 1 == 1

    ii = ii + 1;

    waitbar(ii/n_data,wb,'Parsing atom data...');

    % Skip first line
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end

    % Timestep
    tstep_temp = str2double(fgetl(fid));

    % Skip
    tline = fgetl(fid);

    % Number of atoms in the system
    natoms(ii) = str2double(fgetl(fid));

    % Skip header
    tline = fgetl(fid);

    % Simulation dimensions
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

    if ii/n_data<pcnt_eq
        dt_temp = dt_eq;
    elseif ii/n_data>=pcnt_eq && ii/n_data<pcnt_ld
        dt_temp = dt_load;
    else
        dt_temp = dt_relax;
    end
    if ii==1
        time(ii,1) = 0;
    else
        nsteps = timesteps(ii,1) - timesteps(ii-1,1);
        time(ii,1) = time(ii-1,1)+nsteps*dt_temp;
    end
end
close(wb)
fclose(fid);

tic
rem_indx = find(timesteps(2:end)==0,1,'first');
timesteps(rem_indx:end) = [];
Corners(rem_indx:end,:) = [];
RawAtomsData.timesteps = timesteps;
RawAtomsData.Corners = Corners;
save(raw_data_filename,'-struct','RawAtomsData','-v7.3')
toc

save (raw_data_filename,'time','-append');
save (raw_data_filename,'id','-append');
save (raw_data_filename,'x','-append');
save (raw_data_filename,'y','-append');
save (raw_data_filename,'z','-append');
save (raw_data_filename,'type','-append');
save (raw_data_filename,'mol','-append');

show_fig_debug = 1;
if show_fig_debug==1
    figure(1e3); clf; hold on
    plot(time);
    close
end

clear timesteps Corners id x y z vx vy vz type mol

%% Parsing bonds
fid = fopen(bfile);

wb = waitbar(0,'Parsing bond data...');
nbonds = zeros(n_data,1); tbstep = zeros(n_data,1);
ii = 0; 
while 1 == 1

    ii = ii + 1;
    waitbar(ii/n_data,wb,'Parsing bond data...');

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

save (raw_data_filename,'btype','-append');
save (raw_data_filename,'bid','-append');
save (raw_data_filename,'batom1','-append');
save (raw_data_filename,'batom2','-append');
save (raw_data_filename,'rx','-append');
save (raw_data_filename,'ry','-append');
save (raw_data_filename,'rz','-append');
save (raw_data_filename,'bForce','-append');

clear btype bid batom1  batom2 rx ry rz bForce

end
