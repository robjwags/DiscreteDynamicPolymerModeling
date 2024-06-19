function [dt_eq,dt_ld,dt_rlx,...
    neq,nld,nrlx,nmax,ttot,tdyn,...
    i_dyn_eq,i_dyn_ld,i_dyn_rlx,...
    N_data,...
    iout_eq,iout_ld,iout_rlx,...
    ithermo_eq,ithermo_load,ithermo_rlx,...
    lamdot,taukd,tauka,tauW,teq,tload,trelax] =...
    DefineInputTimescales(tau0,W,kd0,ka0,lambda,dt,...
    BeadSpringOrMeso,Controls)

global omega_max tf omega Np


% Define repulsive properties, if any properties
[lamdot,taukd,tauka,tauW,teq,tload,trelax] = ...
    DefineLoadingTimescale(tau0,W,kd0,ka0,lambda);

% Nominal timestep - that of the equilibrium phases
dt_eq = dt;
dt_rlx = dt;
if Controls.RunLargeDeformation==1 
    dt_ld = dt;
    if W==0.0005
        dt_ld = dt/3;
    elseif W==0.000125
        dt_ld = 2*dt;
    end
elseif Controls.RunOscillatory==1
    dt_ld = dt;
else
    if BeadSpringOrMeso==0
        dt_ld = dt;
    elseif BeadSpringOrMeso==1
        denom = 100*(W-0.01)+1;
        dt_ld = dt/denom;
    end
end

if Controls.RunOscillatory==1
    if omega<0.0002
        tload = tau0/omega;
    elseif omega<0.0003
        tload = 2*tau0/omega;
    elseif omega<0.0004
        tload = 3*tau0/omega;
    else
        tload = tf/omega*tau0;
    end
end

% Number of iterations for each phase
neq = ceil(teq/dt_eq);
nld = ceil(tload/dt_ld);
nrlx = ceil(trelax/dt_rlx);
nmax = neq+nld+nrlx;

if Controls.RunTimingCase==1 && Np>500
    dt_ld = dt_ld/4;
end

ttot = neq*dt_eq + nld*dt_ld + nrlx*dt_rlx;

if Controls.RunOscillatory==1
    nmax = nld;
    ttot = tf;
end

% Define trial frequency for dynamics
tdyn = tau0/20;
i_dyn_eq = round(tdyn/dt_eq);
i_dyn_ld = round(tdyn/dt_ld);
i_dyn_rlx = round(tdyn/dt_rlx);

% Define dump output interval based on desired eq time between points
if BeadSpringOrMeso==0
    tdata = 2*tau0;
elseif BeadSpringOrMeso==1
    tdata = tau0/2;
end
iout_eq = round(tdata/dt_eq);
iout_ld = round(tdata/dt_ld);
iout_rlx = round(tdata/dt_rlx);

% If doing long experiment, set iout to something less frequent
if Controls.RunLargeDeformation==1
    iout_eq = round(nmax/2000);
elseif Controls.RunOscillatory==1
    T_min = 1/omega_max*tau0;
    n_min_per_period = round(T_min/dt_ld);
    iout_ld = n_min_per_period/20;
    iout_eq = 5*iout_eq;
end

% Define total number of output data points
N_data_eq = floor(neq/iout_eq);
N_data_ld = floor(nld/iout_ld);
N_data_rlx = floor(nrlx/iout_rlx);
N_data = N_data_eq + N_data_ld + N_data_rlx;

if Controls.RunOscillatory==1
    N_data = N_data_eq + N_data_ld;
end

% Consol output interval
tthermo = tau0/20;
ithermo_eq = round(tthermo/dt_eq);
ithermo_load = round(tthermo/dt_ld);
ithermo_rlx = round(tthermo/dt_rlx);

end