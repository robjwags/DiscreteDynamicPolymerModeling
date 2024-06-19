function Parameters = DefineparameterPackage(Controls)

% Returns structure "Parameters" containing pertinent constants as well as
% the sweeping parameters, which are stored in an array called "package"
% for which each row stores the sweeping parameter values for a given
% sample, pre-stretch equilibration time, and parameter combination.

global length_conversion damper_conversion...
    force_conversion energy_conversion

% Unit conversions
length_conversion = 3.74e-9;    % meters per unit code length
damper_conversion = 2.89e-4;    % N*s/m per unit damper in code
force_conversion = 1.08e-12;    % N per unit force in code

% Define Input Constants
Np = 60;                % Number of polymers (i.e., molecules)
Nt = 5;                 % Number of side chains per molecule
Ns = 1;                 % Number of free stickers perg grafting site
kbT = 293*1.38e-23;     % Thermal energy in Joules
ea = 0.01*kbT;          % Normalized bond activation energy for association  
b = 0.1667;             % Kuhn length
lambda = 3;             % Extension (200% strain)
energy_conversion = kbT;    % conversion from code units to J;

%% Packae the sweeping parameters
% Define Sweeping Parameters
% Topological parameters
N_Kuhn = [12 36];                   % Number of Kuhn segments in 
phi = [0.2043 0.5];                 % packing volume of polymer assuming 
                                    % v_chain ~ Nb3

% Bond kinetics and loading rate parameters
[~,~,tau0,~] = ...                  % Define timescales
    DefineTimeScale(b,length_conversion,damper_conversion,1,1,1); 
eps_dot = [1/100 1/10];         % Sets loading rate (higher value => 
                                    % higher rate)
ka = 1/tau0*exp(-ea/kbT);
kd = [0 ka/10];

% Sampling
samples = 1:5;                      % Sample IDs
eq_time_factors = [1 1.25 1.5];     % equilibration times as factors of the 
                                    % nominal (1/ka0*160);
% total numer of samples for a given combo. of parameters is therefore
% length(Samples)*length(eq_time_factors) (=9 by default)

if Controls.RunTheLJCase==1
    N_Kuhn = 12;
    kd = 0;
    eps_dot = 0.1;
end

if Controls.RunChainLengthSweep==1
%     N_Kuhn = [12,13,14,15,16,18,24,30,36];
    N_Kuhn = [12,14,16,18,24,30,36];
    phi = phi(1);
    kd = 0;
    eps_dot = 0.01;
elseif Controls.RunLoadingRateSweep==1
    N_Kuhn = 12;
    phi = phi(1);
%     kd = ka/10;
    kd = 0;
    eps_dot = logspace(-2,-1,5);
elseif Controls.RunDetachmentRateSweep==1
    N_Kuhn = 12;
    phi = phi(1);
    kd = [0, [0.001 0.0032 0.01 0.032]/tau0, ka/10];
    eps_dot = 0.01;
elseif Controls.RunOscillatory==1
% 
% % %     batch = 4;
% % %     if batch==1
%         % Batch 1
%         samples = 1:20;
%         lambda = 1.1;      % strain amplitude is lambda-1
%         N_Kuhn = 12;
%         phi = phi(1);
%         kd = 0.002/tau0;
%         eps_dot = 0.001;
% %         omega_min = 0.001; % tau0^(-1)
%         omega_min = 0.0004; % tau0^(-1)
%         omega_max = 0.01;    % tau0^(-1)
% %         tf = 10/0.0008;
%         tf = 10/0.0004;
%         Np = 100;%1e3;
%         eq_time_factors = 1;
% % %     elseif batch==2
% % %         % Batch 2
% % %         samples = 1:6;
% % %         lambda = 1.1;      % strain amplitude is lambda-1
% % %         N_Kuhn = 12;
% % %         phi = phi(1);
% % %         kd = 0.002/tau0;
% % %         eps_dot = 0.001;
% % %         omega_min = 0.0004;  % tau0^(-1)
% % %         omega_max = 0.01;   % tau0^(-1)
% % %         tf = 10/omega_min;
% % %         Np = 100;
% % %         eq_time_factors = 1;
% % %     elseif batch==3
% % %         % Batch 3
% % %         samples = 1:15;
% % %         lambda = 1.1;      % strain amplitude is lambda-1
% % %         N_Kuhn = 12;
% % %         phi = phi(1);
% % %         kd = 0.002/tau0;
% % %         eps_dot = 0.001;
% % %         omega_min = 0.0005;  % tau0^(-1)
% % %         omega_max = 0.01;   % tau0^(-1)
% % %         tf = 10/omega_min;
% % %         Np = 100;
% % %         eq_time_factors = 1;
% % %     elseif batch==4
% % %         % Batch 4
% % %         samples = 1:15;
% % %         lambda = 1.1;      % strain amplitude is lambda-1
% % %         N_Kuhn = 12;
% % %         phi = phi(1);
% % %         kd = 0.002/tau0;
% % %         eps_dot = 0.001;
% % %         omega_min = 0.0002;  % tau0^(-1)
% % %         omega_max = 0.01;   % tau0^(-1)
% % %         tf = 10/omega_min;
% % %         Np = 100;%60;
% % %         eq_time_factors = 1;
% % %     end
% 
% % % % % %     batch = 2;
% % % % % %     samples = 1:12;
% % % % % %     lambda = 1.1;      % strain amplitude is lambda-1
% % % % % %     N_Kuhn = 12;
% % % % % %     phi = phi(1);
% % % % % %     kd = 0.002/tau0;%0.0001/tau0;
% % % % % %     eps_dot = 0.001;
% % % % % %     eq_time_factors = 1;
% % % % % %     Np = 100;%1e3;
% % % % % %     if batch==1
% % % % % %         % Batch 1
% % % % % %         omega_min = 0.0004;     % tau0^(-1)
% % % % % %         omega_max = 0.001;      % tau0^(-1)
% % % % % % % %         tf = 5/omega_min;
% % % % % %     elseif batch==2
% % % % % %         % Batch 2
% % % % % %         omega_min = 0.0007;      % tau0^(-1)
% % % % % %         omega_max = 0.008;      % tau0^(-1)
% % % % % % %         tf = 6/omega_min;
% % % % % %     elseif batch==3
% % % % % %         % Batch 3
% % % % % %         omega_min = 0.005;      % tau0^(-1)
% % % % % %         omega_max = 0.01;       % tau0^(-1)
% % % % % % %         tf = 10/omega_min;
% % % % % %     elseif batch==4
% % % % % % %         % Batch 4
% % % % % % %         omega_min = 0.007;       % tau0^(-1)
% % % % % % %         omega_max = 0.04;       % tau0^(-1)
% % % % % % % %         tf = 15/omega_min;
% % % % % % %     elseif batch==5 % linear ramp in omega
% % % % % % %         omega_min = 0.0001;
% % % % % % %         omega_max = 0.01;
% % % % % % %         kd = 0.001/tau0;
% % % % % % %         % Omega spans from omega_min to omega_max over t_f as
% % % % % % %             % omega = (omega_max-omega_min)*t/t_f + omega_min
% % % % % %     end
% % % % % %     tf = 5/0.0004;

    
    samples = 1:15;
    lambda = 1.1;      % strain amplitude is lambda-1
    N_Kuhn = 12;
    phi = phi(1);
    kd = [1e-4 1e-3 1e-2]/tau0;
%     kd = 1e-2/tau0;
%     kd = [1e-4 1e-3]/tau0;
    eps_dot = 0.001;
    omega_min = 0.0001; % tau0^(-1)
    omega_max = 0.01;    % tau0^(-1)
    omegas = logspace(-4,-2,15);
    % Add a lower frequency
    omegas = [5.18e-5,7.2e-5,omegas];

    tf = 5;  % max run time will be tf/omega
    Np = 100;
    eq_time_factors = 1;

    Parameters.omega_min = omega_min;
    Parameters.omega_max = omega_max;
    Parameters.tf = tf;
elseif Controls.RunLargeDeformation==1
    samples = 1;
	lambda = 20;
    N_Kuhn = 12;
    phi = phi(1);
    kd = 0.001/tau0;
%     eps_dot = 0.00025;
    eps_dot = [0.000125,0.00025,0.0005,0.001];
    Np = 5e3;
    eq_time_factors = [1 1.25 1.5];
elseif Controls.RunTimingCase==1
    N_Kuhn = 12;
    samples = 1;
    Np = 60;
    phi = phi(1);
    kd = 0.01/tau0;
    eps_dot = 0.01;
    eq_time_factors = 1;
end

if Controls.RunOscillatory~=1
    % Package the Input Parameters
    n_params = 8;
    n_perms = length(samples)*length(eq_time_factors)*length(N_Kuhn)...
        *length(phi)*length(kd)*length(eps_dot);
    package = zeros(n_perms,n_params); ct = 0;
    for i=1:length(samples)
        for j=1:length(eq_time_factors)
            for k=1:length(N_Kuhn)
                for l=1:length(phi)
                    for m=1:length(kd)
                        for n=1:length(eps_dot)
                            ct = ct+1;

                            % Assemple the parameter sweep package
                            package(ct,1) = samples(i);
                            package(ct,2) = eq_time_factors(j);
                            package(ct,3) = N_Kuhn(k);
                            package(ct,4) = phi(l);
                            package(ct,5) = kd(m);
                            package(ct,6) = eps_dot(n);

                            % Check number of particles
                            N_backbone_beads = Nt+(Nt-1)*(N_Kuhn(k)-1);
                            N_sidechain_beads = Nt*N_Kuhn(k)*Ns;
                            N_bead = Np*(N_backbone_beads+N_sidechain_beads);
                            N_meso = Np*(Nt*Ns+Nt);

                            package(ct,7) = N_bead;
                            package(ct,8) = N_meso;
                        end
                    end
                end
            end
        end
    end
else
    % Package the Input Parameters
    n_params = 9;
    n_perms = length(samples)*length(eq_time_factors)*length(N_Kuhn)...
        *length(phi)*length(kd)*length(eps_dot)*length(omegas);
    package = zeros(n_perms,n_params); ct = 0;
    for i=1:length(samples)
        for j=1:length(eq_time_factors)
            for k=1:length(N_Kuhn)
                for l=1:length(phi)
                    for m=1:length(kd)
                        for n=1:length(eps_dot)
                            for o=1:length(omegas)
                                ct = ct+1;

                                % Assemple the parameter sweep package
                                package(ct,1) = samples(i);
                                package(ct,2) = eq_time_factors(j);
                                package(ct,3) = N_Kuhn(k);
                                package(ct,4) = phi(l);
                                package(ct,5) = kd(m);
                                package(ct,6) = eps_dot(n);

                                % Check number of particles
                                N_backbone_beads = Nt+(Nt-1)*(N_Kuhn(k)-1);
                                N_sidechain_beads = Nt*N_Kuhn(k)*Ns;
                                N_bead = Np*(N_backbone_beads+N_sidechain_beads);
                                N_meso = Np*(Nt*Ns+Nt);

                                package(ct,7) = N_bead;
                                package(ct,8) = N_meso;

                                package(ct,9) = omegas(o);
                            end
                        end
                    end
                end
            end
        end
    end
end

Parameters.length_conversion = length_conversion;
Parameters.damper_conversion = damper_conversion;
Parameters.force_conversion = force_conversion;
Parameters.energy_conversion = energy_conversion;
Parameters.Np = Np;
Parameters.Nt = Nt;
Parameters.Ns = Ns;
Parameters.kbT = kbT;
Parameters.ea = ea;
Parameters.b = b;
Parameters.lambda = lambda;
Parameters.tau0 = tau0;
Parameters.ka = ka;
Parameters.package = package;


end