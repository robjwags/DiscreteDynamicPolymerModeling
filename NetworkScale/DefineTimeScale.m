function [damps,D,tau0,dtFact] = DefineTimeScale(b,length_conversion,...
    damper_conversion,BeadSpringOrMeso,phi,N_Kuhn)

b_SI = b*length_conversion;
D_nom = 1e-10;          % nominal diffusion coefficient of a monomer 
                        % in [m^2/s]
D = D_nom*10.^[-2 0] ;  % diffusion coefficient [m^2/s]
tau0 = (b_SI^2)./D;     % singler mer diffusion timescale [s]
kbT = 293*1.38e-23;     % Thermal energy [J]
damps = kbT./D;         % Overdamped Langevin eq. of motion damping [kg/s] 
                        % or [N s/m]
damps = damps/damper_conversion; % this is for a single monomer. Effective 
                        % damper for mesoscale also depends on N.

% Define factor by which timestep is smaller than tau0 (i.e., dt =
% tau0/dtFact)
if BeadSpringOrMeso==0
    % dt will be set to tau0/dtFact => larger dtFact means lower dt
    % Ramp dtFact up as function of N_Kuhn based on trial and error and
    % iniital network stability
    if N_Kuhn<=12
        dtFact = 4e4;
    elseif 12<N_Kuhn && N_Kuhn<30
        dtFact = 79083*log(N_Kuhn)-156514;
    elseif N_Kuhn>=30
        dtFact = 12e4;
    end
    if phi>=0.5
        dtFact = dtFact*1.2;  % for higherg packing density, decrease dt by factor of 2
    end
elseif BeadSpringOrMeso==1
    dtFact = 4e3;
end
% Use only one set of the parameters based on reasonable physical values
D = D(2); damps = damps(2); tau0 = min(tau0);

end