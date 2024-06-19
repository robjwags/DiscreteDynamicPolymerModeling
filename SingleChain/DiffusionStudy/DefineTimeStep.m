function [dts,damps,D,dtFact] = DefineTimeStep(b,LengthConversion,DamperConversion,...
    BeadSpringOrMeso)

b_SI = b*LengthConversion;
D_nom = 1e-10;          % nominal diffusion coefficient of a monomer in m^2/s
D = D_nom*10.^[-2 0] ;         % m^2/s
tau0 = (b_SI^2)./D;    % s
kbT = 293*1.38e-23;     % Thermal energy in Joules
damps = kbT./D;                             % kg/s or N s/m
damps = damps/DamperConversion;             % this is for a single monomerq (effective damper will depend on N)
if BeadSpringOrMeso==0
    dtFact = 40;%320;
elseif BeadSpringOrMeso==1
    dtFact = 20;
end
dts = 1/dtFact*min(tau0);             % NOTE THAT Diffusion Coeff, damper, and timestep are all linked
D = D(2); damps = damps(2);

end