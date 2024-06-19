function UnpackParameters(package,Parameters)

global length_conversion damper_conversion energy_conversion...
    Np Nt Ns force_conversion...
    kbT b lambda ka...
    samples eq_time_factors N_Kuhns phis kds Weissenbergs...
    kb T...
    omega_min omega_max tf omegas

force_conversion = Parameters.force_conversion;
length_conversion = Parameters.length_conversion;
damper_conversion = Parameters.damper_conversion;
Np = Parameters.Np;
Nt = Parameters.Nt;
Ns = Parameters.Ns;

kbT = Parameters.kbT;
b = Parameters.b;
lambda = Parameters.lambda;
ka = Parameters.ka;
energy_conversion = 1/kbT;

% select parameters in arbitrary units
kb = 1;
T = 1;

% Swept parameters
samples = unique(package(:,1));
eq_time_factors = unique(package(:,2));
N_Kuhns = unique(package(:,3));
phis = unique(package(:,4));
kds = unique(package(:,5));
Weissenbergs = unique(package(:,6));

if isfield(Parameters,'omega_max')
    omega_max = Parameters.omega_max;
    omega_min = Parameters.omega_min;
    tf = Parameters.tf;
    omegas = unique(package(:,9));
end

end