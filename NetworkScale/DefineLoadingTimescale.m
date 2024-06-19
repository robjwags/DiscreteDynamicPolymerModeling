function [lamdot,tau_d,tau_a,tauW,teq,tload,trelax] = ...
    DefineLoadingTimescale(tau0,Weiss,kd0,ka0,lambda)

global eq_time_factor

lamdot = Weiss/tau0;            % loading rate
tau_d = 1/kd0;              % characteristic dissociation time
tau_a = 1/ka0;              % characteristic association time
tauW = 1/lamdot;

tau_eq_nom = 200*tau_a;     % nominal equilibration time (200x associative 
                            % timescale)

teq = eq_time_factor*tau_eq_nom; % equilibration time before loading
tload = log(lambda)/lamdot; % loading time
trelax = 2.2e2*tau0;        % post-loading relaxationt time

end