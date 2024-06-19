function PlotParameterSpace

close all
clc
clear

% Topological space
phis = [0.204 0.5];
N_Kuhns = [12 18 24 30 36];

figure(1); clf; hold on
set(gca,'FontSize',20/1.5)
xlabel('$N$','FontSize',20,'Interpreter','latex')
ylabel('$\phi$','FontSize',20,'Interpreter','latex')
pbaspect([1 1 1])
xticks(N_Kuhns)
yticks(phis)
xlim([11 37])
ylim([0.1 0.6])


% Kinetic space
eps_dots = logspace(-2,-1,5);
kds = [0 0.001 0.0032 0.01 0.032 0.1];

figure(2); clf; hold on
set(gca,'FontSize',20/1.5)
xlabel('$\dot \varepsilon \tau_0$','FontSize',20,'Interpreter','latex')
ylabel('$k_d \tau_0$','FontSize',20,'Interpreter','latex')
pbaspect([1 1 1])
xticks(eps_dots)
yticks(kds)
set(gca,'yscale','log')
set(gca,'xscale','log')
xlim([min(eps_dots),max(eps_dots)])
ylim([min(kds),max(kds)])

x = (linspace(min(eps_dots),max(eps_dots),100))';
y = 2*x;
plot(x,y,'k--')
xlim([0.01 0.1])
ylim([0.01 0.1])


end