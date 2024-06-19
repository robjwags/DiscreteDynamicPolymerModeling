% Compute number of neighbors and their average distance as functions of
% separation distance inside of a grid

function ComputeNandDbar

global NoFitAttempts PlotDuringFitting plot_ct

PlotDuringFitting = 1;
NoFitAttempts = 50;
plot_ct = 1;
Size = 400;
npts = 100;
d0 = (logspace(-2,0,npts))';

n_neighb = zeros(length(d0),1);
d_bar = zeros(length(d0),1);

for i=1:length(d0)
    d0_temp = d0(i);

    vec = (0:d0_temp:1)';
    vec = unique(sort([-vec;vec]));
    [X,Y,Z] = meshgrid(vec,vec,vec);

    d_all = (X.^2 + Y.^2 + Z.^2).^0.5;
    neighbs = find(d_all<=1);
    n_neighb(i) = length(neighbs);
    d_bar(i) = mean(d_all(neighbs));
end


figure(1); clf; hold on

x = d0;
y = n_neighb;
err = zeros(size(n_neighb));
A_nom = 4.4339; B_nom = -3; C_nom = 0;
[A,B,~,R2,~] = FitFixedPwr(A_nom,B_nom,C_nom,x,y,err,1);


plot(d0,A*x.^B,'r--','LineWidth',1.5);

p = plot(d0,n_neighb);
p.Marker = 'o';
p.Color = 'k';
p.LineStyle = 'none';
ylim([6 Inf])
xlim([1e-2 1])
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'FontSize',20/1.5)

text(1.5e-2,1e1,['$n^{-1/3} \sim $',num2str(A,'%.2f'),...
    '$ {d^*}^{-1}$ $(R^2$ = ',num2str(R2,'%.2f'),'$)$'],...
    'FontSize',15,'Interpreter','latex')

yticks([1e2 1e4 1e6])

xlabel('$d/(Nb)$','FontSize',20,'Interpreter','latex')
ylabel('$n$','FontSize',20,'Interpreter','latex')
set(gcf,'color','w')
pbaspect([1 1 1])

set(gcf,'Position',[100 100 Size Size])


figure(2); clf; hold on

x = d0;
y = n_neighb.^(1/3);
err = zeros(size(n_neighb));
A_nom = 1.632; B_nom = -1; C_nom = 0;
[A,B,~,R2,~] = FitFixedPwr(A_nom,B_nom,C_nom,x,y,err,1);


plot(d0,A*x.^B,'r--','LineWidth',1.5);

p = plot(d0,y);
p.Marker = 'o';
p.Color = 'k';
p.LineStyle = 'none';
ylim([1 Inf])
xlim([1e-2 1])
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'FontSize',20/1.5)

text(1.5e-2,0.15e1,['$n^{1/3} \sim $',num2str(A,'%.2f'),...
    '$ {d^*}^{-1}$ $(R^2$ = ',num2str(R2,'%.2f'),'$)$'],...
    'FontSize',15,'Interpreter','latex')

yticks([1e0 1e1 1e2 1e6])

xlabel('$d/(Nb)$','FontSize',20,'Interpreter','latex')
ylabel('$n^{1/3}$','FontSize',20,'Interpreter','latex')
set(gcf,'color','w')
pbaspect([1 1 1])

set(gcf,'Position',[100 100 Size Size])


figure(3); clf; hold on

plot([1e-2 1],[0.75 0.75],'r--','LineWidth',1.5)
p = plot(d0,d_bar);
p.Marker = 'o';
p.Color = 'k';
p.LineStyle = 'none';
ylim([0 1])
% set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'FontSize',20/1.5)

text(1.5e-2,0.2,'$\bar d \sim 0.75$',...
    'FontSize',15,'Interpreter','latex')

xlabel('$d/(Nb)$','FontSize',20,'Interpreter','latex')
ylabel('$\bar d^*$ ($Nb$)','FontSize',20,'Interpreter','latex')
set(gcf,'color','w')
pbaspect([1 1 1])

set(gcf,'Position',[100 100 Size Size])

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B,C,R2,delta_A]...
    = FitFixedPwr(A_nom,B_nom,C_nom,x,y,err,find_conf_int)

global NoFitAttempts PlotDuringFitting plot_ct

NoPts = 21;
R2 = 0;
ct = 0;
wb5 = waitbar(0,'Fitting Fixed Power...');
while R2<0.995
    ct = ct+1;
    waitbar(ct/NoFitAttempts,wb5,'Fitting Fixed Power...')

    if ct==1
        sig_A = 1/50*A_nom;
        sig_B = 1/8*B_nom;
        sig_C = 1/4*C_nom;
    else
        sig_A = sig_A*0.95;
        sig_B = sig_B*0.95;
        sig_C = sig_C*0.95;
    end
    A_rng = [(linspace(A_nom-sig_A,A_nom+sig_A,NoPts))';A_nom];
    A_rng = unique(A_rng); 
    B_rng = B_nom;
    C_rng = [(linspace(C_nom-sig_C,C_nom+sig_C,NoPts))';C_nom];
    C_rng = unique(C_rng);

    Perms = zeros(length(A_rng)*length(B_rng)*length(C_rng),1);
    i = 0;
    for mi=1:length(A_rng)
        for bi=1:length(B_rng)
            for ci=1:length(C_rng)
                i = i+1;
                Perms(i,1) = A_rng(mi);
                Perms(i,2) = B_rng(bi);
                Perms(i,3) = C_rng(ci);
            end
        end
    end
    R2_all = zeros(size(Perms,1),1);
    A_all = zeros(size(Perms,1),1);
    B_all = zeros(size(Perms,1),1);
    C_all = zeros(size(Perms,1),1);
    Chi2 = zeros(size(Perms,1),1);

    for i=1:size(Perms,1)
        A_all(i) = Perms(i,1);
        B_all(i) = Perms(i,2);
        C_all(i) = Perms(i,3);

        ft = A_all(i)*x.^B_all(i) + C_all(i);

        RSS = sum((y-ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_all(i) = 1-RSS/TSS;
        Chi2(i) = sum(((y-ft)./err).^2);
    end
    diff = abs(Chi2);
    indx = find(diff==min(diff),1,'first');

    A_nom = A_all(indx);
    B_nom = B_all(indx);
    C_nom = C_all(indx);
    R2 = R2_all(indx);

    if PlotDuringFitting==1 && ~mod(ct,plot_ct)
        figure(100); clf; hold on
        scatter(x,y,'k','filled')
        plot(x,A_nom*x.^B_nom + C_nom,'k--')
    end
    if ct>NoFitAttempts
        break;
    end
end
A = A_nom;
B = B_nom;
C = C_nom;
if PlotDuringFitting==1
    figure(100); close
end
close(wb5)

if find_conf_int==1 %Perturb each parameter to find range in which chi2<=1
    R2_temp = R2;
    A_temp = A;
    dA = 0.01*A;
    while R2_temp>0.05
        A_temp = A_temp+dA;
        y_ft = A_temp*x.^B + C;
        RSS = sum((y-y_ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_temp = 1-RSS/TSS;
    end
    delta_A_pos = abs(A_temp-A);

    R2_temp = R2;
    A_temp = A;
    dA = 0.01*A;
    while R2_temp>0.05
        A_temp = A_temp-dA;
        y_ft = A_temp*x.^B + C;
        RSS = sum((y-y_ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_temp = 1-RSS/TSS;
    end
    delta_A_neg = abs(A_temp-A);

    delta_A = (delta_A_neg + delta_A_pos)/2;
end

end
