function CompareData(Package,Override,ToggleDynamics,...
    MakeHistos,MakeMoviesHistos,OverrideHistogramMovies,...
    OverrideMSDFitData,WRTStretchOrLength,HistoUnits,LC,DC,BSOM)

% Computes and plots stress-strain for every iteration of simulation

global LineWidth TurnOnDynamics FontSize...
    LengthConversion DamperConversion DataSize BeadSpringOrMeso

LineWidth = 1.5;
FontSize = 20;
DataSize = 30;
TurnOnDynamics = ToggleDynamics;
LengthConversion = LC;
DamperConversion = DC;
BeadSpringOrMeso = BSOM;

% Samples = unique(Package(:,1));
Np = unique(Package(:,2));        %Number of molecules
Ns = unique(Package(:,3));        %Number of stickeres per tether site
T = unique(Package(:,4));         %Temperature
Nb = unique(Package(:,5));        %Contour length of chains
ka_in = unique(Package(:,6));        %Activation energy of association
kd_in = unique(Package(:,7));        %Activation energy of dissociation
f0 = unique(Package(:,8));       %Force sensitivity to dissociation
dt = unique(Package(:,9));       %timestep size
% DiffCoeffs = unique(Package(:,11));
N_Kuhns = unique(Package(:,12));
b = unique(Package(:,13));
N = Np*Ns;

%% Sweepign parameters are N_Kuhn and b
damps = unique(Package(:,10));    %damping coefficient in units [mass/time]
Separations = unique(Package(:,14));

%% Define filenames for outputs
DefineMSDFileNames;

%% Sweep over all dampers. For each damper, sweep over N
[R2,tau0,taur,MSD_ss,R2_2,tau0_2,taur_2,alpha] = ...
    GeneratePlotsWithNSwept(OverrideMSDFitData,Override,Package,...
    damps,N_Kuhns,Np,Ns,T,Nb,ka_in,kd_in,f0,dt,b,N,Separations);

%% Sweep over all N. For each N, sweep over damper
GeneratePlotsWithDamperSwept(OverrideMSDFitData,...
    tau0,taur,tau0_2,taur_2,alpha,Override,Package,...
    damps,N_Kuhns,Np,Ns,T,Nb,ka_in,kd_in,f0,dt,b,N,Separations);

%% Sweep overa all N and plot fitted parameters wrt damper
PlotDiffusionModelParametersWithNSwept(N_Kuhns,damps,tau0,taur,R2,...
    tau0_2,taur_2,alpha,R2_2,OverrideMSDFitData,b)

%% Sweep overa all dampers and plot fitted parameters wrt N
PlotDiffusionModelParametersWithDSwept(N_Kuhns,damps,tau0,taur,R2,...
    tau0_2,taur_2,alpha,R2_2,OverrideMSDFitData,b)

%% Make surface plots of fitted parameters wrt dampers and N
MakeSurfacePlots(damps,N_Kuhns,b,MSD_ss,tau0,taur,R2,tau0_2,R2_2,taur_2,alpha)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeSurfacePlots(damps,N_Kuhns,b,MSD_ss,tau0,taur,R2,tau0_2,R2_2,taur_2,alpha)

global LengthConversion DamperConversion FontSize AddOn

% MSD_ss
figure(1)
clf; hold on
s = surf(damps*DamperConversion,N_Kuhns,MSD_ss'*LengthConversion^2);
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'FontSize',FontSize/1.5)

xlim([min(damps) max(damps)]*DamperConversion)
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$\langle \textbf{r}^2 (t) \rangle_{ss}$ (nm$^{2}$)');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',FontSize);

xlabel('$\gamma$ (N s m$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
ylabel('$N$','FontSize',FontSize,'Interpreter','latex')
title(['$b$ = ',num2str(b*LengthConversion*1e9,'%.2f'),' nm'],...
    'FontSize',FontSize/1.5','Interpreter','latex')
pbaspect([1 1 1])
set(gcf,'Color','w')

FileName = 'Output Plots/Steady State MSD';
saveas(gcf,[FileName,AddOn,'.png'])
saveas(gcf,[FileName,AddOn,'.fig'])

%tau0
figure(2)
clf; hold on
s = surf(damps*DamperConversion,N_Kuhns,tau0');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'ColorScale','log')
set(gca,'FontSize',FontSize/1.5)

xlim([min(damps) max(damps)]*DamperConversion)
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$\tau_0$ (s)');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',FontSize);

xlabel('$\gamma$ (N s m$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
ylabel('$N$','FontSize',FontSize,'Interpreter','latex')
title(['$b$ = ',num2str(b*LengthConversion*1e9,'%.2f'),' nm'],...
    'FontSize',FontSize/1.5','Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = 'Output Plots/tau0';
saveas(gcf,[FileName,AddOn,'.png'])
saveas(gcf,[FileName,AddOn,'.fig'])

%taur
figure(3)
clf; hold on
s = surf(damps*DamperConversion,N_Kuhns,taur');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'ColorScale','log')
set(gca,'FontSize',FontSize/1.5)

xlim([min(damps) max(damps)]*DamperConversion)
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$\tau_r$ (s)');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',FontSize);

xlabel('$\gamma$ (N s m$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
ylabel('$N$','FontSize',FontSize,'Interpreter','latex')
title(['$b$ = ',num2str(b*LengthConversion*1e9,'%.2f'),' nm'],...
    'FontSize',FontSize/1.5','Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = 'Output Plots/taur';
saveas(gcf,[FileName,AddOn,'.png'])
saveas(gcf,[FileName,AddOn,'.fig'])

%R2
figure(4)
clf; hold on
s = surf(damps*DamperConversion,N_Kuhns,R2');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'FontSize',FontSize/1.5)

xlim([min(damps) max(damps)]*DamperConversion)
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$R^2$');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',FontSize);
caxis([0 1])

xlabel('$\gamma$ (N s m$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
ylabel('$N$','FontSize',FontSize,'Interpreter','latex')
title(['$b$ = ',num2str(b*LengthConversion*1e9,'%.2f'),' nm'],...
    'FontSize',FontSize/1.5','Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = 'Output Plots/R2';
saveas(gcf,[FileName,AddOn,'.png'])
saveas(gcf,[FileName,AddOn,'.fig'])

% tau0 general subdiffusion
figure(5)
clf; hold on
s = surf(damps*DamperConversion,N_Kuhns,tau0_2');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'ColorScale','log')
set(gca,'FontSize',FontSize/1.5)

xlim([min(damps) max(damps)]*DamperConversion)
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$\tau_0^\alpha$ (s)');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',FontSize);

xlabel('$\gamma$ (N s m$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
ylabel('$N$','FontSize',FontSize,'Interpreter','latex')
title(['$b$ = ',num2str(b*LengthConversion*1e9,'%.2f'),' nm'],...
    'FontSize',FontSize/1.5','Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = 'Output Plots/tau0 General Sub-diffusion';
saveas(gcf,[FileName,AddOn,'.png'])
saveas(gcf,[FileName,AddOn,'.fig'])

% taur general subdiffusion
figure(6)
clf; hold on
s = surf(damps*DamperConversion,N_Kuhns,taur_2');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'ColorScale','log')
set(gca,'FontSize',FontSize/1.5)

xlim([min(damps) max(damps)]*DamperConversion)
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$\tau_r^\alpha$ (s)');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',FontSize);

xlabel('$\gamma$ (N s m$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
ylabel('$N$','FontSize',FontSize,'Interpreter','latex')
title(['$b$ = ',num2str(b*LengthConversion*1e9,'%.2f'),' nm'],...
    'FontSize',FontSize/1.5','Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = 'Output Plots/taur General Sub-diffusion';
saveas(gcf,[FileName,AddOn,'.png'])
saveas(gcf,[FileName,AddOn,'.fig'])

% alpha general subdiffusion
figure(7)
clf; hold on
s = surf(damps*DamperConversion,N_Kuhns,alpha');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'FontSize',FontSize/1.5)

xlim([min(damps) max(damps)]*DamperConversion)
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$\alpha$');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',FontSize);
caxis([0.25 0.5])

xlabel('$\gamma$ (N s m$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
ylabel('$N$','FontSize',FontSize,'Interpreter','latex')
title(['$b$ = ',num2str(b*LengthConversion*1e9,'%.2f'),' nm'],...
    'FontSize',FontSize/1.5','Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = 'Output Plots/alpha General Sub-diffusion';
saveas(gcf,[FileName,AddOn,'.png'])
saveas(gcf,[FileName,AddOn,'.fig'])

% R2 general subdiffusion
figure(8)
clf; hold on
s = surf(damps*DamperConversion,N_Kuhns,R2_2');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'FontSize',FontSize/1.5)

xlim([min(damps) max(damps)]*DamperConversion)
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$R^2_\alpha$');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',FontSize);
caxis([0 1])

xlabel('$\gamma$ (N s m$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
ylabel('$N$','FontSize',FontSize,'Interpreter','latex')
title(['$b$ = ',num2str(b*LengthConversion*1e9,'%.2f'),' nm'],...
    'FontSize',FontSize/1.5','Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = 'Output Plots/R2 General Sub-diffusion';
saveas(gcf,[FileName,AddOn,'.png'])
saveas(gcf,[FileName,AddOn,'.fig'])

close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DefineMSDFileNames

global MSD_ss_filename tau0_filename taur_filename R2_filename...
    tau0_SD_filename taur_SD_filename alpha_filename R2_SD_filename...
    BeadSpringOrMeso AddOn

if BeadSpringOrMeso==0
    AddOn = 'Bead.';
else
    AddOn = 'Meso.';
end
MSD_ss_filename = ['Data/Compiled Outputs/',AddOn,'Steady state msd.txt'];
tau0_filename = ['Data/Compiled Outputs/',AddOn,'Rouse tau0.txt'];
taur_filename = ['Data/Compiled Outputs/',AddOn,'Rouse taur.txt'];
R2_filename = ['Data/Compiled Outputs/',AddOn,'Rouse R2.txt'];
tau0_SD_filename = ['Data/Compiled Outputs/',AddOn,'General tau0.txt'];
taur_SD_filename = ['Data/Compiled Outputs/',AddOn,'General taur.txt'];
alpha_filename = ['Data/Compiled Outputs/',AddOn,'General alpha.txt'];
R2_SD_filename = ['Data/Compiled Outputs/',AddOn,'General R2.txt'];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotDiffusionModelParametersWithDSwept(N_Kuhns,damps,tau0,taur,R2,...
    tau0_2,taur_2,alpha,R2_2,OverrideMSDFitData,b)

global FontSize DataSize AddOn

%Sweep over dampers
for dm=1:length(damps)
    damp = damps(dm);
    tau0_temp = tau0(dm,:);
    taur_temp = taur(dm,:);
    R2_temp = R2(dm,:);

    tau02_temp = tau0_2(dm,:);
    taur2_temp = taur_2(dm,:);
    alpha_temp = alpha(dm,:);
    R22_temp = R2_2(dm,:);

    Folder = ['Output Plots/gamma ',num2str(damp,'%.2e')];
    if ~isfolder(Folder)
        mkdir(Folder)
    end
    FileName = [Folder,'/R2 vs N.fig'];

    tau0_eff = zeros(size(tau0_temp));
    tauR_eff = zeros(size(tau0_temp));
    for nk=1:length(N_Kuhns)
        N_Kuhn = N_Kuhns(nk);
        [tau0_eff(nk),tauR_eff(nk),~] = ComputeEffectiveDamper(N_Kuhn,b,damp);
    end
    
    if ~isfile(FileName) || OverrideMSDFitData==1
        figure(1); clf; hold on
        s1 = scatter(N_Kuhns,taur_temp,'k','filled');
        s1.SizeData = DataSize;
        s2 = scatter(N_Kuhns,taur2_temp,'k');
        s2.SizeData = DataSize;
        s3 = scatter(N_Kuhns,tauR_eff,'kx');
        s3.SizeData = DataSize;
        set(gcf,'Color','w')
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        pbaspect([1 1 1])
        xlabel('$N$','FontSize',FontSize,'Interpreter','latex')
        ylabel('$\tau_r$ (s)','FontSize',FontSize,'Interpreter','latex')
        title(['$\gamma =$ ',num2str(damp,'%.2e')],'FontSize',FontSize,'Interpreter','latex')

        l = legend([s1 s2 s3],...
            'Rouse Model','General Subdiff. Model','Set \textit{a priori}');
        l.Location = 'Northwest';
        l.Interpreter = 'latex';
        
        FileName = [Folder,'/',AddOn,'rouse time vs damper.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'rouse time vs damper.fig'];
        saveas(gcf,FileName)

        figure(2); clf; hold on
        s1 = scatter(N_Kuhns,tau0_temp,'k','filled');
        s1.SizeData = DataSize;
        s2 = scatter(N_Kuhns,tau02_temp,'k');
        s2.SizeData = DataSize;
        s3 = scatter(N_Kuhns,tauR_eff,'kx');
        s3.SizeData = DataSize;
        set(gcf,'Color','w')
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        pbaspect([1 1 1])
        xlabel('$N$','FontSize',FontSize,'Interpreter','latex')
        ylabel('$\tau_0$ (s)','FontSize',FontSize,'Interpreter','latex')
        title(['$\gamma =$ ',num2str(damp,'%.2e')],'FontSize',FontSize,'Interpreter','latex')

        x = N_Kuhns;
        y = tau0_temp';

        NoPts = 21;
        A_nom = -tau0_temp(1);
        B_nom = tau0_temp(1);
%         A_nom = 110*tau0_temp(1);
%         B_nom = 1;

        StatusString = 'Fiting tau0 vs N...';
        xlab_str = 'N'; ylab_str = '\tau_0';
        [A,B,R2_tau0vsN] = FitLine(x,y,A_nom,B_nom,NoPts,...
            StatusString,xlab_str,ylab_str);

        figure(2)
        plot(x,A*x+B,'k--')
        string = ['$\tau_0 \approx$',num2str(A,'%.2e'),'$N +$ ',num2str(B,'%.1f'),...
            ' $R^2$ = ',num2str(R2_tau0vsN,'%.3f')];
%         plot(x,A*x.^B,'k--')
%         string = ['$\tau_0 \approx$',num2str(A_nom,'%.2e'),'$N ^{',num2str(B_nom,'%.1f'),'}$',...
%             ' $R^2$ = ',num2str(R2_tau0vsN,'%.3f')];
        text(N_Kuhns(1)+1,1.2*tau0_temp(1),string,'FontSize',FontSize/2,'Interpreter','latex')

        FileName = [Folder,'/',AddOn,'diffusion time vs N.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'diffusion time vs N.fig'];
        saveas(gcf,FileName)

        figure(3); clf; hold on
        s = scatter(N_Kuhns,alpha_temp,'k');
        s.SizeData = DataSize;
        set(gcf,'Color','w')
        set(gca,'XScale','log')
        ylim([0 1])
        pbaspect([1 1 1])
        xlabel('$N$','FontSize',FontSize,'Interpreter','latex')
        ylabel('$\alpha$','FontSize',FontSize,'Interpreter','latex')
        title(['$\gamma =$ ',num2str(damp,'%.2e')],'FontSize',FontSize,'Interpreter','latex')

        x = N_Kuhns;
        y = alpha_temp';

        NoPts = 13;
        A_nom = 0.001;
        B_nom = 0.5;

        StatusString = 'Fiting alpha vs N...';
        xlab_str = 'N'; ylab_str = '\alpha';
        [A,B,R2_alphavsN] = FitLine(x,y,A_nom,B_nom,NoPts,...
            StatusString,xlab_str,ylab_str);

        figure(3)
        plot(x,A*x+B,'k--')
        string = ['$\tau_0 \approx$',num2str(A,'%.2e'),'$N +$ ',num2str(B,'%.1f'),...
            ' $R^2$ = ',num2str(R2_alphavsN,'%.3f')];
%         plot(x,A*x.^B,'k--')
%         string = ['$\alpha \approx$',num2str(A,'%.2e'),'$N ^{',num2str(B,'%.2f'),'}$',...
%             ' $R^2$ = ',num2str(R2_alphavsN,'%.3f')];
        text(N_Kuhns(1)+1,0.9,string,'FontSize',FontSize/2,'Interpreter','latex')

        FileName = [Folder,'/',AddOn,'difussion power vs N.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'difussion power vs N.fig'];
        saveas(gcf,FileName)

        figure(4); clf; hold on
        s = scatter(N_Kuhns,R2_temp,'k','filled');
        s.SizeData = DataSize;
        s = scatter(N_Kuhns,R22_temp,'k');
        s.SizeData = DataSize;
        plot(N_Kuhns,ones(size(N_Kuhns)),'k--')
        set(gcf,'Color','w')
        set(gca,'XScale','log')
        pbaspect([1 1 1])
        ylim([0 1])
        xlabel('$N$','FontSize',FontSize,'Interpreter','latex')
        ylabel('$R^2$','FontSize',FontSize,'Interpreter','latex')
        title(['$\gamma =$ ',num2str(damp,'%.3f')],'FontSize',FontSize,'Interpreter','latex')

        FileName = [Folder,'/',AddOn,'R2 vs N.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'R2 vs N.fig'];
        saveas(gcf,FileName)
    end
end
close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotDiffusionModelParametersWithNSwept(N_Kuhns,damps,tau0,taur,R2,...
    tau0_2,taur_2,alpha,R2_2,OverrideMSDFitData,b)

global FontSize DamperConversion DataSize AddOn

for nk=1:length(N_Kuhns)
    N_Kuhn = N_Kuhns(nk);
    tau0_temp = tau0(:,nk);
    taur_temp = taur(:,nk);
    R2_temp = R2(:,nk);

    tau02_temp = tau0_2(:,nk);
    taur2_temp = taur_2(:,nk);
    alpha_temp = alpha(:,nk);
    R22_temp = R2_2(:,nk);

    Folder = ['Output Plots/N_Kuhn ',num2str(N_Kuhn)];
    if ~isfolder(Folder)
        mkdir(Folder)
    end
    FileName = [Folder,'/R2 vs damper.fig'];

    tau0_eff = zeros(size(tau0_temp));
    tauR_eff = zeros(size(tau0_temp));
    for d=1:length(damps)
        damp = damps(d);
        [tau0_eff(d),tauR_eff(d),R_N] = ComputeEffectiveDamper(N_Kuhn,b,damp);
    end

    if ~isfile(FileName) || OverrideMSDFitData==1
        figure(1); clf; hold on
        s1 = scatter(damps*DamperConversion,taur_temp,'k','filled');
        s1.SizeData = DataSize;
        s2 = scatter(damps*DamperConversion,taur2_temp,'k');
        s2.SizeData = DataSize;
        s3 = scatter(damps*DamperConversion,tauR_eff,'kx');
        s3.SizeData = DataSize;
        set(gcf,'Color','w')
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        pbaspect([1 1 1])
        xlabel('$\gamma$ (N s m$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
        ylabel('$\tau_r$ (s)','FontSize',FontSize,'Interpreter','latex')
        title(['$N =$ ',num2str(N_Kuhn)],'FontSize',FontSize,'Interpreter','latex')

        l = legend([s1 s2 s3],...
            'Fitted Rouse Model','Fitted General Subdiffusion Model','Set \textit{a priori}');
        l.Location = 'Southeast';
        l.Interpreter = 'latex';

        FileName = [Folder,'/',AddOn,'rouse time vs damper.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'rouse time vs damper.fig'];
        saveas(gcf,FileName)


        % Plot tau0 wrt. gamma
        figure(2); clf; hold on
        s1 = scatter(damps*DamperConversion,tau0_temp,'k','filled');
        s1.SizeData = DataSize;
        s2 = scatter(damps*DamperConversion,tau02_temp,'k');
        s2.SizeData = DataSize;
        s3 = scatter(damps*DamperConversion,tau0_eff,'kx');
        s3.SizeData = DataSize;
        set(gcf,'Color','w')
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        pbaspect([1 1 1])
        xlabel('$\gamma$ (N s m$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
        ylabel('$\tau_0$ (s)','FontSize',FontSize,'Interpreter','latex')
        title(['$N =$ ',num2str(N_Kuhn)],'FontSize',FontSize,'Interpreter','latex')

        x = damps*DamperConversion;
        y = tau0_temp;

        NoPts = 21;
        A_nom = 15;
        B_nom = 0;

        StatusString = 'Fiting tau0 vs damper...';
        xlab_str = '\gamma'; ylab_str = '\tau_0';
        [A,B,R2_tau0vsgamma] = FitLine(x,y,A_nom,B_nom,NoPts,...
            StatusString,xlab_str,ylab_str);
       
        plot(x,A*x + B,'k--')
        set(gca,'XScale','log')
        set(gca,'YScale','log')

        figure(2)
        string = ['$\tau_0 \approx$ (',num2str(A,'%.2e'),')$\gamma +$ ',num2str(B,'%.1f'),...
            ' $R^2$ = ',num2str(R2_tau0vsgamma,'%.3f')];
        text(damps(2)*DamperConversion,max(tau0_eff),string,...
            'FontSize',FontSize/2,'Interpreter','latex')


        l = legend([s1 s2 s3],...
            'Rouse Model','General Subdiffusion Model','\textit{A priori}');
        l.Location = 'Southeast';
        l.Interpreter = 'latex';
        
        FileName = [Folder,'/',AddOn,'diffusion time vs damper.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'diffusion time vs damper.fig'];
        saveas(gcf,FileName)

        figure(3); clf; hold on
        s = scatter(damps*DamperConversion,alpha_temp,'k');
        s.SizeData = DataSize;
        set(gcf,'Color','w')
        set(gca,'XScale','log')
        ylim([0 1])
        pbaspect([1 1 1])
        xlabel('$\gamma$ (N s m$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
        ylabel('$\alpha$','FontSize',FontSize,'Interpreter','latex')
        title(['$N =$ ',num2str(N_Kuhn)],'FontSize',FontSize,'Interpreter','latex')

        FileName = [Folder,'/',AddOn,'difussion power vs damper.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'difussion power vs damper.fig'];
        saveas(gcf,FileName)

        figure(4); clf; hold on
        s = scatter(damps*DamperConversion,R2_temp,'k','filled');
        s.SizeData = DataSize;
        s = scatter(damps*DamperConversion,R22_temp,'k');
        s.SizeData = DataSize;
        plot(damps*DamperConversion,ones(size(damps)),'k--')
        set(gcf,'Color','w')
        set(gca,'XScale','log')
        pbaspect([1 1 1])
        ylim([0 1])
        xlabel('$\gamma$ (N s m$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
        ylabel('$R^2$','FontSize',FontSize,'Interpreter','latex')
        title(['$N =$ ',num2str(N_Kuhn)],'FontSize',FontSize,'Interpreter','latex')

        FileName = [Folder,'/',AddOn,'R2 vs damper.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'R2 vs damper.fig'];
        saveas(gcf,FileName)
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B,R2] = FitLine(x,y,A_nom,B_nom,NoPts,...
    StatusString,xlab_str,ylab_str)

R2 = 0;
ct = 0;
wb5 = waitbar(0,StatusString);
figure(100)
while R2<0.999
    ct = ct+1;
    waitbar(R2/1,wb5,StatusString)

    if ct==1
        sigA = 1/2*A_nom;
        sigB = 1/2*B_nom;
    else
        sigA = sigA*0.95;
        sigB = sigB*0.95;
    end
    A_rng = [(linspace(A_nom-sigA,A_nom+sigA,NoPts))';A_nom];
    B_rng = [(linspace(B_nom-sigB,B_nom+sigB,NoPts))';B_nom];
    A_rng = unique(A_rng); %A_rng(A_rng<=0) = [];
    B_rng = unique(B_rng); 

    Perms = zeros(length(A_rng)*length(B_rng),2);
    i = 0;
    for ai=1:length(A_rng)
        for bi=1:length(B_rng)
            i = i+1;
            Perms(i,1) = A_rng(ai);
            Perms(i,2) = B_rng(bi);
        end
    end
    R2_all = zeros(size(Perms,1),1);
    A_all = zeros(size(Perms,1),1);
    B_all = zeros(size(Perms,1),1);

    for i=1:size(Perms,1)
        A_all(i) = Perms(i,1);
        B_all(i) = Perms(i,2);
        ft = A_all(i)*x + B_all(i);

        RSS = sum((y-ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_all(i) = 1-RSS/TSS;
    end
    diff = abs(R2_all-1);
    indx = find(diff==min(diff),1,'first');

    A_nom = A_all(indx);
    B_nom = B_all(indx);
    R2 = R2_all(indx);

    clf; hold on
    scatter(x,y,'k','filled')
    plot(x,A_nom*x + B_nom,'k--')
    xlabel(xlab_str)
    ylabel(ylab_str)

    if ct>75
        break;
    end
end
A = A_nom;
B = B_nom;
figure(100); close
close(wb5)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B,R2] = FitPowerLaw(x,y,A_nom,B_nom,NoPts,...
    StatusString,xlab_str,ylab_str)

R2 = 0;
ct = 0;
wb5 = waitbar(0,StatusString);
figure(100)
while R2<0.999
    ct = ct+1;
    waitbar(R2/1,wb5,StatusString)

    if ct==1
        sigA = 1/2*A_nom;
        sigB = 1/2*B_nom;
    else
        sigA = sigA*0.95;
        sigB = sigB*0.95;
    end
    A_rng = [(linspace(A_nom-sigA,A_nom+sigA,NoPts))';A_nom];
    B_rng = [(linspace(B_nom-sigB,B_nom+sigB,NoPts))';B_nom];
    A_rng = unique(A_rng); %A_rng(A_rng<=0) = [];
    B_rng = unique(B_rng); %B_rng(B_rng<=0) = [];

    Perms = zeros(length(A_rng)*length(B_rng),2);
    i = 0;
    for ai=1:length(A_rng)
        for bi=1:length(B_rng)
            i = i+1;
            Perms(i,1) = A_rng(ai);
            Perms(i,2) = B_rng(bi);
        end
    end
    R2_all = zeros(size(Perms,1),1);
    A_all = zeros(size(Perms,1),1);
    B_all = zeros(size(Perms,1),1);

    for i=1:size(Perms,1)
        A_all(i) = Perms(i,1);
        B_all(i) = Perms(i,2);
        ft = A_all(i)*x.^(B_all(i));

        RSS = sum((y-ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_all(i) = 1-RSS/TSS;
    end
    diff = abs(R2_all-1);
    indx = find(diff==min(diff),1,'first');

    A_nom = A_all(indx);
    B_nom = B_all(indx);
    R2 = R2_all(indx);

    clf; hold on
    scatter(x,y,'k','filled')
    plot(x,A_nom*x.^B_nom,'k--')
    xlabel(xlab_str)
    ylabel(ylab_str)
%     set(gca,'XScale','log')
%     set(gca,'YScale','log')
    if ct>50
        break;
    end
end
A = A_nom;
B = B_nom;
figure(100); close
close(wb5)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GeneratePlotsWithDamperSwept(OverrideMSDFitData,...
    tau0,taur,tau0_2,taur_2,alpha,Override,Package,...
    damps,N_Kuhns,Np,Ns,T,Nb,ka_in,kd_in,f0,dts,b,N,Separations)

global ConstantsFileName FontSize BeadSpringOrMeso AddOn

TotalSims = size(Package,1);

wb1 = waitbar(0,'Post processing for const. N, and swept dampers...');
% MSD wrt time for each damper while sweeping N
SimCt = 0;
for nk = 1:length(N_Kuhns)
    N_Kuhn = N_Kuhns(nk);

    % Figures
    figure(1); clf; hold on % msd
    figure(2); clf; hold on % msd in tau0<t<taur
    figure(3); clf; hold on % Metric tensors of chains
    figure(4); clf; hold on % msd in tau0<t<taur v2

    ColorRange = (linspace(0,1,length(damps)))';
    Colors = [ColorRange flipud(ColorRange) flipud(ColorRange)];
    p0 = [];


    Folder = ['Output Plots/N_Kuhn ',num2str(N_Kuhn)];
    if ~isfolder(Folder)
        mkdir(Folder)
    end
    FileName = [Folder,'/MSD vs t.png'];

    if ~isfile(FileName) || OverrideMSDFitData==1
        for dm = 1:length(damps)
            damp = damps(dm);
            dt = dts(dm);

            Color = Colors(dm,:);
            LegendEntries{dm} = ['$\gamma$ = ',num2str(damp,'%.2e')];

            tau0_temp = tau0(dm,nk);
            taur_temp = taur(dm,nk);
            taur2_temp = taur_2(dm,nk);
            tau02_temp = tau0_2(dm,nk);
            alpha_temp = alpha(dm,nk);

            for sp = 1:length(Separations)
                Separation = Separations(sp);

                PackageTemp = Package(ismember(Package(:,2),Np),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,3),Ns),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,4),T),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,5),Nb),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,6),ka_in),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,7),kd_in),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,8),f0),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,9),dt),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,10),damp),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,12),N_Kuhn),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,13),b),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,14),Separation),:);

                [~,~,R_N] = ComputeEffectiveDamper(N_Kuhn,b,damp);

                if ~isempty(PackageTemp)
                    SimCt = SimCt+1;

                    waitbar(SimCt*size(PackageTemp,1)/TotalSims,wb1,...
                        'Post processing for const. N, and swept dampers...')

                    % Determine allocation size
                    EdgeFactor = 1;
                    [Corners,N_steps] = ...
                        ImportDomainBoundaries(EdgeFactor,Np,N,...
                        ka_in,kd_in,f0,dt,damp,N_Kuhn,b,Separation);
                    N_samples = size(PackageTemp,1);

                    %% Callout Input Script
                    % Initializes non-sweeping parameters
                    InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dt,damp,...
                        N_Kuhn,b,Separation);
                    SetDirAndFileNames;
                    DefineCompiledFileNames;
                    Table = readtable(ConstantsFileName);
                    constants = table2array(Table);
                    N_Kuhn = constants(3); b = constants(4);

                    %% Assemble or read all of the ensemble data
                    if BeadSpringOrMeso==0
                        [time,~,~,~,~,~,msd_st_mean,msd_st_err,...
                            g11_st_mean,g22_st_mean,g33_st_mean,...
                            g11_st_err,g22_st_err,g33_st_err]...
                            = AssembleTheDataBeadSpring(Corners,EdgeFactor,Np,N,...
                            ka_in,kd_in,f0,dt,damp,N_Kuhn,b,Separation,...
                            PackageTemp,Override,N_steps,N_samples);
                    else
                        [time,~,~,~,~,~,msd_st_mean,msd_st_err,...
                            g11_st_mean,g22_st_mean,g33_st_mean,...
                            g11_st_err,g22_st_err,g33_st_err]...
                            = AssembleTheDataMeso(Corners,EdgeFactor,Np,N,...
                            ka_in,kd_in,f0,dt,damp,N_Kuhn,b,Separation,...
                            PackageTemp,Override,N_steps,N_samples);
                    end

                    figure(1) % msd
                    Style = '-';
                    [p0(dm),~] = PlotCurve(time,msd_st_mean,msd_st_err,Color,Style);

                    stopindx = find(time>=taur_temp,1,'first');
                    if isempty(stopindx)
                        stopindx = length(time);
                    end
                    rng = (1:stopindx)';

                    figure(2)
                    Style = '-';
                    [~,~] = PlotCurve(time(rng)/taur_temp,msd_st_mean(rng),...
                        msd_st_err(rng),Color,Style);
                    x = (linspace(0,taur_temp,100))';
                    %             R_N = b*sqrt(N_Kuhn/12);
                    y = (R_N^2)*(x/tau0_temp).^0.5;
                    %                 y = (b^2)*(x/tau0_temp).^0.5;
                    p = plot(x/taur_temp,y);
                    p.Color = Color;
                    p.LineStyle = '--';
                    p.LineWidth = 1.5;

                    figure(4)

                    % For Subdiffusion
                    stopindx = find(time>=taur_temp,1,'first');
                    if isempty(stopindx)
                        stopindx = length(time);
                    end
                    rng = (1:stopindx)';

                    Style = '-';
                    [~,~] = PlotCurve(time(rng)/taur_temp,msd_st_mean(rng),...
                        msd_st_err(rng),Color,Style);
                    x = (linspace(0,taur_temp,100))';
                    %             R_N = b*sqrt(N_Kuhn/12);
                    y = (R_N^2)*(x/tau02_temp).^alpha_temp;
                    %             y = (b^2)*(x/tau02_temp).^alpha_temp;
                    p = plot(x/taur_temp,y);
                    p.Color = Color;
                    p.LineStyle = '--';
                    p.LineWidth = 1.5;


                    figure(3) % g11_st, g22_st, g33_st
                    Style = '-';
                    [~,~] = PlotCurve(time(:,1),g11_st_mean,g11_st_err,Color,Style);
                    Style = '--';
                    [~,~] = PlotCurve(time(:,1),g22_st_mean,g22_st_err,Color,Style);
                    Style = ':';
                    [~,~] = PlotCurve(time(:,1),g33_st_mean,g33_st_err,Color,Style);
                end
            end
        end
        figure(1)
        set(gca,'FontSize',FontSize/1.5)
        set(gcf,'Color','w')
        pbaspect([1 1 1])
        xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
        ylabel('$MSD$ ($\ell^2$)','FontSize',FontSize,'Interpreter','latex')
        title(['$N =$ ',num2str(N_Kuhn)],'FontSize',FontSize,'Interpreter','latex')

        l = legend(p0,LegendEntries);
        l.FontSize = FontSize/2;
        l.Interpreter = 'latex';
        l.Location = 'Southeast';

        FileName = [Folder,'/',AddOn,'MSD vs t.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'MSD vs t.fig'];
        saveas(gcf,FileName)

        figure(2)
        set(gca,'FontSize',FontSize/1.5)
        set(gcf,'Color','w')
        pbaspect([1 1 1])
        xlabel('$t/\tau_r$','FontSize',FontSize,'Interpreter','latex')
        ylabel('$MSD$ ($\ell^2$)','FontSize',FontSize,'Interpreter','latex')
        xlim([0 1])
        title(['$N =$ ',num2str(N_Kuhn)],'FontSize',FontSize,'Interpreter','latex')

        FileName = [Folder,'/',AddOn,'MSD Rouse inset vs t.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'MSD Rouse inset vs t.fig'];
        saveas(gcf,FileName)

        figure(3)
        set(gca,'FontSize',FontSize/1.5)
        set(gcf,'Color','w')
        pbaspect([1 1 1])
        xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
        ylabel('$g_{11}^s,g_{22}^s, g_{33}^s$','FontSize',FontSize,'Interpreter','latex')
        title(['$N =$ ',num2str(N_Kuhn)],'FontSize',FontSize,'Interpreter','latex')

        FileName = [Folder,'/',AddOn,'g-principle inset vs t.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'g-principle inset vs t.fig'];
        saveas(gcf,FileName)

        figure(4)
        set(gca,'FontSize',FontSize/1.5)
        set(gcf,'Color','w')
        pbaspect([1 1 1])
        xlabel('$t/\tau_r$','FontSize',FontSize,'Interpreter','latex')
        ylabel('$MSD$ ($\ell^2$)','FontSize',FontSize,'Interpreter','latex')
        title(['$N =$ ',num2str(N_Kuhn)],'FontSize',FontSize,'Interpreter','latex')
        xlim([0 1])

        FileName = [Folder,'/',AddOn,'MSD Subiffusive inset vs t.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'MSD Subiffusive inset vs t.fig'];
        saveas(gcf,FileName)
    end
end
close(wb1)
close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotAllHistograms(Override,MakeMoviesHistos,OverrideHistogramMovies,...
    WRTStretchOrLength,HistoUnits,Package,...
    damps,N_Kuhns,Np,Ns,T,Nb,ka_in,kd_in,f0,dts,b,N,Separations)

global BeadSpringOrMeso ConstantsFileName FileTagCompiled xlab...
    MovieName1 MovieName2 MovieName3 MovieName4 MovieName5 MovieName6...
    MovieName7 MovieName8 MovieName9 MovieName10 MovieName11 MovieFolder

if WRTStretchOrLength==0
    AddOn = '_stretch';
    xlab = '$\lambda$';
else
    AddOn = '_length';
    xlab = '$r$';
    if HistoUnits==1
        AddOn = '_length_nm';
        xlab = '$r$ (nm)';
    end
end

TotalSims = size(Package,1);
SimCt = 0;
OneDOnly=1;

wb1 = waitbar(0,'Plotting histograms...');

% MSD wrt time for each damper while sweeping N
for dm = 1:length(damps)
    damp = damps(dm);
    dt = dts(dm);

    for nk = 1:length(N_Kuhns)
        N_Kuhn = N_Kuhns(nk);

        for sp=1:length(Separations)
            Separation = Separations(sp);

            PackageTemp = Package(ismember(Package(:,2),Np),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,3),Ns),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,4),T),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,5),Nb),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,6),ka_in),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,7),kd_in),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,8),f0),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,9),dt),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,10),damp),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,12),N_Kuhn),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,13),b),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,14),Separation),:);

            waitbar(SimCt*size(PackageTemp,1)/TotalSims,wb1,...
                'Plotting histograms...')
            if ~isempty(PackageTemp)
                SimCt = SimCt+1;

                %% Callout Input Script
                % Initializes non-sweeping parameters
                EdgeFactor = 1;
                InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dt,damp,...
                    N_Kuhn,b,Separation);
                SetDirAndFileNames;
                DefineCompiledFileNames;
                Table = readtable(ConstantsFileName);
                constants = table2array(Table);
                N_Kuhn = constants(3); b = constants(4);

                %Define movie names
                MovieFolder = ['Movies/damp ',num2str(damp),'.N ',num2str(N_Kuhn),...
                    '.b ',num2str(b,'%.2f')];
                MovieName1 = [MovieFolder,'/End-to-end/Norms',AddOn,FileTagCompiled];
                if OneDOnly~=1
                    MovieName2 = [MovieFolder,'/End-to-end/Norms_xyz',AddOn,FileTagCompiled];

                    MovieName3 = [MovieFolder,'/End-to-end/xy',AddOn,FileTagCompiled];
                    MovieName4 = [MovieFolder,'/End-to-end/yz',AddOn,FileTagCompiled];
                    MovieName5 = [MovieFolder,'/End-to-end/zx',AddOn,FileTagCompiled];

                    if CheckBackBoneVsStickers==1
                        MovieName6 = [MovieFolder,'/End-to-end/Norms_branch',AddOn,FileTagCompiled];
                        MovieName7 = [MovieFolder,'/End-to-end/xy_branch',AddOn,FileTagCompiled];
                        MovieName8 = [MovieFolder,'/End-to-end/yz_branch',AddOn,FileTagCompiled];

                        MovieName9 = [MovieFolder,'/End-to-end/Norms_backbone',AddOn,FileTagCompiled];
                        MovieName10 = [MovieFolder,'/End-to-end/xy_backbone',AddOn,FileTagCompiled];
                        MovieName11 = [MovieFolder,'/End-to-end/yz_backbone',AddOn,FileTagCompiled];
                    end
                end

                if MakeMoviesHistos==1 && (~isfile([MovieName1,'.avi']) || OverrideHistogramMovies==1)
                    % Determine allocation size
                    [Corners,N_steps] = ...
                        ImportDomainBoundaries(EdgeFactor,Np,N,ka_in,kd_in,f0,...
                        dt,damp,N_Kuhn,b,Separation);
                    N_samples = size(PackageTemp,1);

                    %% Assemble or read all of the ensemble data
                    if BeadSpringOrMeso==0
                        [time,stretch,rx,ry,rz,bondtypes,~,~,~,~,~,~,~,~]...
                            = AssembleTheDataBeadSpring(Corners,EdgeFactor,Np,N,...
                            ka_in,kd_in,f0,dt,damp,N_Kuhn,b,Separation,...
                            PackageTemp,Override,N_steps,N_samples);
                    else
                        [time,stretch,rx,ry,rz,bondtypes,~,~,~,~,~,~,~,~]...
                            = AssembleTheDataMeso(Corners,EdgeFactor,Np,N,...
                            ka_in,kd_in,f0,dt,damp,N_Kuhn,b,Separation,...
                            PackageTemp,Override,N_steps,N_samples);
                    end

                    %% Bond statistics
                    % Make histograms of end-to-end distributions (by bond type though)
                    PlotHistograms(rx,ry,rz,bondtypes,time,stretch,...
                        WRTStretchOrLength,HistoUnits,OneDOnly);
                end
            end
        end
    end
end
close(wb1)
close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R2,tau0,tauR,MSD_ss,R2_2,tauR_2,tau0_2,alpha] = ...
    GeneratePlotsWithNSwept(OverrideMSDFitData,Override,Package,...
    damps,N_Kuhns,Np,Ns,T,Nb,ka_in,kd_in,f0,dts,b,N,Separations)

global MSD_ss_filename tau0_filename taur_filename R2_filename...
    tau0_SD_filename taur_SD_filename alpha_filename R2_SD_filename...
    ConstantsFileName FontSize DamperConversion StartMSDMeasurePct...
    BeadSpringOrMeso AddOn

TotalSims = size(Package,1);

tau0 = zeros(length(damps),length(N_Kuhns));
R2 = zeros(length(damps),length(N_Kuhns));
tauR = zeros(length(damps),length(N_Kuhns));
tauR_2 = zeros(length(damps),length(N_Kuhns));
tau0_2 = zeros(length(damps),length(N_Kuhns));
alpha = zeros(length(damps),length(N_Kuhns));
R2_2 = zeros(length(damps),length(N_Kuhns));
MSD_ss = zeros(length(damps),length(N_Kuhns));
SimCt = 0;

wb1 = waitbar(0,'Post processing for const. damper, and swept N...');

if ~isfile(MSD_ss_filename) || OverrideMSDFitData==1
    % MSD wrt time for each damper while sweeping N
    for dm = 1:length(damps)
        damp = damps(dm);
        dt = dts(dm);

        Folder = ['Output Plots/gamma ',num2str(damp,'%.3e')];
        if ~isfolder(Folder)
            mkdir(Folder)
        end

        % Figures
        figure(1); clf; hold on % msd
        figure(2); clf; hold on % msd in tau0<t<taur
        figure(3); clf; hold on % Metric tensors of chains
        figure(4); clf; hold on % msd in tau0<t<taur v2

        ColorRange = (linspace(0,1,length(N_Kuhns)))';
        Colors = [ColorRange flipud(ColorRange) flipud(ColorRange)];
        LegendEntries = [];
        p0 = [];
        for nk = 1:length(N_Kuhns)
            N_Kuhn = N_Kuhns(nk);
            Color = Colors(nk,:);

            LegendEntries{nk} = ['$N =$ ',num2str(N_Kuhn)];

            [tau0_eff,tauR_eff,R_N] = ComputeEffectiveDamper(N_Kuhn,b,damp);

            for sp = 1:length(Separations)
                Separation = Separations(sp);

                PackageTemp = Package(ismember(Package(:,2),Np),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,3),Ns),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,4),T),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,5),Nb),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,6),ka_in),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,7),kd_in),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,8),f0),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,9),dt),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,10),damp),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,12),N_Kuhn),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,13),b),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,14),Separation),:);

                waitbar(SimCt*size(PackageTemp,1)/TotalSims,wb1,...
                    'Post processing for const. damper, and swept N...')

                if ~isempty(PackageTemp)
                    SimCt = SimCt+1;

                    % Determine allocation size
                    EdgeFactor = 1;
                    [Corners,N_steps] = ...
                        ImportDomainBoundaries(EdgeFactor,Np,N,ka_in,kd_in,f0,...
                        dt,damp,N_Kuhn,b,Separation);
                    N_samples = size(PackageTemp,1);

                    %% Callout Input Script
                    % Initializes non-sweeping parameters
                    InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dt,damp,...
                        N_Kuhn,b,Separation);
                    SetDirAndFileNames;
                    DefineCompiledFileNames;
                    Table = readtable(ConstantsFileName);
                    constants = table2array(Table);
                    N_Kuhn = constants(3); b = constants(4);

                    %% Assemble or read all of the ensemble data
                    if BeadSpringOrMeso==0
                        [time,~,~,~,~,~,msd_st_mean,msd_st_err,...
                            g11_st_mean,g22_st_mean,g33_st_mean,...
                            g11_st_err,g22_st_err,g33_st_err]...
                            = AssembleTheDataBeadSpring(Corners,EdgeFactor,Np,N,...
                            ka_in,kd_in,f0,dt,damp,N_Kuhn,b,Separation,...
                            PackageTemp,Override,N_steps,N_samples);
                    else
                        [time,~,~,~,~,~,msd_st_mean,msd_st_err,...
                            g11_st_mean,g22_st_mean,g33_st_mean,...
                            g11_st_err,g22_st_err,g33_st_err]...
                            = AssembleTheDataMeso(Corners,EdgeFactor,Np,N,...
                            ka_in,kd_in,f0,dt,damp,N_Kuhn,b,Separation,...
                            PackageTemp,Override,N_steps,N_samples);
                    end


                    %% Crop MSD Measurement
                    MSDStartIndx = ceil(StartMSDMeasurePct/100*length(time));
                    time = time(MSDStartIndx:end);
                    msd_st_mean = msd_st_mean(MSDStartIndx:end);
                    msd_st_err = msd_st_err(MSDStartIndx:end);
                    g11_st_mean = g11_st_mean(MSDStartIndx:end);
                    g22_st_mean = g22_st_mean(MSDStartIndx:end);
                    g33_st_mean = g33_st_mean(MSDStartIndx:end);
                    g11_st_err = g11_st_err(MSDStartIndx:end);
                    g22_st_err = g22_st_err(MSDStartIndx:end);
                    g33_st_err = g33_st_err(MSDStartIndx:end);
                    time = time-time(1);

                    figure(1) % msd
                    Style = '-';
                    [p0(nk),~] = PlotCurve(time,msd_st_mean,msd_st_err,Color,Style);
                    plot(time,(b^2)*ones(size(time)),'k--')

                    % Fit Rouse diffusion model in time domain tau_0<t<tau_r
                    FitMSD = 1;
                    if FitMSD==1
                        [R2(dm,nk),tau0(dm,nk),tauR(dm,nk),MSD_ss(dm,nk),rng] =...
                            FitRouseModel(time,msd_st_mean,tau0_eff,tauR_eff,R_N,N_Kuhn);

                        %% Fit the Rouse Model
                        figure(2)
                        Style = '-';
                        [~,~] = PlotCurve(time(rng)/tauR(dm,nk),msd_st_mean(rng),...
                            msd_st_err(rng),Color,Style);

                        x = (linspace(0,tauR(dm,nk),100))';
                        y = (R_N^2)*(x/tau0(dm,nk)).^0.5;
                        p = plot(x/tauR(dm,nk),y);
                        p.Color = Color;
                        p.LineStyle = '--';
                        p.LineWidth = 1.5;

%                         [~,~] = PlotCurve(time(rng)/tauR_eff,msd_st_mean(rng),...
%                             msd_st_err(rng),Color,Style);
% 
%                         x = time(rng);
%                         y = (sqrt(N_Kuhn/3)*b)^2*(x/tau0_eff).^0.5;
%                         plot(x/tauR_eff,y,'k')
% 
%                         x = (linspace(0,tauR(dm,nk),100))';
%                         y = (R_N^2)*(x/tau0(dm,nk)).^0.5;
%                         p = plot(x/tauR_eff,y);
%                         p.Color = Color;
%                         p.LineStyle = '--';
%                         p.LineWidth = 1.5;

                        %% Fit the General Sub-diffusion Model
                        figure(4)

                        [R2_2(dm,nk),tauR_2(dm,nk),tau0_2(dm,nk),alpha(dm,nk)] =...
                            FitGeneralSubDiffusionModel(time,msd_st_mean,rng,...
                            N_Kuhn,tau0(dm,nk),R_N);

                        Style = '-';
                        [~,~] = PlotCurve(time(rng)/tauR(dm,nk),msd_st_mean(rng),...
                            msd_st_err(rng),Color,Style);

                        x = (linspace(0,tauR(dm,nk),100))';
                        y = (R_N^2)*(x/tauR_2(dm,nk)).^alpha(dm,nk);
                        p = plot(x/tauR(dm,nk),y);
                        p.Color = Color;
                        p.LineStyle = '--';
                        p.LineWidth = 1.5;
                    end

                    figure(3) % g11_st, g22_st, g33_st
                    Style = '-';
                    [~,~] = PlotCurve(time(:,1),g11_st_mean,g11_st_err,Color,Style);
                    Style = '--';
                    [~,~] = PlotCurve(time(:,1),g22_st_mean,g22_st_err,Color,Style);
                    Style = ':';
                    [~,~] = PlotCurve(time(:,1),g33_st_mean,g33_st_err,Color,Style);
                end
            end
        end
        figure(1)
        set(gca,'FontSize',FontSize/1.5)
        set(gcf,'Color','w')
        pbaspect([1 1 1])
        xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
        ylabel('$MSD$ ($\ell^2$)','FontSize',FontSize,'Interpreter','latex')
        title(['$\gamma$ = ',num2str(damp,'%.2e'),' [force $\cdot$ time/length], ',...
            num2str(damp*DamperConversion,'%.3e'),' [N s m$^{-1}$]'],...
            'FontSize',FontSize/2,'Interpreter','latex')

        l = legend(p0,LegendEntries);
        l.FontSize = FontSize/2;
        l.Interpreter = 'latex';
        l.Location = 'Southeast';

        FileName = [Folder,'/',AddOn,'MSD vs t.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'MSD vs t.fig'];
        saveas(gcf,FileName)

        figure(2)
        set(gca,'FontSize',FontSize/1.5)
        set(gcf,'Color','w')
        pbaspect([1 1 1])
        xlabel('$t/\tau_r$','FontSize',FontSize,'Interpreter','latex')
        ylabel('$MSD$ ($\ell^2$)','FontSize',FontSize,'Interpreter','latex')
        xlim([0 1])
        title(['$\gamma$ = ',num2str(damp,'%.2e'),' [force $\cdot$ time/length], ',...
            num2str(damp*DamperConversion,'%.3e'),' [N s m$^{-1}$]'],...
            'FontSize',FontSize/2,'Interpreter','latex')

        FileName = [Folder,'/',AddOn,'MSD Rouse inset vs t.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'MSD Rouse inset vs t.fig'];
        saveas(gcf,FileName)

        figure(3)
        set(gca,'FontSize',FontSize/1.5)
        set(gcf,'Color','w')
        pbaspect([1 1 1])
        xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
        ylabel('$g_{11}^s,g_{22}^s, g_{33}^s$','FontSize',FontSize,'Interpreter','latex')
        title(['$\gamma$ = ',num2str(damp,'%.2e'),' [force $\cdot$ time/length], ',...
            num2str(damp*DamperConversion,'%.3e'),' [N s m$^{-1}$]'],...
            'FontSize',FontSize/2,'Interpreter','latex')

        FileName = [Folder,'/',AddOn,'g-principle inset vs t.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'g-principle inset vs t.fig'];
        saveas(gcf,FileName)

        figure(4)
        set(gca,'FontSize',FontSize/1.5)
        set(gcf,'Color','w')
        pbaspect([1 1 1])
        xlabel('$t/\tau_r$','FontSize',FontSize,'Interpreter','latex')
        ylabel('$MSD$ ($\ell^2$)','FontSize',FontSize,'Interpreter','latex')
        title(['$\gamma$ = ',num2str(damp,'%.2e'),' [force $\cdot$ time/length], ',...
            num2str(damp*DamperConversion,'%.2e'),' [N s m$^{-1}$]'],...
            'FontSize',FontSize/2,'Interpreter','latex')
        xlim([0 1])

        FileName = [Folder,'/',AddOn,'MSD Subiffusive inset vs t.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'MSD Subiffusive inset vs t.fig'];
        saveas(gcf,FileName)
    end
    writematrix(MSD_ss,MSD_ss_filename)
    writematrix(tau0,tau0_filename)
    writematrix(tauR,taur_filename)
    writematrix(R2,R2_filename)
    writematrix(tau0_2,taur_SD_filename)
    writematrix(tauR_2,tau0_SD_filename)
    writematrix(alpha,alpha_filename)
    writematrix(R2_2,R2_SD_filename)
else
    MSD_ss = readmatrix(MSD_ss_filename);
    tau0 = readmatrix(tau0_filename);
    tauR = readmatrix(taur_filename);
    R2 = readmatrix(R2_filename);
    tau0_2 = readmatrix(taur_SD_filename);
    tauR_2 = readmatrix(tau0_SD_filename);
    alpha = readmatrix(alpha_filename);
    R2_2 = readmatrix(R2_SD_filename);
end
close(wb1)
close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tau0,tauR_eff,R_N] = ComputeEffectiveDamper(N_Kuhn,b,damp)

global DamperConversion LengthConversion

R_N = b;%*sqrt(N_Kuhn/12);
% R_N = b^2;

kbT = 293*1.38e-23;
dampSI = damp*DamperConversion;
D0 = kbT/dampSI;
bSI = b*LengthConversion;
tau0 = (bSI^2)/D0;

% damp_eff = damp;
% tauN_eff = (N_Kuhn^3)*(bSI^2)*dampSI/(24*kbT);
tauR_eff = tau0*N_Kuhn^2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [R2_2,tau0_2,taur_2,alpha] =...
    FitGeneralSubDiffusionModel(time,msd_st_mean,rng,N_Kuhn,tau0,R_N)

% For Subdiffusion
alpha_nom = 0.5;
tau0_nom = tau0;
ct = 0;
NoPts = 13;
x = time(rng);
y = msd_st_mean(rng);
R2_2 = 0;
while R2_2<0.975
    ct = ct+1;
    if ct==1
        sigalpha = 3/4*alpha_nom;
        sigtau0 = 3/4*tau0_nom;
    else
        sigalpha = sigalpha*0.95;
        sigtau0 = sigtau0*0.95;
    end
    alpha_rng = [(linspace(alpha_nom-sigalpha,alpha_nom+sigalpha,NoPts))';alpha_nom];
    tau0_rng = [(linspace(tau0_nom-sigtau0,tau0_nom+sigtau0,NoPts))';tau0_nom];
    alpha_rng = unique(alpha_rng); alpha_rng(alpha_rng<=0) = [];
    tau0_rng = unique(tau0_rng); tau0_rng(tau0_rng<=0) = [];

    Perms = zeros(length(alpha_rng)*length(tau0_rng),2);
    i = 0;
    for ai=1:length(alpha_rng)
        for bi=1:length(tau0_rng)
            i = i+1;
            Perms(i,1) = alpha_rng(ai);
            Perms(i,2) = tau0_rng(bi);
        end
    end
    R2_2_all = zeros(size(Perms,1),1);
    alpha_all = zeros(size(Perms,1),1);
    tau0_all = zeros(size(Perms,1),1);

    for i=1:size(Perms,1)
        alpha_all(i) = Perms(i,1);
        tau0_all(i) = Perms(i,2);
        ft = (R_N^2)*(x/tau0_all(i)).^(alpha_all(i));

        RSS = sum((y-ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_2_all(i) = 1-RSS/TSS;
    end
    diff = abs(R2_2_all-1);
    indx = find(diff==min(diff),1,'first');

    alpha_nom = alpha_all(indx);
    tau0_nom = tau0_all(indx);
    R2_2 = R2_2_all(indx);

    figure(12)
    clf; hold on
    scatter(x,y,'k','filled')
    plot(x,(R_N^2)*(x/tau0_nom).^alpha_nom,'k--')
    xlabel('t')
    ylabel('MSD')
    if ct>100
        break;
    end
end
close(gcf)
taur_2 = tau0_nom*N_Kuhn^2;
tau0_2 = tau0_nom;
alpha = alpha_nom;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R2,tau0,tauR,MSD_ss,rng] = FitRouseModel(time,msd_st_mean,...
    tau0_eff,tauR_eff,R_N,N_Kuhn)  

% global StartMSDMeasurePct
% 
% MSDStartIndx = round(StartMSDMeasurePct/100*length(time));
% time = time(MSDStartIndx:end);
% msd_st_mean = msd_st_mean(MSDStartIndx:end);
% time = time-time(1);

tauR_nom = tauR_eff;
tau0_nom = tau0_eff;
stopindx = find(time>=tauR_nom,1,'first');
if isempty(stopindx)
    stopindx = length(time);
end
rng = (1:stopindx)';

NoPts = 21;
R2 = 0;
ct = 0;
wb5 = waitbar(0,'Fitting Rouse Diffusion Model...');
figure(100)
while R2<0.99
    ct = ct+1;
    waitbar(R2/1,wb5,'Fitting Rouse Diffusion Model...')
    
    if ct==1
        sig_tau0 = 1/8*tau0_nom;
    else
        sig_tau0 = sig_tau0*0.95;
    end
    tau0_rng = [(linspace(tau0_nom-sig_tau0,tau0_nom+sig_tau0,NoPts))';tau0_nom];
    tau0_rng = unique(tau0_rng); tau0_rng(tau0_rng<=0) = [];

    Perms = zeros(length(tau0_rng),1);
    i = 0;
    for ai=1:length(tau0_rng)
        i = i+1;
        Perms(i,1) = tau0_rng(ai);
    end
    R2_all = zeros(size(Perms,1),1);
    tau0_all = zeros(size(Perms,1),1);

    for i=1:size(Perms,1)
        tau0_all(i) = Perms(i,1);

        x = time(rng);
        y = msd_st_mean(rng);

        ft = (R_N^2)*(x/tau0_all(i)).^0.5;

        RSS = sum((y-ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_all(i) = 1-RSS/TSS;
    end
    diff = abs(R2_all-1);
    indx = find(diff==min(diff),1,'first');

    tau0_nom = tau0_all(indx);
    R2 = R2_all(indx);

    clf; hold on
    scatter(x,y,'k','filled')
    plot(x,(R_N^2)*(x/tau0_nom).^0.5,'k--')
    xlabel('t')
    ylabel('MSD')
    if ct>50
        break;
    end
end
tau0 = tau0_nom;
tauR = tau0*N_Kuhn^2;
figure(100); close
close(wb5)

strt = find(time>=tauR_nom,1,'first');
if isempty(strt)
    strt = length(time)-1;
end
MSD_ss = mean(msd_st_mean(strt:end));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotHistograms(rx,ry,rz,bondtypes,time,stretch,...
    WRTStretchOrLength,HistoUnits,OneDOnly)

global binMax binMin binWidth FontSize ConstantsFileName xlab...
    MovieName1 MovieName2 MovieName3 MovieName4 MovieName5 MovieName6...
    MovieName7 MovieName8 MovieName9 MovieName10 MovieName11 MovieFolder

CheckBackBoneVsStickers = 0;
norms = (rx.^2+ry.^2+rz.^2).^(1/2);

% FileName = ['Matlab Topological Data/Constants.S1',FileTagCompiled,'.txt'];
% 
% Table = readtable(FileName);
Table = readtable(ConstantsFileName);
constants = table2array(Table);
N_Kuhn = constants(3); b = constants(4);

%Histo inputs
L = N_Kuhn*b;
binMax = L;   %nm                   %Histogram inputs
binMin = -binMax;
Nbins = 101;
binWidth = (binMax-binMin)/Nbins;

if ~isfolder([MovieFolder,'/End-to-end'])
    mkdir([MovieFolder,'/End-to-end'])
end
v1 = VideoWriter(MovieName1); v1.FrameRate = 10; open(v1)
if OneDOnly~=1
    v2 = VideoWriter(MovieName2); v2.FrameRate = 10; open(v2)
    v3 = VideoWriter(MovieName3); v3.FrameRate = 10; open(v3)
    v4 = VideoWriter(MovieName4); v4.FrameRate = 10; open(v4)
    v5 = VideoWriter(MovieName5); v5.FrameRate = 10; open(v5)

    if CheckBackBoneVsStickers==1
        v6 = VideoWriter(MovieName6); v6.FrameRate = 10; open(v6)
        v7 = VideoWriter(MovieName7); v7.FrameRate = 10; open(v7)
        v8 = VideoWriter(MovieName8); v8.FrameRate = 10; open(v8)
        v9 = VideoWriter(MovieName9); v9.FrameRate = 10; open(v9)
        v10 = VideoWriter(MovieName10); v10.FrameRate = 10; open(v10)
        v11 = VideoWriter(MovieName11); v11.FrameRate = 10; open(v11)
    end
end

wb3 = waitbar(0,'Plotting Histograms...');
i = 0; iout = 10;
stretchindx = find(stretch>1,1,'first');
while i+iout<=size(rx,1)
    i = i+iout;
    if i+iout>=stretchindx
        iout = 1;
        stretchindx = find(stretch>=4,1,'first');
    end
    if i+iout>=stretchindx+20
        iout=10;
    end
    %     for i=1:iout:size(rx,1)
    waitbar(i/size(rx,1),wb3,'Plotting Histograms...')

    % Define temp data
    bondtypes_temp = (bondtypes(i,:))';
    norms_temp = (norms(i,:))';
    rx_temp = (rx(i,:))';
    ry_temp = (ry(i,:))';
    rz_temp = (rz(i,:))';
    rx_temp(isnan(norms_temp)) = [];
    ry_temp(isnan(norms_temp)) = [];
    rz_temp(isnan(norms_temp)) = [];
    norms_temp(isnan(norms_temp)) = [];

    figure(201); clf; hold on
    XLim = sqrt(N_Kuhn);
    YLim = 1.5;
    color_t = 'k';
    SymDist = 0;
    PlotIdeal = 1;
    Plot1DHistogram(norms_temp,color_t,0.6,XLim,YLim,N_Kuhn,b,...
        SymDist,PlotIdeal,WRTStretchOrLength,HistoUnits)

    set(gca,'FontSize',FontSize/1.5)
    xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
    ylabel('$p$','FontSize',FontSize,'Interpreter','latex')

    title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
        ', $t$ = ',num2str(time(i),'%.2e'),' s'],'FontSize',FontSize/2,...
        'Interpreter','latex')

    F1 = getframe(gcf);
    writeVideo(v1,F1);
    mov1(i) = F1;

    if OneDOnly~=1
        figure(202); clf; hold on
        XLim = sqrt(N_Kuhn);
        YLim = 0.2;
        SymDist = 0;
        PlotIdeal = 0;
        color_t = 'k';
        Plot1DHistogram(rx_temp,color_t,0.6,XLim,YLim,N_Kuhn,b,...
            SymDist,PlotIdeal,WRTStretchOrLength,HistoUnits)
        color_t = 'g';
        Plot1DHistogram(ry_temp,color_t,0.6,XLim,YLim,N_Kuhn,b,...
            SymDist,PlotIdeal,WRTStretchOrLength,HistoUnits)
        color_t = 'm';
        Plot1DHistogram(rz_temp,color_t,0.6,XLim,YLim,N_Kuhn,b,...
            SymDist,PlotIdeal,WRTStretchOrLength,HistoUnits)

        xlabel('$\lambda$','FontSize',FontSize,'Interpreter','latex')
        ylabel('$p$','FontSize',FontSize,'Interpreter','latex')

        title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
            ', $t$ = ',num2str(time(i),'%.2e'),' s'],'FontSize',FontSize/2,...
            'Interpreter','latex')

        l = legend('$|r_x|$','$|r_y|$','$|r_z|$');
        l.Interpreter = 'latex'; l.FontSize = FontSize/2;
        l.Location = 'NorthEast';

        F2 = getframe(gcf);
        writeVideo(v2,F2);
        mov2(i) = F2;

        figure(203);

        Plot2DHistogram([rx_temp ry_temp],N_Kuhn,b)

        xlabel('$\lambda_x$','FontSize',FontSize,'Interpreter','latex')
        ylabel('$\lambda_y$','FontSize',FontSize,'Interpreter','latex')

        title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
            ', $t$ = ',num2str(time(i),'%.2e'),' s'],'FontSize',FontSize/2,...
            'Interpreter','latex')

        daspect([1 1 1])
        set(gcf,'Color','w')
        box on

        F3 = getframe(gcf);
        writeVideo(v3,F3);
        mov3(i) = F3;

        figure(204);

        Plot2DHistogram([ry_temp rz_temp],N_Kuhn,b)

        xlabel('$\lambda_y$','FontSize',FontSize,'Interpreter','latex')
        ylabel('$\lambda_z$','FontSize',FontSize,'Interpreter','latex')

        title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
            ', $t$ = ',num2str(time(i),'%.2e'),' s'],'FontSize',FontSize/2,...
            'Interpreter','latex')

        daspect([1 1 1])
        set(gcf,'Color','w')
        box on

        F4 = getframe(gcf);
        writeVideo(v4,F4);
        mov4(i) = F4;

        figure(205);

        Plot2DHistogram([rx_temp rz_temp],N_Kuhn,b)

        xlabel('$\lambda_x$','FontSize',FontSize,'Interpreter','latex')
        ylabel('$\lambda_z$','FontSize',FontSize,'Interpreter','latex')

        title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
            ', $t$ = ',num2str(time(i),'%.2e'),' s'],'FontSize',FontSize/2,...
            'Interpreter','latex')

        daspect([1 1 1])
        set(gcf,'Color','w')
        box on

        F5 = getframe(gcf);
        writeVideo(v5,F5);
        mov5(i) = F5;

        if CheckBackBoneVsStickers==1
            % Branches
            figure(206); clf; hold on
            XLim = sqrt(N_Kuhn);
            YLim = 0.2;
            color_t = 'k';
            SymDist = 0;
            PlotIdeal = 1;
            Plot1DHistogram(norms_temp(bondtypes_temp==2),color_t,0.6,...
                XLim,YLim,N_Kuhn,b,SymDist,PlotIdeal,...
                WRTStretchOrLength,HistoUnits)

            set(gca,'FontSize',FontSize/1.5)
            xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
            ylabel('$p$','FontSize',FontSize,'Interpreter','latex')

            title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
                ', $t$ = ',num2str(time(i),'%.2e'),' s'],'FontSize',FontSize/2,...
                'Interpreter','latex')

            F6 = getframe(gcf);
            writeVideo(v6,F6);
            mov6(i) = F6;

            figure(207);

            Plot2DHistogram([rx_temp(bondtypes_temp==2) ry_temp(bondtypes_temp==2)],N_Kuhn,b)

            xlabel('$\lambda_x$','FontSize',FontSize,'Interpreter','latex')
            ylabel('$\lambda_y$','FontSize',FontSize,'Interpreter','latex')

            title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
                ', $t$ = ',num2str(time(i),'%.2e'),' s'],'FontSize',FontSize/2,...
                'Interpreter','latex')

            daspect([1 1 1])
            set(gcf,'Color','w')
            box on

            F7 = getframe(gcf);
            writeVideo(v7,F7);
            mov7(i) = F7;

            figure(208);

            Plot2DHistogram([ry_temp(bondtypes_temp==2) rz_temp(bondtypes_temp==2)],N_Kuhn,b)

            xlabel('$\lambda_y$','FontSize',FontSize,'Interpreter','latex')
            ylabel('$\lambda_z$','FontSize',FontSize,'Interpreter','latex')

            title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
                ', $t$ = ',num2str(time(i),'%.2e'),' s'],'FontSize',FontSize/2,...
                'Interpreter','latex')

            daspect([1 1 1])
            set(gcf,'Color','w')
            box on

            F8 = getframe(gcf);
            writeVideo(v8,F8);
            mov8(i) = F8;

            % Backbones
            figure(209); clf; hold on
            XLim = sqrt(N_Kuhn);
            YLim = 0.2;
            color_t = 'k';
            SymDist = 0;
            PlotIdeal = 1;
            Plot1DHistogram(norms_temp(bondtypes_temp==1),color_t,0.6,...
                XLim,YLim,N_Kuhn,b,SymDist,PlotIdeal,...
                WRTStretchOrLength,HistoUnits)

            set(gca,'FontSize',FontSize/1.5)
            xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
            ylabel('$p$','FontSize',FontSize,'Interpreter','latex')

            title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
                ', $t$ = ',num2str(time(i),'%.2e'),' s'],'FontSize',FontSize/2,...
                'Interpreter','latex')

            F9 = getframe(gcf);
            writeVideo(v9,F9);
            mov9(i) = F9;

            figure(210)
            Plot2DHistogram([rx_temp(bondtypes_temp==1) ry_temp(bondtypes_temp==1)],N_Kuhn,b)

            xlabel('$\lambda_x$','FontSize',FontSize,'Interpreter','latex')
            ylabel('$\lambda_y$','FontSize',FontSize,'Interpreter','latex')

            title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
                ', $t$ = ',num2str(time(i),'%.2e'),' s'],'FontSize',FontSize/2,...
                'Interpreter','latex')

            daspect([1 1 1])
            set(gcf,'Color','w')
            box on

            F10 = getframe(gcf);
            writeVideo(v10,F10);
            mov10(i) = F10;

            figure(211);

            Plot2DHistogram([ry_temp(bondtypes_temp==1) rz_temp(bondtypes_temp==1)],N_Kuhn,b)

            xlabel('$\lambda_y$','FontSize',FontSize,'Interpreter','latex')
            ylabel('$\lambda_z$','FontSize',FontSize,'Interpreter','latex')

            title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
                ', $t$ = ',num2str(time(i),'%.2e'),' s'],'FontSize',FontSize/2,...
                'Interpreter','latex')

            daspect([1 1 1])
            set(gcf,'Color','w')
            box on

            F11 = getframe(gcf);
            writeVideo(v11,F11);
            mov11(i) = F11;
        end
    end
end
close(wb3)
close(v1);
if OneDOnly~=1

    close(v2);
    close(v3);
    close(v4);
    close(v5);
    if CheckBackBoneVsStickers==1
        close(v6);
        close(v7);
        close(v8);
        close(v9);
        close(v10);
        close(v11);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plot2DHistogram(HistEndtoEnd,N_Kuhn,b)
% Plot2DHistogram(HistEndtoEnd,v_hat,g_hat)

global binMax binMin binWidth FontSize

Edges = {binMin/(sqrt(N_Kuhn)*b):binWidth/(sqrt(N_Kuhn)*b):binMax/(sqrt(N_Kuhn)*b)...
    binMin/(sqrt(N_Kuhn)*b):binWidth/(sqrt(N_Kuhn)*b):binMax/(sqrt(N_Kuhn)*b)};

Lambda = HistEndtoEnd./(ones(size(HistEndtoEnd))*sqrt(N_Kuhn)*b);
hist3(Lambda,'CdataMode','auto',...
    'FaceColor','interp',...
    'Edges',Edges,'LineWidth',0.25);
hold on

% if PlotEigenVectors==1
%     P1 = [0 0 1000]; P2 = [v_hat(1,1)*g_hat(1,1) v_hat(2,1)*g_hat(1,1) 1000];
%     dP = P2-P1;
%     q = quiver3(P1(1),P1(2),P1(3),dP(1),dP(2),dP(3),0);
%     q.LineWidth = 1.5;
%     q.Color = 'k';
%     q.MaxHeadSize = 0.5;
%     
%     P1 = [0 0 1000]; P2 = [v_hat(1,2)*g_hat(2,2) v_hat(2,2)*g_hat(2,2) 1000];
%     dP = P2-P1;
%     q = quiver3(P1(1),P1(2),P1(3),dP(1),dP(2),dP(3),0);
%     q.LineWidth = 1.5;
%     q.Color = 'k';
%     q.MaxHeadSize = 0.5;
% end
% 
% %Plot TNT Data for Full Stress Case
% TNTFileName = ['Stress Relaxation/TNT Data/Full W.',...
%     num2str(W,'%.3f'),'.EF.',...
%     num2str(EdgeFactor,'%.2f'),...
%     '.txt'];
% if isfile(TNTFileName)
%     TNTData = table2array(readtable(TNTFileName));
%     Time = TNTData(:,1);
%     Diff = abs(Time-time);
%     Indx = find(Diff==min(Diff));
%     mu22 = TNTData(Indx,4);
%     mu11 = TNTData(Indx,6);
%     
%     %Plot TNT mu as ellipse
%     x = linspace(-mu11,mu11,1000);
%     yPos = mu22*(1-(x/mu11).^2).^0.5;
%     yNeg = -yPos;
%     p1 = plot3(x,yPos,750*ones(size(x)),'r');
%     p1.LineWidth = 2;
%     p1 = plot3(x,yNeg,750*ones(size(x)),'r');
%     p1.LineWidth = 2;
%     
%     p1 = plot3(0.5*x,0.5*yPos,750*ones(size(x)),'r');
%     p1.LineWidth = 2;
%     p1 = plot3(0.5*x,0.5*yNeg,750*ones(size(x)),'r');
%     p1.LineWidth = 2;
%     
%     p1 = plot3(1.5*x,1.5*yPos,750*ones(size(x)),'r');
%     p1.LineWidth = 2;
%     p1 = plot3(1.5*x,1.5*yNeg,750*ones(size(x)),'r');
%     p1.LineWidth = 2;
% end

% xlim([binMin binMax]);
% ylim([binMin binMax]);
LIM = sqrt(N_Kuhn);
xlim([-LIM LIM]);
ylim([-LIM LIM]);
zlim([1 Inf])

view(0,90)
pbaspect([1 1 1])
box on 
set(gcf,'Color','w')
set(gca,'FontSize',FontSize/1.5)

% LineWidth = 0.75;
% Color = 'k';
% [~] = circle(0,0,binMax,Color,LineWidth);
% Color = 'r';
% [~] = circle(0,0,0.95*binMax,Color,LineWidth);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plot1DHistogram(Norms,Color,Alpha,XLim,YLim,N_Kuhn,b,...
    SymDist,PlotIdeal,WRTStretchOrLength,HistoUnits)

global FontSize LengthConversion %binMax binMin binWidth 

if WRTStretchOrLength==0
    binMax = 4;
    binMin = 0;
    binWidth = 0.08;%0.0396;
    Edges = (binMin:binWidth:binMax);
%     Edges = (binMin:binWidth:binMax)/(sqrt(N_Kuhn)*b);
    XLim = 4;
else
    binMax = 4;
    binMin = 0;
    binWidth = 0.0396;
    XLim = binMax/2;
    Edges = (binMin:binWidth:binMax);
%     XLim = XLim*(sqrt(N_Kuhn)*b);
    if HistoUnits==1
        b = b*LengthConversion;
        Norms = Norms*LengthConversion;
        XLim = XLim*LengthConversion;
        Edges = Edges*LengthConversion;
    end
end

hold on

Lambda = Norms/(sqrt(N_Kuhn)*b);
if WRTStretchOrLength==0
    h = histogram(Lambda,'Normalization','pdf','BinEdges',Edges);
else
    h = histogram(Norms,'Normalization','pdf','BinEdges',Edges);
end
h.FaceColor = Color;
h.FaceAlpha = Alpha;
h.EdgeColor = 'k';
h.EdgeAlpha = 1;
Bins = (h.BinEdges)';

% xlim([0 sqrt(N_Kuhn)])
if SymDist==1
    xlim([-XLim XLim])
else
    xlim([0 XLim])
end
if WRTStretchOrLength==0
    ylim([0 YLim])
else
    ylim([0 2])
    if HistoUnits==1
        ylim([0 1])
    end
end

set(gca,'FontSize',FontSize/1.5);
set(gcf,'Color','w')
pbaspect([3 1 1])

if PlotIdeal==1
    r = (linspace(0,10*N_Kuhn*b,10000))';%*sqrt(N_Kuhn)*b;
    lam = r/(sqrt(N_Kuhn)*b);
    if WRTStretchOrLength==0
        X = lam;
        xDiscrete = (Bins(1:end-1)+Bins(2:end))/2*sqrt(N_Kuhn)*b;
        yshift = 1;
    else
        X = r;
        xDiscrete = (Bins(1:end-1)+Bins(2:end))/2;
        yshift = 1;
        if HistoUnits==1
            yshift = 0.33;
        end
    end
    sigma = sqrt(N_Kuhn/3)*b;
    P = 4*pi*(r.^2)*((sigma*sqrt(2*pi))^(-3)).*exp(-1/2*(r/sigma).^2);

    Area = trapz(lam,P);
    P = P/Area;
    p = plot(X,P,'k:');
    p.LineWidth = 1.5;
    
    P = P*yshift;
    p = plot(X,P,'k-');
    p.LineWidth = 1.5;
    
    yDiscrete = (h.Values)';
    yIdeal = interp1(r,P,xDiscrete);
    yDiscrete(xDiscrete<=0) = [];
    yIdeal(xDiscrete<=0) = [];
    xDiscrete(xDiscrete<=0) = [];

    RSS = sum((yDiscrete-yIdeal).^2);
    TSS = sum((yIdeal-mean(yIdeal)).^2);
    R2 = 1-RSS/TSS;                %For Rouse diffusion

    text(1/2*sqrt(N_Kuhn),0.75,['$R^2$ = ',num2str(R2,'%.2f')],'FontSize',FontSize/2,'Interpreter','latex')
    CheckFig = 0;
    if CheckFig==1
        figure(10); clf; hold on
        scatter(xDiscrete/(sqrt(N_Kuhn)*b),yDiscrete,'k')
        plot(lam,P,'k')
        xlim([0 sqrt(N_Kuhn)])
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,f] = PlotCurve(X,Y,Err,Color,Style)

global LineWidth

ErrX = [X;flipud(X)];
ErrY = [Y+Err;flipud(Y-Err)];
f = fill(ErrX,ErrY,Color);
if ~isempty(f)
    f.FaceAlpha = 0.25;
    f.LineStyle = 'none';
end

p = plot(X,Y);
p.Color = Color;
p.LineWidth = LineWidth;
p.LineStyle = Style;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mean,SE] = CalcMeanAndSE(Data,Dim)

Mean = nanmean(Data,Dim);
SE = nanstd(Data,0,Dim)/sqrt(size(Data,Dim));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time,stretch,rx_ensemble,ry_ensemble,rz_ensemble,bondtypes_ensemble,...
    msd_st_mean,msd_st_err,...
    g11_st_mean,g22_st_mean,g33_st_mean,...
    g11_st_err,g22_st_err,g33_st_err] = ....
    AssembleTheDataMeso(Corners,EdgeFactor,Np,N,ka_in,kd_in,f0,dt,damp,...
    N_Kuhn,b,Separation,PackageTemp,Override,N_steps,N_samples)

global BondsFileName AtomsFileName...
    RxFileName RyFileName RzFileName TimeStretchFileName...
    BondTypesFileName MSDStFileName MSDSt_errFileName...
    g11_stFileName g22_stFileName g33_stFileName...
    g12_stFileName g23_stFileName g31_stFileName...
    g11_st_errFileName g22_st_errFileName g33_st_errFileName...
    g12_st_errFileName g23_st_errFileName g31_st_errFileName...
    StartMSDMeasurePct


% Allocate
dx0 = Corners(1,5)-Corners(1,2);
stretch = (Corners(:,5)-Corners(:,2))/dx0;
time = Corners(:,1)*dt;

StartMSDMeasurePct = 0.001;

if ~isfile(RxFileName) || ~isfile(RyFileName) ||...
        ~isfile(RzFileName) || Override==1

    % Allocate
    msd_st_mean = zeros(N_steps,N_samples);
    msd_st_err = zeros(N_steps,N_samples);

    g11_st_mean = zeros(N_steps,N_samples);
    g22_st_mean = zeros(N_steps,N_samples);
    g33_st_mean = zeros(N_steps,N_samples);
    g12_st_mean = zeros(N_steps,N_samples);
    g23_st_mean = zeros(N_steps,N_samples);
    g31_st_mean = zeros(N_steps,N_samples); 

    g11_st_err = zeros(N_steps,N_samples);
    g22_st_err = zeros(N_steps,N_samples);
    g33_st_err = zeros(N_steps,N_samples);
    g12_st_err = zeros(N_steps,N_samples);
    g23_st_err = zeros(N_steps,N_samples);
    g31_st_err = zeros(N_steps,N_samples); 

    %% PARALLELIZE THIS%%
    for n=1:size(PackageTemp,1)

        %% Unpack Swept Input Parameters
        Sample = PackageTemp(n,1);    %Sample index

        %% Callout Input Script
        % Initializes non-sweeping parameters
        InputScript(EdgeFactor,Sample,Np,N,ka_in,kd_in,f0,dt,damp,...
            N_Kuhn,b,Separation);

        %% Make diriectories and set filenames
        SetDirAndFileNames;
        DefineCompiledFileNames;

        %% Import Bonds Data
        BondData = dlmread(BondsFileName);

        %% Import Atoms Data
        AtomData = dlmread(AtomsFileName);
        % timestep id bond_type node1 node2 rx ry rz norm force
        timesteps = unique(BondData(:,1));
        if size(Corners,1)~=length(timesteps)
            Diff = size(Corners,1)-length(timesteps);
            if Diff<0
                timesteps(end+Diff+1:end) = [];
            end
        end

        MSDStartIndx = ceil(StartMSDMeasurePct/100*length(timesteps));
        MSDStartStep = timesteps(MSDStartIndx);
        %% Initialize dynamics variables
        wb2 = waitbar(0,'Computing stress & dynamics data...');
        for i=1:length(timesteps)
            timestep = timesteps(i);
            waitbar(i/length(timesteps),wb2,'Computing stress & dynamics data...')

            %% Compute virial stress & bond concentration
            rs = BondData(BondData(:,1)==timestep,6:8);

            bondtypes = BondData(BondData(:,1)==timestep,3);
                        
            %% Compute MSDs (by atom type) - can I back out an effective viscosity from this?
            type = 2;
            if timestep<MSDStartStep
                msd_st_mean(i,n) = 0;
                msd_st_err(i,n) = 0;
            else
                [msd_st_mean(i,n),msd_st_err(i,n)] = ...
                    ComputeMSDs(timestep,AtomData,Corners(i,:),MSDStartStep,type);
            end

            %% Store end-to-end components
            rx = rs(:,1); ry = rs(:,2); rz = rs(:,3);
            if i==1 && n==1
                AtomsTemp = AtomData(AtomData(:,1)==timestep,:);
                N_stickers = length(AtomsTemp(AtomsTemp(:,end)==2));
                MaxBonds = size(rs,1)+round(N_stickers/2);
                N_col = MaxBonds*size(PackageTemp,1);
                N_time = length(timesteps);
                rx_ensemble = zeros(length(N_time),N_col);
                ry_ensemble = zeros(length(N_time),N_col);
                rz_ensemble = zeros(length(N_time),N_col);
                bondtypes_ensemble = zeros(length(N_time),N_col); 
            end
            strt = MaxBonds*n-MaxBonds+1;
            fnsh = strt+size(rs,1)-1;
            rx_ensemble(i,strt:fnsh) = rx;
            ry_ensemble(i,strt:fnsh) = ry;
            rz_ensemble(i,strt:fnsh) = rz;
            bondtypes_ensemble(i,strt:fnsh) = bondtypes;

            %% Compute metric tensors by bond type
            [g11_st_mean(i,n),g22_st_mean(i,n),g33_st_mean(i,n),...
                g12_st_mean(i,n),g23_st_mean(i,n),g31_st_mean(i,n),...
                g11_st_err(i,n),g22_st_err(i,n),g33_st_err(i,n),...
                g12_st_err(i,n),g23_st_err(i,n),g31_st_err(i,n)] = ...
                ComputeMetricTensor(rx,ry,rz);
        end
        close(wb2)
    end
    norm_ensemble = (rx_ensemble.^2 + ry_ensemble.^2 ...
        + rz_ensemble.^2).^(1/2);
    rx_ensemble(norm_ensemble==0) = NaN;
    ry_ensemble(norm_ensemble==0) = NaN;
    rz_ensemble(norm_ensemble==0) = NaN;
    bondtypes_ensemble(norm_ensemble==0) = NaN;

    %% Save compiled data
    dlmwrite(RxFileName,rx_ensemble);
    dlmwrite(RyFileName,ry_ensemble);
    dlmwrite(RzFileName,rz_ensemble);

    dlmwrite(BondTypesFileName,bondtypes_ensemble);

    dlmwrite(MSDStFileName,msd_st_mean);
    dlmwrite(MSDSt_errFileName,msd_st_err);

    dlmwrite(g11_stFileName,g11_st_mean);
    dlmwrite(g22_stFileName,g22_st_mean);
    dlmwrite(g33_stFileName,g33_st_mean);
    dlmwrite(g12_stFileName,g12_st_mean);
    dlmwrite(g23_stFileName,g23_st_mean);
    dlmwrite(g31_stFileName,g31_st_mean);

    dlmwrite(g11_st_errFileName,g11_st_err);
    dlmwrite(g22_st_errFileName,g22_st_err);
    dlmwrite(g33_st_errFileName,g33_st_err);
    dlmwrite(g12_st_errFileName,g12_st_err);
    dlmwrite(g23_st_errFileName,g23_st_err);
    dlmwrite(g31_st_errFileName,g31_st_err);

    dlmwrite(TimeStretchFileName,[time stretch])
else
    rx_ensemble=dlmread(RxFileName);
    ry_ensemble=dlmread(RyFileName);
    rz_ensemble=dlmread(RzFileName);

    bondtypes_ensemble=dlmread(BondTypesFileName);

    msd_st_mean=dlmread(MSDStFileName);
    msd_st_err=dlmread(MSDSt_errFileName);

    dat = dlmread(TimeStretchFileName);

    g11_st_mean=dlmread(g11_stFileName);
    g22_st_mean=dlmread(g22_stFileName);
    g33_st_mean=dlmread(g33_stFileName);
    g12_st_mean=dlmread(g12_stFileName);
    g23_st_mean=dlmread(g23_stFileName);
    g31_st_mean=dlmread(g31_stFileName);

    g11_st_err=dlmread(g11_st_errFileName);
    g22_st_err=dlmread(g22_st_errFileName);
    g33_st_err=dlmread(g33_st_errFileName);
    g12_st_err=dlmread(g12_st_errFileName);
    g23_st_err=dlmread(g23_st_errFileName);
    g31_st_err=dlmread(g31_st_errFileName);

    time = dat(:,1); stretch = dat(:,2);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time,stretch,rx_ensemble,ry_ensemble,rz_ensemble,bondtypes_ensemble,...
    msd_st_mean,msd_st_err,...
    g11_st_mean,g22_st_mean,g33_st_mean,...
    g11_st_err,g22_st_err,g33_st_err] = ....
    AssembleTheDataBeadSpring(Corners,EdgeFactor,Np,N,ka_in,kd_in,f0,dt,damp,...
                N_Kuhn,b,Separation,PackageTemp,Override,N_steps,N_samples)

global BondsFileName AtomsFileName...
    RxFileName RyFileName RzFileName TimeStretchFileName...
    BondTypesFileName MSDStFileName MSDSt_errFileName...
    g11_stFileName g22_stFileName g33_stFileName...
    g12_stFileName g23_stFileName g31_stFileName...
    g11_st_errFileName g22_st_errFileName g33_st_errFileName...
    g12_st_errFileName g23_st_errFileName g31_st_errFileName...
    StartMSDMeasurePct


% Allocate
dx0 = Corners(1,5)-Corners(1,2);
stretch = (Corners(:,5)-Corners(:,2))/dx0;
time = Corners(:,1)*dt;

StartMSDMeasurePct = 0.001;

if ~isfile(RxFileName) || ~isfile(RyFileName) ||...
        ~isfile(RzFileName) || Override==1

    % Allocate
%     msd_th_mean = zeros(N_steps,N_samples);
%     msd_th_err = zeros(N_steps,N_samples);
    msd_st_mean = zeros(N_steps,N_samples);
    msd_st_err = zeros(N_steps,N_samples);

    g11_st_mean = zeros(N_steps,N_samples);
    g22_st_mean = zeros(N_steps,N_samples);
    g33_st_mean = zeros(N_steps,N_samples);
    g12_st_mean = zeros(N_steps,N_samples);
    g23_st_mean = zeros(N_steps,N_samples);
    g31_st_mean = zeros(N_steps,N_samples); 

    g11_st_err = zeros(N_steps,N_samples);
    g22_st_err = zeros(N_steps,N_samples);
    g33_st_err = zeros(N_steps,N_samples);
    g12_st_err = zeros(N_steps,N_samples);
    g23_st_err = zeros(N_steps,N_samples);
    g31_st_err = zeros(N_steps,N_samples); 

    %% PARALLELIZE THIS%%
    for n=1:size(PackageTemp,1)

        %% Unpack Swept Input Parameters
        Sample = PackageTemp(n,1);    %Sample index

        %% Callout Input Script
        % Initializes non-sweeping parameters
        InputScript(EdgeFactor,Sample,Np,N,ka_in,kd_in,f0,dt,damp,...
            N_Kuhn,b,Separation);

        %% Make diriectories and set filenames
        SetDirAndFileNames;
        DefineCompiledFileNames;

        %% Import Bonds Data
        BondData = dlmread(BondsFileName);

        %% Import Atoms Data
        AtomData = dlmread(AtomsFileName);
        % timestep id bond_type node1 node2 rx ry rz norm force
        timesteps = unique(BondData(:,1));
        if size(Corners,1)~=length(timesteps)
            Diff = size(Corners,1)-length(timesteps);
            if Diff<0
                timesteps(end+Diff+1:end) = [];
            end
        end

        MSDStartIndx = ceil(StartMSDMeasurePct/100*length(timesteps));
        MSDStartStep = timesteps(MSDStartIndx);
        %% Initialize dynamics variables
        wb2 = waitbar(0,'Computing stress & dynamics data...');
        for i=1:length(timesteps)
            timestep = timesteps(i);
            waitbar(i/length(timesteps),wb2,'Computing stress & dynamics data...')

            %% Isolate important node data
%             rs = BondData(BondData(:,1)==timestep,6:8);

            bondtypes = BondData(BondData(:,1)==timestep,3);
            atomtypes = AtomData(AtomData(:,1)==timestep,end);
            atomnumbers = AtomData(AtomData(:,1)==timestep,2);
            atompositions = AtomData(AtomData(:,1)==timestep,3:5);

            NChains = Np*2;
            rs = zeros(NChains,3);
            for chainnumber=1:NChains
                tether = chainnumber;
                sticker = NChains + N_Kuhn*chainnumber;

                type_t = atomtypes(atomnumbers==tether,1);
                type_s = atomtypes(atomnumbers==sticker,1);
                if type_t~=1 || type_s~=3
                    error('Check atom types')
                end

                pos_t = atompositions(atomnumbers==tether,:);
                pos_s = atompositions(atomnumbers==sticker,:);
                rs(chainnumber,:) = pos_s-pos_t;
            end
                        
            %% Compute MSDs (by atom type) - can I back out an effective viscosity from this?
            if timestep<MSDStartStep
%                 msd_th_mean(i,n) = 0;
%                 msd_th_err(i,n) = 0;
                msd_st_mean(i,n) = 0;
                msd_st_err(i,n) = 0;
            else
%                 type = 1;   %Tether groups (MSD should always be zero)
%                 [msd_th_mean(i,n),msd_th_err(i,n)] = ...
%                     ComputeMSDs(timestep,AtomData,Corners(i,:),MSDStartStep,type);
                type = 3;   %End groups (MSD should follow Rouse diffusion)
                [msd_st_mean(i,n),msd_st_err(i,n)] = ...
                    ComputeMSDs(timestep,AtomData,Corners(i,:),MSDStartStep,type);
            end

            %% Store end-to-end components
            rx = rs(:,1); ry = rs(:,2); rz = rs(:,3);
            if i==1 && n==1
                AtomsTemp = AtomData(AtomData(:,1)==timestep,:);
                N_stickers = length(AtomsTemp(AtomsTemp(:,end)==2));
                MaxBonds = size(rs,1)+round(N_stickers/2);
                N_col = MaxBonds*size(PackageTemp,1);
                N_time = length(timesteps);
                rx_ensemble = zeros(length(N_time),N_col);
                ry_ensemble = zeros(length(N_time),N_col);
                rz_ensemble = zeros(length(N_time),N_col);
                bondtypes_ensemble = zeros(length(N_time),N_col); 
            end
            strt = MaxBonds*n-MaxBonds+1;
            fnsh = strt+size(rs,1)-1;
            rx_ensemble(i,strt:fnsh) = rx;
            ry_ensemble(i,strt:fnsh) = ry;
            rz_ensemble(i,strt:fnsh) = rz;
%             bondtypes_ensemble(i,strt:fnsh) = bondtypes;
            bondtypes_ensemble(i,strt:fnsh) = 4;        %4 for ensemble spring

            %% Compute metric tensors by bond type
            [g11_st_mean(i,n),g22_st_mean(i,n),g33_st_mean(i,n),...
                g12_st_mean(i,n),g23_st_mean(i,n),g31_st_mean(i,n),...
                g11_st_err(i,n),g22_st_err(i,n),g33_st_err(i,n),...
                g12_st_err(i,n),g23_st_err(i,n),g31_st_err(i,n)] = ...
                ComputeMetricTensor(rx,ry,rz);
        end
        close(wb2)
    end
    norm_ensemble = (rx_ensemble.^2 + ry_ensemble.^2 ...
        + rz_ensemble.^2).^(1/2);
    rx_ensemble(norm_ensemble==0) = NaN;
    ry_ensemble(norm_ensemble==0) = NaN;
    rz_ensemble(norm_ensemble==0) = NaN;
    bondtypes_ensemble(norm_ensemble==0) = NaN;

    %% Save compiled data
    dlmwrite(RxFileName,rx_ensemble);
    dlmwrite(RyFileName,ry_ensemble);
    dlmwrite(RzFileName,rz_ensemble);

    dlmwrite(BondTypesFileName,bondtypes_ensemble);

    dlmwrite(MSDStFileName,msd_st_mean);
    dlmwrite(MSDSt_errFileName,msd_st_err);

    dlmwrite(g11_stFileName,g11_st_mean);
    dlmwrite(g22_stFileName,g22_st_mean);
    dlmwrite(g33_stFileName,g33_st_mean);
    dlmwrite(g12_stFileName,g12_st_mean);
    dlmwrite(g23_stFileName,g23_st_mean);
    dlmwrite(g31_stFileName,g31_st_mean);

    dlmwrite(g11_st_errFileName,g11_st_err);
    dlmwrite(g22_st_errFileName,g22_st_err);
    dlmwrite(g33_st_errFileName,g33_st_err);
    dlmwrite(g12_st_errFileName,g12_st_err);
    dlmwrite(g23_st_errFileName,g23_st_err);
    dlmwrite(g31_st_errFileName,g31_st_err);

    dlmwrite(TimeStretchFileName,[time stretch])
else
    rx_ensemble=dlmread(RxFileName);
    ry_ensemble=dlmread(RyFileName);
    rz_ensemble=dlmread(RzFileName);

    bondtypes_ensemble=dlmread(BondTypesFileName);

    msd_st_mean=dlmread(MSDStFileName);
    msd_st_err=dlmread(MSDSt_errFileName);

    dat = dlmread(TimeStretchFileName);

    g11_st_mean=dlmread(g11_stFileName);
    g22_st_mean=dlmread(g22_stFileName);
    g33_st_mean=dlmread(g33_stFileName);
    g12_st_mean=dlmread(g12_stFileName);
    g23_st_mean=dlmread(g23_stFileName);
    g31_st_mean=dlmread(g31_stFileName);

    g11_st_err=dlmread(g11_st_errFileName);
    g22_st_err=dlmread(g22_st_errFileName);
    g33_st_err=dlmread(g33_st_errFileName);
    g12_st_err=dlmread(g12_st_errFileName);
    g23_st_err=dlmread(g23_st_errFileName);
    g31_st_err=dlmread(g31_st_errFileName);

    time = dat(:,1); stretch = dat(:,2);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g11,g22,g33,g12,g23,g31,...
    g11_err,g22_err,g33_err,...
    g12_err,g23_err,g31_err] = ...
    ComputeMetricTensor(rx,ry,rz)

norms = (rx.^2 +ry.^2 +rz.^2).^(1/2);
rx_hat = rx./norms; ry_hat = ry./norms; rz_hat = rz./norms;

g11_temp = rx_hat.*rx_hat;
g22_temp = ry_hat.*ry_hat;
g33_temp = rz_hat.*rz_hat;
g12_temp = rx_hat.*ry_hat;
g23_temp = ry_hat.*rz_hat;
g31_temp = rz_hat.*rx_hat;

g11 = mean(g11_temp);
g22 = mean(g22_temp);
g33 = mean(g33_temp);
g12 = mean(g12_temp);
g23 = mean(g23_temp);
g31 = mean(g31_temp);

g11_err = std(g11_temp)/sqrt(length(g11_temp));
g22_err = std(g22_temp)/sqrt(length(g22_temp));
g33_err = std(g33_temp)/sqrt(length(g33_temp));
g12_err = std(g12_temp)/sqrt(length(g12_temp));
g23_err = std(g23_temp)/sqrt(length(g23_temp));
g31_err = std(g31_temp)/sqrt(length(g31_temp));

% g11_temp = rx_hat(bondtypes==2).*rx_hat(bondtypes==2);
% g22_temp = ry_hat(bondtypes==2).*ry_hat(bondtypes==2);
% g33_temp = rz_hat(bondtypes==2).*rz_hat(bondtypes==2);
% g12_temp = rx_hat(bondtypes==2).*ry_hat(bondtypes==2);
% g23_temp = ry_hat(bondtypes==2).*rz_hat(bondtypes==2);
% g31_temp = rz_hat(bondtypes==2).*rx_hat(bondtypes==2);
% 
% g11_st = mean(g11_temp);
% g22_st = mean(g22_temp);
% g33_st = mean(g33_temp);
% g12_st = mean(g12_temp);
% g23_st = mean(g23_temp);
% g31_st = mean(g31_temp);
% 
% g11_st_err = std(g11_temp)/sqrt(length(g11_temp));
% g22_st_err = std(g22_temp)/sqrt(length(g22_temp));
% g33_st_err = std(g33_temp)/sqrt(length(g33_temp));
% g12_st_err = std(g12_temp)/sqrt(length(g12_temp));
% g23_st_err = std(g23_temp)/sqrt(length(g23_temp));
% g31_st_err = std(g31_temp)/sqrt(length(g31_temp));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [msd,msd_err] = ComputeMSDs(timestep,AtomData,Corners,...
    MSDStartStep,type)

Lx = Corners(5)-Corners(2);
Ly = Corners(6)-Corners(3);
Lz = Corners(7)-Corners(4);

% Sort AtomData by ID number
AtomData = sortrows(AtomData,1);

% Check that these are aligned
AtomData0 = AtomData(AtomData(:,1)==MSDStartStep,:);
AtomDataC = AtomData(AtomData(:,1)==timestep,:);

% Sort AtomDatas by particle ID
AtomData0 = sortrows(AtomData0,2);
AtomDataC = sortrows(AtomDataC,2);

Pos0 = AtomData0(:,3:5);
Pos = AtomDataC(:,3:5);

% Eliminate partilces of wrong type
Pos0(AtomData0(:,end)~=type,:) = [];
Pos(AtomDataC(:,end)~=type,:) = [];

% N = size(Pos,1);
dr = Pos-Pos0;

Maxdr = max(vecnorm(dr,2,2));
% if Maxdr>8
%     disp(Maxdr);
% end

% Adjust for periodic bounds
if ~isempty(dr(abs(dr(:,1))>0.5*Lx,1)) || ~isempty(dr(abs(dr(:,2))>0.5*Ly,2)) ||...
        ~isempty(dr(abs(dr(:,3))>0.5*Lz,3))
    dr(abs(dr(:,1))>0.5*Lx,1);
end
dr(abs(dr(:,1))>0.5*Lx,1) = dr(abs(dr(:,1))>0.5*Lx,1)-Lx*sign(dr(abs(dr(:,1))>0.5*Lx,1));
dr(abs(dr(:,2))>0.5*Ly,2) = dr(abs(dr(:,2))>0.5*Ly,2)-Ly*sign(dr(abs(dr(:,2))>0.5*Ly,2));
dr(abs(dr(:,3))>0.5*Lz,3) = dr(abs(dr(:,3))>0.5*Lz,3)-Lz*sign(dr(abs(dr(:,3))>0.5*Lz,3));

% for type=1:2
%     dr_temp = dr(AtomData0(:,end)==type,:);
MSD_temp = mean((vecnorm(dr,2,2)).^2);
MSD_err = std((vecnorm(dr,2,2)).^2)/sqrt(size(dr,1));

msd = MSD_temp;
msd_err = MSD_err;
% end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kd,ka,dynamic_pairs0,N_bound0,N_free0] = ...
    ComputeBondDynamics(timestep,i,dt,BondData,AtomData,...
    dynamic_pairs0,N_bound0,N_free0,timestep0)

Dt = (timestep-timestep0)*dt;

% Isolate current data
BondsTemp = BondData(BondData(:,1)==timestep,:);
AtomsTemp = AtomData(AtomData(:,1)==timestep,:);

% Isolate dynamic bonds
pairs = BondsTemp(:,4:5);
bondtypes = BondsTemp(:,3);
dynamic_pairs = pairs(bondtypes==3,:);

% Delete redundant dyanmic pairs
dynamic_pairs = sort(dynamic_pairs,2);
dynamic_pairs = sortrows(dynamic_pairs,1);
dynamic_pairs(2:2:end,:) = [];

% compute number of bound and free stickers
N_stickers = length(AtomsTemp(AtomsTemp(:,end)==2));
N_bound = size(dynamic_pairs,1);

N_free = N_stickers/2-N_bound;

% ID previously existing and new dynamic bonds
if isempty(dynamic_pairs0)
%     old_pairs = [];
    new_pairs = dynamic_pairs;
else
    old_nodes1 = ismember(dynamic_pairs(:,1),dynamic_pairs0(:,1));
    old_nodes2 = ismember(dynamic_pairs(:,2),dynamic_pairs0(:,2));
    indx = old_nodes1 + old_nodes2; %new if this sums to 2
%     old_pairs = dynamic_pairs(indx==2,:);
    new_pairs = dynamic_pairs(indx~=2,:);
end

no_attachments = size(new_pairs,1)/2;   %divide by 2 because Sam's code doesn't delete redundant bonds

if i==1
    ka = 0;
else
    ka = no_attachments/N_free0/Dt;
end

% ID lost bonds (i.e., detachments)
if isempty(dynamic_pairs0)
%     surviving_pairs = [];
    deleted_pairs = [];
else
    surviving_nodes1 = ismember(dynamic_pairs0(:,1),dynamic_pairs(:,1));
    surviving_nodes2 = ismember(dynamic_pairs0(:,2),dynamic_pairs(:,2));
    indx = surviving_nodes1 + surviving_nodes2;
%     surviving_pairs = dynamic_pairs0(indx==2);
    deleted_pairs = dynamic_pairs0(indx~=2);
end

no_detachments =  size(deleted_pairs,1);

if i==1 || N_bound0==0
    kd = 0;
else
    kd = no_detachments/N_bound0/Dt;
end

% disp(['ka = ',num2str(ka),', kd = ',num2str(kd)])

% Set old values
dynamic_pairs0 = dynamic_pairs;
N_bound0 = N_bound;
N_free0 = N_free;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigma,conc] = ComputeStressAndBondConc(f,r,Corners)

% Compute volume
LWH = abs(Corners(2:4)-Corners(5:7));
V = prod(LWH);

% Dyad f \otimes r
dyad = [f(:,1).*r(:,1),...
    f(:,2).*r(:,2),...
    f(:,3).*r(:,3),...
    f(:,1).*r(:,2),...
    f(:,2).*r(:,3),...
    f(:,3).*r(:,1)];

N_bonds = size(f,1);
conc = N_bonds/V;

sigma = 1/(2*V)*sum(dyad,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function forces = ComputeForces(types,norms,r)

global ConstantsFileName

% FileName = ['Matlab Topological Data/Constants.Sample_1.Np_',num2str(Np),...
%     '.Nt_',num2str(N),'.Nb_',num2str(Nt),'.Ns_',num2str(Ns),'.txt'];
% FileName = ['Matlab Topological Data/Constants.S1',FileTagCompiled,...
%     '.txt'];

Table = readtable(ConstantsFileName);
constants = table2array(Table);
N_Kuhn = constants(4); b = constants(5); kbT = constants(6);

r_hat = [r(:,1)./norms r(:,2)./norms r(:,3)./norms];

force_mags = zeros(size(types));
Stiffness1 = 3*kbT/(N_Kuhn*b^2);
Stiffness3 = 3*kbT/(b^2)/2; % Divide by 2 because dynamic bonds are double counted in TNT bond package

force_mags(types==1) = Stiffness1*norms(types==1);
force_mags(types==2) = Stiffness1*norms(types==2);
force_mags(types==3) = Stiffness3*norms(types==3);

forces = force_mags.*r_hat;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotIntStressStrain(steps,time,Np,...
    sig11_ensemble, sig22_ensemble,sig33_ensemble)

global FileTag FontSize TurnOnDynamics AddOn

[sig11,sig11_err] = CalcMeanAndSE(sig11_ensemble,2);
[sig22,sig22_err] = CalcMeanAndSE(sig22_ensemble,2);
[sig33,sig33_err] = CalcMeanAndSE(sig33_ensemble,2);

figure(2); clf; hold on
TempStress = sig11; TempErr = sig11_err;
Color = 'k'; Style = '-'; TempStep = 1;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);
Color = 'b'; Style = '-'; TempStep = 2;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);
Color = 'c'; Style = '-'; TempStep = 3;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);

TempStress = sig22; TempErr = sig22_err;
Color = 'k'; Style = '--'; TempStep = 1;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);
Color = 'b'; TempStep = 2;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);
Color = 'c'; TempStep = 3;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);

TempStress = sig33; TempErr = sig33_err;
Color = 'k'; Style = ':'; TempStep = 1;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);
Color = 'b'; TempStep = 2;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);
Color = 'c'; TempStep = 3;
[~,~] = PlotCurve(time(steps==TempStep),...
    TempStress(steps==TempStep),...
    TempErr(steps==TempStep),Color,Style);

set(gca,'FontSize',FontSize/1.5)
xlabel('$t$ [s]','FontSize',FontSize,'Interpreter','latex')
ylabel('$\sigma$ [kPa]','FontSize',FontSize,'Interpreter','latex')
pbaspect([1 1 1])
set(gcf,'Color','w')

if ~isfolder('Output Plots')
    mkdir('Output Plots')
end
if TurnOnDynamics==0
    FileName = 'Control.Stress Check';
    title(['Dynamics OFF, $N_p$ = ',num2str(Np)],'FontSize',FontSize/1.5,'Interpreter','latex')
else
    FileName = 'Stress Check';
    title(['Dynamics ON, $N_p$ = ',num2str(Np)],'FontSize',FontSize/1.5,'Interpreter','latex')
end
CheckStressFileName = ['Output Plots/',FileName,FileTag,'.png'];
saveas(gcf,CheckStressFileName)
CheckStressFileName = ['Output Plots/',FileName,FileTag,'.fig'];
saveas(gcf,CheckStressFileName)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveCompiledFiles(sig11_ensemble,sig22_ensemble,sig33_ensemble,...
    time,stretch,steps,CompiledStress11FileName,CompiledStress22FileName,...
        CompiledStress33FileName)

% global CompiledStress11FileName CompiledStress22FileName...
%     CompiledStress33FileName

writematrix([sig11_ensemble,time,stretch,steps],CompiledStress11FileName)
writematrix([sig22_ensemble,time,stretch,steps],CompiledStress22FileName)
writematrix([sig33_ensemble,time,stretch,steps],CompiledStress33FileName)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveSampleFiles(Corners,Atoms,Bonds,...
    CornersFileName,AtomsFileName,BondsFileName)

% global AtomsFileName BondsFileName CornersFileName

writematrix(Corners,CornersFileName);
writematrix(Atoms,AtomsFileName);
writematrix(Bonds,BondsFileName);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Bonds = UnpackBondsDump(timestep)

global OutputBond

N_output = length(timestep);

% Extract periodic domain boundaries
fmt1=repmat('%f',1,6);             % first data section format
fid=fopen(['Outputs/',OutputBond],'r');     % open file
l=fgetl(fid);                      % read first line
ct = 0;
wb = waitbar(0,'Extracting bonds data...');
while ~feof(fid)                   % loop through file by record
    if strfind(l,'ITEM: ENTRIES index c_Pairs[1] c_Pairs[2] c_Pairs[3] c_Bonds[1] c_Bonds[2] c_Bonds[3] c_Bonds[4] c_Bonds[5]')
        ct = ct+1;
        Bonds_array{:,:,ct} = cell2mat(textscan(fid,'%f'));
        if ~mod(ct,50)
            PcntComplete = ct/N_output;
            waitbar(PcntComplete,wb,'Extracting bonds data...')
        end
    end                                 % break when find first section
    l=fgetl(fid);                       % next record
end
close(wb)

% Unpack bond information
N_var = 9;      % id bond_type node1 node2 rx ry rz norm force
Bonds = zeros(1,N_var+1);
ct = 0;
wb = waitbar(0,'Unpacking bond data...');
StartIndx = 0;
for i=1:N_output
    if ~mod(i,50)
        PcntComplete = i/N_output;
        waitbar(PcntComplete,wb,'Unpacking bond data...')
    end
    timestep_temp = timestep(i);
    Bonds_temp = Bonds_array{:,:,i};
    N_bonds = length(Bonds_temp)/N_var;
    Range = StartIndx+(1:N_bonds)';

    ID_indx = (1:N_var:size(Bonds_temp,1))';
    type_indx = ID_indx+1;
    node1_indx = type_indx+1;
    node2_indx = node1_indx+1;
    rx_indx = node2_indx+1;
    ry_indx = rx_indx+1;
    rz_indx = ry_indx+1;
    norm_indx = rz_indx+1;
    force_indx = norm_indx+1;

    Bonds(Range,1) = timestep_temp*ones(size(Range));
    Bonds(Range,2) = Bonds_temp(ID_indx);
    Bonds(Range,3) = Bonds_temp(type_indx);
    Bonds(Range,4) = Bonds_temp(node1_indx);
    Bonds(Range,5) = Bonds_temp(node2_indx);
    Bonds(Range,6) = Bonds_temp(rx_indx);
    Bonds(Range,7) = Bonds_temp(ry_indx);
    Bonds(Range,8) = Bonds_temp(rz_indx);
    Bonds(Range,9) = Bonds_temp(norm_indx);
    Bonds(Range,10) = Bonds_temp(force_indx);

    StartIndx = size(Bonds,1);
end
close(wb)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [timestep,Corners,Atoms] = UnpackAtomsDump

global OutputAtom

% Extract timesteps
fmt1=repmat('%f',1,6);             % first data section format
fid=fopen(['Outputs/',OutputAtom],'r');     % open file
l=fgetl(fid);                      % read first line
ct = 0;
while ~feof(fid)                   % loop through file by record
    if strfind(l,'ITEM: TIMESTEP')
        ct = ct+1;
        timestep(ct,1) = cell2mat(textscan(fid,'%f'));
    end                                 % break when find first section
    l=fgetl(fid);                       % next record
end
N_output = length(timestep);

% Extract periodic domain boundaries
fmt1=repmat('%f',1,6);             % first data section format
fid=fopen(['Outputs/',OutputAtom],'r');     % open file
l=fgetl(fid);                      % read first line
ct = 0;
wb = waitbar(0,'Extracting domain boundaries...');
while ~feof(fid)                   % loop through file by record
    if strfind(l,'ITEM: BOX BOUNDS pp pp pp')
        ct = ct+1;
        Corners_array{:,:,ct} = cell2mat(textscan(fid,'%f'));
        if ~mod(ct,50)
            PcntComplete = ct/N_output;
            waitbar(PcntComplete,wb,'Extracting domain boundaries...')
        end
    end                                 % break when find first section
    l=fgetl(fid);                       % next record
end
close(wb)

% Unpack corners data
N_var = 6;      % x1 x2 y1 y2 z1 z2
Corners = zeros(N_output,N_var+1);  %+1 for timesteps
Corners(:,1) = timestep;
wb = waitbar(0,'Unpacking boundaries...');
for i=1:N_output
    if ~mod(i,50)
        PcntComplete = i/N_output;
        waitbar(PcntComplete,wb,'Unpacking boundaries...')
    end
    Corners_temp = Corners_array{:,:,i};
    Corners(i,2) = Corners_temp(1);
    Corners(i,3) = Corners_temp(3);
    Corners(i,4) = Corners_temp(5);
    Corners(i,5) = Corners_temp(2);
    Corners(i,6) = Corners_temp(4);
    Corners(i,7) = Corners_temp(6);
end
close(wb)

% Extract atom index, position, velocity, and type
% information
fmt1=repmat('%f',1,6);             % first data section format
fid=fopen(['Outputs/',OutputAtom],'r');     % open file
l=fgetl(fid);                      % read first line
ct = 0;
wb = waitbar(0,'Extracting atom data...');
while ~feof(fid)                   % loop through file by record
    if strfind(l,'ITEM: ATOMS id x y z vx vy vz type')
        ct = ct+1;
        if ~mod(ct,50)
            PcntComplete = ct/N_output;
            waitbar(PcntComplete,wb,'Extracting atom data...')
        end
        Atoms_array{:,:,ct} = cell2mat(textscan(fid,'%f'));
    end                                 % break when find first section
    l=fgetl(fid);                       % next record
end
close(wb)

% Unpack atom information
N_var = 8;      % id x y z vx vy vz type
Atoms = zeros(1,N_var+1);
ct = 0;
wb = waitbar(0,'Unpacking atom data...');
StartIndx = 0;
for i=1:N_output
    if ~mod(i,50)
        PcntComplete = i/N_output;
        waitbar(PcntComplete,wb,'Unpacking atom data...')
    end
    timestep_temp = timestep(i);
    Atoms_temp = Atoms_array{:,:,i};
    N_atoms = length(Atoms_temp)/N_var;
    Range = StartIndx+(1:N_atoms)';
    
    ID_indx = (1:N_var:size(Atoms_temp,1))';
    x_indx = ID_indx+1;
    y_indx = x_indx+1;
    z_indx = y_indx+1;
    vx_indx = z_indx+1;
    vy_indx = vx_indx+1;
    vz_indx = vy_indx+1;
    type_indx = vz_indx+1;
   
    Atoms(Range,1) = timestep_temp*ones(size(Range));
    Atoms(Range,2) = Atoms_temp(ID_indx);
    Atoms(Range,3) = Atoms_temp(x_indx);
    Atoms(Range,4) = Atoms_temp(y_indx);
    Atoms(Range,5) = Atoms_temp(z_indx);
    Atoms(Range,6) = Atoms_temp(vx_indx);
    Atoms(Range,7) = Atoms_temp(vy_indx);
    Atoms(Range,8) = Atoms_temp(vz_indx);
    Atoms(Range,9) = Atoms_temp(type_indx);

    StartIndx = size(Atoms,1);
end
close(wb)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time,stretch,steps,sig11,sig22,sig33]...
    = CompileStress(EqSize,LdSize,RlSize)

global StressEqlFileName StressLdgFileName StressRlxFileName

% Import data
% During equilibration
EquilibrationData = dlmread(['Outputs/',StressEqlFileName],'',1,0);
if size(EquilibrationData,1)~=EqSize
    Diff = size(EquilibrationData,1)-EqSize;
    if Diff>0
        EquilibrationData(end-Diff+1:end,:) = [];
    else
        EquilibrationData(end:end+Diff,:) = NaN;
    end
end
t_eq = EquilibrationData(:,1);
lambda_eq = ones(size(t_eq));
s11_eq = EquilibrationData(:,2)/1e3;
s22_eq = EquilibrationData(:,3)/1e3;
s33_eq = EquilibrationData(:,4)/1e3;
Step_eq = ones(size(t_eq));

% During deformation
LoadingData = dlmread(['Outputs/',StressLdgFileName],'',1,0);
if size(LoadingData,1)~=LdSize
    Diff = size(LoadingData,1)-LdSize;
    if Diff>0
        LoadingData(end-Diff+1:end,:) = [];
    else
        LoadingData(end:end+Diff,:) = NaN;
    end
end
t_ld = LoadingData(:,1) + t_eq(end);
lambda_ld = LoadingData(:,2);
s11_ld = LoadingData(:,3)/1e3;
s22_ld = LoadingData(:,4)/1e3;
s33_ld = LoadingData(:,5)/1e3;
Step_ld = 2*ones(size(t_ld));

% During relaxation
RelaxationData = dlmread(['Outputs/',StressRlxFileName],'',1,0);
if size(RelaxationData,1)~=RlSize
    Diff = size(RelaxationData,1)-RlSize;
    if Diff>0
        RelaxationData(end-Diff+1:end,:) = [];
    else
        RelaxationData(end:end+Diff,:) = NaN;
    end
end
t_rl = RelaxationData(:,1) + t_ld(end);
lambda_rl = RelaxationData(:,2);
s11_rl = RelaxationData(:,3)/1e3;
s22_rl = RelaxationData(:,4)/1e3;
s33_rl = RelaxationData(:,5)/1e3;
Step_rl = 3*ones(size(t_rl));

time = [t_eq;t_ld;t_rl];
stretch = [lambda_eq;lambda_ld;lambda_rl];
steps = [Step_eq;Step_ld;Step_rl];
sig11 = [s11_eq;s11_ld;s11_rl];
sig22 = [s22_eq;s22_ld;s22_rl];
sig33 = [s33_eq;s33_ld;s33_rl];

Check = 0;
if Check==1
    figure(1); clf; hold on
    p = plot(t_eq,s11_eq,'k');
    p.LineWidth = 1.5;
    p = plot(t_ld,s11_ld,'b');
    p.LineWidth = 1.5;
    p = plot(t_rl,s11_rl,'c');
    p.LineWidth = 1.5;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Corners,N_steps] = ImportDomainBoundaries(EdgeFactor,Np,N,...
    ka,kd,f0,dt,damp,N_Kuhn,b,Separation)

global FileTagCompiled 

% Initializes non-sweeping parameters
InputScript(EdgeFactor,1,Np,N,ka,kd,f0,dt,damp,N_Kuhn,b,Separation);
SetDirAndFileNames;

FileName = ['Data/Compiled Outputs/Corners',FileTagCompiled,'.txt'];
Corners = readmatrix(FileName);
N_steps = size(Corners,1);

end
