function PostProcessData(Package,OverridePostProcess,ToggleDynamics,...
    MakeHistos,MakeMoviesHistos,OverrideHistogramMovies,...
    OverrideMSDFitData,WRTStretchOrLength,HistoUnits,LC,DC,BSOM,CM,...
    OverrideCompareData,CF,OF)

% Computes and plots stress-strain for every iteration of simulation

global LineWidth TurnOnDynamics FontSize...
    LengthConversion DamperConversion DataSize BeadSpringOrMeso...
    CompareModels CurrentFolder OutputFolder PlotForFigures

LineWidth = 1.5;
FontSize = 20;
DataSize = 30;
TurnOnDynamics = ToggleDynamics;
LengthConversion = LC;
DamperConversion = DC;
BeadSpringOrMeso = BSOM;
CompareModels = CM;
CurrentFolder = CF;
OutputFolder = OF;
PlotForFigures = 0; %1 for version of code from which figs were plotted

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

%% Plot the histograms
if MakeHistos  %Plots histograms for (meso and bead-spring) vs Gaussian, independently
    PlotAllHistograms(MakeMoviesHistos,OverrideHistogramMovies,...
        WRTStretchOrLength,HistoUnits,Package,...
        N_Kuhns,Ns,T,Nb,ka_in,kd_in,f0,b,N,Separations)
end
if MakeHistos && CompareModels
    PlotAllHistogramsCompare(OverrideCompareData,MakeMoviesHistos,...
        HistoUnits,Package,damps,N_Kuhns,Np,Ns,T,Nb,ka_in,kd_in,f0,dt,b,N,...
        Separations)
end

%% Fit the Rouse Model to the Bead-spring control system
% OverrideMSDFitData = 1;
DetermineRouseFits(OverrideMSDFitData,Package,damps,N_Kuhns,Np,Ns,T,Nb,...
    ka_in,kd_in,f0,dt,b,N,Separations);

%% Measure the steady state MSDs for each model
MeasureSSMSDs(Package,damps,N_Kuhns,Np,Ns,T,Nb,ka_in,kd_in,f0,dt,b,...
    N,Separations,OverridePostProcess)

%% Sweep over all dampers. For each damper, sweep over N and plot MSD data
OverridePostProcess = 1;
GenerateMSDPlotsWithNSwept(OverridePostProcess,Package,damps,N_Kuhns,Np,Ns,...
    T,Nb,ka_in,kd_in,f0,dt,b,N,Separations);

%% Compare MSD behavior accross models
if CompareModels
    OverrideCompareData = 1;
    GenerateMSDPlotsWithNSweptCompare(OverrideCompareData,...
        Package,damps,N_Kuhns,Np,Ns,T,Nb,ka_in,kd_in,f0,dt,b,N,Separations);
end

% OverrideMSDFitData = 0;
%% Sweep over all N. For each N, sweep over damper
GeneratePlotsWithDamperSwept(OverridePostProcess,Package,...
    damps,N_Kuhns,Np,Ns,T,Nb,ka_in,kd_in,f0,dt,b,N,Separations);

%% Sweep overa all N and plot fitted parameters wrt damper
if length(damps)>1
    PlotDiffusionModelParametersWithNSwept(N_Kuhns,damps,OverridePostProcess,b)

    % %% Sweep over all dampers and plot fitted parameters wrt N
    % % THIS PROVES UNNECESSARY AS NO SIGNIFICANT VARIATION IS OBSERVED AS
    % % FUNCITON OF GAMMA - THIS REDUCES THE PARAMETRIC DESIGN SPACE.
    % % PlotDiffusionModelParametersWithDSwept(N_Kuhns,damps,tau0,taur,R2_beadvsrouse,...
    % %     tau0_2,taur_2,alpha,R2_2,OverrideMSDFitData,b)

    %% Make surface plots of fitted parameters wrt dampers and N
    MakeSurfacePlots(damps,N_Kuhns,b)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DetermineRouseFits(OverrideMSDFitData,Package,...
    damps,N_Kuhns,Np,Ns,T,Nb,ka_in,kd_in,f0,dts,b,N,Separations)

global tau0_filename taur_filename R2_beadvsrouse_filename...
    FontSize StartMSDMeasurePct...
    AddOn BeadSpringOrMeso LengthConversion DamperConversion

StartMSDMeasurePct = 0.001; % Leave close to zero - this is only useful if 
                            % chains were initiated anomalously far from 
                            % equilibrium and they weren't as evidenced by
                            % end-to-end distributions and these results'
                            % independence of of StartMSDMeasurePcnt

TotalSims = size(Package,1);

taus = zeros(length(damps),length(N_Kuhns));
R2_beadvsrouse = zeros(length(damps),length(N_Kuhns));
tauR = zeros(length(damps),length(N_Kuhns));
SimCt = 0;

wb1 = waitbar(0,'Post processing for const. damper, and swept N...');

EdgeFactor = 1;
InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dts(1),damps(1),...
    N_Kuhns(1),b,Separations(1));
DefineMSDFileNames;

if ~isfile(tau0_filename) || ~isfile(taur_filename) ||...
        ~isfile(R2_beadvsrouse_filename) || OverrideMSDFitData
    % MSD wrt time for each damper while sweeping N
    for dm = 1:length(damps)
        damp = damps(dm);
        Folder = ['Output Plots/damp ',num2str(damp,'%.3e')];
        if ~isfolder(Folder)
            mkdir(Folder)
        end

        % Figures
        figure(1); clf; hold on % msd in tau0<t<taur with fitted Rouse

        ColorRange = (linspace(0,0.85,length(N_Kuhns)))';

        Colors = [zeros(size(ColorRange)) flipud(ColorRange) zeros(size(ColorRange))]; %Green-to-black

        LegendEntries = [];
        p0 = [];
        for nk = 1:length(N_Kuhns)
            N_Kuhn = N_Kuhns(nk);
            Color = Colors(nk,:);

            LegendEntries{nk} = ['$N =$ ',num2str(N_Kuhn)];

            for sp = 1:length(Separations)
                Separation = Separations(sp);
                dt = dts(dm);

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

                    BSOM_orig = BeadSpringOrMeso;
                    BSOM = 0; Np_temp = 5^3;
                    [dt,damp,~,~] = DefineTimeStep(b,LengthConversion,...
                        DamperConversion,BSOM);

                    [time,msd_bead_mean,msd_err,~,~,~,~,~,~,~,~] = ...
                        ImportMSD(BSOM,Np_temp,N,ka_in,kd_in,f0,dt,damp,...
                        N_Kuhn,b,Separation);
                    BeadSpringOrMeso = BSOM_orig;    %Reset global variable

                    time = time-time(1);    %Set reference time

                    [R2_beadvsrouse(dm,nk),taus(dm,nk),tauR(dm,nk),rng] =...
                        FitRouseModel(time,msd_bead_mean,N_Kuhn);

                    %% Fit the Rouse Model
                    figure(1)
                    Style = '-';
                    [p0(nk),~] = PlotCurve(time(rng)/tauR(dm,nk),...
                        msd_bead_mean(rng)/b^2,...
                        msd_err(rng)/b^2,Color,Style);

                    x = (linspace(0,tauR(dm,nk),100))';
                    y = (b^2)*(x/taus(dm,nk)).^0.5/b^2;
                    p = plot(x/tauR(dm,nk),y);
                    p.Color = Color;
                    p.LineStyle = '--';
                    p.LineWidth = 1.5;
                end
            end
        end
        figure(1)
        set(gca,'FontSize',FontSize/1.5)
        set(gcf,'Color','w')
        pbaspect([1 1 1])
        xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
        ylabel('MSD ($b^2$)','FontSize',FontSize,'Interpreter','latex')
%         title(['$\gamma$ = ',num2str(damp,'%.2e'),' [force $\cdot$ time/length], ',...
%             num2str(damp*DamperConversion,'%.3e'),' [N s m$^{-1}$]'],...
%             'FontSize',FontSize/2,'Interpreter','latex')

        l = legend(p0,LegendEntries);
        l.FontSize = FontSize/2;
        l.Interpreter = 'latex';
        l.Location = 'Northwest';

        FileName = [Folder,'/',AddOn,'Fitted Rouse Model.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'Fitted Rouse Model.fig'];
        saveas(gcf,FileName)
    end
    writematrix(taus,tau0_filename)
    writematrix(tauR,taur_filename)
    writematrix(R2_beadvsrouse,R2_beadvsrouse_filename)
end
close(wb1)
close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeSurfacePlots(damps,N_Kuhns,b)

global DamperConversion FontSize AddOn BeadSpringOrMeso...
    MSD_ss_filename tau0_filename taur_filename R2_beadvsrouse_filename...
    R2_beadvsmeso_filename R2_mesovsrouse_filename...
    R2_msdr_filename R2_msdt_filename...
    MSD_r00_filename MSD_t00_filename %MSD_r0_filename MSD_t0_filename 

xlab = '$\gamma$ (kg s$^{-1}$)';
ylab = '$N$';

%% tau0
tau0 = readmatrix(tau0_filename);
figure(1)
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

xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
ylabel(ylab,'FontSize',FontSize,'Interpreter','latex')
% title(['$b$ = ',num2str(b*LengthConversion*1e9,'%.2f'),' nm'],...
%     'FontSize',FontSize/1.5','Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = 'Output Plots/tau0';
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

%% taur
tauR = readmatrix(taur_filename);
figure(2)
clf; hold on
s = surf(damps*DamperConversion,N_Kuhns,tauR');
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

xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
ylabel(ylab,'FontSize',FontSize,'Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = 'Output Plots/taur';
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

%% R2 bead vs rouse
R2 = readmatrix(R2_beadvsrouse_filename);
figure(3)
clf; hold on
s = surf(damps*DamperConversion,N_Kuhns,R2');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'FontSize',FontSize/1.5)

xlim([min(damps) max(damps)]*DamperConversion)
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$R^2$ (Bead-spring vs. Rouse)');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',FontSize);
caxis([0 1])

xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
ylabel(ylab,'FontSize',FontSize,'Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = 'Output Plots/R2 Bead-spring vs. Rouse';
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

%% R2 bead vs meso
R2 = readmatrix(R2_beadvsmeso_filename);
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
set(get(c,'label'),'string','$R^2$ (Bead-spring vs. Mesoscale)');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',FontSize);
caxis([0 1])

xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
ylabel(ylab,'FontSize',FontSize,'Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = 'Output Plots/R2 Bead-spring vs. mesoscale';
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

%% R2 meso vs Rouse
R2 = readmatrix(R2_mesovsrouse_filename);
figure(5)
clf; hold on
s = surf(damps*DamperConversion,N_Kuhns,R2');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'FontSize',FontSize/1.5)

xlim([min(damps) max(damps)]*DamperConversion)
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$R^2$ (Mesoscale vs. Rouse)');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',FontSize);
caxis([0 1])

xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
ylabel(ylab,'FontSize',FontSize,'Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = 'Output Plots/R2 mesoscale vs. Rouse';
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

%% R2 bead-spring vs meso radial
R2 = readmatrix(R2_msdr_filename);
figure(6)
clf; hold on
s = surf(damps*DamperConversion,N_Kuhns,R2');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'FontSize',FontSize/1.5)

xlim([min(damps) max(damps)]*DamperConversion)
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$R^2$ (Bead-spring vs. Mesoscale - radial MSD)');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',FontSize);
caxis([0 1])

xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
ylabel(ylab,'FontSize',FontSize,'Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = 'Output Plots/R2 mesoscale vs. Rouse - radial';
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

%% R2 bead-spring vs meso radial
R2 = readmatrix(R2_msdt_filename);
figure(7)
clf; hold on
s = surf(damps*DamperConversion,N_Kuhns,R2');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'FontSize',FontSize/1.5)

xlim([min(damps) max(damps)]*DamperConversion)
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$R^2$ (Bead-spring vs. Mesoscale - tang. MSD)');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',FontSize);
caxis([0 1])

xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
ylabel(ylab,'FontSize',FontSize,'Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = 'Output Plots/R2 mesoscale vs. Rouse - tangential';
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

for BeadSpringOrMeso=0:1
    EdgeFactor = 1;
    InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dts(1),damps(1),...
        N_Kuhns(1),b,Separations(1));
    DefineMSDFileNames;

    %% MSD_ss
    MSD_ss = readmatrix(MSD_ss_filename);
    figure(8)
    clf; hold on
    s = surf(damps*DamperConversion,N_Kuhns,MSD_ss'/b^2);
    s.FaceColor = 'interp';
    s.EdgeAlpha = 0.5;
    set(gca,'XScale','log')
    set(gca,'FontSize',FontSize/1.5)

    xlim([min(damps) max(damps)]*DamperConversion)
    ylim([min(N_Kuhns) max(N_Kuhns)])

    c = colorbar;
    set(get(c,'label'),'string','MSD$^{ss}$ (b$^{2}$)');
    set(get(c,'label'),'Interpreter','latex');
    set(get(c,'label'),'FontSize',FontSize);
    caxis([25 85])

    xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
    ylabel(ylab,'FontSize',FontSize,'Interpreter','latex')
    pbaspect([1 1 1])
    set(gcf,'Color','w')

    FileName = 'Output Plots/Steady State MSD';
    saveas(gcf,[FileName,AddOn,'.png'])
    saveas(gcf,[FileName,AddOn,'.fig'])

%%%figures(9-10) (MSDs from times t->t+dt) are better expressed as 2D plots 
%%%and presented as such in the manuscript
%     %% MSD radial
%     figure(9)
%     MSD_r0 = readmatrix(MSD_r0_filename);
%     clf; hold on
%     s = surf(damps*DamperConversion,N_Kuhns,MSD_r0'/b^2);
%     s.FaceColor = 'interp';
%     s.EdgeAlpha = 0.5;
%     set(gca,'XScale','log')
%     set(gca,'FontSize',FontSize/1.5)
% 
%     xlim([min(damps) max(damps)]*DamperConversion)
%     ylim([min(N_Kuhns) max(N_Kuhns)])
% 
%     c = colorbar;
%     set(get(c,'label'),'string','MSD$_r^{ss}$ (b$^{2}$)');
%     set(get(c,'label'),'Interpreter','latex');
%     set(get(c,'label'),'FontSize',FontSize);
%     caxis([-0 1])
% 
%     xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
%     ylabel(ylab,'FontSize',FontSize,'Interpreter','latex')
%     pbaspect([1 1 1])
% 
%     set(gcf,'Color','w')
% 
%     FileName = 'Output Plots/Inst. Radial MSD';
%     saveas(gcf,[FileName,AddOn,'.png'])
%     saveas(gcf,[FileName,AddOn,'.fig'])
% 
%     %% MSD tangential
%     figure(10)
%     MSD_t0 = readmatrix(MSD_t0_filename);
%     clf; hold on
%     s = surf(damps*DamperConversion,N_Kuhns,MSD_t0'/b^2);
%     s.FaceColor = 'interp';
%     s.EdgeAlpha = 0.5;
%     set(gca,'XScale','log')
%     set(gca,'FontSize',FontSize/1.5)
% 
%     xlim([min(damps) max(damps)]*DamperConversion)
%     ylim([min(N_Kuhns) max(N_Kuhns)])
% 
%     c = colorbar;
%     set(get(c,'label'),'string','MSD$_t^{0}$ (b$^{2}$)');
%     set(get(c,'label'),'Interpreter','latex');
%     set(get(c,'label'),'FontSize',FontSize);
%     caxis([0 2])
% 
%     xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
%     ylabel(ylab,'FontSize',FontSize,'Interpreter','latex')
% 
%     pbaspect([1 1 1])
% 
%     set(gcf,'Color','w')
% 
%     FileName = 'Output Plots/Inst. Tang. MSD';
%     saveas(gcf,[FileName,AddOn,'.png'])
%     saveas(gcf,[FileName,AddOn,'.fig'])
    
    %% MSD radial
    figure(11)
    MSD_r00 = readmatrix(MSD_r00_filename);
    clf; hold on
    s = surf(damps*DamperConversion,N_Kuhns,MSD_r00'/b^2);
    s.FaceColor = 'interp';
    s.EdgeAlpha = 0.5;
    set(gca,'XScale','log')
    set(gca,'FontSize',FontSize/1.5)

    xlim([min(damps) max(damps)]*DamperConversion)
    ylim([min(N_Kuhns) max(N_Kuhns)])

    c = colorbar;
    set(get(c,'label'),'string','MSD$_r^{ss}$ (b$^{2}$)');
    set(get(c,'label'),'Interpreter','latex');
    set(get(c,'label'),'FontSize',FontSize);
    caxis([4 17])

    xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
    ylabel(ylab,'FontSize',FontSize,'Interpreter','latex')
    pbaspect([1 1 1])

    set(gcf,'Color','w')

    FileName = 'Output Plots/SS Radial MSD';
    saveas(gcf,[FileName,AddOn,'.png'])
    saveas(gcf,[FileName,AddOn,'.fig'])

    %% MSD tangential
    figure(12)
    MSD_t00 = readmatrix(MSD_t00_filename);
    clf; hold on
    s = surf(damps*DamperConversion,N_Kuhns,MSD_t00'/b^2);
    s.FaceColor = 'interp';
    s.EdgeAlpha = 0.5;
    set(gca,'XScale','log')
    set(gca,'FontSize',FontSize/1.5)

    xlim([min(damps) max(damps)]*DamperConversion)
    ylim([min(N_Kuhns) max(N_Kuhns)])

    c = colorbar;
    set(get(c,'label'),'string','MSD$_t^{ss}$ (b$^{2}$)');
    set(get(c,'label'),'Interpreter','latex');
    set(get(c,'label'),'FontSize',FontSize);
    caxis([35 115])

    xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
    ylabel(ylab,'FontSize',FontSize,'Interpreter','latex')

    pbaspect([1 1 1])

    set(gcf,'Color','w')

    FileName = 'Output Plots/SS Tang. MSD';
    saveas(gcf,[FileName,AddOn,'.png'])
    saveas(gcf,[FileName,AddOn,'.fig'])
end

close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DefineMSDFileNames

global MSD_ss_filename tau0_filename taur_filename R2_beadvsrouse_filename...
    BeadSpringOrMeso AddOn R2_beadvsmeso_filename R2_mesovsrouse_filename...
    MSD_r0_filename MSD_t0_filename R2_msdr_filename R2_msdt_filename...
    OutputFolderName MSD_r00_filename MSD_t00_filename

if BeadSpringOrMeso==0
    AddOn = 'Bead.';
else
    AddOn = 'Meso.';
end

MSD_ss_filename = [OutputFolderName,'Compiled Outputs/',AddOn,'Steady state msd.txt'];
MSD_r0_filename = [OutputFolderName,'Compiled Outputs/',AddOn,'Mean Instant radial MSD.txt'];
MSD_t0_filename = [OutputFolderName,'Compiled Outputs/',AddOn,'Mean Instant tangent MSD.txt'];
MSD_r00_filename = [OutputFolderName,'Compiled Outputs/',AddOn,'Mean SS radial MSD.txt'];
MSD_t00_filename = [OutputFolderName,'Compiled Outputs/',AddOn,'Mean SS tangent MSD.txt'];

tau0_filename = [OutputFolderName,'Compiled Outputs/Rouse tau0.txt'];
taur_filename = [OutputFolderName,'Compiled Outputs/Rouse taur.txt'];

R2_beadvsrouse_filename = [OutputFolderName,'/Compiled Outputs/R2 bead vs rouse.txt'];
R2_beadvsmeso_filename = [OutputFolderName,'Compiled Outputs/R2 bead vs meso.txt'];
R2_mesovsrouse_filename = [OutputFolderName,'Compiled Outputs/R2 meso vs rouse.txt'];

R2_msdr_filename = [OutputFolderName,'Compiled Outputs/',AddOn,'R2 msdr meso vs rouse.txt'];
R2_msdt_filename = [OutputFolderName,'Compiled Outputs/',AddOn,'R2 msdt meso vs rouse.txt'];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotDiffusionModelParametersWithNSwept(N_Kuhns,damps,...
    OverridePostProcess,b)

global FontSize DamperConversion DataSize AddOn tau0_filename...
    taur_filename R2_beadvsrouse_filename LengthConversion

R2_beadvsrouse = readmatrix(R2_beadvsrouse_filename);
tau0 = readmatrix(tau0_filename);
tauR = readmatrix(taur_filename);

kbT = (1.38e-23)*293; %J
b_m = b*LengthConversion;
ExpectedSlope = (b_m)^2/kbT;
disp(['Expected slope, b^2/kbT = ',num2str(ExpectedSlope,'%.2f')])

for nk=1:length(N_Kuhns)
    N_Kuhn = N_Kuhns(nk);
    tau0_temp = tau0(:,nk);
    taur_temp = tauR(:,nk);
    R2_temp = R2_beadvsrouse(:,nk);

    Folder = ['Output Plots/N_Kuhn ',num2str(N_Kuhn)];
    if ~isfolder(Folder)
        mkdir(Folder)
    end
    FileName = [Folder,'/',AddOn,'R2 vs damper.fig'];

    tau0_eff = zeros(size(tau0_temp));
    tauR_eff = zeros(size(tau0_temp));
    for d=1:length(damps)
        damp = damps(d);
        [tau0_eff(d),tauR_eff(d),~] = ComputeEffectiveDamper(N_Kuhn,b,damp);
    end

    if ~isfile(FileName) || OverridePostProcess
        figure(1); clf; hold on
        s1 = scatter(damps*DamperConversion,taur_temp,'k','filled');
        s1.SizeData = DataSize;
        s3 = scatter(damps*DamperConversion,tauR_eff,'kx');
        s3.SizeData = DataSize;
        set(gcf,'Color','w')
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        pbaspect([1 1 1])
        xlabel('$\gamma$ (kg s$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
        ylabel('$\tau_r$ (s)','FontSize',FontSize,'Interpreter','latex')
        title(['$N =$ ',num2str(N_Kuhn)],'FontSize',FontSize,'Interpreter','latex')

        l = legend([s1 s3],...
            'Fitted Rouse Model','Set \textit{a priori} (if untethered)');
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
        s3 = scatter(damps*DamperConversion,tau0_eff,'kx');
        s3.SizeData = DataSize;
        set(gcf,'Color','w')
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        pbaspect([1 1 1])
        xlabel('$\gamma$ (kg s$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
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


        l = legend([s1 s3],...
            'Fitted Rouse Model','Set \textit{a priori} (if untethered)');
        l.Location = 'Southeast';
        l.Interpreter = 'latex';
        
        FileName = [Folder,'/',AddOn,'diffusion time vs damper.png'];
        saveas(gcf,FileName)
        FileName = [Folder,'/',AddOn,'diffusion time vs damper.fig'];
        saveas(gcf,FileName)

        figure(3); clf; hold on
        s = scatter(damps*DamperConversion,R2_temp,'k','filled');
        s.SizeData = DataSize;
        plot(damps*DamperConversion,ones(size(damps)),'k--')
        set(gcf,'Color','w')
        set(gca,'XScale','log')
        pbaspect([1 1 1])
        ylim([0 1])
        xlabel('$\gamma$ (kg s$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
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
function GeneratePlotsWithDamperSwept(OverridePostProcess,Package,...
    damps,N_Kuhns,Np_orig,Ns,T,Nb,ka_in,kd_in,f0,dts,b,N,Separations)

global ConstantsFileName FontSize BeadSpringOrMeso AddOn...
    tau0_filename taur_filename LengthConversion DamperConversion

TotalSims = size(Package,1);
tau0 = readmatrix(tau0_filename);
tauR = readmatrix(taur_filename);

for BeadSpringOrMeso=0:1
    if BeadSpringOrMeso==0
        Np=5^3;
    else
        Np=11^3;
    end

    wb1 = waitbar(0,'Post processing for const. N, and swept dampers...');
    % MSD wrt time for each damper while sweeping N
    SimCt = 0;
    for nk = 1:length(N_Kuhns)
        N_Kuhn = N_Kuhns(nk);

        % Figures
        figure(1); clf; hold on % msd
        figure(2); clf; hold on % msd in tau0<t<taur

        ColorRange = (linspace(0,1,length(damps)))';
        %     Colors = [ColorRange flipud(ColorRange) flipud(ColorRange)];
        Colors = [zeros(size(ColorRange)),...
            flipud(ColorRange),...
            zeros(size(ColorRange))];
        p0 = [];


        Folder = ['Output Plots/N_Kuhn ',num2str(N_Kuhn)];
        if ~isfolder(Folder)
            mkdir(Folder)
        end
        FileName = [Folder,'/',AddOn,'MSD vs t.png'];

        if ~isfile(FileName) || OverridePostProcess==1
            for dm = 1:length(damps)
                damp = damps(dm);
                dt = dts(dm);

                Color = Colors(dm,:);
                LegendEntries{dm} = ['$\gamma$ = ',num2str(damp,'%.2e')];

                tau0_temp = tau0(dm,nk);
                taur_temp = tauR(dm,nk);

                for sp = 1:length(Separations)
                    Separation = Separations(sp);

                    PackageTemp = Package(ismember(Package(:,2),Np_orig),:);
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

                    if ~isempty(PackageTemp)
                        SimCt = SimCt+1;

                        waitbar(SimCt*size(PackageTemp,1)/TotalSims,wb1,...
                            'Post processing for const. N, and swept dampers...')

                        % Determine allocation size
                        EdgeFactor = 1;

                        %% Callout Input Script
                        % Initializes non-sweeping parameters
                        InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dt,damp,...
                            N_Kuhn,b,Separation);
                        SetDirAndFileNames;
                        DefineCompiledFileNames;
%                         Table = readtable(ConstantsFileName);
%                         constants = table2array(Table);
%                         N_Kuhn = constants(3); b = constants(4);

                        BSOM = BeadSpringOrMeso;
                        [dt,damp,~,~] = DefineTimeStep(b,LengthConversion,...
                            DamperConversion,BSOM);

                        [time,msd_st_mean,msd_st_err,~,~,~,~,~,~,~,~] = ...
                            ImportMSD(BeadSpringOrMeso,Np,N,ka_in,kd_in,f0,dt,...
                            damp,N_Kuhn,b,Separation);

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
                        [~,~] = PlotCurve(time(rng)/taur_temp,msd_st_mean(rng)/b^2,...
                            msd_st_err(rng)/b^2,Color,Style);
                        x = (linspace(0,taur_temp,100))';
                        y = (x/tau0_temp).^0.5;
                        p = plot(x/taur_temp,y);
                        p.Color = Color;
                        p.LineStyle = '--';
                        p.LineWidth = 1.5;
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
            xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
            ylabel('MSD ($b^2$)','FontSize',FontSize,'Interpreter','latex')
            xlim([0 1])
            title(['$N =$ ',num2str(N_Kuhn)],'FontSize',FontSize,'Interpreter','latex')

            FileName = [Folder,'/',AddOn,'MSD Rouse inset vs t.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/',AddOn,'MSD Rouse inset vs t.fig'];
            saveas(gcf,FileName)
        end
    end
    close(wb1)
    close all
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotAllHistogramsCompare(OverrideCompare,MakeMoviesHistos,...
    HistoUnits,Package,damps,N_Kuhns,Np,Ns,T,Nb,ka_in,kd_in,f0,dts,b,N,...
    Separations)

global ConstantsFileName xlab FontSize DamperConversion BeadSpringOrMeso...
    LengthConversion

TotalSims = size(Package,1);
SimCt = 0;

wb1 = waitbar(0,'Plotting comparison histograms...');

NpOrig = Np;

% MSD wrt time for each damper while sweeping N
for dm = 1:length(damps)
    damp = damps(dm);
    dt = dts(dm);

    for nk = 1:length(N_Kuhns)
        N_Kuhn = N_Kuhns(nk);

        %Define movie names
        Folder = ['Output Plots/damp ',num2str(damp,'%.2e'),'.N ',num2str(N_Kuhn),...
            '.b ',num2str(b,'%.2f')];
        if ~isfolder(Folder)
            mkdir(Folder)
        end
        FigureName = [Folder,'/Norms - All Models'];

        for sp=1:length(Separations)
            Separation = Separations(sp);

            PackageTemp = Package(ismember(Package(:,2),NpOrig),:);
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
                'Plotting comparison histograms...')
            if ~isempty(PackageTemp) && ((MakeMoviesHistos==1 &&...
                    ~isfile([FigureName,'_stretch.png'])) ||...
                    OverrideCompare)
                SimCt = SimCt+1;

                % Callout Input Script
                EdgeFactor = 1;
                InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dt,damp,...
                    N_Kuhn,b,Separation);
                SetDirAndFileNames;
                DefineCompiledFileNames;
                Table = readtable(ConstantsFileName);
                constants = table2array(Table);
                N_Kuhn = constants(3); b = constants(4);

                %% For Bead-spring 
                % Set Np properly for Bead Model
                Np = 5^3;
                BSOM = 0;   %Toggle bead-spring or meso
                [dt,damp,~,~] = DefineTimeStep(b,LengthConversion,...
                    DamperConversion,BSOM);
                EdgeFactor = 1;
                InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dt,damp,...
                    N_Kuhn,b,Separation);
                SetDirAndFileNames;
                DefineCompiledFileNames;

                [rx_bead,ry_bead,rz_bead] = ImportEndToEnd(BSOM,Np,N,...
                    ka_in,kd_in,f0,dt,damp,N_Kuhn,b,Separation);

                
                % Crop initially stretched chains
                rx_bead(1:20,:) = [];
                ry_bead(1:20,:) = [];
                rz_bead(1:20,:) = [];

                norms_bead = vecnorm([rx_bead(:) ry_bead(:) rz_bead(:)],2,2);
                norms_bead(isnan(norms_bead)) = [];

                %% For Bead-spring 
                % Set Np properly for Bead Model
                Np = 11^3;
                BSOM = 1;       %Toggle bead-spring or meso
                [dt,damp,~,~] = DefineTimeStep(b,LengthConversion,...
                    DamperConversion,BSOM);
                EdgeFactor = 1;
                InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dt,damp,...
                    N_Kuhn,b,Separation);
                SetDirAndFileNames;
                DefineCompiledFileNames;

                [rx_meso,ry_meso,rz_meso] = ImportEndToEnd(BSOM,Np,N,...
                    ka_in,kd_in,f0,dt,damp,N_Kuhn,b,Separation);

                % Crop initially stretched chains
%                 cutoff_indx = round(0.9*size(rx_meso,1));
%                 rx_meso(1:cutoff_indx,:) = [];
%                 ry_meso(1:cutoff_indx,:) = [];
%                 rz_meso(1:cutoff_indx,:) = [];

                norms_meso = vecnorm([rx_meso(:) ry_meso(:) rz_meso(:)],2,2);
                norms_meso(isnan(norms_meso)) = [];

                %% Plot histogram with respect to stretch
                figure(1); clf; hold on
                WRTStretchOrLength = 0;
                AddOn = '_stretch';
                xlab = '$\lambda$';
                XLim = sqrt(N_Kuhn);
                YLim = 1.5;
                SymDist = 0;
                PlotIdeal = 0;
                
                % Mesoscale
                color = [0.75 0.75 0.75]; discret_or_histo = 1;
                [out(1),~] = Plot1DHistogram(norms_meso,color,1,XLim,YLim,N_Kuhn,b,...
                    SymDist,PlotIdeal,WRTStretchOrLength,HistoUnits,discret_or_histo);

                % Bead-spring
                color = 'r';
                PlotIdeal = 1;
                discret_or_histo = 0;
                [out(2),out(3)] = Plot1DHistogram(norms_bead,color,0.6,XLim,YLim,N_Kuhn,b,...
                    SymDist,PlotIdeal,WRTStretchOrLength,HistoUnits,discret_or_histo);

                l = legend(out,{'Mesoscale Model','Bead-spring Model','Gaussian PDF'});
                l.FontSize = FontSize/2;
                l.Interpreter = 'latex';
                l.Location = 'Northeast';

                set(gca,'FontSize',FontSize/1.5)
                xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
                ylabel('$p$','FontSize',FontSize,'Interpreter','latex')

                title(['$N$ =',num2str(N_Kuhn),...
                    ', $\gamma$ = ',num2str(damp*DamperConversion,'%.2e'),' N s m$^{-1}$'],'FontSize',FontSize/2,...
                    'Interpreter','latex')

                xlim([0 3])

                FileName = [FigureName,AddOn,'.png'];
                saveas(gcf,FileName)
                FileName = [FigureName,AddOn,'.fig'];
                saveas(gcf,FileName)
            end
        end
    end
end
close(wb1)
close all

if NpOrig==5^3
    BeadSpringOrMeso = 0;
else
    BeadSpringOrMeso = 1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rx,ry,rz] = ImportEndToEnd(BSOM,Np,N,ka_in,kd_in,f0,dt,damp,...
    N_Kuhn,b,Separation)

global BeadSpringOrMeso EndToEndDataFileName

BeadSpringOrMeso = BSOM;

% Callout Input Script
EdgeFactor = 1;
InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dt,damp,...
    N_Kuhn,b,Separation);
SetDirAndFileNames;
DefineCompiledFileNames;

EndToEndBS = load(EndToEndDataFileName,'-mat');
rx = EndToEndBS.rx;
ry = EndToEndBS.ry;
rz = EndToEndBS.rz;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotAllHistograms(MakeMoviesHistos,OverrideHistogramMovies,...
    WRTStretchOrLength,HistoUnits,Package,...
    N_Kuhns,Ns,T,Nb,ka_in,kd_in,f0,b,N,Separations)

global FileTagCompiled xlab...
    MovieName1 MovieFolder OutputFolderName...
    EndToEndDataFileName TimeStretchDataFileName BeadSpringOrMeso...
    LengthConversion DamperConversion

for BeadSpringOrMeso=0:1
    [dts,damps,~,~] = DefineTimeStep(b,LengthConversion,DamperConversion,...
        1);
    Np = 11^3;

    EdgeFactor = 1;
    InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dts(1),damps(1),...
        N_Kuhns(1),b,Separations(1));
    DefineMSDFileNames;

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

    wb1 = waitbar(0,'Plotting histograms...');

    % MSD wrt time for each damper while sweeping N
    for dm = 1:length(damps)
        for nk = 1:length(N_Kuhns)
%             for sp=1:length(Separations)
                damp = damps(dm);
                dt = dts(dm);
                N_Kuhn = N_Kuhns(nk);
                Separation = Separations(nk);
                Np = 11^3;

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

                    [dt,~,~,~] = DefineTimeStep(b,LengthConversion,DamperConversion,...
                        BeadSpringOrMeso);

                    if BeadSpringOrMeso==0
                        Np = 5^3;
                    else
                        Np = 11^3;
                    end

                    % Callout Input Script
                    % Initializes non-sweeping parameters
                    EdgeFactor = 1;
                    InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dt,damp,...
                        N_Kuhn,b,Separation);
                    SetDirAndFileNames;
                    DefineCompiledFileNames;
%                     Table = readtable(ConstantsFileName);
%                     constants = table2array(Table);
%                     N_Kuhn = constants(3); b = constants(4);

                    MovieFolderName = [OutputFolderName,'Movies'];
                    if ~isfolder(MovieFolderName)
                        mkdir(MovieFolderName)
                    end
                    %Define movie names
                    MovieFolder = [MovieFolderName,'/damp ',num2str(damp,'%.2e'),'.N ',num2str(N_Kuhn),...
                        '.b ',num2str(b,'%.2f')];
                    MovieName1 = [MovieFolder,'/End-to-end/Norms',AddOn,FileTagCompiled];

                    if MakeMoviesHistos==1 && (~isfile([MovieName1,'.avi']) || OverrideHistogramMovies==1)
                        % Determine allocation size

                        TimeStretch = load(TimeStretchDataFileName,'-mat');
                        time = TimeStretch.time;
                        stretch = TimeStretch.stretch;
                        EndToEndData = load(EndToEndDataFileName,'-mat');
                        rx = EndToEndData.rx;
                        ry = EndToEndData.ry;
                        rz = EndToEndData.rz;
                        bondtypes = EndToEndData.bondtypes;

                        %% Bond statistics
                        % Make histograms of end-to-end distributions (by bond type though)
                        PlotHistograms(rx,ry,rz,bondtypes,time,stretch,...
                            WRTStretchOrLength,HistoUnits);
                    end
                end
%             end
        end
    end
    close(wb1)
    close all
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GenerateMSDPlotsWithNSweptCompare(OverrideCompareData,Package,...
    damps,N_Kuhns,Np,Ns,T,Nb,ka_in,kd_in,f0,dts,b,N,Separations)

global FontSize LengthConversion R2_beadvsmeso_filename R2_msdr_filename...
    R2_msdt_filename LineWidth taur_filename PlotForFigures tau0...
    DamperConversion

%% Uncomment below 3 lines to plot figs as presented in manuscript results
%% But COMMENT BACK OFF AND RERUN THIS FUNCTION TO ATTAIN SURFACE PLOTS
% N_Kuhns = N_Kuhns([1 4 7]);
% dts = dts([1 4 7]);
% damps = damps([1 4 7]);

wb1 = waitbar(0,'Plotting comparison of MSDs...');

NpOrig = Np;
dt_orig = dts;
PlotForAppendix = 1;

TotalSims = size(Package,1);

R2_beadvsmeso = zeros(length(damps),length(N_Kuhns));
R2_msdr = zeros(length(damps),length(N_Kuhns));
R2_msdt = zeros(length(damps),length(N_Kuhns));
plateau_MSD_meso = zeros(length(damps),length(N_Kuhns));
plateau_MSD_SE_meso = zeros(length(damps),length(N_Kuhns));
plateau_MSD_bead = zeros(length(damps),length(N_Kuhns));
plateau_MSD_SE_bead = zeros(length(damps),length(N_Kuhns));
SimCt = 0;

ColorRange = (linspace(0,0.85,length(N_Kuhns)))';
% Colors = [ColorRange flipud(ColorRange) flipud(ColorRange)];
Colors = [zeros(size(ColorRange)) flipud(ColorRange) zeros(size(ColorRange))];

RelXMax = 15;       %In units of tauR
SmallSize = 250;
BigSize = 500;
BigXAsp = 2.2;

tauR = readmatrix(taur_filename);

% MSD wrt time for each damper while sweeping N
if ~isfile(R2_beadvsmeso_filename) || ...
        ~isfile(R2_msdr_filename) || ~isfile(R2_msdt_filename) || ...
        OverrideCompareData
    for dm = 1:length(damps)
        damp = damps(dm);
        dt = dts(dm);

        %Define folder name names
        Folder = ['Output Plots/damp ',num2str(damp,'%.2e')];
        if ~isfolder(Folder)
            mkdir(Folder)
        end
        FileName = [Folder,'/Compare Normalized MSD vs Time.png'];

        if  ~isfile(FileName) || OverrideCompareData==1
            figure(1); clf; hold on % msd
            figure(21); clf; hold on % msd in normalized units
            figure(2); clf; hold on % msd in tau0<t<taur
            figure(22); clf; hold on % Inset of fig. (2)
            figure(7); clf; hold on % relative error between bead spring and meso
            figure(8); clf; hold on % absolute error between bead spring and meso
            figure(9); clf; hold on % absolute error between bead spring and meso radial msd
            figure(10); clf; hold on % absolute error between bead spring and meso tang msd
            figure(11); clf; hold on % absolute error between bead spring and meso radial msd
            figure(12); clf; hold on % absolute error between bead spring and meso tang msd
            figure(31); clf; hold on % instantaenous radial MSD
            figure(32); clf; hold on % instantaenous tangent MSD
            figure(33); clf; hold on % radial MSD
            figure(34); clf; hold on % tangent MSD
            nk_ct = 0;
            for nk = 1:length(N_Kuhns)
                N_Kuhn = N_Kuhns(nk);
                Color = Colors(nk,:);

                figure(3); clf; hold on % msd_r instant
                figure(4); clf; hold on % msd_t instant
                figure(5); clf; hold on % msd_r
                figure(6); clf; hold on % msd_t

                for sp=1:length(Separations)
                    Separation = Separations(sp);

                    Np = NpOrig;

                    PackageTemp = Package(ismember(Package(:,2),Np),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,3),Ns),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,4),T),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,5),Nb),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,6),ka_in),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,7),kd_in),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,8),f0),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,9),dt_orig),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,10),damp),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,12),N_Kuhn),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,13),b),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,14),Separation),:);

                    waitbar(SimCt*size(PackageTemp,1)/TotalSims,wb1,...
                        'Plotting comparison of MSDs...')
                    if ~isempty(PackageTemp)
                        SimCt = SimCt+1;

                        %% For Bead-spring
                        % Set Np properly for Bead Model
                        Np = 5^3;
                        BSOM = 0;

                        [dt,damp,~,~] = DefineTimeStep(b,LengthConversion,...
                            DamperConversion,BSOM);

                        [time_bead,msd_bead_mean,msd_bead_err,...
                            msd_bead_r0,msd_bead_r0_err,...
                            msd_bead_t0,msd_bead_t0_err,...
                            msd_bead_r00,msd_bead_r00_err,...
                            msd_bead_t00,msd_bead_t00_err] = ...
                            ImportMSD(BSOM,Np,N,ka_in,kd_in,f0,dt,damp,...
                            N_Kuhn,b,Separation);

                        %% For Mesoscale
                        % Set Np properly for Mesoscale model
                        Np = 11^3;
                        BSOM = 1;

                        [dt,damp,~,~] = DefineTimeStep(b,LengthConversion,...
                            DamperConversion,BSOM);

                        [time_meso,msd_meso_mean,msd_meso_err,...
                            msd_meso_r0,msd_meso_r0_err,...
                            msd_meso_t0,msd_meso_t0_err,...
                            msd_meso_r00,msd_meso_r00_err,...
                            msd_meso_t00,msd_meso_t00_err] = ...
                            ImportMSD(BSOM,Np,N,ka_in,kd_in,f0,dt,damp,...
                            N_Kuhn,b,Separation);

                        D_temp = 1/damp;
                        tau0 = b^2/D_temp;

                        NoPts = 100;
                        Spacing = round(length(time_bead)/(NoPts-1));
                        PlotRangeAbs = (1:Spacing:length(time_bead))';

                        npts = 50;
                        plot_rng = round((linspace(1,length(time_bead),npts))');

                        dtPlot = tauR(dm,nk)/10;
                        dtAvbl = diff(time_bead(1:2));
                        interval = ceil(dtAvbl/dtPlot);
                        PlotRangeRel = (1:interval:length(time_bead))';
                        PlotRangeRel_err = PlotRangeRel(1:10:end);

                        %% Plot both discrete models
                        if (PlotForFigures && (nk==1 || nk==4 || nk==7)) ||...
                                ~PlotForFigures
                            nk_ct = nk_ct+1;
                            LegendEntries{nk_ct} = ['$N =$ ',num2str(N_Kuhn)];

                            figure(1);

                            Style = '-.';
                            p0(nk_ct) = PlotErrorBar(time_bead(PlotRangeAbs),...
                                msd_bead_mean(PlotRangeAbs)*LengthConversion^2,...
                                msd_bead_err(PlotRangeAbs)*LengthConversion^2,...
                                Color);
%                             [p0(nk_ct),~] = PlotCurve(time_bead(PlotRangeAbs),...
%                                 msd_bead_mean(PlotRangeAbs)*LengthConversion^2,...
%                                 msd_bead_err(PlotRangeAbs)*LengthConversion^2,...
%                                 Color,Style);

                            Style = '-';
                            [p1(nk_ct),~] = PlotCurve(time_meso(PlotRangeAbs),...
                                msd_meso_mean(PlotRangeAbs)*LengthConversion^2,...
                                msd_meso_err(PlotRangeAbs)*LengthConversion^2,...
                                Color,Style);

                            figure(21);
                            if PlotForAppendix==1
                                plot_rng_temp = PlotRangeAbs(1:3:end);
                            else
                                plot_rng_temp = PlotRangeAbs;
                            end

                            Style = '-.';
                            p0(nk_ct) = PlotErrorBar(time_bead(plot_rng_temp)/tau0,...
                                msd_bead_mean(plot_rng_temp)/b^2,...
                                msd_bead_err(plot_rng_temp)/b^2,...
                                Color);
%                             [p0(nk_ct),~] = PlotCurve(time_bead(PlotRangeAbs),...
%                                 msd_bead_mean(PlotRangeAbs)*LengthConversion^2,...
%                                 msd_bead_err(PlotRangeAbs)*LengthConversion^2,...
%                                 Color,Style);

                            Style = '-';
                            [p1(nk_ct),~] = PlotCurve(time_meso(plot_rng_temp)/tau0,...
                                msd_meso_mean(plot_rng_temp)/b^2,...
                                msd_meso_err(plot_rng_temp)/b^2,...
                                Color,Style);

                            % Calculate plateau MSD
                            rng = ((round(length(msd_meso_mean)/2)):length(msd_meso_mean))';
                            plateau_MSD_meso(dm,nk) = mean(msd_meso_mean(rng)/b^2);
                            plateau_MSD_SE_meso(dm,nk) = sqrt(sum(msd_meso_err(rng).^2));
%                             plateau_MSD_SE_meso(dm,nk) = std(msd_meso_mean(rng)/b^2,1)/sqrt(length(rng));
                            rng = ((round(length(msd_bead_mean)/2)):length(msd_bead_mean))';
                            plateau_MSD_bead(dm,nk) = mean(msd_bead_mean(rng)/b^2);
                            plateau_MSD_SE_bead(dm,nk) = sqrt(sum(msd_bead_err(rng).^2));
%                             plateau_MSD_SE_bead(dm,nk) = std(msd_bead_mean(rng)/b^2,1)/sqrt(length(rng));


                            %% Plot bead-spring with fitted Rouse model
                            figure(2)

                            max_t = 5*tauR(dm,nk);
                            npts = 50;

                            plot_rng = (linspace(1,length(time_bead(time_bead<=max_t)),npts))';
                            plot_rng = unique(round(plot_rng));
                            p2(nk_ct) = PlotErrorBar(time_bead(plot_rng)/tauR(dm,nk),...
                                msd_bead_mean(plot_rng)/b^2,...
                                msd_bead_err(plot_rng)/b^2,Color);
%                             [p2(nk_ct),~] = PlotCurve(time_bead/tauR(dm,nk),...
%                                 msd_bead_mean/b^2,msd_bead_err/b^2,Color,Style);

                            Style = '-';
                            [p3(nk_ct),~] = PlotCurve(time_meso/tauR(dm,nk),...
                                msd_meso_mean/b^2,msd_meso_err/b^2,Color,Style);

                            %% Plot inset of (2)
                            figure(22)

                            max_t = 5*tauR(dm,nk);
                            npts = 50;

                            t_rouse = (linspace(0,max_t/tauR(dm,nk),100))';
                            MSD_rouse = (t_rouse*N_Kuhn^2).^0.5;
                            p_temp = plot(t_rouse,MSD_rouse);
                            p_temp.LineStyle = '--';
                            p_temp.Color = Color;
                            p_temp.LineWidth = 1.5;

                            plot_rng = (linspace(1,length(time_bead(time_bead<=max_t)),npts))';
                            plot_rng = unique(round(plot_rng));
                            p22a(nk_ct) = PlotErrorBar(time_bead(plot_rng)/tauR(dm,nk),...
                                msd_bead_mean(plot_rng)/b^2,...
                                msd_bead_err(plot_rng)/b^2,Color);


                            Style = '-';
                            [p22b(nk_ct),~] = PlotCurve(time_meso/tauR(dm,nk),...
                                msd_meso_mean/b^2,msd_meso_err/b^2,Color,Style);

                        end

                        %% Plot instantaneous MSD_r for both discrete models
                        figure(3);

                        SubFolder = ['Output Plots/damp ',...
                            num2str(damp,'%.2e'),...
                            '.N ',num2str(N_Kuhn),...
                            '.b ',num2str(b,'%.2f')];
                        if ~isfolder(SubFolder)
                            mkdir(SubFolder)
                        end

                        npts = 50;
                        plot_rng = linspace(1,length(PlotRangeRel(time_bead/tauR(dm,nk)<=RelXMax)),npts)';
                        plot_rng = unique(round(plot_rng'));
                        p4(1) = PlotErrorBar(time_bead(plot_rng)/tauR(dm,nk),...
                            msd_bead_r0(plot_rng)/b^2,...
                            msd_bead_r0_err(plot_rng)/b^2,...
                            Color);

                        Style = '-';
                        [p4(2),~] = PlotCurve(time_meso(PlotRangeRel)/tauR(dm,nk),...
                            msd_meso_r0(PlotRangeRel)/b^2,...
                            msd_meso_r0_err(PlotRangeRel)/b^2,...
                            Color,Style);

                        set(gca,'FontSize',FontSize/1.5)
                        set(gcf,'Color','w')
                        pbaspect([1 1 1])
                        xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
                        ylabel('MSD$_r^0$ ($b^2$)',...
                            'FontSize',FontSize,'Interpreter','latex')

                        xlim([0 RelXMax])
                        set(gcf,'Position',[200 200 SmallSize SmallSize])

                        FileName = [SubFolder,'/MSD_r inst.png'];
                        saveas(gcf,FileName)
                        FileName = [SubFolder,'/MSD_r inst.fig'];
                        saveas(gcf,FileName)

                        %% Plot instantaneous MSD_t for both discrete models
                        figure(4);

                        p5(1) = PlotErrorBar(time_bead(plot_rng)/tauR(dm,nk),...
                            msd_bead_t0(plot_rng)/b^2,...
                            msd_bead_t0_err(plot_rng)/b^2,...
                            Color);

                        Style = '-';
                        [p5(2),~] = PlotCurve(time_meso(PlotRangeRel)/tauR(dm,nk),...
                            msd_meso_t0(PlotRangeRel)/b^2,...
                            msd_meso_t0_err(PlotRangeRel)/b^2,...
                            Color,Style);

                        l_entries = {'Bead-spring','Mesoscale'};
                        l = legend(p5,l_entries);
                        l.FontSize = FontSize/2;
                        l.Interpreter = 'latex';
                        l.Location = 'Best';

                        set(gca,'FontSize',FontSize/1.5)
                        set(gcf,'Color','w')
                        pbaspect([1 1 1])
                        xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
                        ylabel('MSD$_t^0$ ($b^2$)',...
                            'FontSize',FontSize,'Interpreter','latex')

                        xlim([0 RelXMax])
                        set(gcf,'Position',[200 200 SmallSize SmallSize])

                        FileName = [SubFolder,'/MSD_t inst.png'];
                        saveas(gcf,FileName)
                        FileName = [SubFolder,'/MSD_t inst.fig'];
                        saveas(gcf,FileName)

                        %% Plot MSD_r for both discrete models
                        figure(5);

                        p6(1) = PlotErrorBar(time_bead(plot_rng)/tauR(dm,nk),...
                            msd_bead_r00(plot_rng)/b^2,...
                            msd_bead_r00_err(plot_rng)/b^2,...
                            Color);

                        Style = '-';
                        [p6(2),~] = PlotCurve(time_meso(PlotRangeRel)/tauR(dm,nk),...
                            msd_meso_r00(PlotRangeRel)/b^2,...
                            msd_meso_r00_err(PlotRangeRel)/b^2,...
                            Color,Style);

                        l = legend(p6,l_entries);
                        l.FontSize = FontSize/2;
                        l.Interpreter = 'latex';
                        l.Location = 'Best';

                        set(gca,'FontSize',FontSize/1.5)
                        set(gcf,'Color','w')
                        pbaspect([1 1 1])
                        xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
                        ylabel('MSD$_r$ ($b^2$)',...
                            'FontSize',FontSize,'Interpreter','latex')

                        set(gcf,'Position',[200 200 SmallSize SmallSize])
                        xlim([0 RelXMax])
                        ylim([0 Inf])
                        
                        FileName = [SubFolder,'/MSD_r.png'];
                        saveas(gcf,FileName)
                        FileName = [SubFolder,'/MSD_r.fig'];
                        saveas(gcf,FileName)
                        
                        %% Plot MSD_t for both discrete models
                        figure(6);

                        p7(1) = PlotErrorBar(time_bead(plot_rng)/tauR(dm,nk),...
                            msd_bead_t00(plot_rng)/b^2,...
                            msd_bead_t00_err(plot_rng)/b^2,...
                            Color);

                        Style = '-';
                        [p7(2),~] = PlotCurve(time_meso(PlotRangeRel)/tauR(dm,nk),...
                            msd_meso_t00(PlotRangeRel)/b^2,...
                            msd_meso_t00_err(PlotRangeRel)/b^2,...
                            Color,Style);

                        set(gca,'FontSize',FontSize/1.5)
                        set(gcf,'Color','w')
                        pbaspect([1 1 1])
                        xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
                        ylabel('MSD$_t$ ($b^2$)',...
                            'FontSize',FontSize,'Interpreter','latex')

                        set(gcf,'Position',[200 200 SmallSize SmallSize])
                        xlim([0 RelXMax])

                        FileName = [SubFolder,'/MSD_t.png'];
                        saveas(gcf,FileName)
                        FileName = [SubFolder,'/MSD_t.fig'];
                        saveas(gcf,FileName)


                        %% Plot instantaneous MSD_r for both discrete models
                        if nk==1 || nk==length(N_Kuhns)
                            figure(31);

                            npts = 50;
                            plot_rng = linspace(1,length(PlotRangeRel(time_bead/tauR(dm,nk)<=RelXMax)),npts)';
                            plot_rng = unique(round(plot_rng'));
                            p31(1) = PlotErrorBar(time_bead(plot_rng)/tauR(dm,nk),...
                                msd_bead_r0(plot_rng)/b^2,...
                                msd_bead_r0_err(plot_rng)/b^2,...
                                Color);

                            Style = '-';
                            [p31(2),~] = PlotCurve(time_meso(PlotRangeRel)/tauR(dm,nk),...
                                msd_meso_r0(PlotRangeRel)/b^2,...
                                msd_meso_r0_err(PlotRangeRel)/b^2,...
                                Color,Style);

                            set(gca,'FontSize',FontSize/1.5)
                            set(gcf,'Color','w')
                            pbaspect([1 1 1])
                            xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
%                             ylabel('$\Delta r_r^2(\tau_0/10)$ ($b^2$)',...
%                                 'FontSize',FontSize,'Interpreter','latex')
                            ylabel('$\Delta r_r^2(\Delta t)$ ($b^2$)',...
                                'FontSize',FontSize,'Interpreter','latex')
                            xlim([0 RelXMax])
                            set(gcf,'Position',[200 200 SmallSize SmallSize])

                            figure(32);

                            p32(1) = PlotErrorBar(time_bead(plot_rng)/tauR(dm,nk),...
                                msd_bead_t0(plot_rng)/b^2,...
                                msd_bead_t0_err(plot_rng)/b^2,...
                                Color);

                            Style = '-';
                            [p32(2),~] = PlotCurve(time_meso(PlotRangeRel)/tauR(dm,nk),...
                                msd_meso_t0(PlotRangeRel)/b^2,...
                                msd_meso_t0_err(PlotRangeRel)/b^2,...
                                Color,Style);

                            set(gca,'FontSize',FontSize/1.5)
                            set(gcf,'Color','w')
                            pbaspect([1 1 1])
                            xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
%                             ylabel('$\Delta r_t^2(\tau_0/10)$ ($b^2$)',...
%                                 'FontSize',FontSize,'Interpreter','latex')
                            ylabel('$\Delta r_t^2(\Delta t)$ ($b^2$)',...
                                'FontSize',FontSize,'Interpreter','latex')

                            xlim([0 RelXMax])
                            set(gcf,'Position',[200 200 SmallSize SmallSize])

                            figure(33);

                            p33(1) = PlotErrorBar(time_bead(plot_rng)/tauR(dm,nk),...
                                msd_bead_r00(plot_rng)/b^2,...
                                msd_bead_r00_err(plot_rng)/b^2,...
                                Color);

                            Style = '-';
                            [p33(2),~] = PlotCurve(time_meso(PlotRangeRel)/tauR(dm,nk),...
                                msd_meso_r00(PlotRangeRel)/b^2,...
                                msd_meso_r00_err(PlotRangeRel)/b^2,...
                                Color,Style);

                            set(gca,'FontSize',FontSize/1.5)
                            set(gcf,'Color','w')
                            pbaspect([1 1 1])
                            xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
                            ylabel('$\langle \Delta r_r^2 (t) \rangle$ ($b^2$)',...
                                'FontSize',FontSize,'Interpreter','latex')

                            set(gcf,'Position',[200 200 SmallSize SmallSize])
                            xlim([0 RelXMax])
                            ylim([0 Inf])

                            figure(34);

                            p34(1) = PlotErrorBar(time_bead(plot_rng)/tauR(dm,nk),...
                                msd_bead_t00(plot_rng)/b^2,...
                                msd_bead_t00_err(plot_rng)/b^2,...
                                Color);

                            Style = '-';
                            [p34(2),~] = PlotCurve(time_meso(PlotRangeRel)/tauR(dm,nk),...
                                msd_meso_t00(PlotRangeRel)/b^2,...
                                msd_meso_t00_err(PlotRangeRel)/b^2,...
                                Color,Style);

                            set(gca,'FontSize',FontSize/1.5)
                            set(gcf,'Color','w')
                            pbaspect([1 1 1])
                            xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
                            ylabel('$\langle \Delta r_t^2 (t) \rangle$ ($b^2$)',...
                                'FontSize',FontSize,'Interpreter','latex')

                            set(gcf,'Position',[200 200 SmallSize SmallSize])
                            xlim([0 RelXMax])
                        end

                        if (PlotForFigures && (nk==1 || nk==4 || nk==7)) ||...
                                ~PlotForFigures
                            % Compute relative error
                            figure(7)
                            if nk==1
                                plot(time_meso/tauR(dm,nk),...
                                    zeros(size(time_meso)),'k:')
                            end
                            fact = 10;
                            mmmeso = movmean(msd_meso_mean,1);
                            mmbead = movmean(msd_bead_mean,1);
                            relative_err_mean = (mmmeso-mmbead)./mmbead;
                            p = plot(time_meso(PlotRangeRel)/tauR(dm,nk),...
                                relative_err_mean(PlotRangeRel));
                            p.Color = Color;
                            p.LineWidth = LineWidth;

                            % Compute absolute error
                            figure(8)
                            if nk==1
                                plot(time_meso/tauR(dm,nk),...
                                    zeros(size(time_meso)),'k:')
                            end
                            mmmeso = movmean(msd_meso_mean,1);
                            mmbead = movmean(msd_bead_mean,1);
                            absolute_err_mean = (mmmeso-mmbead)./b^2;
                            p = plot(time_meso(PlotRangeRel_err)/tauR(dm,nk),...
                                absolute_err_mean(PlotRangeRel_err));
                            p.Color = Color;
                            p.LineWidth = LineWidth;

                            % Compute absolute error
                            figure(9)
                            if nk==1
                                plot(time_meso/tauR(dm,nk),...
                                    zeros(size(time_meso)),'k:')
                            end
                            mmmeso = movmean(msd_meso_r0,1);
                            mmbead = movmean(msd_bead_r0,1);
                            absolute_err_r0 = (mmmeso-mmbead)./b^2;
                            p = plot(time_meso(PlotRangeRel)/tauR(dm,nk),...
                                absolute_err_r0(PlotRangeRel));
                            p.Color = Color;
                            p.LineWidth = LineWidth;

                            % Compute absolute error
                            figure(10)
                            if nk==1
                                plot(time_meso/tauR(dm,nk),...
                                    zeros(size(time_meso)),'k:')
                            end
                            mmmeso = movmean(msd_meso_t0,11);
                            mmbead = movmean(msd_bead_t0,fact);
                            absolute_err_t0 = (mmmeso-mmbead)./b^2;
                            p = plot(time_meso(PlotRangeRel)/tauR(dm,nk),...
                                absolute_err_t0(PlotRangeRel));
                            p.Color = Color;
                            p.LineWidth = LineWidth;

                            % Compute absolute error
                            figure(11)
                            if nk==1
                                plot(time_meso/tauR(dm,nk),...
                                    zeros(size(time_meso)),'k:')
                            end
                            mmmeso = movmean(msd_meso_r00,1);
                            mmbead = movmean(msd_bead_r00,1);
%                             absolute_err_r00 = (mmmeso-mmbead)./b^2;
                            absolute_err_r00 = (mmmeso-mmbead)./mmbead;
                            p = plot(time_meso(PlotRangeRel)/tauR(dm,nk),...
                                absolute_err_r00(PlotRangeRel));
                            p.Color = Color;
                            p.LineWidth = LineWidth;

                            % Compute absolute error
                            figure(12)
                            if nk==1
                                plot(time_meso/tauR(dm,nk),...
                                    zeros(size(time_meso)),'k:')
                            end
                            mmmeso = movmean(msd_meso_t00,11);
                            mmbead = movmean(msd_bead_t00,fact);
%                             absolute_err_t00 = (mmmeso-mmbead)./b^2;
                            absolute_err_t00 = (mmmeso-mmbead)./mmbead;
                            p = plot(time_meso(PlotRangeRel)/tauR(dm,nk),...
                                absolute_err_t00(PlotRangeRel));
                            p.Color = Color;
                            p.LineWidth = LineWidth;
                        end

                        %% Compute R2 between models
                        mmint = 1;
                        rednoise_msd_bead = movmean(msd_bead_mean,mmint);
                        rednoise_msd_meso = movmean(msd_meso_mean,mmint);
                        rednoise_msdr_bead = movmean(msd_bead_r00,mmint);
                        rednoise_msdr_meso = movmean(msd_meso_r00,mmint);
                        rednoise_msdt_bead = movmean(msd_bead_t00,mmint);
                        rednoise_msdt_meso = movmean(msd_meso_t00,mmint);
                        Nsamp = 25;
                        samprange = round(linspace(1,length(time_meso),Nsamp));

                        y = rednoise_msd_bead(samprange);
                        f = rednoise_msd_meso(samprange);
                        R2_beadvsmeso(dm,nk) = ComputeR2(y,f);

                        y = rednoise_msdr_bead(samprange);
                        f = rednoise_msdr_meso(samprange);
                        R2_msdr(dm,nk) = ComputeR2(y,f);

                        y = rednoise_msdt_bead(samprange);
                        f = rednoise_msdt_meso(samprange);
                        R2_msdt(dm,nk) = ComputeR2(y,f);

                    end
                end
            end
            figure(1)

            l = legend(p1,LegendEntries);
            l.FontSize = FontSize/2;
            l.Interpreter = 'latex';
            l.Location = 'Southeast';

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([BigXAsp 1 1])
            xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
            ylabel('$\langle \Delta r^2 (t) \rangle$ (nm$^2$)','FontSize',FontSize,'Interpreter','latex')

            set(gcf,'Position',[200 200 BigSize BigSize])

            FileName = [Folder,'/Compare MSD vs Time.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/Compare MSD vs Time.fig'];
            saveas(gcf,FileName)

            figure(21)

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            xlabel('$t$ ($\tau_0$)','FontSize',FontSize,'Interpreter','latex')
            ylabel('$\langle \Delta  r^2 (t) \rangle$ ($b^2$)','FontSize',FontSize,'Interpreter','latex')

            if PlotForAppendix==1
                pbaspect([1 1 1])
                set(gcf,'Position',[200 200 SmallSize SmallSize])

                FileName = [Folder,'/Appendix.Compare MSD vs Time.png'];
                saveas(gcf,FileName)
                FileName = [Folder,'/Appendix.Compare MSD vs Time.fig'];
                saveas(gcf,FileName)
            else
                pbaspect([BigXAsp 1 1])
                set(gcf,'Position',[200 200 BigSize BigSize])

                FileName = [Folder,'/Compare MSD vs Time.png'];
                saveas(gcf,FileName)
                FileName = [Folder,'/Compare MSD vs Time.fig'];
                saveas(gcf,FileName)
            end

            figure(2)
            l = legend(p3,LegendEntries);
            l.FontSize = FontSize/2;
            l.Interpreter = 'latex';
            l.Location = 'Southeast';

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([BigXAsp 1 1])
            xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
            ylabel('$\langle \Delta r^2 (t) \rangle$ ($b^2$)','FontSize',FontSize,'Interpreter','latex')

            pbaspect([1.4 1 1])
            ylim([0 120])
            set(gcf,'Position',[200 200 BigSize BigSize])

            FileName = [Folder,'/Compare Normalized MSD vs Time.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/Compare Normalized MSD vs Time.fig'];
            saveas(gcf,FileName)
            xlim([0 RelXMax/3])


            figure(22)
            l = legend(p5,LegendEntries);
            l.FontSize = FontSize/2;
            l.Interpreter = 'latex';
            l.Location = 'Southeast';

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([1 1 1])

            set(gcf,'Position',[200 200 150 150])

            xlim([0 1])
            ylim([0 40])

            FileName = [Folder,'/Compare Normalized MSD vs Time.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/Compare Normalized MSD vs Time.fig'];
            saveas(gcf,FileName)


            figure(7)

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
            ylabel('Rel. err.','FontSize',FontSize,'Interpreter','latex')

            xlim([0 RelXMax/3])

            pbaspect([4 1 1])
            set(gcf,'Position',[200 200 BigSize BigSize])


            FileName = [Folder,'/Relative error of Meso vs BeadSpring MSD vs Time.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/Relative error of Meso vs BeadSpring MSD vs Time.fig'];
            saveas(gcf,FileName)


            figure(8)

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([1 1 1])
            xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
            ylabel('Abs. err. ($b^2$)','FontSize',FontSize,'Interpreter','latex')

            set(gcf,'Position',[200 200 BigSize BigSize])

            FileName = [Folder,'/Absolute error of Meso vs BeadSpring MSD vs Time.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/Absolute error of Meso vs BeadSpring MSD vs Time.fig'];
            saveas(gcf,FileName)
            
            figure(9)

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([3 1 1])
            xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
%             ylabel('$\Delta_r^0$ ($b^2$)','FontSize',FontSize,'Interpreter','latex')
            ylabel('Err. ($b^2$)','FontSize',FontSize,'Interpreter','latex')

            val = mean(absolute_err_r0(end-100:end));
            plot([0 RelXMax],[val val],'k--')
            
            xlim([0 RelXMax])
            ylim([-2 0])
            set(gcf,'Position',[200 200 1.05*SmallSize SmallSize])

            FileName = [Folder,'/Absolute error of Meso vs BeadSpring Rad Inst. MSD vs Time.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/Absolute error of Meso vs BeadSpring Rad Inst. MSD vs Time.fig'];
            saveas(gcf,FileName)
                    
            figure(10)

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([3 1 1])
            xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
%             ylabel('$\Delta_t^0$ ($b^2$)','FontSize',FontSize,'Interpreter','latex')
            ylabel('Err. ($b^2$)','FontSize',FontSize,'Interpreter','latex')

            val = mean(absolute_err_t0(end-100:end));
            plot([0 RelXMax],[val val],'k--')
            
            xlim([0 RelXMax])
            ylim([-3 0])
            set(gcf,'Position',[200 200 1.05*SmallSize SmallSize])

            FileName = [Folder,'/Absolute error of Meso vs BeadSpring Tan Inst. MSD vs Time.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/Absolute error of Meso vs BeadSpring Tan Inst. MSD vs Time.fig'];
            saveas(gcf,FileName)
            
            figure(11)

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([3 1 1])
            xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
%             ylabel('$\Delta_r$ ($b^2$)','FontSize',FontSize,'Interpreter','latex')
%             ylabel('Err. ($b^2$)','FontSize',FontSize,'Interpreter','latex')
            ylabel('Rel. err.','FontSize',FontSize,'Interpreter','latex')

            xlim([0 RelXMax])
            set(gcf,'Position',[200 200 SmallSize SmallSize])

            FileName = [Folder,'/Absolute error of Meso vs BeadSpring Rad MSD vs Time.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/Absolute error of Meso vs BeadSpring Rad MSD vs Time.fig'];
            saveas(gcf,FileName)
                    
            figure(12)

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([3 1 1])
            xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
%             ylabel('$\Delta_t$ ($b^2$)','FontSize',FontSize,'Interpreter','latex')
            ylabel('Rel. err.','FontSize',FontSize,'Interpreter','latex')

            xlim([0 RelXMax])
            set(gcf,'Position',[200 200 SmallSize SmallSize])

            FileName = [Folder,'/Absolute error of Meso vs BeadSpring Tan MSD vs Time.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/Absolute error of Meso vs BeadSpring Tan MSD vs Time.fig'];
            saveas(gcf,FileName)

            figure(13); clf; hold on

            e = errorbar(N_Kuhns,plateau_MSD_bead(dm,:),plateau_MSD_SE_bead(dm,:));
            e.Marker = 'o';
            e.MarkerFaceColor = 'k';
            e.MarkerEdgeColor = 'k';
            e.Color = 'k';
            e.LineStyle = 'none';

            e = errorbar(N_Kuhns,plateau_MSD_meso(dm,:),plateau_MSD_SE_meso(dm,:));
            e.Marker = '^';
            e.MarkerFaceColor = 'k';
            e.MarkerEdgeColor = 'k';
            e.Color = 'k';
            e.LineStyle = 'none';            


            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([1 1 1])

            xlabel('$N$','FontSize',FontSize,'Interpreter','latex')
            ylabel('$\langle \Delta r^2 (t) \rangle_{ss}$ ($b^2$)','FontSize',FontSize,'Interpreter','latex')
            set(gcf,'Position',[200 200 SmallSize SmallSize])

            x = [N_Kuhns;N_Kuhns];
            y = [plateau_MSD_bead(dm,:)';plateau_MSD_meso(dm,:)'];
            A_nom = 1; B_nom = 0;
            StatusString = 'Fitting line';
            xlab_str = '$N$';
            ylab_str = 'MSD';
            [A,B,R2] = FitLine(x,y,A_nom,B_nom,NoPts,...
                StatusString,xlab_str,ylab_str);
            y = A*N_Kuhns+B;
            plot(N_Kuhns,y,'k--')

            FileName = [Folder,'/SS MSD vs N.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/SS MSD vs N.fig'];
            saveas(gcf,FileName)

            figure(31)
            FileName = [Folder,'/MSD0_r.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/MSD0_r.fig'];
            saveas(gcf,FileName)

            figure(31)
            FileName = [Folder,'/MSD0_t.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/MSD0_t.fig'];
            saveas(gcf,FileName)

        end
    end
    writematrix(R2_beadvsmeso,R2_beadvsmeso_filename)
    writematrix(R2_msdr,R2_msdr_filename)
    writematrix(R2_msdt,R2_msdt_filename)
end

% Plot rouse timescale as function of N and D
D = 1./damps;
tau0_vals = b^2./D;

if length(tau0_vals)>1
    figure(22); clf; hold on

    [tau0_vals_grid,N_Kuhns_grid] = meshgrid(tau0_vals,N_Kuhns);
    s = surface(tau0_vals_grid,N_Kuhns_grid,tauR'./tau0_vals_grid);
    s.FaceColor = 'interp';

    set(gca,'FontSize',FontSize/1.5)
    set(gca,'xscale','log')
    xlabel('$\tau_0$ (s)','FontSize',FontSize,'Interpreter','latex')
    ylabel('$N$','FontSize',FontSize,'Interpreter','latex')

    xlim([min(tau0_vals),max(tau0_vals)])
    ylim([min(N_Kuhns),max(N_Kuhns)])
    yticks(N_Kuhns)

    c = colorbar;
    clab = c.Label;
    clab.String = '$\tau_r^*$ ($\tau_0$)';
    clab.FontSize = FontSize;
    clab.Interpreter = 'latex';

    set(gcf,'color','w')
    set(gcf,'Position',[100 100 500 500])
    pbaspect([1 1 1])


    figure(23); clf; hold on

    [tau0_vals_grid,N_Kuhns_grid] = meshgrid(tau0_vals,N_Kuhns);
    s = surface(tau0_vals_grid,N_Kuhns_grid,tauR'./tau0_vals_grid./N_Kuhns_grid.^2);
    s.FaceColor = 'interp';

    set(gca,'FontSize',FontSize/1.5)
    set(gca,'xscale','log')
    xlabel('$\tau_0$ (s)','FontSize',FontSize,'Interpreter','latex')
    ylabel('$N$','FontSize',FontSize,'Interpreter','latex')

    xlim([min(tau0_vals),max(tau0_vals)])
    ylim([min(N_Kuhns),max(N_Kuhns)])
    yticks(N_Kuhns)

    c = colorbar;
    clab = c.Label;
    clab.String = '$\tau_s^* (\tau_0)$';
    clab.FontSize = FontSize;
    clab.Interpreter = 'latex';
    caxis([0 0.1])

    set(gcf,'color','w')
    set(gcf,'Position',[100 100 500 500])
    pbaspect([1 1 1])

    all_taus_star = tauR'./tau0_vals_grid./N_Kuhns_grid.^2;
    avg_taus_star = mean(all_taus_star(:));
    se_taus_star = std(all_taus_star(:))/sqrt(length(all_taus_star(:)));
end

close(wb1)
close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R2 = ComputeR2(y,f)

RSS = sum((y-f).^2);
TSS = sum((y-mean(y)).^2);
R2 = 1-RSS/TSS;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time,msd_mean,msd_err,...
    msd_r0,msd_r0_err,...
    msd_t0,msd_t0_err,...
    msd_r00,msd_r00_err,...
    msd_t00,msd_t00_err] = ...
    ImportMSD(BSOM,Np,N,ka_in,kd_in,f0,dt,damp,...
    N_Kuhn,b,Separation)

global BeadSpringOrMeso TimeStretchDataFileName MSDDataFileName

BeadSpringOrMeso = BSOM;

% Callout Input Script
EdgeFactor = 1;
InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dt,damp,...
    N_Kuhn,b,Separation);
SetDirAndFileNames;
DefineCompiledFileNames;

TimeStretch = load(TimeStretchDataFileName,'-mat');
time = TimeStretch.time;
MSD = load(MSDDataFileName,'-mat');

msd_mean = MSD.msd;
msd_err = MSD.msd_err;
msd_r0 = MSD.msd_st_r0;
msd_r0_err = MSD.msd_st_r0_err;
msd_t0 = MSD.msd_st_t0;
msd_t0_err = MSD.msd_st_t0_err;
msd_r00 = MSD.msd_st_r00;
msd_r00_err = MSD.msd_st_r00_err;
msd_t00 = MSD.msd_st_t00;
msd_t00_err = MSD.msd_st_t00_err;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R2_beadvsmeso,R2_mesovsrouse,R2_msdr,R2_msdt] = ...
    GeneratePlotsWithNSweptCompare(OverrideCompareData,Override,Package,...
    damps,N_Kuhns,Np,Ns,T,Nb,ka_in,kd_in,f0,dts,b,N,Separations)

global FontSize DamperConversion BeadSpringOrMeso LengthConversion...
    R2_beadvsmeso_filename R2_mesovsrouse_filename...
    R2_msdr_filename R2_msdt_filename 

% if length(N_Kuhns)>2
%     N_Kuhns(2:2:end) = [];
% end

wb1 = waitbar(0,'Plotting comparison of MSDs...');

NpOrig = Np;

TotalSims = size(Package,1);

tau0 = zeros(length(damps),length(N_Kuhns));
tau0_eff = zeros(length(damps),length(N_Kuhns));
R2 = zeros(length(damps),length(N_Kuhns));
tauR = zeros(length(damps),length(N_Kuhns));
MSD_ss = zeros(length(damps),length(N_Kuhns));
R2_beadvsmeso = zeros(size(tau0));
R2_mesovsrouse = zeros(size(tau0));
R2_msdr = zeros(length(damps),length(N_Kuhns));
R2_msdt = zeros(length(damps),length(N_Kuhns));
SimCt = 0;

ColorRange = (linspace(0,1,length(N_Kuhns)))';
Colors = [ColorRange flipud(ColorRange) flipud(ColorRange)];


% MSD wrt time for each damper while sweeping N
if ~isfile(R2_beadvsmeso_filename) || ...
        ~isfile(R2_msdr_filename) || ~isfile(R2_msdt_filename) || ...
        OverrideCompareData==1
    for dm = 1:length(damps)
        damp = damps(dm);
        dt = dts(dm);

        %Define movie names
        Folder = ['Output Plots/damp ',num2str(damp,'%.2e')];
        if ~isfolder(Folder)
            mkdir(Folder)
        end
        FileName = [Folder,'/Compare Normalized MSD vs Time.png'];
%         Folder = 'Output Plots';
%         FigureTag = ['Output Plots/damp ',num2str(damp,'%.2e'),'.b ',num2str(b,'%.2f')];
%         if ~isfolder(Folder)
%             mkdir(Folder)
%         end
%         FileName = [FigureTag,'.Normalized MSD vs Time.png'];

        if  ~isfile(FileName) || OverrideCompareData==1
            figure(1); clf; hold on % msd
            figure(2); clf; hold on % msd in tau0<t<taur
            for nk = 1:length(N_Kuhns)
                N_Kuhn = N_Kuhns(nk);
                Color = Colors(nk,:);

                LegendEntries{nk} = ['$N =$ ',num2str(N_Kuhn)];

                [tau0_eff(dm,nk),tauR_eff,R_N] = ComputeEffectiveDamper(N_Kuhn,b,damp);

                figure(3); clf; hold on % msd_r instant
                figure(4); clf; hold on % msd_t instant
                figure(5); clf; hold on % msd_r
                figure(6); clf; hold on % msd_t

                for sp=1:length(Separations)
                    Separation = Separations(sp);

                    Np = NpOrig;

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
                        'Plotting comparison of MSDs...')
                    if ~isempty(PackageTemp)
                        SimCt = SimCt+1;

                        %% For Bead-spring
                        % Set Np properly for Bead Model
                        Np = 5^3;
                        BeadSpringOrMeso = 0;

                        % Callout Input Script
                        EdgeFactor = 1;
                        InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dt,damp,...
                            N_Kuhn,b,Separation);
                        SetDirAndFileNames;
                        DefineCompiledFileNames;

                        % Determine allocation size
                        [Corners,N_steps] = ...
                            ImportDomainBoundaries(EdgeFactor,Np,N,ka_in,kd_in,f0,...
                            dt,damp,N_Kuhn,b,Separation);
                        N_samples = size(PackageTemp,1);

                        %% Assemble or read all of the ensemble data
                        [time_bead,~,~,~,~,~,msd_bead_mean,msd_bead_err,...
                            msd_bead_r0,msd_bead_r0_err,...
                            msd_bead_t0,msd_bead_t0_err,...
                            msd_bead_r00,msd_bead_r00_err,...
                            msd_bead_t00,msd_bead_t00_err,...
                            ~,~,~,~,~,~]...
                            = AssembleTheDataBeadSpring(Corners,EdgeFactor,Np,N,...
                            ka_in,kd_in,f0,dt,damp,N_Kuhn,b,Separation,...
                            PackageTemp,Override,N_steps,N_samples);

                        %% For Mesoscale
                        % Set Np properly for Mesoscale model
                        Np = 11^3;
                        BeadSpringOrMeso = 1;

                        % Callout Input Script
                        EdgeFactor = 1;
                        InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dt,damp,...
                            N_Kuhn,b,Separation);
                        SetDirAndFileNames;
                        DefineCompiledFileNames;

                        % Determine allocation size
                        [Corners,N_steps] = ...
                            ImportDomainBoundaries(EdgeFactor,Np,N,ka_in,kd_in,f0,...
                            dt,damp,N_Kuhn,b,Separation);
                        N_samples = size(PackageTemp,1);
                        [time_meso,~,~,~,~,~,msd_meso_mean,msd_meso_err,...
                            msd_meso_r0,msd_meso_r0_err,...
                            msd_meso_t0,msd_meso_t0_err,...
                            msd_meso_r00,msd_meso_r00_err,...
                            msd_meso_t00,msd_meso_t00_err,...
                            ~,~,~,~,~,~]...
                            = AssembleTheDataMeso(Corners,EdgeFactor,Np,N,...
                            ka_in,kd_in,f0,dt,damp,N_Kuhn,b,Separation,...
                            PackageTemp,Override,N_steps,N_samples);

                        [R2(dm,nk),tau0(dm,nk),tauR(dm,nk),MSD_ss(dm,nk)] =...
                            FitRouseModel(time_bead,msd_bead_mean,...
                            tau0_eff(dm,nk),tauR_eff,R_N,N_Kuhn);

                        %% Plot both discrete models
                        figure(1);

                        Style = '-.';
                        [p0(nk),~] = PlotCurve(time_bead,msd_bead_mean*LengthConversion^2,...
                            msd_bead_err*LengthConversion^2,Color,Style);

                        Style = '-';
                        [p1(nk),~] = PlotCurve(time_meso,msd_meso_mean*LengthConversion^2,...
                            msd_meso_err*LengthConversion^2,Color,Style);

                        x = time_bead;
%                         x(x>tauR(dm,nk)) = [];
                        y = (R_N^2)*(x/tau0(dm,nk)).^0.5;
                        p = plot(x,y*LengthConversion^2);
                        p.Color = Color;
                        p.LineStyle = '--';
                        p.LineWidth = 1.5;


                        %% Plot bead-spring with fitted Rouse model
                        figure(2)

                        Style = '-.';
                        [p2(nk),~] = PlotCurve(time_bead/tauR(dm,nk),...
                            msd_bead_mean/b^2,msd_bead_err/b^2,Color,Style);

                        Style = '-';
                        [p3(nk),~] = PlotCurve(time_meso/tauR(dm,nk),...
                            msd_meso_mean/b^2,msd_meso_err/b^2,Color,Style);

                        x = time_bead;
                        x(x>tauR(dm,nk)) = [];
                        y = (R_N^2)*(x/tau0(dm,nk)).^0.5;
                        p = plot(x/tauR(dm,nk),y/b^2);
                        p.Color = Color;
                        p.LineStyle = '--';
                        p.LineWidth = 1.5;

                        y = msd_bead_mean;
                        f = msd_meso_mean;
                        RSS = sum((y-f).^2);
                        TSS = sum((y-mean(y)).^2);
                        R2_beadvsmeso(dm,nk) = 1-RSS/TSS;

                        y = (R_N^2)*(x/tau0(dm,nk)).^0.5;
                        f = msd_meso_mean;
                        f(time_meso>tauR(dm,nk)) = [];
                        RSS = sum((y-f).^2);
                        TSS = sum((y-mean(y)).^2);
                        R2_mesovsrouse(dm,nk) = 1-RSS/TSS;

                        %% Plot instantaneous MSD_r for both discrete models
                        figure(3);

                        SubFolder = ['Output Plots/damp ',num2str(damp,'%.2e'),'.N ',num2str(N_Kuhn),...
                            '.b ',num2str(b,'%.2f')];
                        if ~isfolder(SubFolder)
                            mkdir(SubFolder)
                        end
%                         FigureTag = [SubFolder,'/damp ',num2str(damp,'%.2e'),'.b ',num2str(b,'%.2f')];
%                         ImgTag = [FigureTag,'.N ',num2str(N_Kuhn)];

                        Style = '-.';
                        [p4(1),~] = PlotCurve(time_bead,msd_bead_r0*LengthConversion^2,...
                            msd_bead_r0_err*LengthConversion^2,Color,Style);

                        Style = '-';
                        [p4(2),~] = PlotCurve(time_meso,msd_meso_r0*LengthConversion^2,...
                            msd_meso_r0_err*LengthConversion^2,Color,Style);

                        l_entries = {'Bead-spring','Mesoscale'};
                        l = legend(p4,l_entries);
                        l.FontSize = FontSize/2;
                        l.Interpreter = 'latex';
                        l.Location = 'Best';

                        set(gca,'FontSize',FontSize/1.5)
                        set(gcf,'Color','w')
                        pbaspect([1 1 1])
                        xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
                        ylabel('MSD$_r(\tau)$ (nm$^2$)','FontSize',FontSize,'Interpreter','latex')
%                         title(['$\gamma$ = ',num2str(damp,'%.2e'),' [force $\cdot$ time/length], ',...
%                             num2str(damp*DamperConversion,'%.2e'),' [N s m$^{-1}$]'],...
%                             'FontSize',FontSize/3,'Interpreter','latex')
                        %     xlim([0 1])

                        set(gcf,'Position',[200 200 275 275])

%                         FileName = [ImgTag,'.MSD_r inst.png'];
%                         saveas(gcf,FileName)
%                         FileName = [ImgTag,'.MSD_r inst.fig'];
%                         saveas(gcf,FileName)
                        FileName = [SubFolder,'/MSD_r inst.png'];
                        saveas(gcf,FileName)
                        FileName = [SubFolder,'/MSD_r inst.fig'];
                        saveas(gcf,FileName)

                        %% Plot instantaneous MSD_t for both discrete models
                        figure(4);

                        Style = '-.';
                        [p5(1),~] = PlotCurve(time_bead,msd_bead_t0*LengthConversion^2,...
                            msd_bead_t0_err*LengthConversion^2,Color,Style);

                        Style = '-';
                        [p5(2),~] = PlotCurve(time_meso,msd_meso_t0*LengthConversion^2,...
                            msd_meso_t0_err*LengthConversion^2,Color,Style);

                        l = legend(p4,l_entries);
                        l.FontSize = FontSize/2;
                        l.Interpreter = 'latex';
                        l.Location = 'Best';

                        set(gca,'FontSize',FontSize/1.5)
                        set(gcf,'Color','w')
                        pbaspect([1 1 1])
                        xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
                        ylabel('MSD$_t(\tau)$ (nm$^2$)','FontSize',FontSize,'Interpreter','latex')
%                         title(['$\gamma$ = ',num2str(damp,'%.2e'),' [force $\cdot$ time/length], ',...
%                             num2str(damp*DamperConversion,'%.2e'),' [N s m$^{-1}$]'],...
%                             'FontSize',FontSize/3,'Interpreter','latex')
                        %     xlim([0 1])

                        set(gcf,'Position',[200 200 275 275])

%                         FileName = [ImgTag,'.MSD_t inst.png'];
%                         saveas(gcf,FileName)
%                         FileName = [ImgTag,'.MSD_t inst.fig'];
%                         saveas(gcf,FileName)
                        FileName = [SubFolder,'/MSD_t inst.png'];
                        saveas(gcf,FileName)
                        FileName = [SubFolder,'/MSD_t inst.fig'];
                        saveas(gcf,FileName)

                        %% Plot MSD_r for both discrete models
                        figure(5);

                        Style = '-.';
                        [p6(1),~] = PlotCurve(time_bead,msd_bead_r00*LengthConversion^2,...
                            msd_bead_r00_err*LengthConversion^2,Color,Style);

                        Style = '-';
                        [p6(2),~] = PlotCurve(time_meso,msd_meso_r00*LengthConversion^2,...
                            msd_meso_r00_err*LengthConversion^2,Color,Style);

                        l = legend(p6,l_entries);
                        l.FontSize = FontSize/2;
                        l.Interpreter = 'latex';
                        l.Location = 'Best';

                        set(gca,'FontSize',FontSize/1.5)
                        set(gcf,'Color','w')
                        pbaspect([1 1 1])
                        xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
                        ylabel('MSD$_r$ (nm$^2$)','FontSize',FontSize,'Interpreter','latex')

                        set(gcf,'Position',[200 200 275 275])
                        
                        npt = 30;
                        y = msd_bead_r00;
                        f = msd_meso_r00;
%                         RSS = sum((y-f).^2);
%                         TSS = sum((y-mean(y)).^2);
                        
                        R2_rng = round(linspace(1,length(msd_meso_r00),npt));
                        RSS = sum((y(R2_rng)-f(R2_rng)).^2);
                        TSS = sum((y(R2_rng)-mean(y(R2_rng))).^2);
                        R2_msdr(dm,nk) = 1-RSS/TSS;

%                         FileName = [ImgTag,'.MSD_r.png'];
%                         saveas(gcf,FileName)
%                         FileName = [ImgTag,'.MSD_r.fig'];
%                         saveas(gcf,FileName)
                        FileName = [SubFolder,'/MSD_r.png'];
                        saveas(gcf,FileName)
                        FileName = [SubFolder,'/MSD_r.fig'];
                        saveas(gcf,FileName)
                        %% Plot MSD_t for both discrete models
                        figure(6);

                        Style = '-.';
                        [p7(1),~] = PlotCurve(time_bead,msd_bead_t00*LengthConversion^2,...
                            msd_bead_t00_err*LengthConversion^2,Color,Style);

                        Style = '-';
                        [p7(2),~] = PlotCurve(time_meso,msd_meso_t00*LengthConversion^2,...
                            msd_meso_t00_err*LengthConversion^2,Color,Style);

                        l = legend(p7,l_entries);
                        l.FontSize = FontSize/2;
                        l.Interpreter = 'latex';
                        l.Location = 'Best';

                        set(gca,'FontSize',FontSize/1.5)
                        set(gcf,'Color','w')
                        pbaspect([1 1 1])
                        xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
                        ylabel('MSD$_t$ (nm$^2$)','FontSize',FontSize,'Interpreter','latex')

                        set(gcf,'Position',[200 200 275 275])

                        y = msd_bead_t00;
                        f = msd_meso_t00;
                        RSS = sum((y(R2_rng)-f(R2_rng)).^2);
                        TSS = sum((y(R2_rng)-mean(y(R2_rng))).^2);
                        R2_msdt(dm,nk) = 1-RSS/TSS;

%                         FileName = [ImgTag,'.MSD_t.png'];
%                         saveas(gcf,FileName)
%                         FileName = [ImgTag,'.MSD_t.fig'];
%                         saveas(gcf,FileName)
                        FileName = [SubFolder,'/MSD_t.png'];
                        saveas(gcf,FileName)
                        FileName = [SubFolder,'/MSD_t.fig'];
                        saveas(gcf,FileName)
                    end
                end
            end
            figure(1)

            l = legend(p1,LegendEntries);
            l.FontSize = FontSize/2;
            l.Interpreter = 'latex';
            l.Location = 'Southeast';

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([1 1 1])
            xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
            ylabel('MSD (nm$^2$)','FontSize',FontSize,'Interpreter','latex')
            title(['$\gamma$ = ',num2str(damp,'%.2e'),' [force $\cdot$ time/length], ',...
                num2str(damp*DamperConversion,'%.2e'),' [N s m$^{-1}$]'],...
                'FontSize',FontSize/2,'Interpreter','latex')
            %     xlim([0 1])

            FileName = [Folder,'/Compare MSD vs Time.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/Compare MSD vs Time.fig'];
            saveas(gcf,FileName)

            figure(2)
            l = legend(p3,LegendEntries);
            l.FontSize = FontSize/2;
            l.Interpreter = 'latex';
            l.Location = 'Southeast';

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([1 1 1])
            xlabel('$t/\tau_r$','FontSize',FontSize,'Interpreter','latex')
            ylabel('MSD/$b^2$','FontSize',FontSize,'Interpreter','latex')
            title(['$\gamma$ = ',num2str(damp,'%.2e'),' [force $\cdot$ time/length], ',...
                num2str(damp*DamperConversion,'%.2e'),' [N s m$^{-1}$]'],...
                'FontSize',FontSize/2,'Interpreter','latex')

            FileName = [Folder,'/Compare Normalized MSD vs Time.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/Compare Normalized MSD vs Time.fig'];
            saveas(gcf,FileName)
        end
    end
    writematrix(R2_beadvsmeso,R2_beadvsmeso_filename)
    writematrix(R2_mesovsrouse,R2_mesovsrouse_filename)
    writematrix(R2_msdr,R2_msdr_filename)
    writematrix(R2_msdt,R2_msdt_filename)
else
    R2_beadvsmeso = readmatrix(R2_beadvsmeso_filename);
    R2_mesovsrouse = readmatrix(R2_mesovsrouse_filename);
    R2_msdr = readmatrix(R2_msdr_filename);
    R2_msdt = readmatrix(R2_msdt_filename);
end

close(wb1)
close all

if NpOrig==5^3
    BeadSpringOrMeso = 0;
else
    BeadSpringOrMeso = 1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MeasureSSMSDs(Package,damps,N_Kuhns,Np_orig,Ns,T,Nb,ka_in,kd_in,f0,dts,b,...
    N,Separations,OverridePostProcess)

global MSD_ss_filename StartMSDMeasurePct...
    MSD_r0_filename MSD_t0_filename BeadSpringOrMeso...
    MSD_r00_filename MSD_t00_filename LengthConversion DamperConversion

StartMSDMeasurePct = 0.001; % Leave close to zero - this is only useful if 
                            % chains were initiated anomalously far from 
                            % equilibrium and they weren't as evidenced by
                            % end-to-end distributions and these results'
                            % independence of of StartMSDMeasurePcnt

TotalSims = size(Package,1);

MSD_ss = zeros(length(damps),length(N_Kuhns));
MSD_t0 = zeros(length(damps),length(N_Kuhns));
MSD_r0 = zeros(length(damps),length(N_Kuhns));
MSD_r00 = zeros(length(damps),length(N_Kuhns));
MSD_t00 = zeros(length(damps),length(N_Kuhns));
SimCt = 0;
for BeadSpringOrMeso=0:1
    if BeadSpringOrMeso==0
        Np=5^3;
    else
        Np=11^3;
    end
    EdgeFactor = 1;
    InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dts(1),damps(1),...
        N_Kuhns(1),b,Separations(1));
    DefineMSDFileNames;
    wb1 = waitbar(0,'Post processing for const. damper, and swept N...');

    if ~isfile(MSD_ss_filename) || ~isfile(MSD_r0_filename) ||...
        ~isfile(MSD_t0_filename) || ~isfile(MSD_r00_filename) ||...
        ~isfile(MSD_t00_filename) ||OverridePostProcess
        % MSD wrt time for each damper while sweeping N
        for dm = 1:length(damps)
            damp = damps(dm);
            dt = dts(dm);

            for nk = 1:length(N_Kuhns)
                N_Kuhn = N_Kuhns(nk);
                for sp = 1:length(Separations)
                    Separation = Separations(sp);

                    PackageTemp = Package(ismember(Package(:,2),Np_orig),:);
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

                        [dt,damp,~,~] = DefineTimeStep(b,LengthConversion,...
                            DamperConversion,BeadSpringOrMeso);

                        [~,msd_mean,~,msd_r0,~,msd_t0,~,msd_r00,~,msd_t00,~] = ...
                            ImportMSD(BeadSpringOrMeso,Np,N,ka_in,kd_in,...
                            f0,dt,damp,N_Kuhn,b,Separation);

                        %% Crop MSD Measurement
%                         MSDStartIndx = ceil(StartMSDMeasurePct/100*length(time));
%                         msd_mean = msd_mean(MSDStartIndx:end);
%                         msd_r0 = msd_r0(MSDStartIndx:end);
%                         msd_t0 = msd_t0(MSDStartIndx:end);

                        MSD_t0(dm,nk) = mean(msd_t0(end-100:end));
                        MSD_r0(dm,nk) = mean(msd_r0(end-100:end));
                        MSD_t00(dm,nk) = mean(msd_t00(end-100:end));
                        MSD_r00(dm,nk) = mean(msd_r00(end-100:end));
                        MSD_ss(dm,nk) = mean(msd_mean(end-100:end));
                    end
                end
            end
        end
        writematrix(MSD_ss,MSD_ss_filename)
        writematrix(MSD_r0,MSD_r0_filename)
        writematrix(MSD_t0,MSD_t0_filename)
        writematrix(MSD_r00,MSD_r00_filename)
        writematrix(MSD_t00,MSD_t00_filename)
    end
    close(wb1)
    close all
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GenerateMSDPlotsWithNSwept(OverridePostProcess,Package,...
    damps,N_Kuhns,Np_orig,Ns,T,Nb,ka_in,kd_in,f0,dts,b,N,Separations)

global R2_mesovsrouse_filename FontSize DamperConversion...
    StartMSDMeasurePct AddOn BeadSpringOrMeso tau0_filename taur_filename...
    LengthConversion

StartMSDMeasurePct = 0.001; % Leave close to zero - this is only useful if 
                            % chains were initiated anomalously far from 
                            % equilibrium and they weren't as evidenced by
                            % end-to-end distributions and these results'
                            % independence of of StartMSDMeasurePcnt

dt_orig = dts;  
for BeadSpringOrMeso=0:1
    if BeadSpringOrMeso==0
        Np=5^3;
    else
        Np=11^3;
    end
    % Callout Input Script
    EdgeFactor = 1; 
    InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dts(1),damps(1),...
        N_Kuhns(1),b,Separations(1));
    DefineMSDFileNames;

    taus = readmatrix(tau0_filename);
    tauR = readmatrix(taur_filename);

    R2_mesovsrouse = zeros(length(damps),length(N_Kuhns));

    TotalSims = size(Package,1);

    SimCt = 0;

    wb1 = waitbar(0,'Post processing for const. damper, and swept N...');

    if ~isfile(R2_mesovsrouse_filename) || OverridePostProcess
        % MSD wrt time for each damper while sweeping N
        for dm = 1:length(damps)
            damp = damps(dm);
            dt = dts(dm);

            D_temp = 1/damp;
            tau0 = b^2/D_temp;

            Folder = ['Output Plots/damp ',num2str(damp,'%.3e')];
            if ~isfolder(Folder)
                mkdir(Folder)
            end

            % Figures
            figure(1); clf; hold on % msd
            figure(2); clf; hold on % msd in tau0<t<taur
            figure(3); clf; hold on % radial and circumferential msd over delta t vs t
            figure(4); clf; hold on % cumulative radial and circumferential msd vs t

            ColorRange = (linspace(0,0.85,length(N_Kuhns)))';
%         Colors = [ColorRange flipud(ColorRange) flipud(ColorRange)]; %Blue-to-red
            Colors = [zeros(size(ColorRange)),...
                flipud(ColorRange),...
                zeros(size(ColorRange))]; %Green-to-black
            LegendEntries = [];
            p0 = [];
            for nk = 1:length(N_Kuhns)
                N_Kuhn = N_Kuhns(nk);
                Color = Colors(nk,:);

                LegendEntries{nk} = ['$N =$ ',num2str(N_Kuhn)];

                for sp = 1:length(Separations)
                    Separation = Separations(sp);

                    PackageTemp = Package(ismember(Package(:,2),Np_orig),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,3),Ns),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,4),T),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,5),Nb),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,6),ka_in),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,7),kd_in),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,8),f0),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,9),dt_orig),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,10),damp),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,12),N_Kuhn),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,13),b),:);
                    PackageTemp = PackageTemp(ismember(PackageTemp(:,14),Separation),:);

                    waitbar(SimCt*size(PackageTemp,1)/TotalSims,wb1,...
                        'Post processing for const. damper, and swept N...')

                    if ~isempty(PackageTemp)
                        SimCt = SimCt+1;

                        % Callout Input Script
                        EdgeFactor = 1;
                        InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dt,damp,...
                            N_Kuhn,b,Separation);
                        SetDirAndFileNames;
                        DefineCompiledFileNames;
                        DefineMSDFileNames; % Define filenames for outputs
                        
                        BSOM = BeadSpringOrMeso;
                        [dt,damp,~,~] = DefineTimeStep(b,LengthConversion,...
                            DamperConversion,BSOM);

                        [time,msd_mean,msd_err,...
                            msd_r0,msd_r0_err,...
                            msd_t0,msd_t0_err,...
                            msd_r00,msd_r00_err,...
                            msd_t00,msd_t00_err] = ...
                            ImportMSD(BSOM,Np,N,ka_in,kd_in,f0,dt,damp,...
                            N_Kuhn,b,Separation);

                        %% Crop MSD Measurement
                        MSDStartIndx = ceil(StartMSDMeasurePct/100*length(time));
                        time = time(MSDStartIndx:end);
                        msd_mean = msd_mean(MSDStartIndx:end);
                        msd_err = msd_err(MSDStartIndx:end);
                        msd_r0 = msd_r0(MSDStartIndx:end);
                        msd_r0_err = msd_r0_err(MSDStartIndx:end);
                        msd_t0 = msd_t0(MSDStartIndx:end);
                        msd_t0_err = msd_t0_err(MSDStartIndx:end);
                        msd_r00 = msd_r00(MSDStartIndx:end);
                        msd_r00_err = msd_r00_err(MSDStartIndx:end);
                        msd_t00 = msd_t00(MSDStartIndx:end);
                        msd_t00_err = msd_t00_err(MSDStartIndx:end);

                        time = time-time(1);

                        npts = 50;
                        plot_rng = round((linspace(1,length(time),npts))');
                        figure(1) % msd
%                         if BeadSpringOrMeso==0
%                             p0(nk) = PlotErrorBar(time(plot_rng),...
%                                 msd_mean(plot_rng),...
%                                 msd_err(plot_rng),Color);
%                         else
                            Style = '-';
                            [p0(nk),~] = PlotCurve(time(plot_rng)/tau0,...
                                msd_mean(plot_rng)/b^2,...
                                msd_err(plot_rng)/b^2,Color,Style);
                            plot(time/tau0,ones(size(time)),'k--')
%                         end

                        %% Plot the rouse model
                        endindx = find(time>=tauR(dm,nk),1,'first');
                        rng = (1:endindx)';
                        figure(2)
%                         if BeadSpringOrMeso==0
%                             [~] = PlotErrorBar(time(rng)/tauR(dm,nk),...
%                                 msd_mean(rng)/b^2,...
%                                  msd_err(rng)/b^2,Color);
%                         else
                            Style = '-';
                            [~,~] = PlotCurve(time(rng)/tauR(dm,nk),...
                                msd_mean(rng)/b^2,...
                                msd_err(rng)/b^2,Color,Style);
%                         end

                        x = (linspace(0,tauR(dm,nk),100))';
                        y = (b^2)*(x/taus(dm,nk)).^0.5/b^2;
                        p = plot(x/tauR(dm,nk),y);
                        p.Color = Color;
                        p.LineStyle = '--';
                        p.LineWidth = 1.5;

                        f = msd_mean(rng);
                        y = (b^2)*(time(rng)/taus(dm,nk)).^0.5;
                        RSS = sum((y-f).^2);
                        TSS = sum((y-mean(y)).^2);
                        R2_mesovsrouse(dm,nk) = 1-RSS/TSS;

                        figure(3) % smoothed, instantaneous radial vs circumferential MSD in delta t
                        window = 30;
%                         if BeadSpringOrMeso==0
%                             [~] = PlotErrorBar(time(plot_rng,1),...
%                                 movmean(msd_r0(plot_rng),window),...
%                                 movmean(msd_r0_err(plot_rng),window),Color);
% 
%                             [~] = PlotErrorBar(time(plot_rng,1),...
%                                 movmean(msd_t0(plot_rng),window),...
%                                 movmean(msd_t0_err(plot_rng),window),Color);
%                         else
                            Style = '-';
                            [~,~] = PlotCurve(time(:,1)/tau0,...
                                movmean(msd_r0,window),...
                                movmean(msd_r0_err,window),Color,Style);
                            Style = '--';
                            [~,~] = PlotCurve(time(:,1)/tau0,...
                                movmean(msd_t0,window),...
                                movmean(msd_t0_err,window),Color,Style);
%                         end

                        figure(4) % radial vs circumferential MSD in delta t
%                         if BeadSpringOrMeso==0
%                             [~] = PlotErrorBar(time(plot_rng,1),...
%                                 movmean(msd_r0(plot_rng),window),...
%                                 movmean(msd_r0_err(plot_rng),window),Color);
% 
%                             [~] = PlotErrorBar(time(plot_rng,1),...
%                                 movmean(msd_t0(plot_rng),window),...
%                                 movmean(msd_t0_err(plot_rng),window),Color);
%                         else
                            Style = '-';
                            [~,~] = PlotCurve(time(plot_rng,1)/tau0,...
                                msd_r00(plot_rng),...
                                msd_r00_err(plot_rng),Color,Style);
                            Style = '--';
                            [~,~] = PlotCurve(time(plot_rng,1)/tau0,...
                                msd_t00(plot_rng),...
                                msd_t00_err(plot_rng),Color,Style);
%                         end
                    end
                end
            end
            figure(1)
            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([1 1 1])
            xlabel('$t$ ($\tau_0$)','FontSize',FontSize,'Interpreter','latex')
            ylabel('$\langle r^2(t) \rangle$ ($b^2$)','FontSize',FontSize,'Interpreter','latex')
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
            xlabel('$t$ ($\tau_r$)','FontSize',FontSize,'Interpreter','latex')
            ylabel('$\langle r^2(t) \rangle$ ($b^2$)','FontSize',FontSize,'Interpreter','latex')
            xlim([0 1])

            l = legend(p0,LegendEntries);
            l.FontSize = FontSize/2;
            l.Interpreter = 'latex';
            l.Location = 'Southeast';

            %         title(['$\gamma$ = ',num2str(damp,'%.2e'),' [force $\cdot$ time/length], ',...
            %             num2str(damp*DamperConversion,'%.3e'),' [N s m$^{-1}$]'],...
            %             'FontSize',FontSize/2,'Interpreter','latex')

            FileName = [Folder,'/',AddOn,'MSD Rouse inset vs t.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/',AddOn,'MSD Rouse inset vs t.fig'];
            saveas(gcf,FileName)

            figure(3)
            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([1 1 1])
            xlabel('$t$ ($\tau_0$)','FontSize',FontSize,'Interpreter','latex')
            ylabel('$\langle r^2(t) \rangle$ ($b^2$)','FontSize',FontSize,'Interpreter','latex')
            title(['$\gamma$ = ',num2str(damp,'%.2e'),' [force $\cdot$ time/length], ',...
                num2str(damp*DamperConversion,'%.3e'),' [N s m$^{-1}$]'],...
                'FontSize',FontSize/2,'Interpreter','latex')

            FileName = [Folder,'/',AddOn,'Instant Radial & Circ. MSD vs t.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/',AddOn,'Instant Radial & Circ. MSD vs t.fig'];
            saveas(gcf,FileName)

            figure(4)
            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([1 1 1])
            xlabel('$t$ ($\tau_0$)','FontSize',FontSize,'Interpreter','latex')
            ylabel('$\langle r^2(t) \rangle$ ($b^2$)','FontSize',FontSize,'Interpreter','latex')
            title(['$\gamma$ = ',num2str(damp,'%.2e'),' [force $\cdot$ time/length], ',...
                num2str(damp*DamperConversion,'%.3e'),' [N s m$^{-1}$]'],...
                'FontSize',FontSize/2,'Interpreter','latex')

            FileName = [Folder,'/',AddOn,'Radial & Circ. MSD vs t.png'];
            saveas(gcf,FileName)
            FileName = [Folder,'/',AddOn,'Radial & Circ. MSD vs t.fig'];
            saveas(gcf,FileName)
        end
        if BeadSpringOrMeso==1
            writematrix(R2_mesovsrouse,R2_mesovsrouse_filename)
        end
    end
    close(wb1)
    close all
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = PlotErrorBar(X,Y,Err,Color)

e = errorbar(X,Y,Err);
e.Marker = 'o';
e.Color = Color;
e.MarkerEdgeColor = Color;
e.MarkerFaceColor = Color;
e.LineStyle = 'none';
e.MarkerSize = 2;

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
function [R2,taus,tauR,rng] = FitRouseModel(time,msd_st_mean,N_Kuhn)  

global b

% global StartMSDMeasurePct
% 
% MSDStartIndx = round(StartMSDMeasurePct/100*length(time));
% time = time(MSDStartIndx:end);
% msd_st_mean = msd_st_mean(MSDStartIndx:end);
% time = time-time(1);

ARouse = N_Kuhn*b^2;
indx_temp = find(msd_st_mean>=ARouse,1,'first');
tauR = time(indx_temp);
if ~isempty(tauR)
    taus = tauR/N_Kuhn^2;

    % Initial guess based on data
    taus_nom = taus;

    stopindx = find(time>=tauR,1,'first');
else
    taus_nom = 3e-8;
    stopindx = [];
    taus = [];
end

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
        sig_taus = 1/8*taus_nom;
    else
        sig_taus = sig_taus*0.95;
    end
    taus_rng = [(linspace(taus_nom-sig_taus,taus_nom+sig_taus,NoPts))';taus_nom];
    taus_rng = unique(taus_rng); taus_rng(taus_rng<=0) = [];

    Perms = zeros(length(taus_rng),1);
    i = 0;
    for ai=1:length(taus_rng)
        i = i+1;
        Perms(i,1) = taus_rng(ai);
    end
    R2_all = zeros(size(Perms,1),1);
    taus_all = zeros(size(Perms,1),1);

    for i=1:size(Perms,1)
        taus_all(i) = Perms(i,1);

        x = time(rng);
        y = msd_st_mean(rng);

        ft = (b^2)*(x/taus_all(i)).^0.5;

        RSS = sum((y-ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_all(i) = 1-RSS/TSS;
    end
    diff = abs(R2_all-1);
    indx = find(diff==min(diff),1,'first');

    taus_nom = taus_all(indx);
    R2 = R2_all(indx);

    clf; hold on
    scatter(x,y,'k','filled')
    plot(x,(b^2)*(x/taus_nom).^0.5,'k--')

    if ~isempty(taus)
        plot(x,b^2*(x/taus).^0.5,'k')
    end

    xlabel('t')
    ylabel('MSD')
    if ct>50
        break;
    end
end
taus = taus_nom;
tauR = taus*N_Kuhn^2;
figure(100); close
close(wb5)
stopindx = find(time>=tauR,1,'first');
if isempty(stopindx)
    stopindx = length(time);
end
rng = (1:stopindx)';


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotHistograms(rx,ry,rz,~,time,stretch,...
    WRTStretchOrLength,HistoUnits)

global binMax binMin binWidth FontSize ConstantsFileName xlab...
    MovieName1 MovieFolder

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
%     bondtypes_temp = (bondtypes(i,:))';
    norms_temp = (norms(i,:))';
%     rx_temp = (rx(i,:))';
%     ry_temp = (ry(i,:))';
%     rz_temp = (rz(i,:))';
%     rx_temp(isnan(norms_temp)) = [];
%     ry_temp(isnan(norms_temp)) = [];
%     rz_temp(isnan(norms_temp)) = [];
    norms_temp(isnan(norms_temp)) = [];

    figure(201); clf; hold on
    XLim = sqrt(N_Kuhn);
    YLim = 1.5;
    color_t = 'k';
    SymDist = 0;
    PlotIdeal = 1;
    discrete_or_histo = 1;
    Plot1DHistogram(norms_temp,color_t,0.6,XLim,YLim,N_Kuhn,b,...
        SymDist,PlotIdeal,WRTStretchOrLength,HistoUnits,discrete_or_histo)

    set(gca,'FontSize',FontSize/1.5)
    xlabel(xlab,'FontSize',FontSize,'Interpreter','latex')
    ylabel('$p$','FontSize',FontSize,'Interpreter','latex')

    title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
        ', $t$ = ',num2str(time(i),'%.2e'),' s'],'FontSize',FontSize/2,...
        'Interpreter','latex')

    F1 = getframe(gcf);
    writeVideo(v1,F1);
    mov1(i) = F1;
end
close(wb3)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,theory] = Plot1DHistogram(Norms,Color,Alpha,~,YLim,N_Kuhn,b,...
    SymDist,PlotIdeal,WRTStretchOrLength,HistoUnits,discret_or_histo)

global FontSize LengthConversion %binMax binMin binWidth 

if WRTStretchOrLength==0
    binMax = 4;
    binMin = 0;
    binWidth = 0.08;%0.0396;
    Edges = (binMin:binWidth:binMax);
    XLim = 4;
else
    binMax = 4;
    binMin = 0;
    binWidth = 0.0396;
    XLim = binMax/2;
    Edges = (binMin:binWidth:binMax);
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
    h_temp = histogram(Lambda,'Normalization','pdf','BinEdges',Edges);
else
    if HistoUnits==0
        h_temp = histogram(Norms,'Normalization','pdf','BinEdges',Edges);
    else 
        h_temp = histogram(Norms,'Normalization','probability','BinEdges',Edges);
    end
end
if discret_or_histo==1
    h_temp.FaceColor = Color;
    h_temp.FaceAlpha = Alpha;
    h_temp.EdgeColor = 'k';
    h_temp.EdgeAlpha = 1;
    Bins = (h_temp.BinEdges)';
    model = h_temp;
else
    h_temp.FaceColor = 'none';
    h_temp.EdgeColor = 'none';
    Bins = (h_temp.BinEdges);
    yDiscrete = (h_temp.Values)';
    xDiscrete = (mean([Bins(1:end-1);Bins(2:end)]))';
    model = scatter(xDiscrete,yDiscrete);
    model.Marker = 'd';
    model.MarkerEdgeColor = 'k';
    model.MarkerFaceColor = 'k';
    model.SizeData = 10;
end

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
        ylim([0 0.1])
    end
end

set(gca,'FontSize',FontSize/1.5);
set(gcf,'Color','w')
pbaspect([3 1 1])

theory = [];
if PlotIdeal==1
    r = (linspace(0,10*N_Kuhn*b,10000))';
    lam = r/(sqrt(N_Kuhn)*b);
    if WRTStretchOrLength==0
        X = lam;
        xDiscrete = (Bins(1:end-1)+Bins(2:end))/2*sqrt(N_Kuhn)*b;
        yshift = 1;
    else
        X = r;
        xDiscrete = (Bins(1:end-1)+Bins(2:end))/2;
        yshift = 1.75;
        if HistoUnits==1
            yshift = 1/15;
        end
    end
    sigma = sqrt((N_Kuhn)/3)*b;
    P = 4*pi*(r.^2)*((sigma*sqrt(2*pi))^(-3)).*exp(-1/2*(r/sigma).^2);

    Area = trapz(lam,P);
    P = P/Area;
    
    P = P*yshift;
    theory = plot(X,P,'k-');
    theory.LineWidth = 1.5;
    
    yDiscrete = (h_temp.Values)';
    yIdeal = interp1(r,P,xDiscrete);
    yDiscrete(xDiscrete<=0) = [];
    yIdeal(xDiscrete<=0) = [];
    xDiscrete(xDiscrete<=0) = [];

    RSS = sum((yDiscrete-yIdeal).^2);
    TSS = sum((yIdeal-mean(yIdeal)).^2);
    R2 = 1-RSS/TSS;                %For Rouse diffusion

    text(2/3*sqrt(N_Kuhn),0.725,['$R^2$ = ',num2str(R2,'%.2f')],'FontSize',FontSize/2,'Interpreter','latex')
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
function [Corners,N_steps] = ImportDomainBoundaries(EdgeFactor,Np,N,...
    ka,kd,f0,dt,damp,N_Kuhn,b,Separation)

global RawDataFileName 

% Initializes non-sweeping parameters
InputScript(EdgeFactor,1,Np,N,ka,kd,f0,dt,damp,N_Kuhn,b,Separation);
SetDirAndFileNames;

RawData = load(RawDataFileName,'-mat');
Corners = RawData.Corners;
% FileName = ['Data/Compiled Outputs/Corners',FileTagCompiled,'.txt'];
% Corners = readmatrix(FileName);
N_steps = size(Corners,1);

end



%% Unused functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mean,SE] = CalcMeanAndSE(Data,Dim)

Mean = nanmean(Data,Dim);
SE = nanstd(Data,0,Dim)/sqrt(size(Data,Dim));

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

    Folder = ['Output Plots/damp ',num2str(damp,'%.2e')];
    if ~isfolder(Folder)
        mkdir(Folder)
    end
    FileName = [Folder,'/',AddOn,'R2 vs N.fig'];

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
function MakeSurfacePlotsCompare(damps,N_Kuhns,b)

global LengthConversion DamperConversion FontSize AddOn BeadSpringOrMeso...
    R2_beadvsrouse_filename R2_beadvsmeso_filename R2_mesovsrouse_filename...
    MSD_ss_filename R2_msdr_filename R2_msdt_filename
   
if isfile(R2_beadvsrouse_filename)
    R2_beadvsrouse = readmatrix(R2_beadvsrouse_filename);
end
if isfile(R2_beadvsmeso_filename)
    R2_beadvsmeso = readmatrix(R2_beadvsmeso_filename);
end
if isfile(R2_mesovsrouse_filename)
    R2_mesovsrouse = readmatrix(R2_mesovsrouse_filename);
end
if isfile(R2_msdr_filename)
    R2_msdr = readmatrix(R2_msdr_filename);
end
if isfile(R2_msdt_filename)
    R2_msdt = readmatrix(R2_msdt_filename);
end

BeadSpringOrMeso=0;
if isfile(MSD_ss_filename)
    MSD_ss_bead = readmatrix(MSD_ss_filename);
end
BeadSpringOrMeso=1;
DefineMSDFileNames;
if isfile(MSD_ss_filename)
    MSD_ss_meso= readmatrix(MSD_ss_filename);
end

if ~isempty(R2_beadvsrouse)
    % R2 bead vs meso
    figure(1)
    clf; hold on
    s = surf(damps*DamperConversion,N_Kuhns,R2_beadvsrouse');
    s.FaceColor = 'interp';
    s.EdgeAlpha = 0.5;
    set(gca,'XScale','log')
    set(gca,'FontSize',FontSize/1.5)

    xlim([min(damps) max(damps)]*DamperConversion)
    ylim([min(N_Kuhns) max(N_Kuhns)])

    c = colorbar;
    set(get(c,'label'),'string','$R^2$ (Bead-spring vs. Rouse)');
    set(get(c,'label'),'Interpreter','latex');
    set(get(c,'label'),'FontSize',FontSize);
    caxis([0 1])

    xlabel('$\gamma$ (kg s$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
    ylabel('$N$','FontSize',FontSize,'Interpreter','latex')
    title(['$b$ = ',num2str(b*LengthConversion*1e9,'%.2f'),' nm'],...
        'FontSize',FontSize/1.5','Interpreter','latex')
    pbaspect([1 1 1])

    set(gcf,'Color','w')

    FileName = 'Output Plots/R2 Bead-spring vs Rouse';
    saveas(gcf,[FileName,AddOn,'.png'])
    saveas(gcf,[FileName,AddOn,'.fig'])
end

if ~isempty(R2_beadvsmeso)
    % R2 bead vs meso
    figure(2)
    clf; hold on
    s = surf(damps*DamperConversion,N_Kuhns,R2_beadvsmeso');
    s.FaceColor = 'interp';
    s.EdgeAlpha = 0.5;
    set(gca,'XScale','log')
    set(gca,'FontSize',FontSize/1.5)

    xlim([min(damps) max(damps)]*DamperConversion)
    ylim([min(N_Kuhns) max(N_Kuhns)])

    c = colorbar;
    set(get(c,'label'),'string','$R^2$ (Bead-spring vs. Meso.)');
    set(get(c,'label'),'Interpreter','latex');
    set(get(c,'label'),'FontSize',FontSize);
    caxis([0 1])

    xlabel('$\gamma$ (kg s$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
    ylabel('$N$','FontSize',FontSize,'Interpreter','latex')
    title(['$b$ = ',num2str(b*LengthConversion*1e9,'%.2f'),' nm'],...
        'FontSize',FontSize/1.5','Interpreter','latex')
    pbaspect([1 1 1])

    set(gcf,'Color','w')

    FileName = 'Output Plots/R2 Bead-spring vs meso';
    saveas(gcf,[FileName,AddOn,'.png'])
    saveas(gcf,[FileName,AddOn,'.fig'])
end

if ~isempty(R2_mesovsrouse)
    % R2 bead vs meso
    figure(3)
    clf; hold on
    s = surf(damps*DamperConversion,N_Kuhns,R2_mesovsrouse');
    s.FaceColor = 'interp';
    s.EdgeAlpha = 0.5;
    set(gca,'XScale','log')
    set(gca,'FontSize',FontSize/1.5)

    xlim([min(damps) max(damps)]*DamperConversion)
    ylim([min(N_Kuhns) max(N_Kuhns)])

    c = colorbar;
    set(get(c,'label'),'string','$R^2$ (Meso. vs. Rouse)');
    set(get(c,'label'),'Interpreter','latex');
    set(get(c,'label'),'FontSize',FontSize);
    caxis([0 1])

    xlabel('$\gamma$ (kg s$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
    ylabel('$N$','FontSize',FontSize,'Interpreter','latex')
    title(['$b$ = ',num2str(b*LengthConversion*1e9,'%.2f'),' nm'],...
        'FontSize',FontSize/1.5','Interpreter','latex')
    pbaspect([1 1 1])

    set(gcf,'Color','w')

    FileName = 'Output Plots/R2 Meso vs rouse';
    saveas(gcf,[FileName,AddOn,'.png'])
    saveas(gcf,[FileName,AddOn,'.fig'])
end

if ~isempty(MSD_ss_meso) && ~isempty(MSD_ss_bead)
    dMSD = MSD_ss_meso-MSD_ss_bead;
    RMSD_MSD = mean(mean(dMSD.^2))^0.5;
    PcntDev = abs(dMSD)./MSD_ss_bead*100;

    % Relative deviation
    figure(4)
    clf; hold on
    s = surf(damps*DamperConversion,N_Kuhns,PcntDev');
    s.FaceColor = 'interp';
    s.EdgeAlpha = 0.5;
    set(gca,'XScale','log')
    set(gca,'FontSize',FontSize/1.5)

    xlim([min(damps) max(damps)]*DamperConversion)
    ylim([min(N_Kuhns) max(N_Kuhns)])

    c = colorbar;
    set(get(c,'label'),'string','$\%$ Dev. in Steady-state MSD');
    set(get(c,'label'),'Interpreter','latex');
    set(get(c,'label'),'FontSize',FontSize);
    caxis([0 50])

    xlabel('$\gamma$ (kg s$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
    ylabel('$N$','FontSize',FontSize,'Interpreter','latex')
    title(['$b$ = ',num2str(b*LengthConversion*1e9,'%.2f'),' nm'],...
        'FontSize',FontSize/1.5','Interpreter','latex')
    pbaspect([1 1 1])

    set(gcf,'Color','w')

    FileName = 'Output Plots/Pct Deviation in SS MSD Bead-spring vs meso';
    saveas(gcf,[FileName,AddOn,'.png'])
    saveas(gcf,[FileName,AddOn,'.fig'])
end

if ~isempty(R2_msdr) 
    figure(5)
    clf; hold on
    s = surf(damps*DamperConversion,N_Kuhns,R2_msdr');
    s.FaceColor = 'interp';
    s.EdgeAlpha = 0.5;
    set(gca,'XScale','log')
    set(gca,'FontSize',FontSize/1.5)

    xlim([min(damps) max(damps)]*DamperConversion)
    ylim([min(N_Kuhns) max(N_Kuhns)])

    c = colorbar;
    set(get(c,'label'),'string','$R^2$ of radial MSD');
    set(get(c,'label'),'Interpreter','latex');
    set(get(c,'label'),'FontSize',FontSize);
    caxis([0 1])

    xlabel('$\gamma$ (kg s$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
    ylabel('$N$','FontSize',FontSize,'Interpreter','latex')
    title(['$b$ = ',num2str(b*LengthConversion*1e9,'%.2f'),' nm'],...
        'FontSize',FontSize/1.5','Interpreter','latex')
    pbaspect([1 1 1])

    set(gcf,'Color','w')

    FileName = 'Output Plots/R2 MSD_r Bead-spring vs meso';
    saveas(gcf,[FileName,AddOn,'.png'])
    saveas(gcf,[FileName,AddOn,'.fig'])
end

if ~isempty(R2_msdt) 
    figure(6)
    clf; hold on
    s = surf(damps*DamperConversion,N_Kuhns,R2_msdt');
    s.FaceColor = 'interp';
    s.EdgeAlpha = 0.5;
    set(gca,'XScale','log')
    set(gca,'FontSize',FontSize/1.5)

    xlim([min(damps) max(damps)]*DamperConversion)
    ylim([min(N_Kuhns) max(N_Kuhns)])

    c = colorbar;
    set(get(c,'label'),'string','$R^2$ of tangential MSD');
    set(get(c,'label'),'Interpreter','latex');
    set(get(c,'label'),'FontSize',FontSize);
    caxis([0 1])

    xlabel('$\gamma$ (kg s$^{-1}$)','FontSize',FontSize,'Interpreter','latex')
    ylabel('$N$','FontSize',FontSize,'Interpreter','latex')
    title(['$b$ = ',num2str(b*LengthConversion*1e9,'%.2f'),' nm'],...
        'FontSize',FontSize/1.5','Interpreter','latex')
    pbaspect([1 1 1])

    set(gcf,'Color','w')

    FileName = 'Output Plots/R2 MSD_t Bead-spring vs meso';
    saveas(gcf,[FileName,AddOn,'.png'])
    saveas(gcf,[FileName,AddOn,'.fig'])
end

close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time,stretch,rx_ensemble,ry_ensemble,rz_ensemble,bondtypes_ensemble,...
    msd_st_mean,msd_st_err,...
    msd_st_r0_mean,msd_st_r0_err,msd_st_t0_mean,msd_st_t0_err,...
    msd_st_r00_mean,msd_st_r00_err,msd_st_t00_mean,msd_st_t00_err,...
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
    StartMSDMeasurePct BeadSpringOrMeso...
    MSDSt_r0FileName MSDSt_r0errFileName MSDSt_t0FileName MSDSt_t0errFileName...
    MSDSt_r00FileName MSDSt_r00errFileName MSDSt_t00FileName MSDSt_t00errFileName

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

    msd_st_r0_mean = zeros(N_steps,N_samples);
    msd_st_r0_err = zeros(N_steps,N_samples);
    msd_st_t0_mean = zeros(N_steps,N_samples);
    msd_st_t0_err = zeros(N_steps,N_samples);
    msd_st_r00_mean = zeros(N_steps,N_samples);
    msd_st_r00_err = zeros(N_steps,N_samples);
    msd_st_t00_mean = zeros(N_steps,N_samples);
    msd_st_t00_err = zeros(N_steps,N_samples);

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
        BondData = readmatrix(BondsFileName);

        %% Import Atoms Data
        AtomData = readmatrix(AtomsFileName);
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
            if BeadSpringOrMeso==0
                stickertype = 3; tethertype = 1;
            else
                stickertype = 2; tethertype = 1;
            end
            if timestep<=MSDStartStep
                msd_st_mean(i,n) = 0;
                msd_st_err(i,n) = 0;
                PrevStep = MSDStartStep;
            else
                [msd_st_mean(i,n),msd_st_err(i,n),...                
                msd_st_r0_mean(i,n),msd_st_r0_err(i,n),...
                msd_st_t0_mean(i,n),msd_st_t0_err(i,n),...
                msd_st_r00_mean(i,n),msd_st_r00_err(i,n),...
                msd_st_t00_mean(i,n),msd_st_t00_err(i,n),PrevStep] = ...
                    ComputeMSDs(timestep,AtomData,Corners(i,:),...
                    MSDStartStep,PrevStep,stickertype,tethertype);
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
    writematrix(rx_ensemble,RxFileName);
    writematrix(ry_ensemble,RyFileName);
    writematrix(rz_ensemble,RzFileName);

    writematrix(bondtypes_ensemble,BondTypesFileName);

    writematrix(msd_st_mean,MSDStFileName);
    writematrix(msd_st_err,MSDSt_errFileName);


    writematrix(msd_st_r0_mean,MSDSt_r0FileName);
    writematrix(msd_st_r0_err,MSDSt_r0errFileName);
    writematrix(msd_st_t0_mean,MSDSt_t0FileName);
    writematrix(msd_st_t0_err,MSDSt_t0errFileName);

    writematrix(msd_st_r00_mean,MSDSt_r00FileName);
    writematrix(msd_st_r00_err,MSDSt_r00errFileName);
    writematrix(msd_st_t00_mean,MSDSt_t00FileName);
    writematrix(msd_st_t00_err,MSDSt_t00errFileName);
    
    writematrix(g11_st_mean,g11_stFileName);
    writematrix(g22_st_mean,g22_stFileName);
    writematrix(g33_st_mean,g33_stFileName);
    writematrix(g12_st_mean,g12_stFileName);
    writematrix(g23_st_mean,g23_stFileName);
    writematrix(g31_st_mean,g31_stFileName);

    writematrix(g11_st_err,g11_st_errFileName);
    writematrix(g22_st_err,g22_st_errFileName);
    writematrix(g33_st_err,g33_st_errFileName);
    writematrix(g12_st_err,g12_st_errFileName);
    writematrix(g23_st_err,g23_st_errFileName);
    writematrix(g31_st_err,g31_st_errFileName);

    writematrix([time stretch],TimeStretchFileName)
else
    rx_ensemble=readmatrix(RxFileName);
    ry_ensemble=readmatrix(RyFileName);
    rz_ensemble=readmatrix(RzFileName);

    bondtypes_ensemble=readmatrix(BondTypesFileName);

    msd_st_mean=readmatrix(MSDStFileName);
    msd_st_err=readmatrix(MSDSt_errFileName);

    msd_st_r0_mean=readmatrix(MSDSt_r0FileName);
    msd_st_r0_err=readmatrix(MSDSt_r0errFileName);
    msd_st_t0_mean=readmatrix(MSDSt_t0FileName);
    msd_st_t0_err=readmatrix(MSDSt_t0errFileName);
    
    msd_st_r00_mean=readmatrix(MSDSt_r00FileName);
    msd_st_r00_err=readmatrix(MSDSt_r00errFileName);
    msd_st_t00_mean=readmatrix(MSDSt_t00FileName);
    msd_st_t00_err=readmatrix(MSDSt_t00errFileName);

    dat = readmatrix(TimeStretchFileName);

    g11_st_mean=readmatrix(g11_stFileName);
    g22_st_mean=readmatrix(g22_stFileName);
    g33_st_mean=readmatrix(g33_stFileName);
    g12_st_mean=readmatrix(g12_stFileName);
    g23_st_mean=readmatrix(g23_stFileName);
    g31_st_mean=readmatrix(g31_stFileName);

    g11_st_err=readmatrix(g11_st_errFileName);
    g22_st_err=readmatrix(g22_st_errFileName);
    g33_st_err=readmatrix(g33_st_errFileName);
    g12_st_err=readmatrix(g12_st_errFileName);
    g23_st_err=readmatrix(g23_st_errFileName);
    g31_st_err=readmatrix(g31_st_errFileName);

    time = dat(:,1); stretch = dat(:,2);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time,stretch,rx_ensemble,ry_ensemble,rz_ensemble,bondtypes_ensemble,...
    msd_st_mean,msd_st_err,...
    msd_st_r0_mean,msd_st_r0_err,msd_st_t0_mean,msd_st_t0_err,...
    msd_st_r00_mean,msd_st_r00_err,msd_st_t00_mean,msd_st_t00_err,...
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
    StartMSDMeasurePct BeadSpringOrMeso...
    MSDSt_r0FileName MSDSt_r0errFileName MSDSt_t0FileName MSDSt_t0errFileName...
    MSDSt_r00FileName MSDSt_r00errFileName MSDSt_t00FileName MSDSt_t00errFileName

% Allocate
dx0 = Corners(1,5)-Corners(1,2);
stretch = (Corners(:,5)-Corners(:,2))/dx0;
time = Corners(:,1)*dt;


if ~isfile(RxFileName) || ~isfile(RyFileName) ||...
        ~isfile(RzFileName) || Override==1

    % Allocate
%     msd_th_mean = zeros(N_steps,N_samples);
%     msd_th_err = zeros(N_steps,N_samples);
    msd_st_mean = zeros(N_steps,N_samples);
    msd_st_err = zeros(N_steps,N_samples);

    msd_st_r0_mean = zeros(N_steps,N_samples);
    msd_st_r0_err = zeros(N_steps,N_samples);
    msd_st_t0_mean = zeros(N_steps,N_samples);
    msd_st_t0_err = zeros(N_steps,N_samples);
    msd_st_r00_mean = zeros(N_steps,N_samples);
    msd_st_r00_err = zeros(N_steps,N_samples);
    msd_st_t00_mean = zeros(N_steps,N_samples);
    msd_st_t00_err = zeros(N_steps,N_samples);

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
        BondData = readmatrix(BondsFileName);

        %% Import Atoms Data
        AtomData = readmatrix(AtomsFileName);
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
            atomtypes = AtomData(AtomData(:,1)==timestep,end-1);
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
            if timestep<=MSDStartStep
                msd_st_mean(i,n) = 0;
                msd_st_err(i,n) = 0;
                PrevStep = MSDStartStep;
            else
                if BeadSpringOrMeso==0
                    stickertype = 3;   %End groups (MSD should follow Rouse diffusion)
                    tethertype = 1;
                else
                    stickertype = 2;
                    tethertype = 1;
                end
                [msd_st_mean(i,n),msd_st_err(i,n),...
                    msd_st_r0_mean(i,n),msd_st_r0_err(i,n),...
                    msd_st_t0_mean(i,n),msd_st_t0_err(i,n),...
                    msd_st_r00_mean(i,n),msd_st_r00_err(i,n),...
                    msd_st_t00_mean(i,n),msd_st_t00_err(i,n),PrevStep] = ...
                    ComputeMSDs(timestep,AtomData,Corners(i,:),...
                    MSDStartStep,PrevStep,stickertype,tethertype);
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
    writematrix(rx_ensemble,RxFileName);
    writematrix(ry_ensemble,RyFileName);
    writematrix(rz_ensemble,RzFileName);

    writematrix(bondtypes_ensemble,BondTypesFileName);

    writematrix(msd_st_mean,MSDStFileName);
    writematrix(msd_st_err,MSDSt_errFileName);


    writematrix(msd_st_r0_mean,MSDSt_r0FileName);
    writematrix(msd_st_r0_err,MSDSt_r0errFileName);
    writematrix(msd_st_t0_mean,MSDSt_t0FileName);
    writematrix(msd_st_t0_err,MSDSt_t0errFileName);

    writematrix(msd_st_r00_mean,MSDSt_r00FileName);
    writematrix(msd_st_r00_err,MSDSt_r00errFileName);
    writematrix(msd_st_t00_mean,MSDSt_t00FileName);
    writematrix(msd_st_t00_err,MSDSt_t00errFileName);
    
    writematrix(g11_st_mean,g11_stFileName);
    writematrix(g22_st_mean,g22_stFileName);
    writematrix(g33_st_mean,g33_stFileName);
    writematrix(g12_st_mean,g12_stFileName);
    writematrix(g23_st_mean,g23_stFileName);
    writematrix(g31_st_mean,g31_stFileName);

    writematrix(g11_st_err,g11_st_errFileName);
    writematrix(g22_st_err,g22_st_errFileName);
    writematrix(g33_st_err,g33_st_errFileName);
    writematrix(g12_st_err,g12_st_errFileName);
    writematrix(g23_st_err,g23_st_errFileName);
    writematrix(g31_st_err,g31_st_errFileName);

    writematrix([time stretch],TimeStretchFileName)
else
    rx_ensemble=readmatrix(RxFileName);
    ry_ensemble=readmatrix(RyFileName);
    rz_ensemble=readmatrix(RzFileName);

    bondtypes_ensemble=readmatrix(BondTypesFileName);

    msd_st_mean=readmatrix(MSDStFileName);
    msd_st_err=readmatrix(MSDSt_errFileName);

    msd_st_r0_mean=readmatrix(MSDSt_r0FileName);
    msd_st_r0_err=readmatrix(MSDSt_r0errFileName);
    msd_st_t0_mean=readmatrix(MSDSt_t0FileName);
    msd_st_t0_err=readmatrix(MSDSt_t0errFileName);
    
    msd_st_r00_mean=readmatrix(MSDSt_r00FileName);
    msd_st_r00_err=readmatrix(MSDSt_r00errFileName);
    msd_st_t00_mean=readmatrix(MSDSt_t00FileName);
    msd_st_t00_err=readmatrix(MSDSt_t00errFileName);
    
    dat = readmatrix(TimeStretchFileName);

    g11_st_mean=readmatrix(g11_stFileName);
    g22_st_mean=readmatrix(g22_stFileName);
    g33_st_mean=readmatrix(g33_stFileName);
    g12_st_mean=readmatrix(g12_stFileName);
    g23_st_mean=readmatrix(g23_stFileName);
    g31_st_mean=readmatrix(g31_stFileName);

    g11_st_err=readmatrix(g11_st_errFileName);
    g22_st_err=readmatrix(g22_st_errFileName);
    g33_st_err=readmatrix(g33_st_errFileName);
    g12_st_err=readmatrix(g12_st_errFileName);
    g23_st_err=readmatrix(g23_st_errFileName);
    g31_st_err=readmatrix(g31_st_errFileName);

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
function [msd,msd_err,msd0_r,msd0_r_err,msd0_t,msd0_t_err,...
    msd00_r,msd00_r_err,msd00_t,msd00_t_err,PrevStep]...
    = ComputeMSDs(timestep,AtomData,Corners,MSDStartStep,PrevStep,...
    stickertype,tethertype)

if timestep==30 || timestep==35
    timestep;
end

Lx = Corners(5)-Corners(2);
Ly = Corners(6)-Corners(3);
Lz = Corners(7)-Corners(4);

% Sort AtomData by ID number
AtomData = sortrows(AtomData,1);

% Check that these are aligned
AtomData00 = AtomData(AtomData(:,1)==MSDStartStep,:);   %First step
AtomData0 = AtomData(AtomData(:,1)==PrevStep,:);   %Previous step
AtomDataC = AtomData(AtomData(:,1)==timestep,:);

% Sort AtomDatas by particle ID
AtomData00 = sortrows(AtomData00,2);
AtomData0 = sortrows(AtomData0,2);
AtomDataC = sortrows(AtomDataC,2);

Pos00 = AtomData00(:,3:5);
Pos = AtomDataC(:,3:5);

% Find tether particle of corresponding molecule
AtomData_st = AtomDataC(AtomDataC(:,end-1)==tethertype,:);
AtomData_th = AtomDataC(AtomDataC(:,end-1)==stickertype,:);
AtomData_st = sortrows(AtomData_st,10);
AtomData_th = sortrows(AtomData_th,10);
R = AtomData_st(:,3:5)-AtomData_th(:,3:5);

if sum(ismember(AtomData_st(:,end),AtomData_th(:,end),'rows'))~=size(AtomData_st,1)
    [AtomData_st(:,end),AtomData_th(:,end)];
end

AtomData0_st = AtomData0(AtomData0(:,end-1)==tethertype,:);
AtomData0_th = AtomData0(AtomData0(:,end-1)==stickertype,:);
AtomData0_st = sortrows(AtomData0_st,10);
AtomData0_th = sortrows(AtomData0_th,10);
R0 = AtomData0_st(:,3:5)-AtomData0_th(:,3:5);

if sum(ismember(AtomData0_st(:,end),AtomData0_th(:,end),'rows'))~=size(AtomData0_st,1)
    [AtomData0_st(:,end),AtomData0_th(:,end)];
end

AtomData00_st = AtomData00(AtomData00(:,end-1)==tethertype,:);
AtomData00_th = AtomData00(AtomData00(:,end-1)==stickertype,:);
AtomData00_st = sortrows(AtomData00_st,10);
AtomData00_th = sortrows(AtomData00_th,10);
R00 = AtomData00_st(:,3:5)-AtomData00_th(:,3:5);

if sum(ismember(AtomData00_st(:,end),AtomData00_th(:,end),'rows'))~=size(AtomData00_st,1)
    [AtomData00_st(:,end),AtomData00_th(:,end)];
end

%Correct for periodic bounds by moving R0 back, as needed
dr0_temp = R-R0;
R0(dr0_temp(:,1)>Lx/2,1) = R0(dr0_temp(:,1)>Lx/2,1)+Lx;
R0(dr0_temp(:,1)<-Lx/2,1) = R0(dr0_temp(:,1)<-Lx/2,1)-Lx;
R0(dr0_temp(:,2)>Lx/2,2) = R0(dr0_temp(:,2)>Lx/2,2)+Lx;
R0(dr0_temp(:,2)<-Lx/2,2) = R0(dr0_temp(:,2)<-Lx/2,2)-Lx;
R0(dr0_temp(:,3)>Lx/2,3) = R0(dr0_temp(:,3)>Lx/2,3)+Lx;
R0(dr0_temp(:,3)<-Lx/2,3) = R0(dr0_temp(:,3)<-Lx/2,3)-Lx;

theta0 = acos(dot(R0,R,2)./(vecnorm(R0,2,2).*vecnorm(R,2,2)));
dR0_r = vecnorm(R,2,2)-vecnorm(R0,2,2);
dR0_t = vecnorm(R0,2,2).*theta0;
dR0_r(isnan(dR0_r)) = [];
dR0_t(isnan(dR0_t)) = [];
msd0_r = mean(dR0_r.^2);
msd0_r_err = std(dR0_r.^2)/sqrt(size(dR0_r,1));
msd0_t = mean(dR0_t.^2);
msd0_t_err = std(dR0_t.^2)/sqrt(size(dR0_t,1));

if max(theta0)>pi
    theta0;
end

%Correct for periodic bounds by moving R00 back, as needed
dr00_temp = R-R00;
R00(dr00_temp(:,1)>Lx/2,1) = R00(dr00_temp(:,1)>Lx/2,1)+Lx;
R00(dr00_temp(:,1)<-Lx/2,1) = R00(dr00_temp(:,1)<-Lx/2,1)-Lx;
R00(dr00_temp(:,2)>Lx/2,2) = R00(dr00_temp(:,2)>Lx/2,2)+Lx;
R00(dr00_temp(:,2)<-Lx/2,2) = R00(dr00_temp(:,2)<-Lx/2,2)-Lx;
R00(dr00_temp(:,3)>Lx/2,3) = R00(dr00_temp(:,3)>Lx/2,3)+Lx;
R00(dr00_temp(:,3)<-Lx/2,3) = R00(dr00_temp(:,3)<-Lx/2,3)-Lx;

theta00 = acos(dot(R00,R,2)./(vecnorm(R00,2,2).*vecnorm(R,2,2)));
dR00_r = vecnorm(R,2,2)-vecnorm(R00,2,2);
dR00_t = vecnorm(R00,2,2).*theta00;
dR00_r(isnan(dR00_r)) = [];
dR00_t(isnan(dR00_t)) = [];
msd00_r = mean(dR00_r.^2);
msd00_r_err = std(dR00_r.^2)/sqrt(size(dR00_r,1));
msd00_t = mean(dR00_t.^2);
msd00_t_err = std(dR00_t.^2)/sqrt(size(dR00_t,1));

if max(theta00)>pi
    theta00;
end

% Eliminate partilces of wrong type
Pos00(AtomData00(:,end-1)~=stickertype,:) = [];
Pos(AtomDataC(:,end-1)~=stickertype,:) = [];

% N = size(Pos,1);
dr = Pos-Pos00;

% Maxdr = max(vecnorm(dr,2,2));
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

PrevStep = timestep;

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
