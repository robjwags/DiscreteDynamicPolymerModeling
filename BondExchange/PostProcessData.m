function PostProcessData(Package,Override,ToggleDynamics,...
    LC,DC,BSOM,ASD,CF,LOL,TAB)

% Computes and plots stress-strain for every iteration of simulation

global LineWidth TurnOnDynamics font_size NoFitAttempts...
    LengthConversion DamperConversion DataSize BeadSpringOrMeso...
    OutputFolder PlotScalingTheory tau0... 
    PlotBondKineticsAndLifetimes PlotBondKineticsError CurrentFolder...
    NormalizeRates plot_ct AppendScalingData LinearOrLangevin TraceABond

OutputFolder = 'Output Plots';
if ~isfolder(OutputFolder)
    mkdir(OutputFolder)
end

NoFitAttempts = 100;
LineWidth = 1.5;
font_size = 20;
DataSize = 30;
TurnOnDynamics = ToggleDynamics;
LengthConversion = LC;
DamperConversion = DC;
BeadSpringOrMeso = BSOM;
LinearOrLangevin = LOL;
PlotScalingTheory = 1;
CurrentFolder = CF;
AppendScalingData = ASD;
TraceABond = TAB;

PlotBondKineticsAndLifetimes = 1;
PlotBondKineticsError = 1;

NormalizeRates = 1;     % Set 1 to normalize rates by tau0^{-1} and times by tau0
plot_ct = 20;           % How often to show fitting of parameters

plot_wrtc = 1;
plot_wrtd = 0;
plot_wrtphi = 1;
plot_wrtlambda = 1;
plot_wrtcopen = 0;
PlotDynwrtTime = 0;

% End-to-end histo plotting options
MakeHistos = 1;                 % Set 1 to make end-to-end histos
OverrideCompare = 0;            % Set 1 to override comparison between bead-spring and mesoscale end-to-end
WRTStretchOrLength = 0;         % Set 0 to make histos wrt. lambda, and 1 for length
HistoUnits = 0;                 % If WRTStretchOrLength==1, set 0 to make length in arbitrary units and 1 in nm

% MSD plotting options
MakeMSDPlots = 1;
OverrideCompareMSDs = 0;

Sample = unique(Package(:,1));
Np = unique(Package(:,2));        %Number of molecules
Ns = unique(Package(:,3));        %Number of stickeres per tether site
ka_in = unique(Package(:,4));         %Temperature
kd = unique(Package(:,5));        %Contour length of chains
f0 = unique(Package(:,6));        %Activation energy of association
N_Kuhns = unique(Package(:,7));        %Activation energy of dissociation
b = unique(Package(:,8));
phis = unique(Package(:,10));
kbT = unique(Package(:,11));
N = Np*Ns;

%% Sweepign parameters are N_Kuhn and phi
for i=1:length(ka_in)
    ka = ka_in(i);

    [damp,D,tau0,dtFact] = DefineTimeScale(b,LengthConversion,DamperConversion,...
        BeadSpringOrMeso);
    DiffCoeffs = D;
    
    dt = tau0/dtFact;
    
    %% Trace a bond through time and make a movie of its bonding history
    if TraceABond
        eaStar = -log(ka*tau0);
        edStar = -log(kd*tau0);
        TraceBondThroughTime(Override,Package,...
        N_Kuhns,Np,Ns,kbT,phis,f0,b,N,eaStar,edStar,D,dt)
    end
    
    %% Assemble the Ensemble Data
    Override = 1;
    AssembleEnsembleDynamicsData(Override,Package,damp,...
        N_Kuhns,Np,Ns,kbT,phis,ka,kd,f0,b,N,PlotDynwrtTime);
    close all

    %% Plot Outputs
    
    % Plot histograms of various bond lifetimes
    fig_nom =0;
    no_figs = 5;
    Override = 0;
    PlotLifetimesHistograms(fig_nom,no_figs,N_Kuhns,D,b,phis,Override)
    
    % Plot bond kinetic rates with respect to crosslink concentration
    % open crosslink concentration, crosslink separation distance, and/or 
    % overall polymer packing fraction
    Override = 1;
    PlotKineticRates(N_Kuhns,b,phis,kd,DiffCoeffs,plot_wrtd,plot_wrtc,...
        plot_wrtphi,plot_wrtlambda,plot_wrtcopen,Override);

    
    %% Plot end-to-end distribution comparison
    PlotAllHistogramsCompare(MakeHistos,OverrideCompare,...
        WRTStretchOrLength,HistoUnits,Package,ka_in,N,tau0,D,dt)

    %% Plot MSD comparison
    PlotMSDsCompare(MakeMSDPlots,OverrideCompareMSDs,...
        Package,N_Kuhns,Np,Ns,ka_in,kd,f0,dt,b,N,phis,D,kbT,Sample)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotAllHistogramsCompare(MakeHistos,OverrideCompare,...
        WRTStretchOrLength,HistoUnits,Package,ka_in,N,tau0,D,dt)

global xlab font_size 

SimCt = 0;

wb1 = waitbar(0,'Plotting comparison histograms...');

DefineAssembledFileNames(ka_in);


for p=1:size(Package,1)
    Sample = Package(1,2);  
    Np = Package(p,2);  
    Ns = Package(p,3);      
    ka_in = Package(p,4);       
    kd_in = Package(p,5);      
    N_Kuhn = Package(p,7);   
    b = Package(p,8);   
    dist = Package(p,9);
    phi_fn = Package(p,10);
    kbT = Package(p,11);
           
    phi = round(pi*N_Kuhn/6.*(b./dist).^3,4);
    
    eaStar = -log(ka_in*tau0);
    edStar = -log(kd_in*tau0);

    waitbar(p/size(Package,1),wb1,...
        'Plotting comparison histograms...')

    %Define movie names
    Folder = 'Output Plots/End-to-End';
    if ~isfolder(Folder)
        mkdir(Folder)
    end
    FigureName = [Folder,'/N ',num2str(N_Kuhn),'.phi ',num2str(phi,'%.2f')];

    if MakeHistos==1 && (~isfile([FigureName,'_stretch.png']) ||...
            OverrideCompare)

        SimCt = SimCt+1;

        % Callout Input Script
        InputScript(Sample,Np,Ns,N,N_Kuhn,b,phi_fn,eaStar,edStar,D,dt,kbT);
        SetDirAndFileNames;
        DefineCompiledFileNames;
%         Table = readtable(ConstantsFileName);
%         constants = table2array(Table);
%         N_Kuhn = constants(3); b = constants(4);

        %% For Bead-spring
        % Set Np properly for Bead Model
        BSOM = 0;   %Toggle bead-spring or meso

        [rx_bead,ry_bead,rz_bead] = ImportEndToEnd(BSOM,Sample,...
            Np,Ns,N,N_Kuhn,b,phi_fn,eaStar,edStar,D,dt,kbT);

        % Crop initially stretched chains
        rx_bead(1:20,:) = [];
        ry_bead(1:20,:) = [];
        rz_bead(1:20,:) = [];

        norms_bead = vecnorm([rx_bead(:) ry_bead(:) rz_bead(:)],2,2);
        norms_bead(isnan(norms_bead)) = [];

        %% For Bead-spring
        % Set Np properly for Bead Model
        BSOM = 1;       %Toggle bead-spring or meso

        [rx_meso,ry_meso,rz_meso] = ImportEndToEnd(BSOM,Sample,...
            Np,Ns,N,N_Kuhn,b,phi_fn,eaStar,edStar,D,dt,kbT);

        % Crop initially stretched chains
        rx_meso(1:20,:) = [];
        ry_meso(1:20,:) = [];
        rz_meso(1:20,:) = [];

        norms_meso = vecnorm([rx_meso(:) ry_meso(:) rz_meso(:)],2,2);
        norms_meso(isnan(norms_meso)) = [];

        %% Plot histogram with respect to stretch
        figure(1); clf; hold on
        AddOn = '_stretch';
        xlab = '$\lambda$';
        XLim = sqrt(N_Kuhn);
        YLim = 1.5;
        SymDist = 0;
        PlotIdeal = 0;

        % Mesoscale
        color = 'c';
        Plot1DHistogram(norms_meso,color,1,XLim,YLim,N_Kuhn,b,...
            SymDist,PlotIdeal,WRTStretchOrLength,HistoUnits)

        % Bead-spring
        color = 'r';
        PlotIdeal = 1;
        Plot1DHistogram(norms_bead,color,0.6,XLim,YLim,N_Kuhn,b,...
            SymDist,PlotIdeal,WRTStretchOrLength,HistoUnits)

        l = legend('Mesoscale Model','Bead-spring Model','Gaussian PDF');
        l.FontSize = font_size/2;
        l.Interpreter = 'latex';
        l.Location = 'Northeast';

        set(gca,'FontSize',font_size/1.5)
        xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
        ylabel('$p$','FontSize',font_size,'Interpreter','latex')

        dstar = dist/N_Kuhn/b;
        cstar = (N_Kuhn*dstar)^(-3);
        
        title(['$N$ =',num2str(N_Kuhn),...
            ', $\bar{d}^*$ = ',num2str(dstar,'%.2f'),...
            ', $\bar{c}^*$ = ',num2str(cstar,'%.2e'),...
            ', $\phi$ = ',num2str(phi,'%.2f')],'FontSize',font_size/2,...
            'Interpreter','latex')

        FileName = [FigureName,AddOn,'.png'];
        saveas(gcf,FileName)
        FileName = [FigureName,AddOn,'.fig'];
        saveas(gcf,FileName)
    end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotMSDsCompare(MakeMSDPlots,OverrideCompareMSDs,...
        Package,N_Kuhns,Np,Ns,ka_in,kd_in,f0,dt,b,N,phis,D,kbT,Sample)

global ConstantsFileName tau0 font_size

TotalSims = size(Package,1);
SimCt = 0;

wb1 = waitbar(0,'Plotting comparison of MSDs...');

DefineAssembledFileNames(ka_in);

% Replace phis
phi_pckg = phis;
for nk = 1:size(Package,1)/length(phis)-1
    phi_pckg = cat(1,phi_pckg,phis); 
end
Package(:,15) = phi_pckg;

R2_msd = zeros(length(phis),length(N_Kuhns));
R2_msdr0 = zeros(length(phis),length(N_Kuhns));
R2_msdt0 = zeros(length(phis),length(N_Kuhns));
R2_msdr00 = zeros(length(phis),length(N_Kuhns));
R2_msdt00 = zeros(length(phis),length(N_Kuhns));

MSE_msd = zeros(length(phis),length(N_Kuhns));
MSE_msdr0 = zeros(length(phis),length(N_Kuhns));
MSE_msdt0 = zeros(length(phis),length(N_Kuhns));
MSE_msdr00 = zeros(length(phis),length(N_Kuhns));
MSE_msdt00 = zeros(length(phis),length(N_Kuhns));

ColorRange = (linspace(0,1,length(N_Kuhns)))';
Colors = [zeros(size(ColorRange)) flipud(ColorRange) zeros(size(ColorRange))];

eaStar = -log(ka_in*tau0);
edStar = -log(kd_in*tau0);

marker_type = 'o';
marker_edge_color = 'k';
marker_size = 4;

% MSD wrt time for each damper while sweeping N
%Define movie names
Folder = 'Output Plots/MSDs';
if ~isfolder(Folder)
    mkdir(Folder)
end
FigureName = [Folder,'/MSD - phi.',num2str(phis(1),'%.3f'),'.png'];
FileName = [Folder,'/R2_msd.txt'];
if ~isfile(FigureName) || ~isfile(FileName) || OverrideCompareMSDs==1
    for ph = 1:length(phis)
        phi = phis(ph);
        
        figure(1); clf; hold on
        figure(2); clf; hold on
        figure(3); clf; hold on
        figure(4); clf; hold on
        figure(5); clf; hold on
        figure(6); clf; hold on
        
        tit = ['$\phi = $',num2str(phi,'%.3f')];
        
        for nk = 1:length(N_Kuhns)
            N_Kuhn = N_Kuhns(nk);
            Color = Colors(nk,:);
            
            PackageTemp = Package(ismember(Package(:,3),Ns),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,4),ka_in),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,5),kd_in),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,6),f0),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,7),N_Kuhn),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,8),b),:);
            PackageTemp = PackageTemp(ismember(PackageTemp(:,10),phi),:);
           
            waitbar(SimCt*size(PackageTemp,1)/TotalSims,wb1,...
                'Plotting comparison of MSDs...')
            if ~isempty(PackageTemp) && MakeMSDPlots==1
                SimCt = SimCt+1;
                
                % Callout Input Script
                InputScript(Sample,Np,Ns,N,N_Kuhn,b,phi,eaStar,edStar,D,dt,kbT);
                SetDirAndFileNames;
                DefineCompiledFileNames;
                Table = readtable(ConstantsFileName);
                constants = table2array(Table);
                N_Kuhn = constants(3); b = constants(4);
                
                % For Bead-spring
                BSOM = 0;   %Toggle bead-spring or meso
                
                [time_bead,msd_bead_mean,msd_bead_err,...
                    msd_bead_r0,msd_bead_r0_err,...
                    msd_bead_t0,msd_bead_t0_err,...
                    msd_bead_r00,msd_bead_r00_err,...
                    msd_bead_t00,msd_bead_t00_err] = ...
                    ImportMSD(BSOM,Sample,Np,Ns,N,N_Kuhn,b,phi,eaStar,edStar,D,dt,kbT);
                
                % For Mesoscale
                BSOM = 1;   %Toggle bead-spring or meso
                
                [time_meso,msd_meso_mean,msd_meso_err,...
                    msd_meso_r0,msd_meso_r0_err,...
                    msd_meso_t0,msd_meso_t0_err,...
                    msd_meso_r00,msd_meso_r00_err,...
                    msd_meso_t00,msd_meso_t00_err] = ...
                    ImportMSD(BSOM,Sample,Np,Ns,N,N_Kuhn,b,phi,eaStar,edStar,D,dt,kbT);
                
                if length(time_meso)>=length(time_bead)
                    npts = 50;
                    plot_rng = (round(linspace(1,length(time_bead),npts)))';
                    n_fit_pts = length(time_bead);
                    R2_rng = (round(linspace(1,length(time_bead),n_fit_pts)))';
                else
                    npts = 50;
                    plot_rng = (round(linspace(1,length(time_meso),npts)))';
                    n_fit_pts = length(time_meso);
                    R2_rng = (round(linspace(1,length(time_meso),n_fit_pts)))';
                end
                plot_rng(end) = [];
                
                l_entries = {'Bead-spring','Mesoscale'};
                
                % Plot overall MSD
                figure(1)
                
                p1(1) = PlotErrorBar(time_bead(plot_rng)/tau0,msd_bead_mean(plot_rng)/b^2,...
                    msd_bead_err(plot_rng)/b^2,Color,marker_type,marker_edge_color,...
                    marker_size);
                
                Style = '-';
                [p1(2),~] = PlotCurve(time_meso(plot_rng)/tau0,msd_meso_mean(plot_rng)/b^2,...
                    msd_meso_err(plot_rng)/b^2,Color,Style);
                
                %             l = legend(p1,l_entries);
                %             l.FontSize = font_size/2;
                %             l.Interpreter = 'latex';
                %             l.Location = 'Best';
                
                xticks([0 1e3])
                ylim([0 100])
                
                set(gca,'FontSize',font_size/1.5)
                set(gcf,'Color','w')
                pbaspect([1 1 1])
                xlabel('$t$ ($\tau_0$)','FontSize',font_size,'Interpreter','latex')
                ylabel('MSD ($b^2$)','FontSize',font_size,'Interpreter','latex')
                title(tit,'FontSize',font_size/2,'Interpreter','latex')
                
                set(gcf,'Position',[200 200 275 275])
                
                y = msd_bead_mean;
                f = msd_meso_mean;
                RSS = sum((y(R2_rng)-f(R2_rng)).^2);
                TSS = sum((y(R2_rng)-mean(y(R2_rng))).^2);
                MSE_msd(ph,nk) = mean((y(R2_rng)-f(R2_rng)).^2);
                R2_msd(ph,nk) = 1-RSS/TSS;
                
                % Plot incremental radial MSD
                figure(2)
                
                p2(1) = PlotErrorBar(time_bead(plot_rng)/tau0,msd_bead_r0(plot_rng)/b^2,...
                    msd_bead_r0_err(plot_rng)/b^2,Color,marker_type,marker_edge_color,...
                    marker_size);
                
                Style = '-';
                [p2(2),~] = PlotCurve(time_meso(plot_rng)/tau0,msd_meso_r0(plot_rng)/b^2,...
                    msd_meso_r0_err(plot_rng)/b^2,Color,Style);
                %
                %             l = legend(p2,l_entries);
                %             l.FontSize = font_size/2;
                %             l.Interpreter = 'latex';
                %             l.Location = 'Best';
                
                xticks([0 1e3])
                ylim([0 0.3])
                
                set(gca,'FontSize',font_size/1.5)
                set(gcf,'Color','w')
                pbaspect([1 1 1])
                xlabel('$t$ ($\tau_0$)','FontSize',font_size,'Interpreter','latex')
                ylabel('MSD$^0_r(\tau)$ ($b^2$)','FontSize',font_size,'Interpreter','latex')
                title(tit,'FontSize',font_size/2,'Interpreter','latex')
                
                set(gcf,'Position',[200 200 275 275])
                
                y = msd_bead_r0;
                f = msd_meso_r0;
                RSS = sum((y(R2_rng)-f(R2_rng)).^2);
                TSS = sum((y(R2_rng)-mean(y(R2_rng))).^2);
                MSE_msdr0(ph,nk) = mean((y(R2_rng)-f(R2_rng)).^2);
                R2_msdr0(ph,nk) = 1-RSS/TSS;
                
                % Plot incremental tangential MSD
                figure(3)
                
                p3(1) = PlotErrorBar(time_bead(plot_rng)/tau0,msd_bead_t0(plot_rng)/b^2,...
                    msd_bead_t0_err(plot_rng)/b^2,Color,marker_type,marker_edge_color,...
                    marker_size);
                
                Style = '-';
                [p3(2),~] = PlotCurve(time_meso(plot_rng)/tau0,msd_meso_t0(plot_rng)/b^2,...
                    msd_meso_t0_err(plot_rng)/b^2,Color,Style);
                
                l = legend(p3,l_entries);
                l.FontSize = font_size/2;
                l.Interpreter = 'latex';
                l.Location = 'Best';
                
                xticks([0 1e3])
                ylim([0 0.6])
                
                set(gca,'FontSize',font_size/1.5)
                set(gcf,'Color','w')
                pbaspect([1 1 1])
                xlabel('$t$ ($\tau_0$)','FontSize',font_size,'Interpreter','latex')
                ylabel('MSD$^0_t(\tau)$ ($b^2$)','FontSize',font_size,'Interpreter','latex')
                title(tit,'FontSize',font_size/2,'Interpreter','latex')
                
                set(gcf,'Position',[200 200 275 275])
                
                y = msd_bead_t0;
                f = msd_meso_t0;
                RSS = sum((y(R2_rng)-f(R2_rng)).^2);
                TSS = sum((y(R2_rng)-mean(y(R2_rng))).^2);
                MSE_msdt0(ph,nk) = mean((y(R2_rng)-f(R2_rng)).^2);
                R2_msdt0(ph,nk) = 1-RSS/TSS;
                
                % Plot overall radial MSD
                figure(4)
                
                p4(1) = PlotErrorBar(time_bead(plot_rng)/tau0,msd_bead_r00(plot_rng)/b^2,...
                    msd_bead_r00_err(plot_rng)/b^2,Color,marker_type,marker_edge_color,...
                    marker_size);
                
                Style = '-';
                [p4(2),~] = PlotCurve(time_meso(plot_rng)/tau0,msd_meso_r00(plot_rng)/b^2,...
                    msd_meso_r00_err(plot_rng)/b^2,Color,Style);
                
                %             l = legend(p4,l_entries);
                %             l.FontSize = font_size/2;
                %             l.Interpreter = 'latex';
                %             l.Location = 'Best';
                
                xticks([0 1e3])
                ylim([0 50])
                
                set(gca,'FontSize',font_size/1.5)
                set(gcf,'Color','w')
                pbaspect([1 1 1])
                xlabel('$t$ ($\tau_0$)','FontSize',font_size,'Interpreter','latex')
                ylabel('MSD$_r(\tau)$ ($b^2$)','FontSize',font_size,'Interpreter','latex')
                title(tit,'FontSize',font_size/2,'Interpreter','latex')
                
                set(gcf,'Position',[200 200 275 275])
                
                y = msd_bead_r00;
                f = msd_meso_r00;
                RSS = sum((y(R2_rng)-f(R2_rng)).^2);
                TSS = sum((y(R2_rng)-mean(y(R2_rng))).^2);
                MSE_msdr00(ph,nk) = mean((y(R2_rng)-f(R2_rng)).^2);
                R2_msdr00(ph,nk) = 1-RSS/TSS;
                
                % Plot overall tangential MSD
                figure(5)
                
                p5(1) = PlotErrorBar(time_bead(plot_rng)/tau0,msd_bead_t00(plot_rng)/b^2,...
                    msd_bead_t00_err(plot_rng)/b^2,Color,marker_type,marker_edge_color,...
                    marker_size);
                
                Style = '-';
                [p5(2),~] = PlotCurve(time_meso(plot_rng)/tau0,msd_meso_t00(plot_rng)/b^2,...
                    msd_meso_t00_err(plot_rng)/b^2,Color,Style);
                
                %             l = legend(p5,l_entries);
                %             l.FontSize = font_size/2;
                %             l.Interpreter = 'latex';
                %             l.Location = 'Best';
                
                xticks([0 1e3])
                ylim([0 150])
                
                set(gca,'FontSize',font_size/1.5)
                set(gcf,'Color','w')
                pbaspect([1 1 1])
                xlabel('$t$ ($\tau_0$)','FontSize',font_size,'Interpreter','latex')
                ylabel('MSD$_t(\tau)$ ($b^2$)','FontSize',font_size,'Interpreter','latex')
                title(tit,'FontSize',font_size/2,'Interpreter','latex')
                
                set(gcf,'Position',[200 200 275 275])
                
                y = msd_bead_t00;
                f = msd_meso_t00;
                RSS = sum((y(R2_rng)-f(R2_rng)).^2);
                TSS = sum((y(R2_rng)-mean(y(R2_rng))).^2);
                MSE_msdt00(ph,nk) = mean((y(R2_rng)-f(R2_rng)).^2);
                R2_msdt00(ph,nk) = 1-RSS/TSS;
                
                %% FIT ROUSE TIME
                ARouse = N_Kuhn*b^2;
                indx_temp = find(msd_bead_mean>=ARouse,1,'first');
                tauR_bead = time_bead(indx_temp)/tau0;
                taus_bead = tauR_bead/N_Kuhn^2;
                
                tauR = tauR_bead;
                taus = taus_bead;
                
                % Plot overall MSD in Rouse timescale
                figure(6)
                
                npts = round(2e4*144/N_Kuhn^2);
                if length(time_meso)>=length(time_bead)
                    plot_rng = (round(linspace(1,length(time_bead),npts)))';
                else
                    plot_rng = (round(linspace(1,length(time_meso),npts)))';
                end
                plot_rng(end) = [];
                
                [~,~] = PlotCurve(time_bead(plot_rng)/tau0/tauR,msd_bead_mean(plot_rng)/b^2,...
                    msd_bead_err(plot_rng)/b^2,Color,Style);
                
                p = plot(time_bead(plot_rng)/tau0/tauR,(time_bead(plot_rng)/taus/tau0).^0.5);
                p.Color = Color;
                p.LineStyle = '--';
                p.LineWidth = 1.5;
                
                xlim([0 1])
                ylim([0 40])
                
                set(gca,'FontSize',font_size/1.5)
                set(gcf,'Color','w')
                pbaspect([1 1 1])
                xlabel('$t$ ($\tau_r$)','FontSize',font_size,'Interpreter','latex')
                ylabel('MSD ($b^2$)','FontSize',font_size,'Interpreter','latex')
                title(tit,'FontSize',font_size/2,'Interpreter','latex')
                
                set(gcf,'Position',[200 200 300 300])
                
            end
        end
        figure(1)
        file_name = [Folder,'/MSD - phi.',num2str(phi,'%.3f')];
        saveas(gcf,[file_name,'.png'])
        saveas(gcf,[file_name,'.fig'])
        
        figure(2)
        file_name = [Folder,'/MSD_r0 - phi.',num2str(phi,'%.3f')];
        saveas(gcf,[file_name,'.png'])
        saveas(gcf,[file_name,'.fig'])
        
        figure(3)
        file_name = [Folder,'/MSD_t0 - phi.',num2str(phi,'%.3f')];
        saveas(gcf,[file_name,'.png'])
        saveas(gcf,[file_name,'.fig'])
        
        figure(4)
        file_name = [Folder,'/MSD_r00 - phi.',num2str(phi,'%.3f')];
        saveas(gcf,[file_name,'.png'])
        saveas(gcf,[file_name,'.fig'])
        
        figure(5)
        file_name = [Folder,'/MSD_t00 - phi.',num2str(phi,'%.3f')];
        saveas(gcf,[file_name,'.png'])
        saveas(gcf,[file_name,'.fig'])
        
        figure(6)
        file_name = [Folder,'/Rouse MSD - phi.',num2str(phi,'%.3f')];
        saveas(gcf,[file_name,'.png'])
        saveas(gcf,[file_name,'.fig'])
    end
    dlmwrite([Folder,'/R2_msd.txt'],R2_msd)
    dlmwrite([Folder,'/R2_msdr0.txt'],R2_msdr0)
    dlmwrite([Folder,'/R2_msdr00.txt'],R2_msdr00)
    dlmwrite([Folder,'/R2_msdt0.txt'],R2_msdt0)
    dlmwrite([Folder,'/R2_msdt00.txt'],R2_msdt00)
    
    dlmwrite([Folder,'/MSE_msd.txt'],MSE_msd)
    dlmwrite([Folder,'/MSE_msdr0.txt'],MSE_msdr0)
    dlmwrite([Folder,'/MSE_msdr00.txt'],MSE_msdr00)
    dlmwrite([Folder,'/MSE_msdt0.txt'],MSE_msdt0)
    dlmwrite([Folder,'/MSE_msdt00.txt'],MSE_msdt00)
else
    R2_msd = dlmread([Folder,'/R2_msd.txt']);
    R2_msdr0 = dlmread([Folder,'/R2_msdr0.txt']);
    R2_msdr00 = dlmread([Folder,'/R2_msdr00.txt']);
    R2_msdt0 = dlmread([Folder,'/R2_msdt0.txt']);
    R2_msdt00 = dlmread([Folder,'/R2_msdt00.txt']);
    
    MSE_msd = dlmread([Folder,'/MSE_msd.txt']);
    MSE_msdr0 = dlmread([Folder,'/MSE_msdr0.txt']);
    MSE_msdr00 = dlmread([Folder,'/MSE_msdr00.txt']);
    MSE_msdt0 = dlmread([Folder,'/MSE_msdt0.txt']);
    MSE_msdt00 = dlmread([Folder,'/MSE_msdt00.txt']);
end
close(wb1)
close all


% MSE of overall MSD
figure(1)
clf; hold on
s = surf(phis,N_Kuhns,MSE_msd');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'YScale','log')
xtickangle(45)
set(gca,'FontSize',font_size/1.5)

xlim([min(phis) max(phis)])
xticks(phis(2:2:end))
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','MSE ($b^2$) for MSD');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',font_size/1.5);

xlabel('$\phi$','FontSize',font_size,'Interpreter','latex')
ylabel('$N$','FontSize',font_size,'Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = [Folder,'/MSE_msd'];
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

% MSE of overall radial MSD
figure(2)
clf; hold on
s = surf(phis,N_Kuhns,MSE_msdr00');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'YScale','log')
xtickangle(45)
set(gca,'FontSize',font_size/1.5)

xlim([min(phis) max(phis)])
xticks(phis(2:2:end))
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','MSE ($b^2$) for radial MSD');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',font_size/1.5);

xlabel('$\phi$','FontSize',font_size,'Interpreter','latex')
ylabel('$N$','FontSize',font_size,'Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = [Folder,'/MSE_msd r'];
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

% MSE of overall tangential MSD
figure(3)
clf; hold on
s = surf(phis,N_Kuhns,MSE_msdt00');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'YScale','log')
xtickangle(45)
set(gca,'FontSize',font_size/1.5)

xlim([min(phis) max(phis)])
xticks(phis(2:2:end))
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','MSE ($b^2$) for tangential MSD');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',font_size/1.5);

xlabel('$\phi$','FontSize',font_size,'Interpreter','latex')
ylabel('$N$','FontSize',font_size,'Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = [Folder,'/MSE_msd t'];
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

% MSE of instant radial MSD
figure(4)
clf; hold on
s = surf(phis,N_Kuhns,MSE_msdr0');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'YScale','log')
xtickangle(45)
set(gca,'FontSize',font_size/1.5)

xlim([min(phis) max(phis)])
xticks(phis(2:2:end))
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','MSE ($b^2$) for radial MSD($\Delta t$)');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',font_size/1.5);

xlabel('$\phi$','FontSize',font_size,'Interpreter','latex')
ylabel('$N$','FontSize',font_size,'Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = [Folder,'/MSE_msd inst. r'];
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

% MSE of instant tangential MSD
figure(5)
clf; hold on
s = surf(phis,N_Kuhns,MSE_msdt0');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'YScale','log')
xtickangle(45)
set(gca,'FontSize',font_size/1.5)

xlim([min(phis) max(phis)])
xticks(phis(2:2:end))
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','MSE ($b^2$) for tangential MSD($\Delta t$)');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',font_size/1.5);

xlabel('$\phi$','FontSize',font_size,'Interpreter','latex')
ylabel('$N$','FontSize',font_size,'Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = [Folder,'/MSE_msd inst. t'];
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

close all

% R2 of overall MSD
figure(1)
clf; hold on
s = surf(phis,N_Kuhns,R2_msd');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'YScale','log')
% set(gca,'ColorScale','default')
xtickangle(45)
set(gca,'FontSize',font_size/1.5)

xlim([min(phis) max(phis)])
xticks(phis(2:2:end))
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$R^2$ for MSD');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',font_size/1.5);
c.Limits = [0 1]; caxis([0 1])

xlabel('$\phi$','FontSize',font_size,'Interpreter','latex')
ylabel('$N$','FontSize',font_size,'Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = [Folder,'/R2_msd'];
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

% R2 of overall radial MSD
figure(2)
clf; hold on
s = surf(phis,N_Kuhns,R2_msdr00');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'YScale','log')
% set(gca,'ColorScale','default')
xtickangle(45)
set(gca,'FontSize',font_size/1.5)

xlim([min(phis) max(phis)])
xticks(phis(2:2:end))
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$R^2$ for radial MSD');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',font_size/1.5);
c.Limits = [0 1]; caxis([0 1])

xlabel('$\phi$','FontSize',font_size,'Interpreter','latex')
ylabel('$N$','FontSize',font_size,'Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = [Folder,'/R2_msd r'];
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

% R2 of overall tangential MSD
figure(3)
clf; hold on
s = surf(phis,N_Kuhns,R2_msdt00');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'YScale','log')
% set(gca,'ColorScale','default')
xtickangle(45)
set(gca,'FontSize',font_size/1.5)

xlim([min(phis) max(phis)])
xticks(phis(2:2:end))
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$R^2$ for tangential MSD');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',font_size/1.5);
c.Limits = [0 1]; caxis([0 1])

xlabel('$\phi$','FontSize',font_size,'Interpreter','latex')
ylabel('$N$','FontSize',font_size,'Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = [Folder,'/R2_msd t'];
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

% R2 of instant radial MSD
figure(4)
clf; hold on
s = surf(phis,N_Kuhns,R2_msdr0');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'YScale','log')
% set(gca,'ColorScale','default')
xtickangle(45)
set(gca,'FontSize',font_size/1.5)

xlim([min(phis) max(phis)])
xticks(phis(2:2:end))
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$R^2$ for radial MSD($\Delta t$)');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',font_size/1.5);
c.Limits = [0 1]; caxis([0 1])

xlabel('$\phi$','FontSize',font_size,'Interpreter','latex')
ylabel('$N$','FontSize',font_size,'Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = [Folder,'/R2_msd inst. r'];
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

% R2 of instant tangential MSD
figure(5)
clf; hold on
s = surf(phis,N_Kuhns,R2_msdt0');
s.FaceColor = 'interp';
s.EdgeAlpha = 0.5;
set(gca,'XScale','log')
set(gca,'YScale','log')
% set(gca,'ColorScale','default')
xtickangle(45)
set(gca,'FontSize',font_size/1.5)

xlim([min(phis) max(phis)])
xticks(phis(2:2:end))
ylim([min(N_Kuhns) max(N_Kuhns)])

c = colorbar;
set(get(c,'label'),'string','$R^2$ for tangential MSD($\Delta t$)');
set(get(c,'label'),'Interpreter','latex');
set(get(c,'label'),'FontSize',font_size/1.5);
c.Limits = [0 1]; caxis([0 1])

xlabel('$\phi$','FontSize',font_size,'Interpreter','latex')
ylabel('$N$','FontSize',font_size,'Interpreter','latex')
pbaspect([1 1 1])

set(gcf,'Color','w')

FileName = [Folder,'/R2_msd inst. t'];
saveas(gcf,[FileName,'.png'])
saveas(gcf,[FileName,'.fig'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = PlotErrorBar(x,y,err,color,marker_type,marker_edge_color,...
    marker_size)

e = errorbar(x,y,err);
e.Marker = marker_type;
e.Color = color;
e.MarkerEdgeColor = marker_edge_color;
e.MarkerFaceColor = color;
e.LineStyle = 'none';
e.MarkerSize = marker_size;

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
function [time,msd_mean,msd_err,...
    msd_r0,msd_r0_err,...
    msd_t0,msd_t0_err,...
    msd_r00,msd_r00_err,...
    msd_t00,msd_t00_err] = ...
    ImportMSD(BSOM,Sample,Np,Ns,N,N_Kuhn,b,phi,eaStar,edStar,D,dt,kbT)

global BeadSpringOrMeso TimeStretchDataFileName MSDDataFileName

BeadSpringOrMeso = BSOM;

% Callout Input Script
InputScript(Sample,Np,Ns,N,N_Kuhn,b,phi,eaStar,edStar,D,dt,kbT);
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
function Plot1DHistogram(Norms,Color,Alpha,~,YLim,N_Kuhn,b,...
    SymDist,PlotIdeal,WRTStretchOrLength,HistoUnits)

global font_size LengthConversion %binMax binMin binWidth 

if WRTStretchOrLength==0
    binMax = 10;
    binMin = 0;
    binWidth = 0.08;%0.0396;
    Edges = (binMin:binWidth:binMax);
    XLim = 4;
else
    binMax = 10;
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
    h = histogram(Lambda,'Normalization','pdf','BinEdges',Edges);
else
    if HistoUnits==0
        h = histogram(Norms,'Normalization','pdf','BinEdges',Edges);
    else 
        h = histogram(Norms,'Normalization','probability','BinEdges',Edges);
    end
end
h.FaceColor = Color;
h.FaceAlpha = Alpha;
h.EdgeColor = 'k';
h.EdgeAlpha = 1;
Bins = (h.BinEdges)';

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

set(gca,'FontSize',font_size/1.5);
set(gcf,'Color','w')
pbaspect([3 1 1])

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

    text(1/2*sqrt(N_Kuhn),0.75,['$R^2$ = ',num2str(R2,'%.2f')],'FontSize',font_size/2,'Interpreter','latex')
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
function [rx,ry,rz] = ImportEndToEnd(BSOM,Sample,Np,Ns,N,N_Kuhn,b,phi,...
    eaStar,edStar,D,dt,kbT)

global BeadSpringOrMeso EndToEndDataFileName

BeadSpringOrMeso = BSOM;

% Callout Input Script
InputScript(Sample,Np,Ns,N,N_Kuhn,b,phi,eaStar,edStar,D,dt,kbT);
SetDirAndFileNames;
DefineCompiledFileNames;

EndToEndBS = load(EndToEndDataFileName,'-mat');
rx = EndToEndBS.rx;
ry = EndToEndBS.ry;
rz = EndToEndBS.rz;

norms = (rx.^2 + ry.^2 + rz.^2).^0.5;

mean_stretch = sqrt(N_Kuhn)*b;

% Adjust for periodic boundaries (some chains are measured as longer than
% they really are)
rx(norms/mean_stretch>3) = NaN;
ry(norms/mean_stretch>3) = NaN;
rz(norms/mean_stretch>3) = NaN;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotKineticRates(N_Kuhns,b,phis,kd_in,DiffCoeffs,...
    plot_wrtd,plot_wrtc,plot_wrtphi,plot_wrtlambda,plot_wrtcopen,Override)

global conc_all_bead tau_a_ens_bead LengthConversion...
    PlotBondKineticsAndLifetimes

% define number of clusters for bar chart
norm_conc = conc_all_bead*LengthConversion^3;
norm_conc = norm_conc*b^3;
no_figs = 13;

if plot_wrtc && (~isfile('Output Plots/Wrt. c/fa vs. c.png') || Override==1)
    x = norm_conc; 
    xlab = '$\bar c b^3$'; 
    save_tag = 'c';
    fig_nom = 0; 
    xlimits = [4e-4 0.1]; 
    xtickers = logspace(-3,-1,3);
    xstyle = 'log'; 
    ystyle = 'default';
    if PlotBondKineticsAndLifetimes==1
        MakeDynamicsPlotsWrtX(x,xlab,fig_nom,no_figs,N_Kuhns,kd_in,save_tag,...
            xlimits,xtickers,xstyle,ystyle,DiffCoeffs,b);
        close all
    end
end

% Plot dynamics wrt. nominal separation distance
if plot_wrtd && (~isfile('Output Plots/Wrt. d/fa vs. d.png') || Override==1)
    nom_sep = zeros(size(tau_a_ens_bead));
    for nk=1:length(N_Kuhns) % Define the nominal separation (it is a function of N_kuhn)
        nom_sep(:,nk) = (conc_all_bead(:,nk).^(-1/3))/LengthConversion/(N_Kuhns(nk)*b);
    end
    x = nom_sep; 
    xlab = '$\bar d/(Nb)$'; 
    save_tag = 'd';
    fig_nom = 0; 
    xlimits = [0 0.8]; 
    xtickers = min(xlimits):0.1:max(xlimits);
    xstyle = 'default'; 
    ystyle = 'default';
    if PlotBondKineticsAndLifetimes==1
        MakeDynamicsPlotsWrtX(x,xlab,fig_nom,no_figs,N_Kuhns,kd_in,save_tag,...
            xlimits,xtickers,xstyle,ystyle,DiffCoeffs,b);
        close all
    end
end

if plot_wrtphi && (~isfile('Output Plots/Wrt. phi/fa vs. phi.png') ||...
        Override==1)
    x = [phis phis phis]; 
    xlab = '$\bar \phi$'; 
    save_tag = 'phi';
    fig_nom = 0; 
    xlimits = [1e-2 1]; 
    xtickers = [1e-2 1e-1 1];
    xstyle = 'log'; 
    ystyle = 'default';
    if PlotBondKineticsAndLifetimes==1
        MakeDynamicsPlotsWrtX(x,xlab,fig_nom,no_figs,N_Kuhns,kd_in,save_tag,...
            xlimits,xtickers,xstyle,ystyle,DiffCoeffs,b);
        close all
    end
end

% Plot dynamics wrt. nominal chain stretch
if plot_wrtlambda && (~isfile('Output Plots/Wrt. lambda/fa vs. d.png') || Override==1)
    nom_stretch = zeros(size(tau_a_ens_bead));
    for nk=1:length(N_Kuhns) % Define the nominal separation (it is a function of N_kuhn)
        nom_stretch(:,nk) = (conc_all_bead(:,nk).^(-1/3))/LengthConversion/((2*N_Kuhns(nk))^0.5*b);
    end
    x = nom_stretch; 
    xlab = '$\bar \lambda$'; 
    save_tag = 'lambda';
    fig_nom = 0; 
    xlimits = [0 2]; 
    xtickers = min(xlimits):0.5:max(xlimits);
    xstyle = 'default'; 
    ystyle = 'default';
    if PlotBondKineticsAndLifetimes==1
        MakeDynamicsPlotsWrtX(x,xlab,fig_nom,no_figs,N_Kuhns,kd_in,save_tag,...
            xlimits,xtickers,xstyle,ystyle,DiffCoeffs,b);
        close all
    end
end

if plot_wrtcopen && (~isfile('Output Plots/Wrt. c_open/fa vs. c_open.png')...
        || Override==1)
    x = norm_conc; 
    xlab = '$\bar c_o b^3$'; 
    save_tag = 'c_open';
    fig_nom = 0; 
    xlimits = [4e-4 0.1]; 
    xtickers = logspace(-3,-1,3);
    xstyle = 'log'; 
    ystyle = 'default';
    if PlotBondKineticsAndLifetimes==1
        MakeDynamicsPlotsWrtX(x,xlab,fig_nom,no_figs,N_Kuhns,kd_in,save_tag,...
            xlimits,xtickers,xstyle,ystyle,DiffCoeffs,b);
        close all
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotLifetimesHistograms(fig_nom,no_figs,N_Kuhns,D,b,phis,Override)

global tau_a_all_bead tau_d_all_bead...
    tau_exc_all_bead tau_rpt_all_bead tau_adj_all_bead...
    tau_a_all_meso tau_d_all_meso...
    tau_exc_all_meso tau_rpt_all_meso tau_adj_all_meso...
    tau0 N_Kuhn...
    LengthConversion...
    tau_a_med_bs tau_a_range_bs tau_a_sample_ct_bs...
    tau_a_med_ms tau_a_range_ms tau_a_sample_ct_ms...
    tau_d_med_bs tau_d_range_bs tau_d_sample_ct_bs...
    tau_d_med_ms tau_d_range_ms tau_d_sample_ct_ms...    
    tau_rpt_med_bs tau_rpt_range_bs tau_rpt_sample_ct_bs...
    tau_rpt_med_ms tau_rpt_range_ms tau_rpt_sample_ct_ms...    
    tau_exc_med_bs tau_exc_range_bs tau_exc_sample_ct_bs...
    tau_exc_med_ms tau_exc_range_ms tau_exc_sample_ct_ms...
    tau_adj_med_bs tau_adj_range_bs tau_adj_sample_ct_bs...
    tau_adj_med_ms tau_adj_range_ms tau_adj_sample_ct_ms...
    
NormLength = b*LengthConversion;
tau0 = (NormLength^2)/D;


% Plot sampling counts of each event type
folder_name = 'Output Plots/Event Frequencies';
if ~isfolder(folder_name)
    mkdir(folder_name);
end

if ~isfile([folder_name,'/repeat-to-exchange ratio.png']) || Override==1
    
    % Calculate the ensemble statistics
    y_bs = tau_a_all_bead;
    y_ms = tau_a_all_meso;
    [tau_a_med_bs,tau_a_range_bs,tau_a_sample_ct_bs,...
        tau_a_med_ms,tau_a_range_ms,tau_a_sample_ct_ms] = ...
        CalculateLifetimeStatistics(y_bs,y_ms);
    
    y_bs = tau_d_all_bead;
    y_ms = tau_d_all_meso;
    [tau_d_med_bs,tau_d_range_bs,tau_d_sample_ct_bs,...
        tau_d_med_ms,tau_d_range_ms,tau_d_sample_ct_ms] = ...
        CalculateLifetimeStatistics(y_bs,y_ms);
    
    y_bs = tau_rpt_all_bead;
    y_ms = tau_rpt_all_meso;
    [tau_rpt_med_bs,tau_rpt_range_bs,tau_rpt_sample_ct_bs,...
        tau_rpt_med_ms,tau_rpt_range_ms,tau_rpt_sample_ct_ms] = ...
        CalculateLifetimeStatistics(y_bs,y_ms);
    
    y_bs = tau_exc_all_bead;
    y_ms = tau_exc_all_meso;
    [tau_exc_med_bs,tau_exc_range_bs,tau_exc_sample_ct_bs,...
        tau_exc_med_ms,tau_exc_range_ms,tau_exc_sample_ct_ms] = ...
        CalculateLifetimeStatistics(y_bs,y_ms);
    
    y_bs = tau_adj_all_bead;
    y_ms = tau_adj_all_meso;
    [tau_adj_med_bs,tau_adj_range_bs,tau_adj_sample_ct_bs,...
        tau_adj_med_ms,tau_adj_range_ms,tau_adj_sample_ct_ms] = ...
        CalculateLifetimeStatistics(y_bs,y_ms);
    
    
    ratio_of_rpt_to_exchange_bs = tau_rpt_sample_ct_bs./tau_exc_sample_ct_bs;
    ratio_of_rpt_to_exchange_ms = tau_rpt_sample_ct_ms./tau_exc_sample_ct_ms;
    
    fraction_of_repeats_bs = tau_rpt_sample_ct_bs./tau_a_sample_ct_bs;
    fraction_of_repeats_ms = tau_rpt_sample_ct_ms./tau_a_sample_ct_ms;
    
    fraction_of_exchanges_bs = tau_exc_sample_ct_bs./tau_a_sample_ct_bs;
    fraction_of_exchanges_ms = tau_exc_sample_ct_ms./tau_a_sample_ct_ms;
    
    first_times_bs = tau_a_sample_ct_bs-tau_exc_sample_ct_bs-tau_rpt_sample_ct_bs;
    first_times_ms = tau_a_sample_ct_ms-tau_exc_sample_ct_ms-tau_rpt_sample_ct_ms;
    fraction_of_virgin_bs = first_times_bs./tau_a_sample_ct_bs;
    fraction_of_virgin_ms = first_times_ms./tau_a_sample_ct_ms;
    
    % x-limits for the following bar charts
    xlims = [0.5 length(phis)+0.5];
    
    %%%% RAW SAMPLING EVENT COUNTS
      
    % Plot raw number of total detachment events
    figno = 1;
    y_bs = tau_d_sample_ct_bs;
    y_ms = tau_d_sample_ct_ms;
    ylab = '$n_{d}$';
    ylims = [0 400];
    leg_pos = 'northwest';
    file_name = [folder_name,'/n_d'];
    PlotStatisticalCounts(figno,N_Kuhns,phis,y_bs,y_ms,ylab,ylims,xlims,leg_pos,...
        file_name);
         
    % Plot raw number of total attachment events
    figno = figno + 1;
    y_bs = tau_a_sample_ct_bs;
    y_ms = tau_a_sample_ct_ms;
    ylab = '$n_{a}$';
    ylims = [0 400];
    leg_pos = 'northwest';
    file_name = [folder_name,'/n_total'];
    PlotStatisticalCounts(figno,N_Kuhns,phis,y_bs,y_ms,ylab,ylims,xlims,leg_pos,...
        file_name);   
    
    % Plot raw number of total repeat attachment events
    figno = figno + 1;
    y_bs = tau_rpt_sample_ct_bs;
    y_ms = tau_rpt_sample_ct_ms;
    ylab = '$n_{rpt}$';
    ylims = [0 200];
    leg_pos = 'northwest';
    file_name = [folder_name,'/n_repeat'];
    PlotStatisticalCounts(figno,N_Kuhns,phis,y_bs,y_ms,ylab,ylims,xlims,leg_pos,...
        file_name);
    
    % Plot raw number of total exchange events
    figno = figno + 1;
    y_bs = tau_exc_sample_ct_bs;
    y_ms = tau_exc_sample_ct_ms;
    ylab = '$n_{exc}$';
    ylims = [0 200];
    leg_pos = 'northwest';
    file_name = [folder_name,'/n_exchange'];
    PlotStatisticalCounts(figno,N_Kuhns,phis,y_bs,y_ms,ylab,ylims,xlims,leg_pos,...
        file_name);    
    
    % Raw number of total adjusted bond lifetimes events is same as raw
    % number of exchange events by definition
      
    %%%% FRACTIONS OF ATTACHMENT EVENT TYPES
    % Plot fraction of attachment events that are repeats
    figno = figno + 1;
    y_bs = fraction_of_repeats_bs;
    y_ms = fraction_of_repeats_ms;
    ylab = '$n_{rpt}/n_{tot}$';
    ylims = [0 1];
    leg_pos = 'northwest';
    file_name = [folder_name,'/repeat-to-total ratio'];
    PlotStatisticalCounts(figno,N_Kuhns,phis,y_bs,y_ms,ylab,ylims,xlims,leg_pos,...
        file_name);
    
    % Plot fraction of attachment events that are exchanges
    figno = figno + 1;
    y_bs = fraction_of_exchanges_bs;
    y_ms = fraction_of_exchanges_ms;
    ylab = '$n_{exc}/n_{tot}$';
    ylims = [0 1];
    leg_pos = 'northwest';
    file_name = [folder_name,'/exchange-to-total ratio'];
    PlotStatisticalCounts(figno,N_Kuhns,phis,y_bs,y_ms,ylab,ylims,xlims,leg_pos,...
        file_name);
    
    % Plot fraction of attachment events that are first time attachments
    figno = figno + 1;
    y_bs = fraction_of_virgin_bs;
    y_ms = fraction_of_virgin_ms;
    ylab = '$n_{first}/n_{tot}$';
    ylims = [0 1];
    leg_pos = 'northeast';
    file_name = [folder_name,'/first-to-total ratio'];
    PlotStatisticalCounts(figno,N_Kuhns,phis,y_bs,y_ms,ylab,ylims,xlims,leg_pos,...
        file_name);
    
    
    %%%% RATIO OF REPEAT TO EXCHANGE MEASURES
    
    % Plot ratio of number of repeat to exchange events
    figno = figno + 1;
    y_bs = ratio_of_rpt_to_exchange_bs;
    y_ms = ratio_of_rpt_to_exchange_ms;
    ylab = '$n_{rpt}/n_{exc}$';
    ylims = [0 20];
    leg_pos = 'northwest';
    file_name = [folder_name,'/repeat-to-exchange ratio'];
    PlotStatisticalCounts(figno,N_Kuhns,phis,y_bs,y_ms,ylab,ylims,xlims,leg_pos,...
        file_name);
  
    
    %%%% MEDIAN BOND STATUS LIFETIMES    
    
    % Plot median time of attached bond
    figno = figno + 1;
    y_bs = tau_a_med_bs;
    y_ms = tau_a_med_ms;
    ylab = '$\tilde{\tau}_a/\tau_0$';
    ylims = [0 500];
    leg_pos = 'northeast';
    file_name = [folder_name,'/median tau a'];
    PlotStatisticalCounts(figno,N_Kuhns,phis,y_bs,y_ms,ylab,ylims,xlims,leg_pos,...
        file_name);
    
    % Plot median time of detached bond
    figno = figno + 1;
    y_bs = tau_d_med_bs;
    y_ms = tau_d_med_ms;
    ylab = '$\tilde{\tau}_d/\tau_0$';
    ylims = [0 110];
    leg_pos = 'northeast';
    file_name = [folder_name,'/median tau d'];
    PlotStatisticalCounts(figno,N_Kuhns,phis,y_bs,y_ms,ylab,ylims,xlims,leg_pos,...
        file_name);
    
    % Plot median time of detached bond before repeat
    figno = figno + 1;
    y_bs = tau_rpt_med_bs;
    y_ms = tau_rpt_med_ms;
    ylab = '$\tilde{\tau}_{rpt}/\tau_0$';
    ylims = [0 60];
    leg_pos = 'northeast';
    file_name = [folder_name,'/median tau rpt'];
    PlotStatisticalCounts(figno,N_Kuhns,phis,y_bs,y_ms,ylab,ylims,xlims,leg_pos,...
        file_name);
    
    
    % Plot median time of detached bond before exchange
    figno = figno + 1;
    y_bs = tau_exc_med_bs;
    y_ms = tau_exc_med_ms;
    ylab = '$\tilde{\tau}_{exc}/\tau_0$';
    ylims = [0 600];
    leg_pos = 'northeast';
    file_name = [folder_name,'/median tau exc'];
    PlotStatisticalCounts(figno,N_Kuhns,phis,y_bs,y_ms,ylab,ylims,xlims,leg_pos,...
        file_name);
    
    
    % Plot median time of detached bond before exchange
    figno = figno + 1;
    y_bs = tau_adj_med_bs;
    y_ms = tau_adj_med_ms;
    ylab = '$\tilde{\tau}_{adj}/\tau_0$';
    ylims = [0 1000];
    leg_pos = 'northeast';
    file_name = [folder_name,'/median tau exc'];
    PlotStatisticalCounts(figno,N_Kuhns,phis,y_bs,y_ms,ylab,ylims,xlims,leg_pos,...
        file_name);
    
    close all
end


output_dir = 'Output Plots/Bond Lifetimes';
if ~isfolder(output_dir)
    mkdir(output_dir)
end

if ~isfile([output_dir,'/tau_a_hist.png']) || Override==1
    % Plot histograms of various bond status lifetimes
    InitializeFigures(fig_nom,no_figs);
    
    n_clr = length(phis);
    colors_phi = [(linspace(0.15,1,n_clr))',...
        (linspace(0.15,0.75,n_clr))',...
        flipud((linspace(0,0.15,n_clr))')];
    n_models = 2;
    
    for nk=1:length(N_Kuhns)
        N_Kuhn = N_Kuhns(nk);
        
        % 1) Attached bond lifetime
        figno = fig_nom+1;
        figure(figno)
        
        y_bs = tau_a_all_bead;
        y_ms = tau_a_all_meso;
        xlab_hist = '$\tau_a/\tau_0$';
        bin_width = 150;
        xlims = [0 1500];
%         ylims = [0 112];
        ylims = [0 1];
        PlotLifetimeHistograms(n_models,N_Kuhns,xlab_hist,y_bs,y_ms,phis,nk,...
            colors_phi,bin_width,xlims,ylims)
        
        
        % 2) Ajusted bond lifetime (bond lifetimes including detached phases
        % prior to repeat attachments)
        figno = fig_nom+2;
        figure(figno)
        
        y_bs = tau_adj_all_bead;
        y_ms = tau_adj_all_meso;
        xlab_hist = '$\tau_{adj}/\tau_0$';
        bin_width = 150;
        xlims = [0 1500];
        %     ylims = [0 50];
%         ylims = [0 112];
        ylims = [0 1];
        PlotLifetimeHistograms(n_models,N_Kuhns,xlab_hist,y_bs,y_ms,phis,nk,...
            colors_phi,bin_width,xlims,ylims)
        
        
        
        % 3) Detached bond lifetime
        figno = fig_nom+3;
        figure(figno)
        
        y_bs = tau_d_all_bead;
        y_ms = tau_d_all_meso;
        xlab_hist = '$\tau_d/\tau_0$';
        bin_width = 50;
        %     xlims = [0 1000];
        xlims = [0 500];
%         ylims = [0 200];
        ylims = [0 1];
        PlotLifetimeHistograms(n_models,N_Kuhns,xlab_hist,y_bs,y_ms,phis,nk,...
            colors_phi,bin_width,xlims,ylims)
        
        
        % 4) Detached bond lifetime before repeat
        figno = fig_nom+4;
        figure(figno)
        
        y_bs = tau_rpt_all_bead;
        y_ms = tau_rpt_all_meso;
        xlab_hist = '$\tau_{rpt}/\tau_0$';
        bin_width = 50;
        %     xlims = [0 1000];
        xlims = [0 500];
%         ylims = [0 150];
        ylims = [0 1];
        PlotLifetimeHistograms(n_models,N_Kuhns,xlab_hist,y_bs,y_ms,phis,nk,...
            colors_phi,bin_width,xlims,ylims)
        
        
        % 5) Detached bond lifetime before exchange
        figno = fig_nom+5;
        figure(figno)
        
        y_bs = tau_exc_all_bead;
        y_ms = tau_exc_all_meso;
        xlab_hist = '$\tau_{exc}/\tau_0$';
        bin_width = 50;
        %     xlims = [0 1000];
        xlims = [0 500];
%         ylims = [0 71];
        ylims = [0 1];
        PlotLifetimeHistograms(n_models,N_Kuhns,xlab_hist,y_bs,y_ms,phis,nk,...
            colors_phi,bin_width,xlims,ylims)
    end
    
    XSize = 1080;
    YSize = 720;
    
    SaveTheBondLifetimeFigures(fig_nom+1,XSize,YSize,'tau_a_hist',output_dir);
    SaveTheBondLifetimeFigures(fig_nom+2,XSize,YSize,'tau_a_adj_hist',output_dir);
    SaveTheBondLifetimeFigures(fig_nom+3,XSize,YSize,'tau_d_hist',output_dir);
    SaveTheBondLifetimeFigures(fig_nom+4,XSize,YSize,'tau_d_rpt_hist',output_dir);
    SaveTheBondLifetimeFigures(fig_nom+5,XSize,YSize,'tau_d_exc_hist',output_dir);
    
    close all
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotStatisticalCounts(figno,N_Kuhns,phis,y_bs,y_ms,ylab,ylims,xlims,...
    leg_pos,file_name)

global font_size


colors_nk = [zeros(length(N_Kuhns),1),...
    flipud((linspace(0,1,length(N_Kuhns)))'),...
    zeros(length(N_Kuhns),1)];

figure(figno); clf; hold on

subplot(1,2,1); hold on
b = bar(y_bs);
legend_entries = cell(length(N_Kuhns),1);
for nk=1:length(N_Kuhns)
    b(nk).FaceColor = colors_nk(nk,:);
    N_Kuhn = N_Kuhns(nk);
    legend_entries{nk} = ['$N=$ ',num2str(N_Kuhn)];
end
l = legend(legend_entries);
l.FontSize = font_size/2;
l.Interpreter = 'latex';
l.Location = leg_pos;
xticks(1:length(phis))
xticklabels(num2str(phis,'%.2f'));
xtickangle(45)
set(gca,'FontSize',font_size/1.75)
xlabel('$\phi$','FontSize',font_size,'Interpreter','latex')
ylabel(ylab,'FontSize',font_size,'Interpreter','latex')
ylim(ylims)
xlim(xlims)
pbaspect([1.5 1 1])
title('Bead-spring model','FontSize',font_size,'Interpreter','latex')

subplot(1,2,2); hold on
b = bar(y_ms);
for nk=1:length(N_Kuhns)
    b(nk).FaceColor = colors_nk(nk,:);
end
xticks(1:length(phis))
xticklabels(num2str(phis,'%.2f'));
xtickangle(45)
set(gca,'FontSize',font_size/1.75)
xlabel('$\phi$','FontSize',font_size,'Interpreter','latex')
ylabel(ylab,'FontSize',font_size,'Interpreter','latex')
ylim(ylims)
xlim(xlims)
pbaspect([1.5 1 1])
title('Mesoscale model','FontSize',font_size,'Interpreter','latex')

set(gcf,'color','w')
set(gcf,'Position',[100 100 800 400])

saveas(gcf,[file_name,'.png'])
saveas(gcf,[file_name,'.fig'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [med_bs,rng_bs,count_bs,med_ms,rng_ms,count_ms] =...
    CalculateLifetimeStatistics(y_bs,y_ms)

y_bs(y_bs==0) = NaN;
y_ms(y_ms==0) = NaN;

med_bs = zeros(size(y_bs,2),size(y_bs,3));
rng_bs = zeros(size(y_bs,2),size(y_bs,3));
count_bs = zeros(size(y_bs,2),size(y_bs,3));

med_ms = zeros(size(y_ms,2),size(y_ms,3));
rng_ms = zeros(size(y_ms,2),size(y_ms,3));
count_ms = zeros(size(y_ms,2),size(y_ms,3));
for sp=1:size(y_bs,2)   
    for nk=1:size(y_bs,3)
        y_tmp = y_bs;
        med_bs(sp,nk) = nanmedian(y_tmp(:,sp,nk));
        rng_bs(sp,nk) = max(y_tmp(:,sp,nk))-min(y_tmp(:,sp,nk));
        count_bs(sp,nk) = sum(~isnan(y_tmp(:,sp,nk)));

        y_tmp = y_ms;
        med_ms(sp,nk) = nanmedian(y_tmp(:,sp,nk));
        rng_ms(sp,nk) = max(y_tmp(:,sp,nk))-min(y_tmp(:,sp,nk));
        count_ms(sp,nk) = sum(~isnan(y_tmp(:,sp,nk)));       
    end
end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeDynamicsPlotsWrtX(x,xlab,fig_nom,no_figs,N_Kuhns,kd_in,...
    save_tag,xlimits,xtickers,xstyle,ystyle,D,b)

global font_size...
    ka_ens_meso kd_ens_meso kaSE_ens_meso kdSE_ens_meso...
    ka1_ens_meso ka1SE_ens_meso...
    fa_ens_meso fd_ens_meso faSE_ens_meso fdSE_ens_meso...
    kexc_ens_meso kexcSE_ens_meso ...
    krpt_ens_meso krptSE_ens_meso...
    ka_ens_bead kd_ens_bead kaSE_ens_bead kdSE_ens_bead...
    ka1_ens_bead ka1SE_ens_bead...
    fa_ens_bead fd_ens_bead faSE_ens_bead fdSE_ens_bead...
    kexc_ens_bead kexcSE_ens_bead ...
    krpt_ens_bead krptSE_ens_bead...
    tau0 N_Kuhn...
    NormalizeRates LengthConversion AppendScalingData...
    get_ka_from_adding...
    tau_a_ens_meso tau_aSE_ens_meso...
    tau_d_ens_meso tau_dSE_ens_meso...
    tau_exc_ens_meso tau_excSE_ens_meso ...
    tau_rpt_ens_meso tau_rptSE_ens_meso ...
    tau_adj_ens_meso tau_adjSE_ens_meso ...
    tau_rnm_ens_meso tau_rnmSE_ens_meso...
    tau_a_ens_bead tau_aSE_ens_bead...
    tau_d_ens_bead tau_dSE_ens_bead...
    tau_exc_ens_bead tau_excSE_ens_bead ...
    tau_rpt_ens_bead tau_rptSE_ens_bead ...
    tau_adj_ens_bead tau_adjSE_ens_bead ...
    tau_rnm_ens_bead tau_rnmSE_ens_bead...
    
% Remove the outlier for N=18, phi is second to lowest for bead-spring
x_out = 2; y_out = 2; mts = 0;
RemoveOutlier(x_out,y_out,mts);

get_ka_from_adding = 1;
if get_ka_from_adding==1
    ka_ens_bead = krpt_ens_bead +kexc_ens_bead;
    ka_ens_meso = krpt_ens_meso +kexc_ens_meso;
end

if AppendScalingData
    save_tag = [save_tag,'_with scaling data'];
end

NormLength = b*LengthConversion;
tau0 = (NormLength^2)/D;

marker_bs = 'o';
marker_ms = '^';
marker_edge_color = 'k';
line_style_dat = 'none';

InitializeFigures(fig_nom,no_figs);
color_range = (linspace(0,1,length(N_Kuhns)))';
colors = [zeros(size(color_range)) flipud(color_range) zeros(size(color_range))];

if contains(save_tag,'c')
    ft_ka = '$k_a^{semi} (1-c_{min}/c)$';
    ft_ka1 = '$k_{a,1}^{semi} [1-(c_{min}/c)^{1/3}]$';
    ft_kexc = '$k_{exc}^{semi} [1-(c_{min}/c)^{1/3}]$';
    ft_krpt = '$k_{rpt}^{semi} [1-(c_{min}/c)^{1/3}]$';
elseif contains(save_tag,'d')
    ft_ka = '$k_a^{semi} [1-(d/d_{max})^3]$';
    ft_ka1 = '$k_{a,1}^{semi} (1-d/d_{max})$';
    ft_kexc = '$k_{exc}^{semi} (1-d/d_{max})$';
    ft_krpt = '$k_{rpt}^{semi} (1-d/d_{max})$';
elseif contains(save_tag,'phi')
    ft_ka = '$k_a^{semi} (1-\phi_{min}/\phi)$';
    ft_ka1 = '$k_{a,1}^{semi} [1-(\phi_{min}/\phi)^{1/3}]$';
    ft_kexc = '$k_{exc}^{semi} [1-(\phi_{min}/\phi)^{1/3}]$';
    ft_krpt = '$k_{rpt}^{semi} [1-(\phi_{min}/\phi)^{1/3}]$';
else
    ft_ka = '$k_a^{semi} (1-c_{min}/c_o)$';
end

n_param = 2;  % two fitted parameters for every model type
coeff_ka_bs = zeros(length(N_Kuhns),n_param);
coeff_ka_ms = zeros(length(N_Kuhns),n_param);
coeff_ka1_bs = zeros(length(N_Kuhns),n_param);
coeff_ka1_ms = zeros(length(N_Kuhns),n_param);
coeff_kexc_bs = zeros(length(N_Kuhns),n_param);
coeff_kexc_ms = zeros(length(N_Kuhns),n_param);
coeff_krpt_bs = zeros(length(N_Kuhns),n_param);
coeff_krpt_ms = zeros(length(N_Kuhns),n_param);
ci_ka_bs = zeros(length(N_Kuhns),n_param);
ci_ka_ms = zeros(length(N_Kuhns),n_param);
ci_ka1_bs = zeros(length(N_Kuhns),n_param);
ci_ka1_ms = zeros(length(N_Kuhns),n_param);
ci_kexc_bs = zeros(length(N_Kuhns),n_param);
ci_kexc_ms = zeros(length(N_Kuhns),n_param);
ci_krpt_bs = zeros(length(N_Kuhns),n_param);
ci_krpt_ms = zeros(length(N_Kuhns),n_param);
R2_ka_bs = zeros(length(N_Kuhns),1);
R2_ka_ms = zeros(length(N_Kuhns),1);
R2_ka1_bs = zeros(length(N_Kuhns),1);
R2_ka1_ms = zeros(length(N_Kuhns),1);
R2_kexc_bs = zeros(length(N_Kuhns),1);
R2_kexc_ms = zeros(length(N_Kuhns),1);
R2_krpt_bs = zeros(length(N_Kuhns),1);
R2_krpt_ms = zeros(length(N_Kuhns),1);

for nk=1:length(N_Kuhns)
    
    color = colors(nk,:);
    x_temp = x(:,nk);
%     x_plt = (linspace(x_temp(1),x_temp(end),100))';
    N_Kuhn = N_Kuhns(nk);
 
    if contains(save_tag,'lambda')
        angle = 0;
        x_plt = (linspace(x_temp(1),x_temp(end),100))';
    elseif contains(save_tag,'d')
        angle = 45;
        x_plt = (linspace(x_temp(1),x_temp(end),100))';
    elseif contains(save_tag,'c')
        angle = 0;
        x_plt = (logspace(-4,-1,100))';
        x_plt(x_plt<x_temp(1)) = [];
        x_plt(x_plt>x_temp(end)) = [];
    else
        angle = 0;
        x_plt = (logspace(-2,0,100))';
        x_plt(x_plt<x_temp(1)) = [];
        x_plt(x_plt>x_temp(end)) = [];
    end

    
    if contains(save_tag,'c_open')
        x_temp = x(:,nk);
        N_Kuhn = N_Kuhns(nk);
        
        % 1) Detached bond lifetime wrt. open sticker concentration
        fd = mean([fd_ens_meso(:,nk),fd_ens_bead(:,nk)],2);
        x_temp = fd.*x_temp;
    end
    
    [ft_ka_bs,gof_ka_bs,...
        ft_ka_ms,gof_ka_ms,...
        ft_ka1_bs,gof_ka1_bs,...
        ft_ka1_ms,gof_ka1_ms,...
        ft_kexc_bs,gof_kexc_bs,...
        ft_kexc_ms,gof_kexc_ms,...
        ft_krpt_bs,gof_krpt_bs,...
        ft_krpt_ms,gof_krpt_ms] = ...
        IntakeFittedParameters(x_temp,save_tag,nk,N_Kuhn);
    
    % 1) Detachment rate
    figno = 1;
    figure(figno)
    if NormalizeRates==1
        plot(xlimits,[kd_in kd_in]*tau0*1e3,'k--');
    else
        plot(xlimits,[kd_in kd_in]/1e3,'k--');
    end
    
    y = kd_ens_bead(:,nk)*1e3;
    err = kdSE_ens_bead(:,nk)*1e3;
%     if get_ka_from_adding==1
%         err = NaN*ones(size(y));
%     end
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
        marker_edge_color)
    
    y = kd_ens_meso(:,nk)*1e3;
    err = kdSE_ens_meso(:,nk)*1e3;
%     if get_ka_from_adding==1
%         err = NaN*ones(size(y));
%     end
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
        marker_edge_color)
    
    PolishTheFigure(xlab,xlimits,xtickers,xstyle,ystyle,angle);
    if NormalizeRates==1
        ylim([0 1.5])
        ylabel('$\bar{k}_d \tau_0$ $(10^{-3})$','FontSize',font_size,'Interpreter','latex')
    else
        ylim([0 400])
        ylabel('$\bar{k}_d$ (kHz)','FontSize',font_size,'Interpreter','latex')
    end

    % 2) Attachment rate
    figno = figno+1;
    figure(figno)    
    if ~isempty(ft_ka_bs)
        [coeff_ka_bs(nk,:),ci_ka_bs(nk,:),coeff_ka_ms(nk,:),ci_ka_ms(nk,:)] = ...
            PlotTheEmpericalFits(x_plt,ft_ka_bs,ft_ka_ms,color,'-','--',1e3);
        R2_ka_bs(nk,:) = gof_ka_bs.rsquare;
        R2_ka_ms(nk,:) = gof_ka_ms.rsquare;
%         if nk==1
%             text(x_temp(1),1.5,ft_ka,'FontSize',font_size/1.75,'Interpreter','latex')
%         end
    end
    
    y = ka_ens_bead(:,nk)*1e3;
    err = kaSE_ens_bead(:,nk)*1e3;
%     if get_ka_from_adding==1
%         err = NaN*ones(size(y));
%     end
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
        marker_edge_color)
    
    y = ka_ens_meso(:,nk)*1e3;
    err = kaSE_ens_meso(:,nk)*1e3;
%     if get_ka_from_adding==1
%         err = NaN*ones(size(y));
%     end
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
        marker_edge_color)
    
    PolishTheFigure(xlab,xlimits,xtickers,xstyle,ystyle,angle);
    if NormalizeRates==1
        ylim([0 15])
        ylabel('$\bar{k}_a \tau_0$ $(10^{-3})$','FontSize',font_size,'Interpreter','latex')
    else
        ylim([0 400])
        ylabel('$\bar{k}_a$ (kHz)','FontSize',font_size,'Interpreter','latex')
    end
   
    
%     % 3) First-time attachment rate
%     figno = figno+1;
%     figure(figno)
%     if ~isempty(ft_ka1_bs)
%         [coeff_ka1_bs(nk,:),ci_ka1_bs(nk,:),coeff_ka1_ms(nk,:),ci_ka1_ms(nk,:)] = ...
%             PlotTheEmpericalFits(x_plt,ft_ka1_bs,ft_ka1_ms,color,'-','--',1e3);
%         R2_ka1_bs(nk,:) = gof_ka1_bs.rsquare;
%         R2_ka1_ms(nk,:) = gof_ka1_ms.rsquare;
%         if nk==1
%             text(x_temp(1),1.5,ft_ka1,'FontSize',font_size/1.75,'Interpreter','latex')
%         end
%     end
% 
%     y = ka1_ens_bead(:,nk)*1e3;
%     err = ka1SE_ens_bead(:,nk)*1e3;
%     MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
%         marker_edge_color)
%     
%     y = ka1_ens_meso(:,nk)*1e3;
%     err = ka1SE_ens_meso(:,nk)*1e3;
%     MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
%         marker_edge_color)
%     
%     PolishTheFigure(xlab,xlimits,xtickers,xstyle,ystyle,angle);
%     if NormalizeRates==1
%         ylim([0 1])
%         ylabel('$\bar{k}_{a,1} \tau_0$ $(10^{-3})$','FontSize',font_size,'Interpreter','latex')
%     else
%         ylim([0 400])
%         ylabel('$\bar{k}_{a,1}$ (kHz)','FontSize',font_size,'Interpreter','latex')
%     end
    
    
    % 3) Exchange rate
    figno = figno+1;
    figure(figno)
    if ~isempty(ft_kexc_bs)
        [coeff_kexc_bs(nk,:),ci_kexc_bs(nk,:),coeff_kexc_ms(nk,:),ci_kexc_ms(nk,:)] = ...
            PlotTheEmpericalFits(x_plt,ft_kexc_bs,ft_kexc_ms,color,'-','--',1e3);
        R2_kexc_bs(nk,:) = gof_kexc_bs.rsquare;
        R2_kexc_ms(nk,:) = gof_kexc_ms.rsquare;
%         if nk==1
%             text(x_temp(1),1.5,ft_kexc,'FontSize',font_size/1.75,'Interpreter','latex')
%         end
    end
    
    y = kexc_ens_bead(:,nk)*1e3;
    err = kexcSE_ens_bead(:,nk)*1e3;
%     if get_ka_from_adding==1
%         err = NaN*ones(size(y));
%     end
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
        marker_edge_color)
    
    y = kexc_ens_meso(:,nk)*1e3;
    err = kexcSE_ens_meso(:,nk)*1e3;
%     if get_ka_from_adding==1
%         err = NaN*ones(size(y));
%     end
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
        marker_edge_color)
    
    PolishTheFigure(xlab,xlimits,xtickers,xstyle,ystyle,angle);
    if NormalizeRates==1
        ylabel('$\bar{k}_{exc} \tau_0$ $(10^{-3})$','FontSize',font_size,'Interpreter','latex')
        ylim([0 15])
    else
        ylabel('$\bar{k}_{exc}$ (kHz)','FontSize',font_size,'Interpreter','latex')
        ylim([0 400])
    end
    
    
    % 4) Repeat rate
    figno = figno+1;
    figure(figno)
    if ~isempty(ft_krpt_bs)
        [coeff_krpt_bs(nk,:),ci_krpt_bs(nk,:),coeff_krpt_ms(nk,:),ci_krpt_ms(nk,:)] = ...
            PlotTheEmpericalFits(x_plt,ft_krpt_bs,ft_krpt_ms,color,'-','--',1e3);
        R2_krpt_bs(nk,:) = gof_krpt_bs.rsquare;
        R2_krpt_ms(nk,:) = gof_krpt_ms.rsquare;
%         if nk==1
%             text(x_temp(1),1.5,ft_krpt,'FontSize',font_size/1.75,'Interpreter','latex')
%         end
    end
    
    y = krpt_ens_bead(:,nk)*1e3;
    err = krptSE_ens_bead(:,nk)*1e3;
%     if get_ka_from_adding==1
%         err = NaN*ones(size(y));
%     end
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
        marker_edge_color)
    
    y = krpt_ens_meso(:,nk)*1e3;
    err = krptSE_ens_meso(:,nk)*1e3;
%     if get_ka_from_adding==1
%         err = NaN*ones(size(y));
%     end
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
        marker_edge_color)
    
    PolishTheFigure(xlab,xlimits,xtickers,xstyle,ystyle,angle);
    if NormalizeRates==1
        ylim([0 15])
        ylabel('$\bar{k}_{rpt} \tau_0$ $(10^{-3})$','FontSize',font_size,'Interpreter','latex')
    else
        ylim([0 400])
        ylabel('$\bar{k}_{rpt}$ (kHz)','FontSize',font_size,'Interpreter','latex')
    end
    
    
    % 5) Fraction of attached bonds
    figno = figno+1;
    figure(figno)

    y = fa_ens_bead(:,nk);
    err = faSE_ens_bead(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
        marker_edge_color)
    
    y = fa_ens_meso(:,nk);
    err = faSE_ens_meso(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
        marker_edge_color)
    
    PolishTheFigure(xlab,xlimits,xtickers,xstyle,ystyle,angle);
    ylabel('$\bar{f}_a$','FontSize',font_size,'Interpreter','latex')
    ylim([0 1])

    
    % 6) Fraction of detached bonds
    figno = figno+1;
    figure(figno)

    y = fd_ens_bead(:,nk);
    err = fdSE_ens_bead(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
        marker_edge_color)
    
    y = fd_ens_meso(:,nk);
    err = fdSE_ens_meso(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
        marker_edge_color)
    
    PolishTheFigure(xlab,xlimits,xtickers,xstyle,ystyle,angle);
    ylabel('$\bar{f}_d$','FontSize',font_size,'Interpreter','latex')
    ylim([0 1])


    % 7) Attached bond lifetime (tau_a)
    figno = figno+1;
    figure(figno)
%     if NormalizeRates==1
%         plot(xlimits,1./([kd_in kd_in]*tau0),'k--');
%     else
%         plot(xlimits,[kd_in kd_in]/1e3,'k--');
%     end
    
%     [y,err] = InvertTheRate(kd_ens_bead,kdSE_ens_bead,nk);
%     MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
%         marker_edge_color)
    y = tau_a_ens_bead(:,nk);
    err = tau_aSE_ens_bead(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
        marker_edge_color)
    
%     [y,err] = InvertTheRate(kd_ens_meso,kdSE_ens_meso,nk);
%     MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
%         marker_edge_color)
    y = tau_a_ens_meso(:,nk);
    err = tau_aSE_ens_meso(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
        marker_edge_color)
    
    PolishTheFigure(xlab,xlimits,xtickers,xstyle,ystyle,angle);
    if NormalizeRates==1
        ylim([0 1500])
        ylabel('$\bar \tau_{a}/\tau_0$','FontSize',font_size,'Interpreter','latex')
    else
        ylim([0 0.8])
        ylabel('$\bar \tau_{a}$ (s)','FontSize',font_size,'Interpreter','latex')
    end

    
    % 8) Detached bond lifetime (tau_d)
    figno = figno+1;
    figure(figno)
%     if ~isempty(ft_ka_bs)
%         PlotTheInverseFit(ft_ka_bs,ft_ka_ms,x_plt,color,1.5)
%     end
    
    y = tau_d_ens_bead(:,nk);
    err = tau_dSE_ens_bead(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
        marker_edge_color)
%     [y,err] = InvertTheRate(ka_ens_bead,kaSE_ens_bead,nk);
%     MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
%         marker_edge_color)
    
    y = tau_d_ens_meso(:,nk);
    err = tau_dSE_ens_meso(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
        marker_edge_color)
%     [y,err] = InvertTheRate(ka_ens_meso,kaSE_ens_meso,nk);
%     MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
%         marker_edge_color)
    
    PolishTheFigure(xlab,xlimits,xtickers,xstyle,ystyle,angle);
    if NormalizeRates==1
        ylim([0 500])
        ylabel('$\bar \tau_{d}/\tau_0$','FontSize',font_size,'Interpreter','latex')
    else
        ylim([0 0.8])
        ylabel('$\bar \tau_{d}$ (s)','FontSize',font_size,'Interpreter','latex')
    end
%     set(gca,'yscale','log')

    
    % 9) Detached bond lifetime prior to repeat (tau_rpt)
    figno = figno+1;
    figure(figno)
%     if ~isempty(ft_ka_bs)
%         PlotTheInverseFit(ft_krpt_bs,ft_krpt_ms,x_plt,color,1.5)
%     end
    
    y = tau_rpt_ens_bead(:,nk);
    err = tau_rptSE_ens_bead(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
        marker_edge_color)
%     [y,err] = InvertTheRate(krpt_ens_bead,krptSE_ens_bead,nk);
%     MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
%         marker_edge_color)
    
    y = tau_rpt_ens_meso(:,nk);
    err = tau_rptSE_ens_meso(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
        marker_edge_color)
%     [y,err] = InvertTheRate(krpt_ens_meso,krptSE_ens_meso,nk);
%     MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
%         marker_edge_color)
    
    PolishTheFigure(xlab,xlimits,xtickers,xstyle,ystyle,angle);
    if NormalizeRates==1
        ylim([0 500])
        ylabel('$\bar \tau_{rpt}/\tau_0$','FontSize',font_size,'Interpreter','latex')
    else
        ylim([0 0.8])
        ylabel('$\bar \tau_{rpt}$ (s)','FontSize',font_size,'Interpreter','latex')
    end
       
    
    % 10) Detached bond lifetime prior to exchange (tau_exc)
    figno = figno+1;
    figure(figno)
%     if ~isempty(ft_ka_bs)
%         PlotTheInverseFit(ft_kexc_bs,ft_kexc_ms,x_plt,color,1.5)
%     end
    
    y = tau_exc_ens_bead(:,nk);
    err = tau_excSE_ens_bead(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
        marker_edge_color)    
%     [y,err] = InvertTheRate(kexc_ens_bead,kexcSE_ens_bead,nk);
%     MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
%         marker_edge_color)
    
    y = tau_exc_ens_meso(:,nk);
    err = tau_excSE_ens_meso(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
        marker_edge_color)
%     [y,err] = InvertTheRate(kexc_ens_meso,kexcSE_ens_meso,nk);
%     MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
%         marker_edge_color)
    
    PolishTheFigure(xlab,xlimits,xtickers,xstyle,ystyle,angle);
    if NormalizeRates==1
        ylim([0 500])
        ylabel('$\bar \tau_{exc}/\tau_0$','FontSize',font_size,'Interpreter','latex')
    else
        ylim([0 0.8])
        ylabel('$\bar \tau_{exc}$ (s)','FontSize',font_size,'Interpreter','latex')
    end
    
    
    % 11) Adjusted bond lifetime (tau_adj)
    figno = figno+1;
    figure(figno)
    
    y = tau_adj_ens_bead(:,nk);
    err = tau_adjSE_ens_bead(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
        marker_edge_color)
    
    y = tau_adj_ens_meso(:,nk);
    err = tau_adjSE_ens_meso(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
        marker_edge_color)
    
    PolishTheFigure(xlab,xlimits,xtickers,xstyle,ystyle,angle);
    if NormalizeRates==1
        ylim([0 1500])
        ylabel('$\bar \tau_{adj}/\tau_0$','FontSize',font_size,'Interpreter','latex')
    else
        ylim([0 0.8])
        ylabel('$\bar \tau_{adj}$ (s)','FontSize',font_size,'Interpreter','latex')
    end
    
    % 12) Renormalized bond lifetime (tau_exc + tau_adj)
    figno = figno+1;
    figure(figno)
    
    y = tau_rnm_ens_bead(:,nk);
    err = tau_rnmSE_ens_bead(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
        marker_edge_color)
    
    y = tau_rnm_ens_meso(:,nk);
    err = tau_rnmSE_ens_meso(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
        marker_edge_color)
    
    PolishTheFigure(xlab,xlimits,xtickers,xstyle,ystyle,angle);
    if NormalizeRates==1
        ylim([0 1500])
        ylabel('$\bar \tau_{rnm}/\tau_0$','FontSize',font_size,'Interpreter','latex')
    else
        ylim([0 0.8])
        ylabel('$\bar \tau_{rnm}$ (s)','FontSize',font_size,'Interpreter','latex')
    end
    
    % 13) Avg attached, adjusted, and renormalized bond lifetimes 
    figno = figno+1;
    figure(figno)
    
    y = tau_a_ens_bead(:,nk);
    err = tau_aSE_ens_bead(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
        marker_edge_color)

    y = tau_a_ens_meso(:,nk);
    err = tau_aSE_ens_meso(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
        marker_edge_color)
    
    y = tau_adj_ens_bead(:,nk);
    err = tau_adjSE_ens_bead(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
        marker_edge_color)
    
    y = tau_adj_ens_meso(:,nk);
    err = tau_adjSE_ens_meso(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
        marker_edge_color)
       
    y = tau_rnm_ens_bead(:,nk);
    err = tau_rnmSE_ens_bead(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_bs,...
        marker_edge_color)
    
    y = tau_rnm_ens_meso(:,nk);
    err = tau_rnmSE_ens_meso(:,nk);
    MakeErrorPlot(figno,x_temp,y,err,color,line_style_dat,marker_ms,...
        marker_edge_color)
    
    PolishTheFigure(xlab,xlimits,xtickers,xstyle,ystyle,angle);
    ylim([0 1500])
    ylabel('$\bar \tau_{a}/\tau_0$','FontSize',font_size,'Interpreter','latex')

end

output_dir = ['Output Plots/Wrt. ',save_tag];
if ~isfolder(output_dir)
    mkdir(output_dir)
end

% Save the figures
XSize = 300;
YSize = 300;

SaveTheKineticsFigures(fig_nom+1,XSize,YSize,'kd',output_dir,save_tag);
pbaspect([2 1 1])
set(gcf,'Position',[100 400 310 310])

SaveTheKineticsFigures(fig_nom+2,XSize,YSize,'ka',output_dir,save_tag);
SaveTheKineticsFigures(fig_nom+3,XSize,YSize,'kexc',output_dir,save_tag);
SaveTheKineticsFigures(fig_nom+4,XSize,YSize,'krpt',output_dir,save_tag);

SaveTheKineticsFigures(fig_nom+5,XSize,YSize,'fa',output_dir,save_tag);
pbaspect([2 1 1])
set(gcf,'Position',[100 400 305 305])
SaveTheKineticsFigures(fig_nom+6,XSize,YSize,'fd',output_dir,save_tag);
pbaspect([2 1 1])
set(gcf,'Position',[100 400 305 305])

SaveTheKineticsFigures(fig_nom+8,XSize,YSize,'tau_d',output_dir,save_tag);
set(gca,'FontSize',20/1.6)
set(gcf,'Position',[100 400 190 190])
xticks([min(xtickers) max(xtickers)])
yticks([0 500])
SaveTheKineticsFigures(fig_nom+9,XSize,YSize,'tau_rpt',output_dir,save_tag);
set(gca,'FontSize',20/1.6)
set(gcf,'Position',[100 400 190 190])
xticks([min(xtickers) max(xtickers)])
yticks([0 500])
SaveTheKineticsFigures(fig_nom+10,XSize,YSize,'tau_exc',output_dir,save_tag);
set(gca,'FontSize',20/1.6)
set(gcf,'Position',[100 400 190 190])
xticks([min(xtickers) max(xtickers)])
yticks([0 500])

sz_2 = 322;
SaveTheKineticsFigures(fig_nom+7,XSize,YSize,'tau_a',output_dir,save_tag);
set(gcf,'Position',[100 400 sz_2 sz_2])
ylim([0 1500])
pbaspect([1.65 1 1])
SaveTheKineticsFigures(fig_nom+11,XSize,YSize,'tau_adj',output_dir,save_tag);
set(gcf,'Position',[100 400 sz_2 sz_2])
pbaspect([1.65 1 1])
ylim([0 1500])
SaveTheKineticsFigures(fig_nom+12,XSize,YSize,'tau_rnm',output_dir,save_tag);
set(gcf,'Position',[100 400 sz_2 sz_2])
pbaspect([1.65 1 1])
ylim([0 1500])

SaveTheKineticsFigures(fig_nom+13,XSize,YSize,'tau_a_comb',output_dir,save_tag);

PlotTheParameterFits(coeff_ka_bs,coeff_ka_ms,...
    coeff_kexc_bs,coeff_kexc_ms,...
    coeff_krpt_bs,coeff_krpt_ms,...
    ci_ka_bs,ci_ka_ms,...
    ci_kexc_bs,ci_kexc_ms,...
    ci_krpt_bs,ci_krpt_ms,...
    save_tag,N_Kuhns)

% Write the fitted parameters to tables
folder_temp = 'Fitted Parameter Tables';
if ~isfolder(folder_temp)
    mkdir(folder_temp)
end
if contains(save_tag,'c')
    N = N_Kuhns;
    ka_max = coeff_ka_bs(:,2);
    pm_ka = ci_ka_bs(:,2);
%     if get_ka_from_adding==0
%         c_min = coeff_ka_bs(:,1);
%         pm_c = ci_ka_bs(:,1);
%     else
        c_min = (coeff_ka_bs(:,1)).^(-3);
%         c_min = (coeff_ka_bs(:,1)).^3;
%         pm_c = 1/3*(coeff_ka_bs(:,1).^2).*ci_ka_bs(:,1)*sqrt(3); % from variance formula
        pm_c = VarianceFunction(coeff_ka_bs(:,1),ci_ka_bs(:,1));
%     end
    R2 = R2_ka_bs;
    tab = table(N,ka_max,pm_ka,c_min,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.ka_bs.txt'];
    writetable(tab,file_name);
    
    ka_max = coeff_ka_ms(:,2);
    pm_ka = ci_ka_ms(:,2);
%     if get_ka_from_adding==0
%         c_min = coeff_ka_bs(:,1);
%         pm_c = ci_ka_ms(:,1);
%     else
        c_min = (coeff_ka_bs(:,1)).^(-3);
%         c_min = (coeff_ka_bs(:,1)).^3;
%         pm_c = 1/3*(coeff_ka_ms(:,1).^2).*ci_ka_ms(:,1)*sqrt(3); % from variance formula
        pm_c = VarianceFunction(coeff_ka_ms(:,1),ci_ka_ms(:,1));
%     end
    R2 = R2_ka_ms;
    tab = table(N,ka_max,pm_ka,c_min,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.ka_ms.txt'];
    writetable(tab,file_name);
    
%     ka1_max = coeff_ka1_bs(:,2);
%     pm_ka = ci_ka1_bs(:,2);
%     c_min = (coeff_ka1_bs(:,1)).^(-3);
% %     c_min = (coeff_ka1_bs(:,1)).^3;
% %  	pm_c = 1/3*(coeff_ka1_bs(:,1).^2).*ci_ka1_bs(:,1)*sqrt(3); % from variance formula
%     pm_c = VarianceFunction(coeff_ka1_bs(:,1),ci_ka1_bs(:,1));
%     R2 = R2_ka1_bs;
%     tab = table(N,ka1_max,pm_ka,c_min,pm_c,R2);
%     file_name = [folder_temp,'/',save_tag,'.ka1_bs.txt'];
%     writetable(tab,file_name);
%     
%     ka1_max = coeff_ka1_ms(:,2);
%     pm_ka = ci_ka1_ms(:,2);
%     c_min = (coeff_ka1_ms(:,1)).^(-3);
% %     c_min = (coeff_ka1_ms(:,1)).^3;
% %  	pm_c = 1/3*(coeff_ka1_ms(:,1).^2).*ci_ka1_ms(:,1)*sqrt(3); % from variance formula
%     pm_c = VarianceFunction(coeff_ka1_ms(:,1),ci_ka1_ms(:,1));
%     R2 = R2_ka1_ms;
%     tab = table(N,ka1_max,pm_ka,c_min,pm_c,R2);
%     file_name = [folder_temp,'/',save_tag,'.ka1_ms.txt'];
%     writetable(tab,file_name);
    
    krpt_max = coeff_krpt_bs(:,2);
    pm_ka = ci_krpt_bs(:,2);
    c_min = (coeff_krpt_bs(:,1)).^(-3);
%     c_min = (coeff_krpt_bs(:,1)).^3;
%  	pm_c = 1/3*(coeff_krpt_bs(:,1).^2).*ci_krpt_bs(:,1)*sqrt(3); % from variance formula
    pm_c = VarianceFunction(coeff_krpt_bs(:,1),ci_krpt_bs(:,1));
    R2 = R2_krpt_bs;
    tab = table(N,krpt_max,pm_ka,c_min,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.krpt_bs.txt'];
    writetable(tab,file_name);
    
    krpt_max = coeff_krpt_ms(:,2);
    pm_ka = ci_krpt_ms(:,2);
    c_min = (coeff_krpt_ms(:,1)).^(-3);
%     c_min = (coeff_krpt_ms(:,1)).^3;
%  	pm_c = 1/3*(coeff_krpt_ms(:,1).^2).*ci_krpt_ms(:,1)*sqrt(3); % from variance formula
    pm_c = VarianceFunction(coeff_krpt_ms(:,1),ci_krpt_ms(:,1));
    R2 = R2_krpt_ms;
    tab = table(N,krpt_max,pm_ka,c_min,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.krpt_ms.txt'];
    writetable(tab,file_name);
    
    kexc_max = coeff_kexc_bs(:,2);
    pm_ka = ci_kexc_bs(:,2);
    c_min = (coeff_kexc_bs(:,1)).^(-3);
%     c_min = (coeff_kexc_bs(:,1)).^3;
%  	pm_c = 1/3*(coeff_kexc_bs(:,1).^2).*ci_kexc_bs(:,1)*sqrt(3); % from variance formula
    pm_c = VarianceFunction(coeff_kexc_bs(:,1),ci_kexc_bs(:,1));
    R2 = R2_kexc_bs;
    tab = table(N,kexc_max,pm_ka,c_min,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.kexc_bs.txt'];
    writetable(tab,file_name);
    
    kexc_max = coeff_kexc_ms(:,2);
    pm_ka = ci_kexc_ms(:,2);
    c_min = (coeff_kexc_ms(:,1)).^(-3);
%     c_min = (coeff_kexc_ms(:,1)).^3;
%  	pm_c = 1/3*(coeff_kexc_ms(:,1).^2).*ci_kexc_ms(:,1)*sqrt(3); % from variance formula
    pm_c = VarianceFunction(coeff_kexc_ms(:,1),ci_kexc_ms(:,1));
    R2 = R2_kexc_ms;
    tab = table(N,kexc_max,pm_ka,c_min,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.kexc_ms.txt'];
    writetable(tab,file_name);
elseif contains(save_tag,'d')
    N = N_Kuhns;
    ka_max = coeff_ka_bs(:,2);
    pm_ka = ci_ka_bs(:,2);
    d_max = coeff_ka_bs(:,1);
    pm_c = ci_ka_bs(:,1);
    R2 = R2_ka_bs;
    tab = table(N,ka_max,pm_ka,d_max,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.ka_bs.txt'];
    writetable(tab,file_name);
    
    ka_max = coeff_ka_ms(:,2);
    pm_ka = ci_ka_ms(:,2);
    d_max = coeff_ka_ms(:,1);
    pm_c = ci_ka_ms(:,1);
    R2 = R2_ka_ms;
    tab = table(N,ka_max,pm_ka,d_max,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.ka_ms.txt'];
    writetable(tab,file_name);
    
    ka1_max = coeff_ka1_bs(:,2);
    pm_ka = ci_ka1_bs(:,2);
    d_max = coeff_ka1_bs(:,1);
    pm_c = ci_ka1_bs(:,1);
    R2 = R2_ka1_bs;
    tab = table(N,ka1_max,pm_ka,d_max,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.ka1_bs.txt'];
    writetable(tab,file_name);
    
    ka1_max = coeff_ka1_ms(:,2);
    pm_ka = ci_ka1_ms(:,2);
    d_max = coeff_ka1_ms(:,1);
    pm_c = ci_ka1_ms(:,1);
    R2 = R2_ka1_ms;
    tab = table(N,ka1_max,pm_ka,d_max,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.ka1_ms.txt'];
    writetable(tab,file_name);
    
    krpt_max = coeff_krpt_bs(:,2);
    pm_ka = ci_krpt_bs(:,2);
    d_max = coeff_krpt_bs(:,1);
    pm_c = ci_krpt_bs(:,1);
    R2 = R2_krpt_bs;
    tab = table(N,krpt_max,pm_ka,d_max,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.krpt_bs.txt'];
    writetable(tab,file_name);
    
    krpt_max = coeff_krpt_ms(:,2);
    pm_ka = ci_krpt_ms(:,2);
    d_max = coeff_krpt_ms(:,1);
    pm_c = ci_krpt_ms(:,1);
    R2 = R2_krpt_ms;
    tab = table(N,krpt_max,pm_ka,d_max,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.krpt_ms.txt'];
    writetable(tab,file_name);
    
    kexc_max = coeff_kexc_bs(:,2);
    pm_ka = ci_kexc_bs(:,2);
    d_max = coeff_kexc_bs(:,1);
    pm_c = ci_kexc_bs(:,1);
    R2 = R2_kexc_bs;
    tab = table(N,kexc_max,pm_ka,d_max,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.kexc_bs.txt'];
    writetable(tab,file_name);
    
    kexc_max = coeff_kexc_ms(:,2);
    pm_ka = ci_kexc_ms(:,2);
    d_max = coeff_kexc_ms(:,1);
    pm_c = ci_kexc_ms(:,1);
    R2 = R2_kexc_ms;
    tab = table(N,kexc_max,pm_ka,d_max,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.kexc_ms.txt'];
    writetable(tab,file_name);
elseif contains(save_tag,'phi')
    N = N_Kuhns;
    ka_max = coeff_ka_bs(:,2);
    pm_ka = ci_ka_bs(:,2);
%     if get_ka_from_adding==0
%         phi_min = coeff_ka_bs(:,1);
%         pm_c = ci_ka_bs(:,1);
%     else
        phi_min = (coeff_ka_bs(:,1)).^(-3);
%         pm_c = 1/3*(coeff_ka_bs(:,1).^2).*ci_ka_bs(:,1)*sqrt(3); % from variance formula
        pm_c = VarianceFunction(coeff_ka_bs(:,1),ci_ka_bs(:,1));
%     end
    R2 = R2_ka_bs;
    tab = table(N,ka_max,pm_ka,phi_min,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.ka_bs.txt'];
    writetable(tab,file_name);
    
    ka_max = coeff_ka_ms(:,2);
    pm_ka = ci_ka_ms(:,2);
%     if get_ka_from_adding==0
%         phi_min = coeff_ka_bs(:,1);
%         pm_c = ci_ka_ms(:,1);
%     else
        phi_min = (coeff_ka_ms(:,1)).^(-3);
%         pm_c = 1/3*(coeff_ka_ms(:,1).^2).*ci_ka_ms(:,1)*sqrt(3); % from variance formula
        pm_c = VarianceFunction(coeff_ka_ms(:,1),ci_ka_ms(:,1));
%     end
    R2 = R2_ka_ms;
    tab = table(N,ka_max,pm_ka,phi_min,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.ka_ms.txt'];
    writetable(tab,file_name);
    
%     ka1_max = coeff_ka1_bs(:,2);
%     pm_ka = ci_ka1_bs(:,2);
%     phi_min = (coeff_ka1_bs(:,1)).^(-3);
% %     pm_c = 1/3*(coeff_ka1_bs(:,1).^2).*ci_ka1_bs(:,1)*sqrt(3); % from variance formula
%     pm_c = VarianceFunction(coeff_ka1_bs(:,1),ci_ka1_bs(:,1));
%     R2 = R2_ka1_bs;
%     tab = table(N,ka1_max,pm_ka,phi_min,pm_c,R2);
%     file_name = [folder_temp,'/',save_tag,'.ka1_bs.txt'];
%     writetable(tab,file_name);
%     
%     ka1_max = coeff_ka1_ms(:,2);
%     pm_ka = ci_ka1_ms(:,2);
%     phi_min = (coeff_ka1_ms(:,1)).^(-3);
% %     pm_c = 1/3*(coeff_ka1_ms(:,1).^2).*ci_ka1_ms(:,1)*sqrt(3); % from variance formula
%     pm_c = VarianceFunction(coeff_ka1_ms(:,1),ci_ka1_ms(:,1));
%     R2 = R2_ka1_ms;
%     tab = table(N,ka1_max,pm_ka,phi_min,pm_c,R2);
%     file_name = [folder_temp,'/',save_tag,'.ka1_ms.txt'];
%     writetable(tab,file_name);
    
    krpt_max = coeff_krpt_bs(:,2);
    pm_ka = ci_krpt_bs(:,2);
    phi_min = (coeff_krpt_bs(:,1)).^(-3);
%     pm_c = 1/3*(coeff_krpt_bs(:,1).^2).*ci_krpt_bs(:,1)*sqrt(3); % from variance formula
    pm_c = VarianceFunction(coeff_krpt_bs(:,1),ci_krpt_bs(:,1));
    R2 = R2_krpt_bs;
    tab = table(N,krpt_max,pm_ka,phi_min,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.krpt_bs.txt'];
    writetable(tab,file_name);
    
    krpt_max = coeff_krpt_ms(:,2);
    pm_ka = ci_krpt_ms(:,2);
    phi_min = (coeff_krpt_ms(:,1)).^(-3);
%     pm_c = 1/3*(coeff_krpt_ms(:,1).^2).*ci_krpt_ms(:,1)*sqrt(3); % from variance formula
    pm_c = VarianceFunction(coeff_krpt_ms(:,1),ci_krpt_ms(:,1));
    R2 = R2_krpt_ms;
    tab = table(N,krpt_max,pm_ka,phi_min,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.krpt_ms.txt'];
    writetable(tab,file_name);
    
    kexc_max = coeff_kexc_bs(:,2);
    pm_ka = ci_kexc_bs(:,2);
    phi_min = (coeff_kexc_bs(:,1)).^(-3);
%     pm_c = 1/3*(coeff_kexc_bs(:,1).^2).*ci_kexc_bs(:,1)*sqrt(3); % from variance formula
    pm_c = VarianceFunction(coeff_kexc_bs(:,1),ci_kexc_bs(:,1));
    R2 = R2_kexc_bs;
    tab = table(N,kexc_max,pm_ka,phi_min,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.kexc_bs.txt'];
    writetable(tab,file_name);
    
    kexc_max = coeff_kexc_ms(:,2);
    pm_ka = ci_kexc_ms(:,2);
    phi_min = (coeff_kexc_ms(:,1)).^(-3);
%     pm_c = 1/3*(coeff_kexc_ms(:,1).^2).*ci_kexc_ms(:,1)*sqrt(3); % from variance formula
    pm_c = VarianceFunction(coeff_kexc_ms(:,1),ci_kexc_ms(:,1));
    R2 = R2_kexc_ms;
    tab = table(N,kexc_max,pm_ka,phi_min,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.kexc_ms.txt'];
    writetable(tab,file_name);
else
    N = N_Kuhns;
    ka_max = coeff_ka_bs(:,2);
    pm_ka = ci_ka_bs(:,2);
    co_min = coeff_ka_bs(:,1);
    pm_c = ci_ka_bs(:,1);
    R2 = R2_ka_bs;
    tab = table(N,ka_max,pm_ka,co_min,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.ka_bs.txt'];
    writetable(tab,file_name);
    
    ka_max = coeff_ka_ms(:,2);
    pm_ka = ci_ka_ms(:,2);
    co_min = coeff_ka_ms(:,1);
    pm_c = ci_ka_ms(:,1);
    R2 = R2_ka_ms;
    tab = table(N,ka_max,pm_ka,co_min,pm_c,R2);
    file_name = [folder_temp,'/',save_tag,'.ka_ms.txt'];
    writetable(tab,file_name);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotTheInverseFit(ft_bs,ft_ms,x,color,line_width)

fit_orig = ft_bs{1};
y_orig = fit_orig(x);
y_new = 1./y_orig;
plot(x,y_new,'LineWidth',line_width,'Color',color);

fit_orig = ft_ms{1};
y_orig = fit_orig(x);
y_new = 1./y_orig;

plot(x,y_new,'LineWidth',line_width,'LineStyle','--','Color',color);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RemoveOutlier(nk_indices,ph_indices,mt)

global ka_ens_meso kd_ens_meso...
    ka1_ens_meso...
    fa_ens_meso fd_ens_meso...
    kexc_ens_meso ...
    krpt_ens_meso...
    ka_ens_bead kd_ens_bead...
    ka1_ens_bead...
    fa_ens_bead fd_ens_bead...
    kexc_ens_bead ...
    krpt_ens_bead...
    tau_a_ens_meso...
    tau_d_ens_meso...
    tau_exc_ens_meso...
    tau_rpt_ens_meso ...
    tau_adj_ens_meso...
    tau_rnm_ens_meso...
    tau_a_ens_bead...
    tau_d_ens_bead...
    tau_exc_ens_bead...
    tau_rpt_ens_bead...
    tau_adj_ens_bead...
    tau_rnm_ens_bead...
    
for mt_temp=0:1
    nk_temp = nk_indices(mt==mt_temp);
    ph_temp = ph_indices(mt==mt_temp);
    for i=1:length(nk_temp)
        nk = nk_temp(i);
        ph = ph_temp(i);
        
        if mt_temp==0
            ka_ens_bead(ph,nk) = NaN;
            ka1_ens_bead(ph,nk) = NaN;
            kd_ens_bead(ph,nk) = NaN;
            kexc_ens_bead(ph,nk) = NaN;
            krpt_ens_bead(ph,nk) = NaN;
            fa_ens_bead(ph,nk) = NaN;
            fd_ens_bead(ph,nk) = NaN;
            tau_a_ens_bead(ph,nk) = NaN;
            tau_d_ens_bead(ph,nk) = NaN;
            tau_exc_ens_bead(ph,nk) = NaN;
            tau_rpt_ens_bead(ph,nk) = NaN;
            tau_adj_ens_bead(ph,nk) = NaN;
            tau_rnm_ens_bead(ph,nk) = NaN;
        elseif mt_temp==1
            ka_ens_meso(ph,nk) = NaN;
            ka1_ens_meso(ph,nk) = NaN;
            kd_ens_meso(ph,nk) = NaN;
            kexc_ens_meso(ph,nk) = NaN;
            krpt_ens_meso(ph,nk) = NaN;
            fa_ens_meso(ph,nk) = NaN;
            fd_ens_meso(ph,nk) = NaN;
            tau_a_ens_meso(ph,nk) = NaN;
            tau_d_ens_meso(ph,nk) = NaN;
            tau_exc_ens_meso(ph,nk) = NaN;
            tau_rpt_ens_meso(ph,nk) = NaN;
            tau_adj_ens_meso(ph,nk) = NaN;
            tau_rnm_ens_meso(ph,nk) = NaN;            
        end
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tau,err] = InvertTheRate(k,k_se,i)

tau = 1./k(:,i);
err = (k_se(:,i).*(k(:,i)).^(-2));
err(err>0.5*tau)= NaN;    % This is just so that the plot is not 
                        % obscured by the error bars. At very low chain
                        % concentrations or high chain separations, the
                        % error will diverge along with the values of tau

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var = VarianceFunction(x,ci)

% var = sqrt(3)*3*(x.^2).*ci(:,1); % from variance formula for x^3
var = 3*(x.^(-4)).*ci; % from variance formula for x^3
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotTheParameterFits(coeff_ka_bs,coeff_ka_ms,...
    coeff_kexc_bs,coeff_kexc_ms,...
    coeff_krpt_bs,coeff_krpt_ms,...
    ci_ka_bs,ci_ka_ms,...
    ci_kexc_bs,ci_kexc_ms,...
    ci_krpt_bs,ci_krpt_ms,...
    save_tag,N_Kuhns)

global get_ka_from_adding

folder = 'Output Plots/Fitted Param';
if ~isfolder(folder)
    mkdir(folder)
end

if contains(save_tag,'lambda')
    y_lab = '$\lambda_{max}$';
elseif contains(save_tag,'d')
    y_lab = '$d_{max}/(Nb)$';
elseif contains(save_tag,'phi')
    y_lab = '$\phi_{min}$';
elseif contains(save_tag,'c_open')
    y_lab = '$c_{open,min}b^3$ $(10^{-3})$';
else
    y_lab = '$c_{min}b^3$ $(10^{-3})$';
end

tit = 'Fitted to $k_a$';

% min or max int from fitting ka
fig_no = 201;

y_bs = coeff_ka_bs(:,1);
err_bs = ci_ka_bs(:,1);
y_ms = coeff_ka_ms(:,1);
err_ms = ci_ka_ms(:,1);
ylims = [0 Inf];
if get_ka_from_adding && contains(save_tag,'c')
%     err_bs = 1/3*(y_bs.^2).*err_bs*sqrt(3)*1e3; % from variance formula
%     err_ms = 1/3*(y_ms.^2).*err_ms*sqrt(3)*1e3; % from variance formula
    err_bs = VarianceFunction(y_bs,err_bs);
    err_ms = VarianceFunction(y_ms,err_ms);
    err_bs = err_bs*1e3;
    err_ms = err_ms*1e3;
    y_bs = (y_bs.^(-3))*1e3;
    y_ms = (y_ms.^(-3))*1e3;
    ylims = [0 Inf];
    add_leg = 0;
elseif get_ka_from_adding && contains(save_tag,'d')
%     y_bs = -coeff_ka_bs(:,2)./coeff_ka_bs(:,1);
%     y_ms = -coeff_ka_ms(:,2)./coeff_ka_ms(:,1);
%     dfdd = coeff_ka_bs(:,2).*coeff_ka_bs(:,1).^(-2);
%     dfdk = 1./coeff_ka_bs(:,1);
%     err_bs = ((dfdd.^2).*(ci_ka_bs(:,2)).^2 + ...
%         (dfdk.^2).*(ci_ka_bs(:,1).^2)).^0.5; % from variance formula
%     err_ms = ((dfdd.^2).*(ci_ka_bs(:,2)).^2 + ...
%         (dfdk.^2).*(ci_ka_bs(:,1).^2)).^0.5; % from variance formula
    y_bs = coeff_ka_bs(:,1);
    y_ms = coeff_ka_ms(:,1);
    err_bs = ci_ka_bs(:,1); % from variance formula
    err_ms = ci_ka_bs(:,1); % from variance formula
    add_leg = 0;
elseif get_ka_from_adding && contains(save_tag,'phi')
%     err_bs = 1/3*(y_bs.^2).*err_bs*sqrt(3); % from variance formula
%     err_ms = 1/3*(y_ms.^2).*err_ms*sqrt(3); % from variance formula
    err_bs = VarianceFunction(y_bs,err_bs);
    err_ms = VarianceFunction(y_ms,err_ms);
    y_bs = (y_bs.^(-3));
    y_ms = (y_ms.^(-3));
    add_leg = 1;
end
PlotThePolishedFitValues(fig_no,y_bs,y_ms,err_bs,err_ms,y_lab,tit,N_Kuhns,...
    add_leg,ylims)

saveas(gcf,[folder,'/max_or_min.',save_tag,'.png'])
saveas(gcf,[folder,'/max_or_min.',save_tag,'.fig'])

% plot the maximum rates predicted by each fit
ylims = [0 10];

% min or max k from fitting ka
fig_no = 202;
tit = 'Fitted to $k_a$';

y_bs = coeff_ka_bs(:,2)*1000;
err_bs = ci_ka_bs(:,2)*1000;
y_ms = coeff_ka_ms(:,2)*1000;
err_ms = ci_ka_ms(:,2)*1000;
add_leg = 0;

PlotThePolishedFitValues(fig_no,y_bs,y_ms,err_bs,err_ms,...
    '$k_a^{semi}\tau_0$ ($10^{-3}$)',tit,N_Kuhns,add_leg,ylims)

saveas(gcf,[folder,'/ka_max.',save_tag,'.png'])
saveas(gcf,[folder,'/ka_max.',save_tag,'.fig'])

% min or max krpt from fitting krpt
fig_no = 203;
tit = 'Fitted to $k_{rpt}$';

y_bs = coeff_krpt_bs(:,2)*1000;
err_bs = ci_krpt_bs(:,2)*1000;
y_ms = coeff_krpt_ms(:,2)*1000;
err_ms = ci_krpt_ms(:,2)*1000;
add_leg = 0;

PlotThePolishedFitValues(fig_no,y_bs,y_ms,err_bs,err_ms,...
    '$k_{rpt}^{semi}\tau_0$ ($10^{-3}$)',tit,N_Kuhns,add_leg,ylims)

saveas(gcf,[folder,'/krpt_max.',save_tag,'.png'])
saveas(gcf,[folder,'/krpt_max.',save_tag,'.fig'])

% min or max kexc from fitting krpt
fig_no = 204;
tit = 'Fitted to $k_{exc}$';

y_bs = coeff_kexc_bs(:,2)*1000;
err_bs = ci_kexc_bs(:,2)*1000;
y_ms = coeff_kexc_ms(:,2)*1000;
err_ms = ci_kexc_ms(:,2)*1000;
add_leg = 0;

PlotThePolishedFitValues(fig_no,y_bs,y_ms,err_bs,err_ms,...
    '$k_{exc}^{semi}\tau_0$ ($10^{-3}$)',tit,N_Kuhns,add_leg,ylims)

saveas(gcf,[folder,'/kexc_max.',save_tag,'.png'])
saveas(gcf,[folder,'/kexc_max.',save_tag,'.fig'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotThePolishedFitValues(fig_no,y_bs,y_ms,err_bs,err_ms,y_lab,tit,...
    N_Kuhns,add_leg,ylims)

global font_size

marg = 0.145;
figure(fig_no); clf; hold on
marker_type = 'none';
marker_edge_color = 'k';
marker_size = 1;

b = bar([y_bs,y_ms]);
b(1).FaceColor = [0.35 0.35 0.35];
b(2).FaceColor = [0.85 0.85 0.85];
[~] = PlotErrorBar((1:length(N_Kuhns))'-marg,y_bs,err_bs,'k',marker_type,...
    marker_edge_color,marker_size);
[~] = PlotErrorBar((1:length(N_Kuhns))'+marg,y_ms,err_ms,'k',marker_type,...
    marker_edge_color,marker_size);

set(gca,'FontSize',font_size/1.5);

if add_leg==1
    l = legend(b,{'Bead-spring','Mesoscale'});
    l.FontSize = font_size/2;
    l.Interpreter = 'latex';
    l.Location = 'northeast';
end

ylabel(y_lab,'FontSize',font_size,'Interpreter','latex')
xlabel('$N$','FontSize',font_size,'Interpreter','latex')
xticks([1 2 3])
xticklabels({'12','18','36'})
xlim([0.5 3.5])
set(gcf,'color','w')
pbaspect([1 1 1])
set(gcf,'Position',[100 100 275 275])
title(tit,'FontSize',font_size/2,'Interpreter','latex')

ylim(ylims)

% if contains(y_lab,'$d')
%     marg = 0.145;
%     figure(fig_no+1); clf; hold on
%     marker_type = 'none';
%     marker_edge_color = 'k';
%     marker_size = 1;
%     
%     y_bs_tmp = y_bs.*(N_Kuhns/2).^0.5;
%     err_bs_tmp = err_bs.*(N_Kuhns/2).^0.5;
%     y_ms_tmp = y_ms.*(N_Kuhns/2).^0.5;
%     err_ms_tmp = err_ms.*(N_Kuhns/2).^0.5;
% 
%     br = bar([y_bs_tmp,y_ms_tmp]);
%     br(1).FaceColor = [0.35 0.35 0.35];
%     br(2).FaceColor = [0.85 0.85 0.85];
%     [~] = PlotErrorBar((1:length(N_Kuhns))'-marg,y_bs_tmp,...
%         err_bs_tmp,'k',marker_type,...
%         marker_edge_color,marker_size);
%     [~] = PlotErrorBar((1:length(N_Kuhns))'+marg,y_ms_tmp,...
%         err_ms_tmp,'k',marker_type,...
%         marker_edge_color,marker_size);
%     
%     set(gca,'FontSize',font_size/1.5);
%     
%     if add_leg==1
%         l = legend(br,{'Bead-spring','Mesoscale'});
%         l.FontSize = font_size/2;
%         l.Interpreter = 'latex';
%         l.Location = 'northeast';
%     end
%     
%     ylabel('$\lambda_{max}$','FontSize',font_size,'Interpreter','latex')
%     xlabel('$N$','FontSize',font_size,'Interpreter','latex')
%     xticks([1 2 3])
%     xticklabels({'12','18','36'})
%     xlim([0.5 3.5])
%     set(gcf,'color','w')
%     pbaspect([1 1 1])
%     set(gcf,'Position',[100 100 275 275])
%     title(tit,'FontSize',font_size/2,'Interpreter','latex')
%     
%     ylim(ylims)
%     saveas(gcf,['Output Plots/Fitted Param','/krpt_max.lambda.png'])
%     saveas(gcf,['Output Plots/Fitted Param','/krpt_max.lambda.fig'])
%     close(fig_no+1)
% end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotLifetimeHistograms(n_mod,N_Kuhns,x_lab,y_bs,y_ms,phis,nk,...
    colors,bin_width,xlims,ylims)

N_Kuhn = N_Kuhns(nk);

if nk==length(N_Kuhns)
    x_lab_temp = x_lab;
else
    x_lab_temp = '';
end

ylab = ['$n, N=$',num2str(N_Kuhn)];
if nk==1
    tit = 'Bead-spring';
else
    tit = '';
end
subplot(length(N_Kuhns),n_mod,nk*n_mod-1); hold on
MakeSubplotHisto(phis,y_bs,nk,bin_width,...
    x_lab_temp,ylab,tit,colors,0,xlims,ylims)

ylab = '';
if nk==1
    tit = 'Mesoscale';
else
    tit = '';
end
subplot(length(N_Kuhns),n_mod,nk*n_mod); hold on
MakeSubplotHisto(phis,y_ms,nk,bin_width,...
    x_lab_temp,ylab,tit,colors,1,xlims,ylims)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeSubplotHisto(phis,y,nk,bin_width,...
        xlab,ylab,tit,colors,model_no,xlims,ylims)

global font_size LineWidth tau0

sweep_rng = [12 8 4 1];
legend_entries = cell(length(sweep_rng),1);
% p_leg = cell(length(sweep_rng),1);
ct = 0;
for sp=sweep_rng
    ct = ct+1;
    phi = phis(sp);

    y_temp = y(:,sp,nk);
    y_temp(y_temp==0) = [];
    h = histogram(y_temp);
    h.FaceColor = colors(sp,:);
    max_x = max(h.BinEdges)+bin_width;
    h.BinEdges = (0:bin_width:max_x);
%     h.Normalization = 'count';
%     h.FaceAlpha = 0.8;
%     
%     x_ft = h.BinEdges;
%     x_ft = x_ft(1:end-1)+bin_width/2;
%     y_ft = h.Values;
%     if length(y_ft)>5
%         [fitresult, gof] = CreatePowerLawFit(x_ft, y_ft);
%         x_plot = logspace(-5,4,1e3);
%         plot(x_plot,fitresult(x_plot),...
%             'Color','k','LineWidth',LineWidth*1.5);
%         p_leg(ct) = plot(x_plot,fitresult(x_plot),...
%             'Color',colors(sp,:),'LineWidth',LineWidth);
%         
%         coeffs = coeffvalues(fitresult);
%     else
%         coeffs = [NaN NaN];
%         gof.rsquare = NaN;
%         p_leg(ct) = plot(NaN,NaN,...
%             'Color',colors(sp,:),'LineWidth',LineWidth);
%     end
    h.Normalization = 'probability';
    h.FaceAlpha = 0.8;
    
    x_ft = h.BinEdges;
    x_ft = x_ft(1:end-1)+bin_width/2;
    y_ft = h.Values;
    if length(y_ft)>5
        [fitresult, gof] = CreateExponentialDecayFits(x_ft, y_ft);
        x_plot = logspace(-5,4,1e3);
        plot(x_plot,fitresult(x_plot),...
            'Color','k','LineWidth',LineWidth*1.5);
        p_leg(ct) = plot(x_plot,fitresult(x_plot),...
            'Color',colors(sp,:),'LineWidth',LineWidth);
        
        coeffs = coeffvalues(fitresult);
    else
        coeffs = [NaN NaN];
        gof.rsquare = NaN;
        p_leg(ct) = plot(NaN,NaN,...
            'Color',colors(sp,:),'LineWidth',LineWidth);
    end
    
    if model_no==0 && nk==1
        legend_entries{ct} = ['$\phi = $ ',num2str(phi,'%.2f'),...
            ', $k^*=$ ',num2str(coeffs(2),'%.2e'),...
            ', $R^2=$ ',num2str(gof.rsquare,'%.2f')];
    else
        legend_entries{ct} = ['$k^*=$ ',num2str(coeffs(2),'%.2e'),...
            ', $R^2=$ ',num2str(gof.rsquare,'%.2f')];        
    end
end
xlim(xlims)
ylim(ylims)
set(gca,'FontSize',font_size/1.5)
xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
ylabel(ylab,'FontSize',font_size,'Interpreter','latex')
title(tit,'FontSize',font_size,'Interpreter','latex')

l = legend(p_leg,legend_entries);
l.FontSize = font_size/1.8;
l.Interpreter = 'latex';
l.Location = 'northeast';

set(gcf,'Color','w')
pbaspect([2.25 1 1])
set(gcf,'Position',[100 100,1080,720])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fitresult, gof] = CreateExponentialDecayFits(x_temp, y_temp)
%CREATEFIT(X_TEMP,Y_TEMP)
%  Create a fit.
%
%  Data for 'untitled fit 5' fit:
%      X Input : x_temp
%      Y Output: y_temp
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 16-Mar-2024 21:32:37


%% Fit: 'untitled fit 5'.
[xData, yData] = prepareCurveData( x_temp, y_temp );

% Set up fittype and options.
ft = fittype( 'a*exp(-b*x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.4 0.01];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 5' );
% h = plot( fitresult, xData, yData );
% legend( h, 'y_temp vs. x_temp', 'untitled fit 5', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'x_temp', 'Interpreter', 'none' );
% ylabel( 'y_temp', 'Interpreter', 'none' );
% grid on

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fitresult, gof] = CreatePowerLawFit(x_temp, y_temp)
%CREATEFIT(X_TEMP,Y_TEMP)
%  Create a fit.
%
%  Data for 'untitled fit 5' fit:
%      X Input : x_temp
%      Y Output: y_temp
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 16-Mar-2024 21:32:37


%% Fit: 'untitled fit 5'.
[xData, yData] = prepareCurveData( x_temp, y_temp );

% Set up fittype and options.
ft = fittype( 'a*x^b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.688170666082093 0.725985709942296];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 5' );
% h = plot( fitresult, xData, yData );
% legend( h, 'y_temp vs. x_temp', 'untitled fit 5', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'x_temp', 'Interpreter', 'none' );
% ylabel( 'y_temp', 'Interpreter', 'none' );
% grid on

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coeff_bs,ci_bs,coeff_ms,ci_ms] = ...
    PlotTheEmpericalFits(x_plt,ft_bs,ft_ms,color,style_bs,style_ms,y_mult)

global LineWidth

func = ft_bs{1};
y_plt = func(x_plt)*y_mult;
p = plot(x_plt,y_plt);
p.Color = color;
p.LineStyle = style_bs;
p.LineWidth = LineWidth;

coeff_bs = coeffvalues(func);
ci_bs = confint(func);
ci_bs = diff(ci_bs,1);

func = ft_ms{1};
y_plt = func(x_plt)*y_mult;
p = plot(x_plt,y_plt);
p.Color = color;
p.LineStyle = style_ms;
p.LineWidth = LineWidth;

coeff_ms = coeffvalues(func);
ci_ms = confint(func);
ci_ms = diff(ci_ms,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PolishTheFigure(xlab,xlimits,xtickers,xstyle,ystyle,angle)

global font_size

set(gca,'FontSize',font_size/1.75)
xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
set(gcf,'Color','w')
pbaspect([1 1 1])
xlim(xlimits)
xticks(xtickers)
xtickangle(angle) 
set(gca,'XScale',xstyle)
set(gca,'YScale',ystyle)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ft_ka_bs,gof_ka_bs,...
    ft_ka_ms,gof_ka_ms,...
    ft_ka1_bs,gof_ka1_bs,...
    ft_ka1_ms,gof_ka1_ms,...
    ft_kexc_bs,gof_kexc_bs,...
    ft_kexc_ms,gof_kexc_ms,...
    ft_krpt_bs,gof_krpt_bs,...
    ft_krpt_ms,gof_krpt_ms] = IntakeFittedParameters(x,save_tag,nk,N_Kuhn)

global ka_ens_bead ka_ens_meso...
    ka1_ens_bead ka1_ens_meso...
    kexc_ens_bead kexc_ens_meso...
    krpt_ens_bead krpt_ens_meso...
    get_ka_from_adding
    

ft_ka_bs = [];
gof_ka_bs = [];
ft_ka_ms = [];
gof_ka_ms = [];
ft_ka1_bs = [];
gof_ka1_bs = [];
ft_ka1_ms = [];
gof_ka1_ms = [];
ft_kexc_bs = [];
gof_kexc_bs = [];
ft_kexc_ms = [];
gof_kexc_ms = [];
ft_krpt_bs = [];
gof_krpt_bs = [];
ft_krpt_ms = [];
gof_krpt_ms = [];


ka_ens_bead_ft = ka_ens_bead;
ka_ens_meso_ft = ka_ens_meso;
kexc_ens_bead_ft = kexc_ens_bead;
kexc_ens_meso_ft = kexc_ens_meso;
krpt_ens_bead_ft = krpt_ens_bead;
krpt_ens_meso_ft = krpt_ens_meso;

k_thresh = 0.0003;

ka_ens_bead_ft(ka_ens_bead_ft<k_thresh) = NaN;
ka_ens_meso_ft(ka_ens_meso_ft<k_thresh) = NaN;
kexc_ens_bead_ft(kexc_ens_bead_ft<k_thresh) = NaN;
kexc_ens_meso_ft(kexc_ens_meso_ft<k_thresh) = NaN;
krpt_ens_bead_ft(krpt_ens_bead_ft<k_thresh) = NaN;
krpt_ens_meso_ft(krpt_ens_meso_ft<k_thresh) = NaN;

switch save_tag
    case 'c'
        [ft_ka_bs,gof_ka_bs] = ...
            createFits_ka_wrtc(x,ka_ens_bead_ft(:,nk),N_Kuhn,0);
        
        [ft_ka_ms,gof_ka_ms] = ...
            createFits_ka_wrtc(x,ka_ens_meso_ft(:,nk),N_Kuhn,1);
%         
%         [ft_ka1_bs,gof_ka1_bs] = ...
%             createFits_ka1_wrtc(x,ka1_ens_bead(:,nk),N_Kuhn,0);
%         
%         [ft_ka1_ms,gof_ka1_ms] = ...
%             createFits_ka1_wrtc(x,ka1_ens_meso(:,nk),N_Kuhn,1);

        [ft_kexc_bs,gof_kexc_bs] = ...
            createFits_kexc_wrtc(x,kexc_ens_bead_ft(:,nk),N_Kuhn,0);

        [ft_kexc_ms,gof_kexc_ms] = ...
            createFits_kexc_wrtc(x,kexc_ens_meso_ft(:,nk),N_Kuhn,1);

        [ft_krpt_bs,gof_krpt_bs] = ...
            createFits_krpt_wrtc(x,krpt_ens_bead_ft(:,nk),N_Kuhn,0);

        [ft_krpt_ms,gof_krpt_ms] = ...
            createFits_krpt_wrtc(x,krpt_ens_meso_ft(:,nk),N_Kuhn,1);

    case {'d','lambda'}
        [ft_ka_bs,gof_ka_bs] = ...
            createFits_ka_wrtd(x,ka_ens_bead_ft(:,nk),N_Kuhn,0);
        
        [ft_ka_ms,gof_ka_ms] = ...
            createFits_ka_wrtd(x,ka_ens_meso_ft(:,nk),N_Kuhn,1);
%         
%         [ft_ka1_bs,gof_ka1_bs] = ...
%             createFits_ka1_wrtd(x,ka1_ens_bead(:,nk),N_Kuhn,0);
% 
%         [ft_ka1_ms,gof_ka1_ms] = ...
%             createFits_ka1_wrtd(x,ka1_ens_meso(:,nk),N_Kuhn,1);

        [ft_kexc_bs,gof_kexc_bs] = ...
            createFits_kexc_wrtd(x,kexc_ens_bead_ft(:,nk),N_Kuhn,0);

        [ft_kexc_ms,gof_kexc_ms] = ...
            createFits_kexc_wrtd(x,kexc_ens_meso_ft(:,nk),N_Kuhn,1);

        [ft_krpt_bs,gof_krpt_bs] = ...
            createFits_krpt_wrtd(x,krpt_ens_bead_ft(:,nk),N_Kuhn,0);

        [ft_krpt_ms,gof_krpt_ms] = ...
            createFits_krpt_wrtd(x,krpt_ens_meso_ft(:,nk),N_Kuhn,1);

    case 'phi'
        [ft_ka_bs,gof_ka_bs] = ...
            createFits_ka_wrtphi(x,ka_ens_bead_ft(:,nk),N_Kuhn,0);
        
        [ft_ka_ms,gof_ka_ms] = ...
            createFits_ka_wrtphi(x,ka_ens_meso_ft(:,nk),N_Kuhn,1);
%         
%         [ft_ka1_bs,gof_ka1_bs] = ...
%             createFits_ka1_wrtphi(x,ka1_ens_bead(:,nk),N_Kuhn,0);
%         
%         [ft_ka1_ms,gof_ka1_ms] = ...
%             createFits_ka1_wrtphi(x,ka1_ens_meso(:,nk),N_Kuhn,1);
        
        [ft_kexc_bs,gof_kexc_bs] = ...
            createFits_kexc_wrtphi(x,kexc_ens_bead_ft(:,nk),N_Kuhn,0);

        [ft_kexc_ms,gof_kexc_ms] = ...
            createFits_kexc_wrtphi(x,kexc_ens_meso_ft(:,nk),N_Kuhn,1);

        [ft_krpt_bs,gof_krpt_bs] = ...
            createFits_krpt_wrtphi(x,krpt_ens_bead_ft(:,nk),N_Kuhn,0);

        [ft_krpt_ms,gof_krpt_ms] = ...
            createFits_krpt_wrtphi(x,krpt_ens_meso_ft(:,nk),N_Kuhn,1);

    case 'c_open'
        [ft_ka_bs,gof_ka_bs] = ...
            createFits_ka_wrtcopen(x,ka_ens_bead_ft(:,nk),N_Kuhn,0);
        
        [ft_ka_ms,gof_ka_ms] = ...
            createFits_ka_wrtcopen(x,ka_ens_meso_ft(:,nk),N_Kuhn,1);
        
        [ft_kexc_bs,gof_kexc_bs] = ...
            createFits_kexc_wrtcopen(x,kexc_ens_bead_ft(:,nk),N_Kuhn,0);
        
        [ft_kexc_ms,gof_kexc_ms] = ...
            createFits_kexc_wrtcopen(x,kexc_ens_meso_ft(:,nk),N_Kuhn,1);
        
        [ft_krpt_bs,gof_krpt_bs] = ...
            createFits_krpt_wrtcopen(x,krpt_ens_bead_ft(:,nk),N_Kuhn,0);

        [ft_krpt_ms,gof_krpt_ms] = ...
            createFits_krpt_wrtcopen(x,krpt_ens_meso_ft(:,nk),N_Kuhn,1);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveTheBondLifetimeFigures(figno,XSize,YSize,label,output_dir)

figure(figno)
set(gcf,'Position',[100 100 XSize YSize])
file_name = [output_dir,'/',label];
saveas(gcf,[file_name,'.png'])
saveas(gcf,[file_name,'.fig'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveTheKineticsFigures(figno,XSize,YSize,label,output_dir,save_tag)

global get_ka_from_adding

figure(figno)
set(gcf,'Position',[100 100 XSize YSize])
file_name = [output_dir,'/',label,' vs. ',save_tag];
if get_ka_from_adding==1
    file_name = [file_name,'.add'];
end
saveas(gcf,[file_name,'.png'])
saveas(gcf,[file_name,'.fig'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitializeFigures(fig_nom,no_figs)

for i = 1:no_figs
    fig_no = fig_nom+i;
    figure(fig_no); clf; hold on
end

% for i = 1:no_figs+20
%     fig_no = fig_nom+i;
%     figure(fig_no); clf; hold on
% end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = MakeErrorPlot(figno,x,y,err,color,style,marker_style,marker_edge_color)

figure(figno)
e(1) = errorbar(x,y,err);
e(1).Color = 'k';
e(1).LineStyle = style;
e(1).Marker = marker_style;
e(1).MarkerFaceColor = color;
e(1).MarkerEdgeColor = marker_edge_color;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DefineTracedBondNames(N_Kuhn,phi)

global OutputDir EnsembleDynamicsFileName AllDynamicsFileName...
    FittedDynamicsFileName R2DynamicsFileName

EnsembleDynamicsFileName = [OutputDir,'EnsembleDynamics.ka-',...
    num2str(ka_in,'%.2e'),'.m'];
AllDynamicsFileName = [OutputDir,'AllDynamicswrtTime.ka-',...
    num2str(ka_in,'%.2e'),'.m'];
R2DynamicsFileName = [OutputDir,'R2Dynamics.ka-',...
    num2str(ka_in,'%.2e'),'.m'];
FittedDynamicsFileName = [OutputDir,'FittedDynamics.ka-',...
    num2str(ka_in,'%.2e'),'.m'];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mols_edge,mols_center] = ...
    FindMoleculesNearBorders(x,y,z,types,mols,Lx,Ly,Lz,N_Kuhn,b,tether_type)

check_figs = 0;
cutoff = N_Kuhn*b/5;

x = x{1}; y = y{1}; z = z{1};
types = types{1}; mols = mols{1};

x(types~=tether_type) = [];
y(types~=tether_type) = [];
z(types~=tether_type) = [];
mols(types~=tether_type) = [];

mols_center = mols;

if check_figs==1
    figure(1); clf; hold on
    scatter3(x,y,z,'k','filled')
    xlim([-Lx/2 Lx/2])
    ylim([-Lx/2 Lx/2])
    zlim([-Lx/2 Lx/2])
    daspect([1 1 1])
    view(45,30)
end

% Remove those too close to x-bounds
mols_center((abs(abs(x)-Lx/2))<=cutoff) = [];
y((abs(abs(x)-Lx/2))<=cutoff) = [];
z((abs(abs(x)-Lx/2))<=cutoff) = [];
x((abs(abs(x)-Lx/2))<=cutoff) = [];

if check_figs==1
    figure(1); clf; hold on
    scatter3(x,y,z,'r','filled')
    xlim([-Lx/2 Lx/2])
    ylim([-Lx/2 Lx/2])
    zlim([-Lx/2 Lx/2])
    daspect([1 1 1])
    view(45,30)
end
    

% Remove those too close to y-bounds
mols_center((abs(abs(y)-Ly/2))<=cutoff) = [];
x((abs(abs(y)-Ly/2))<=cutoff) = [];
z((abs(abs(y)-Ly/2))<=cutoff) = [];
y((abs(abs(y)-Ly/2))<=cutoff) = [];

if check_figs==1
    figure(1); clf; hold on
    scatter3(x,y,z,'g','filled')
    xlim([-Lx/2 Lx/2])
    ylim([-Lx/2 Lx/2])
    zlim([-Lx/2 Lx/2])
    daspect([1 1 1])
    view(45,30)
end

% Remove those too close to z-bounds
mols_center((abs(abs(z)-Lz/2))<=cutoff) = [];
x((abs(abs(z)-Lz/2))<=cutoff) = [];
y((abs(abs(z)-Lz/2))<=cutoff) = [];
z((abs(abs(z)-Lz/2))<=cutoff) = [];

if check_figs==1
    figure(1); clf; hold on
    scatter3(x,y,z,'g','filled')
    xlim([-Lx/2 Lx/2])
    ylim([-Lx/2 Lx/2])
    zlim([-Lx/2 Lx/2])
    daspect([1 1 1])
    view(45,30)
end

mols_edge = mols;
mols_edge(ismember(mols_edge,mols_center)) = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mols_ngb = FindPotentialBondNeighbors(x,y,z,types,mols,...
    N_Kuhn,b,tether_type,mol_of_interest)

cutoff = N_Kuhn*b;
x = x{1}; y = y{1}; z = z{1};
types = types{1}; mols = mols{1};

% Examine only the tethers of each molecule
x(types~=tether_type) = [];
y(types~=tether_type) = [];
z(types~=tether_type) = [];
mols(types~=tether_type) = [];

% Find position of mol_of_interest
x_int = x(mols==mol_of_interest);
y_int = y(mols==mol_of_interest);
z_int = z(mols==mol_of_interest);

% Find positions of all other molecules
x_ngb = x(mols~=mol_of_interest);
y_ngb = y(mols~=mol_of_interest);
z_ngb = z(mols~=mol_of_interest);
mols_ngb = mols(mols~=mol_of_interest);

% Find sepatration distances between each
dist = ((x_ngb-x_int).^2 + (y_ngb-y_int).^2 + (z_ngb-z_int).^2).^0.5;

% Remove molecules for which distance is greater than cutoff
mols_ngb(dist>=cutoff) = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_st,y_st,z_st,id_st] = FindStickerID(x_int,y_int,z_int,ids_int,...
    types_int,sticker_type)

x_st = x_int(types_int==sticker_type);
y_st = y_int(types_int==sticker_type);
z_st = z_int(types_int==sticker_type);
id_st = ids_int(types_int==sticker_type);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_tmp,y_tmp,z_tmp,x_int,y_int,z_int,...
    x_th,y_th,z_th,id_th] = CenterAboutTether(x_int,y_int,z_int,ids_int,...
    types_int,tether_type,x_tmp,y_tmp,z_tmp)

x_th = x_int(types_int==tether_type);
y_th = y_int(types_int==tether_type);
z_th = z_int(types_int==tether_type);
id_th = ids_int(types_int==tether_type);
x_tmp = x_tmp-x_th;
y_tmp = y_tmp-y_th;
z_tmp = z_tmp-z_th;
x_int = x_int-x_th;
y_int = y_int-y_th;
z_int = z_int-z_th;
x_th = 0;
y_th = 0;
z_th = 0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TraceBondThroughTime(Override,Package,N_Kuhns,Np,Ns,kbT,phis,f0,...
    b,N,eaStar,edStar,D,dt)

global LengthConversion DamperConversion RawDataFileName
 
Folder = 'Traced Bonds';
if ~isfolder(Folder)
    mkdir(Folder)
end

line_width = 2;
marker_size = 6;

for BeadSpringOrMeso=0:1
    % Loop over all sims
    
    if BeadSpringOrMeso==0
        tether_type = 1;
        intermed_type = 2;
        sticker_type = 3;
    else
        tether_type = 1;
        sticker_type = 2;
    end

    
    % Samples = unique(Package(:,1));
    Nps = unique(Package(:,2));        %Number of molecules
    Nss = unique(Package(:,3));        %Number of stickeres per tether site
    kas = unique(Package(:,4));         %Temperature
    kds = unique(Package(:,5));        %Contour length of chains
    f0s = unique(Package(:,6));        %Activation energy of association
    N_Kuhns = unique(Package(:,7));        %Activation energy of dissociation
    bs = unique(Package(:,8));
    Separations = unique(Package(:,9));
    phis = unique(Package(:,10));
    kbTs = unique(Package(:,11));
    
    for j=1:length(Nps)
     for k=1:length(Nss)
      for l=1:length(kas)
       for m=1:length(kds)
        for n=1:length(f0s)
         for o=1:length(N_Kuhns)
          for p=1:length(bs)
           for q=1:length(Separations)
            for r=1:length(phis)
             for s=1:length(kbTs)
                Np = Nps(j);
                Ns = Nss(k);
                N = Np*Ns;
                ka = kas(l);
                kd = kds(m);
                f0 = f0s(n);
                N_Kuhn = N_Kuhns(o);
                b = bs(p);
                Separation = Separations(q);
                phi = phis(r);
                kbT = kbTs(s);
                
                [~,D,tau0,dtFact] = DefineTimeScale(b,LengthConversion,DamperConversion,...
                    BeadSpringOrMeso);
                
                dt = tau0/dtFact;
                
                eaStar = -log(ka*(b*LengthConversion)^2/D);
                edStar = -log(kd*(b*LengthConversion)^2/D);
                
                PackageTemp = Package(ismember(Package(:,2),Np),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,3),Ns),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,4),ka),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,5),kd),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,6),f0),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,7),N_Kuhn),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,8),b),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,9),Separation),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,10),phi),:);
                PackageTemp = PackageTemp(ismember(PackageTemp(:,11),kbT),:);
                
                
                %% Compile the Raw Data
                for pkg=1:2:size(PackageTemp,1)
                    file_tag = ['N_',num2str(N_Kuhn),'.phi_',num2str(phi,'%.3f')];
                    
                    % Only proceed if there isn't a video tracing this
                    % particular molecule in this particular sim
                    file_name = [Folder,'/',file_tag,'.avi'];
                    if ~isfile(file_name) || Override==1
                        InputScript(pkg,Np,Ns,N,N_Kuhn,b,phi,eaStar,edStar,D,dt,kbT);
                        
                        %% Make diriectories and set filenames
                        SetDirAndFileNames;
                        DefineCompiledFileNames;
                        SizeTheDomain(Np,Separation);
                        
                        dat = load(RawDataFileName,'-mat','timesteps','Corners',...
                            'type','id','mol',...
                            'x','y','z',...
                            'btype','batom1','batom2','rx','ry','rz');
                        
                        timesteps = dat.timesteps;
                        Corners = dat.Corners;
                        Corners = Corners(1,:);
                        Lx = abs(Corners(2))*2;
                        Ly = abs(Corners(3))*2;
                        Lz = abs(Corners(4))*2;
                        types = dat.type;
                        ids = dat.id;
                        x = dat.x; y = dat.y; z = dat.z;
                        %                     rx_all = dat.rx; ry_all = dat.ry; rz_all = dat.rz;
                        bond_types = dat.btype;
                        partner1 = dat.batom1;
                        partner2 = dat.batom2;
                        rx = dat.rx;
                        ry = dat.ry;
                        rz = dat.rz;
                        mols = dat.mol;
                        
                        num_fade_steps = round(length(timesteps/20));
                        
                        clear dat
                        
                        % 1. Pick a molecule at random from those far from
                        % periodic bounds (simpler)
%                         x_tmp = x{1}; y_tmp = y{1}; z_tmp = z{1};
%                         ids_tmp = ids{1};
%                         types_tmp = types{1};
%                         mols_tmp = mols{1};
                        
                        [~,mols_center] = ...
                            FindMoleculesNearBorders(x,y,z,types,mols,...
                            Lx,Ly,Lz,N_Kuhn,b,tether_type);
                        
%                         mol_of_interest = randi(numel(unique(mols_center)));
                        
                        mol_of_interest = 21;
                        file_tag = ['N_',num2str(N_Kuhn),'.phi_',num2str(phi,'%.3f'),...
                            '.mol_',num2str(mol_of_interest)];
                        %                     file_tag = ['N_',num2str(N_Kuhn),'.phi_',num2str(phi,'%.3f')];
                        %
                        %                     % Only proceed if there isn't a video tracing this
                        %                     % particular molecule in this particular sim
                        %                     file_name = [Folder,'/',file_tag,'.avi'];
                        %                     if ~isfile(file_name) || Override==1
                        
                        % 2. ID the potential bonding neighbors of this molecule
%                         mols_neighbs = FindPotentialBondNeighbors(x,y,z,...
%                             types,mols,N_Kuhn,b,tether_type,mol_of_interest);
                        
                        v1 = VideoWriter([Folder,'/',file_tag]);
                        v1.FrameRate = 10;
                        open(v1)
                        
                        wb = waitbar(0,'Making bond exchange video');
                        
                        bonded_neighbs = [];
                        for tm=1:2:length(timesteps)
                            waitbar(tm/length(timesteps),wb);
                            
                            % 3. Isolate Data
                            
                            % Unpack the data
                            x_tmp = x{tm}; y_tmp = y{tm}; z_tmp = z{tm};
                            types_tmp = types{tm}; ids_tmp = ids{tm};
                            mols_tmp = mols{tm};
                            bondtypes_tmp = bond_types{tm};
                            p1_tmp = partner1{tm};
                            p2_tmp = partner2{tm};
                            rx_tmp = rx{tm};
                            ry_tmp = ry{tm};
                            rz_tmp = rz{tm};
                            all_pairs = [p1_tmp p2_tmp rx_tmp ry_tmp rz_tmp bondtypes_tmp];
                            
                            % Find out if the molecule of interest is
                            % bonded to another molecule
                           
                                                     
                            % Isolate data for the molecule of interest
                            x_int = x_tmp(mols_tmp==mol_of_interest);
                            y_int = y_tmp(mols_tmp==mol_of_interest);
                            z_int = z_tmp(mols_tmp==mol_of_interest);
                            types_int = types_tmp(mols_tmp==mol_of_interest);
                            ids_int = ids_tmp(mols_tmp==mol_of_interest);
                            
                            % Identify the tether of the molecule of
                            % interest and center everything about it
                            [x_tmp,y_tmp,z_tmp,x_int,y_int,z_int,...
                                x_th,y_th,z_th,id_th] = CenterAboutTether(x_int,y_int,z_int,ids_int,...
                                types_int,tether_type,x_tmp,y_tmp,z_tmp);
                            
%                             % Center everything about the tether of the
%                             % molecule of interest
%                             [x_tmp,y_tmp,z_tmp,x_int,y_int,z_int,...
%                                 x_th,y_th,z_th,id_th] = CenterAboutTether(x_int,y_int,z_int,ids_int,...
%                                 types_int,tether_type,x_tmp,y_tmp,z_tmp);
                            
                            % Find its sticker ID
                            [~,~,~,~] = FindStickerID(x_int,y_int,z_int,ids_int,...
                                types_int,sticker_type);
                            
%                             % Determine whether sticker is open or closed
%                             indx1 = find(p1_tmp==id_st);
%                             indx2 = find(p2_tmp==id_st);
%                             n_bnds_st = length(indx1)+length(indx2);
%                             
%                             if n_bnds_st>1
%                                 color = [22 71 167]/255;
%                                 color_st = 'b';
%                             else
%                                 color = [119 3 3]/255;
%                                 color_st = 'r';
%                             end
                            
                            % Starting from the tether position, propogate
                            % through the bond pairs [p1 p2] until the
                            % chain terminates. Store the data in a matrix
                            bond_indices = [find(all_pairs(:,1)==id_th);...
                                find(all_pairs(:,2)==id_th)]; % find where tether is in a bond pair
                            temp_pair = all_pairs(bond_indices,:); % define temp pair of interest
                            all_pairs(bond_indices,:) = []; % remove pair from pool
                            id_from = id_th;
                            x_chain = x_th; y_chain = y_th; z_chain = z_th;
                            ptype_chain = 1;
                            x_from = x_th; y_from = y_th; z_from = z_th;
                            rx_chain = []; ry_chain = []; rz_chain = [];
                            btype_chain = [];
                            ct = 0;
                            add_sign = 1;
                            while ~isempty(bond_indices)
                                ct = ct+1;
                                if size(temp_pair,1)==2
                                    temp_pair = temp_pair(temp_pair(:,1)==id_from,:);
                                    add_sign = 1; % will be moving backwards along chain now so subtract [rx,ry rz]
                                end
                                id_to = temp_pair(1:2); % find particle the current particle is bonded to
                                id_to(id_to==id_from) = [];
                                
                                x_to = x_tmp(ids_tmp==id_to); % find position of next particle in line
                                y_to = y_tmp(ids_tmp==id_to);
                                z_to = z_tmp(ids_tmp==id_to);
                                rx_to = x_to-x_from;
                                ry_to = y_to-y_from;
                                rz_to = z_to-z_from;
                                if rx_to<=-Lx/2
                                    rx_to = rx_to+Lx;
                                end
                                if rx_to>=Lx/2
                                    rx_to = rx_to-Lx;
                                end
                                if ry_to<=-Ly/2
                                    ry_to = ry_to+Ly;
                                end
                                if ry_to>=Ly/2
                                    ry_to = ry_to-Ly;
                                end                                
                                if rz_to<=-Lz/2
                                    rz_to = rz_to+Lz;
                                end
                                if rz_to>=Lz/2
                                    rz_to = rz_to-Lz;
                                end
                                %                                 rx_to = temp_pair(3); % find bond vector for from-to pair
%                                 ry_to = temp_pair(4);
%                                 rz_to = temp_pair(5);
                                btype_to = temp_pair(end);
                                x_to = x_from + add_sign*rx_to; % find position of next particle in line
                                y_to = y_from + add_sign*ry_to;
                                z_to = z_from + add_sign*rz_to;
                                ptype_to = types_tmp(ids_tmp==id_to);
                                
                                x_chain = cat(1,x_chain,x_to); % append particle positions to the list
                                y_chain = cat(1,y_chain,y_to);
                                z_chain = cat(1,z_chain,z_to);
                                ptype_chain = cat(1,ptype_chain,ptype_to);
%                                 if abs(rx_to)>3*b || abs(ry_to)>3*b...
%                                         || abs(rz_to)>3*b
%                                     temp_pair;
%                                 end
                                rx_chain = cat(1,rx_chain,rx_to); % append rx to list
                                ry_chain = cat(1,ry_chain,ry_to); 
                                rz_chain = cat(1,rz_chain,rz_to); 
                                btype_chain = cat(1,btype_chain,btype_to);
                                
                                bond_indices = [find(all_pairs(:,1)==id_to);...
                                    find(all_pairs(:,2)==id_to)]; % find where next particle in-line is in another bond pair
                                if ~isempty(bond_indices)
                                    temp_pair = all_pairs(bond_indices,:);
                                    all_pairs(bond_indices,:) = []; % remove pair from pool so nothing is double counted
                                end
                                id_from = id_to;
                                x_from = x_to;
                                y_from = y_to;
                                z_from = z_to;
                            end
                            
                            % Plot the chain
                            
                            % ID Color of chain (if bonded, use blue)
                            color = [119 3 3]/255;
                            color_st = 'r';
                            if ~isempty(btype_chain(btype_chain==2))
                                color = [22 71 167]/255;
                                color_st = 'b';
                            end
                            figure(1); clf; hold on
                            % Identify which chain has been bonded to
                            if length(x_chain(ptype_chain==tether_type))>1
                                id_th_to = id_to;
                                x_bonded = x_chain(ptype_chain==tether_type);
                                y_bonded = y_chain(ptype_chain==tether_type);
                                z_bonded = z_chain(ptype_chain==tether_type);
                                rem0 = find(x_bonded==0 & y_bonded==0 & z_bonded==0);
                                x_bonded(rem0) = [];
                                y_bonded(rem0) = [];
                                z_bonded(rem0) = [];
                                bonded_neighbs_tmp = ...
                                    [id_th_to x_bonded y_bonded z_bonded tm];
                                bonded_neighbs = cat(1,bonded_neighbs,...
                                    bonded_neighbs_tmp);
                                if size(bonded_neighbs,1)>1
                                    bonded_neighbs = sortrows(bonded_neighbs,5);
                                end
                                % clean up bonded neighbs
                                unique_neighbs = unique(bonded_neighbs(:,1));
                                rem_indx_all = [];
                                for un=1:length(unique_neighbs)
                                    neighb_tmp = unique_neighbs(un);
                                    indx1 = find(bonded_neighbs(:,1)==neighb_tmp,1,'last'); % save most recent instance
                                    rem_indices = find(bonded_neighbs(:,1)==neighb_tmp);
                                    rem_indices = rem_indices(rem_indices~=indx1);
                                    rem_indx_all = cat(1,rem_indx_all,rem_indices);
                                end
                                if ~isempty(rem_indx_all)
                                    bonded_neighbs(rem_indx_all,:) = [];
                                end
                            end
                            sct = scatter3(x_chain(ptype_chain==tether_type),...
                                y_chain(ptype_chain==tether_type),...
                                z_chain(ptype_chain==tether_type),'filled');
                            sct.MarkerFaceColor = [0.5 0.5 0.5];
                            sct.MarkerEdgeColor = 'k';
                            sct.SizeData = marker_size*2;
                            if BeadSpringOrMeso==0
                                sct = scatter3(x_chain(ptype_chain==intermed_type),...
                                    y_chain(ptype_chain==intermed_type),...
                                    z_chain(ptype_chain==intermed_type),'filled');
                                sct.MarkerFaceColor = color;
                                sct.SizeData = marker_size;
                            end
                            sct = scatter3(x_chain(ptype_chain==sticker_type),...
                                y_chain(ptype_chain==sticker_type),...
                                z_chain(ptype_chain==sticker_type),'filled');
                            sct.MarkerFaceColor = color_st;
                            sct.MarkerEdgeColor = 'k';
                            sct.SizeData = marker_size*2;
                          
                            % Plot the bonds
                            x_from = x_th;
                            y_from = y_th;
                            z_from = z_th;
                            for bnd_no=1:length(rx_chain)
                                if ~isempty(btype_chain(btype_chain==2))
                                    color_tmp = color;
                                else
                                    color_tmp = color*bnd_no/length(rx_chain);
                                end
                                x_to = x_from+rx_chain(bnd_no);
                                y_to = y_from+ry_chain(bnd_no);
                                z_to = z_from+rz_chain(bnd_no);
                                if btype_chain(bnd_no)==1
                                    line_style = '-';
                                else
                                    color_tmp = 'r';
                                    line_style = ':';
                                end
                                
                                pl = plot3([x_from x_to],[y_from y_to],[z_from z_to]);
                                pl.LineStyle = line_style;
                                pl.Color = color_tmp; 
                                pl.LineWidth = line_width;
                                x_from = x_to;
                                y_from = y_to;
                                z_from = z_to;
                            end
                            
                            % Plot previously bonded_neighbs
                            for bn=1:size(bonded_neighbs,1)
                                sct = scatter3(bonded_neighbs(bn,2),...
                                    bonded_neighbs(bn,3),...
                                    bonded_neighbs(bn,4),...
                                    'filled');
                                alpha = (tm-bonded_neighbs(bn,end))/num_fade_steps;
                                alpha(alpha>1) = 1;
                                alpha = 1-alpha;
                                sct.MarkerFaceColor = [0.5 0.5 0.5];
                                sct.MarkerEdgeColor = 'k';
                                sct.MarkerFaceAlpha = alpha;
                                sct.MarkerEdgeAlpha = alpha;
                                sct.SizeData = marker_size*1.5;
                            end
                            
                            % plot box
                            plot3([-Lx/2 Lx/2],[-Ly/2 -Ly/2],[-Lz/2 -Lz/2],'k--')
                            plot3([-Lx/2 Lx/2],[+Ly/2 +Ly/2],[+Lz/2 +Lz/2],'k--')
                            plot3([-Lx/2 Lx/2],[+Ly/2 +Ly/2],[-Lz/2 -Lz/2],'k--')
                            plot3([-Lx/2 Lx/2],[-Ly/2 -Ly/2],[+Lz/2 +Lz/2],'k--')
                            
                            plot3([-Lx/2 -Lx/2],[-Ly/2 +Ly/2],[-Lz/2 -Lz/2],'k--')
                            plot3([-Lx/2 -Lx/2],[-Ly/2 +Ly/2],[+Lz/2 +Lz/2],'k--')
                            plot3([+Lx/2 +Lx/2],[-Ly/2 +Ly/2],[-Lz/2 -Lz/2],'k--')
                            plot3([+Lx/2 +Lx/2],[-Ly/2 +Ly/2],[+Lz/2 +Lz/2],'k--')
                            
                            plot3([-Lx/2 -Lx/2],[-Ly/2 -Ly/2],[-Lz/2 +Lz/2],'k--')
                            plot3([-Lx/2 -Lx/2],[+Ly/2 +Ly/2],[-Lz/2 +Lz/2],'k--')
                            plot3([+Lx/2 +Lx/2],[-Ly/2 -Ly/2],[-Lz/2 +Lz/2],'k--')
                            plot3([+Lx/2 +Lx/2],[+Ly/2 +Ly/2],[-Lz/2 +Lz/2],'k--')
                            
                            view(60,30)
                            pbaspect([1 1 1])
                            xlim([-Lx/2 Lx/2])
                            ylim([-Lx/2 Lx/2])
                            zlim([-Lx/2 Lx/2])
                            daspect([1 1 1]) 
                            set(gcf,'color','w')
                            xticks([])
                            yticks([])
                            zticks([])
                            set(gcf,'Position',[1000 50 600 600])
                            
                            title(['$t = ',num2str(timesteps(tm)*dt/1e-9,'%.2f'),...
                                '$ ns'],'FontSize',10,'Interpreter','latex')
                            
                            F1 = getframe(gcf);
                            writeVideo(v1,F1);
                            mov(tm) = F1;
                        end
                        close(v1)
                        close(wb)
                        
                        clear timesteps Corners Lx Ly Lz types ids x y z...
                            bond_types partner1 partner2 rx ry rz mols
                    end

                    
                    

                    
                    % 4. Plot the central molecule in color
                    
                    % If chain is detachd then it is maroon
                    
                    % If chain is attached then it is cyan
                    
                    % 5. Plot the neighbors in grayscale
                    
                end
             end
            end
           end
          end
         end
        end
       end
      end
     end
    end
    
end


    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rs = ExtractEndToEnd(types,id,mol,x,y,z,bond_types,batom1,batom2,...
    rx,ry,rz,Np,i,samp,Corners)

global BeadSpringOrMeso bondtypes Nt Ns meso_pairs N_Kuhn b...
%      Lx Ly Lz meso_pairs0
 
CheckFig = 0;

L = abs(Corners(3)-Corners(2));
W = abs(Corners(5)-Corners(4));
H = abs(Corners(7)-Corners(6));

if BeadSpringOrMeso==1
    rx = rx(bond_types~=2);
    ry = ry(bond_types~=2);
    rz = rz(bond_types~=2);
    rs = [rx ry rz];
    
    meso_pairs = [batom1(bond_types==1), batom2(bond_types==1)];
else
    if isempty(meso_pairs) 
        % Initialize all vectors
        n_tethers = Nt-1;
        n_stickers = Nt*Ns; % per molecule
        Nsegments = Np*(n_tethers+n_stickers);
        rx_out = zeros(Nsegments,1);
        ry_out = zeros(Nsegments,1);
        rz_out = zeros(Nsegments,1);

        batom1_out = zeros(Nsegments,1);
        batom2_out = zeros(Nsegments,1);

        % Define segments of side-branches
        ct = 0;
        all_pairs = [batom1 batom2];
        all_pairs(bond_types==2,:) = [];
        mols = unique(mol);
        for i=1:length(mols)
            molecule = mols(i);     % define the molecule
            temp_types = types(mol==molecule);      % define the particle types belonging to the molecule
            temp_id = id(mol==molecule);            % find particles belonging to the molecule
            [members,sort_indx] = sort(temp_id);    % sort particles belonging to the molecule
            temp_types = temp_types(sort_indx);     % sort the particle types

            % Isolate pairs belonging to the molecule
            indx1 = find(ismember(all_pairs(:,1),members));     % find all bound pairs belonging to the molecule
            indx2 = find(ismember(all_pairs(:,2),members));
            indx_mps = intersect(indx1,indx2);
            mol_pairs = all_pairs(indx_mps,:);      % bonds belonging to the molecule (i.e., Kuhn segments)
            mol_rx = rx(indx_mps,:);
            mol_ry = ry(indx_mps,:);
            mol_rz = rz(indx_mps,:);

            tethers = members(1:Nt);
            sticker_rng = Nt+(1:n_stickers);
            stickers = members(sticker_rng);
            others = members;
            others(ismember(others,tethers)) = [];
            others(ismember(others,stickers)) = [];
            if isempty(temp_types(sticker_rng)==1)
                warning('Incorrect stickers selected')
            end

            % Debugging fig
            if CheckFig==1
                figure(100); clf; hold on
                tether_pos = [x(ismember(id,tethers)),...
                    y(ismember(id,tethers)),...
                    z(ismember(id,tethers))];
                sticker_pos = [x(ismember(id,stickers)),...
                    y(ismember(id,stickers)),...
                    z(ismember(id,stickers))];
                other_pos = [x(ismember(id,others)),...
                    y(ismember(id,others)),...
                    z(ismember(id,others))];
                s = scatter3(tether_pos(:,1),tether_pos(:,2),tether_pos(:,3),'filled');
                s.MarkerFaceColor = 'c'; s.SizeData = 10;
                s = scatter3(sticker_pos(:,1),sticker_pos(:,2),sticker_pos(:,3),'filled');
                s.MarkerFaceColor = 'r'; s.SizeData = 10;
                s = scatter3(other_pos(:,1),other_pos(:,2),other_pos(:,3),'filled');
                s.MarkerFaceColor = [0.5 0.5 0.5]; s.SizeData = 10;
                view([30,45])
            end
            % Define summated vector of bead-springs between adjacent tethers
            mol_pairs = sort(mol_pairs,2);  % sort the Kuhn segments' end particles from left-to-right
            [mol_pairs,unq_indx] = unique(mol_pairs,'rows');  % This should be unnecessary
            mol_rx = mol_rx(unq_indx);
            mol_ry = mol_ry(unq_indx);
            mol_rz = mol_rz(unq_indx);
            flag = 1;
            while ~isempty(mol_pairs)  %while there remains molecule pairs
                % ID first tether in line
                if flag  %if its the first particle of the chain segment
                    indx = find(ismember(mol_pairs(:,1),tethers),1,'first'); % find the first instance of a tether in the Kuhn segments list
                    from = mol_pairs(indx,1);   % Going from...
                    to = mol_pairs(indx,2);     % ... to
                    from0 = from;

                    [rx_temp,ry_temp,rz_temp] =  ...
                        Calcdr(from,to,x,y,z,id,mol_rx,mol_ry,mol_rz,indx,L,W,H);

                    flag = 0; add_vecs = 0;
                    rxs = []; rys = []; rzs = [];   % Reset the rxs, rys and rzs for this chain segment
                    ct = ct+1;
                else
                    indx3 = find(ismember(mol_pairs(:,1),from),1);
                    indx4 = find(ismember(mol_pairs(:,2),from),1);
                    indx = [indx3;indx4];
                    to  = mol_pairs(indx,mol_pairs(indx,:)~=from);

                    [rx_temp,ry_temp,rz_temp] =  ...
                        Calcdr(from,to,x,y,z,id,mol_rx,mol_ry,mol_rz,indx,L,W,H);

                    if ismember(to,tethers) || ismember(to,stickers)
                        flag = 1; add_vecs = 1;
                    end
                end

                % Append vector components to stored
                rxs = cat(1,rxs,rx_temp);
                rys = cat(1,rys,ry_temp);
                rzs = cat(1,rzs,rz_temp);

                if add_vecs
                    rx_out(ct) = sum(rxs);
                    ry_out(ct) = sum(rys);
                    rz_out(ct) = sum(rzs);

                    rx_check = x(id==to)-x(id==from0);
                    ry_check = y(id==to)-y(id==from0);
                    rz_check = z(id==to)-z(id==from0);

                    bound_fact = 2; % Bigger bound factor means more 
                                    % sensitively going to treat springs as
                                    % if they span the periodic boundaries

                    if rx_check<-L/bound_fact
                        rx_check = rx_check+L;
                    elseif rx_check>L/bound_fact
                        rx_check = rx_check-L;
                    end
%                     rx_check(rx_check<-L/bound_fact) = rx_check(rx_check<-L/bound_fact)+L;
%                     rx_check(rx_check>L/bound_fact) = rx_check(rx_check>L/bound_fact)-L;

                    if ry_check<-W/bound_fact
                        ry_check = ry_check+W;
                    elseif ry_check>W/bound_fact
                        ry_check = ry_check-W;
                    end
%                     ry_check(ry_check<-W/bound_fact) = ry_check(ry_check<-W/bound_fact)+W;
%                     ry_check(ry_check>W/bound_fact) = ry_check(ry_check>W/bound_fact)-W;

                    if rz_check<-H/bound_fact
                        rz_check = rz_check+H;
                    elseif rz_check>H/bound_fact
                        rz_check = rz_check-H;
                    end
%                      rz_check(rz_check<-H/bound_fact) = rz_check(rz_check<-H/bound_fact)+H;
%                     rz_check(rz_check>H/bound_fact) = rz_check(rz_check>H/bound_fact)-H;

                    if round(rx_check,2)~=round(rx_out(ct),2) || ...
                            round(ry_check,2)~=round(ry_out(ct),2) || ...
                            round(rz_check,2)~=round(rz_out(ct),2)
                        warning('Mismatch in values for two methods - check round off')
                    end

                    batom1_out(ct) = from0;
                    batom2_out(ct) = to;
                end

                if CheckFig==1
                    x1 = x(id==from);
                    y1 = y(id==from);
                    z1 = z(id==from);
                    x2 = x1+rxs(end);
                    y2 = y1+rys(end);
                    z2 = z1+rzs(end);
                    plot3([x1 x2],[y1 y2],[z1 z2],'k-')
                    if add_vecs
                        x4 = x(id==to);
                        y4 = y(id==to);
                        z4 = z(id==to);
                        x3 = x4-rx_out(ct);
                        y3 = y4-ry_out(ct);
                        z3 = z4-rz_out(ct);
                        plot3([x4 x3],[y4 y3],[z4 z3],'r')
                    end
                end

                % Initialize reference bead for next iteration
                from = to;      % the new particle 'going from' is the old 'going to'
                mol_pairs(indx,:) = [];
                mol_rx(indx,:) = [];
                mol_ry(indx,:) = [];
                mol_rz(indx,:) = [];
            end
        end
        rs = [rx_out ry_out rz_out];
        meso_pairs = [batom1_out batom2_out];
    else
        [~,indx] = sort(id,1);
        x_sort = x(indx); 
        y_sort = y(indx); 
        z_sort = z(indx); 
        x_from = x_sort(meso_pairs(:,1));
        y_from = y_sort(meso_pairs(:,1));
        z_from = z_sort(meso_pairs(:,1));
        x_to = x_sort(meso_pairs(:,2));
        y_to = y_sort(meso_pairs(:,2));
        z_to = z_sort(meso_pairs(:,2));
        rx_out = x_to-x_from;
        ry_out = y_to-y_from;
        rz_out = z_to-z_from;

        % adjust for periodic boundaries
        indx_less = find(rx_out<-L/2);
        indx_more = find(rx_out>L/2);
        if ~isempty(indx_less)
            rx_out(indx_less) = rx_out(indx_less)+L;
        end
        if ~isempty(indx_more)
            rx_out(indx_more) = rx_out(indx_more)-L;
        end

        indx_less = find(ry_out<-W/2);
        indx_more = find(ry_out>W/2);
        if ~isempty(indx_less)
            ry_out(indx_less) = ry_out(indx_less)+W;
        end
        if ~isempty(indx_more)
            ry_out(indx_more) = ry_out(indx_more)-W;
        end

        indx_less = find(rz_out<-H/2);
        indx_more = find(rz_out>H/2);
        if ~isempty(indx_less)
            rz_out(indx_less) = rz_out(indx_less)+H;
        end
        if ~isempty(indx_more)
            rz_out(indx_more) = rz_out(indx_more)-H;
        end

        if ~isempty(rx_out(abs(rx_out)>N_Kuhn*b))
            warning('spring too long in x-direction');
        end
        if ~isempty(ry_out(abs(ry_out)>N_Kuhn*b))
            warning('spring too long in y-direction');
        end
        if ~isempty(rz_out(abs(rz_out)>N_Kuhn*b))
            warning('spring too long in z-direction');
        end
% 
%         rx_out(rx_out<-L/2) = rx_out(rx_out<-L/2)+L;
%         rx_out(rx_out>L/2) = rx_out(rx_out>L/2)-L;
% 
%         ry_out(ry_out<-W/2) = ry_out(ry_out<-W/2)+W;
%         ry_out(ry_out>W/2) = ry_out(ry_out>W/2)-W;    
% 
%         rz_out(rz_out<-H/2) = rz_out(rz_out<-H/2)+H;
%         rz_out(rz_out>H/2) = rz_out(rz_out>H/2)-H;
% 
%         rx_check = rx_out(abs(rx_out)>1.1*N_Kuhn*b);
%         ry_check = ry_out(abs(ry_out)>1.1*N_Kuhn*b);
%         
%         
%         rz_check = rz_out(abs(rz_out)>1.1*N_Kuhn*b);
%         if ~isempty(rx_check) || ~isempty(ry_check) || ~isempty(rz_check) 
%             warning('One of the springs is too long')
%         end
        rs = [rx_out ry_out rz_out];
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SizeTheDomain(Np,Separation)

global Lx Ly Lz PairsPerEdge c0

PairsPerEdge = round(Np^(1/3));

Lx = PairsPerEdge*Separation;
Ly = Lx;
Lz = Lx;

c0 = Np/(Lx*Ly*Lz);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AssembleEnsembleDynamicsData(Override,Package,...
    damps,N_Kuhns,Np,Ns,kbT,phis,ka_in,kd_in,f0,b,N,PlotDynwrtTime)

global N_Kuhn EnsembleDynamicsFileName AllDynamicsFileName...
    conc_all_meso...
    ka_ens_meso kd_ens_meso kaSE_ens_meso kdSE_ens_meso...
    ka1_ens_meso ka1SE_ens_meso...
    Na_ens_meso Nd_ens_meso NaSE_ens_meso NdSE_ens_meso...
    fa_ens_meso fd_ens_meso faSE_ens_meso fdSE_ens_meso...
    tau_a_all_meso tau_d_all_meso...
    tau_exc_all_meso tau_rpt_all_meso tau_adj_all_meso...
    ka_all_meso ka1_all_meso kd_all_meso...
    Na_all_meso Nd_all_meso...
    fa_all_meso fd_all_meso...
    tau_a_ens_meso tau_aSE_ens_meso...
    tau_d_ens_meso tau_dSE_ens_meso...
    tau_exc_ens_meso tau_excSE_ens_meso ...
    tau_rpt_ens_meso tau_rptSE_ens_meso ...
    tau_adj_ens_meso tau_adjSE_ens_meso ...
    tau_rnm_ens_meso tau_rnmSE_ens_meso...
    kexc_ens_meso kexcSE_ens_meso ...
    krpt_ens_meso krptSE_ens_meso...
    conc_all_bead...
    ka_ens_bead kd_ens_bead kaSE_ens_bead kdSE_ens_bead...
    ka1_ens_bead ka1SE_ens_bead...
    Na_ens_bead Nd_ens_bead NaSE_ens_bead NdSE_ens_bead...
    fa_ens_bead fd_ens_bead faSE_ens_bead fdSE_ens_bead...
    tau_a_all_bead tau_d_all_bead...
    tau_exc_all_bead tau_rpt_all_bead tau_adj_all_bead...
    ka_all_bead ka1_all_bead kd_all_bead...
    Na_all_bead Nd_all_bead...
    fa_all_bead fd_all_bead...
    tau_a_ens_bead tau_aSE_ens_bead...
    tau_d_ens_bead tau_dSE_ens_bead...
    tau_exc_ens_bead tau_excSE_ens_bead ...
    tau_rpt_ens_bead tau_rptSE_ens_bead ...
    tau_adj_ens_bead tau_adjSE_ens_bead ...
    tau_rnm_ens_bead tau_rnmSE_ens_bead...
    kexc_ens_bead kexcSE_ens_bead ...
    krpt_ens_bead krptSE_ens_bead...
    conc_all_thry...
    ka_ens_thry kd_ens_thry kaSE_ens_thry kdSE_ens_thry...
    ka1_ens_thry ka1SE_ens_thry...
    Na_ens_thry Nd_ens_thry NaSE_ens_thry NdSE_ens_thry...
    fa_ens_thry fd_ens_thry faSE_ens_thry fdSE_ens_thry...
    tau_a_all_thry tau_d_all_thry...
    tau_exc_all_thry tau_rpt_all_thry tau_adj_all_thry...
    ka_all_thry ka1_all_thry kd_all_thry...
    Na_all_thry Nd_all_thry...
    fa_all_thry fd_all_thry...
    tau_a_ens_thry tau_aSE_ens_thry...
    tau_d_ens_thry tau_dSE_ens_thry...
    tau_exc_ens_thry tau_excSE_ens_thry ...
    tau_rpt_ens_thry tau_rptSE_ens_thry ...
    tau_adj_ens_thry tau_adjSE_ens_thry ...
    tau_rnm_ens_thry tau_rnmSE_ens_thry...
    kexc_ens_thry kexcSE_ens_thry ...
    krpt_ens_thry krptSE_ens_thry...
    AppendScalingData...
    LengthConversion DamperConversion...
%     tau_a_bead tau_d_bead tau_exc_bead tau_adj_bead tau_rpt_bead...
%     tau_a_meso tau_d_meso tau_exc_meso tau_adj_meso tau_rpt_meso...

Folder = 'Output Plots';
if ~isfolder(Folder)
    mkdir(Folder)
end

DefineAssembledFileNames(ka_in);

% MSD wrt time for each damper while sweeping N
for dm = 1:length(damps)
    for nk = 1:length(N_Kuhns)
        N_Kuhn = N_Kuhns(nk);

        PackageTemp = Package(ismember(Package(:,2),Np),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,3),Ns),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,4),ka_in),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,5),kd_in),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,6),f0),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,7),N_Kuhn),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,8),b),:);
        PackageTemp = PackageTemp(ismember(PackageTemp(:,11),kbT),:);
        
        if PlotDynwrtTime==1
            figure(1); clf; hold on;
            figure(2); clf; hold on;
            figure(3); clf; hold on;
            figure(4); clf; hold on;
            
            n_clr = length(phis);
            colors = [(linspace(0.15,1,n_clr))',...
                (linspace(0.15,0.75,n_clr))',...
                flipud((linspace(0,0.15,n_clr))')];
            
            fig_folder = 'Output Plots/Kinetics vs Time';
            if ~isfolder(fig_folder)
                mkdir(fig_folder)
            end
        end
        for sp = 1:size(PackageTemp,1) % loop over the separation distances
            Sample = PackageTemp(sp,1);
            phi = phis(sp);
            
            [~,D,~,~] = DefineTimeScale(b,LengthConversion,...
                DamperConversion,1);
             
            eaStar = -log(ka_in*(b*LengthConversion)^2/D);
            edStar = -log(kd_in*(b*LengthConversion)^2/D);
            
            if sp==1 && nk==1
                InitializeEnsembleDynamics(N_Kuhns,size(PackageTemp,1));
            end

            %% Assemble the kinetics data by separation distance
            if ~isempty(PackageTemp) && (~isfile(EnsembleDynamicsFileName) ||...
                    Override)

                % Pull in the Bead-spring data
                BSOM = 0;
                [~,D,tau0,dtFact] = DefineTimeScale(b,LengthConversion,...
                    DamperConversion,BSOM);
                dt = tau0/dtFact;
                [time_bead,conc_all_bead(sp,nk),ka_bead,ka1_bead,kd_bead,...
                    kexc_bead,krpt_bead,....
                    Na_bead,Nd_bead,fa_bead,fd_bead,...
                    tau_a_bead,tau_d_bead,...
                    tau_exc_bead,tau_rpt_bead,tau_adj_bead,...
                    ka_ens_bead(sp,nk),kaSE_ens_bead(sp,nk),...
                    ka1_ens_bead(sp,nk),ka1SE_ens_bead(sp,nk),...
                    kd_ens_bead(sp,nk),kdSE_ens_bead(sp,nk),...
                    Na_ens_bead(sp,nk),NaSE_ens_bead(sp,nk),...
                    Nd_ens_bead(sp,nk),NdSE_ens_bead(sp,nk),...
                    fa_ens_bead(sp,nk),faSE_ens_bead(sp,nk),...
                    fd_ens_bead(sp,nk),fdSE_ens_bead(sp,nk),...
                    tau_a_ens_bead(sp,nk),tau_aSE_ens_bead(sp,nk),...
                    tau_d_ens_bead(sp,nk),tau_dSE_ens_bead(sp,nk),...
                    tau_exc_ens_bead(sp,nk),tau_excSE_ens_bead(sp,nk),...
                    tau_rpt_ens_bead(sp,nk),tau_rptSE_ens_bead(sp,nk),...
                    tau_adj_ens_bead(sp,nk),tau_adjSE_ens_bead(sp,nk),...
                    kexc_ens_bead(sp,nk),kexcSE_ens_bead(sp,nk),...
                    krpt_ens_bead(sp,nk),krptSE_ens_bead(sp,nk)] =...
                    CalculateEnsembleDynamics(BSOM,Sample,Np,Ns,N,phi,eaStar,edStar,...
                        dt,N_Kuhn,b,D,kbT);
                    

                % Pull in the Mesoscale data
                BSOM = 1;
                [~,D,tau0,dtFact] = DefineTimeScale(b,LengthConversion,...
                    DamperConversion,BSOM);
                dt = tau0/dtFact;
                [time_meso,conc_all_meso(sp,nk),ka_meso,ka1_meso,kd_meso,...
                    kexc_meso,krpt_meso,....
                    Na_meso,Nd_meso,fa_meso,fd_meso,...
                    tau_a_meso,tau_d_meso,...
                    tau_exc_meso,tau_rpt_meso,tau_adj_meso,...
                    ka_ens_meso(sp,nk),kaSE_ens_meso(sp,nk),...
                    ka1_ens_meso(sp,nk),ka1SE_ens_meso(sp,nk),...
                    kd_ens_meso(sp,nk),kdSE_ens_meso(sp,nk),...
                    Na_ens_meso(sp,nk),NaSE_ens_meso(sp,nk),...
                    Nd_ens_meso(sp,nk),NdSE_ens_meso(sp,nk),...
                    fa_ens_meso(sp,nk),faSE_ens_meso(sp,nk),...
                    fd_ens_meso(sp,nk),fdSE_ens_meso(sp,nk),...
                    tau_a_ens_meso(sp,nk),tau_aSE_ens_meso(sp,nk),...
                    tau_d_ens_meso(sp,nk),tau_dSE_ens_meso(sp,nk),...
                    tau_exc_ens_meso(sp,nk),tau_excSE_ens_meso(sp,nk),...
                    tau_rpt_ens_meso(sp,nk),tau_rptSE_ens_meso(sp,nk),...
                    tau_adj_ens_meso(sp,nk),tau_adjSE_ens_meso(sp,nk),...
                    kexc_ens_meso(sp,nk),kexcSE_ens_meso(sp,nk),...
                    krpt_ens_meso(sp,nk),krptSE_ens_meso(sp,nk)] =...
                    CalculateEnsembleDynamics(BSOM,Sample,Np,Ns,N,phi,eaStar,edStar,...
                        dt,N_Kuhn,b,D,kbT);

                % Plot dynamics wrt time
                if PlotDynwrtTime==1
                    color = colors(sp,:);
                    
                    figure(1)
                    save_tag = 'k_d';
                    y_bs = kd_bead;
                    y_ms = kd_meso; 
                    mm_pcnt = 20;
                    PlotKineticswrtTime(time_bead,time_meso,y_bs,y_ms,...
                        color,tau0,fig_folder,phis,N_Kuhn,sp,save_tag,mm_pcnt)
                    ylim([0 2e-3])

                    figure(2)
                    save_tag = 'k_a';
                    y_bs = ka_bead;
                    y_ms = ka_meso;                    
                    PlotKineticswrtTime(time_bead,time_meso,y_bs,y_ms,...
                        color,tau0,fig_folder,phis,N_Kuhn,sp,save_tag,mm_pcnt)
                    ylim([0 2e-2])
                    
                    figure(3)
                    save_tag = 'k_{exc}';
                    y_bs = kexc_bead;
                    y_ms = kexc_meso;                    
                    PlotKineticswrtTime(time_bead,time_meso,y_bs,y_ms,...
                        color,tau0,fig_folder,phis,N_Kuhn,sp,save_tag,mm_pcnt)
                    ylim([0 2e-2])

                    figure(4)
                    save_tag = 'k_{rpt}';
                    y_bs = krpt_bead;
                    y_ms = krpt_meso;                    
                    PlotKineticswrtTime(time_bead,time_meso,y_bs,y_ms,...
                        color,tau0,fig_folder,phis,N_Kuhn,sp,save_tag,mm_pcnt)                    
                    ylim([0 2e-2])
                end
                
                if AppendScalingData
                    % Pull in the scaling theory data
                    BSOM = 2;
                    [time_thry,conc_all_thry(sp,nk),ka_thry,ka1_thry,kd_thry,...
                        kexc_thry,krpt_thry,....
                        Na_thry,Nd_thry,fa_thry,fd_thry,...
                        tau_a_thry,tau_d_thry,...
                        tau_exc_thry,tau_rpt_thry,tau_adj_thry,...
                        ka_ens_thry(sp,nk),kaSE_ens_thry(sp,nk),...
                        ka1_ens_thry(sp,nk),ka1SE_ens_thry(sp,nk),...
                        kd_ens_thry(sp,nk),kdSE_ens_thry(sp,nk),...
                        Na_ens_thry(sp,nk),NaSE_ens_thry(sp,nk),...
                        Nd_ens_thry(sp,nk),NdSE_ens_thry(sp,nk),...
                        fa_ens_thry(sp,nk),faSE_ens_thry(sp,nk),...
                        fd_ens_thry(sp,nk),fdSE_ens_thry(sp,nk),...
                        tau_a_ens_thry(sp,nk),tau_aSE_ens_thry(sp,nk),...
                        tau_d_ens_thry(sp,nk),tau_dSE_ens_thry(sp,nk),...
                        tau_exc_ens_thry(sp,nk),tau_excSE_ens_thry(sp,nk),...
                        tau_rpt_ens_thry(sp,nk),tau_rptSE_ens_thry(sp,nk),...
                        tau_adj_ens_thry(sp,nk),tau_adjSE_ens_thry(sp,nk),...
                        kexc_ens_thry(sp,nk),kexcSE_ens_thry(sp,nk),...
                        krpt_ens_thry(sp,nk),krptSE_ens_thry(sp,nk)] =...
                        CalculateEnsembleDynamics(BSOM,Sample,Np,Ns,N,phi,eaStar,edStar,...
                            dt,N_Kuhn,b,D,kbT);
                end

                if sp==1 && nk==1
                    InitializeAllDynamics(time_bead,time_meso,size(PackageTemp,1));
                end

                % make them all 10x longer than initial case
                
                % fill empties with NaNs
                rng = (1:length(time_bead))';

                ka_all_bead(rng,sp,nk) = ka_bead;
                ka1_all_bead(rng,sp,nk) = ka1_bead;
                kd_all_bead(rng,sp,nk) = kd_bead;
                Na_all_bead(rng,sp,nk) = Na_bead;
                Nd_all_bead(rng,sp,nk) = Nd_bead;
                fa_all_bead(rng,sp,nk) = fa_bead;
                fd_all_bead(rng,sp,nk) = fd_bead;
                
                if length(rng)>size(ka_all_bead,1)
                    ka_all_bead(rng(end)+1:end,sp,nk) = NaN;
                    ka1_all_bead(rng(end)+1:end,sp,nk) = NaN;
                    kd_all_bead(rng(end)+1:end,sp,nk) = NaN;
                    Na_all_bead(rng(end)+1:end,sp,nk) = NaN;
                    Nd_all_bead(rng(end)+1:end,sp,nk) = NaN;
                    fa_all_bead(rng(end)+1:end,sp,nk) = NaN;
                    fd_all_bead(rng(end)+1:end,sp,nk) = NaN;
                end
                
                tau_a_all_bead(1:length(tau_a_bead),sp,nk) = tau_a_bead;
                tau_d_all_bead(1:length(tau_d_bead),sp,nk) = tau_d_bead;
                tau_exc_all_bead(1:length(tau_exc_bead),sp,nk) = tau_exc_bead;
                tau_rpt_all_bead(1:length(tau_rpt_bead),sp,nk) = tau_rpt_bead;
                tau_adj_all_bead(1:length(tau_adj_bead),sp,nk) = tau_adj_bead;


                rng = (1:length(time_meso))';

                ka_all_meso(rng,sp,nk) = ka_meso;
                ka1_all_meso(rng,sp,nk) = ka1_meso;
                kd_all_meso(rng,sp,nk) = kd_meso;
                Na_all_meso(rng,sp,nk) = Na_meso;
                Nd_all_meso(rng,sp,nk) = Nd_meso;
                fa_all_meso(rng,sp,nk) = fa_meso;
                fd_all_meso(rng,sp,nk) = fd_meso;  
                 
                if length(rng)>size(ka_all_meso,1)
                    ka_all_meso(rng(end)+1:end,sp,nk) = NaN;
                    ka1_all_meso(rng(end)+1:end,sp,nk) = NaN;
                    kd_all_meso(rng(end)+1:end,sp,nk) = NaN;
                    Na_all_meso(rng(end)+1:end,sp,nk) = NaN;
                    Nd_all_meso(rng(end)+1:end,sp,nk) = NaN;
                    fa_all_meso(rng(end)+1:end,sp,nk) = NaN;
                    fd_all_meso(rng(end)+1:end,sp,nk) = NaN;
                end
                
                tau_a_all_meso(1:length(tau_a_meso),sp,nk) = tau_a_meso;
                tau_d_all_meso(1:length(tau_d_meso),sp,nk) = tau_d_meso;
                tau_exc_all_meso(1:length(tau_exc_meso),sp,nk) = tau_exc_meso;
                tau_rpt_all_meso(1:length(tau_rpt_meso),sp,nk) = tau_rpt_meso;
                tau_adj_all_meso(1:length(tau_adj_meso),sp,nk) = tau_adj_meso;

                if AppendScalingData
                    rng = (1:length(time_thry))';
                    ka_all_thry(rng,sp,nk) = ka_thry;
                    ka1_all_thry(rng,sp,nk) = ka1_thry;
                    kd_all_thry(rng,sp,nk) = kd_thry;
                    Na_all_thry(rng,sp,nk) = Na_thry;
                    Nd_all_thry(rng,sp,nk) = Nd_thry;
                    fa_all_thry(rng,sp,nk) = fa_thry;
                    fd_all_thry(rng,sp,nk) = fd_thry;

                    tau_a_all_thry(1:length(tau_a_thry),sp,nk) = tau_a_thry;
                    tau_d_all_thry(1:length(tau_d_thry),sp,nk) = tau_d_thry;
                    tau_exc_all_thry(1:length(tau_exc_thry),sp,nk) = tau_exc_thry;
                    tau_rpt_all_thry(1:length(tau_rpt_thry),sp,nk) = tau_rpt_thry;
                    tau_adj_all_thry(1:length(tau_adj_thry),sp,nk) = tau_adj_thry;
                end
            end
        end
    end
end

if ~isfile(EnsembleDynamicsFileName) || Override
    EnsembleDynamics.ka_bead = ka_ens_bead;
    EnsembleDynamics.kaSE_bead = kaSE_ens_bead;
    EnsembleDynamics.ka1_bead = ka1_ens_bead;
    EnsembleDynamics.ka1SE_bead = ka1SE_ens_bead;
    EnsembleDynamics.kd_bead = kd_ens_bead;
    EnsembleDynamics.kdSE_bead = kdSE_ens_bead;
    EnsembleDynamics.Na_bead = Na_ens_bead;
    EnsembleDynamics.NaSE_bead = NaSE_ens_bead;
    EnsembleDynamics.Nd_bead = Nd_ens_bead;
    EnsembleDynamics.NdSE_bead = NdSE_ens_bead;
    EnsembleDynamics.fa_bead = fa_ens_bead;
    EnsembleDynamics.faSE_bead = faSE_ens_bead;
    EnsembleDynamics.fd_bead = fd_ens_bead;
    EnsembleDynamics.fdSE_bead = fdSE_ens_bead;

    EnsembleDynamics.tau_a_ens_bead = tau_a_ens_bead;
    EnsembleDynamics.tau_d_ens_bead = tau_d_ens_bead;
    EnsembleDynamics.tau_exc_ens_bead = tau_exc_ens_bead;
    EnsembleDynamics.tau_rpt_ens_bead = tau_rpt_ens_bead;
    EnsembleDynamics.tau_adj_ens_bead = tau_adj_ens_bead;
    
    EnsembleDynamics.tau_aSE_ens_bead = tau_aSE_ens_bead;
    EnsembleDynamics.tau_dSE_ens_bead = tau_dSE_ens_bead;
    EnsembleDynamics.tau_excSE_ens_bead = tau_excSE_ens_bead;
    EnsembleDynamics.tau_rptSE_ens_bead = tau_rptSE_ens_bead;
    EnsembleDynamics.tau_adjSE_ens_bead = tau_adjSE_ens_bead;
    
    EnsembleDynamics.kexc_ens_bead = kexc_ens_bead;
    EnsembleDynamics.kexcSE_ens_bead = kexcSE_ens_bead;
    EnsembleDynamics.krpt_ens_bead = krpt_ens_bead;
    EnsembleDynamics.krptSE_ens_bead = krptSE_ens_bead;


    EnsembleDynamics.tau_a_ens_meso = tau_a_ens_meso;
    EnsembleDynamics.tau_d_ens_meso = tau_d_ens_meso;
    EnsembleDynamics.tau_exc_ens_meso = tau_exc_ens_meso;
    EnsembleDynamics.tau_rpt_ens_meso = tau_rpt_ens_meso;
    EnsembleDynamics.tau_adj_ens_meso = tau_adj_ens_meso;
    
    EnsembleDynamics.tau_aSE_ens_meso = tau_aSE_ens_meso;
    EnsembleDynamics.tau_dSE_ens_meso = tau_dSE_ens_meso;
    EnsembleDynamics.tau_excSE_ens_meso = tau_excSE_ens_meso;
    EnsembleDynamics.tau_rptSE_ens_meso = tau_rptSE_ens_meso;
    EnsembleDynamics.tau_adjSE_ens_meso = tau_adjSE_ens_meso;
    
    EnsembleDynamics.ka_meso = ka_ens_meso;
    EnsembleDynamics.kaSE_meso = kaSE_ens_meso;
    EnsembleDynamics.ka1_meso = ka1_ens_meso;
    EnsembleDynamics.ka1SE_meso = ka1SE_ens_meso;
    EnsembleDynamics.kd_meso = kd_ens_meso;
    EnsembleDynamics.kdSE_meso = kdSE_ens_meso;
    EnsembleDynamics.Na_meso = Na_ens_meso;
    EnsembleDynamics.NaSE_meso = NaSE_ens_meso;
    EnsembleDynamics.Nd_meso = Nd_ens_meso;
    EnsembleDynamics.NdSE_meso = NdSE_ens_meso;
    EnsembleDynamics.fa_meso = fa_ens_meso;
    EnsembleDynamics.faSE_meso = faSE_ens_meso;
    EnsembleDynamics.fd_meso = fd_ens_meso;
    EnsembleDynamics.fdSE_meso = fdSE_ens_meso;

    EnsembleDynamics.kexc_ens_meso = kexc_ens_meso;
    EnsembleDynamics.kexcSE_ens_meso = kexcSE_ens_meso;
    EnsembleDynamics.krpt_ens_meso = krpt_ens_meso;
    EnsembleDynamics.krptSE_ens_meso = krptSE_ens_meso;

    if AppendScalingData
        EnsembleDynamics.tau_a_ens_thry = tau_a_ens_thry;
        EnsembleDynamics.tau_d_ens_thry = tau_d_ens_thry;
        EnsembleDynamics.tau_exc_ens_thry = tau_exc_ens_thry;
        EnsembleDynamics.tau_rpt_ens_thry = tau_rpt_ens_thry;
        EnsembleDynamics.tau_adj_ens_thry = tau_adj_ens_thry;

        EnsembleDynamics.tau_aSE_ens_thry = tau_aSE_ens_thry;
        EnsembleDynamics.tau_dSE_ens_thry = tau_dSE_ens_thry;
        EnsembleDynamics.tau_excSE_ens_thry = tau_excSE_ens_thry;
        EnsembleDynamics.tau_rptSE_ens_thry = tau_rptSE_ens_thry;
        EnsembleDynamics.tau_adjSE_ens_thry = tau_adjSE_ens_thry;
        
        EnsembleDynamics.ka_thry = ka_ens_thry;
        EnsembleDynamics.kaSE_thry = kaSE_ens_thry;
        EnsembleDynamics.ka1_thry = ka1_ens_thry;
        EnsembleDynamics.ka1SE_thry = ka1SE_ens_thry;
        EnsembleDynamics.kd_thry = kd_ens_thry;
        EnsembleDynamics.kdSE_thry = kdSE_ens_thry;
        EnsembleDynamics.Na_thry = Na_ens_thry;
        EnsembleDynamics.NaSE_thry = NaSE_ens_thry;
        EnsembleDynamics.Nd_thry = Nd_ens_thry;
        EnsembleDynamics.NdSE_thry = NdSE_ens_thry;
        EnsembleDynamics.fa_thry = fa_ens_thry;
        EnsembleDynamics.faSE_thry = faSE_ens_thry;
        EnsembleDynamics.fd_thry = fd_ens_thry;
        EnsembleDynamics.fdSE_thry = fdSE_ens_thry;

        EnsembleDynamics.kexc_ens_thry = kexc_ens_thry;
        EnsembleDynamics.kexcSE_ens_thry = kexcSE_ens_thry;
        EnsembleDynamics.krpt_ens_thry = krpt_ens_thry;
        EnsembleDynamics.krptSE_ens_thry = krptSE_ens_thry;
    end

    save(EnsembleDynamicsFileName,'-struct','EnsembleDynamics');
else
    EnsembleDynamics = load(EnsembleDynamicsFileName,'-mat');

    ka_ens_bead = EnsembleDynamics.ka_bead;
    kaSE_ens_bead = EnsembleDynamics.kaSE_bead;
    ka1_ens_bead = EnsembleDynamics.ka1_bead;
    ka1SE_ens_bead = EnsembleDynamics.ka1SE_bead;
    kd_ens_bead = EnsembleDynamics.kd_bead;
    kdSE_ens_bead =  EnsembleDynamics.kdSE_bead;
    Na_ens_bead = EnsembleDynamics.Na_bead;
    NaSE_ens_bead = EnsembleDynamics.NaSE_bead;
    Nd_ens_bead = EnsembleDynamics.Nd_bead;
    NdSE_ens_bead = EnsembleDynamics.NdSE_bead;
    fa_ens_bead = EnsembleDynamics.fa_bead;
    faSE_ens_bead = EnsembleDynamics.faSE_bead;
    fd_ens_bead = EnsembleDynamics.fd_bead;
    fdSE_ens_bead = EnsembleDynamics.fdSE_bead;

    tau_a_ens_bead = EnsembleDynamics.tau_a_ens_bead;
    tau_d_ens_bead = EnsembleDynamics.tau_d_ens_bead;
    tau_exc_ens_bead = EnsembleDynamics.tau_exc_ens_bead;
    tau_rpt_ens_bead = EnsembleDynamics.tau_rpt_ens_bead;
    tau_adj_ens_bead = EnsembleDynamics.tau_adj_ens_bead;

    tau_aSE_ens_bead = EnsembleDynamics.tau_aSE_ens_bead;
    tau_dSE_ens_bead = EnsembleDynamics.tau_dSE_ens_bead;
    tau_excSE_ens_bead = EnsembleDynamics.tau_excSE_ens_bead;
    tau_rptSE_ens_bead = EnsembleDynamics.tau_rptSE_ens_bead;
    tau_adjSE_ens_bead = EnsembleDynamics.tau_adjSE_ens_bead;
    
    kexc_ens_bead = EnsembleDynamics.kexc_ens_bead;
    kexcSE_ens_bead = EnsembleDynamics.kexc_ens_bead;
    krpt_ens_bead = EnsembleDynamics.krpt_ens_bead;
    krptSE_ens_bead = EnsembleDynamics.krpt_ens_bead;


    tau_a_ens_meso = EnsembleDynamics.tau_a_ens_meso;
    tau_d_ens_meso = EnsembleDynamics.tau_d_ens_meso;
    tau_exc_ens_meso = EnsembleDynamics.tau_exc_ens_meso;
    tau_rpt_ens_meso = EnsembleDynamics.tau_rpt_ens_meso;
    tau_adj_ens_meso = EnsembleDynamics.tau_rpt_ens_meso;

    tau_aSE_ens_meso = EnsembleDynamics.tau_aSE_ens_meso;
    tau_dSE_ens_meso = EnsembleDynamics.tau_dSE_ens_meso;
    tau_excSE_ens_meso = EnsembleDynamics.tau_excSE_ens_meso;
    tau_rptSE_ens_meso = EnsembleDynamics.tau_rptSE_ens_meso;
    tau_adjSE_ens_meso = EnsembleDynamics.tau_adjSE_ens_meso;
    
    ka_ens_meso = EnsembleDynamics.ka_meso;
    kaSE_ens_meso = EnsembleDynamics.kaSE_meso;
    ka1_ens_meso = EnsembleDynamics.ka1_meso;
    ka1SE_ens_meso = EnsembleDynamics.ka1SE_meso;
    kd_ens_meso = EnsembleDynamics.kd_meso;
    kdSE_ens_meso = EnsembleDynamics.kdSE_meso;
    Na_ens_meso = EnsembleDynamics.Na_meso;
    NaSE_ens_meso = EnsembleDynamics.NaSE_meso;
    Nd_ens_meso = EnsembleDynamics.Nd_meso;
    NdSE_ens_meso = EnsembleDynamics.NdSE_meso;
    fa_ens_meso = EnsembleDynamics.fa_meso;
    faSE_ens_meso = EnsembleDynamics.faSE_meso;
    fd_ens_meso = EnsembleDynamics.fd_meso;
    fdSE_ens_meso = EnsembleDynamics.fdSE_meso;

    kexc_ens_meso = EnsembleDynamics.kexc_ens_meso;
    kexcSE_ens_meso = EnsembleDynamics.kexc_ens_meso;
    krpt_ens_meso = EnsembleDynamics.krpt_ens_meso;
    krptSE_ens_meso = EnsembleDynamics.krpt_ens_meso;

    if AppendScalingData
        tau_a_ens_thry = EnsembleDynamics.tau_a_ens_thry;
        tau_d_ens_thry = EnsembleDynamics.tau_d_ens_thry;
        tau_exc_ens_thry = EnsembleDynamics.tau_exc_ens_thry;
        tau_rpt_ens_thry = EnsembleDynamics.tau_rpt_ens_thry;
        tau_adj_ens_thry = EnsembleDynamics.tau_adj_ens_thry;

        tau_aSE_ens_thry = EnsembleDynamics.tau_aSE_ens_thry;
        tau_dSE_ens_thry = EnsembleDynamics.tau_dSE_ens_thry;
        tau_excSE_ens_thry = EnsembleDynamics.tau_excSE_ens_thry;
        tau_rptSE_ens_thry = EnsembleDynamics.tau_rptSE_ens_thry;
        tau_adjSE_ens_thry = EnsembleDynamics.tau_adjSE_ens_thry;
        
        ka_ens_thry = EnsembleDynamics.ka_thry;
        kaSE_ens_thry = EnsembleDynamics.kaSE_thry;
        ka1_ens_thry = EnsembleDynamics.ka1_thry;
        ka1SE_ens_thry = EnsembleDynamics.ka1SE_thry;
        kd_ens_thry = EnsembleDynamics.kd_thry;
        kdSE_ens_thry = EnsembleDynamics.kdSE_thry;
        Na_ens_thry = EnsembleDynamics.Na_thry;
        NaSE_ens_thry = EnsembleDynamics.NaSE_thry;
        Nd_ens_thry = EnsembleDynamics.Nd_thry;
        NdSE_ens_thry = EnsembleDynamics.NdSE_thry;
        fa_ens_thry = EnsembleDynamics.fa_thry;
        faSE_ens_thry = EnsembleDynamics.faSE_thry;
        fd_ens_thry = EnsembleDynamics.fd_thry;
        fdSE_ens_thry = EnsembleDynamics.fdSE_thry;

        kexc_ens_thry = EnsembleDynamics.kexc_ens_thry;
        kexcSE_ens_thry = EnsembleDynamics.kexcSE_ens_thry;
        krpt_ens_thry = EnsembleDynamics.krpt_ens_thry;
        krptSE_ens_thry = EnsembleDynamics.krptSE_ens_thry;
    end
end

% Per Stukalin, et al.'s definition
tau_rnm_ens_bead = tau_exc_ens_bead+tau_adj_ens_bead;
tau_rnmSE_ens_bead = (tau_excSE_ens_bead.^2 + tau_adjSE_ens_bead.^2).^0.5;

tau_rnm_ens_meso = tau_exc_ens_meso+tau_adj_ens_meso;
tau_rnmSE_ens_meso = (tau_excSE_ens_meso.^2 + tau_adjSE_ens_meso.^2).^0.5;

if AppendScalingData
    tau_rnm_ens_thry = tau_exc_ens_thry+tau_adj_ens_thry;
    tau_rnmSE_ens_thry = (tau_excSE_ens_thry.^2 + tau_adjSE_ens_thry.^2).^0.5;
end

if ~isfile(AllDynamicsFileName) || Override
    AllDynamics.conc_all_bead = conc_all_bead;
    AllDynamics.time_bead = time_bead;
    AllDynamics.ka_bead = ka_all_bead;
    AllDynamics.kd_bead = kd_all_bead;
    AllDynamics.Na_bead = Na_all_bead;
    AllDynamics.Nd_bead = Nd_all_bead;
    AllDynamics.fa_bead = fa_all_bead;
    AllDynamics.fd_bead = fd_all_bead;
    AllDynamics.tau_a_bead = tau_a_all_bead;
    AllDynamics.tau_d_bead = tau_d_all_bead;
    AllDynamics.tau_exc_bead = tau_exc_all_bead;
    AllDynamics.tau_rpt_bead = tau_rpt_all_bead;
    AllDynamics.tau_adj_bead = tau_adj_all_bead;

    AllDynamics.conc_all_meso = conc_all_meso;
    AllDynamics.time_meso = time_meso;
    AllDynamics.ka_meso = ka_all_meso;
    AllDynamics.kd_meso = kd_all_meso;
    AllDynamics.Na_meso = Na_all_meso;
    AllDynamics.Nd_meso = Nd_all_meso;
    AllDynamics.fa_meso = fa_all_meso;
    AllDynamics.fd_meso = fd_all_meso;
    AllDynamics.tau_a_meso = tau_a_all_meso;
    AllDynamics.tau_d_meso = tau_d_all_meso;
    AllDynamics.tau_exc_meso = tau_exc_all_meso;
    AllDynamics.tau_rpt_meso = tau_rpt_all_meso;
    AllDynamics.tau_adj_meso = tau_adj_all_meso;

    if AppendScalingData
        AllDynamics.conc_all_thry = conc_all_thry;
        AllDynamics.time_thry = time_thry;
        AllDynamics.ka_thry = ka_all_thry;
        AllDynamics.kd_thry = kd_all_thry;
        AllDynamics.Na_thry = Na_all_thry;
        AllDynamics.Nd_thry = Nd_all_thry;
        AllDynamics.fa_thry = fa_all_thry;
        AllDynamics.fd_thry = fd_all_thry;
        AllDynamics.tau_a_thry = tau_a_all_thry;
        AllDynamics.tau_d_thry = tau_d_all_thry;
        AllDynamics.tau_exc_thry = tau_exc_all_thry;
        AllDynamics.tau_rpt_thry = tau_rpt_all_thry;
        AllDynamics.tau_adj_thry = tau_adj_all_thry;
    end

    save(AllDynamicsFileName,'-struct','AllDynamics');
else
    AllDynamics = load(AllDynamicsFileName,'-mat');

    conc_all_bead = AllDynamics.conc_all_bead;

    ka_all_bead = AllDynamics.ka_bead;
    kd_all_bead = AllDynamics.kd_bead;
    Na_all_bead = AllDynamics.Na_bead;
    Nd_all_bead = AllDynamics.Nd_bead;
    fa_all_bead = AllDynamics.fa_bead;
    fd_all_bead = AllDynamics.fd_bead;
    tau_a_all_bead = AllDynamics.tau_a_bead;
    tau_d_all_bead = AllDynamics.tau_d_bead;
    tau_exc_all_bead = AllDynamics.tau_exc_bead;
    tau_rpt_all_bead = AllDynamics.tau_rpt_bead;
    tau_adj_all_bead = AllDynamics.tau_adj_bead;

    conc_all_meso = AllDynamics.conc_all_meso;

    ka_all_meso = AllDynamics.ka_meso;
    kd_all_meso = AllDynamics.kd_meso;
    Na_all_meso = AllDynamics.Na_meso;
    Nd_all_meso = AllDynamics.Nd_meso;
    fa_all_meso = AllDynamics.fa_meso;
    fd_all_meso = AllDynamics.fd_meso;
    tau_a_all_meso = AllDynamics.tau_a_meso;
    tau_d_all_meso = AllDynamics.tau_d_meso;
    tau_exc_all_meso = AllDynamics.tau_exc_meso;
    tau_rpt_all_meso = AllDynamics.tau_rpt_meso;
    tau_adj_all_meso = AllDynamics.tau_adj_meso;

    if AppendScalingData
        conc_all_thry = AllDynamics.conc_all_thry;

        ka_all_thry = AllDynamics.ka_thry;
        kd_all_thry = AllDynamics.kd_thry;
        Na_all_thry = AllDynamics.Na_thry;
        Nd_all_thry = AllDynamics.Nd_thry;
        fa_all_thry = AllDynamics.fa_thry;
        fd_all_thry = AllDynamics.fd_thry;
        tau_a_all_thry = AllDynamics.tau_a_thry;
        tau_d_all_thry = AllDynamics.tau_d_thry;
        tau_exc_all_thry = AllDynamics.tau_exc_thry;
        tau_rpt_all_thry = AllDynamics.tau_rpt_thry;
        tau_adj_all_thry = AllDynamics.tau_adj_thry;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotKineticswrtTime(time_bs,time_ms,k_bs,k_ms,color,tau0,...
    fig_folder,phis,N_Kuhn,sp,save_tag,mm_pcnt)

k_bs(isnan(k_bs)) = 0;
k_ms(isnan(k_ms)) = 0;
mm_divisor = 1/(mm_pcnt/100);
mm_win = round(length(time_bs)/mm_divisor);
rng = (1:500:length(time_bs)-1);
k_bs_smth = movmean(k_bs,mm_win);
s = scatter(time_bs(rng)/tau0,k_bs_smth(rng));
s.MarkerFaceColor = 'none';
s.MarkerEdgeColor = color;
s.SizeData = 10;

mm_win = round(length(time_ms)/mm_divisor);
rng = (1:length(time_ms)-1);
k_ms_smth = movmean(k_ms,mm_win);
p = plot(time_ms(rng)/tau0,k_ms_smth(rng));
p.Color = color;
p.LineWidth = 1.5;

if sp==length(phis)
    ylab = ['$',save_tag,'\tau_0$'];
    set(gca,'FontSize',20/1.5)
    xlabel('$t/\tau_0$','FontSize',20,'Interpreter','latex')
    ylabel(ylab,'FontSize',20,'Interpreter','latex')
    set(gcf,'color','w')
    pbaspect([1 1 1])
    set(gcf,'Position',[100 400 275 275])
    saveas(gcf,[fig_folder,'/ka.N',num2str(N_Kuhn),'.png'])
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DefineAssembledFileNames(ka_in)

global OutputDir EnsembleDynamicsFileName AllDynamicsFileName...
    FittedDynamicsFileName R2DynamicsFileName

EnsembleDynamicsFileName = [OutputDir,'EnsembleDynamics.ka-',...
    num2str(ka_in,'%.2e'),'.m'];
AllDynamicsFileName = [OutputDir,'AllDynamicswrtTime.ka-',...
    num2str(ka_in,'%.2e'),'.m'];
R2DynamicsFileName = [OutputDir,'R2Dynamics.ka-',...
    num2str(ka_in,'%.2e'),'.m'];
FittedDynamicsFileName = [OutputDir,'FittedDynamics.ka-',...
    num2str(ka_in,'%.2e'),'.m'];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitializeAllDynamics(time_bead,time_meso,NoSep)

global ka_all_meso kd_all_meso...
    Na_all_meso Nd_all_meso...
    fa_all_meso fd_all_meso...
    ka_all_bead kd_all_bead...
    Na_all_bead Nd_all_bead...
    fa_all_bead fd_all_bead...
    tau_a_all_bead tau_d_all_bead...
    tau_exc_all_bead tau_rpt_all_bead tau_adj_all_bead...
    tau_a_all_meso tau_d_all_meso...
    tau_exc_all_meso tau_rpt_all_meso tau_adj_all_meso...

ka_all_bead = zeros(3*size(time_bead,1),NoSep);
kd_all_bead = zeros(3*size(time_bead,1),NoSep);
Na_all_bead = zeros(3*size(time_bead,1),NoSep);
Nd_all_bead = zeros(3*size(time_bead,1),NoSep);
fa_all_bead = zeros(3*size(time_bead,1),NoSep);
fd_all_bead = zeros(3*size(time_bead,1),NoSep);

ka_all_meso = zeros(3*size(time_meso,1),NoSep);
kd_all_meso = zeros(3*size(time_meso,1),NoSep);
Na_all_meso = zeros(3*size(time_meso,1),NoSep);
Nd_all_meso = zeros(3*size(time_meso,1),NoSep);
fa_all_meso = zeros(3*size(time_meso,1),NoSep);
fd_all_meso = zeros(3*size(time_meso,1),NoSep);

tau_a_all_bead = zeros(5e3,NoSep);
tau_d_all_bead = zeros(5e3,NoSep);
tau_exc_all_bead = zeros(5e3,NoSep);
tau_rpt_all_bead = zeros(5e3,NoSep);
tau_adj_all_bead = zeros(5e3,NoSep);

tau_a_all_meso = zeros(5e3,NoSep);
tau_d_all_meso = zeros(5e3,NoSep);
tau_exc_all_meso = zeros(5e3,NoSep);
tau_rpt_all_meso = zeros(5e3,NoSep);
tau_adj_all_meso = zeros(5e3,NoSep);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitializeEnsembleDynamics(N_Kuhns,NoSep)

global ka_ens_meso kd_ens_meso kaSE_ens_meso kdSE_ens_meso...
    Na_ens_meso Nd_ens_meso NaSE_ens_meso NdSE_ens_meso...
    fa_ens_meso fd_ens_meso faSE_ens_meso fdSE_ens_meso...
    ka_ens_bead kd_ens_bead kaSE_ens_bead kdSE_ens_bead...
    Na_ens_bead Nd_ens_bead NaSE_ens_bead NdSE_ens_bead...
    fa_ens_bead fd_ens_bead faSE_ens_bead fdSE_ens_bead...
    tau_a1_ens_bead tau_a1SE_ens_bead...
    tau_a_ens_bead tau_aSE_ens_bead...
    tau_d_ens_bead tau_dSE_ens_bead...
    tau_a1_ens_meso tau_a1SE_ens_meso...
    tau_a_ens_meso tau_aSE_ens_meso...
    tau_d_ens_meso tau_dSE_ens_meso...
    conc_all_meso conc_all_bead...
    tau_exc_ens_meso tau_excSE_ens_meso...
    tau_rpt_ens_meso tau_rptSE_ens_meso...
    tau_adj_ens_meso tau_adjSE_ens_meso...
    kexc_ens_meso kexcSE_ens_meso...
    krpt_ens_meso krptSE_ens_meso

conc_all_meso = zeros(NoSep,length(N_Kuhns));
conc_all_bead = zeros(NoSep,length(N_Kuhns));

ka_ens_meso = zeros(NoSep,length(N_Kuhns));
kd_ens_meso = zeros(NoSep,length(N_Kuhns));
kaSE_ens_meso = zeros(NoSep,length(N_Kuhns));
kdSE_ens_meso = zeros(NoSep,length(N_Kuhns));
Na_ens_meso = zeros(NoSep,length(N_Kuhns));
Nd_ens_meso = zeros(NoSep,length(N_Kuhns));
NaSE_ens_meso = zeros(NoSep,length(N_Kuhns));
NdSE_ens_meso = zeros(NoSep,length(N_Kuhns));
fa_ens_meso = zeros(NoSep,length(N_Kuhns));
fd_ens_meso = zeros(NoSep,length(N_Kuhns));
faSE_ens_meso = zeros(NoSep,length(N_Kuhns));
fdSE_ens_meso = zeros(NoSep,length(N_Kuhns));

ka_ens_bead = zeros(NoSep,length(N_Kuhns));
kd_ens_bead = zeros(NoSep,length(N_Kuhns));
kaSE_ens_bead = zeros(NoSep,length(N_Kuhns));
kdSE_ens_bead = zeros(NoSep,length(N_Kuhns));
Na_ens_bead = zeros(NoSep,length(N_Kuhns));
Nd_ens_bead = zeros(NoSep,length(N_Kuhns));
NaSE_ens_bead = zeros(NoSep,length(N_Kuhns));
NdSE_ens_bead = zeros(NoSep,length(N_Kuhns));
fa_ens_bead = zeros(NoSep,length(N_Kuhns));
fd_ens_bead = zeros(NoSep,length(N_Kuhns));
faSE_ens_bead = zeros(NoSep,length(N_Kuhns));
fdSE_ens_bead = zeros(NoSep,length(N_Kuhns));

tau_a1_ens_bead = zeros(NoSep,length(N_Kuhns));
tau_a1SE_ens_bead = zeros(NoSep,length(N_Kuhns));
tau_a_ens_bead = zeros(NoSep,length(N_Kuhns));
tau_aSE_ens_bead = zeros(NoSep,length(N_Kuhns));
tau_d_ens_bead = zeros(NoSep,length(N_Kuhns));
tau_dSE_ens_bead = zeros(NoSep,length(N_Kuhns));

tau_a1_ens_meso = zeros(NoSep,length(N_Kuhns));
tau_a1SE_ens_meso = zeros(NoSep,length(N_Kuhns));
tau_a_ens_meso = zeros(NoSep,length(N_Kuhns));
tau_aSE_ens_meso = zeros(NoSep,length(N_Kuhns));
tau_d_ens_meso = zeros(NoSep,length(N_Kuhns));
tau_dSE_ens_meso = zeros(NoSep,length(N_Kuhns));

tau_exc_ens_meso = zeros(NoSep,length(N_Kuhns));
tau_excSE_ens_meso = zeros(NoSep,length(N_Kuhns));
tau_rpt_ens_meso = zeros(NoSep,length(N_Kuhns));
tau_rptSE_ens_meso = zeros(NoSep,length(N_Kuhns));
tau_adj_ens_meso = zeros(NoSep,length(N_Kuhns));
tau_adjSE_ens_meso = zeros(NoSep,length(N_Kuhns));
kexc_ens_meso = zeros(NoSep,length(N_Kuhns));
kexcSE_ens_meso = zeros(NoSep,length(N_Kuhns));
krpt_ens_meso = zeros(NoSep,length(N_Kuhns));
krptSE_ens_meso = zeros(NoSep,length(N_Kuhns));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time,conc,ka,ka1,kd,...
    k_exc,k_rpt,...
    Na,Nd,fa,fd,...
    tau_a,tau_d,tau_exchange,tau_repeat,tau_renorm,...
    mean_ka,se_ka,mean_ka1,se_ka1,mean_kd,se_kd,...
    mean_Na,se_Na,mean_Nd,se_Nd,...
    mean_fa,se_fa,mean_fd,se_fd,...
    mean_tau_a,se_tau_a,...
    mean_tau_d,se_tau_d,...
    mean_tau_exchange,se_tau_exchange,...
    mean_tau_repeat,se_tau_repeat,...
    mean_tau_renorm,se_tau_renorm,...
    mean_k_exc,se_k_exc,mean_k_rpt,se_k_rpt] =...
    CalculateEnsembleDynamics(BSOM,Sample,Np,Ns,N,phi,eaStar,edStar,...
    dt,N_Kuhn,b,D,kbT)

global BeadSpringOrMeso TimeStretchDataFileName BondKineticsDataFileName...
    NormalizeRates tau0 LengthConversion

% Switch for comparing kinetic rates measured over entire time domain or
% steady sate
SampleAllOrSS = 0;

% Define model type and inputs, and set directories/file names
BeadSpringOrMeso = BSOM;
InputScript(Sample,Np,Ns,N,N_Kuhn,b,phi,eaStar,edStar,D,dt,kbT);
SetDirAndFileNames;
DefineCompiledFileNames;
tau0 = (b*LengthConversion)^2/D;

% Import the time and stretch data
TimeAndStretch = load(TimeStretchDataFileName,'-mat');
time= TimeAndStretch.time;
Delta_t = time(2)-time(1);

% Import the bond kinetics data
BondKinetics = load(BondKineticsDataFileName,'-mat');
ka = BondKinetics.ka;       % raw attachment rates
ka1 = BondKinetics.ka1;     % raw first-time attachment rates
kd = BondKinetics.kd;       % raw detachment rates
Na = BondKinetics.Na;       % raw number of attached chains
Nd = BondKinetics.Nd;       % raw number of detached chains
fa = Na./(Na+Nd);           % computed attached bond fraction
fd = Nd./(Na+Nd);           % computed detached bond fraction
conc = BondKinetics.Concentration;  % chain concentration
tau_a = BondKinetics.AttachedLifetimes; % ensemble of attached bond lifetimes
tau_d = BondKinetics.DetachedLifetimes; % ensemble of detached bond lifetimes
tau_exchange = BondKinetics.OpenExchangeLifetimes;  % ensemble of detached bond lifetimes prior to exchange
tau_repeat = BondKinetics.OpenRepeatLifetimes;      % ensemble of detached bond lifetiems prior to repeat attachment
tau_renorm = BondKinetics.RenormalizedLifetimes;    % ensemble of renormalized bond lifetimes as defined by Stukalin, et al. (2013)
exchange_events = BondKinetics.ExchangeEvents;  % number of exchange events at each timestep
repeat_events = BondKinetics.RepeatEvents;      % number of repeat events at each timestep

% ID the steady state regime as when 
n_smooth = round(length(time)/5);   % smooth data over a 20% window to mitigate 
                                    % noise and smooth gradients can be taken
                                    % for rate change of attached bond
                                    % fractions
smooth_Na = smooth(Na,n_smooth);
smooth_Nd = smooth(Nd,n_smooth);
smooth_fa = smooth(fa,n_smooth);
n_pts = 100;
rng = round(linspace(1,length(time),n_pts));
dfa = gradient(smooth_fa(rng))./smooth_fa(rng)*100;

temp_time = time(rng);
SSIndx = find(dfa<0.1,1,'first'); % ID index at which change in attached fraction with respect to 1/100th of total time is less than 1%
ss_time = temp_time(SSIndx);
SSIndx = find(ss_time<=time,1,'first');

SSIndx_pcnt = round(SSIndx)/length(time);

showssfigs = 0; % set 1 to make sure that index is reasonably capture "steady state" conditions
if showssfigs==1
    figure(1); clf; hold on

    subplot(2,1,1); hold on
    plot(time,Na,'c-'); plot(time,Nd,'r-');
    plot(time,smooth_Na,'c--'); plot(time,smooth_Nd,'r--');
    plot([time(SSIndx),time(SSIndx)],[0 343],'k--')
    xlabel('time')
    ylabel('No. attchd/dtchd bonds')

    subplot(2,1,2); hold on

    plot([time(SSIndx),time(SSIndx)],[-100 100],'k--')
    plot(temp_time,dNadt,'c-.'); plot(temp_time,dNddt,'r-.');
    ylim([-100 100])
    xlabel('time')
    ylabel('dN/dt')

    close
end

% If normalizing the output rates
if NormalizeRates==1
    ka = ka*tau0;
    ka1 = ka1*tau0;
    kd = kd*tau0;
end

% Compute the mean and standard error of the mean for attached and detached
% fractions only at steady state
avg_rng = SSIndx:length(fa);
mean_fa = mean(fa(avg_rng,:),1);
se_fa = ComputeSEM(fa, 100);
mean_fd = mean(fd(avg_rng,:),1);
se_fd = ComputeSEM(fd, 100);


% Compute means and standard error of the mean for virgin attachment rate,
% total attachment rate, detachment rate, number of attachment events,
% number of detachment events, 
if SampleAllOrSS==0 % if wanting to sample all time for kinetic rates
    avg_rng = (1:avg_rng(end)); % reset averaging range to all samples
end
if NormalizeRates==1
    mean_ka1 = mean(ka1,1);         %
    se_ka1 = ComputeSEM(ka1,100);
    mean_ka = mean(ka(avg_rng),1); 
    se_ka = ComputeSEM(ka, 100); 
    mean_kd = nanmean(kd(avg_rng),1);     
    se_kd = ComputeSEM(kd, 100); 
else
    mean_ka1 = mean(ka1,1)/1e3;
    se_ka1 = ComputeSEM(ka1(avg_rng)/1e3,100);
    mean_ka = mean(ka(avg_rng),1)/1e3; 
    se_ka = ComputeSEM(ka(avg_rng)/1e3,100);
    mean_kd = mean(kd(avg_rng),1)/1e3; 
    se_kd = ComputeSEM(kd(avg_rng)/1e3,100);
end
mean_Na = mean(Na(avg_rng),1); 
se_Na = ComputeSEM(Na, 100);
mean_Nd = mean(Nd(avg_rng),1); 
se_Nd = ComputeSEM(Nd, 100);

% Compute averaged lifetimes of events
if round(length(tau_a)*SSIndx_pcnt)~=0
    tau_a = tau_a(ceil(length(tau_a)*SSIndx_pcnt):end);
    tau_d = tau_d(ceil(length(tau_d)*SSIndx_pcnt):end);
    tau_exchange = tau_exchange(ceil(length(tau_exchange)*SSIndx_pcnt):end);
    tau_repeat = tau_repeat(ceil(length(tau_repeat)*SSIndx_pcnt):end);
    tau_renorm = tau_renorm(ceil(length(tau_renorm)*SSIndx_pcnt):end);
end

if NormalizeRates==1
    tau_a = tau_a/tau0;
    tau_d = tau_d/tau0;
    tau_exchange = tau_exchange/tau0;
    tau_repeat = tau_repeat/tau0;
    tau_renorm = tau_renorm/tau0;
end

% Compute exchange and repeat rates using same methodology as ka, ka1, and
% kd
k_exc = exchange_events./Nd/Delta_t;	% k_exc(t) = N_exc(t)/Nd(t)/dt/2
                        % where N_exc(t) is the number of exchange events 
                        % at time, t, Nd(t) is the number of detached bonds
                        % at time, t, dt is the time duration between
                        % observations, and it is multiplied by two because
                        % attachment rates should be normalized by no. of 
                        % open pairs, not open chains. 
k_rpt = repeat_events./Nd/Delta_t;         
if NormalizeRates==1
    k_exc = k_exc*tau0;
    k_rpt = k_rpt*tau0;
end
mean_k_exc = nanmean(k_exc(avg_rng,:));
se_k_exc = ComputeSEM(k_exc, 100);
mean_k_rpt = nanmean(k_rpt(avg_rng,:));
se_k_rpt = ComputeSEM(k_rpt, 100);


% total_repeats = sum(repeat_events);
% total_chains = Nd(1);
% repeats_per_chain = total_repeats/total_chains;
% 
% total_repeats = length(tau_a);


if mean_k_exc>mean_k_rpt
   disp(['No. repeats is ',num2str(sum(repeat_events))])
   disp(['No. exchanges is ',num2str(sum(exchange_events))])
   disp(['Avg. repeat rate ',num2str(sum(mean_k_rpt))])
   disp(['Avg. exchange rate ',num2str(sum(mean_k_exc))])
end

% Compute mean and standard error of mean of lifetimes if at least
% 'mi_sample_size' observations were made
min_sample_size = 10; 
[mean_tau_a,se_tau_a] = DetermineMeanAndSEIfNGeq(tau_a,min_sample_size);
[mean_tau_d,se_tau_d] = DetermineMeanAndSEIfNGeq(tau_d,min_sample_size);
[mean_tau_exchange,se_tau_exchange] = DetermineMeanAndSEIfNGeq(tau_exchange,min_sample_size);
[mean_tau_repeat,se_tau_repeat] = DetermineMeanAndSEIfNGeq(tau_repeat,min_sample_size);
[mean_tau_renorm,se_tau_renorm] = DetermineMeanAndSEIfNGeq(tau_renorm,min_sample_size);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SEM = ComputeSEM(data, num_bins)
    % Check if input data is a column vector
    if ~iscolumn(data)
        error('Input data must be a column vector');
    end
    
    % Check if number of bins is valid
    if num_bins < 1
        error('Number of bins must be at least 1');
    end
    
    % Compute the number of elements in each bin
    bin_size = ceil(length(data) / num_bins);
    
    % Randomly shuffle the data
    shuffled_data = data(randperm(length(data)));
    
    ct = 0;
    means = [];
    rng = 0;
    while ~isempty(shuffled_data)
        ct = ct+1;
        if size(shuffled_data,1)>=bin_size
            rng = 1:bin_size;
        else
            rng = 1:size(shuffled_data,1);
        end
        dat_temp = shuffled_data(rng);
        mean_temp = nanmean(dat_temp);
        means = cat(1,means,mean_temp);
        shuffled_data(rng) = [];
    end
    
%     % Bin the data into approximately evenly sized arrays
%     bins = mat2cell(shuffled_data, bin_size * ones(1, num_bins-1), 1);
%     bins{num_bins} = shuffled_data((num_bins-1)*bin_size+1:end);
%     
%     % Compute averages of each bin
%     bin_means = cellfun(@mean, bins);
    
    % Compute standard error of the averages
    SEM = nanstd(means) / sqrt(ct);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mean_tau,se_tau] = DetermineMeanAndSEIfNGeq(tau,min_sample_size)

if length(tau)>=min_sample_size
    mean_tau = mean(tau);
%     se_tau = NaN;
%     se_tau = std(tau);
%     se_tau = std(tau)/sqrt(length(tau));
    se_tau = ComputeSEM(tau, 100);
else
    mean_tau = NaN;
    se_tau = NaN;
end

end
