function PostProcessData(Package,OverridePostProccess,ToggleDynamics,...
    LC,DC,BSOM,FPE)

% Computes and plots stress-strain for every iteration of simulation

global LineWidth TurnOnDynamics FontSize NoFitAttempts...
    LengthConversion DamperConversion DataSize BeadSpringOrMeso...
    OutputFolder FitPowerOrExponent PlotScalingTheory NormalizeRates

OutputFolder = 'Output Plots';
if ~isfolder(OutputFolder)
    mkdir(OutputFolder)
end

NoFitAttempts = 75;
LineWidth = 1.5;
FontSize = 20;
DataSize = 30;
TurnOnDynamics = ToggleDynamics;
LengthConversion = LC;
DamperConversion = DC;
BeadSpringOrMeso = BSOM;
FitPowerOrExponent = FPE;
PlotScalingTheory = 0;
NormalizeRates = 1;

% Samples = unique(Package(:,1));
Np = unique(Package(:,2));        %Number of molecules
Ns = unique(Package(:,3));        %Number of stickeres per tether site
T = unique(Package(:,4));         %Temperature
Nb = unique(Package(:,5));        %Contour length of chains
ka_in = unique(Package(:,6));        %Activation energy of association
kd_in = unique(Package(:,7));        %Activation energy of dissociation
f0 = unique(Package(:,8));       %Force sensitivity to dissociation
dts = unique(Package(:,9));       %timestep size
DiffCoeffs = unique(Package(:,11));
N_Kuhns = unique(Package(:,12));
b = unique(Package(:,13));
N = Np*Ns;
dtFacts = unique(Package(:,15));

%% Sweepign parameters are N_Kuhn and b
damps = unique(Package(:,10));    %damping coefficient in units [mass/time]

%% Assemble the Ensemble Data
AssembleEnsembleDynamicsData(OverridePostProccess,Package,damps,N_Kuhns,Np,Ns,T,Nb,...
    ka_in,kd_in,f0,dts,b,N);

%% Plot Outputs
% Dynamics wrt. time
PlotDynamicsWrtTime(OverridePostProccess,Package,N_Kuhns,b,dts);

% Ensemble average dynamics
PlotEnsembleAvgDynamics(OverridePostProccess,Package,N_Kuhns,b,...
    dts,ka_in,kd_in,DiffCoeffs,dtFacts);

close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotEnsembleAvgDynamics(Override,Package,N_Kuhns,b,dts,...
    ka_in,kd_in,D,dtFacts)

global EnsembleDynamicsFileName FontSize OutputFolder...
    LengthConversion FittedDynamicsFileName NormalizeRates

% For plotting outputs wrt time
figure(1); clf; hold on %ka
figure(2); clf; hold on %kd
figure(3); clf; hold on %Na
figure(4); clf; hold on %Nd

A = zeros(size(N_Kuhns));
R2_meso_vs_bead = zeros(size(N_Kuhns));

ColorRange = (linspace(0,1,length(dts)))';
ColorFactn = [1 165/255 0];
Colorsn = [ColorRange ColorRange ColorRange].*ColorFactn;

NormLength = b*LengthConversion;
tau0 = (NormLength^2)/D;

dts = flipud(dts);
plotct = 0;
for dti=1:length(dts)
    dt = dts(dti);
    dtFact = dtFacts(dti);
    AddOn = ['Ensemble.dt-',num2str(dt,'%.2e'),'.'];
    plotct = plotct+1;

    DefineAssembledFileNames(dt);

    if ~isfile(FittedDynamicsFileName) || Override

        %Import the data
        EnsembleDynamics = load(EnsembleDynamicsFileName,'-mat');

        ka_bead = EnsembleDynamics.ka_bead*1e3;
        kaSE_bead = EnsembleDynamics.kaSE_bead*1e3;
        kd_bead = EnsembleDynamics.kd_bead*1e3;
        kdSE_bead = EnsembleDynamics.kdSE_bead*1e3;

        ka_meso = EnsembleDynamics.ka_meso*1e3;
        kaSE_meso = EnsembleDynamics.kaSE_meso*1e3;
        kd_meso = EnsembleDynamics.kd_meso*1e3;
        kdSE_meso = EnsembleDynamics.kdSE_meso*1e3;

        fa_bead = EnsembleDynamics.fa_bead;
        faSE_bead = EnsembleDynamics.faSE_bead;
        fd_bead = EnsembleDynamics.fd_bead;
        fdSE_bead = EnsembleDynamics.fdSE_bead;

        fa_meso = EnsembleDynamics.fa_meso;
        faSE_meso = EnsembleDynamics.faSE_meso;
        fd_meso = EnsembleDynamics.fd_meso;
        fdSE_meso = EnsembleDynamics.fdSE_meso;

        % Loop over N_Kuhn
        for nk = 1:length(N_Kuhns)
            N_Kuhn = N_Kuhns(nk);
            PackageTemp = Package(Package(:,12)==N_Kuhn,:);
            Separations = unique(PackageTemp(:,14));

            % Isolte data for N_Kuhn
            if NormalizeRates==1
                ka_meso_temp = ka_meso(:,nk)*tau0;
                ka_bead_temp = ka_bead(:,nk)*tau0;
                kaSE_meso_temp = kaSE_meso(:,nk)*tau0;
                kaSE_bead_temp = kaSE_bead(:,nk)*tau0;

                kd_meso_temp = kd_meso(:,nk)*tau0;
                kd_bead_temp = kd_bead(:,nk)*tau0;
                kdSE_meso_temp = kdSE_meso(:,nk)*tau0;
                kdSE_bead_temp = kdSE_bead(:,nk)*tau0;

                kd_plot = kd_in*ones(size(Separations))*tau0;
                yhi = 6e-3;
            else
                ka_meso_temp = ka_meso(:,nk);
                ka_bead_temp = ka_bead(:,nk);
                kaSE_meso_temp = kaSE_meso(:,nk);
                kaSE_bead_temp = kaSE_bead(:,nk);

                kd_meso_temp = kd_meso(:,nk);
                kd_bead_temp = kd_bead(:,nk);
                kdSE_meso_temp = kdSE_meso(:,nk);
                kdSE_bead_temp = kdSE_bead(:,nk);

                kd_plot = kd_in*ones(size(Separations))/1e3;
                yhi = 1.25e3;
            end

            fa_bead_temp = fa_bead(:,nk);
            faSE_bead_temp = faSE_bead(:,nk);
            fa_meso_temp = fa_meso(:,nk);
            faSE_meso_temp = faSE_meso(:,nk);

            fd_bead_temp = fd_bead(:,nk);
            fdSE_bead_temp = fdSE_bead(:,nk);
            fd_meso_temp = fd_meso(:,nk);
            fdSE_meso_temp = fdSE_meso(:,nk);

            % Plot attachment rate wrt. sep
            A_nom = 0.9;
            x = Separations;
            y = ka_bead_temp;

            % Regression analysis between models
            RSS = sum((ka_bead_temp-ka_meso_temp).^2);
            TSS = sum((ka_meso_temp-mean(ka_meso_temp)).^2);
            R2_meso_vs_bead(dti) = 1-RSS/TSS;

            NormLength = b*LengthConversion;
            if NormalizeRates==1
                eaOverkbT = -log(ka_in*tau0*NormLength^2/D);
            else
                eaOverkbT = -log(ka_in*NormLength^2/D);
            end
            sigma_a_theor = (3/4)^(-1/2)*sqrt(N_Kuhn)*b;
            A(dti) = FitThePrefactor(x,y,A_nom,...
                D,NormLength,N_Kuhn,eaOverkbT,sigma_a_theor);

            LegendEntries{plotct} = ['$\Delta t^*$ = $\tau_0/$',num2str(dtFact)];

            figure(1)
            xplot = (linspace(Separations(1),Separations(end),100))';
            ka0_theor = A(dti)*3/4/pi*D/(NormLength^2)*N_Kuhn^(-3/2)*exp(-eaOverkbT);
            ka_theor = ka0_theor*exp(-(xplot.^2)/(sigma_a_theor^2));
            
            Ptheor(plotct) = plot(xplot/(N_Kuhn*b),ka_theor);
            Ptheor(plotct).LineStyle = '-';
            Ptheor(plotct).Color = Colorsn(plotct,:);
            Ptheor(plotct).LineWidth = 1.5;

            e =  errorbar(Separations/(N_Kuhn*b),ka_bead_temp,kaSE_bead_temp);
            e.LineStyle = 'none';
            e.Color = 'k';
            e.Marker = 'o';
            e.MarkerFaceColor = Colorsn(plotct,:);

            e =  errorbar(Separations/(N_Kuhn*b),ka_meso_temp,kaSE_meso_temp);
            e.LineStyle = 'none';
            e.Color  = 'k';
            e.Marker = '^';
            e.MarkerFaceColor = Colorsn(plotct,:);
            

            % Plot detachment rate wrt. sep
            figure(2)
            plot(Separations/(N_Kuhn*b),kd_plot,'k--')

            e =  errorbar(Separations/(N_Kuhn*b),kd_bead_temp,kdSE_bead_temp);
            e.LineStyle = 'none';
            e.Color = 'k';
            e.Marker = 'o';
            e.MarkerFaceColor = Colorsn(plotct,:);

            e =  errorbar(Separations/(N_Kuhn*b),kd_meso_temp,kdSE_meso_temp');
            e.LineStyle = 'none';
            e.Color = 'k';
            e.Marker = '^';
            e.MarkerFaceColor = Colorsn(plotct,:);
            % Plot just for 1 Nk

            if nk==1
                figure(10); clf; hold on

                e =  errorbar(Separations/(N_Kuhn*b),ka_bead_temp,kaSE_bead_temp);
                e.LineStyle = 'none';
                e.Color = 'k';
                e.Marker = 'o';
                e.MarkerFaceColor = 'k';

                e =  errorbar(Separations/(N_Kuhn*b),ka_meso_temp,kaSE_meso_temp);
                e.LineStyle = 'none';
                e.Color  = 'k';
                e.Marker = '^';
                e.MarkerFaceColor = 'none';

                set(gca,'FontSize',FontSize/1.5)
                set(gcf,'Color','w')
                pbaspect([1 1 1])

                l = legend('Bead-spring','Mesoscale');
                l.FontSize = FontSize/2;
                l.Interpreter = 'latex';

                xlabel('$d^*$ ($Nb$)','FontSize',FontSize,'Interpreter','latex')
                if NormalizeRates==1
                    ylabel('$\bar{k}_a$ ($\tau_0^{-1}$)','FontSize',FontSize,'Interpreter','latex')
                else
                    ylabel('$\bar{k}_a$ (kHz)','FontSize',FontSize,'Interpreter','latex')
                end

                ylim([0 yhi])

                FileTag = ['dt.',num2str(dt,'%.2e'),'.ka vs d'];
                saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
                saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])

                figure(20); clf; hold on

                e =  errorbar(Separations/(N_Kuhn*b),kd_bead_temp,kdSE_bead_temp);
                e.LineStyle = 'none';
                e.Color = 'k';
                e.Marker = 'o';
                e.MarkerFaceColor = 'k';

                e =  errorbar(Separations/(N_Kuhn*b),kd_meso_temp,kdSE_meso_temp);
                e.LineStyle = 'none';
                e.Color  = 'k';
                e.Marker = '^';
                e.MarkerFaceColor = 'none';

                plot(Separations/(N_Kuhn*b),kd_plot,'k--')

                set(gca,'FontSize',FontSize/1.5)
                set(gcf,'Color','w')
                pbaspect([1 1 1])

                l = legend('Bead-spring','Mesoscale');
                l.FontSize = FontSize/2;
                l.Interpreter = 'latex';

                ylim([0 yhi])

                xlabel('$d^*$ ($Nb$)','FontSize',FontSize,'Interpreter','latex')
                if NormalizeRates==1
                    ylabel('$\bar{k}_d$ ($\tau_0^{-1}$)','FontSize',FontSize,'Interpreter','latex')
                else
                    ylabel('$\bar{k}_d$ (kHz)','FontSize',FontSize,'Interpreter','latex')
                end

                FileTag = ['dt.',num2str(dt,'%.2e'),'.kd vs d'];
                saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
                saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])
            end

            figure(3)

            fa_theor = ka_theor./(ka_theor+unique(kd_plot));
            p = plot(xplot/(N_Kuhn*b),fa_theor);
            p.Color = Colorsn(plotct,:);
            p.LineWidth = 1.5;
                        
            e =  errorbar(Separations/(N_Kuhn*b),fa_bead_temp,faSE_bead_temp);
            e.LineStyle = 'none';
            e.Color = 'k';
            e.Marker = 'o';
            e.MarkerFaceColor = Colorsn(plotct,:);

            e =  errorbar(Separations/(N_Kuhn*b),fa_meso_temp,faSE_meso_temp);
            e.LineStyle = 'none';
            e.Color  = 'k';
            e.Marker = '^';
            e.MarkerFaceColor = Colorsn(plotct,:);

            figure(4)

            fd_theor = 1-ka_theor./(ka_theor+unique(kd_plot));
            p = plot(xplot/(N_Kuhn*b),fd_theor);
            p.Color = Colorsn(plotct,:);
            p.LineWidth = 1.5;

            e =  errorbar(Separations/(N_Kuhn*b),fd_bead_temp,fdSE_bead_temp);
            e.LineStyle = 'none';
            e.Color = 'k';
            e.Marker = 'o';
            e.MarkerFaceColor = Colorsn(plotct,:);

            e =  errorbar(Separations/(N_Kuhn*b),fd_meso_temp,fdSE_meso_temp);
            e.LineStyle = 'none';
            e.Color  = 'k';
            e.Marker = '^';
            e.MarkerFaceColor = Colorsn(plotct,:);

            % Plot attached and detached fractions wrt. sep
            figure(5); clf; hold on
            e1 =  errorbar(Separations/(N_Kuhn*b),fa_bead_temp,faSE_bead_temp);
            e1.LineStyle = 'none';
            e1.Color = 'k';
            e1.Marker = 'o';
            e1.MarkerFaceColor = 'c';
            e1.MarkerSize = 8;

            e2 =  errorbar(Separations/(N_Kuhn*b),fd_bead_temp,fdSE_bead_temp);
            e2.LineStyle = 'none';
            e2.Color = 'k';
            e2.Marker = '^';
            e2.MarkerFaceColor = 'r';
            e2.MarkerSize = 8;

            e =  errorbar(Separations/(N_Kuhn*b),fa_meso_temp,faSE_meso_temp);
            e.LineStyle = 'none';
            e.Color = 'c';
            e.Marker = 'o';
            e.MarkerFaceColor = 'none';
            e.MarkerSize = 8;

            e =  errorbar(Separations/(N_Kuhn*b),fd_meso_temp,fdSE_meso_temp);
            e.LineStyle = 'none';
            e.Color = 'r';
            e.Marker = '^';
            e.MarkerFaceColor = 'none';
            e.MarkerSize = 8;

            set(gcf,'Position',[200 200 210 210])

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([1 1 1])

            xlabel('$d^*$ ($Nb$)','FontSize',FontSize,'Interpreter','latex')
            ylabel('$f$','FontSize',FontSize,'Interpreter','latex')

            FileTag = ['dt.',num2str(dt,'%.2e'),'.f vs d'];
            saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
            saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])
        end
    end
end

Xsize = 350;
Ysize = 350;

%% Polish emergent ka
figure(1)
set(gca,'FontSize',FontSize/1.5)
set(gcf,'Color','w')
pbaspect([1 1 1])

l = legend(Ptheor,LegendEntries);
l.FontSize = FontSize/2;
l.Interpreter = 'latex';

xlabel('$d^*$ ($Nb$)','FontSize',FontSize,'Interpreter','latex')
if NormalizeRates==1
    ylabel('$\bar{k}_a$ ($\tau_0^{-1}$)','FontSize',FontSize,'Interpreter','latex')
else
    ylabel('$\bar{k}_a$ (kHz)','FontSize',FontSize,'Interpreter','latex')
end

set(gcf,'Position',[100 100 Xsize Ysize])

FileTag = 'ka ';
saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])

%% Polish emergent kd
figure(2)
set(gca,'FontSize',FontSize/1.5)
set(gcf,'Color','w')
pbaspect([1 1 1])

xlabel('$d^*$ ($Nb$)','FontSize',FontSize,'Interpreter','latex')
if NormalizeRates==1
    ylabel('$\bar{k}_d$ ($\tau_0^{-1}$)','FontSize',FontSize,'Interpreter','latex')
else
    ylabel('$\bar{k}_d$ (kHz)','FontSize',FontSize,'Interpreter','latex')
end

ylim([0 yhi])

set(gcf,'Position',[100 100 Xsize Ysize])

FileTag = 'kd ';
saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])

%% Polish emergent fa
figure(3)
set(gca,'FontSize',FontSize/1.5)
set(gcf,'Color','w')
pbaspect([1 1 1])

xlabel('$d^*$ ($Nb$)','FontSize',FontSize,'Interpreter','latex')
ylabel('$f_a$','FontSize',FontSize,'Interpreter','latex')

ylim([0 1])

set(gcf,'Position',[100 100 Xsize Ysize])

FileTag = 'fa ';
saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])

%% Polish emergent fd
figure(4)
set(gca,'FontSize',FontSize/1.5)
set(gcf,'Color','w')
pbaspect([1 1 1])

xlabel('$d^*$ ($Nb$)','FontSize',FontSize,'Interpreter','latex')
ylabel('$f_d$','FontSize',FontSize,'Interpreter','latex')

ylim([0 1])

set(gcf,'Position',[100 100 Xsize Ysize])

FileTag = 'fd ';
saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])

%% Save fa and fd
figure(5)
FileTag = 'fa and fd ';
saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotDynamicsWrtTime(Override,Package,N_Kuhns,b,dts)

global AllDynamicsFileName FontSize OutputFolder

AddOn = 'Wrt. time.';

% MSD wrt time for each damper while sweeping dt
for dti = 1:length(dts)
    dt = dts(dti);

    DefineAssembledFileNames(dt);

    DynamicsWrtTime = load(AllDynamicsFileName,'-mat');

    time_bead = DynamicsWrtTime.time_bead;
    ka_bead = DynamicsWrtTime.ka_bead;
    kd_bead = DynamicsWrtTime.kd_bead;
    fa_bead = DynamicsWrtTime.fa_bead;
    fd_bead = DynamicsWrtTime.fd_bead;
    time_meso = DynamicsWrtTime.time_meso;
    ka_meso = DynamicsWrtTime.ka_meso;
    kd_meso = DynamicsWrtTime.kd_meso;
    fa_meso = DynamicsWrtTime.fa_meso;
    fd_meso = DynamicsWrtTime.fd_meso;

    window = 100;
    smoothka_bead = movmean(ka_bead/1e3,window,1);
    smoothkd_bead = movmean(kd_bead/1e3,window,1);
    smoothka_meso = movmean(ka_meso/1e3,window,1);
    smoothkd_meso = movmean(kd_meso/1e3,window,1);

    Npt = 25;
    ScatterSpacing = round((size(time_bead,1)-1)/Npt);
    rng = (1:ScatterSpacing:size(time_bead,1));
    DataSize = 10;

    for nk = 1:length(N_Kuhns)
        N_Kuhn = N_Kuhns(nk);

        FileTag = ['dt.',num2str(dt,'%.2e'),'.ka vs t'];
        SampName = [OutputFolder,'/',AddOn,FileTag,'.png'];
        if ~isfile(SampName) || Override

            PackageTemp = Package(Package(:,12)==N_Kuhn,:);
            Separations = unique(PackageTemp(:,14));

            ColorRange = (linspace(0,1,length(Separations)))';
            ColorFactd = [1 0 0];%[1 0.647 0];      %Orange-to-black as function of separation distance for detachment
            Colorsd = [flipud(ColorRange) flipud(ColorRange) flipud(ColorRange)].*ColorFactd;
            ColorFacta = [0 1 1];      %Cyan-to-black as function of separation distance for attachment
            Colorsa = [flipud(ColorRange) flipud(ColorRange) flipud(ColorRange)].*ColorFacta;

            % For plotting outputs wrt time
            figure(1); clf; hold on %ka
            figure(2); clf; hold on %kd
            figure(3); clf; hold on %Na
            figure(4); clf; hold on %Nd
            for sp = 1:length(Separations)
                Separation = Separations(sp);
                LegendEntries{sp} = ['$d^*$ = ',num2str(Separation/N_Kuhn/b,'%.2f')];

                %% Plot attachment rates wrt. time
                figure(1)
                s1 = scatter(time_bead(rng),smoothka_bead(rng,sp,nk));
                s1.MarkerEdgeColor = Colorsa(sp,:);
                s1.Marker = 'o';
                s1.SizeData = DataSize;

                p1(sp) = plot(time_meso,smoothka_meso(:,sp,nk));
                p1(sp).Color = Colorsa(sp,:);
                p1(sp).LineWidth = 1.5;
                p1(sp).LineStyle = '-';
                set(gca,'FontSize',FontSize)

                ylim([0 Inf])

                set(gca,'FontSize',FontSize/1.5)
                set(gcf,'Color','w')
                pbaspect([2 1 1])

                set(gcf,'Position',[200 200 250 250])

                xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
                ylabel('$k_a$ (kHz)','FontSize',FontSize,'Interpreter','latex')

                %% Plot detachment rates wrt. time
                figure(2)
                s2 = scatter(time_bead(rng),smoothkd_bead(rng,sp,nk));
                s2.MarkerEdgeColor = Colorsd(sp,:);
                s2.Marker = 'o';
                s2.SizeData = DataSize;

                p2(sp) = plot(time_meso,smoothkd_meso(:,sp,nk));
                p2(sp).Color = Colorsd(sp,:);
                p2(sp).LineWidth = 1.5;
                p2(sp).LineStyle = '-';

                ylim([0 Inf])

                set(gca,'FontSize',FontSize/1.5)
                set(gcf,'Color','w')
                pbaspect([2 1 1])

                set(gcf,'Position',[200 200 250 250])

                xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
                ylabel('$k_d$ (kHz)','FontSize',FontSize,'Interpreter','latex')

                %% Plot attached fraction of bonds wrt. time
                figure(3)
                s3 = scatter(time_bead(rng),fa_bead(rng,sp,nk)*100);
                s3.MarkerEdgeColor = Colorsa(sp,:);
                s3.Marker = 'o';
                s3.SizeData = DataSize;

                p3 = plot(time_meso,fa_meso(:,sp,nk)*100);
                p3.Color = Colorsa(sp,:);
                p3.LineWidth = 1.5;
                p3.LineStyle = '-';

                ylim([0 100])

                set(gca,'FontSize',FontSize/1.5)
                set(gcf,'Color','w')
                pbaspect([2 1 1])

                set(gcf,'Position',[200 200 250 250])

                xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
                ylabel('$f_a$ ($\%$)','FontSize',FontSize,'Interpreter','latex')

                %% Plot detached fraction of bonds wrt. time
                figure(4)
                s4 = scatter(time_bead(rng),fd_bead(rng,sp,nk)*100);
                s4.MarkerEdgeColor = Colorsd(sp,:);
                s4.Marker = 'o';
                s4.SizeData = DataSize;

                p4 = plot(time_meso,fd_meso(:,sp,nk)*100);
                p4.Color = Colorsd(sp,:);
                p4.LineWidth = 1.5;
                p4.LineStyle = '-';

                ylim([0 100])

                set(gca,'FontSize',FontSize/1.5)
                set(gcf,'Color','w')
                pbaspect([2 1 1])

                set(gcf,'Position',[200 200 250 250])

                xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
                ylabel('$f_d$ ($\%$)','FontSize',FontSize,'Interpreter','latex')
            end

            fig = figure(1);
            set(fig,'name',['dt = ',num2str(dt,'%.2e'),',ka vs t'])
            FileTag = ['dt.',num2str(dt,'%.2e'),'.ka vs t'];
            saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
            saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])
            l = legend(p1,LegendEntries);
            l.FontSize = FontSize/2;
            l.Interpreter = 'latex';

            fig = figure(2);
            set(fig,'name',['dt = ',num2str(dt,'%.2e'),',kd vs t'])
            FileTag = ['dt.',num2str(dt,'%.2e'),'.kd vs t'];
            saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
            saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])
            l = legend(p2,LegendEntries);
            l.FontSize = FontSize/2;
            l.Interpreter = 'latex';

            fig = figure(3);
            set(fig,'name',['dt = ',num2str(dt,'%.2e'),',fa vs t'])
            FileTag = ['dt.',num2str(dt,'%.2e'),'.fa vs t'];
            saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
            saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])

            fig = figure(4);
            set(fig,'name',['dt = ',num2str(dt,'%.2e'),',fd vs t'])
            FileTag = ['dt.',num2str(dt,'%.2e'),'.fd vs t'];
            saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
            saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])
        end
    end
end
close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AssembleEnsembleDynamicsData(Override,Package,...
    damps,N_Kuhns,Np,Ns,T,Nb,ka_in,kd_in,f0,dts,b,N)

global N_Kuhn...
    ka_ens_meso kd_ens_meso kaSE_ens_meso kdSE_ens_meso...
    Na_ens_meso Nd_ens_meso NaSE_ens_meso NdSE_ens_meso...
    fa_ens_meso fd_ens_meso faSE_ens_meso fdSE_ens_meso...
    ka_ens_bead kd_ens_bead kaSE_ens_bead kdSE_ens_bead...
    Na_ens_bead Nd_ens_bead NaSE_ens_bead NdSE_ens_bead...
    fa_ens_bead fd_ens_bead faSE_ens_bead fdSE_ens_bead...
    tau_a1_all_bead tau_a_all_bead tau_d_all_bead...
    tau_a1_all_meso tau_a_all_meso tau_d_all_meso...
    ka_all_meso kd_all_meso...
    Na_all_meso Nd_all_meso...
    fa_all_meso fd_all_meso...
    ka_all_bead kd_all_bead...
    Na_all_bead Nd_all_bead...
    fa_all_bead fd_all_bead...
    tau_a1_ens_meso tau_a1SE_ens_meso...
    tau_a_ens_meso tau_aSE_ens_meso...
    tau_d_ens_meso tau_dSE_ens_meso...
    tau_a1_ens_bead tau_a1SE_ens_bead...
    tau_a_ens_bead tau_aSE_ens_bead...
    tau_d_ens_bead tau_dSE_ens_bead...
    EnsembleDynamicsFileName AllDynamicsFileName

Folder = 'Output Plots';
if ~isfolder(Folder)
    mkdir(Folder)
end

% MSD wrt time for each damper while sweeping N
for dti = 1:length(dts)
    dt = dts(dti);

    DefineAssembledFileNames(dt);
    for dm = 1:length(damps)
        damp = damps(dm);

        for nk = 1:length(N_Kuhns)
            N_Kuhn = N_Kuhns(nk);


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

            for sp = 1:size(PackageTemp,1)
                Separation = PackageTemp(sp,14);

                if sp==1 && nk==1
                    InitializeEnsembleDynamics(N_Kuhns,size(PackageTemp,1));
                end

                %% Assemble the kinetics data by separation distance
                if ~isempty(PackageTemp) && (~isfile(EnsembleDynamicsFileName) ||...
                        Override)

                    % Pull in the Bead-spring data
                    BSOM = 0;
                    [time_bead,ka_bead,kd_bead,Na_bead,Nd_bead,fa_bead,fd_bead,...
                        tau_a1_bead,tau_a_bead,tau_d_bead,...
                        ka_ens_bead(sp,nk),kaSE_ens_bead(sp,nk),...
                        kd_ens_bead(sp,nk),kdSE_ens_bead(sp,nk),...
                        Na_ens_bead(sp,nk),NaSE_ens_bead(sp,nk),...
                        Nd_ens_bead(sp,nk),NdSE_ens_bead(sp,nk),...
                        fa_ens_bead(sp,nk),faSE_ens_bead(sp,nk),...
                        fd_ens_bead(sp,nk),fdSE_ens_bead(sp,nk),...
                        tau_a1_ens_bead(sp,nk),tau_a1SE_ens_bead(sp,nk),...
                        tau_a_ens_bead(sp,nk),tau_aSE_ens_bead(sp,nk),...
                        tau_d_ens_bead(sp,nk),tau_dSE_ens_bead(sp,nk)] =...
                        CalculateEnsembleDynamics(BSOM,...
                        Np,N,ka_in,kd_in,f0,dt,damp,N_Kuhn,b,Separation);

                    % Pull in the Mesoscale data
                    BSOM = 1;
                    [time_meso,ka_meso,kd_meso,Na_meso,Nd_meso,fa_meso,fd_meso,...
                        tau_a1_meso,tau_a_meso,tau_d_meso,...
                        ka_ens_meso(sp,nk),kaSE_ens_meso(sp,nk),...
                        kd_ens_meso(sp,nk),kdSE_ens_meso(sp,nk),...
                        Na_ens_meso(sp,nk),NaSE_ens_meso(sp,nk),...
                        Nd_ens_meso(sp,nk),NdSE_ens_meso(sp,nk),...
                        fa_ens_meso(sp,nk),faSE_ens_meso(sp,nk),...
                        fd_ens_meso(sp,nk),fdSE_ens_meso(sp,nk),...
                        tau_a1_ens_meso(sp,nk),tau_a1SE_ens_meso(sp,nk),...
                        tau_a_ens_meso(sp,nk),tau_aSE_ens_meso(sp,nk),...
                        tau_d_ens_meso(sp,nk),tau_dSE_ens_meso(sp,nk)] =...
                        CalculateEnsembleDynamics(BSOM,...
                        Np,N,ka_in,kd_in,f0,dt,damp,N_Kuhn,b,Separation);

                    if sp==1 && nk==1
                        InitializeAllDynamics(time_bead,time_meso,size(PackageTemp,1));
                    end

                    ka_all_bead(:,sp,nk) = ka_bead;
                    kd_all_bead(:,sp,nk) = kd_bead;
                    Na_all_bead(:,sp,nk) = Na_bead;
                    Nd_all_bead(:,sp,nk) = Nd_bead;
                    fa_all_bead(:,sp,nk) = fa_bead;
                    fd_all_bead(:,sp,nk) = fd_bead;

                    ka_all_meso(:,sp,nk) = ka_meso;
                    kd_all_meso(:,sp,nk) = kd_meso;
                    Na_all_meso(:,sp,nk) = Na_meso;
                    Nd_all_meso(:,sp,nk) = Nd_meso;
                    fa_all_meso(:,sp,nk) = fa_meso;
                    fd_all_meso(:,sp,nk) = fd_meso;

                    tau_a1_all_bead(1:length(tau_a1_bead),sp,nk) = tau_a1_bead;
                    tau_a_all_bead(1:length(tau_a_bead),sp,nk) = tau_a_bead;
                    tau_d_all_bead(1:length(tau_d_bead),sp,nk) = tau_d_bead;
                    tau_a1_all_meso(1:length(tau_a1_meso),sp,nk) = tau_a1_meso;
                    tau_a_all_meso(1:length(tau_a_meso),sp,nk) = tau_a_meso;
                    tau_d_all_meso(1:length(tau_d_meso),sp,nk) = tau_d_meso;
                end
            end
        end
    end
    if ~isfile(EnsembleDynamicsFileName) || Override
        EnsembleDynamics.ka_bead = ka_ens_bead;
        EnsembleDynamics.kaSE_bead = kaSE_ens_bead;
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

        EnsembleDynamics.ka_meso = ka_ens_meso;
        EnsembleDynamics.kaSE_meso = kaSE_ens_meso;
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

        EnsembleDynamics.tau_a1_meso = tau_a1_ens_meso;
        EnsembleDynamics.tau_a1SE_meso = tau_a1SE_ens_meso;
        EnsembleDynamics.tau_a_meso = tau_a_ens_meso;
        EnsembleDynamics.tau_aSE_meso = tau_aSE_ens_meso;
        EnsembleDynamics.tau_d_meso = tau_d_ens_meso;
        EnsembleDynamics.tau_dSE_meso = tau_dSE_ens_meso;
        EnsembleDynamics.tau_a1_bead = tau_a1_ens_bead;
        EnsembleDynamics.tau_a1SE_bead = tau_a1SE_ens_bead;
        EnsembleDynamics.tau_a_bead = tau_a_ens_bead;
        EnsembleDynamics.tau_aSE_bead = tau_aSE_ens_bead;
        EnsembleDynamics.tau_d_bead = tau_d_ens_bead;
        EnsembleDynamics.tau_dSE_bead = tau_dSE_ens_bead;

        save(EnsembleDynamicsFileName,'-struct','EnsembleDynamics');
    end

    if ~isfile(AllDynamicsFileName) || Override
        AllDynamics.time_bead = time_bead;
        AllDynamics.ka_bead = ka_all_bead;
        AllDynamics.kd_bead = kd_all_bead;
        AllDynamics.Na_bead = Na_all_bead;
        AllDynamics.Nd_bead = Nd_all_bead;
        AllDynamics.fa_bead = fa_all_bead;
        AllDynamics.fd_bead = fd_all_bead;
        AllDynamics.tau_a1_bead = tau_a1_all_bead;
        AllDynamics.tau_a_bead = tau_a_all_bead;
        AllDynamics.tau_d_bead = tau_d_all_bead;

        AllDynamics.time_meso = time_meso;
        AllDynamics.ka_meso = ka_all_meso;
        AllDynamics.kd_meso = kd_all_meso;
        AllDynamics.Na_meso = Na_all_meso;
        AllDynamics.Nd_meso = Nd_all_meso;
        AllDynamics.fa_meso = fa_all_meso;
        AllDynamics.fd_meso = fd_all_meso;
        AllDynamics.tau_a1_meso = tau_a1_all_meso;
        AllDynamics.tau_a_meso = tau_a_all_meso;
        AllDynamics.tau_d_meso = tau_d_all_meso;

        save(AllDynamicsFileName,'-struct','AllDynamics');
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DefineAssembledFileNames(dt)

global OutputDir EnsembleDynamicsFileName AllDynamicsFileName...
    FittedDynamicsFileName

EnsembleDynamicsFileName = [OutputDir,'EnsembleDynamics.dt-',...
    num2str(dt,'%.2e'),'.m'];
AllDynamicsFileName = [OutputDir,'AllDynamicswrtTime.dt-',...
    num2str(dt,'%.2e'),'.m'];
FittedDynamicsFileName = [OutputDir,'FittedDynamics.dt-',...
    num2str(dt,'%.2e'),'.m'];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitializeAllDynamics(time_bead,time_meso,NoSep)

global ka_all_meso kd_all_meso...
    Na_all_meso Nd_all_meso...
    fa_all_meso fd_all_meso...
    ka_all_bead kd_all_bead...
    Na_all_bead Nd_all_bead...
    fa_all_bead fd_all_bead...
    tau_a1_all_bead tau_a_all_bead...
    tau_d_all_bead tau_a1_all_meso...
    tau_a_all_meso tau_d_all_meso

ka_all_bead = zeros(size(time_bead,1),NoSep);
kd_all_bead = zeros(size(time_bead,1),NoSep);
Na_all_bead = zeros(size(time_bead,1),NoSep);
Nd_all_bead = zeros(size(time_bead,1),NoSep);
fa_all_bead = zeros(size(time_bead,1),NoSep);
fd_all_bead = zeros(size(time_bead,1),NoSep);

ka_all_meso = zeros(size(time_meso,1),NoSep);
kd_all_meso = zeros(size(time_meso,1),NoSep);
Na_all_meso = zeros(size(time_meso,1),NoSep);
Nd_all_meso = zeros(size(time_meso,1),NoSep);
fa_all_meso = zeros(size(time_meso,1),NoSep);
fd_all_meso = zeros(size(time_meso,1),NoSep);

tau_a1_all_bead = zeros(5e3,NoSep);
tau_a_all_bead = zeros(5e3,NoSep);
tau_d_all_bead = zeros(5e3,NoSep);
tau_a1_all_meso = zeros(5e3,NoSep);
tau_a_all_meso = zeros(5e3,NoSep);
tau_d_all_meso = zeros(5e3,NoSep);

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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time,ka,kd,Na,Nd,fa,fd,tau_a1,tau_a,tau_d,...
    mean_ka,se_ka,mean_kd,se_kd,...
    mean_Na,se_Na,mean_Nd,se_Nd,...
    mean_fa,se_fa,mean_fd,se_fd,...
    mean_tau_a1,se_tau_a1,...
    mean_tau_a,se_tau_a,...
    mean_tau_d,se_tau_d] =...
    CalculateEnsembleDynamics(BSOM,Np,N,ka_in,kd_in,f0,...
    dt,damp,N_Kuhn,b,Separation)

global BeadSpringOrMeso TimeStretchDataFileName BondKineticsDataFileName

BeadSpringOrMeso = BSOM;
EdgeFactor = 1;
InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dt,damp,...
    N_Kuhn,b,Separation);
SetDirAndFileNames;
DefineCompiledFileNames;

TimeAndStretch = load(TimeStretchDataFileName,'-mat');
time= TimeAndStretch.time;

BondKinetics = load(BondKineticsDataFileName,'-mat');
ka = BondKinetics.ka;
kd = BondKinetics.kd;
Na = BondKinetics.Na;
Nd = BondKinetics.Nd;

fa = Na./(Na+Nd);
fd = Nd./(Na+Nd);

SSPcnt = 20; SSIndx = round(SSPcnt/100*size(fa,1));

mean_fa = mean(fa(SSIndx:end,:),1);
se_fa = std(fa(SSIndx:end,:),1,1)/sqrt(size(fa(SSIndx:end,:),1));
mean_fd = mean(fd(SSIndx:end,:),1);
se_fd = std(fd(SSIndx:end,:),1,1)/sqrt(size(fd(SSIndx:end,:),1));
mean_ka = mean(ka,1)/1e3; se_ka = std(ka,1,1)/sqrt(size(ka,1))/1e3;
mean_kd = mean(kd,1)/1e3; se_kd = std(kd,1,1)/sqrt(size(kd,1))/1e3;
mean_Na = mean(Na,1); se_Na = std(Na,1,1)/sqrt(size(Na,1));
mean_Nd = mean(Nd,1); se_Nd = std(Nd,1,1)/sqrt(size(Nd,1));

% Attached and detached bond lifetimes
tau_a1 = BondKinetics.FirstAttachment;
tau_a = BondKinetics.AttachedLifetimes;
tau_d = BondKinetics.DetachedLifetimes;

mean_tau_a1 = mean(tau_a1);
se_tau_a1 = std(tau_a1)/sqrt(length(tau_a1));

mean_tau_a = mean(tau_a);
se_tau_a = std(tau_a)/sqrt(length(tau_a));

mean_tau_d = mean(tau_d);
se_tau_d = std(tau_d)/sqrt(length(tau_d));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = FitThePrefactor(x,y,A_nom,...
    D,NormLength,N_Kuhn,eaOverkbT,sigma_a_theor)

% global NoFitAttempts
NoFitAttempts = 50;
ReductionFactor = 0.98;
NoPts = 21;
R2 = 0;
ct = 0;
wb5 = waitbar(0,'Fitting Prefactor...');
figure(100)
while R2<0.999
    ct = ct+1;
    waitbar(ct/NoFitAttempts,wb5,'Fitting Prefactor...')

    if ct==1
        sig_A = 1/8*A_nom;
    else
        sig_A = sig_A*ReductionFactor;
    end
    A_rng = [(linspace(A_nom-sig_A,A_nom+sig_A,NoPts))';A_nom];
    A_rng = unique(A_rng); 
    A_rng(A_rng<=0) = [];

    Perms = zeros(length(A_rng),1);
    i = 0;
    for mi=1:length(A_rng)
        i = i+1;
        Perms(i,1) = A_rng(mi);
    end
    R2_all = zeros(size(Perms,1),1);
    A_all = zeros(size(Perms,1),1);

    for i=1:size(Perms,1)
        A_all(i) = Perms(i,1);

        ka0_theor = A_all(i)*3/4/pi*D/(NormLength^2)*N_Kuhn^(-3/2)*exp(-eaOverkbT);
        ft = ka0_theor*exp(-(x.^2)/(sigma_a_theor^2));

        RSS = sum((y-ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_all(i) = 1-RSS/TSS;
    end
    diff = abs(R2_all-1);
    indx = find(diff==min(diff),1,'first');

    A_nom = A_all(indx);
    R2 = R2_all(indx);

    figure(100); clf; hold on
    scatter(x,y,'k','filled')
    ka0_theor = A_nom*3/4/pi*D/(NormLength^2)*N_Kuhn^(-3/2)*exp(-eaOverkbT);
    yplot = ka0_theor*exp(-(x.^2)/(sigma_a_theor^2));
    plot(x,yplot,'k--')
    if ct>NoFitAttempts
        break;
    end
end
A = A_nom;
figure(100); close
close(wb5)

end