function PostProcessData(Package,OverridePostProccess,ToggleDynamics,...
    LC,DC,FPE)

% Computes and plots stress-strain for every iteration of simulation

global LineWidth TurnOnDynamics FontSize NoFitAttempts...
    LengthConversion DamperConversion DataSize...
    OutputFolder FitPowerOrExponent PlotScalingTheory

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
% BeadSpringOrMeso = BSOM;
FitPowerOrExponent = FPE;
PlotScalingTheory = 1;
NormalizeRates = 1;

% Samples = unique(Package(:,1));
Np = unique(Package(:,2));        %Number of molecules
Ns = unique(Package(:,3));        %Number of stickeres per tether site
T = unique(Package(:,4));         %Temperature
Nb = unique(Package(:,5));        %Contour length of chains
ka_in = unique(Package(:,6));        %Activation energy of association
kd_in = unique(Package(:,7));        %Activation energy of dissociation
f0 = unique(Package(:,8));       %Force sensitivity to dissociation
dt = unique(Package(:,9));       %timestep size
DiffCoeffs = unique(Package(:,11));
N_Kuhns = unique(Package(:,12));
b = unique(Package(:,13));
N = Np*Ns;

%% Sweepign parameters are N_Kuhn and b
damps = unique(Package(:,10));    %damping coefficient in units [mass/time]

ka_bar =  zeros(length(ka_in),1);
kd_bar =  zeros(length(ka_in),1);
ka_SE =  zeros(length(ka_in),1);
kd_SE =  zeros(length(ka_in),1);

for i=1:length(ka_in)
    ka = ka_in(i);

    NormLength = b*LengthConversion;
    tau0 = (NormLength^2)/DiffCoeffs;
    eaStar = -log(tau0*ka);

    %% Assemble the Ensemble Data
    AssembleEnsembleDynamicsData(OverridePostProccess,Package,damps,N_Kuhns,Np,Ns,T,Nb,...
        ka,kd_in,f0,dt,b,N);

    % Dynamics wrt. time
    [ka_bar(i,1),ka_SE(i,1),kd_bar(i,1),kd_SE(i,1)] = ...
        PlotDynamicsWrtTime(OverridePostProccess,Package,N_Kuhns,b,eaStar,...
        DiffCoeffs,kd_in);
end
close all

%% Plot outputs swept over activation energy
NormLength = b*LengthConversion;
tau0 = (NormLength^2)/DiffCoeffs;
eaStar = -log(tau0*ka_in);

ka_bar(eaStar==10) = [];
kd_bar(eaStar==10) = [];
ka_SE(eaStar==10) = [];
kd_SE(eaStar==10) = [];
eaStar(eaStar==10) = [];

figure(1); clf; hold on

yyaxis left
if NormalizeRates==1
    sa = errorbar(eaStar,ka_bar*tau0,ka_SE*tau0);
else
    sa = errorbar(eaStar,ka_bar/1e6,ka_SE/1e6);
end
sa.Marker = 'o';
sa.LineStyle = 'none';
sa.MarkerFaceColor = 'auto';

eaPlot = (linspace(min(eaStar),10,100))';
kaPlot = DiffCoeffs/(NormLength^2)*exp(-eaPlot);

if NormalizeRates==1
    p = plot(eaPlot,kaPlot*tau0);
else
    p = plot(eaPlot,kaPlot/1e6);
end
p.Color = sa.Color;

yyaxis right

if NormalizeRates==1
    sd = errorbar(eaStar,kd_bar*tau0*1e3,kd_SE*tau0*1e3);
else
    sd = errorbar(eaStar,kd_bar/1e3,kd_SE/1e3);
end
sd.Marker = '^';
sd.LineStyle = 'none';
sd.MarkerFaceColor = 'auto';

if NormalizeRates==1
    p = plot([min(eaStar),10],[kd_in kd_in]*tau0);
else
    p = plot([min(eaStar),10],[kd_in kd_in]/1e3);
end
p.Color = sd.Color;

set(gca,'XScale','log')
set(gca,'FontSize',FontSize/1.5)

xlabel('$\varepsilon_a/k_b T$','FontSize',FontSize,'Interpreter','latex')
yyaxis left
if NormalizeRates==1
    ylabel('$\bar k_a^*$ ($\tau_0^{-1}$)','FontSize',FontSize,'Interpreter','latex')
    ylim([0 1.1])
else
    ylabel('$\bar k_a$ (MHz)','FontSize',FontSize,'Interpreter','latex')
    ylim([0 300])
end
yyaxis right
if NormalizeRates==1
    ylabel('$\bar k_d^*$ ($\tau_0^{-1}$)','FontSize',FontSize,'Interpreter','latex')
    ylim([0 1.1])
%     ylim([0 2*kd_in*tau0])
else
    ylabel('$\bar k_d$ (kHz)','FontSize',FontSize,'Interpreter','latex')
    ylim([0 2*kd_in/1e3])
end

set(gcf,'Color','w')
pbaspect([1 1 1])

l = legend('$k_a$','$k_a^{ap}$','$k_d$','$k_d^{ap}$');
l.FontSize = FontSize/2;
l.Location = 'NorthEast';
l.Interpreter = 'latex';

if NormalizeRates==1
    FileTag = 'Output Plots/ka_star and kd_star.png';
else
    FileTag = 'Output Plots/ka and kd.png';
end
saveas(gcf,[FileTag,'.png'])
saveas(gcf,[FileTag,'.fig'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ka_bar,ka_SE,kd_bar,kd_SE] = ...
    PlotDynamicsWrtTime(Override,Package,N_Kuhns,b,eaStar,D,kd_ap)

global AllDynamicsFileName FontSize OutputFolder LengthConversion

AddOn = 'Wrt. time.';
DynamicsWrtTime = load(AllDynamicsFileName,'-mat');

ka_ap = D/(LengthConversion*b)^2*exp(-eaStar)/1e3;

time = DynamicsWrtTime.time;
ka = DynamicsWrtTime.ka;
kd = DynamicsWrtTime.kd;
fa = DynamicsWrtTime.fa;
fd = DynamicsWrtTime.fd;

ka_bar = mean(ka);
kd_bar = mean(kd);
ka_SE = std(ka,1)/sqrt(length(ka));
kd_SE = std(kd,1)/sqrt(length(kd));

window = 200;
smoothka = movmean(ka/1e3,window,1);
smoothkd = movmean(kd/1e3,window,1);

Npt = 50;
ScatterSpacing = round((size(time,1)-1)/Npt);
rng = (1:ScatterSpacing:size(time,1));
DataSize = 10;

% MSD wrt time for each damper while sweeping N
for nk = 1:length(N_Kuhns)
    N_Kuhn = N_Kuhns(nk);

    FileTag = ['Nk.',num2str(N_Kuhn),'.ka vs t'];
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
        figure(1); clf; hold on %ka and kd
        figure(2); clf; hold on %fa and fd
        for sp = 1%1:length(Separations)    %It is identical across all 4 separations GOOD
            Separation = Separations(sp);
%             LegendEntries{sp} = ['$d^*$ = ',num2str(Separation/N_Kuhn/b,'%.2f')];

            %% Plot attachment rates wrt. time
            figure(1)
            yyaxis left
            plot([0 4e-7],[ka_ap ka_ap],'k--');

            s1 = scatter(time(rng),smoothka(rng,sp,nk),'filled');
            s1.MarkerEdgeColor = 'k';
%             s1.MarkerFaceColor = Colorsa(sp,:);
            s1.Marker = 'o';
            s1.SizeData = 4*DataSize;

            ylim([0 ka_ap*1.5])

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([1 1 1])

%             set(gcf,'Position',[200 200 250 250])

            %% Plot detachment rates wrt. time
%             figure(2)
            yyaxis right
%             plot([0 4e-7],[kd_ap kd_ap]/1e3,'k--');

            s2 = scatter(time(rng),smoothkd(rng,sp,nk),'filled');
            s2.MarkerEdgeColor = 'k';
%             s2.MarkerFaceColor = Colorsd(sp,:);
            s2.Marker = 'diamond';
            s2.SizeData = 4*DataSize;

            ylim([0 kd_ap/1e3*1.5])

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([1 1 1])

%             set(gcf,'Position',[200 200 250 250])

            xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
            yyaxis left
            ylabel('$k_a$ (kHz)','FontSize',FontSize,'Interpreter','latex')
            yyaxis right
            ylabel('$k_d$ (kHz)','FontSize',FontSize,'Interpreter','latex')

            title(['$\varepsilon_a^*$ = ',num2str(eaStar,'%.2f')],...
                'FontSize',FontSize/1.5,'Interpreter','latex')

            %% Plot attached fraction of bonds wrt. time
            figure(2)
            yyaxis left
            plot([0 4e-7],[100 100],'k--');

            s3 = scatter(time(rng),fa(rng,sp,nk)*100,'filled');
            s3.MarkerEdgeColor = 'k';%Colorsa(sp,:);
            s3.Marker = 'o';
            s3.SizeData = 4*DataSize;

            ylim([0 105])


            %% Plot detached fraction of bonds wrt. time
            yyaxis right

            s4 = scatter(time(rng),fd(rng,sp,nk)*100,'filled');
            s4.MarkerEdgeColor = 'k';%Colorsd(sp,:);
            s4.Marker = 'diamond';
            s4.SizeData = 4*DataSize;

            ylim([0 105])

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([2 1 1])

            set(gcf,'Position',[200 200 400 400])

            xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
            yyaxis left
            ylabel('$f_a$ ($\%$)','FontSize',FontSize,'Interpreter','latex')
            yyaxis right
            ylabel('$f_d$ ($\%$)','FontSize',FontSize,'Interpreter','latex')

            title(['$\varepsilon_a^*$ = ',num2str(eaStar,'%.2f')],...
                'FontSize',FontSize/1.5,'Interpreter','latex')
        end

        fig = figure(1);
        set(fig,'name',['N = ',num2str(N_Kuhn),',ka and kd vs t'])
        FileTag = ['Nk.',num2str(N_Kuhn),'.ka and kd vs t'];
        saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
        saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])

        fig = figure(2);
        set(fig,'name',['N = ',num2str(N_Kuhn),',fa and fd vs t'])
        FileTag = ['Nk.',num2str(N_Kuhn),'.fa and fd vs t'];
        saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
        saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])
    end
end
close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AssembleEnsembleDynamicsData(Override,Package,...
    damps,N_Kuhns,Np,Ns,T,Nb,ka_in,kd_in,f0,dts,b,N)

global N_Kuhn...
    ka_ens kd_ens kaSE_ens kdSE_ens...
    Na_ens Nd_ens NaSE_ens NdSE_ens...
    fa_ens fd_ens faSE_ens fdSE_ens...
    tau_a1_all tau_a_all tau_d_all...
    ka_all kd_all...
    Na_all Nd_all...
    fa_all fd_all...
    tau_a1_ens tau_a1SE_ens...
    tau_a_ens tau_aSE_ens...
    tau_d_ens tau_dSE_ens...
    EnsembleDynamicsFileName AllDynamicsFileName

Folder = 'Output Plots';
if ~isfolder(Folder)
    mkdir(Folder)
end

DefineAssembledFileNames(ka_in);

% MSD wrt time for each damper while sweeping N
for dm = 1:length(damps)
    damp = damps(dm);
    dt = dts(dm);

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
                [time,ka,kd,Na,Nd,fa,fd,tau_a1,tau_a,tau_d,...
                    ka_ens(sp,nk),kaSE_ens(sp,nk),...
                    kd_ens(sp,nk),kdSE_ens(sp,nk),...
                    Na_ens(sp,nk),NaSE_ens(sp,nk),...
                    Nd_ens(sp,nk),NdSE_ens(sp,nk),...
                    fa_ens(sp,nk),faSE_ens(sp,nk),...
                    fd_ens(sp,nk),fdSE_ens(sp,nk),...
                    tau_a1_ens(sp,nk),tau_a1SE_ens(sp,nk),...
                    tau_a_ens(sp,nk),tau_aSE_ens(sp,nk),...
                    tau_d_ens(sp,nk),tau_dSE_ens(sp,nk)] =...
                    CalculateEnsembleDynamics(Np,N,ka_in,kd_in,f0,dt,damp,...
                    N_Kuhn,b,Separation);

                if sp==1 && nk==1
                    InitializeAllDynamics(time,size(PackageTemp,1));
                end

                ka_all(:,sp,nk) = ka;
                kd_all(:,sp,nk) = kd;
                Na_all(:,sp,nk) = Na;
                Nd_all(:,sp,nk) = Nd;
                fa_all(:,sp,nk) = fa;
                fd_all(:,sp,nk) = fd;

                tau_a1_all(1:length(tau_a1),sp,nk) = tau_a1;
                tau_a_all(1:length(tau_a),sp,nk) = tau_a;
                tau_d_all(1:length(tau_d),sp,nk) = tau_d;
            end
        end
    end
end

if ~isfile(EnsembleDynamicsFileName) || Override
    EnsembleDynamics.ka = ka_ens;
    EnsembleDynamics.kaSE = kaSE_ens;
    EnsembleDynamics.kd = kd_ens;
    EnsembleDynamics.kdSE = kdSE_ens;
    EnsembleDynamics.Na = Na_ens;
    EnsembleDynamics.NaSE = NaSE_ens;
    EnsembleDynamics.Nd = Nd_ens;
    EnsembleDynamics.NdSE = NdSE_ens;
    EnsembleDynamics.fa = fa_ens;
    EnsembleDynamics.faSE = faSE_ens;
    EnsembleDynamics.fd = fd_ens;
    EnsembleDynamics.fdSE = fdSE_ens;

    EnsembleDynamics.tau_a1 = tau_a1_ens;
    EnsembleDynamics.tau_a1SE = tau_a1SE_ens;
    EnsembleDynamics.tau_a = tau_a_ens;
    EnsembleDynamics.tau_aSE = tau_aSE_ens;
    EnsembleDynamics.tau_d = tau_d_ens;
    EnsembleDynamics.tau_dSE = tau_dSE_ens;

    save(EnsembleDynamicsFileName,'-struct','EnsembleDynamics');
else
    EnsembleDynamics = load(EnsembleDynamicsFileName,'-mat');

    ka_ens = EnsembleDynamics.ka;
    kaSE_ens = EnsembleDynamics.kaSE;
    kd_ens = EnsembleDynamics.kd;
    kdSE_ens =  EnsembleDynamics.kdSE;
    Na_ens = EnsembleDynamics.Na;
    NaSE_ens = EnsembleDynamics.NaSE;
    Nd_ens = EnsembleDynamics.Nd;
    NdSE_ens = EnsembleDynamics.NdSE;
    fa_ens = EnsembleDynamics.fa;
    faSE_ens = EnsembleDynamics.faSE;
    fd_ens = EnsembleDynamics.fd;
    fdSE_ens = EnsembleDynamics.fdSE;

    tau_a1_ens = EnsembleDynamics.tau_a1;
    tau_a1SE_ens = EnsembleDynamics.tau_a1SE;
    tau_a_ens = EnsembleDynamics.tau_a;
    tau_aSE_ens = EnsembleDynamics.tau_aSE;
    tau_d_ens = EnsembleDynamics.tau_d;
    tau_dSE_ens = EnsembleDynamics.tau_dSE;
end

if ~isfile(AllDynamicsFileName) || Override
    AllDynamics.time = time;
    AllDynamics.ka = ka_all;
    AllDynamics.kd = kd_all;
    AllDynamics.Na = Na_all;
    AllDynamics.Nd = Nd_all;
    AllDynamics.fa = fa_all;
    AllDynamics.fd = fd_all;
    AllDynamics.tau_a1 = tau_a1_all;
    AllDynamics.tau_a = tau_a_all;
    AllDynamics.tau_d = tau_d_all;

    save(AllDynamicsFileName,'-struct','AllDynamics');
else
    AllDynamics = load(AllDynamicsFileName,'-mat');

    ka_all = AllDynamics.ka;
    kd_all = AllDynamics.kd;
    Na_all = AllDynamics.Na;
    Nd_all = AllDynamics.Nd;
    fa_all = AllDynamics.fa;
    fd_all = AllDynamics.fd;
    tau_a1_all = AllDynamics.tau_a1;
    tau_a_all = AllDynamics.tau_a;
    tau_d_all = AllDynamics.tau_d;
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
function InitializeAllDynamics(time,NoSep)

global ka_all kd_all...
    Na_all Nd_all...
    fa_all fd_all...
    tau_a1_all tau_a_all tau_d_all 

ka_all = zeros(size(time,1),NoSep);
kd_all = zeros(size(time,1),NoSep);
Na_all = zeros(size(time,1),NoSep);
Nd_all = zeros(size(time,1),NoSep);
fa_all = zeros(size(time,1),NoSep);
fd_all = zeros(size(time,1),NoSep);

tau_a1_all = zeros(5e3,NoSep);
tau_a_all = zeros(5e3,NoSep);
tau_d_all = zeros(5e3,NoSep);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitializeEnsembleDynamics(N_Kuhns,NoSep)

global ka_ens_meso kd_ens_meso kaSE_ens_meso kdSE_ens_meso...
    Na_ens_meso Nd_ens_meso NaSE_ens_meso NdSE_ens_meso...
    fa_ens_meso fd_ens_meso faSE_ens_meso fdSE_ens_meso...
    ka_ens kd_ens kaSE_ens kdSE_ens...
    Na_ens Nd_ens NaSE_ens NdSE_ens...
    fa_ens fd_ens faSE_ens fdSE_ens...
    tau_a1_ens tau_a1SE_ens...
    tau_a_ens tau_aSE_ens...
    tau_d_ens tau_dSE_ens...
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

ka_ens = zeros(NoSep,length(N_Kuhns));
kd_ens = zeros(NoSep,length(N_Kuhns));
kaSE_ens = zeros(NoSep,length(N_Kuhns));
kdSE_ens = zeros(NoSep,length(N_Kuhns));
Na_ens = zeros(NoSep,length(N_Kuhns));
Nd_ens = zeros(NoSep,length(N_Kuhns));
NaSE_ens = zeros(NoSep,length(N_Kuhns));
NdSE_ens = zeros(NoSep,length(N_Kuhns));
fa_ens = zeros(NoSep,length(N_Kuhns));
fd_ens = zeros(NoSep,length(N_Kuhns));
faSE_ens = zeros(NoSep,length(N_Kuhns));
fdSE_ens = zeros(NoSep,length(N_Kuhns));

tau_a1_ens = zeros(NoSep,length(N_Kuhns));
tau_a1SE_ens = zeros(NoSep,length(N_Kuhns));
tau_a_ens = zeros(NoSep,length(N_Kuhns));
tau_aSE_ens = zeros(NoSep,length(N_Kuhns));
tau_d_ens = zeros(NoSep,length(N_Kuhns));
tau_dSE_ens = zeros(NoSep,length(N_Kuhns));

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
    CalculateEnsembleDynamics(Np,N,ka_in,kd_in,f0,...
    dt,damp,N_Kuhn,b,Separation)

global BondKineticsDataFileName


% BeadSpringOrMeso = BSOM;
EdgeFactor = 1;
InputScript(EdgeFactor,1,Np,N,ka_in,kd_in,f0,dt,damp,...
    N_Kuhn,b,Separation);
SetDirAndFileNames;
DefineCompiledFileNames;

BondKinetics = load(BondKineticsDataFileName,'-mat');
time = BondKinetics.time;
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
