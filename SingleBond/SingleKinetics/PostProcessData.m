function PostProcessData(Package,OverridePostProcess,ToggleDynamics,...
    LC,DC,BSOM,FPE,ASD,CF,OF)

% Computes and plots stress-strain for every iteration of simulation

global LineWidth TurnOnDynamics FontSize NoFitAttempts...
    LengthConversion DamperConversion DataSize BeadSpringOrMeso...
    OutputFolder FitPowerOrExponent PlotScalingTheory PlotOnlyScalingTheory...
    timenorm NormalizeRates AppendScalingData CurrentFolder

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
PlotScalingTheory = 1;
PlotOnlyScalingTheory = 1;  % 1 to show only the continuous scaling theory curve with A treated as a fitting parameter
AppendScalingData = ASD;
CurrentFolder = CF;
OutputFolder = OF;

NormalizeRates = 1;         % Set 1 to normalize rates by 1/tau0

% Samples = unique(Package(:,1));
Np = unique(Package(:,2));        %Number of molecules
Ns = unique(Package(:,3));        %Number of stickeres per tether site
T = unique(Package(:,4));         %Temperature
Nb = unique(Package(:,5));        %Contour length of chains
ka_in = unique(Package(:,6));        %Activation energy of association
kd_in = unique(Package(:,7));        %Activation energy of dissociation
f0 = unique(Package(:,8));       %Force sensitivity to dissociation
dt = unique(Package(:,9));       %timestep size
damp = unique(Package(:,10));    %damping coefficient in units [mass/time]
DiffCoeffs = unique(Package(:,11));
N_Kuhns = unique(Package(:,12));
b = unique(Package(:,13));
Separations = unique(Package(:,14));
N = Np*Ns;

%% Sweepign parameters are N_Kuhn and ka
R2_meso_vs_bead = zeros(length(N_Kuhns),length(ka_in));

A_bead = zeros(length(N_Kuhns),length(ka_in));
A_meso = zeros(length(N_Kuhns),length(ka_in));
dA_bead = zeros(length(N_Kuhns),length(ka_in));
dA_meso = zeros(length(N_Kuhns),length(ka_in));
R2_bead = zeros(length(N_Kuhns),length(ka_in));
R2_meso = zeros(length(N_Kuhns),length(ka_in));

R2a1 = zeros(length(N_Kuhns),length(ka_in));
R2a = zeros(length(N_Kuhns),length(ka_in));
R2d = zeros(length(N_Kuhns),length(ka_in));

for i=1:length(ka_in)
    ka = ka_in(i);

    NormLength = b*LengthConversion;
    tau0 = (NormLength^2)/DiffCoeffs;

    if NormalizeRates==1
        timenorm = 1/tau0;
    else
        timenorm = 1e3;
    end

    EdgeFactor = Nb;
    InputScript(EdgeFactor,1,Np,N,ka,kd_in,f0,dt,damp,...
        N_Kuhns(1),b,Separations(1));
    SetDirAndFileNames;

    %% Assemble the Ensemble Data
    AssembleEnsembleDynamicsData(OverridePostProcess,Package,damp,N_Kuhns,Np,Ns,T,Nb,...
        ka,kd_in,f0,dt,b,N,AppendScalingData);

    %% Plot Outputs
    % Histograms of bond lifetimes  
    [R2a1(:,i),R2a(:,i),R2d(:,i)] = PlotBondLifetimes(OverridePostProcess,...
        Package,N_Kuhns,b,DiffCoeffs,ka);
    close all

    % Dynamics wrt. time
    PlotDynamicsWrtTime(OverridePostProcess,Package,N_Kuhns,b);
    close all

    % Ensemble average dynamics
    OverridePostProcess = 1;
    [A_bead(:,i),dA_bead(:,i),A_meso(:,i),dA_meso(:,i),...
    R2_bead(:,i),R2_meso(:,i),R2_meso_vs_bead(:,i)] = ...
        PlotEnsembleAvgDynamics(OverridePostProcess,Package,N_Kuhns,b,ka,...
        kd_in,DiffCoeffs);
    close all
    OverridePostProcess = 0;
end
close all

%% Plot Theory on top of existing bead-spring and mesoscale data
OverridePostProcess = 1;

PlotScalingTheoryPhaseSpace(OverridePostProcess,b,DiffCoeffs)

%% Plot outputs swept over activation energy
PlotOutputsOverActivationEnergy(N_Kuhns,b,DiffCoeffs,ka_in,...
    A_bead,dA_bead,R2_bead,A_meso,dA_meso,R2_meso,R2_meso_vs_bead);
close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotScalingTheoryPhaseSpace(OverridePostProcess,b,D)

global FontSize LengthConversion 

NormLength = b*LengthConversion;
tau0 = (NormLength^2)/D;

ea_rng = (logspace(-2,1,15))';
lambda_rng = (linspace(0,1,15))';
[ea_grid,lambda_grid] = meshgrid(ea_rng,lambda_rng);

N_Kuhns = 12:36;

Folder = 'Output Plots/Scaling Theory';
if ~isfolder(Folder)
    mkdir(Folder)
end

MovieName = [Folder,'/Scaling Theory Phase Diagram'];
if ~isfile([MovieName,'.avi']) || OverridePostProcess==1
    v1 = VideoWriter(MovieName); v1.FrameRate = 10; open(v1)

    for nk_i = 1:length(N_Kuhns)
        figure(4); clf; hold on

        N_Kuhn = N_Kuhns(nk_i);
        dist_grid = lambda_grid*(sqrt(2*N_Kuhn)*b);

        sigma_a_theor = (3/2)^(-1/2)*sqrt(2*N_Kuhn)*b;
        k_encounter = ((3/2/pi)^(3/2))*D/(NormLength^2)*N_Kuhn^(-3/2);
        ka0_theor = k_encounter*exp(-ea_grid);
        ka_plot = ka0_theor.*exp(-(dist_grid.^2)/(sigma_a_theor^2));

        cmax = ((3/2/pi)^(3/2))*D/(NormLength^2)*(min(N_Kuhns))^(-3/2)...
            *exp(-min(min(ea_grid)))*exp(-(0.^2)/(((3/2)^(-1/2)*sqrt(min(N_Kuhns)*b)^2)));

        s = surf(ea_grid,lambda_grid,ka_plot*tau0);
        s.FaceColor = 'interp';

        c = colorbar;
        c.Ticks = [1e-3 2e-3 4e-3 8e-3];
        set(gca,'FontSize',FontSize/1.5)
        set(gca,'colorscale','log')
        clab = c.Label;
        clab.String = '$\bar{k}_a\tau_0$';
        clab.Interpreter = 'latex';
        clab.FontSize = FontSize;

        caxis([1e-3 cmax*tau0])
        %     set(gca,'Colorscale','log')
        set(gca,'xscale','log')

        xlabel('$\varepsilon_a/k_b T$','FontSize',FontSize,'Interpreter','latex')
        ylabel('$d/(\sqrt N b)$','FontSize',FontSize,'Interpreter','latex')

        title(['$2\times$ Chains, $N = $ ',num2str(N_Kuhn)],'FontSize',FontSize,'Interpreter','latex')

        pbaspect([1 1 1])
        set(gcf,'color','w')

        F1 = getframe(gcf);
        writeVideo(v1,F1);
        mov(nk_i) = F1;
    end
    N_Kuhns = fliplr(N_Kuhns);
    for nk_i = 1:length(N_Kuhns)
        figure(4); clf; hold on

        N_Kuhn = N_Kuhns(nk_i);
        dist_grid = lambda_grid*(sqrt(2*N_Kuhn)*b);

        sigma_a_theor = (3/2)^(-1/2)*sqrt(2*N_Kuhn)*b;
        k_encounter = ((3/2/pi)^(3/2))*D/(NormLength^2)*N_Kuhn^(-3/2);
        ka0_theor = k_encounter*exp(-ea_grid);
        ka_plot = ka0_theor.*exp(-(dist_grid.^2)/(sigma_a_theor^2));

        cmax = ((3/2/pi)^(3/2))*D/(NormLength^2)*(min(N_Kuhns))^(-3/2)...
            *exp(-min(min(ea_grid)))*exp(-(0.^2)/(((3/2)^(-1/2)*sqrt(min(N_Kuhns)*b)^2)));

        s = surf(ea_grid,lambda_grid,ka_plot*tau0);
        s.FaceColor = 'interp';

        c = colorbar;
        c.Ticks = [1e-3 2e-3 4e-3 8e-3];
        set(gca,'FontSize',FontSize/1.5)
        set(gca,'colorscale','log')
        clab = c.Label;
        clab.String = '$\bar{k}_a \tau_0$';
        clab.Interpreter = 'latex';
        clab.FontSize = FontSize;

        caxis([1e-3 cmax*tau0])
        %     set(gca,'Colorscale','log')
        set(gca,'xscale','log')

        xlabel('$\varepsilon_a/k_b T$','FontSize',FontSize,'Interpreter','latex')
        ylabel('$d/(\sqrt N b)$','FontSize',FontSize,'Interpreter','latex')

        title(['$2\times$ Chains, $N = $ ',num2str(N_Kuhn)],'FontSize',FontSize,'Interpreter','latex')

        pbaspect([1 1 1])
        set(gcf,'color','w')

        F1 = getframe(gcf);
        writeVideo(v1,F1);
        mov1(nk_i+length(N_Kuhns)) = F1;
    end
    close(v1)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotOutputsOverActivationEnergy(N_Kuhns,b,DiffCoeffs,ka_in,...
    A_bead,dA_bead,R2_bead,A_meso,dA_meso,R2_meso,R2_meso_vs_bead)

global LengthConversion FitPowerOrExponent FontSize...
    OutputFolder
 
NormLength = b*LengthConversion;
tau0 = (NormLength^2)/DiffCoeffs;
eaStar = -log(tau0*ka_in);

if FitPowerOrExponent==0
    AddOn = 'Power.';
else
    AddOn = 'Exp.';
end

%% Plot pre-factor
figure(5); clf; hold on
plot([1e-3 1e1],[1 1],'k--')
% styles = {'o','^','sq','v','d'};
% LegendEntries = [];
ColorRng = (linspace(0,1,length(N_Kuhns)))';
Colors = [zeros(size(ColorRng)) flipud(ColorRng) zeros(size(ColorRng))];

width = 0.125;
for i=1:size(A_bead,1)
    Color = Colors(i,:);
    if i==1
        dx1 = -2*width;
    elseif i==2
        dx1 = 0;
    elseif i==3
        dx1 = 2*width;
    end
   
    dx2 = -width/2;    
    dx = dx1+dx2;
    x_bead = [1 2 3]+dx;
    y_bead = fliplr(A_bead(i,:));
    dy_bead = fliplr(dA_bead(i,:));

    dx2 = +width/2;    
    dx = dx1+dx2;
    x_meso = [1 2 3]+dx;
    y_meso = fliplr(A_meso(i,:));
    dy_meso = fliplr(dA_meso(i,:));
        
    b_bead(i) = bar(x_bead,y_bead);
    b_bead(i).FaceColor = Color;
    b_bead(i).BarWidth = width;
    b_bead(i).FaceAlpha = 0.25;

    b_meso(i) = bar(x_meso,y_meso);
    b_meso(i).FaceColor = Color;
    b_meso(i).BarWidth = width;
    b_meso(i).FaceAlpha = 0.25;

    e_bead(i) = errorbar(x_bead,y_bead,dy_bead);
    e_bead(i).Marker = 'o';
    e_bead(i).MarkerEdgeColor = 'k';
    e_bead(i).MarkerFaceColor = Color;    
    e_bead(i).Color = 'k';
    e_bead(i).LineStyle = 'none';
    e_bead(i).MarkerSize = 4;

    e_meso(i) = errorbar(x_meso,y_meso,dy_meso);
    e_meso(i).Marker = '^';
    e_meso(i).MarkerEdgeColor = 'k';
    e_meso(i).MarkerFaceColor = Color;    
    e_meso(i).Color = 'k';
    e_meso(i).LineStyle = 'none';
    e_meso(i).MarkerSize = 4;
end

set(gca,'FontSize',FontSize/1.5)
% set(gca,'XScale','log')
set(gca,'YScale','default')
ylim([0 Inf])
% xlim([5e-3 2e0])
% xticks([1e-2 1e-1 1e0])
xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'10^{-2}','10^{-1}','10^0'})

% l = legend([s,e],LegendEntries);
l = legend([e_bead(1),e_meso(1)],{'Bead-spring','Mesoscale'});
l.FontSize = FontSize/1.75;
l.Location = 'Northwest';
l.Interpreter = 'latex';

xlabel('$\varepsilon_a^*$ ($k_b T$)','FontSize',FontSize,'Interpreter','latex')
ylabel('$A$','FontSize',FontSize,'Interpreter','latex')

set(gcf,'Color','w')
pbaspect([1 1.5 1])
set(gcf,'Position',[100 100 325 400])

FileTag = 'Attempt frequency pre-factor vs ea';
saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])

%% Plot goodness of fit between models
figure(6); clf; hold on
plot([1e-3 1e1],[1 1],'k--')
width = 0.125*2/3;
for i=1:size(A_bead,1)
    Color = Colors(i,:);
    if i==1
        dx1 = -3*width;
    elseif i==2
        dx1 = 0;
    elseif i==3
        dx1 = 3*width;
    end
   
    dx2 = -width;    
    dx = dx1+dx2;
    x_bead = [1 2 3]+dx;
    y_bead = fliplr(R2_bead(i,:));

    dx2 = 0;    
    dx = dx1+dx2;
    x_meso = [1 2 3]+dx;
    y_meso = fliplr(R2_meso(i,:));

    dx2 = +width;    
    dx = dx1+dx2;
    x_comp = [1 2 3]+dx;
    y_comp = fliplr(R2_meso_vs_bead(i,:));

    b_bead(i) = bar(x_bead,y_bead);
    b_bead(i).FaceColor = Color;
    b_bead(i).BarWidth = width;
    b_bead(i).FaceAlpha = 0.25;

    s_bead(i) = scatter(x_bead,y_bead);
    s_bead(i).Marker = 'o';
    s_bead(i).MarkerFaceColor = Color;
    s_bead(i).MarkerEdgeColor = 'k';
    s_bead(i).SizeData = 20;

    b_meso(i) = bar(x_meso,y_meso);
    b_meso(i).FaceColor = Color;
    b_meso(i).BarWidth = width;
    b_meso(i).FaceAlpha = 0.25;

    s_meso(i) = scatter(x_meso,y_meso);
    s_meso(i).Marker = '^';
    s_meso(i).MarkerFaceColor = Color;
    s_meso(i).MarkerEdgeColor = 'k';
    s_meso(i).SizeData = 20;

    b_comp(i) = bar(x_comp,y_comp);
    b_comp(i).FaceColor = Color;
    b_comp(i).BarWidth = width;
    b_comp(i).FaceAlpha = 0.25;

    s_comp(i) = scatter(x_comp,y_comp);
    s_comp(i).Marker = 'sq';
    s_comp(i).MarkerFaceColor = Color;
    s_comp(i).MarkerEdgeColor = 'k';
    s_comp(i).SizeData = 20;
end

set(gca,'FontSize',FontSize/1.5)
set(gca,'YScale','default')
ylim([0.5 1.5])

xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'10^{-2}','10^{-1}','10^0'})

% l = legend([s,e],LegendEntries);
l = legend([s_bead(1),s_meso(1),s_comp(1)],{'Scaling vs bead-spring','Scaling vs. mesoscale','Mesoscale vs. bead-spring'});
l.FontSize = FontSize/1.75;
l.Location = 'Northwest';
l.Interpreter = 'latex';

xlabel('$\varepsilon_a^*$ ($k_b T$)','FontSize',FontSize,'Interpreter','latex')
ylabel('$R^2$','FontSize',FontSize,'Interpreter','latex')

set(gcf,'Color','w')
pbaspect([2 1 1])
set(gcf,'Position',[100 100 275 400])

FileTag = 'R2 between models - including meso vs bead-spring';
saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])

%% Plot goodness of fit between models
figure(7); clf; hold on
plot([1e-3 1e1],[1 1],'k--')
width = 0.125;
for i=1:size(A_bead,1)
    Color = Colors(i,:);
    if i==1
        dx1 = -2*width;
    elseif i==2
        dx1 = 0;
    elseif i==3
        dx1 = 2*width;
    end
   
    dx2 = -width/2;    
    dx = dx1+dx2;
    x_bead = [1 2 3]+dx;
    y_bead = fliplr(R2_bead(i,:));

    dx2 = +width/2;    
    dx = dx1+dx2;
    x_meso = [1 2 3]+dx;
    y_meso = fliplr(R2_meso(i,:));

    b_bead(i) = bar(x_bead,y_bead);
    b_bead(i).FaceColor = Color;
    b_bead(i).BarWidth = width;
    b_bead(i).FaceAlpha = 0.25;

    s_bead(i) = scatter(x_bead,y_bead);
    s_bead(i).Marker = 'o';
    s_bead(i).MarkerFaceColor = Color;
    s_bead(i).MarkerEdgeColor = 'k';
    s_bead(i).SizeData = 20;

    b_meso(i) = bar(x_meso,y_meso);
    b_meso(i).FaceColor = Color;
    b_meso(i).BarWidth = width;
    b_meso(i).FaceAlpha = 0.25;

    s_meso(i) = scatter(x_meso,y_meso);
    s_meso(i).Marker = '^';
    s_meso(i).MarkerFaceColor = Color;
    s_meso(i).MarkerEdgeColor = 'k';
    s_meso(i).SizeData = 20;
end

set(gca,'FontSize',FontSize/1.5)
set(gca,'YScale','default')
ylim([0.5 1.5])

xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'10^{-2}','10^{-1}','10^0'})

% l = legend([s,e],LegendEntries);
l = legend([s_bead(1),s_meso(1)],{'Bead-spring','Mesoscale'});
l.FontSize = FontSize/1.75;
l.Location = 'Northwest';
l.Interpreter = 'latex';

xlabel('$\varepsilon_a^*$ ($k_b T$)','FontSize',FontSize,'Interpreter','latex')
ylabel('$R^2$','FontSize',FontSize,'Interpreter','latex')

set(gcf,'Color','w')
pbaspect([2 1 1])
set(gcf,'Position',[100 100 275 400])

FileTag = 'R2 between models';
saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A_bead,dA_bead,A_meso,dA_meso,...
    R2_bead,R2_meso,R2_meso_vs_bead] = ...
    PlotEnsembleAvgDynamics(Override,Package,N_Kuhns,b,ka_in,kd_in,D)

global EnsembleDynamicsFileName FontSize OutputFolder...
    LengthConversion PlotScalingTheory FittedDynamicsFileName...
    PlotOnlyScalingTheory timenorm NormalizeRates AppendScalingData

NormLength = b*LengthConversion;
tau0 = (NormLength^2)/D;
eaStar = -log(tau0*ka_in);

PlotwrtStretchOrNormalizedDist = 0;     %0 to plot wrt stretch
if PlotwrtStretchOrNormalizedDist==1
    denoms = N_Kuhns*b;
    AddOn = ['Ensemble.ea.',num2str(eaStar,'%.1e'),'.'];
else
    denoms = sqrt(N_Kuhns)*b;
    AddOn = ['Ensemble.ea.',num2str(eaStar,'%.1e'),'.stretch.'];
end

if ~isfile(FittedDynamicsFileName) || Override

    %Import the data
    EnsembleDynamics = load(EnsembleDynamicsFileName,'-mat');

    ka_bead = EnsembleDynamics.ka_bead;
    kaSE_bead = EnsembleDynamics.kaSE_bead;
    ka1_bead = EnsembleDynamics.ka1_bead;
    ka1SE_bead = EnsembleDynamics.ka1SE_bead;
    kd_bead = EnsembleDynamics.kd_bead;
    kdSE_bead = EnsembleDynamics.kdSE_bead;

    fa_bead = EnsembleDynamics.fa_bead;
    faSE_bead = EnsembleDynamics.faSE_bead;
    fd_bead = EnsembleDynamics.fd_bead;
    fdSE_bead = EnsembleDynamics.fdSE_bead;

    ka_meso = EnsembleDynamics.ka_meso;
    kaSE_meso = EnsembleDynamics.kaSE_meso;
    ka1_meso = EnsembleDynamics.ka1_meso;
    ka1SE_meso = EnsembleDynamics.ka1SE_meso;
    kd_meso = EnsembleDynamics.kd_meso;
    kdSE_meso = EnsembleDynamics.kdSE_meso;

    fa_meso = EnsembleDynamics.fa_meso;
    faSE_meso = EnsembleDynamics.faSE_meso;
    fd_meso = EnsembleDynamics.fd_meso;
    fdSE_meso = EnsembleDynamics.fdSE_meso;

    if AppendScalingData
        ka_thry = EnsembleDynamics.ka_thry;
        kaSE_thry = EnsembleDynamics.kaSE_thry;
        ka1_thry = EnsembleDynamics.ka1_thry;
        ka1SE_thry = EnsembleDynamics.ka1SE_thry;
        kd_thry = EnsembleDynamics.kd_thry;
        kdSE_thry = EnsembleDynamics.kdSE_thry;

        fa_thry = EnsembleDynamics.fa_thry;
        faSE_thry = EnsembleDynamics.faSE_thry;
        fd_thry = EnsembleDynamics.fd_thry;
        fdSE_thry = EnsembleDynamics.fdSE_thry;
    end

    % For plotting outputs wrt dist
    figure(1); clf; hold on %ka
    figure(2); clf; hold on %kd
    figure(3); clf; hold on %ka1 vs ka

    R2_bead = zeros(size(N_Kuhns));
    R2_meso = zeros(size(N_Kuhns));
    A_bead = zeros(size(N_Kuhns));
    A_meso = zeros(size(N_Kuhns));
    dA_bead = zeros(size(N_Kuhns));
    dA_meso = zeros(size(N_Kuhns));
    R2_meso_vs_bead = zeros(size(N_Kuhns));

    R21_bead = zeros(size(N_Kuhns));
    R21_meso = zeros(size(N_Kuhns));
    R21_thry = zeros(size(N_Kuhns));
    A1_bead = zeros(size(N_Kuhns));
    A1_meso = zeros(size(N_Kuhns));
    A1_thry = zeros(size(N_Kuhns));
    dA1_bead = zeros(size(N_Kuhns));
    dA1_meso = zeros(size(N_Kuhns));
    dA1_thry = zeros(size(N_Kuhns));

    ForPaper = 1; %Set 1 to plot versions of figures in paper
    if ForPaper==1
        PlotCases = [12 18 36];
    else
        PlotCases = N_Kuhns;
    end

    % Loop over N_Kuhn
    plotct = 0;
    for nk = 1:length(N_Kuhns)
        N_Kuhn = N_Kuhns(nk);
        PackageTemp = Package(Package(:,12)==N_Kuhn,:);
        Separations = unique(PackageTemp(:,14));

        if NormalizeRates==1
            ka_meso_temp = ka_meso(:,nk);
            ka_bead_temp = ka_bead(:,nk);
            kaSE_meso_temp = kaSE_meso(:,nk);
            kaSE_bead_temp = kaSE_bead(:,nk);

            ka1_meso_temp = ka1_meso(:,nk);
            ka1_bead_temp = ka1_bead(:,nk);
            ka1SE_meso_temp = ka1SE_meso(:,nk);
            ka1SE_bead_temp = ka1SE_bead(:,nk);
            
            kd_meso_temp = kd_meso(:,nk);
            kd_bead_temp = kd_bead(:,nk);
            kdSE_meso_temp = kdSE_meso(:,nk);
            kdSE_bead_temp = kdSE_bead(:,nk);
            if AppendScalingData
                ka_thry_temp = ka_thry(:,nk);
                kaSE_thry_temp = kaSE_thry(:,nk);

                ka1_thry_temp = ka1_thry(:,nk);
                ka1SE_thry_temp = ka1SE_thry(:,nk);
                
                kd_thry_temp = kd_thry(:,nk);
                kdSE_thry_temp = kdSE_thry(:,nk);
            end
        else
            % Isolte data for N_Kuhn
            ka_meso_temp = ka_meso(:,nk)*1e3;
            ka_bead_temp = ka_bead(:,nk)*1e3;
            kaSE_meso_temp = kaSE_meso(:,nk)*1e3;
            kaSE_bead_temp = kaSE_bead(:,nk)*1e3;

            ka1_meso_temp = ka1_meso(:,nk)*1e3;
            ka1_bead_temp = ka1_bead(:,nk)*1e3;
            ka1SE_meso_temp = ka1SE_meso(:,nk)*1e3;
            ka1SE_bead_temp = ka1SE_bead(:,nk)*1e3;
            
            kd_meso_temp = kd_meso(:,nk)*1e3;
            kd_bead_temp = kd_bead(:,nk)*1e3;
            kdSE_meso_temp = kdSE_meso(:,nk)*1e3;
            kdSE_bead_temp = kdSE_bead(:,nk)*1e3;
            
            if AppendScalingData
                ka_thry_temp = ka_thry(:,nk)*1e3;
                kaSE_thry_temp = kaSE_thry(:,nk)*1e3;

                ka1_thry_temp = ka1_thry(:,nk)*1e3;
                ka1SE_thry_temp = ka1SE_thry(:,nk)*1e3;

                kd_thry_temp = kd_thry(:,nk)*1e3;
                kdSE_thry_temp = kdSE_thry(:,nk)*1e3;
            end
        end

        fa_bead_temp = fa_bead(:,nk);
        faSE_bead_temp = faSE_bead(:,nk);
        fa_meso_temp = fa_meso(:,nk);
        faSE_meso_temp = faSE_meso(:,nk);

        fd_bead_temp = fd_bead(:,nk);
        fdSE_bead_temp = fdSE_bead(:,nk);
        fd_meso_temp = fd_meso(:,nk);
        fdSE_meso_temp = fdSE_meso(:,nk);

        if AppendScalingData
            fa_thry_temp = fa_thry(:,nk);
            faSE_thry_temp = faSE_thry(:,nk);
            fd_thry_temp = fd_thry(:,nk);
            fdSE_thry_temp = fdSE_thry(:,nk);
        end

        ColorRange = (linspace(0,1,length(PlotCases)))';
        ColorFactn = [0 0.85 0];
        Colorsn = [flipud(ColorRange) flipud(ColorRange) flipud(ColorRange)].*ColorFactn;

        ColorRange2 = (linspace(0.4,1,length(PlotCases)))';
        ColorFactka = [0 0.85 0];
        Colorsnka = [flipud(ColorRange2) flipud(ColorRange2) flipud(ColorRange2)].*ColorFactka;
        ColorFactka1 = [1 0.8 0];
        % ColorFactka1 = [0 0.5 0.85];
        Colorsnka1 = [flipud(ColorRange2) flipud(ColorRange2) flipud(ColorRange2)].*ColorFactka1;

        % Regression analysis between models
        RSS = sum((ka_bead_temp-ka_meso_temp).^2);
        TSS = sum((ka_meso_temp-mean(ka_meso_temp)).^2);
        R2_meso_vs_bead(nk) = 1-RSS/TSS;

        % Regression analysis between models for rate to first attachment
        NormLength = b*LengthConversion;
        if NormalizeRates==1
            eaOverkbT = -log(ka_in*NormLength^2/D);
            D_in = tau0*D;
        else
            eaOverkbT = -log(ka_in*NormLength^2/D);
            D_in = D;
        end
%         sigma_a_theor = (3/2)^(-1/2)*sqrt(2*N_Kuhn)*b;

        A_nom = 0.9;
        [A_bead(nk),dA_bead(nk),R2_bead(nk)] =...
            FitThePrefactor(Separations,ka_bead_temp,A_nom,...
            N_Kuhn,b,eaOverkbT);
        [A_meso(nk),dA_meso(nk),R2_meso(nk)]...
            = FitThePrefactor(Separations,ka_meso_temp,A_nom,...
            N_Kuhn,b,eaOverkbT);

        [A1_bead(nk),dA1_bead(nk),R21_bead(nk)] =...
            FitThePrefactor(Separations,ka1_bead_temp,A_nom,...
            N_Kuhn,b,eaOverkbT);
        [A1_meso(nk),dA1_meso(nk),R21_meso(nk)]...
            = FitThePrefactor(Separations,ka1_meso_temp,A_nom,...
            N_Kuhn,b,eaOverkbT);
        [A1_thry(nk),dA1_thry(nk),R21_thry(nk)]...
            = FitThePrefactor(Separations,ka1_thry_temp,A_nom,...
            N_Kuhn,b,eaOverkbT);

        if ismember(N_Kuhn,PlotCases)
            plotct = plotct+1;
            LegendEntries{plotct} = ['$N$ = ',num2str(N_Kuhn)];

            %% Plot attachment rate
            figure(1)
            xplot = (linspace(Separations(1),Separations(end),100))';


            if PlotScalingTheory
                P_encounter_bead = A_bead(nk)*((3/2/pi/N_Kuhn)^(3/2))*D_in/(NormLength^2);
                ka0_theor_bead = P_encounter_bead*exp(-eaOverkbT);
                lambda = xplot/sqrt(N_Kuhn)/b;
                yplot_bead = ka0_theor_bead*exp(-3/2*lambda.^2);

                Ptheor(plotct) = plot(xplot/denoms(nk),yplot_bead);
                Ptheor(plotct).LineStyle = '-';
                Ptheor(plotct).Color = Colorsn(plotct,:);
                Ptheor(plotct).LineWidth = 1;

                P_encounter_meso = A_meso(nk)*((3/2/pi/N_Kuhn)^(3/2))*D_in/(NormLength^2);
                ka0_theor_meso = P_encounter_meso*exp(-eaOverkbT);
                yplot_meso = ka0_theor_meso*exp(-3/2*lambda.^2);


                Ptheor_meso(plotct) = plot(xplot/denoms(nk),yplot_meso);
                Ptheor_meso(plotct).LineStyle = '-.';
                Ptheor_meso(plotct).Color = Colorsn(plotct,:);
                Ptheor_meso(plotct).LineWidth = 1;
            end

            e =  errorbar(Separations/denoms(nk),ka_bead_temp,kaSE_bead_temp);
            e.LineStyle = 'none';
            e.Color = 'k';
            e.Marker = 'o';
            if ~PlotOnlyScalingTheory
                e.MarkerFaceColor = 'k';
            else
                e.MarkerFaceColor = Colorsn(plotct,:);
            end

            e =  errorbar(Separations/denoms(nk),ka_meso_temp,kaSE_meso_temp);
            e.LineStyle = 'none';
            e.Color  = 'k';
            e.Marker = '^';
            if ~PlotOnlyScalingTheory
                e.MarkerFaceColor = 'none';
            else
                e.MarkerFaceColor = Colorsn(plotct,:);
            end

            if AppendScalingData
                e =  errorbar(Separations/denoms(nk),ka_thry_temp,kaSE_thry_temp);
                e.LineStyle = 'none';
                e.Color = 'k';
                e.Marker = 'd';
                if ~PlotOnlyScalingTheory
                    e.MarkerFaceColor = 'k';
                else
                    e.MarkerFaceColor = Colorsn(plotct,:);
                end
            end

            %% Plot virgin vs overall attachment rates
            figure(3)
            subplot(1,3,1); hold on
            if PlotScalingTheory
                P_encounter_bead = A_bead(nk)*((3/2/pi/N_Kuhn)^(3/2))*D_in/(NormLength^2);
                ka0_theor_bead = P_encounter_bead*exp(-eaOverkbT);
                lambda = xplot/sqrt(N_Kuhn)/b;
                yplot_bead = ka0_theor_bead*exp(-3/2*lambda.^2);

                p1 = plot(xplot/denoms(nk),yplot_bead);
                p1.LineStyle = '-';
                p1.Color = Colorsnka(plotct,:);
                p1.LineWidth = 1;

                P_encounter_bead = A1_bead(nk)*((3/2/pi/N_Kuhn)^(3/2))*D_in/(NormLength^2);
                ka0_theor_bead = P_encounter_bead*exp(-eaOverkbT);
                lambda = xplot/sqrt(N_Kuhn)/b;
                yplot_bead = ka0_theor_bead*exp(-3/2*lambda.^2);

                p1 = plot(xplot/denoms(nk),yplot_bead);
                p1.LineStyle = '-';
                p1.Color = Colorsnka1(plotct,:);
                p1.LineWidth = 1;
            end

            e =  errorbar(Separations/denoms(nk),ka_bead_temp,kaSE_bead_temp);
            e.LineStyle = 'none';
            e.Color = 'k';
            e.Marker = 'o';
            e.MarkerFaceColor = Colorsnka(plotct,:);

            e =  errorbar(Separations/denoms(nk),ka1_bead_temp,ka1SE_bead_temp);
            e.LineStyle = 'none';
            e.Color = 'k';
            e.Marker = 'o';
            e.MarkerFaceColor = Colorsnka1(plotct,:);
            
            subplot(1,3,2); hold on
            if PlotScalingTheory
                P_encounter_meso = A_meso(nk)*((3/2/pi/N_Kuhn)^(3/2))*D_in/(NormLength^2);
                ka0_theor_meso = P_encounter_meso*exp(-eaOverkbT);
                yplot_meso = ka0_theor_meso*exp(-3/2*lambda.^2);

                p1 = plot(xplot/denoms(nk),yplot_meso);
                p1.LineStyle = '-';
                p1.Color = Colorsnka(plotct,:);
                p1.LineWidth = 1;

                P_encounter_meso = A1_meso(nk)*((3/2/pi/N_Kuhn)^(3/2))*D_in/(NormLength^2);
                ka0_theor_meso = P_encounter_meso*exp(-eaOverkbT);
                yplot_meso = ka0_theor_meso*exp(-3/2*lambda.^2);

                p1 = plot(xplot/denoms(nk),yplot_meso);
                p1.LineStyle = '-';
                p1.Color = Colorsnka1(plotct,:);
                p1.LineWidth = 1;
            end

            e =  errorbar(Separations/denoms(nk),ka_meso_temp,kaSE_meso_temp);
            e.LineStyle = 'none';
            e.Color  = 'k';
            e.Marker = '^';
            e.MarkerFaceColor = Colorsnka(plotct,:);

            e =  errorbar(Separations/denoms(nk),ka1_meso_temp,ka1SE_meso_temp);
            e.LineStyle = 'none';
            e.Color  = 'k';
            e.Marker = '^';
            e.MarkerFaceColor = Colorsnka1(plotct,:);

            if AppendScalingData
                subplot(1,3,3); hold on
                if PlotScalingTheory
                    P_encounter_bead = A_bead(nk)*((3/2/pi/N_Kuhn)^(3/2))*D_in/(NormLength^2);
                    ka0_theor_bead = P_encounter_bead*exp(-eaOverkbT);
                    lambda = xplot/sqrt(N_Kuhn)/b;
                    yplot_bead = ka0_theor_bead*exp(-3/2*lambda.^2);

                    p1 = plot(xplot/denoms(nk),yplot_bead);
                    p1.LineStyle = '-';
                    p1.Color = Colorsnka(plotct,:);
                    p1.LineWidth = 1;


                    P_encounter_bead = A1_bead(nk)*((3/2/pi/N_Kuhn)^(3/2))*D_in/(NormLength^2);
                    ka0_theor_bead = P_encounter_bead*exp(-eaOverkbT);
                    lambda = xplot/sqrt(N_Kuhn)/b;
                    yplot_bead = ka0_theor_bead*exp(-3/2*lambda.^2);

                    p1 = plot(xplot/denoms(nk),yplot_bead);
                    p1.LineStyle = '-';
                    p1.Color = Colorsnka1(plotct,:);
                    p1.LineWidth = 1;
                end

                e =  errorbar(Separations/denoms(nk),ka_thry_temp,kaSE_thry_temp);
                e.LineStyle = 'none';
                e.Color = 'k';
                e.Marker = 'd';
                e.MarkerFaceColor = Colorsnka(plotct,:);

                e =  errorbar(Separations/denoms(nk),ka1_thry_temp,ka1SE_thry_temp);
                e.LineStyle = 'none';
                e.Color = 'k';
                e.Marker = 'd';
                e.MarkerFaceColor = Colorsnka1(plotct,:);
            end

            %% Plot detachment rate wrt. sep
            figure(2)
            plot(Separations/denoms(nk),kd_in/timenorm*ones(size(Separations)),'k--')

            e =  errorbar(Separations/denoms(nk),kd_bead_temp,kdSE_bead_temp);
            e.LineStyle = 'none';
            e.Color = 'k';
            e.Marker = 'o';
            if ~PlotOnlyScalingTheory
                e.MarkerFaceColor = 'none';
            else
                e.MarkerFaceColor = Colorsn(plotct,:);
            end

            e =  errorbar(Separations/denoms(nk),kd_meso_temp,kdSE_meso_temp');
            e.LineStyle = 'none';
            e.Color = 'k';
            e.Marker = 'd';
            if ~PlotOnlyScalingTheory
                e.MarkerFaceColor = 'none';
            else
                e.MarkerFaceColor = Colorsn(plotct,:);
            end           
        end

        if nk==1
            figure(10); clf; hold on

            e =  errorbar(Separations/denoms(nk),ka_bead_temp,kaSE_bead_temp);
            e.LineStyle = 'none';
            e.Color = 'k';
            e.Marker = 'o';
            e.MarkerFaceColor = 'k';

            e =  errorbar(Separations/denoms(nk),ka_meso_temp,kaSE_meso_temp);
            e.LineStyle = 'none';
            e.Color  = 'k';
            e.Marker = '^';
            e.MarkerFaceColor = 'none';

            if AppendScalingData
                e =  errorbar(Separations/denoms(nk),ka_thry_temp,kaSE_thry_temp);
                e.LineStyle = 'none';
                e.Color  = 'k';
                e.Marker = 'd';
                e.MarkerFaceColor = 'k';
            end

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([1 1 1])

            l = legend('Bead-spring','Mesoscale','Scaling Thry.');
            l.FontSize = FontSize/2;
            l.Interpreter = 'latex';

            if PlotwrtStretchOrNormalizedDist==0
                xlabel('$\lambda$','FontSize',FontSize,'Interpreter','latex')
                xlim([0 0.6])
            else
                xlabel('$d^*$ ($Nb$)','FontSize',FontSize,'Interpreter','latex')
            end
            if NormalizeRates==1
                ylabel('$\bar{k}_a^*$ ($\tau_0^{-1}$)','FontSize',FontSize,'Interpreter','latex')
                ylim([0 7.5e-3])
            else
                ylabel('$\bar{k}_a$ (kHz)','FontSize',FontSize,'Interpreter','latex')
                ylim([0 1.25e3])
            end


            FileTag = ['Nk.',num2str(N_Kuhn),'.ka vs d'];
            saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
            saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])

            figure(20); clf; hold on

            e =  errorbar(Separations/denoms(nk),kd_bead_temp,kdSE_bead_temp);
            e.LineStyle = 'none';
            e.Color = 'k';
            e.Marker = 'o';
            e.MarkerFaceColor = 'k';

            e =  errorbar(Separations/denoms(nk),kd_meso_temp,kdSE_meso_temp);
            e.LineStyle = 'none';
            e.Color  = 'k';
            e.Marker = '^';
            e.MarkerFaceColor = 'none';

            if AppendScalingData
                e =  errorbar(Separations/denoms(nk),kd_thry_temp,kdSE_thry_temp);
                e.LineStyle = 'none';
                e.Color  = 'k';
                e.Marker = 'd';
                e.MarkerFaceColor = 'none';
            end

            plot(Separations/denoms(nk),kd_in/timenorm*ones(size(Separations)),'r-')

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([1 1 1])

            l = legend('Bead-spring','Mesoscale','\textit{A priori}');
            l.FontSize = FontSize/2;
            l.Interpreter = 'latex';
            l.Location = 'Northwest';


            if PlotwrtStretchOrNormalizedDist==0
                xlabel('$\lambda$','FontSize',FontSize,'Interpreter','latex')
                xlim([0 0.6])
            else
                xlabel('$d^*$ ($Nb$)','FontSize',FontSize,'Interpreter','latex')
            end
            if NormalizeRates==1
                ylabel('$\bar{k}_d^*$ ($\tau_0^{-1}$)','FontSize',FontSize,'Interpreter','latex')
                ylim([0 7.5e-3])
                yticks([0 7.5e-3])
            else
                ylabel('$\bar{k}_d$ (kHz)','FontSize',FontSize,'Interpreter','latex')
                ylim([0 1.25e3])
            end

            if PlotOnlyScalingTheory
                if PlotwrtStretchOrNormalizedDist==1
                    set(gcf,'Position',[100 100 225 225])
                else
                    set(gcf,'Position',[100 100 200 200])
                    xticks([0 7.5e-3])
                end
            end
            FileTag = ['Nk.',num2str(N_Kuhn),'.kd vs d'];
            saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
            saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])
        end

        % Plot attached and detached fractions wrt. sep
        figure(5); clf; hold on
        e1 =  errorbar(Separations/denoms(nk),fa_bead_temp,faSE_bead_temp);
        e1.LineStyle = 'none';
        e1.Color = 'k';
        e1.Marker = 'o';
        e1.MarkerFaceColor = 'c';
        e1.MarkerSize = 8;

        e2 =  errorbar(Separations/denoms(nk),fd_bead_temp,fdSE_bead_temp);
        e2.LineStyle = 'none';
        e2.Color = 'k';
        e2.Marker = 'o';
        e2.MarkerFaceColor = 'r';
        e2.MarkerSize = 8;

        if AppendScalingData
            e3 =  errorbar(Separations/denoms(nk),fd_thry_temp,fdSE_thry_temp);
            e3.LineStyle = 'none';
            e3.Color = 'r';
            e3.Marker = 'd';
            e3.MarkerFaceColor = 'r';
            e3.MarkerSize = 12;
        end

        e1 =  errorbar(Separations/denoms(nk),fa_meso_temp,faSE_meso_temp);
        e1.LineStyle = 'none';
        e1.Color = 'c';
        e1.Marker = '^';
        e1.MarkerFaceColor = 'none';
        e1.MarkerSize = 8;

        e2 =  errorbar(Separations/denoms(nk),fd_meso_temp,fdSE_meso_temp);
        e2.LineStyle = 'none';
        e2.Color = 'r';
        e2.Marker = '^';
        e2.MarkerFaceColor = 'none';
        e2.MarkerSize = 8;

        if AppendScalingData
            e3 =  errorbar(Separations/denoms(nk),fa_thry_temp,faSE_thry_temp);
            e3.LineStyle = 'none';
            e3.Color = 'c';
            e3.Marker = 'd';
            e3.MarkerFaceColor = 'c';
            e3.MarkerSize = 12;
        end

        set(gcf,'Position',[200 200 210 210])

        set(gca,'FontSize',FontSize/1.5)
        set(gcf,'Color','w')
        pbaspect([1 1 1])

        if PlotwrtStretchOrNormalizedDist==0
            xlabel('$\lambda$','FontSize',FontSize,'Interpreter','latex')
            xlim([0 0.6])
        else
            xlabel('$d^*$ ($Nb$)','FontSize',FontSize,'Interpreter','latex')
        end
        ylabel('$f$','FontSize',FontSize,'Interpreter','latex')

        FileTag = ['Nk.',num2str(N_Kuhn),'.f vs d'];
        saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
        saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])
    end

    %% Polish emergent ka
    figure(1)
    set(gca,'FontSize',FontSize/1.5)
    set(gcf,'Color','w')
    pbaspect([1 1 1])

    if PlotOnlyScalingTheory
        l = legend(Ptheor,LegendEntries);
    else
        l = legend(P,LegendEntries);
    end
    l.FontSize = FontSize/2;
    l.Interpreter = 'latex';


    if PlotwrtStretchOrNormalizedDist==0
        xlabel('$\lambda$','FontSize',FontSize,'Interpreter','latex')
        xlim([0 1.5])
    else
        xlabel('$d^*$ ($Nb$)','FontSize',FontSize,'Interpreter','latex')
    end
    if NormalizeRates==1
        ylabel('$\bar{k}_a^*$ ($\tau_0^{-1}$)','FontSize',FontSize,'Interpreter','latex')
        ylim([0 6e-3])
    else
        ylabel('$\bar{k}_a$ (kHz)','FontSize',FontSize,'Interpreter','latex')
        ylim([0 1.25e3])
    end

    FileTag = 'ka ';
    saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
    saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])



    %% Polish emergent ka1 vs ka
    % if eaStar==1
    %     ylims = [0 4e-3];
    % elseif (eaStar-0.1)<0.001
    %     ylims = [0 5e-3];
    % else
    %     ylims = [0 Inf];
    % end
    ylims = [0 6e-3];
    figure(3)
    set(gcf,'Color','w')
    set(gcf,'Position',[1000 100 1080 400])
    subplot(1,3,1)
    set(gca,'FontSize',FontSize/1.5)
    if eaStar<0.099
        xlabel('$\lambda$','FontSize',FontSize,'Interpreter','latex')
    else
        xticks([0 0.5 1 1.5])
        xticklabels([])
    end
    % ylabel('$\bar{k}_a \tau_0$ (Bead-spring)','FontSize',FontSize,'Interpreter','latex')
    ylabel('$\bar{k}_a \tau_0$','FontSize',FontSize,'Interpreter','latex')
    pbaspect([1 1 1])
    xlim([0 1.5])
    ylim(ylims)

    subplot(1,3,2)
    set(gca,'FontSize',FontSize/1.5)
    if eaStar<0.099
        xlabel('$\lambda$','FontSize',FontSize,'Interpreter','latex')
    else
        xticks([0 0.5 1 1.5])
        xticklabels([])
    end
    % ylabel('$\bar{k}_a \tau_0$ (Mesoscale)','FontSize',FontSize,'Interpreter','latex')
    % ylabel('$\bar{k}_a \tau_0$','FontSize',FontSize,'Interpreter','latex')
    yticks([0 2e-3 4e-3 6e-3])
    yticklabels([])
    pbaspect([1 1 1])
    xlim([0 1.5])
    ylim(ylims)

    subplot(1,3,3)
    set(gca,'FontSize',FontSize/1.5)
    if eaStar<0.099
        xlabel('$\lambda$','FontSize',FontSize,'Interpreter','latex')
    else
        xticks([0 0.5 1 1.5])
        xticklabels([])
    end
    % ylabel('$\bar{k}_a \tau_0$ (Scaling Thry.)','FontSize',FontSize,'Interpreter','latex')
    % ylabel('$\bar{k}_a \tau_0$','FontSize',FontSize,'Interpreter','latex')
    yticks([0 2e-3 4e-3 6e-3])
    yticklabels([])
    pbaspect([1 1 1])
    xlim([0 1.5])
    ylim(ylims)

    FileTag = 'ka vs ka1';
    saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
    saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])



    %% Polish emergent kd
    figure(2)
    set(gca,'FontSize',FontSize/1.5)
    set(gcf,'Color','w')
    pbaspect([1 1 1])

    if PlotwrtStretchOrNormalizedDist==0
        xlabel('$\lambda$','FontSize',FontSize,'Interpreter','latex')
        xlim([0 0.6])
    else
        xlabel('$d^*$ ($Nb$)','FontSize',FontSize,'Interpreter','latex')
    end
    if NormalizeRates==1
        ylabel('$\bar{k}_d^*$ ($\tau_0^{-1}$)','FontSize',FontSize,'Interpreter','latex')
        ylim([0 7.5e-3])
    else
        ylabel('$\bar{k}_d$ (kHz)','FontSize',FontSize,'Interpreter','latex')
        ylim([0 1.25e3])
    end

    if PlotOnlyScalingTheory
        set(gcf,'Position',[100 100 250 250])
    end
    FileTag = 'kd ';
    saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
    saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])

    %% Plot the fitted pre-factor
    figure(9); clf; hold on
    plot([0 sqrt(max(N_Kuhns))],[1 1],'k--')
    x = sqrt(N_Kuhns);
    y = A_meso;         % unitless
    s = scatter(x,y);
    s.MarkerEdgeColor = 'k';
    s.MarkerFaceColor = 'k';
    set(gca,'FontSize',FontSize/1.5)

    set(gcf,'Color','w')
    pbaspect([2 1 1])
    xlim([3 6])
    ylim([0 1.5])

    set(gcf,'Position',[200 200 260 260])

    xlabel('$\sqrt N$','FontSize',FontSize,'Interpreter','latex')
    ylabel('$A$','FontSize',FontSize,'Interpreter','latex')

    FileTag = 'Pre-factor vs sqrt N b';
    saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
    saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])

    FittedDynamics.A_bead = A_bead;
    FittedDynamics.dA_bead = dA_bead;
    FittedDynamics.R2_bead = R2_bead;
    FittedDynamics.A_meso = A_meso;
    FittedDynamics.dA_meso = dA_meso;
    FittedDynamics.R2_meso = R2_meso;
    FittedDynamics.R2_meso_vs_bead = R2_meso_vs_bead;

    save(FittedDynamicsFileName,'-struct','FittedDynamics');
else
    FittedDynamics = load(FittedDynamicsFileName,'-mat');

    A_bead = FittedDynamics.A_bead;
    dA_bead = FittedDynamics.dA_bead;
    R2_bead = FittedDynamics.R2_bead;
    A_meso = FittedDynamics.A_meso;
    dA_meso = FittedDynamics.dA_meso;
    R2_meso = FittedDynamics.R2_meso;
    R2_meso_vs_bead = FittedDynamics.R2_meso_vs_bead;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R2a1,R2a,R2d,...
    k_a1_meso_fit,k_a_meso_fit,k_d_meso_fit,...
    k_a1_bead_fit,k_a_bead_fit,k_d_bead_fit]...
    = PlotBondLifetimes(Override,Package,N_Kuhns,b,DiffCoeffs,ka)

global AllDynamicsFileName FontSize OutputFolder LengthConversion...
    R2DynamicsFileName NormalizeRates AppendScalingData

NormLength = b*LengthConversion;
tau0 = (NormLength^2)/DiffCoeffs;
eaStar = -log(tau0*ka);

AddOn = [];
if ~isfile(R2DynamicsFileName) || Override
    AllDynamics = load(AllDynamicsFileName,'-mat');

    tau_d1_bead = AllDynamics.tau_d1_bead;
    tau_a_bead = AllDynamics.tau_a_bead;
    tau_d_bead = AllDynamics.tau_d_bead;

    tau_d1_meso = AllDynamics.tau_d1_meso;
    tau_a_meso = AllDynamics.tau_a_meso;
    tau_d_meso = AllDynamics.tau_d_meso;

    tau_d1_thry = AllDynamics.tau_d1_thry;
    tau_a_thry = AllDynamics.tau_a_thry;
    tau_d_thry = AllDynamics.tau_d_thry;

    R2a1 = zeros(length(N_Kuhns),1);
    R2a = zeros(length(N_Kuhns),1);
    R2d = zeros(length(N_Kuhns),1);

    k_a1_meso_fit = zeros(length(N_Kuhns),1);
    k_a_meso_fit = zeros(length(N_Kuhns),1);
    k_d_meso_fit = zeros(length(N_Kuhns),1);
    k_a1_bead_fit = zeros(length(N_Kuhns),1);
    k_a_bead_fit = zeros(length(N_Kuhns),1);
    k_d_bead_fit = zeros(length(N_Kuhns),1);
    k_a1_thry_fit = zeros(length(N_Kuhns),1);
    k_a_thry_fit = zeros(length(N_Kuhns),1);
    k_d_thry_fit = zeros(length(N_Kuhns),1);

    se_k_a1_meso_fit = zeros(length(N_Kuhns),1);
    se_k_a_meso_fit = zeros(length(N_Kuhns),1);
    se_k_d_meso_fit = zeros(length(N_Kuhns),1);
    se_k_a1_bead_fit = zeros(length(N_Kuhns),1);
    se_k_a_bead_fit = zeros(length(N_Kuhns),1);
    se_k_d_bead_fit = zeros(length(N_Kuhns),1);
    se_k_a1_thry_fit = zeros(length(N_Kuhns),1);
    se_k_a_thry_fit = zeros(length(N_Kuhns),1);
    se_k_d_thry_fit = zeros(length(N_Kuhns),1);

    % OverwriteR2s = 0;

    for nk = 1:length(N_Kuhns)
        N_Kuhn = N_Kuhns(nk);

        FileTag = ['Nk.',num2str(N_Kuhn),'.ea.',num2str(eaStar,'%.2f'),...
            '.bond lifetimes hist'];
        SampName = [OutputFolder,'/',AddOn,FileTag,'.png'];
        if ~isfile(SampName) || Override || ~isfile(R2DynamicsFileName)

            PackageTemp = Package(Package(:,12)==N_Kuhn,:);
            Separations = unique(PackageTemp(:,14));

            for sp = 1:length(Separations)
                Separation = Separations(sp);
                if ~AppendScalingData
                    FileTag = ['Nk.',num2str(N_Kuhn),'.ea.',num2str(eaStar,'%.2f'),...
                        '.d.',num2str(Separation,'%.3f'),...
                        '.bond lifetimes hist'];
                else
                    FileTag = ['Nk.',num2str(N_Kuhn),'.ea.',num2str(eaStar,'%.2f'),...
                        '.d.',num2str(Separation,'%.3f'),...
                        '.bond lifetimes hist_thry'];
                end
                SampName = [OutputFolder,'/',AddOn,FileTag,'.png'];
                if ~isfile(SampName) || Override || ~isfile(R2DynamicsFileName)

                    if NormalizeRates==1
                        tau_d1_meso_temp = tau_d1_meso(:,sp,nk)/tau0;
                        tau_a_meso_temp = tau_a_meso(:,sp,nk)/tau0;
                        tau_d_meso_temp = tau_d_meso(:,sp,nk)/tau0;
                        tau_d1_bead_temp = tau_d1_bead(:,sp,nk)/tau0;
                        tau_a_bead_temp = tau_a_bead(:,sp,nk)/tau0;
                        tau_d_bead_temp = tau_d_bead(:,sp,nk)/tau0;
                        tau_d1_thry_temp = tau_d1_thry(:,sp,nk)/tau0;
                        tau_a_thry_temp = tau_a_thry(:,sp,nk)/tau0;
                        tau_d_thry_temp = tau_d_thry(:,sp,nk)/tau0;
                        xlab_units = ' ($\tau_0$)';
                        High = 0.5e3;
                    else
                        tau_d1_meso_temp = tau_d1_meso(:,sp,nk)*1e6;
                        tau_a_meso_temp = tau_a_meso(:,sp,nk)*1e6;
                        tau_d_meso_temp = tau_d_meso(:,sp,nk)*1e6;
                        tau_d1_bead_temp = tau_d1_bead(:,sp,nk)*1e6;
                        tau_a_bead_temp = tau_a_bead(:,sp,nk)*1e6;
                        tau_d_bead_temp = tau_d_bead(:,sp,nk)*1e6;
                        tau_d1_thry_temp = tau_d1_thry(:,sp,nk)*1e6;
                        tau_a_thry_temp = tau_a_thry(:,sp,nk)*1e6;
                        tau_d_thry_temp = tau_d_thry(:,sp,nk)*1e6;
                        xlab_units = ' ($\mu$s)';
                        High = 5;
                    end
                    tau_d1_meso_temp(tau_d1_meso_temp==0) = [];
                    tau_a_meso_temp(tau_a_meso_temp==0) = [];
                    tau_d_meso_temp(tau_d_meso_temp==0) = [];
                    tau_d1_bead_temp(tau_d1_bead_temp==0) = [];
                    tau_a_bead_temp(tau_a_bead_temp==0) = [];
                    tau_d_bead_temp(tau_d_bead_temp==0) = [];
                    tau_d1_thry_temp(tau_d1_thry_temp==0) = [];
                    tau_a_thry_temp(tau_a_thry_temp==0) = [];
                    tau_d_thry_temp(tau_d_thry_temp==0) = [];

                    Nbins = 8;
                    Low = 0;
                    Edges = linspace(Low,High,Nbins+1);

                    % Meso tau_d1
                    figure(1); clf; hold on
                    if ~AppendScalingData
                        subplot(2,3,1)
                    else
                        subplot(3,3,1)
                    end

                    hold on

                    h1 = histogram(tau_d1_meso_temp);
                    h1.FaceColor = 'g';
                    h1.Normalization = 'probability';
                    h1.BinEdges = Edges;

                    set(gca,'FontSize',FontSize/1.5)
                    %             pbaspect([4 1 1])
                    pbaspect([1 1 1])
                    ylim([0 1])
                    xlim([Low High])

                    ylabel('$P$','FontSize',FontSize,'Interpreter','latex')

                    if ~AppendScalingData
                        subplot(2,3,4)
                        xlabel(['$\tau_{d,1}$',xlab_units],'FontSize',FontSize,'Interpreter','latex')
                    else
                        subplot(3,3,4)
                    end
                    hold on

                    h1b = histogram(tau_d1_bead_temp);
                    h1b.FaceColor = 0.67*[0 1 0];
                    h1b.Normalization = 'probability';
                    h1b.BinEdges = Edges;

                    set(gca,'FontSize',FontSize/1.5)
                    set(gcf,'Color','w')

                    %             pbaspect([4 1 1])
                    pbaspect([1 1 1])
                    ylim([0 1])
                    xlim([Low High])

                    ylabel('$P$','FontSize',FontSize,'Interpreter','latex')

                    % Meso tau_a
                    %             figure(2); clf; hold on
                    if ~AppendScalingData
                        subplot(2,3,3)
                    else
                        subplot(3,3,3)
                    end
                    hold on

                    h2 = histogram(tau_a_meso_temp);
                    h2.FaceColor = 'c';
                    h2.Normalization = 'probability';
                    h2.BinEdges = Edges;

                    set(gca,'FontSize',FontSize/1.5)
                    %             pbaspect([4 1 1])
                    pbaspect([1 1 1])
                    ylim([0 1])
                    xlim([Low High])

                    %             ylabel('$P$','FontSize',FontSize,'Interpreter','latex')

                    if ~AppendScalingData
                        subplot(2,3,6)
                        xlabel(['$\tau_a$',xlab_units],'FontSize',FontSize,'Interpreter','latex')
                    else
                        subplot(3,3,6)
                    end
                    hold on

                    h2b = histogram(tau_a_bead_temp);
                    h2b.FaceColor = 0.67*[0 1 1];
                    h2b.Normalization = 'probability';
                    h2b.BinEdges = Edges;

                    set(gca,'FontSize',FontSize/1.5)
                    set(gcf,'Color','w')

                    %             pbaspect([4 1 1])
                    pbaspect([1 1 1])
                    ylim([0 1])
                    xlim([Low High])

                    %             ylabel('$P$','FontSize',FontSize,'Interpreter','latex')

                    % Meso tau_d
                    %             figure(3); clf; hold on
                    if ~AppendScalingData
                        subplot(2,3,2)
                        xlabel(['$\tau_a$',xlab_units],'FontSize',FontSize,'Interpreter','latex')
                    else
                        subplot(3,3,2)
                    end
                    hold on

                    h3 = histogram(tau_d_meso_temp);
                    h3.FaceColor = 'r';
                    h3.Normalization = 'probability';
                    h3.BinEdges = Edges;

                    set(gca,'FontSize',FontSize/1.5)
                    %             pbaspect([4 1 1])
                    pbaspect([1 1 1])
                    ylim([0 1])
                    xlim([Low High])

                    title(['$N$ = ',num2str(N_Kuhn),', $d^*$ = ',...
                        num2str(Separation/N_Kuhn/b,'%.2f')],...
                        'FontSize',FontSize/2,'Interpreter','latex')

                    %             ylabel('$P$','FontSize',FontSize,'Interpreter','latex')

                    if ~AppendScalingData
                        subplot(2,3,5)
                        xlabel(['$\tau_d$',xlab_units],'FontSize',FontSize,'Interpreter','latex')
                    else
                        subplot(3,3,5)
                    end
                    hold on

                    h3b = histogram(tau_d_bead_temp);
                    h3b.FaceColor = 0.67*[1 0 0];
                    h3b.Normalization = 'probability';
                    h3b.BinEdges = Edges;

                    set(gca,'FontSize',FontSize/1.5)
                    set(gcf,'Color','w')

                    %             pbaspect([4 1 1])
                    pbaspect([1 1 1])
                    ylim([0 1])
                    xlim([Low High])

                    %             ylabel('$P$','FontSize',FontSize,'Interpreter','latex')

                    set(gcf,'Position',[100 100 700 450])


                    if AppendScalingData
                        subplot(3,3,7)
                        h1c = histogram(tau_d1_thry_temp);
                        h1c.FaceColor = 0.33*[0 1 0];
                        h1c.Normalization = 'probability';
                        h1c.BinEdges = Edges;

                        set(gca,'FontSize',FontSize/1.5)
                        pbaspect([1 1 1])
                        ylim([0 1])
                        xlim([Low High])
                        box off

                        xlabel(['$\tau_{d,1}$',xlab_units],'FontSize',FontSize,'Interpreter','latex')
                        ylabel('$P$','FontSize',FontSize,'Interpreter','latex')


                        subplot(3,3,9)
                        h2c = histogram(tau_a_thry_temp);
                        h2c.FaceColor = 0.33*[0 1 1];
                        h2c.Normalization = 'probability';
                        h2c.BinEdges = Edges;

                        set(gca,'FontSize',FontSize/1.5)
                        pbaspect([1 1 1])
                        ylim([0 1])
                        xlim([Low High])
                        box off

                        xlabel(['$\tau_a$',xlab_units],'FontSize',FontSize,'Interpreter','latex')

                        
                        subplot(3,3,8)
                        h3 = histogram(tau_d_thry_temp);
                        h3.FaceColor = 0.33*[1 0 0];
                        h3.Normalization = 'probability';
                        h3.BinEdges = Edges;

                        set(gca,'FontSize',FontSize/1.5)
                        pbaspect([1 1 1])
                        ylim([0 1])
                        xlim([Low High])
                        box off

                        xlabel(['$\tau_d$',xlab_units],'FontSize',FontSize,'Interpreter','latex')
                    end

                    SampName = [OutputFolder,'/',AddOn,FileTag,'.png'];
                    saveas(gcf,SampName)
                    SampName = [OutputFolder,'/',AddOn,FileTag,'.fig'];
                    saveas(gcf,SampName)
                end
            end
            %% Same thing for all separations
            if NormalizeRates==1
                tau_d1_meso_temp = tau_d1_meso(:,:,nk)/tau0;
                tau_a_meso_temp = tau_a_meso(:,:,nk)/tau0;
                tau_d_meso_temp = tau_d_meso(:,:,nk)/tau0;
                tau_d1_bead_temp = tau_d1_bead(:,:,nk)/tau0;
                tau_a_bead_temp = tau_a_bead(:,:,nk)/tau0;
                tau_d_bead_temp = tau_d_bead(:,:,nk)/tau0;
                xlab_units = ' ($\tau_0$)';
                High = 0.5e3;
                tau_nom = -3e3;
                k_nom = -3e-3;
                p_nom = 0.25;
            else
                tau_d1_meso_temp = tau_d1_meso(:,:,nk)*1e6;
                tau_a_meso_temp = tau_a_meso(:,:,nk)*1e6;
                tau_d_meso_temp = tau_d_meso(:,:,nk)*1e6;
                tau_d1_bead_temp = tau_d1_bead(:,:,nk)*1e6;
                tau_a_bead_temp = tau_a_bead(:,:,nk)*1e6;
                tau_d_bead_temp = tau_d_bead(:,:,nk)*1e6;
                xlab_units = ' ($\mu$s)';
                High = 1;
                tau_nom = -3;
                p_nom = 0.25;
            end
            tau_d1_meso_temp = tau_d1_meso_temp(:);
            tau_a_meso_temp = tau_a_meso_temp(:);
            tau_d_meso_temp = tau_d_meso_temp(:);
            tau_d1_bead_temp = tau_d1_bead_temp(:);
            tau_a_bead_temp = tau_a_bead_temp(:);
            tau_d_bead_temp = tau_d_bead_temp(:);
            tau_d1_meso_temp(tau_d1_meso_temp==0) = [];
            tau_a_meso_temp(tau_a_meso_temp==0) = [];
            tau_d_meso_temp(tau_d_meso_temp==0) = [];
            tau_d1_bead_temp(tau_d1_bead_temp==0) = [];
            tau_a_bead_temp(tau_a_bead_temp==0) = [];
            tau_d_bead_temp(tau_d_bead_temp==0) = [];

            Nbins = 10;
            Low = 0;
            yhi = 0.75;
            Edges = linspace(Low,High,Nbins+1);

            % Meso tau_d1
            figure(2); clf; hold on
            if ~AppendScalingData
                subplot(2,3,1)
            else
                subplot(3,3,1)
            end
            hold on

            h1 = histogram(tau_d1_meso_temp);
            h1.FaceColor = 'g';
            h1.Normalization = 'probability';
            h1.BinEdges = Edges;

            set(gca,'FontSize',FontSize/1.5)
            pbaspect([1 1 1])
            ylim([0 yhi])
            xlim([Low High])

            ylabel('$P$','FontSize',FontSize,'Interpreter','latex')

            if ~AppendScalingData
                subplot(2,3,4)
            xlabel(['$\tau_{d,1}$',xlab_units],'FontSize',FontSize,'Interpreter','latex')
            else
                subplot(3,3,4)
            end
            hold on

            h1b = histogram(tau_d1_bead_temp);
            h1b.FaceColor = 0.67*[0 1 0];
            h1b.Normalization = 'probability';
            h1b.BinEdges = Edges;

            RSS = sum((h1.Values-h1b.Values).^2);
            TSS = sum((h1.Values-(mean(h1.Values))).^2);
            R2a1(nk) = 1-RSS/TSS;

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')

            pbaspect([1 1 1])
            ylim([0 yhi])
            xlim([Low High])

            ylabel('$P$','FontSize',FontSize,'Interpreter','latex')

            % Meso tau_a
            if ~AppendScalingData
                subplot(2,3,3)
            else
                subplot(3,3,3)
            end
            hold on

            h2 = histogram(tau_a_meso_temp);
            h2.FaceColor = 'r';
            h2.Normalization = 'probability';
            h2.BinEdges = Edges;

            set(gca,'FontSize',FontSize/1.5)
            pbaspect([1 1 1])
            ylim([0 yhi])
            xlim([Low High])

            if ~AppendScalingData
                subplot(2,3,6)
            xlabel(['$\tau_a$',xlab_units],'FontSize',FontSize,'Interpreter','latex')
            else
                subplot(3,3,6)
            end
            hold on

            h2b = histogram(tau_a_bead_temp);
            h2b.FaceColor = 0.67*[1 0 0];
            h2b.Normalization = 'probability';
            h2b.BinEdges = Edges;

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')

            pbaspect([1 1 1])
            ylim([0 yhi])
            xlim([Low High])


            RSS = sum((h2.Values-h2b.Values).^2);
            TSS = sum((h2.Values-(mean(h2.Values))).^2);
            R2a(nk) = 1-RSS/TSS;

            % Meso tau_d
            if ~AppendScalingData
                subplot(2,3,2)
            else
                subplot(3,3,2)
            end
            hold on

            h3 = histogram(tau_d_meso_temp);
            h3.FaceColor = 'c';
            h3.Normalization = 'probability';
            h3.BinEdges = Edges;

            set(gca,'FontSize',FontSize/1.5)
            pbaspect([1 1 1])
            ylim([0 1])
            xlim([Low High])

            if ~AppendScalingData
                subplot(2,3,5)
                xlabel(['$\tau_d$',xlab_units],'FontSize',FontSize,'Interpreter','latex')
            else
                subplot(3,3,5)
            end
            hold on

            h3b = histogram(tau_d_bead_temp);
            h3b.FaceColor = 0.67*[0 1 1];
            h3b.Normalization = 'probability';
            h3b.BinEdges = Edges;

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')

            pbaspect([1 1 1])
            ylim([0 1])
            xlim([Low High])


            RSS = sum((h3.Values-h3b.Values).^2);
            TSS = sum((h3.Values-(mean(h3.Values))).^2);
            R2d(nk) = 1-RSS/TSS;

            set(gcf,'Position',[100 100 700 450])

            % Estimate kinetic rates through fitting of exponential
            if ~AppendScalingData
                subplot(2,3,1)
            else
                subplot(3,3,1)
            end
            x = h1.BinEdges;
            xplot = linspace(min(h1.BinEdges),max(h1.BinEdges),100);
            x = diff(x)/2+h1.BinEdges(1:end-1);
            y = h1.Values;
            if NormalizeRates==1
                [p_out,k_out,R2Check,~,dk] = FitExp(p_nom,k_nom,x,y);
                fplot = p_out*exp(xplot*k_out);
                k_a1_meso_fit(nk) = -k_out;
                se_k_a1_meso_fit(nk) = dk;
            else
                ft = fittype('a*exp(-b*t)','indep','t');

                f = fit(x',y',ft,'start',[p_nom,tau_nom]);
                fplot = f(xplot);
                coeffvals = coeffvalues(f);
                k_a1_meso_fit(nk) = coeffvals(2);
                conf = confint(f);
                se_k_a1_meso_fit(nk) = abs(conf(2,2)-conf(1,2));

                RSS = sum((y'-f(x)).^2);
                TSS = sum((y-mean(y)).^2);
                R2Check = 1-RSS/TSS;
            end

            p = plot(xplot,fplot);
            p.Color = 'k';
            p.LineWidth = 2;
            p = plot(xplot,fplot);
            p.Color = 'g';
            p.LineWidth = 1;

            if NormalizeRates==1
                text(50,0.4,['$\bar k_{a,1} \sim$ ',num2str(k_a1_meso_fit(nk),'%.1e'),...
                    '$\pm$',num2str(se_k_a1_meso_fit(nk),'%.1e'),...
                    ' $\tau_0^{-1}$, $R^2$=',num2str(R2Check,'%.2f')],...
                    'FontSize',FontSize/3,'Interpreter','latex')
            else
                text(0.25,0.4,['$\bar k_{a,1} \sim$ ',num2str(k_a1_meso_fit(nk),'%.2f'),...
                    '$\pm$',num2str(se_k_a1_meso_fit(nk),'%.2f'),' MHz, ',...
                    '$R^2$=',num2str(R2Check,'%.2f')],'FontSize',FontSize/2.5,...
                    'Interpreter','latex')
            end

            if ~AppendScalingData
                subplot(2,3,4)
            else
                subplot(3,3,4)
            end
            x = h1b.BinEdges;
            xplot = linspace(min(h1b.BinEdges),max(h1b.BinEdges),100);
            x = diff(x)/2+h1b.BinEdges(1:end-1);
            y = h1b.Values;

            if NormalizeRates==1
                [p_out,k_out,R2Check,~,dk] = FitExp(p_nom,k_nom,x,y);
                fplot = p_out*exp(xplot*k_out);
                k_a1_bead_fit(nk) = -k_out;
                se_k_a1_bead_fit(nk) = dk;
            else
                ft = fittype('a*exp(-b*t)','indep','t');

                f = fit(x',y',ft,'start',[p_nom,tau_nom]);
                fplot = f(xplot);
                coeffvals = coeffvalues(f);
                k_a1_bead_fit(nk) = coeffvals(2);
                conf = confint(f);
                se_k_a1_bead_fit(nk) = abs(conf(2,2)-conf(1,2));

                RSS = sum((y'-f(x)).^2);
                TSS = sum((y-mean(y)).^2);
                R2Check = 1-RSS/TSS;
            end

            p = plot(xplot,fplot);
            p.Color = 'k';
            p.LineWidth = 2;
            p = plot(xplot,fplot);
            p.Color = [0 0.5 0];
            p.LineWidth = 1;

            if NormalizeRates==1
                text(50,0.4,['$\bar k_{a,1} \sim$ ',num2str(k_a1_bead_fit(nk),'%.1e'),...
                    '$\pm$',num2str(se_k_a1_bead_fit(nk),'%.1e'),...
                    ' $\tau_0^{-1}$, $R^2$=',num2str(R2Check,'%.2f')],...
                    'FontSize',FontSize/3,'Interpreter','latex')
            else
                text(0.25,0.4,['$\bar k_{a,1} \sim$ ',num2str(k_a1_bead_fit(nk),'%.2f'),...
                    '$\pm$',num2str(se_k_a1_bead_fit(nk),'%.2f'),' MHz, ',...
                    '$R^2$=',num2str(R2Check,'%.2f')],'FontSize',FontSize/2.5,...
                    'Interpreter','latex')
            end

            if ~AppendScalingData
                subplot(2,3,3)
            else
                subplot(3,3,3)
            end
            x = h2.BinEdges;
            xplot = linspace(min(h2.BinEdges),max(h2.BinEdges),100);
            x = diff(x)/2+h2.BinEdges(1:end-1);
            y = h2.Values;

            if NormalizeRates==1
                [p_out,k_out,R2Check,~,dk] = FitExp(p_nom,k_nom,x,y);
                fplot = p_out*exp(xplot*k_out);
                k_d_meso_fit(nk) = -k_out;
                se_k_d_meso_fit(nk) = dk;
            else
                ft = fittype('a*exp(-b*t)','indep','t');

                f = fit(x',y',ft,'start',[p_nom,tau_nom]);
                fplot = f(xplot);
                coeffvals = coeffvalues(f);
                k_d_meso_fit(nk) = coeffvals(2);
                conf = confint(f);
                se_k_d_meso_fit(nk) = abs(conf(2,2)-conf(1,2));

                RSS = sum((y'-f(x)).^2);
                TSS = sum((y-mean(y)).^2);
                R2Check = 1-RSS/TSS;
            end

            p = plot(xplot,fplot);
            p.Color = 'k';
            p.LineWidth = 2;
            p = plot(xplot,fplot);
            p.Color = [1 0 0];
            p.LineWidth = 1;

            if NormalizeRates==1
                text(50,0.4,['$\bar k_d \sim$ ',num2str(k_d_meso_fit(nk),'%.1e'),...
                    '$\pm$',num2str(se_k_d_meso_fit(nk),'%.1e'),...
                    ' $\tau_0^{-1}$, $R^2$=',num2str(R2Check,'%.2f')],...
                    'FontSize',FontSize/3,'Interpreter','latex')
            else
                text(0.25,0.4,['$\bar k_d \sim$ ',num2str(k_d_meso_fit(nk),'%.2f'),...
                    '$\pm$',num2str(se_k_d_meso_fit(nk),'%.2f'),' MHz, ',...
                    '$R^2$=',num2str(R2Check,'%.2f')],'FontSize',FontSize/2.5,...
                    'Interpreter','latex')
            end

            if ~AppendScalingData
                subplot(2,3,6)
            else
                subplot(3,3,6)
            end
            x = h2b.BinEdges;
            xplot = linspace(min(h2b.BinEdges),max(h2b.BinEdges),100);
            x = diff(x)/2+h2b.BinEdges(1:end-1);
            y = h2b.Values;

            if NormalizeRates==1
                [p_out,k_out,R2Check,~,dk] = FitExp(p_nom,k_nom,x,y);
                fplot = p_out*exp(xplot*k_out);
                k_d_bead_fit(nk) = -k_out;
                se_k_d_bead_fit(nk) = dk;
            else
                ft = fittype('a*exp(-b*t)','indep','t');

                f = fit(x',y',ft,'start',[p_nom,tau_nom]);
                fplot = f(xplot);
                coeffvals = coeffvalues(f);
                k_d_bead_fit(nk) = coeffvals(2);
                conf = confint(f);
                se_k_d_bead_fit(nk) = abs(conf(2,2)-conf(1,2));

                RSS = sum((y'-f(x)).^2);
                TSS = sum((y-mean(y)).^2);
                R2Check = 1-RSS/TSS;
            end

            p = plot(xplot,fplot);
            p.Color = 'k';
            p.LineWidth = 2;
            p = plot(xplot,fplot);
            p.Color = [0.5 0 0];
            p.LineWidth = 1;

            if NormalizeRates==1
                text(50,0.4,['$\bar k_d \sim$ ',num2str(k_d_bead_fit(nk),'%.1e'),...
                    '$\pm$',num2str(se_k_d_bead_fit(nk),'%.1e'),...
                    ' $\tau_0^{-1}$, $R^2$=',num2str(R2Check,'%.2f')],...
                    'FontSize',FontSize/3,'Interpreter','latex')
            else
                text(0.25,0.4,['$\bar k_d \sim$ ',num2str(k_d_bead_fit(nk),'%.2f'),...
                    '$\pm$',num2str(se_k_d_bead_fit(nk),'%.2f'),' MHz, ',...
                    '$R^2$=',num2str(R2Check,'%.2f')],'FontSize',FontSize/2.5,...
                    'Interpreter','latex')
            end

            if ~AppendScalingData
                subplot(2,3,2)
            else
                subplot(3,3,2)
            end
            x = h3.BinEdges;
            xplot = linspace(min(h3.BinEdges),max(h3.BinEdges),100);
            x = diff(x)/2+h3.BinEdges(1:end-1);
            y = h3.Values;

            if NormalizeRates==1
                [p_out,k_out,R2Check,~,dk] = FitExp(p_nom,k_nom,x,y);
                fplot = p_out*exp(xplot*k_out);
                k_a_meso_fit(nk) = -k_out;
                se_k_d_meso_fit(nk) = dk;
            else
                ft = fittype('a*exp(-b*t)','indep','t');

                f = fit(x',y',ft,'start',[p_nom,tau_nom]);
                fplot = f(xplot);
                coeffvals = coeffvalues(f);
                k_a_meso_fit(nk) = coeffvals(2);
                conf = confint(f);
                se_k_d_meso_fit(nk) = abs(conf(2,2)-conf(1,2));

                RSS = sum((y'-f(x)).^2);
                TSS = sum((y-mean(y)).^2);
                R2Check = 1-RSS/TSS;
            end

            p = plot(xplot,fplot);
            p.Color = 'k';
            p.LineWidth = 2;
            p = plot(xplot,fplot);
            p.Color = [0 1 1];
            p.LineWidth = 1;

            if NormalizeRates==1
                text(50,0.4,['$\bar k_d \sim$ ',num2str(k_a_meso_fit(nk),'%.1e'),...
                    '$\pm$',num2str(se_k_a_meso_fit(nk),'%.1e'),...
                    ' $\tau_0^{-1}$, $R^2$=',num2str(R2Check,'%.2f')],...
                    'FontSize',FontSize/3,'Interpreter','latex')
            else
                text(0.25,0.4,['$\bar k_a \sim$ ',num2str(k_a_meso_fit(nk),'%.2f'),...
                    '$\pm$',num2str(se_k_a_meso_fit(nk),'%.2f'),' MHz, ',...
                    '$R^2$=',num2str(R2Check,'%.2f')],'FontSize',FontSize/2.5,...
                    'Interpreter','latex')
            end

            if ~AppendScalingData
                subplot(2,3,5)
            else
                subplot(3,3,5)
            end
            x = h3b.BinEdges;
            xplot = linspace(min(h3b.BinEdges),max(h3b.BinEdges),100);
            x = diff(x)/2+h3b.BinEdges(1:end-1);
            y = h3b.Values;

            if NormalizeRates==1
                [p_out,k_out,R2Check,~,dk] = FitExp(p_nom,k_nom,x,y);
                fplot = p_out*exp(xplot*k_out);
                k_a_bead_fit(nk) = -k_out;
                se_k_a_bead_fit(nk) = dk;
            else
                ft = fittype('a*exp(-b*t)','indep','t');

                f = fit(x',y',ft,'start',[p_nom,tau_nom]);
                fplot = f(xplot);
                coeffvals = coeffvalues(f);
                k_a_bead_fit(nk) = coeffvals(2);
                conf = confint(f);
                se_k_a_bead_fit(nk) = abs(conf(2,2)-conf(1,2));

                RSS = sum((y'-f(x)).^2);
                TSS = sum((y-mean(y)).^2);
                R2Check = 1-RSS/TSS;
            end

            p = plot(xplot,fplot);
            p.Color = 'k';
            p.LineWidth = 2;
            p = plot(xplot,fplot);
            p.Color = [0 0.5 0.5];
            p.LineWidth = 1;

            if NormalizeRates==1
                text(50,0.4,['$\bar k_d \sim$ ',num2str(k_a_bead_fit(nk),'%.1e'),...
                    '$\pm$',num2str(se_k_a_bead_fit(nk),'%.1e'),...
                    ' $\tau_0^{-1}$, $R^2$=',num2str(R2Check,'%.2f')],...
                    'FontSize',FontSize/3,'Interpreter','latex')
            else
                text(0.25,0.4,['$\bar k_a \sim$ ',num2str(k_a_bead_fit(nk),'%.2f'),...
                    '$\pm$',num2str(se_k_a_bead_fit(nk),'%.2f'),' MHz, ',...
                    '$R^2$=',num2str(R2Check,'%.2f')],'FontSize',FontSize/2.5,...
                    'Interpreter','latex')
            end

            if ~AppendScalingData
                subplot(2,3,2)
            else
                subplot(3,3,2)
            end
            title(['$N$ = ',num2str(N_Kuhn)],...
                'FontSize',FontSize/2,'Interpreter','latex')

            if AppendScalingData
                subplot(3,3,7)
                hold on

                h1c = histogram(tau_d1_meso_temp);
                h1c.FaceColor = 0.33*[0 1 0];
                h1c.Normalization = 'probability';
                h1c.BinEdges = Edges;

                set(gca,'FontSize',FontSize/1.5)
                pbaspect([1 1 1])
                ylim([0 yhi])
                xlim([Low High])

                ylabel('$P$','FontSize',FontSize,'Interpreter','latex')
                xlabel(['$\tau_{d,1}$',xlab_units],'FontSize',FontSize,'Interpreter','latex')


                % Meso tau_a
                subplot(3,3,9)
                hold on

                h2c = histogram(tau_a_meso_temp);
                h2c.FaceColor = 0.33*[1 0 0];
                h2c.Normalization = 'probability';
                h2c.BinEdges = Edges;

                set(gca,'FontSize',FontSize/1.5)
                pbaspect([1 1 1])
                ylim([0 yhi])
                xlim([Low High])

                xlabel(['$\tau_a$',xlab_units],'FontSize',FontSize,'Interpreter','latex')


                % Meso tau_d
                subplot(3,3,8)
                hold on

                h3c = histogram(tau_d_meso_temp);
                h3c.FaceColor = 0.33*[0 1 1];
                h3c.Normalization = 'probability';
                h3c.BinEdges = Edges;

                set(gca,'FontSize',FontSize/1.5)
                pbaspect([1 1 1])
                ylim([0 1])
                xlim([Low High])

                xlabel(['$\tau_d$',xlab_units],'FontSize',FontSize,'Interpreter','latex')

                set(gcf,'Position',[100 100 700 450])



                % Estimate kinetic rates through fitting of exponential
                subplot(3,3,7)
                x = h1c.BinEdges;
                xplot = linspace(min(h1c.BinEdges),max(h1c.BinEdges),100);
                x = diff(x)/2+h1c.BinEdges(1:end-1);
                y = h1c.Values;
                if NormalizeRates==1
                    [p_out,k_out,R2Check,~,dk] = FitExp(p_nom,k_nom,x,y);
                    fplot = p_out*exp(xplot*k_out);
                    k_a1_thry_fit(nk) = -k_out;
                    se_k_a1_thry_fit(nk) = dk;
                else
                    ft = fittype('a*exp(-b*t)','indep','t');

                    f = fit(x',y',ft,'start',[p_nom,tau_nom]);
                    fplot = f(xplot);
                    coeffvals = coeffvalues(f);
                    k_a1_thry_fit(nk) = coeffvals(2);
                    conf = confint(f);
                    se_k_a1_thry_fit(nk) = abs(conf(2,2)-conf(1,2));

                    RSS = sum((y'-f(x)).^2);
                    TSS = sum((y-mean(y)).^2);
                    R2Check = 1-RSS/TSS;
                end

                p = plot(xplot,fplot);
                p.Color = 'k';
                p.LineWidth = 2;
                p = plot(xplot,fplot);
                p.Color = 0.33*[0 1 0];
                p.LineWidth = 1;

                if NormalizeRates==1
                    text(50,0.4,['$\bar k_{a,1} \sim$ ',num2str(k_a1_thry_fit(nk),'%.1e'),...
                        '$\pm$',num2str(se_k_a1_thry_fit(nk),'%.1e'),...
                        ' $\tau_0^{-1}$, $R^2$=',num2str(R2Check,'%.2f')],...
                        'FontSize',FontSize/3,'Interpreter','latex')
                else
                    text(0.25,0.4,['$\bar k_{a,1} \sim$ ',num2str(k_a1_thry_fit(nk),'%.2f'),...
                        '$\pm$',num2str(se_k_a1_thry_fit(nk),'%.2f'),' MHz, ',...
                        '$R^2$=',num2str(R2Check,'%.2f')],'FontSize',FontSize/2.5,...
                        'Interpreter','latex')
                end

                subplot(3,3,9)
                x = h2c.BinEdges;
                xplot = linspace(min(h2c.BinEdges),max(h2c.BinEdges),100);
                x = diff(x)/2+h2c.BinEdges(1:end-1);
                y = h2c.Values;
                if NormalizeRates==1
                    [p_out,k_out,R2Check,~,dk] = FitExp(p_nom,k_nom,x,y);
                    fplot = p_out*exp(xplot*k_out);
                    k_d_thry_fit(nk) = -k_out;
                    se_k_d_thry_fit(nk) = dk;
                else
                    ft = fittype('a*exp(-b*t)','indep','t');

                    f = fit(x',y',ft,'start',[p_nom,tau_nom]);
                    fplot = f(xplot);
                    coeffvals = coeffvalues(f);
                    k_d_thry_fit(nk) = coeffvals(2);
                    conf = confint(f);
                    se_k_d_thry_fit(nk) = abs(conf(2,2)-conf(1,2));

                    RSS = sum((y'-f(x)).^2);
                    TSS = sum((y-mean(y)).^2);
                    R2Check = 1-RSS/TSS;
                end

                p = plot(xplot,fplot);
                p.Color = 'k';
                p.LineWidth = 2;
                p = plot(xplot,fplot);
                p.Color = 0.33*[1 0 0];
                p.LineWidth = 1;

                if NormalizeRates==1
                    text(50,0.4,['$\bar k_{a,1} \sim$ ',num2str(k_d_thry_fit(nk),'%.1e'),...
                        '$\pm$',num2str(se_k_d_thry_fit(nk),'%.1e'),...
                        ' $\tau_0^{-1}$, $R^2$=',num2str(R2Check,'%.2f')],...
                        'FontSize',FontSize/3,'Interpreter','latex')
                else
                    text(0.25,0.4,['$\bar k_{a,1} \sim$ ',num2str(k_d_thry_fit(nk),'%.2f'),...
                        '$\pm$',num2str(se_k_d_thry_fit(nk),'%.2f'),' MHz, ',...
                        '$R^2$=',num2str(R2Check,'%.2f')],'FontSize',FontSize/2.5,...
                        'Interpreter','latex')
                end

                subplot(3,3,8)
                x = h3c.BinEdges;
                xplot = linspace(min(h3c.BinEdges),max(h3c.BinEdges),100);
                x = diff(x)/2+h3c.BinEdges(1:end-1);
                y = h3c.Values;
                if NormalizeRates==1
                    [p_out,k_out,R2Check,~,dk] = FitExp(p_nom,k_nom,x,y);
                    fplot = p_out*exp(xplot*k_out);
                    k_a_thry_fit(nk) = -k_out;
                    se_k_a_thry_fit(nk) = dk;
                else
                    ft = fittype('a*exp(-b*t)','indep','t');

                    f = fit(x',y',ft,'start',[p_nom,tau_nom]);
                    fplot = f(xplot);
                    coeffvals = coeffvalues(f);
                    k_a_thry_fit(nk) = coeffvals(2);
                    conf = confint(f);
                    se_k_a_thry_fit(nk) = abs(conf(2,2)-conf(1,2));

                    RSS = sum((y'-f(x)).^2);
                    TSS = sum((y-mean(y)).^2);
                    R2Check = 1-RSS/TSS;
                end

                p = plot(xplot,fplot);
                p.Color = 'k';
                p.LineWidth = 2;
                p = plot(xplot,fplot);
                p.Color = 0.33*[0 1 1];
                p.LineWidth = 1;

                if NormalizeRates==1
                    text(50,0.4,['$\bar k_{a,1} \sim$ ',num2str(k_a_thry_fit(nk),'%.1e'),...
                        '$\pm$',num2str(se_k_a_thry_fit(nk),'%.1e'),...
                        ' $\tau_0^{-1}$, $R^2$=',num2str(R2Check,'%.2f')],...
                        'FontSize',FontSize/3,'Interpreter','latex')
                else
                    text(0.25,0.4,['$\bar k_{a,1} \sim$ ',num2str(k_a_thry_fit(nk),'%.2f'),...
                        '$\pm$',num2str(se_k_a_thry_fit(nk),'%.2f'),' MHz, ',...
                        '$R^2$=',num2str(R2Check,'%.2f')],'FontSize',FontSize/2.5,...
                        'Interpreter','latex')
                end

            end


            if ~AppendScalingData
                FileTag = ['Nk.',num2str(N_Kuhn),'.ea.',num2str(eaStar,'%.2f'),...
                    '.bond lifetimes hist'];
            else
                FileTag = ['Nk.',num2str(N_Kuhn),'.ea.',num2str(eaStar,'%.2f'),...
                    '.bond lifetimes hist_thry'];
            end
            SampName = [OutputFolder,'/',AddOn,FileTag,'.png'];
            saveas(gcf,SampName)
            SampName = [OutputFolder,'/',AddOn,FileTag,'.fig'];
            saveas(gcf,SampName)
        end
    end


    R2Dynamics.R2a1 = R2a1;
    R2Dynamics.R2a = R2a;
    R2Dynamics.R2d = R2d;

    R2Dynamics.k_a1_meso_fit = k_a1_meso_fit;
    R2Dynamics.k_a_meso_fit = k_a_meso_fit;
    R2Dynamics.k_d_meso_fit = k_d_meso_fit;
    R2Dynamics.k_a1_bead_fit = k_a1_bead_fit;
    R2Dynamics.k_a_bead_fit = k_a_bead_fit;
    R2Dynamics.k_d_bead_fit = k_d_bead_fit;

    R2Dynamics.se_k_a1_meso_fit = se_k_a1_meso_fit;
    R2Dynamics.se_k_a_meso_fit = se_k_a_meso_fit;
    R2Dynamics.se_k_d_meso_fit = se_k_d_meso_fit;
    R2Dynamics.se_k_a1_bead_fit = se_k_a1_bead_fit;
    R2Dynamics.se_k_a_bead_fit = se_k_a_bead_fit;
    R2Dynamics.se_k_d_bead_fit = se_k_d_bead_fit;

    save(R2DynamicsFileName,'-struct','R2Dynamics');
else
    R2Dynamics = load(R2DynamicsFileName,'-mat');
    R2a1 = R2Dynamics.R2a1;
    R2a = R2Dynamics.R2a;
    R2d = R2Dynamics.R2d;

    k_a1_meso_fit = R2Dynamics.k_a1_meso_fit;
    k_a_meso_fit = R2Dynamics.k_a_meso_fit;
    k_d_meso_fit = R2Dynamics.k_d_meso_fit;
    k_a1_bead_fit = R2Dynamics.k_a1_bead_fit;
    k_a_bead_fit = R2Dynamics.k_a_bead_fit;
    k_d_bead_fit = R2Dynamics.k_d_bead_fit;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotDynamicsWrtTime(Override,Package,N_Kuhns,b)

global AllDynamicsFileName FontSize OutputFolder timenorm NormalizeRates...
    AppendScalingData

AddOn = 'Wrt. time.';
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

if AppendScalingData==1
    time_thry = DynamicsWrtTime.time_thry;
    ka_thry = DynamicsWrtTime.ka_thry;
    kd_thry = DynamicsWrtTime.kd_thry;
    fa_thry = DynamicsWrtTime.fa_thry;
    fd_thry = DynamicsWrtTime.fd_thry;
end

window = 100;
smoothka_bead = movmean(ka_bead/timenorm,window,1);
smoothkd_bead = movmean(kd_bead/timenorm,window,1);

smoothka_meso = movmean(ka_meso/timenorm,window,1);
smoothkd_meso = movmean(kd_meso/timenorm,window,1);

if AppendScalingData==1
    smoothka_thry = movmean(ka_thry/timenorm,window,1);
    smoothkd_thry = movmean(kd_thry/timenorm,window,1);
end

Npt = 25;
ScatterSpacing = round((size(time_bead,1)-1)/Npt);
rng = (1:ScatterSpacing:size(time_bead,1));
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

            if AppendScalingData
                s1a = scatter(time_thry(rng),smoothka_thry(rng,sp,nk));
                s1a.MarkerEdgeColor = Colorsa(sp,:);
                s1a.Marker = 'd';
                s1a.SizeData = DataSize;
            end

            ylim([0 Inf])

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([2 1 1])

            set(gcf,'Position',[200 200 250 250])

            xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
            if NormalizeRates==1
                ylabel('$k_a^*$ ($\tau_0^{-1}$)','FontSize',FontSize,'Interpreter','latex')
                ylim([0 7.5e-3])
            else
                ylabel('$k_a$ (kHz)','FontSize',FontSize,'Interpreter','latex')
            end

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

            if AppendScalingData
                s2a = scatter(time_thry(rng),smoothkd_thry(rng,sp,nk));
                s2a.MarkerEdgeColor = Colorsd(sp,:);
                s2a.Marker = 'd';
                s2a.SizeData = DataSize;
            end
            
            ylim([0 Inf])

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([2 1 1])

            set(gcf,'Position',[200 200 250 250])

            xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
            if NormalizeRates==1
                ylabel('$k_d^*$ ($\tau_0^{-1}$)','FontSize',FontSize,'Interpreter','latex')
                ylim([0 7.5e-3])
            else
                ylabel('$k_d$ (kHz)','FontSize',FontSize,'Interpreter','latex')
            end

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

            if AppendScalingData
                s3a = scatter(time_thry(rng),fa_thry(rng,sp,nk)*100);
                s3a.MarkerEdgeColor = Colorsa(sp,:);
                s3a.Marker = 'd';
                s3a.SizeData = DataSize;
            end

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

            if AppendScalingData
                s4a = scatter(time_thry(rng),fd_thry(rng,sp,nk)*100);
                s4a.MarkerEdgeColor = Colorsd(sp,:);
                s4a.Marker = 'd';
                s4a.SizeData = DataSize;
            end

            ylim([0 100])

            set(gca,'FontSize',FontSize/1.5)
            set(gcf,'Color','w')
            pbaspect([2 1 1])

            set(gcf,'Position',[200 200 250 250])

            xlabel('$t$ (s)','FontSize',FontSize,'Interpreter','latex')
            ylabel('$f_d$ ($\%$)','FontSize',FontSize,'Interpreter','latex')
        end

        fig = figure(1);
        set(fig,'name',['N = ',num2str(N_Kuhn),',ka vs t'])
        FileTag = ['Nk.',num2str(N_Kuhn),'.ka vs t'];
        saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
        saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])
        l = legend(p1,LegendEntries);
        l.FontSize = FontSize/2;
        l.Interpreter = 'latex';

        fig = figure(2);
        set(fig,'name',['N = ',num2str(N_Kuhn),',kd vs t'])
        FileTag = ['Nk.',num2str(N_Kuhn),'.kd vs t'];
        saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
        saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])
        l = legend(p2,LegendEntries);
        l.FontSize = FontSize/2;
        l.Interpreter = 'latex';

        fig = figure(3);
        set(fig,'name',['N = ',num2str(N_Kuhn),',fa vs t'])
        FileTag = ['Nk.',num2str(N_Kuhn),'.fa vs t'];
        saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
        saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])

        fig = figure(4);
        set(fig,'name',['N = ',num2str(N_Kuhn),',fd vs t'])
        FileTag = ['Nk.',num2str(N_Kuhn),'.fd vs t'];
        saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.png'])
        saveas(gcf,[OutputFolder,'/',AddOn,FileTag,'.fig'])
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AssembleEnsembleDynamicsData(Override,Package,...
    damps,N_Kuhns,Np,Ns,T,Nb,ka_in,kd_in,f0,dts,b,N,AppendScalingData)

global N_Kuhn...
    ka_ens_meso ka1_ens_meso kd_ens_meso...
    kaSE_ens_meso ka1SE_ens_meso kdSE_ens_meso...
    Na_ens_meso Nd_ens_meso NaSE_ens_meso NdSE_ens_meso...
    fa_ens_meso fd_ens_meso faSE_ens_meso fdSE_ens_meso...
    ka_ens_bead ka1_ens_bead kd_ens_bead...
    kaSE_ens_bead ka1SE_ens_bead kdSE_ens_bead...
    Na_ens_bead Nd_ens_bead NaSE_ens_bead NdSE_ens_bead...
    fa_ens_bead fd_ens_bead faSE_ens_bead fdSE_ens_bead...
    ka_ens_thry ka1_ens_thry kd_ens_thry...
    kaSE_ens_thry ka1SE_ens_thry kdSE_ens_thry...
    Na_ens_thry Nd_ens_thry NaSE_ens_thry NdSE_ens_thry...
    fa_ens_thry fd_ens_thry faSE_ens_thry fdSE_ens_thry...
    tau_d1_all_bead tau_a_all_bead tau_d_all_bead...
    tau_d1_all_meso tau_a_all_meso tau_d_all_meso...
    tau_d1_all_thry tau_a_all_thry tau_d_all_thry...
    ka_all_meso ka1_all_meso kd_all_meso...
    Na_all_meso Nd_all_meso...
    fa_all_meso fd_all_meso...
    ka_all_bead ka1_all_bead kd_all_bead...
    Na_all_bead Nd_all_bead...
    fa_all_bead fd_all_bead...
    ka_all_thry ka1_all_thry kd_all_thry...
    Na_all_thry Nd_all_thry...
    fa_all_thry fd_all_thry...
    tau_d1_ens_meso tau_d1SE_ens_meso...
    tau_a_ens_meso tau_aSE_ens_meso...
    tau_d_ens_meso tau_dSE_ens_meso...
    tau_d1_ens_bead tau_d1SE_ens_bead...
    tau_a_ens_bead tau_aSE_ens_bead...
    tau_d_ens_bead tau_dSE_ens_bead...
    tau_d1_ens_thry tau_d1SE_ens_thry...
    tau_a_ens_thry tau_aSE_ens_thry...
    tau_d_ens_thry tau_dSE_ens_thry...
    EnsembleDynamicsFileName AllDynamicsFileName

Folder = 'Output Plots';
if ~isfolder(Folder)
    mkdir(Folder)
end

DefineAssembledFileNames(ka_in);

% MSD wrt time for each damper while sweeping N
Override = 1;
if ~isfile(EnsembleDynamicsFileName) || ~isfile(AllDynamicsFileName) || Override
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
                    BSOM = 0;
                    [time_bead,ka_bead,ka1_bead,kd_bead,Na_bead,Nd_bead,fa_bead,fd_bead,...
                        tau_d1_bead,tau_a_bead,tau_d_bead,...
                        ka_ens_bead(sp,nk),kaSE_ens_bead(sp,nk),...
                        ka1_ens_bead(sp,nk),ka1SE_ens_bead(sp,nk),...
                        kd_ens_bead(sp,nk),kdSE_ens_bead(sp,nk),...
                        Na_ens_bead(sp,nk),NaSE_ens_bead(sp,nk),...
                        Nd_ens_bead(sp,nk),NdSE_ens_bead(sp,nk),...
                        fa_ens_bead(sp,nk),faSE_ens_bead(sp,nk),...
                        fd_ens_bead(sp,nk),fdSE_ens_bead(sp,nk),...
                        tau_d1_ens_bead(sp,nk),tau_d1SE_ens_bead(sp,nk),...
                        tau_a_ens_bead(sp,nk),tau_aSE_ens_bead(sp,nk),...
                        tau_d_ens_bead(sp,nk),tau_dSE_ens_bead(sp,nk)] =...
                        CalculateEnsembleDynamics(BSOM,...
                        Np,N,ka_in,kd_in,f0,dt,damp,N_Kuhn,b,Separation);

                    % Pull in the Mesoscale data
                    BSOM = 1;
                    [time_meso,ka_meso,ka1_meso,kd_meso,Na_meso,Nd_meso,fa_meso,fd_meso,...
                        tau_d1_meso,tau_a_meso,tau_d_meso,...
                        ka_ens_meso(sp,nk),kaSE_ens_meso(sp,nk),...
                        ka1_ens_meso(sp,nk),ka1SE_ens_meso(sp,nk),...
                        kd_ens_meso(sp,nk),kdSE_ens_meso(sp,nk),...
                        Na_ens_meso(sp,nk),NaSE_ens_meso(sp,nk),...
                        Nd_ens_meso(sp,nk),NdSE_ens_meso(sp,nk),...
                        fa_ens_meso(sp,nk),faSE_ens_meso(sp,nk),...
                        fd_ens_meso(sp,nk),fdSE_ens_meso(sp,nk),...
                        tau_d1_ens_meso(sp,nk),tau_d1SE_ens_meso(sp,nk),...
                        tau_a_ens_meso(sp,nk),tau_aSE_ens_meso(sp,nk),...
                        tau_d_ens_meso(sp,nk),tau_dSE_ens_meso(sp,nk)] =...
                        CalculateEnsembleDynamics(BSOM,...
                        Np,N,ka_in,kd_in,f0,dt,damp,N_Kuhn,b,Separation);

                    if AppendScalingData
                        % Pull in the Scaling theory data data
                        BSOM = 2;
                        [time_thry,ka_thry,ka1_thry,kd_thry,Na_thry,Nd_thry,fa_thry,fd_thry,...
                            tau_d1_thry,tau_a_thry,tau_d_thry,...
                            ka_ens_thry(sp,nk),kaSE_ens_thry(sp,nk),...
                            ka1_ens_thry(sp,nk),ka1SE_ens_thry(sp,nk),...
                            kd_ens_thry(sp,nk),kdSE_ens_thry(sp,nk),...
                            Na_ens_thry(sp,nk),NaSE_ens_thry(sp,nk),...
                            Nd_ens_thry(sp,nk),NdSE_ens_thry(sp,nk),...
                            fa_ens_thry(sp,nk),faSE_ens_thry(sp,nk),...
                            fd_ens_thry(sp,nk),fdSE_ens_thry(sp,nk),...
                            tau_d1_ens_thry(sp,nk),tau_d1SE_ens_thry(sp,nk),...
                            tau_a_ens_thry(sp,nk),tau_aSE_ens_thry(sp,nk),...
                            tau_d_ens_thry(sp,nk),tau_dSE_ens_thry(sp,nk)] =...
                            CalculateEnsembleDynamics(BSOM,...
                            Np,N,ka_in,kd_in,f0,dt,damp,N_Kuhn,b,Separation);
                    else
                        time_thry = [];
                    end

                    if sp==1 && nk==1
                        InitializeAllDynamics(time_bead,time_meso,time_thry,size(PackageTemp,1));
                    end

                    ka_all_bead(:,sp,nk) = ka_bead;
                    ka1_all_bead(:,sp,nk) = ka1_bead;
                    kd_all_bead(:,sp,nk) = kd_bead;
                    Na_all_bead(:,sp,nk) = Na_bead;
                    Nd_all_bead(:,sp,nk) = Nd_bead;
                    fa_all_bead(:,sp,nk) = fa_bead;
                    fd_all_bead(:,sp,nk) = fd_bead;
                    
                    tau_d1_all_bead(1:length(tau_d1_bead),sp,nk) = tau_d1_bead;
                    tau_a_all_bead(1:length(tau_a_bead),sp,nk) = tau_a_bead;
                    tau_d_all_bead(1:length(tau_d_bead),sp,nk) = tau_d_bead;
                    
                    
                    ka_all_meso(:,sp,nk) = ka_meso;
                    ka1_all_meso(:,sp,nk) = ka1_meso;
                    kd_all_meso(:,sp,nk) = kd_meso;
                    Na_all_meso(:,sp,nk) = Na_meso;
                    Nd_all_meso(:,sp,nk) = Nd_meso;
                    fa_all_meso(:,sp,nk) = fa_meso;
                    fd_all_meso(:,sp,nk) = fd_meso;
                    
                    tau_d1_all_meso(1:length(tau_d1_meso),sp,nk) = tau_d1_meso;
                    tau_a_all_meso(1:length(tau_a_meso),sp,nk) = tau_a_meso;
                    tau_d_all_meso(1:length(tau_d_meso),sp,nk) = tau_d_meso;

                    if AppendScalingData
                        ka_all_thry(:,sp,nk) = ka_thry;
                        ka1_all_thry(:,sp,nk) = ka1_thry;
                        kd_all_thry(:,sp,nk) = kd_thry;
                        Na_all_thry(:,sp,nk) = Na_thry;
                        Nd_all_thry(:,sp,nk) = Nd_thry;
                        fa_all_thry(:,sp,nk) = fa_thry;
                        fd_all_thry(:,sp,nk) = fd_thry;

                        tau_d1_all_thry(1:length(tau_d1_thry),sp,nk) = tau_d1_thry;
                        tau_a_all_thry(1:length(tau_a_thry),sp,nk) = tau_a_thry;
                        tau_d_all_thry(1:length(tau_d_thry),sp,nk) = tau_d_thry;
                    end
                end
            end
        end
    end

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

    EnsembleDynamics.tau_d1_bead = tau_d1_ens_bead;
    EnsembleDynamics.tau_d1SE_bead = tau_d1SE_ens_bead;
    EnsembleDynamics.tau_a_bead = tau_a_ens_bead;
    EnsembleDynamics.tau_aSE_bead = tau_aSE_ens_bead;
    EnsembleDynamics.tau_d_bead = tau_d_ens_bead;
    EnsembleDynamics.tau_dSE_bead = tau_dSE_ens_bead;


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
    
    EnsembleDynamics.tau_d1_meso = tau_d1_ens_meso;
    EnsembleDynamics.tau_d1SE_meso = tau_d1SE_ens_meso;
    EnsembleDynamics.tau_a_meso = tau_a_ens_meso;
    EnsembleDynamics.tau_aSE_meso = tau_aSE_ens_meso;
    EnsembleDynamics.tau_d_meso = tau_d_ens_meso;
    EnsembleDynamics.tau_dSE_meso = tau_dSE_ens_meso;


    if AppendScalingData
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

        EnsembleDynamics.tau_d1_thry = tau_d1_ens_thry;
        EnsembleDynamics.tau_d1SE_thry = tau_d1SE_ens_thry;
        EnsembleDynamics.tau_a_thry = tau_a_ens_thry;
        EnsembleDynamics.tau_aSE_thry = tau_aSE_ens_thry;
        EnsembleDynamics.tau_d_thry = tau_d_ens_thry;
        EnsembleDynamics.tau_dSE_thry = tau_dSE_ens_thry;
    end
    
    save(EnsembleDynamicsFileName,'-struct','EnsembleDynamics');


    AllDynamics.time_bead = time_bead;
    AllDynamics.ka_bead = ka_all_bead;
    AllDynamics.ka1_bead = ka1_all_bead;
    AllDynamics.kd_bead = kd_all_bead;
    AllDynamics.Na_bead = Na_all_bead;
    AllDynamics.Nd_bead = Nd_all_bead;
    AllDynamics.fa_bead = fa_all_bead;
    AllDynamics.fd_bead = fd_all_bead;
    AllDynamics.tau_d1_bead = tau_d1_all_bead;
    AllDynamics.tau_a_bead = tau_a_all_bead;
    AllDynamics.tau_d_bead = tau_d_all_bead;

    AllDynamics.time_meso = time_meso;
    AllDynamics.ka_meso = ka_all_meso;
    AllDynamics.ka1_meso = ka1_all_meso;
    AllDynamics.kd_meso = kd_all_meso;
    AllDynamics.Na_meso = Na_all_meso;
    AllDynamics.Nd_meso = Nd_all_meso;
    AllDynamics.fa_meso = fa_all_meso;
    AllDynamics.fd_meso = fd_all_meso;
    AllDynamics.tau_d1_meso = tau_d1_all_meso;
    AllDynamics.tau_a_meso = tau_a_all_meso;
    AllDynamics.tau_d_meso = tau_d_all_meso;

    if AppendScalingData
        AllDynamics.time_thry = time_thry;
        AllDynamics.ka_thry = ka_all_thry;
        AllDynamics.ka1_thry = ka1_all_thry;
        AllDynamics.kd_thry = kd_all_thry;
        AllDynamics.Na_thry = Na_all_thry;
        AllDynamics.Nd_thry = Nd_all_thry;
        AllDynamics.fa_thry = fa_all_thry;
        AllDynamics.fd_thry = fd_all_thry;
        AllDynamics.tau_d1_thry = tau_d1_all_thry;
        AllDynamics.tau_a_thry = tau_a_all_thry;
        AllDynamics.tau_d_thry = tau_d_all_thry;
    end

    save(AllDynamicsFileName,'-struct','AllDynamics');
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

    tau_d1_ens_bead = EnsembleDynamics.tau_d1_bead;
    tau_d1SE_ens_bead = EnsembleDynamics.tau_d1SE_bead;
    tau_a_ens_bead = EnsembleDynamics.tau_a_bead;
    tau_aSE_ens_bead = EnsembleDynamics.tau_aSE_bead;
    tau_d_ens_bead = EnsembleDynamics.tau_d_bead;
    tau_dSE_ens_bead = EnsembleDynamics.tau_dSE_bead;

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

    tau_d1_ens_meso = EnsembleDynamics.tau_d1_meso;
    tau_d1SE_ens_meso = EnsembleDynamics.tau_d1SE_meso;
    tau_a_ens_meso = EnsembleDynamics.tau_a_meso;
    tau_aSE_ens_meso = EnsembleDynamics.tau_aSE_meso;
    tau_d_ens_meso = EnsembleDynamics.tau_d_meso;
    tau_dSE_ens_meso = EnsembleDynamics.tau_dSE_meso;


    if AppendScalingData
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

        tau_d1_ens_thry = EnsembleDynamics.tau_d1_thry;
        tau_d1SE_ens_thry = EnsembleDynamics.tau_d1SE_thry;
        tau_a_ens_thry = EnsembleDynamics.tau_a_thry;
        tau_aSE_ens_thry = EnsembleDynamics.tau_aSE_thry;
        tau_d_ens_thry = EnsembleDynamics.tau_d_thry;
        tau_dSE_ens_thry = EnsembleDynamics.tau_dSE_thry;
    end

    AllDynamics = load(AllDynamicsFileName,'-mat');

    ka_all_bead = AllDynamics.ka_bead;
    ka1_all_bead = AllDynamics.ka1_bead;
    kd_all_bead = AllDynamics.kd_bead;
    Na_all_bead = AllDynamics.Na_bead;
    Nd_all_bead = AllDynamics.Nd_bead;
    fa_all_bead = AllDynamics.fa_bead;
    fd_all_bead = AllDynamics.fd_bead;
    tau_d1_all_bead = AllDynamics.tau_d1_bead;
    tau_a_all_bead = AllDynamics.tau_a_bead;
    tau_d_all_bead = AllDynamics.tau_d_bead;

    ka_all_meso = AllDynamics.ka_meso;
    ka1_all_meso = AllDynamics.ka1_meso;
    kd_all_meso = AllDynamics.kd_meso;
    Na_all_meso = AllDynamics.Na_meso;
    Nd_all_meso = AllDynamics.Nd_meso;
    fa_all_meso = AllDynamics.fa_meso;
    fd_all_meso = AllDynamics.fd_meso;
    tau_d1_all_meso = AllDynamics.tau_d1_meso;
    tau_a_all_meso = AllDynamics.tau_a_meso;
    tau_d_all_meso = AllDynamics.tau_d_meso;

    if AppendScalingData
        ka_all_thry = AllDynamics.ka_thry;
        ka1_all_thry = AllDynamics.ka1_thry;
        kd_all_thry = AllDynamics.kd_thry;
        Na_all_thry = AllDynamics.Na_thry;
        Nd_all_thry = AllDynamics.Nd_thry;
        fa_all_thry = AllDynamics.fa_thry;
        fd_all_thry = AllDynamics.fd_thry;
        tau_d1_all_thry = AllDynamics.tau_d1_thry;
        tau_a_all_thry = AllDynamics.tau_a_thry;
        tau_d_all_thry = AllDynamics.tau_d_thry;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DefineAssembledFileNames(ka_in)

global OutputDir EnsembleDynamicsFileName AllDynamicsFileName...
    FittedDynamicsFileName R2DynamicsFileName R2ThryFileName

EnsembleDynamicsFileName = [OutputDir,'EnsembleDynamics.ka-',...
    num2str(ka_in,'%.2e'),'.m'];
AllDynamicsFileName = [OutputDir,'AllDynamicswrtTime.ka-',...
    num2str(ka_in,'%.2e'),'.m'];
R2DynamicsFileName = [OutputDir,'R2Dynamics.ka-',...
    num2str(ka_in,'%.2e'),'.m'];
FittedDynamicsFileName = [OutputDir,'FittedDynamics.ka-',...
    num2str(ka_in,'%.2e'),'.m'];
R2ThryFileName = [OutputDir,'R2ScalingDynamics.ka-',...
    num2str(ka_in,'%.2e'),'.m'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitializeAllDynamics(time_bead,time_meso,time_thry,NoSep)

global ka_all_meso kd_all_meso...
    Na_all_meso Nd_all_meso...
    fa_all_meso fd_all_meso...
    tau_d1_all_meso tau_a_all_meso tau_d_all_meso...
    ka_all_bead kd_all_bead...
    Na_all_bead Nd_all_bead...
    fa_all_bead fd_all_bead...
    tau_d1_all_bead tau_a_all_bead tau_d_all_bead ...
    ka_all_thry kd_all_thry...
    Na_all_thry Nd_all_thry...
    fa_all_thry fd_all_thry...
    tau_d1_all_thry tau_a_all_thry tau_d_all_thry ...
    AppendScalingData

ka_all_bead = zeros(size(time_bead,1),NoSep);
kd_all_bead = zeros(size(time_bead,1),NoSep);
Na_all_bead = zeros(size(time_bead,1),NoSep);
Nd_all_bead = zeros(size(time_bead,1),NoSep);
fa_all_bead = zeros(size(time_bead,1),NoSep);
fd_all_bead = zeros(size(time_bead,1),NoSep);

tau_d1_all_bead = zeros(5e3,NoSep);
tau_a_all_bead = zeros(5e3,NoSep);
tau_d_all_bead = zeros(5e3,NoSep);


ka_all_meso = zeros(size(time_meso,1),NoSep);
kd_all_meso = zeros(size(time_meso,1),NoSep);
Na_all_meso = zeros(size(time_meso,1),NoSep);
Nd_all_meso = zeros(size(time_meso,1),NoSep);
fa_all_meso = zeros(size(time_meso,1),NoSep);
fd_all_meso = zeros(size(time_meso,1),NoSep);

tau_d1_all_meso = zeros(5e3,NoSep);
tau_a_all_meso = zeros(5e3,NoSep);
tau_d_all_meso = zeros(5e3,NoSep);

if AppendScalingData
    ka_all_thry = zeros(size(time_thry,1),NoSep);
    kd_all_thry = zeros(size(time_thry,1),NoSep);
    Na_all_thry = zeros(size(time_thry,1),NoSep);
    Nd_all_thry = zeros(size(time_thry,1),NoSep);
    fa_all_thry = zeros(size(time_thry,1),NoSep);
    fd_all_thry = zeros(size(time_thry,1),NoSep);

    tau_d1_all_thry = zeros(5e3,NoSep);
    tau_a_all_thry = zeros(5e3,NoSep);
    tau_d_all_thry = zeros(5e3,NoSep);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitializeEnsembleDynamics(N_Kuhns,NoSep)

global ka_ens_meso kd_ens_meso kaSE_ens_meso kdSE_ens_meso...
    Na_ens_meso Nd_ens_meso NaSE_ens_meso NdSE_ens_meso...
    fa_ens_meso fd_ens_meso faSE_ens_meso fdSE_ens_meso...
    ka_ens_bead kd_ens_bead kaSE_ens_bead kdSE_ens_bead...
    Na_ens_bead Nd_ens_bead NaSE_ens_bead NdSE_ens_bead...
    fa_ens_bead fd_ens_bead faSE_ens_bead fdSE_ens_bead...
    tau_d1_ens_bead tau_d1SE_ens_bead...
    tau_a_ens_bead tau_aSE_ens_bead...
    tau_d_ens_bead tau_dSE_ens_bead...
    tau_d1_ens_meso tau_d1SE_ens_meso...
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

tau_d1_ens_bead = zeros(NoSep,length(N_Kuhns));
tau_d1SE_ens_bead = zeros(NoSep,length(N_Kuhns));
tau_a_ens_bead = zeros(NoSep,length(N_Kuhns));
tau_aSE_ens_bead = zeros(NoSep,length(N_Kuhns));
tau_d_ens_bead = zeros(NoSep,length(N_Kuhns));
tau_dSE_ens_bead = zeros(NoSep,length(N_Kuhns));

tau_d1_ens_meso = zeros(NoSep,length(N_Kuhns));
tau_d1SE_ens_meso = zeros(NoSep,length(N_Kuhns));
tau_a_ens_meso = zeros(NoSep,length(N_Kuhns));
tau_aSE_ens_meso = zeros(NoSep,length(N_Kuhns));
tau_d_ens_meso = zeros(NoSep,length(N_Kuhns));
tau_dSE_ens_meso = zeros(NoSep,length(N_Kuhns));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time,ka,ka1,kd,Na,Nd,fa,fd,tau_d1,tau_a,tau_d,...
    mean_ka,se_ka,mean_ka1,se_ka1,mean_kd,se_kd,...
    mean_Na,se_Na,mean_Nd,se_Nd,...
    mean_fa,se_fa,mean_fd,se_fd,...
    mean_tau_d1,se_tau_d1,...
    mean_tau_a,se_tau_a,...
    mean_tau_d,se_tau_d] =...
    CalculateEnsembleDynamics(BSOM,Np,N,ka_in,kd_in,f0,...
    dt,damp,N_Kuhn,b,Separation)

global BeadSpringOrMeso TimeStretchDataFileName BondKineticsDataFileName...
    timenorm


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
ka1 = BondKinetics.ka1;
kd = BondKinetics.kd;
Na = BondKinetics.Na;
Nd = BondKinetics.Nd;

fa = Na./(Na+Nd);
fd = Nd./(Na+Nd);

SSPcnt = 80; SSIndx = round(SSPcnt/100*size(fa,1));

mean_fa = nanmean(fa(SSIndx:end,:),1);
se_fa = nanstd(fa(SSIndx:end,:),1,1)/sqrt(size(fa(SSIndx:end,:),1));
mean_fd = nanmean(fd(SSIndx:end,:),1);
se_fd = nanstd(fd(SSIndx:end,:),1,1)/sqrt(size(fd(SSIndx:end,:),1));

%% FIND THAT I NEED A BETTER WAY TO COMPUTE THE MEAN KA AND KD - FIRST MUST SMOOTH SINCE THE DISTRIBUTION IS SKEWED

mean_ka = nanmean(ka,1)/timenorm; se_ka = nanstd(ka,1,1)/sqrt(size(ka,1))/timenorm;
mean_ka1 = nanmean(ka1,1)/timenorm; se_ka1 = nanstd(ka1,1,1)/sqrt(size(ka1,1))/timenorm;
mean_kd = nanmean(kd,1)/timenorm; se_kd = nanstd(kd,1,1)/sqrt(size(kd,1))/timenorm;
mean_Na = nanmean(Na,1); se_Na = nanstd(Na,1,1)/sqrt(size(Na,1));
mean_Nd = nanmean(Nd,1); se_Nd = nanstd(Nd,1,1)/sqrt(size(Nd,1));

% Attached and detached bond lifetimes
tau_d1 = BondKinetics.FirstAttachment;
tau_a = BondKinetics.AttachedLifetimes;
tau_d = BondKinetics.DetachedLifetimes;

mean_tau_d1 = mean(tau_d1);
se_tau_d1 = std(tau_d1)/sqrt(length(tau_d1));

mean_tau_a = mean(tau_a);
se_tau_a = std(tau_a)/sqrt(length(tau_a));

mean_tau_d = mean(tau_d);
se_tau_d = std(tau_d)/sqrt(length(tau_d));

check_ka_kd = 0;
if check_ka_kd==1
    figure(1); clf; hold on
    x = time;
    smooth_ka = smooth(ka,1e3);
    mean_smooth_ka = mean(smooth_ka(round(length(x)/2):end));
    scatter(x,ka)
    scatter(x,smooth_ka)

    plot(x,timenorm*mean_ka*ones(size(x)),'k--')
    plot(x,mean_smooth_ka*ones(size(x)),'g--')
    ylim([0 ka_in/7])

    figure(2); clf; hold on
    x = time;
    smooth_kd = smooth(kd,1e3);
    mean_smooth_kd = mean(smooth_kd(round(length(x)/2):end));
    scatter(x,kd)
    scatter(x,smooth_kd)

    plot(x,timenorm*mean_kd*ones(size(x)),'k--')
    plot(x,mean_smooth_kd*ones(size(x)),'g--')
    plot(x,kd_in*ones(size(x)),'c--')
    ylim([0 kd_in*1.5])
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,delta_A,R2] = FitThePrefactor(x,y,A_nom,N_Kuhn,b,eaOverkbT)

% global NoFitAttempts
find_conf_int = 1;
[A,R2,delta_A] = RunPrefactorFittingLoop(x,y,A_nom,N_Kuhn,b,eaOverkbT,...
    find_conf_int);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,R2,delta_A] = RunPrefactorFittingLoop(x,y,A_nom,N_Kuhn,b,...
    eaOverkbT,find_conf_int)

NoFitAttempts = 50;
ReductionFactor = 0.98;
NoPts = 21;
R2 = 0;
ct = 0;
lambda = x/sqrt(N_Kuhn)/b;

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

        ft = A_all(i)*(3/2/pi/N_Kuhn)^(3/2)*exp(-eaOverkbT-3/2*lambda.^2);

        RSS = sum((y-ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_all(i) = 1-RSS/TSS;
    end
    diff = abs(R2_all-1);
    indx = find(diff==min(diff),1,'first');

    A_nom = A_all(indx);
    R2 = R2_all(indx);

    if ~mod(ct,20)
        figure(100); clf; hold on
        scatter(x,y,'k','filled')

        yplot = A_nom*(3/2/pi/N_Kuhn)^(3/2)*exp(-eaOverkbT-3/2*lambda.^2);

        plot(x,yplot,'k--')
    end
    if ct>NoFitAttempts
        break;
    end
end
A = A_nom;
figure(100); close
close(wb5)


if find_conf_int==1 %Perturb each parameter to find range in which chi2<=1
    R2_temp = R2;
    A_temp = A;
    dA = 0.01*A;
    while R2_temp>0.05
        A_temp = A_temp+dA;
        y_ft = A_temp*(3/2/pi/N_Kuhn)^(3/2)*exp(-eaOverkbT-3/2*lambda.^2);
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
        y_ft = A_temp*(3/2/pi/N_Kuhn)^(3/2)*exp(-eaOverkbT-3/2*lambda.^2);
        RSS = sum((y-y_ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_temp = 1-RSS/TSS;
    end
    delta_A_neg = abs(A_temp-A);

    delta_A = (delta_A_neg + delta_A_pos)/2;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,k,R2,delta_A,delta_k] = FitExp(A_nom,k_nom,x,y)

global NoFitAttempts

NoPts = 21;
R2 = 0;
ct = 0;
wb5 = waitbar(0,'Fitting Exponential...');
figure(100)
while R2<0.995
    ct = ct+1;
    waitbar(ct/NoFitAttempts,wb5,'Fitting Exponential...')

    if ct==1
        sig_A = 1/8*A_nom;
        sig_k = 1/8*k_nom;
    else
        sig_A = sig_A*0.95;
        sig_k = sig_k*0.95;
    end
    A_rng = [(linspace(A_nom-sig_A,A_nom+sig_A,NoPts))';A_nom];
    A_rng = unique(A_rng); 
    k_rng = [(linspace(k_nom-sig_k,k_nom+sig_k,NoPts))';k_nom];
    k_rng = unique(k_rng);

    Perms = zeros(length(A_rng)*length(k_rng),1);
    i = 0;
    for mi=1:length(A_rng)
        for bi=1:length(k_rng)
            i = i+1;
            Perms(i,1) = A_rng(mi);
            Perms(i,2) = k_rng(bi);
        end
    end
    R2_all = zeros(size(Perms,1),1);
    A_all = zeros(size(Perms,1),1);
    k_all = zeros(size(Perms,1),1);

    for i=1:size(Perms,1)
        A_all(i) = Perms(i,1);
        k_all(i) = Perms(i,2);

        ft = A_all(i)*exp(x*k_all(i));

        RSS = sum((y-ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_all(i) = 1-RSS/TSS;
    end
    diff = abs(R2_all-1);
    indx = find(diff==min(diff),1,'first');

    A_nom = A_all(indx);
    k_nom = k_all(indx);
    R2 = R2_all(indx);

    if ~mod(ct,20)
        figure(100); clf; hold on
        scatter(x,y,'k','filled')
        plot(x,A_nom*exp(x*k_nom),'k--')
    end
    if ct>NoFitAttempts
        break;
    end
end
A = A_nom;
k = k_nom;
figure(100); close
close(wb5)

find_conf_int = 1;
if find_conf_int==1 %Perturb each parameter to find range in which chi2<=1
    R2_temp = R2;
    A_temp = A;
    dA = 0.01*A;
    while R2_temp>0.05
        A_temp = A_temp+dA;
        y_ft = A_temp*exp(x*k);
        RSS = sum((y-y_ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_temp = 1-RSS/TSS;
    end
    delta_A_pos = abs(A_temp-A);
    R2_temp = R2;
    k_temp = k;
    dtau = 0.01*k;
    while R2_temp>0.05
        k_temp = k_temp+dtau;
        y_ft = A*exp(x*k_temp);
        RSS = sum((y-y_ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_temp = 1-RSS/TSS;
    end
    delta_k_pos = abs(k_temp-k);

    R2_temp = R2;
    A_temp = A;
    dA = 0.01*A;
    while R2_temp>0.05
        A_temp = A_temp-dA;
        y_ft = A_temp*exp(x/k);
        RSS = sum((y-y_ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_temp = 1-RSS/TSS;
    end
    delta_A_neg = abs(A_temp-A);
    R2_temp = R2;
    k_temp = k;
    dtau = 0.01*k;
    while R2_temp>0.05
        k_temp = k_temp-dtau;
        y_ft = A*exp(x/k_temp);
        RSS = sum((y-y_ft).^2);
        TSS = sum((y-mean(y)).^2);
        R2_temp = 1-RSS/TSS;
    end
    delta_k_neg = abs(k_temp-k);

    delta_A = (delta_A_neg + delta_A_pos)/2;
    delta_k = (delta_k_neg + delta_k_pos)/2;
end

end
