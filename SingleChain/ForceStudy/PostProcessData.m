function PostProcessData(Package,ToggleDynamics,LC,DC,FC,BSOM,CM,...
    CF,OF,BT)

% Computes and plots stress-strain for every iteration of simulation

global LineWidth TurnOnDynamics FontSize...
    LengthConversion DamperConversion DataSize BeadSpringOrMeso...
    CompareModels CurrentFolder OutputFolder PlotForFigures ForceConversion...
    BondType

LineWidth = 1.5;
FontSize = 20;
DataSize = 30;
TurnOnDynamics = ToggleDynamics;
LengthConversion = LC;
DamperConversion = DC;
ForceConversion = FC;
BeadSpringOrMeso = BSOM;
CompareModels = CM;
CurrentFolder = CF;
OutputFolder = OF;
PlotForFigures = 0; %1 for version of code from which figs were plotted
BondType = BT;

%% Unpack Swept Input Parameters
Nps = unique(Package(:,2));        %Number of molecules
Ds = unique(Package(:,3));         %Diffusion coefficient [m2/s]
N_Kuhns = unique(Package(:,4));    %Number of Kuhn segments in chain
stiffnesses = unique(Package(:,5)); %Stiffness of single harmonic bond

kbTs = unique(Package(:,6));       %Thermal energy
bs = unique(Package(:,7));         %Kuhn length
dts = unique(Package(:,8));        %timestep

PlotForceVersusDist(Nps,Ds,N_Kuhns,stiffnesses,...
    kbTs,bs,dts);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotForceVersusDist(Nps,Ds,N_Kuhns,stiffnesses,...
    kbTs,bs,dts)

global FontSize LengthConversion ForceDataFileName EndToEndDataFileName...
    AlignmentDataFileName TimeStretchDataFileName ForceConversion...
    BeadSpringOrMeso

if ~isfolder('Output Plots')
    mkdir('Output Plots')
end

Np = unique(Nps);
D = unique(Ds);
kbT = unique(kbTs);
b = unique(bs);
dt = unique(dts);

color1 = [0 0.85 0];
color2 = [0 0 0];
% color2 = [0 0.85 0];
colors_nk = DefinePlotColors(color1,color2,length(N_Kuhns));
figure(4); clf; hold on
for nk=1:length(N_Kuhns)
    N_Kuhn = N_Kuhns(nk);
    color_nk = colors_nk(nk,:);

    figure(1); clf; hold on
    figure(2); clf; hold on
    figure(3); clf; hold on

    color1 = [0 0 0];
    color2 = [1 0 0];
    colors_stf = DefinePlotColors(color1,color2,length(stiffnesses));

    marker_styles = {'o','^','sq','>','d','v'};
    legend_entries_stf = [];
    for stf = 1:length(stiffnesses)
        stiffness = stiffnesses(stf);
        color_stf = colors_stf(stf,:);
        marker_style = marker_styles{stf};

        % Import the bead-spring data
        BeadSpringOrMeso = 0;
        InputScript(1,Np,D,N_Kuhn,stiffness,kbT,b,dt);
        SetDirAndFileNames;

        dat_bs = load(TimeStretchDataFileName,'-mat');
        stretch_bs = dat_bs.stretch;

        stretches_bs = unique(round(stretch_bs,2));

        % Note, distances are given in normalized units
        dat_bs = load(EndToEndDataFileName,'-mat');
        rx_bs = dat_bs.rx;
        ry_bs = dat_bs.ry;
        rz_bs = dat_bs.rz;
        r_bs = (rx_bs.^2 + ry_bs.^2 + rz_bs.^2).^0.5;

        dat_bs = load(AlignmentDataFileName,'-mat');
        rxr11_bs = dat_bs.rxr11;
        rxr22_bs = dat_bs.rxr22;
        rxr33_bs = dat_bs.rxr33;
        rxr12_bs = dat_bs.rxr12;
        rxr23_bs = dat_bs.rxr23;
        rxr31_bs = dat_bs.rxr31;

        % Note, forces are given in SI units
        dat_bs = load(ForceDataFileName,'-mat');
        fx_bs = dat_bs.fx;
        fy_bs = dat_bs.fy;
        fz_bs = dat_bs.fz;
        
        % Histogram of end-to-ends
        figure(1) 
        hoa(stf) = histogram(r_bs/b);
        bins = (0.5:0.01:1.5);
        hoa(stf).BinEdges = bins;
        hoa(stf).EdgeColor = 'k';
        hoa(stf).FaceColor = color_stf;
        hoa(stf).FaceAlpha = 1;
        hoa(stf).Normalization = 'probability';

        % Import the mesoscale data
        BeadSpringOrMeso = 1;
        InputScript(1,Np,D,N_Kuhn,stiffness,kbT,b,dt);
        SetDirAndFileNames;

        dat_ms = load(TimeStretchDataFileName,'-mat');
        stretch_ms = dat_ms.stretch;

        stretches_ms = unique(round(stretch_ms,2));

        % Note, distances are given in normalized units
        dat_ms = load(EndToEndDataFileName,'-mat');
        rx_ms = dat_ms.rx;
        ry_ms = dat_ms.ry;
        rz_ms = dat_ms.rz;
        r_ms = (rx_ms.^2 + ry_ms.^2 + rz_ms.^2).^0.5;

        dat_ms = load(AlignmentDataFileName,'-mat');
        rxr11_ms = dat_ms.rxr11;
        rxr22_ms = dat_ms.rxr22;
        rxr33_ms = dat_ms.rxr33;
        rxr12_ms = dat_ms.rxr12;
        rxr23_ms = dat_ms.rxr23;
        rxr31_ms = dat_ms.rxr31;

        % Note, forces are given in SI units
        dat_ms = load(ForceDataFileName,'-mat');
        fx_ms = dat_ms.fx;
        fy_ms = dat_ms.fy;
        fz_ms = dat_ms.fz;
        
        % Plot accross stiffnesses
        figure(101); clf; hold on

        % Initialize the bead-spring data
        mean_rxr11_bs = zeros(length(stretches_bs),1);
        mean_rxr22_bs = zeros(length(stretches_bs),1);
        mean_rxr33_bs = zeros(length(stretches_bs),1);
        mean_rxr12_bs = zeros(length(stretches_bs),1);
        mean_rxr23_bs = zeros(length(stretches_bs),1);
        mean_rxr31_bs = zeros(length(stretches_bs),1);
        se_rxr11_bs = zeros(length(stretches_bs),1);
        se_rxr22_bs = zeros(length(stretches_bs),1);
        se_rxr33_bs = zeros(length(stretches_bs),1);
        se_rxr12_bs = zeros(length(stretches_bs),1);
        se_rxr23_bs = zeros(length(stretches_bs),1);
        se_rxr31_bs = zeros(length(stretches_bs),1);

        mean_fx_bs = zeros(length(stretches_bs),1);
        mean_fy_bs = zeros(length(stretches_bs),1);
        mean_fz_bs = zeros(length(stretches_bs),1);
        se_fx_bs = zeros(length(stretches_bs),1);
        se_fy_bs = zeros(length(stretches_bs),1);
        se_fz_bs = zeros(length(stretches_bs),1);

        mean_end_fx_bs = zeros(length(stretches_bs),1);
        mean_end_fy_bs = zeros(length(stretches_bs),1);
        mean_end_fz_bs = zeros(length(stretches_bs),1);
        se_end_fx_bs = zeros(length(stretches_bs),1);
        se_end_fy_bs = zeros(length(stretches_bs),1);
        se_end_fz_bs = zeros(length(stretches_bs),1);

 
        % Initialize the mesoscale data
        mean_rxr11_ms = zeros(length(stretches_ms),1);
        mean_rxr22_ms = zeros(length(stretches_ms),1);
        mean_rxr33_ms = zeros(length(stretches_ms),1);
        mean_rxr12_ms = zeros(length(stretches_ms),1);
        mean_rxr23_ms = zeros(length(stretches_ms),1);
        mean_rxr31_ms = zeros(length(stretches_ms),1);
        se_rxr11_ms = zeros(length(stretches_ms),1);
        se_rxr22_ms = zeros(length(stretches_ms),1);
        se_rxr33_ms = zeros(length(stretches_ms),1);
        se_rxr12_ms = zeros(length(stretches_ms),1);
        se_rxr23_ms = zeros(length(stretches_ms),1);
        se_rxr31_ms = zeros(length(stretches_ms),1);

        mean_fx_ms = zeros(length(stretches_ms),1);
        mean_fy_ms = zeros(length(stretches_ms),1);
        mean_fz_ms = zeros(length(stretches_ms),1);
        se_fx_ms = zeros(length(stretches_ms),1);
        se_fy_ms = zeros(length(stretches_ms),1);
        se_fz_ms = zeros(length(stretches_ms),1);

        mean_end_fx_ms = zeros(length(stretches_ms),1);
        mean_end_fy_ms = zeros(length(stretches_ms),1);
        mean_end_fz_ms = zeros(length(stretches_ms),1);
        se_end_fx_ms = zeros(length(stretches_ms),1);
        se_end_fy_ms = zeros(length(stretches_ms),1);
        se_end_fz_ms = zeros(length(stretches_ms),1);


        legend_entries_str = []; lgnd_ct = 0;
        color1 = [0 1 1];
        color2 = [1 0 1];
        colors_str = DefinePlotColors(color1,color2,length(stretches_bs));    
        legend_entries_str = []; clear h
        for str=1:length(stretches_bs)
            stretch_temp = stretches_bs(str);
            indices = find(abs(stretch_temp-stretch_bs)<0.01);
            color_str = colors_str(str,:);

            % Compile the bead-spring data
            rx_temp_bs = rx_bs(indices,:,:);
            ry_temp_bs = ry_bs(indices,:,:);
            rz_temp_bs = rz_bs(indices,:,:);
            r_temp_bs = (rx_temp_bs.^2 + ry_temp_bs.^2 + rz_temp_bs.^2).^0.5;

            rxr11_temp_bs = rxr11_bs(indices,:);
            rxr22_temp_bs = rxr22_bs(indices,:);
            rxr33_temp_bs = rxr33_bs(indices,:);
            rxr12_temp_bs = rxr12_bs(indices,:);
            rxr23_temp_bs = rxr23_bs(indices,:);
            rxr31_temp_bs = rxr31_bs(indices,:);

            mean_rxr11_bs(str) = mean(mean(rxr11_temp_bs));
            se_rxr11_bs(str) = std(rxr11_temp_bs(:))/sqrt(length(rxr11_temp_bs(:)));
            mean_rxr22_bs(str) = mean(mean(rxr22_temp_bs));
            se_rxr22_bs(str) = std(rxr22_temp_bs(:))/sqrt(length(rxr22_temp_bs(:)));
            mean_rxr33_bs(str) = mean(mean(rxr33_temp_bs));
            se_rxr33_bs(str) = std(rxr33_temp_bs(:))/sqrt(length(rxr33_temp_bs(:)));

            mean_rxr12_bs(str) = mean(mean(rxr12_temp_bs));
            se_rxr12_bs(str) = std(rxr12_temp_bs(:))/sqrt(length(rxr12_temp_bs(:)));
            mean_rxr23_bs(str) = mean(mean(rxr23_temp_bs));
            se_rxr23_bs(str) = std(rxr23_temp_bs(:))/sqrt(length(rxr23_temp_bs(:)));
            mean_rxr31_bs(str) = mean(mean(rxr31_temp_bs));
            se_rxr31_bs(str) = std(rxr31_temp_bs(:))/sqrt(length(rxr31_temp_bs(:)));

            fx_temp_bs = -(fx_bs(indices,:,:));
            fy_temp_bs = -(fz_bs(indices,:,:));
            fz_temp_bs = -(fy_bs(indices,:,:));

            fx_tot_bs = squeeze(mean(fx_temp_bs,2));
            fy_tot_bs = squeeze(mean(fy_temp_bs,2));
            fz_tot_bs = squeeze(mean(fz_temp_bs,2));

            rng = (round(size(fx_tot_bs,1)/2):size(fx_tot_bs,1))';
            fx_mean_bs = mean(fx_tot_bs(rng,:),1);
            fy_mean_bs = mean(fy_tot_bs(rng,:),1);
            fz_mean_bs = mean(fz_tot_bs(rng,:),1);

            fx_mean_bs(isnan(fx_mean_bs)) = [];
            fy_mean_bs(isnan(fy_mean_bs)) = [];
            fz_mean_bs(isnan(fz_mean_bs)) = [];

            mean_fx_bs(str) = mean(fx_mean_bs);
            mean_fy_bs(str) = mean(fy_mean_bs);
            mean_fz_bs(str) = mean(fz_mean_bs);

            se_fx_bs(str) = std(fx_mean_bs,1)/sqrt(length(fx_mean_bs));
            se_fy_bs(str) = std(fy_mean_bs,1)/sqrt(length(fy_mean_bs));
            se_fz_bs(str) = std(fz_mean_bs,1)/sqrt(length(fz_mean_bs));
           
            fx_tot_end_bs = squeeze(fx_temp_bs(:,end,:));
            fy_tot_end_bs = squeeze(fy_temp_bs(:,end,:));
            fz_tot_end_bs = squeeze(fz_temp_bs(:,end,:));

            rng = (round(size(fx_tot_end_bs,1)/2):size(fx_tot_end_bs,1))';
            fx_mean_end_bs = mean(fx_tot_end_bs(rng,:),1);
            fy_mean_end_bs = mean(fy_tot_end_bs(rng,:),1);
            fz_mean_end_bs = mean(fz_tot_end_bs(rng,:),1);

            fx_mean_end_bs(isnan(fx_mean_end_bs)) = [];
            fy_mean_end_bs(isnan(fy_mean_end_bs)) = [];
            fz_mean_end_bs(isnan(fz_mean_end_bs)) = [];
            
            mean_end_fx_bs(str) = mean(fx_mean_end_bs);
            mean_end_fy_bs(str) = mean(fy_mean_end_bs);
            mean_end_fz_bs(str) = mean(fz_mean_end_bs);

            se_end_fx_bs(str) = std(fx_mean_end_bs)/sqrt(length(fx_mean_end_bs));
            se_end_fy_bs(str) = std(fy_mean_end_bs)/sqrt(length(fy_mean_end_bs));
            se_end_fz_bs(str) = std(fz_mean_end_bs)/sqrt(length(fz_mean_end_bs));
                       

            % Compile the bead-spring data
            rx_temp_ms = rx_ms(indices,:,:);
            ry_temp_ms = ry_ms(indices,:,:);
            rz_temp_ms = rz_ms(indices,:,:);
            r_temp_ms = (rx_temp_ms.^2 + ry_temp_ms.^2 + rz_temp_ms.^2).^0.5;

            rxr11_temp_ms = rxr11_ms(indices,:);
            rxr22_temp_ms = rxr22_ms(indices,:);
            rxr33_temp_ms = rxr33_ms(indices,:);
            rxr12_temp_ms = rxr12_ms(indices,:);
            rxr23_temp_ms = rxr23_ms(indices,:);
            rxr31_temp_ms = rxr31_ms(indices,:);

            mean_rxr11_ms(str) = mean(mean(rxr11_temp_ms));
            se_rxr11_ms(str) = std(rxr11_temp_ms(:))/sqrt(length(rxr11_temp_ms(:)));
            mean_rxr22_ms(str) = mean(mean(rxr22_temp_ms));
            se_rxr22_ms(str) = std(rxr22_temp_ms(:))/sqrt(length(rxr22_temp_ms(:)));
            mean_rxr33_ms(str) = mean(mean(rxr33_temp_ms));
            se_rxr33_ms(str) = std(rxr33_temp_ms(:))/sqrt(length(rxr33_temp_ms(:)));

            mean_rxr12_ms(str) = mean(mean(rxr12_temp_ms));
            se_rxr12_ms(str) = std(rxr12_temp_ms(:))/sqrt(length(rxr12_temp_ms(:)));
            mean_rxr23_ms(str) = mean(mean(rxr23_temp_ms));
            se_rxr23_ms(str) = std(rxr23_temp_ms(:))/sqrt(length(rxr23_temp_ms(:)));
            mean_rxr31_ms(str) = mean(mean(rxr31_temp_ms));
            se_rxr31_ms(str) = std(rxr31_temp_ms(:))/sqrt(length(rxr31_temp_ms(:)));

            fx_temp_ms = -(fx_ms(indices,:,:));
            fy_temp_ms = -(fz_ms(indices,:,:));
            fz_temp_ms = -(fy_ms(indices,:,:));

            fx_tot_ms = squeeze(mean(fx_temp_ms,2));
            fy_tot_ms = squeeze(mean(fy_temp_ms,2));
            fz_tot_ms = squeeze(mean(fz_temp_ms,2));

            rng = (round(size(fx_tot_ms,1)/2):size(fx_tot_ms,1))';
            fx_mean_ms = mean(fx_tot_ms(rng,:),1);
            fy_mean_ms = mean(fy_tot_ms(rng,:),1);
            fz_mean_ms = mean(fz_tot_ms(rng,:),1);

            fx_mean_ms(isnan(fx_mean_ms)) = [];
            fy_mean_ms(isnan(fy_mean_ms)) = [];
            fz_mean_ms(isnan(fz_mean_ms)) = [];

            mean_fx_ms(str) = mean(fx_mean_ms);
            mean_fy_ms(str) = mean(fy_mean_ms);
            mean_fz_ms(str) = mean(fz_mean_ms);

            se_fx_ms(str) = std(fx_mean_ms,1)/sqrt(length(fx_mean_ms));
            se_fy_ms(str) = std(fy_mean_ms,1)/sqrt(length(fy_mean_ms));
            se_fz_ms(str) = std(fz_mean_ms,1)/sqrt(length(fz_mean_ms));
           
            fx_tot_end_ms = squeeze(fx_temp_ms(:,end,:));
            fy_tot_end_ms = squeeze(fy_temp_ms(:,end,:));
            fz_tot_end_ms = squeeze(fz_temp_ms(:,end,:));

            rng = (round(size(fx_tot_end_ms,1)/2):size(fx_tot_end_ms,1))';
            fx_mean_end_ms = mean(fx_tot_end_ms(rng,:),1);
            fy_mean_end_ms = mean(fy_tot_end_ms(rng,:),1);
            fz_mean_end_ms = mean(fz_tot_end_ms(rng,:),1);

            fx_mean_end_ms(isnan(fx_mean_end_ms)) = [];
            fy_mean_end_ms(isnan(fy_mean_end_ms)) = [];
            fz_mean_end_ms(isnan(fz_mean_end_ms)) = [];
            
            mean_end_fx_ms(str) = mean(fx_mean_end_ms);
            mean_end_fy_ms(str) = mean(fy_mean_end_ms);
            mean_end_fz_ms(str) = mean(fz_mean_end_ms);

            se_end_fx_ms(str) = std(fx_mean_end_ms)/sqrt(length(fx_mean_end_ms));
            se_end_fy_ms(str) = std(fy_mean_end_ms)/sqrt(length(fy_mean_end_ms));
            se_end_fz_ms(str) = std(fz_mean_end_ms)/sqrt(length(fz_mean_end_ms));
            

            % Histogram of end-to-ends
            figure(101)
            if ~mod(str+1,2)
                lgnd_ct = lgnd_ct+1;
                h(lgnd_ct) = histogram(r_temp_bs/b);
                h(lgnd_ct).BinEdges = bins;
                h(lgnd_ct).EdgeColor = 'k';
                h(lgnd_ct).FaceColor = color_str;
                h(lgnd_ct).Normalization = 'probability';
                h(lgnd_ct).FaceAlpha = 1;
                legend_entries_str{lgnd_ct} = ...
                    ['$\lambda^* = $ ',num2str(stretch_temp/sqrt(N_Kuhn),'%.2f')];
            else
                temp_h = histogram(r_temp_bs/b);
                temp_h.BinEdges = bins;
                temp_h.EdgeColor = 'k';
                temp_h.FaceColor = color_str;
                temp_h.Normalization = 'probability';
                temp_h.FaceAlpha = 1;
            end           
        end
        figure(101)
        plot([1 1],[0 0.12],'k--')
        xlim([0.75 1.25])
        ylim([0 Inf])
        set(gca,'FontSize',FontSize/1.75)
        set(gcf,'color','w')
        pbaspect([1 1 1])
        l = legend(h,legend_entries_str);
        l.FontSize = FontSize/2;
        l.Location = 'Northeast';
        l.Interpreter = 'latex';
        xlabel('$\bar{r}_{\alpha \beta}/b$','FontSize',FontSize,'Interpreter','latex')
        ylabel('$p$','FontSize',FontSize,'Interpreter','latex')
        title(['$K=$',num2str(stiffness),' $k_b T/b^2$'],'FontSize',FontSize/2,'Interpreter','latex')
        saveas(gcf,['Output Plots/Histos_stretch.stiffness',num2str(stiffness),'.png'])
        saveas(gcf,['Output Plots/Histos_stretch.stiffness',num2str(stiffness),'.fig'])

        figure(2)
        subplot(2,2,1); hold on
        e2(stf) = errorbar(stretches_bs,mean_rxr11_bs,se_rxr11_bs);
        e2(stf).Marker = marker_style;
        e2(stf).MarkerFaceColor = color_stf;
        e2(stf).MarkerEdgeColor = 'k';
        e2(stf).Color = 'none';
        e2(stf).LineStyle = 'none';
        set(gca,'FontSize',FontSize/1.5)
        set(gcf,'color','w')
        ylabel('$\langle \hat r_1 \hat r_1 \rangle$','FontSize',FontSize,'Interpreter','latex')
        ylim([0 1])

        subplot(2,2,2); hold on
        e = errorbar(stretches_bs,mean_rxr22_bs,se_rxr22_bs);
        e.Marker = marker_style;
        e.MarkerFaceColor = color_stf;
        e.MarkerEdgeColor = 'k';
        e.Color = 'none';
        e.LineStyle = 'none';
        set(gca,'FontSize',FontSize/1.5)
        set(gcf,'color','w')
        ylabel('$\langle \hat r_2 \hat r_2 \rangle$','FontSize',FontSize,'Interpreter','latex')
        ylim([0 1])

        subplot(2,2,3); hold on
        e = errorbar(stretches_bs,mean_rxr33_bs,se_rxr33_bs);
        e.Marker = marker_style;
        e.MarkerFaceColor = color_stf;
        e.MarkerEdgeColor = 'k';
        e.Color = 'none';
        e.LineStyle = 'none';
        set(gca,'FontSize',FontSize/1.5)
        set(gcf,'color','w')
        ylabel('$\langle \hat r_3 \hat r_3 \rangle$','FontSize',FontSize,'Interpreter','latex')
        xlabel('$\lambda$','FontSize',FontSize,'Interpreter','latex')
        ylim([0 1])

        subplot(2,2,4); hold on
        e = errorbar(stretches_bs,mean_rxr12_bs,se_rxr12_bs);
        e.Marker = marker_style;
        e.MarkerFaceColor = color_stf;
        e.MarkerEdgeColor = 'k';
        e.Color = 'none';
        e.LineStyle = 'none';
        set(gca,'FontSize',FontSize/1.5)
        set(gcf,'color','w')
        ylabel('$\langle \hat r_1 \hat r_2 \rangle$','FontSize',FontSize,'Interpreter','latex')
        xlabel('$\lambda$','FontSize',FontSize,'Interpreter','latex')
        ylim([-0.5 0.5])

        figure(3)
        if stf==1
            % plot entropic tension from ideal spring model
            npts = 50;
            lambda = (linspace(0,max(stretches_bs),npts))';
            r_bs = lambda*sqrt(N_Kuhn)*b*LengthConversion;

            f_lin = 3*kbT*r_bs/(N_Kuhn*(b*LengthConversion)^2);
            plot(lambda,f_lin/ForceConversion,'k--','LineWidth',1.5)

            f_lang = kbT*lambda/(sqrt(N_Kuhn)*b*LengthConversion).*...
                ((lambda.^2-3*N_Kuhn)./(lambda.^2-N_Kuhn));
            plot(lambda,f_lang/ForceConversion,'k-.','LineWidth',1.5)
        end
        ebs3(stf) = errorbar(stretches_bs,mean_fx_bs/ForceConversion,se_fx_bs/ForceConversion); % Normalized by kbT/b
        ebs3(stf).Marker = marker_style;
        ebs3(stf).MarkerFaceColor = color_stf;
        ebs3(stf).MarkerEdgeColor = 'k';
        ebs3(stf).Color = 'k';
        ebs3(stf).LineStyle = 'none';

        ems3(stf) = errorbar(stretches_ms,mean_fx_ms/ForceConversion,se_fx_ms/ForceConversion); % Normalized by kbT/b
        ems3(stf).Marker = marker_style;
        ems3(stf).MarkerFaceColor = 'none';
        ems3(stf).MarkerEdgeColor = color_stf;
        ems3(stf).Color = color_stf;
        ems3(stf).LineStyle = 'none';
        
        legend_entries_stf{stf} = ['$K=$ ',num2str(stiffness),' $k_b T/b^2$'];

        set(gca,'FontSize',FontSize/1.5)
        set(gcf,'color','w')
        xlabel('$\lambda$','FontSize',FontSize,'Interpreter','latex')
        ylabel('$\bar{f}/(k_b T/b)$','FontSize',FontSize,'Interpreter','latex')
        ylim([0 50])

    end

    figure(1)
    plot([1 1],[0 0.12],'k--')
    xlim([0.75 1.25])
    ylim([0 Inf])
    set(gca,'FontSize',FontSize/1.75)
    set(gcf,'color','w')
    pbaspect([1 1 1])
    l = legend(hoa,legend_entries_stf);
    l.FontSize = FontSize/2;
    l.Location = 'Northeast';
    l.Interpreter = 'latex';
    xlabel('$\bar{r}_{\alpha \beta}/b$','FontSize',FontSize,'Interpreter','latex')
    ylabel('$p$','FontSize',FontSize,'Interpreter','latex')
    title(['$N=$',num2str(N_Kuhn)],'FontSize',FontSize/2,'Interpreter','latex')
    saveas(gcf,['Output Plots/Histos_chain_length.N',num2str(N_Kuhn),'.png'])
    saveas(gcf,['Output Plots/Histos_chain_length.N',num2str(N_Kuhn),'.fig'])


    figure(2)
    subplot(2,2,1)
    l = legend(e2,legend_entries_stf);
    l.FontSize = FontSize/4;
    l.Location = 'Northwest';
    l.Interpreter = 'latex';
    saveas(gcf,['Output Plots/Chain_alignment.N',num2str(N_Kuhn),'.png'])
    saveas(gcf,['Output Plots/Chain_alignment.N',num2str(N_Kuhn),'.fig'])
    
    
    figure(3)
    l = legend(ebs3,legend_entries_stf);
    l.FontSize = FontSize/1.75;
    l.Location = 'Northwest';
    l.Interpreter = 'latex';
    pbaspect([1 1 1])
    xlim([0 3.5])
    ylim([0 35])
    title(['$N=$',num2str(N_Kuhn)],'FontSize',FontSize/2,'Interpreter','latex')
    saveas(gcf,['Output Plots/Force_stretch.N',num2str(N_Kuhn),'.png'])
    saveas(gcf,['Output Plots/Force_stretch.N',num2str(N_Kuhn),'.fig'])

    figure(4)
    % plot entropic tension from ideal spring model
    npts = 50;
    lambda = (linspace(0,max(stretches_bs),npts))';

    f_lang = kbT*lambda/(sqrt(N_Kuhn)*b*LengthConversion).*...
        ((lambda.^2-3*N_Kuhn)./(lambda.^2-N_Kuhn));
    p = plot(lambda,f_lang/ForceConversion,'LineWidth',1.5);
    p.Color = color_nk;
    p.LineStyle = '-.';
    
    ebs4(nk) = errorbar(stretches_bs,mean_fx_bs/ForceConversion,se_fx_bs/ForceConversion); % Normalized by kbT/b
    ebs4(nk).Marker = 'o';
    ebs4(nk).MarkerFaceColor = color_nk;
    ebs4(nk).MarkerEdgeColor = 'k';
    ebs4(nk).Color = 'k';
    ebs4(nk).LineStyle = 'none';
    
%     ems4(nk) = errorbar(stretches_ms,-mean_fx_ms/ForceConversion,se_fx_ms/ForceConversion); % Normalized by kbT/b
%     ems4(nk).Marker = 'o';
%     ems4(nk).MarkerFaceColor = 'none';
%     ems4(nk).MarkerEdgeColor = color_nk;
%     ems4(nk).Color = color_nk;
%     ems4(nk).LineStyle = 'none';

    legend_entries_nk{nk} = ['$N=$ ',num2str(N_Kuhn)];

    set(gca,'FontSize',FontSize/1.5)
    set(gcf,'color','w')
    xlabel('$\lambda$','FontSize',FontSize,'Interpreter','latex')
    ylabel('$\bar{f}/(k_b T/b)$','FontSize',FontSize,'Interpreter','latex')
    ylim([0 50])
end

figure(4)
l = legend(ebs4,legend_entries_nk);
l.FontSize = FontSize/1.75;
l.Location = 'Northwest';
l.Interpreter = 'latex';
pbaspect([1 1 1])
xlim([0 3.5])
ylim([0 35])
saveas(gcf,'Output Plots/Force_stretch.png')
saveas(gcf,'Output Plots/Force_stretch.fig')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function colors = DefinePlotColors(color1,color2,ncolors)

rinterp = (linspace(color1(1),color2(1),ncolors))';
ginterp = (linspace(color1(2),color2(2),ncolors))';
binterp = (linspace(color1(3),color2(3),ncolors))';
colors = [rinterp ginterp binterp];

end
