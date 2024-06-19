function PostProcessData(Parameters,Controls,Directories,MT)

% Computes and plots stress-strain for every iteration of simulation

global line_width TurnOnDynamics font_size samples Np Nt eq_time_factors...
    LinearOrLangevin current_folder n_samps n_eqtimes n_pts N...
    length_conversion damper_conversion force_conversion energy_conversion...
    ModelTypes

% Unpack the parameters
package = Parameters.package;
UnpackParameters(package,Parameters);

length_conversion = Parameters.length_conversion;
damper_conversion = Parameters.damper_conversion;
force_conversion = Parameters.force_conversion;
energy_conversion = Parameters.energy_conversion;
ModelTypes = MT;

% Define plot settings
line_width = 1.5;   % default curve line width
font_size = 20;     % default figure font size
n_pts = 500;        % default number of points for interpolation and plotting

% Define controls
TurnOnDynamics = Controls.ToggleDynamics;
LinearOrLangevin = Controls.LinearOrLangevin;

% Define directories
current_folder = Directories.current_folder;

%% Unpack Swept Input Parameters
N = Np*Nt;
n_samps = length(samples);
n_eqtimes = length(eq_time_factors);

% Consolidate the data for all three equlibrium times, interpolate it, and
% package it into one array for both bead-spring and mesoscale models
ConsolidateTheData(Controls.OverrideConsolidate,Controls.RunTheLJCase,Controls);

%% MAKE POLISHED PLOTS
PlotPolishedFigs(Controls);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotPolishedFigs(Controls)

global N_Kuhns phis kds Weissenbergs...
    N_Kuhn b phi kd Weissenberg dt D damp...
    length_conversion damper_conversion n_pts tau0

Override = Controls.OverridePostProcess;
RunTheLJCase = Controls.RunTheLJCase;
RunChainLengthSweep = Controls.RunChainLengthSweep;
RunLoadingRateSweep = Controls.RunLoadingRateSweep;
RunDetachmentRateSweep = Controls.RunDetachmentRateSweep;
RunOscillatory = Controls.RunOscillatory;
RunLargeDeformation = Controls.RunLargeDeformation;

n_pts = 45;
folder_name = 'Polished Plots';
if RunTheLJCase==1
    folder_name = [folder_name,'/LJ Study'];
elseif RunChainLengthSweep==1
    folder_name = [folder_name,'/Length Sweep'];
elseif RunLoadingRateSweep==1
    folder_name = [folder_name,'/Loading Rate Sweep'];
elseif RunDetachmentRateSweep==1
    folder_name = [folder_name,'/Detachment Rate Sweep'];
elseif RunOscillatory==1
    folder_name = [folder_name,'/Frequency Sweep'];
elseif RunLargeDeformation==1
    folder_name = [folder_name,'/Large Deformation'];
else
    folder_name = [folder_name,'/Overall Sweep'];
end
if ~isfolder(folder_name)
    mkdir(folder_name)
end

if RunChainLengthSweep
    for i=1:length(phis)
        phi = phis(i);
        PlotChainLengthSweep(phi,N_Kuhns,kd,Weissenberg,folder_name,n_pts,Controls)
    end
elseif RunLoadingRateSweep
    for i=1:length(phis)
        phi = phis(i);
        PlotLoadingRateSweep(phi,N_Kuhns,kd,Weissenbergs,folder_name,n_pts,Controls)
    end    
elseif RunDetachmentRateSweep
    for i=1:length(phis)
        phi = phis(i);
        PlotDetachmentRateSweep(phi,N_Kuhns,kds,Weissenberg,folder_name,n_pts,Controls)
    end   
elseif RunOscillatory
    PlotStorageLoss(kds,Controls);
elseif RunLargeDeformation
    PlotLargeDeformationResponse(phi,N_Kuhns,kds,Weissenbergs,folder_name,n_pts,Controls)
else
    for i=1:length(N_Kuhns)
     N_Kuhn = N_Kuhns(i);
     for ii = 1:length(kds)
      kd = kds(ii);
      for iii=1:length(Weissenbergs)
        Weissenberg = Weissenbergs(iii);

        [damp,D,tau0,dtFact] = DefineTimeScale(b,length_conversion,damper_conversion,...
            1,phi,N_Kuhn);
        dt = tau0/dtFact;

        if RunTheLJCase
            PlotTheLJStudyOutputs(phis,N_Kuhn,kd,Weissenberg,folder_name,n_pts,Controls)
        else
            PlotTheFullParameterSweep(phis,N_Kuhn,kd,Weissenberg,folder_name,n_pts,Controls)
        end

      end
     end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotStorageLoss(kds,Controls)

global storeage_loss_filename tau0 omega_min omega_max lambda...
    b length_conversion D edStar

color_rng = (linspace(0,1,6))';
color_rng = color_rng([1,4,6]);
colors = [color_rng,zeros(size(color_rng)),zeros(size(color_rng))];

figure(1); clf; hold on
figno = 0;
freq_xover = zeros(size(kds));
for i=1:length(kds)
    color = colors(i,:);
    kd = kds(i);
    figno = figno+1;

    % Import data
    edStar = -log(kd*(b*length_conversion)^2/D);
    DefineFileNames(Controls);
    dat = load(storeage_loss_filename,'-mat');

    E_prime = dat.E_prime_mean;
    E_prime_err = dat.E_prime_se;
    E_prime = E_prime/1e3;
    E_prime_err = E_prime_err/1e3;

    E_dbl_prime = dat.E_dbl_prime_mean;
    E_dbl_prime_err = dat.E_dbl_prime_se;
    E_dbl_prime = E_dbl_prime/1e3;
    E_dbl_prime_err = E_dbl_prime_err/1e3;

    freq = dat.freq_mean;

    delta = dat.delta_mean;
    delta_err = dat.delta_se;

    sig0 = dat.sig0_mean;
    sig0_err = dat.sig0_se;

    figure(figno); clf; hold on

    % plot([kd kd]*tau0,[0.0001 3e5]/1e3,'k--')
    % ylim([0 Inf])
%     plot([kd kd]*tau0,[0.0001 1e6],'k--')

    % plot kd and cross-over frequency
%     p = plot([kd kd]*tau0,[0.0001 1e6]);
    p = plot([kd kd]*tau0,[0.0001 1e3]);
    p.Color = color;
    p.LineStyle = '--';
    p.LineWidth = 1.5;
%     ylim([1e3 1e6])
    ylim([8 3e2])
    yticks([10 100])

    freq_interp = (linspace(min(freq),max(freq),1e6))';
    E_prime_interp = interp1(freq,E_prime,freq_interp);
    E_dbl_prime_interp = interp1(freq,E_dbl_prime,freq_interp);
    freq_xover_tmp = freq_interp(find(E_prime_interp<=E_dbl_prime_interp,1,'last'));
    if ~isempty(freq_xover_tmp)
        freq_xover(i) = freq_xover_tmp;
    end
    p = plot([freq_xover(i) freq_xover(i)],[0.0001 1e6]);
    p.Color = color;
    p.LineStyle = '-.';
    p.LineWidth = 1;

    e1 = errorbar(freq,E_prime,E_prime_err);
    e1.MarkerSize = 8;
    e1.Marker = 'o';
    e1.MarkerFaceColor = color;
    e1.MarkerEdgeColor = 'k';
    e1.Color = 'k';


    e2 = errorbar(freq,E_dbl_prime,E_dbl_prime_err);
    e2.MarkerSize = 8;
    e2.Marker = '^';
    e2.MarkerFaceColor = color;
    e2.MarkerEdgeColor = 'k';
    e2.Color = 'k';
    e2.LineStyle = '--';


    set(gca,'FontSize',20/1.5)
    set(gca,'xscale','log')
    set(gca,'yscale','log')

    set(gcf,'color','w')

    if i==length(kds)
        xlabel('$f \tau_0$','FontSize',20,'Interpreter','latex')
    else
        xticklabels({})
    end
%     ylabel('$G^{\prime}$, $G^{\prime \prime}$ (Pa)','FontSize',20,'Interpreter','latex')
    ylabel('$G^{\prime}$, $G^{\prime \prime}$ (kPa)','FontSize',20,'Interpreter','latex')
    xlim([4e-5 1.5e-2])

    pbaspect([3 1 1])
    set(gcf,'position',[1000 100 600 300])

    yyaxis right
%     delta_color = [0.3 0 0.9];
    delta_color = [0.5 0.5 0.5];
    ax = gca;
    ax.YAxis(2).Color = delta_color*0.5;

    e3 = errorbar(freq,E_dbl_prime./E_prime,delta_err);
    e3.Marker = 'sq';
    e3.MarkerSize = 6;
    e3.MarkerFaceColor = 'none';
    e3.MarkerEdgeColor = delta_color;
    e3.Color = delta_color;
    e3.LineStyle = '-.';
    ylabel('$\tan (\delta)$','FontSize',20,'Interpreter','latex')
    ylim([0 3])
    yticks([1 2])
    
end
ratio = kds*tau0./freq_xover;


    % yyaxis right
    %
    % e = errorbar(freq,delta,delta_err);
    % e.Marker = '^';
    % ylabel('tan$(\delta)$','FontSize',20,'Interpreter','latex')
% 
% 
%     figure(2); clf
% 
%     subplot(2,1,1); hold on
%     tf = 10/omega_min*tau0*1e6;
%     time = linspace(0,tf,1e5);
%     omega_max_temp = omega_max/10;
%     omega = omega_min*exp(time/tf*log(omega_max_temp/omega_min));
% 
%     p = plot(time,omega,'k');
%     p.LineWidth = 1.5;
%     set(gca,'Fontsize',20/1.5)
% 
%     xticks((0:50:200))
%     xticklabels({})
% 
%     ylabel('$\omega(t) \tau_0$','FontSize',20,'Interpreter','latex')
%     set(gca,'yscale','log')
% 
%     subplot(2,1,2);  hold on
%     eps0 = lambda-1;
% 
%     eps = eps0*sin(2*pi*omega.*time/tau0*1e-6);
% 
%     p = plot(time,eps,'k');
%     p.LineWidth = 1.5;
%     p = plot(time,-eps,'k:');
%     p.LineWidth = 1.5;
%     set(gca,'FontSize',20/1.5)
% 
%     ylabel('$\epsilon(t)$','FontSize',20,'Interpreter','latex')
%     xlabel('$t$ ($\mu$s)','FontSize',20,'Interpreter','latex')
% 
%     set(gcf,'color','w')
%     ylim([-0.15 0.15])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotDetachmentRateSweep(phi,N_Kuhn,kds,loading_rate,folder_name,n_pts,Controls)

global raw_data_folder_name tau0 sig11_lmp_temp

n_fig = 15;
InitiateTheFigures(n_fig,Controls)
xsize = 400;
ysize = 400;

param_rng = kds([1 3 6]);
% param_rng = kds([2 3 6]);
% param_rng = kds; param_rng(5) = [];
c1 = [0 0 0];
c2 = [1 0 0];
n_colors = length(kds);
colors = [(linspace(c1(1),c2(1),n_colors))',...
    (linspace(c1(2),c2(2),n_colors))',...
    (linspace(c1(3),c2(3),n_colors))'];

peak_abs_err_sig11 = zeros(length(kds),1);
mean_rel_err_sig11 = zeros(length(kds),1);
mean_rel_err_ka = zeros(length(kds),1);
mean_rel_err_kd = zeros(length(kds),1);
mean_rel_err_fa = zeros(length(kds),1);
mean_rel_err_fd = zeros(length(kds),1);
mean_rel_err_g11 = zeros(length(kds),1);
mean_rel_err_g12 = zeros(length(kds),1);
mean_rel_err_g22 = zeros(length(kds),1);
mean_rel_err_MSDst = zeros(length(kds),1);
mean_rel_err_MSDth = zeros(length(kds),1);
mean_rel_err_ca = zeros(length(kds),1);
mean_rel_err_nbonds = zeros(length(kds),1);
mean_rel_err_nothers = zeros(length(kds),1);
mean_rel_err_nself = zeros(length(kds),1);

falpha_bead = zeros(length(kds),1);
kalpha_bead = zeros(length(kds),1);
chi_bead = zeros(length(kds),1);
SEfalpha_bead = zeros(length(kds),1);
SEkalpha_bead = zeros(length(kds),1);
SEchi_bead = zeros(length(kds),1);
falpha_meso = zeros(length(kds),1);
kalpha_meso = zeros(length(kds),1);
chi_meso = zeros(length(kds),1);
SEfalpha_meso = zeros(length(kds),1);
SEkalpha_meso = zeros(length(kds),1);
SEchi_meso = zeros(length(kds),1);

kreconfig_bead = zeros(length(kds),1);
SEkreconfig_bead = zeros(length(kds),1);
kreconfig_meso = zeros(length(kds),1);
SEkreconfig_meso = zeros(length(kds),1);

for iv=1:length(kds)
    kd = kds(iv);
    color = colors(iv,:);
    file_tag = ['N_',num2str(N_Kuhn),...
        '.kd_',num2str(kd,'%.2e'),...
        '.W_',num2str(loading_rate,'%.3f'),...
        '.phi_',num2str(phi,'%.3f')];

    % Import the Conslidated Data

    % Time and stretch
    file_name = [raw_data_folder_name,'/','timestretch.',file_tag,'.mat'];
    timestretch = load(file_name,'-mat');
    time_ld = timestretch.time_ld;
    time_eq = timestretch.time_eq;

    % Stress
    file_name = [raw_data_folder_name,'/','stress.',file_tag,'.mat'];
    stress = load(file_name,'-mat');
    sig11 = stress.sig11;
    sig11_lmp = stress.sig11_lmp;

    % Bond kinetics and topological data
    file_name = [raw_data_folder_name,'/','kinetics.',file_tag,'.mat'];
    kinetics = load(file_name,'-mat');
    ka_out = kinetics.ka;
    kd_out = kinetics.kd;
    fa = kinetics.fa;
    fd = kinetics.fd;
    n_attach = kinetics.n_attach;
    n_detach = kinetics.n_detach;
    n_bonds = kinetics.n_bonds;
    n_self = kinetics.n_self;
    n_others = kinetics.n_others;
    cluster_coeff = kinetics.cluster_coeff;

    % Chain alignment
    file_name = [raw_data_folder_name,'/','covariance.',file_tag,'.mat'];
    covariance = load(file_name,'-mat');
    rxr11 = covariance.rxr11_all;
    rxr22 = covariance.rxr22_all;
    rxr12 = covariance.rxr12_all;

    % MSD
    file_name = [raw_data_folder_name,'/','msd.',file_tag,'.mat'];
    MSD = load(file_name,'-mat');
    msd_sticker = MSD.msd_st;
    msd_tether = MSD.msd_th;

    sig11_lmp_temp = (sig11_lmp-sig11_lmp(1,:,:))/1e6;
    sig11_lmp_temp = squeeze(nanmean(sig11_lmp_temp,2));

    % Stress in loading direction
    figno = 1;
    figure(figno)
    x = time_ld/1e-6;
    y_tmp = (sig11_lmp-sig11_lmp(1,:,:))/1e6;
    y_lmp = squeeze(nanmean(y_tmp,2));
    err_lmp = squeeze(nanstd(y_tmp,1,2));
    y_tmp = 2*sig11/1e6;
    y_ntwrk = squeeze(nanmean(y_tmp,2));
    err_ntwrk = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$\sigma_{11}$ (MPa)';
    tag = 'sig11';
    ylims = [0 2.5];
    xlims = [0 Inf];
    PlotErr = 1;

    if ismember(kd,param_rng)
        PlotPolishedOutputLAMMPS(kd,kds,figno,x,y_lmp(:,1),err_lmp(:,1),...
            y_ntwrk(:,1),err_ntwrk(:,1),xlab,ylab,xlims,ylims,color,...
            file_tag,folder_name,tag,PlotErr);
    end


    % stress in loading direction from post-processing virial formulation
    figno = figno+1;
    figure(figno)
    x = time_ld/1e-6;
    y_tmp = 2*sig11/1e6;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$\sigma_{11}$ (MPa)';
    tag = 'sig_11_endtoend';
    ylims = [0 2.5];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 1;

    [~,mean_rel_err_sig11(iv),~,~,peak_abs_err_sig11(iv)] = ...
        PlotPolishedOutputDetachmentRates(kd,kds,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);
    
    % Compute relaxation timescale
    figno = figno+1;
    figure(figno); hold on
    rlx_indx = find(y(:,2)==max(y(:,2)),1,'first');
    x = time_ld(rlx_indx:end)/(1e3*tau0);
    x = x-x(1);
    sig11_bead = y(rlx_indx:end,1)/y(rlx_indx,1);
    sig11_meso = y(rlx_indx:end,2)/y(rlx_indx,2);
    rng = (ceil(linspace(1,length(x),n_pts)))';

    % Fit the relaxation timescale
    ft = fittype(['A*exp(-x^C*B) + (1-A)*exp(-x*D*',num2str(kd*1e3*tau0),')']);
%     options = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 15 0.1],'Upper',[1 Inf 1]);
%     options = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 15 0.1],'Upper',[1 Inf 1]);
%     options = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 15 0.9999],'Upper',[1 Inf 1]);
    chi_lo = 0.8; chi_hi = 1;
    lowers = [0 15 chi_lo];
    uppers = [1 Inf chi_hi];
    options = fitoptions('Method', 'NonlinearLeastSquares','Lower',lowers,'Upper',uppers);
    ft_meso = ft;
    ft_bead = ft;

    sig11_ft = sig11_bead(~isnan(sig11_bead));
    x_ft = x(~isnan(sig11_bead));
    func_bead = fit(x_ft,sig11_ft,ft_bead,options);

    sig11_ft = sig11_meso(~isnan(sig11_meso));
    x_ft = x(~isnan(sig11_meso));
    func_meso = fit(x_ft,sig11_ft,ft_meso,options);



    if kd==kds(end)
        lowers = [func_meso.A func_meso.B func_meso.C func_meso.D]*0.85;
        uppers = [func_meso.A func_meso.B func_meso.C func_meso.D]*1.15;
        lowers(2) = 15;
        uppers(2) = Inf;
        lowers(3) = chi_lo; 
        uppers(3) = chi_hi;
        options = fitoptions('Method', 'NonlinearLeastSquares','Lower',lowers,'Upper',uppers);
        func_bead = fit(x_ft,sig11_ft,ft_bead,options);
    end
    
    val_bead = coeffvalues(func_bead);
    pm_bead = mean(abs(confint(func_bead)-val_bead),1);
    val_meso = coeffvalues(func_meso);
    pm_meso = mean(abs(confint(func_meso)-val_meso),1);

    falpha_bead(iv) = val_bead(1);
    SEfalpha_bead(iv) = pm_bead(1);

    falpha_meso(iv) = val_meso(1);
    SEfalpha_meso(iv) = pm_meso(1);

    kalpha_bead(iv) = val_bead(2)^(1/val_bead(3));
    SEkalpha_bead(iv) = pm_bead(2);

    kalpha_meso(iv) = val_meso(2)^(1/val_meso(3));
    SEkalpha_meso(iv) = pm_meso(2);

    chi_bead(iv) = val_bead(3);
    SEchi_bead(iv) = pm_bead(3);

    chi_meso(iv) = val_meso(3);
    SEchi_meso(iv) = pm_meso(3);

    kr_bead(iv) = val_bead(4)*kd*1e3*tau0;
    SEkr_bead(iv) = pm_bead(4)*kd*1e3*tau0;
    kr_meso(iv) = val_meso(4)*kd*1e3*tau0;
    SEkr_meso(iv) = pm_meso(4)*kd*1e3*tau0;

    kreconfig_bead(iv) = kd;
    SEkreconfig_bead(iv) = 0;
    kreconfig_meso(iv) = kd;
    SEkreconfig_meso(iv) = 0;

    if ismember(kd,param_rng)
        [~] = PlotErrorBar(x(rng),sig11_bead(rng),NaN*ones(size(rng)),color,'filled');
        e = PlotErrorBar(x(rng),sig11_meso(rng),NaN*ones(size(rng)),color,'');
        e.Marker = '^';
        p = plot(x,func_bead(x));
        p.LineWidth = 1;
        p.LineStyle = '-';
        p.Color = color;

        p = plot(x,func_meso(x));
        p.LineWidth = 1;
        p.LineStyle = '--';
        p.Color = color;
    end

    if kd==kds(end)
        set(gca,'FontSize',20/1.5)
        set(gcf,'color','w')
        set(gcf,'Position',[1000 100 400 400])
        xlabel('$t/(\tau_0 \cdot 10^3)$','FontSize',20,'Interpreter','latex')
        ylabel('$\sigma_{11}^*$','FontSize',20,'Interpreter','latex')
    end

    % Attachment rate
    figno = figno+1;
    figure(figno)
    x = time_eq/tau0;
    y_tmp = ka_out*tau0;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t/\tau_0$';
    ylab = '$k_a\tau_0$';
    tag = 'ka';
    ylims = [0 Inf];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 1;

    [~,mean_rel_err_ka(iv),~,~,~] = ...
        PlotPolishedOutputDetachmentRates(kd,kds,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);
    set(gca,'yscale','log')
    ylim([1e-4 1e-1])

    % Detachment rate
    figno = figno+1;
    figure(figno)
    x = time_eq/tau0;
    y_tmp = kd_out*tau0;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$k_d\tau_0$';
    tag = 'kd';
    ylims = [0 Inf];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 1;

    [~,mean_rel_err_kd(iv),~,~,~] = ...
        PlotPolishedOutputDetachmentRates(kd,kds,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);


    % Attached and detached chain fractions
    figno = figno+1;
    figure(figno)
    x = time_eq/tau0;
    y_tmp = fa;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t/\tau_0$';
    ylab = '$f_a$, $f_d$';
    tag = 'f';
    ylims = [0 Inf];
    xlims = [0 Inf];
    PlotErr = 1;

    [~,mean_rel_err_fa(iv),~,~,~] = ...
        PlotPolishedOutputDetachmentRates(kd,kds,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);

    y_tmp = fd;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    [~,mean_rel_err_fd(iv),~,~,~] = ...
        PlotPolishedOutputDetachmentRates(kd,kds,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);

    % Alignment tensor g11
    figno = figno+1;
    figure(figno)
    x = time_ld/1e-6;
    y_tmp = rxr11;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t/\tau_0$';
    ylab = '$g_{11}$';
    tag = 'g11';
    ylims = [-0.1 1];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 1;

    plot(x,0.33*ones(size(x)),'k--')

    [~,mean_rel_err_g11(iv),~,~,~] = ...
        PlotPolishedOutputDetachmentRates(kd,kds,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);


    % Alignment tensor g22
    figno = figno+1;
    figure(figno)
    x = time_ld/tau0;
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$g_{22}$';
    tag = 'g22';
    ylims = [-0.1 1];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 1;

    plot(x,0.33*ones(size(x)),'k--')
    y_tmp = rxr22;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    [~,mean_rel_err_g22(iv),~,~,~] = ...
        PlotPolishedOutputDetachmentRates(kd,kds,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);


    % Alignemnt tensor g12
    figno = figno+1;
    figure(figno)
    x = time_ld/tau0;
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$g_{12}$';
    tag = 'g12';
    ylims = [-0.1 1];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 1;

    y_tmp = rxr12;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    [~,mean_rel_err_g12(iv),~,~,~] = ...
        PlotPolishedOutputDetachmentRates(kd,kds,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);

    % MSD st
    figno = figno+1;
    figure(figno)
    x = time_ld/tau0;
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$\langle \Delta x^2 \rangle_{st}/b^2$';
    tag = 'MSDst';
    ylims = [-0.1 1];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 1;

    y_tmp = msd_sticker;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    [~,mean_rel_err_MSDst(iv),~,~,~] = ...
        PlotPolishedOutputDetachmentRates(kd,kds,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);

    % MSD th
    figno = figno+1;
    figure(figno)
    x = time_ld/tau0;
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$\langle \Delta x^2 \rangle_{th}/b^2$';
    tag = 'MSDst';
    ylims = [-0.1 1];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 1;

    y_tmp = msd_tether;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    [~,mean_rel_err_MSDth(iv),~,~,~] = ...
        PlotPolishedOutputDetachmentRates(kd,kds,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);
    
    % clustering coefficient
    figno = figno+1;
    figure(figno)
    x = time_ld/tau0;
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$\langle c_\alpha \rangle$';
    tag = 'ca';
    ylims = [-0.1 1];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 1;

    y_tmp = cluster_coeff;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    [~,mean_rel_err_ca(iv),~,~,~] = ...
        PlotPolishedOutputDetachmentRates(kd,kds,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);

    % number of bonds 
    figno = figno+1;
    figure(figno)
    x = time_ld/tau0;
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$n_{bonds}$';
    tag = 'nbonds';
    ylims = [-0.1 1];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 1;

    y_tmp = n_bonds;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    [~,mean_rel_err_nbonds(iv),~,~,~] = ...
        PlotPolishedOutputDetachmentRates(kd,kds,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);
    
    % number of bonds to others
    figno = figno+1;
    figure(figno)
    x = time_ld/tau0;
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$n_{others}$';
    tag = 'nothers';
    ylims = [-0.1 1];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 1;

    y_tmp = n_others;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    [~,mean_rel_err_nothers(iv),~,~,~] = ...
        PlotPolishedOutputDetachmentRates(kd,kds,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);

    % number of bonds to self
    figno = figno+1;
    figure(figno)
    x = time_ld/tau0;
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$n_{self}$';
    tag = 'nself';
    ylims = [-0.1 1];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 1;

    y_tmp = n_self;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    [~,mean_rel_err_nself(iv),~,~,~] = ...
        PlotPolishedOutputDetachmentRates(kd,kds,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);
end

close all
figno = 0;
font_size = 20;
x = log10(kds(2:end)*tau0);
xlims = [-3.5 -0.5];
xlab = 'log$(k_d \tau_0)$';
xsize = 225; 
ysize = 225;
fig_aspect = [1 1 1];

% outputs = {'sig11','g11','g22','ka','fa','nbonds','nother','nself','ca','MSDst','MSDth'}';
outputs = {'sig11_abs','sig11_rel','g11','g22','ka','kd','fa','MSDst','MSDth','nbonds','nothers','nself','ca'}';
slope = zeros(length(outputs),1);
SEm = zeros(length(outputs),1);
intercept = zeros(length(outputs),1);
SEb = zeros(length(outputs),1);
R2 = zeros(length(outputs),1);
p_val = zeros(length(outputs),1);

figno = figno + 1;
y = peak_abs_err_sig11(2:end);
ylims = [-Inf Inf];
ylab = '$| Err. |_{max}$';
[~,f2_sig11,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','sig11_error.png'])
saveas(gcf,[folder_name,'/','sig11_error.fig'])
intercept(figno) = f2_sig11.Coefficients.Estimate(1);
slope(figno) = f2_sig11.Coefficients.Estimate(2);
SEm(figno) = f2_sig11.Coefficients.SE(2);
SEb(figno) = f2_sig11.Coefficients.SE(1);
R2(figno) = f2_sig11.Rsquared.Ordinary;
p_val(figno) = f2_sig11.anova.pValue(1);
% set(gca,'xscale','log')

figno = figno + 1;
y = mean_rel_err_sig11(2:end);
ylims = [-Inf Inf];
ylab = '$\langle \% Err. \rangle$';
[~,f2_sig11_rel,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','sig11_error.png'])
saveas(gcf,[folder_name,'/','sig11_error.fig'])
intercept(figno) = f2_sig11_rel.Coefficients.Estimate(1);
slope(figno) = f2_sig11_rel.Coefficients.Estimate(2);
SEm(figno) = f2_sig11_rel.Coefficients.SE(1);
SEb(figno) = f2_sig11_rel.Coefficients.SE(2);
R2(figno) = f2_sig11_rel.Rsquared.Ordinary;
p_val(figno) = f2_sig11_rel.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_g11(2:end);
ylims = [-10 10];
ylab = '$\langle \% Err. \rangle$';
[~,f2_g11,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','g11_error.png'])
saveas(gcf,[folder_name,'/','g11_error.fig'])
intercept(figno) = f2_g11.Coefficients.Estimate(1);
slope(figno) = f2_g11.Coefficients.Estimate(2);
SEm(figno) = f2_g11.Coefficients.SE(1);
SEb(figno) = f2_g11.Coefficients.SE(2);
R2(figno) = f2_g11.Rsquared.Ordinary;
p_val(figno) = f2_g11.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_g22(2:end);
ylims = [-10 10];
ylab = '$\langle \% Err. \rangle$';
[~,f2_g22,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','g22_error.png'])
saveas(gcf,[folder_name,'/','g22_error.fig'])
intercept(figno) = f2_g22.Coefficients.Estimate(1);
slope(figno) = f2_g22.Coefficients.Estimate(2);
SEm(figno) = f2_g22.Coefficients.SE(1);
SEb(figno) = f2_g22.Coefficients.SE(2);
R2(figno) = f2_g22.Rsquared.Ordinary;
p_val(figno) = f2_g22.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_ka(2:end);
ylims = [-Inf Inf];
ylab = '$\langle \% Err. \rangle$';
[~,f2_ka_rel,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','ka_error.png'])
saveas(gcf,[folder_name,'/','ka_error.fig'])
intercept(figno) = f2_ka_rel.Coefficients.Estimate(1);
slope(figno) = f2_ka_rel.Coefficients.Estimate(2);
SEm(figno) = f2_ka_rel.Coefficients.SE(1);
SEb(figno) = f2_ka_rel.Coefficients.SE(2);
R2(figno) = f2_ka_rel.Rsquared.Ordinary;
p_val(figno) = f2_ka_rel.anova.pValue(1);

if nansum(mean_rel_err_kd)~=0
    figno = figno + 1;
    y = mean_rel_err_kd(2:end);
    ylab = '$\langle \% Err. \rangle$';
    ylims = [-Inf Inf];
    [~,f2_kd,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
        xlab,ylab,font_size,fig_aspect);
    saveas(gcf,[folder_name,'/','kd_error.png'])
    saveas(gcf,[folder_name,'/','kd_error.fig'])
    intercept(figno) = f2_kd.Coefficients.Estimate(1);
    slope(figno) = f2_kd.Coefficients.Estimate(2);
    SEm(figno) = f2_kd.Coefficients.SE(1);
    SEb(figno) = f2_kd.Coefficients.SE(2);
    R2(figno) = f2_kd.Rsquared.Ordinary;
    p_val(figno) = f2_kd.anova.pValue(1);
end

figno = figno + 1;
y = mean_rel_err_fa(2:end);
ylims = [-Inf Inf];
ylab = '$\langle \% Err. \rangle$';
[~,f2_fa_rel,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','fa_error.png'])
saveas(gcf,[folder_name,'/','fa_error.fig'])
intercept(figno) = f2_fa_rel.Coefficients.Estimate(1);
slope(figno) = f2_fa_rel.Coefficients.Estimate(2);
SEm(figno) = f2_fa_rel.Coefficients.SE(1);
SEb(figno) = f2_fa_rel.Coefficients.SE(2);
R2(figno) = f2_fa_rel.Rsquared.Ordinary;
p_val(figno) = f2_fa_rel.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_MSDst(2:end);
ylims = [-Inf Inf];
ylab = '$\langle \% Err. \rangle$';
[~,f2_MSDst_rel,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','MSDst_error.png'])
saveas(gcf,[folder_name,'/','MSDst_error.fig'])
intercept(figno) = f2_MSDst_rel.Coefficients.Estimate(1);
slope(figno) = f2_MSDst_rel.Coefficients.Estimate(2);
SEm(figno) = f2_MSDst_rel.Coefficients.SE(1);
SEb(figno) = f2_MSDst_rel.Coefficients.SE(2);
R2(figno) = f2_MSDst_rel.Rsquared.Ordinary;
p_val(figno) = f2_MSDst_rel.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_MSDth(2:end);
ylims = [-Inf Inf];
ylab = '$\langle \% Err. \rangle$';
[~,f2_MSDth_rel,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','MSDth_error.png'])
saveas(gcf,[folder_name,'/','MSDth_error.fig'])
intercept(figno) = f2_MSDth_rel.Coefficients.Estimate(1);
slope(figno) = f2_MSDth_rel.Coefficients.Estimate(2);
SEm(figno) = f2_MSDth_rel.Coefficients.SE(1);
SEb(figno) = f2_MSDth_rel.Coefficients.SE(2);
R2(figno) = f2_MSDth_rel.Rsquared.Ordinary;
p_val(figno) = f2_MSDth_rel.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_nbonds(2:end);
ylims = [-Inf Inf];
ylab = '$\langle \% Err. \rangle$';
[~,f2_nbonds_rel,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','nbonds_error.png'])
saveas(gcf,[folder_name,'/','nbonds_error.fig'])
intercept(figno) = f2_nbonds_rel.Coefficients.Estimate(1);
slope(figno) = f2_nbonds_rel.Coefficients.Estimate(2);
SEm(figno) = f2_nbonds_rel.Coefficients.SE(1);
SEb(figno) = f2_nbonds_rel.Coefficients.SE(2);
R2(figno) = f2_nbonds_rel.Rsquared.Ordinary;
p_val(figno) = f2_nbonds_rel.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_nothers(2:end);
ylims = [-Inf Inf];
ylab = '$\langle \% Err. \rangle$';
[~,f2_nothers_rel,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','nothers_error.png'])
saveas(gcf,[folder_name,'/','nothers_error.fig'])
intercept(figno) = f2_nothers_rel.Coefficients.Estimate(1);
slope(figno) = f2_nothers_rel.Coefficients.Estimate(2);
SEm(figno) = f2_nothers_rel.Coefficients.SE(1);
SEb(figno) = f2_nothers_rel.Coefficients.SE(2);
R2(figno) = f2_nothers_rel.Rsquared.Ordinary;
p_val(figno) = f2_nothers_rel.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_nself(2:end);
ylims = [-Inf Inf];
ylab = '$\langle \% Err. \rangle$';
[~,f2_nself_rel,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','nself_error.png'])
saveas(gcf,[folder_name,'/','nself_error.fig'])
intercept(figno) = f2_nself_rel.Coefficients.Estimate(1);
slope(figno) = f2_nself_rel.Coefficients.Estimate(2);
SEm(figno) = f2_nself_rel.Coefficients.SE(1);
SEb(figno) = f2_nself_rel.Coefficients.SE(2);
R2(figno) = f2_nself_rel.Rsquared.Ordinary;
p_val(figno) = f2_nself_rel.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_ca(2:end);
ylims = [-Inf Inf];
ylab = '$\langle \% Err. \rangle$';
[~,f2_ca_rel,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','ca_error.png'])
saveas(gcf,[folder_name,'/','ca_error.fig'])
intercept(figno) = f2_ca_rel.Coefficients.Estimate(1);
slope(figno) = f2_ca_rel.Coefficients.Estimate(2);
SEm(figno) = f2_ca_rel.Coefficients.SE(1);
SEb(figno) = f2_ca_rel.Coefficients.SE(2);
R2(figno) = f2_ca_rel.Rsquared.Ordinary;
p_val(figno) = f2_ca_rel.anova.pValue(1);

tab = table(outputs,slope,SEm,intercept,SEb,R2,p_val);
writetable(tab,[folder_name,'/Error Analysis Table.txt'])

tab = table(kds,falpha_bead,SEfalpha_bead,kreconfig_bead,SEkreconfig_bead,...
    falpha_meso,SEfalpha_meso,kreconfig_meso,SEkreconfig_meso);
writetable(tab,[folder_name,'/Relaxation Timescale Analysis Table.txt'])


kds(1) = kds(3)/10; %to plot in lieu of 0

figno = figno + 1; 
c1 = [0.2 0 0.5];
c2 = [0 0.5 0];
c3 = [0.8 0.4 0.1];
c4 = [0.1 0.2 0.8];
figure(figno); clf; hold on
yyaxis left
set(gca,'FontSize',20/1.5)
e = errorbar(kds,falpha_bead,SEfalpha_bead);
e.Marker = 'o';
e.LineWidth = 1.5;
% e.LineStyle = 'none';
e.Color = c1;
e.MarkerEdgeColor = c1;
% e.MarkerFaceColor = 'k';

e = errorbar(kds,falpha_meso,SEfalpha_meso);
e.Marker = '^';
e.LineWidth = 1.5;
% e.LineStyle = 'none';
e.Color = c1;
e.MarkerEdgeColor = c1;
% e.MarkerFaceColor = 'none';
ylabel('$f_\alpha$','FontSize',20,'Interpreter','latex')
ylim([0 1])
set(gca, 'YColor',c1); % Set left y-axis color

yyaxis right
e = errorbar(kds,chi_bead,SEchi_bead);
e.Marker = 'o';
e.LineWidth = 1.5;
% e.LineStyle = 'none';
e.Color = c3;
e.MarkerEdgeColor = c3;
% e.MarkerFaceColor = 'k';

e = errorbar(kds,chi_meso,SEchi_meso);
e.Marker = '^';
e.LineWidth = 1.5;
% e.LineStyle = 'none';
e.Color = c3;
e.MarkerEdgeColor = c3;
% e.MarkerFaceColor = 'none';
ylabel('$\chi$','FontSize',20,'Interpreter','latex')
ylim([0 1.05])
set(gca, 'YColor',c3); % Set left y-axis color

xlabel('$k_d \tau_0$','FontSize',20,'Interpreter','latex')
set(gcf,'Color','w')
set(gcf,'Position',[1000 100 400 400])
set(gca,'xscale','log')
xticks(flipud(loading_rate))
% xlim([0.02 3.1]*10^7)
xlim([0.007 3.1]*10^7)
xticks(kds)
xticklabels({'0','0.001','0.003','0.01','0.03','0.1'})
pbaspect([2 1 1])


figno = figno + 1; 
figure(figno); clf; hold on
yyaxis left
set(gca,'FontSize',20/1.5)
e = errorbar(kds,kalpha_bead,SEfalpha_bead);
e.Marker = 'o';
e.LineWidth = 1.5;
% e.LineStyle = 'none';
e.Color = c2;
e.MarkerEdgeColor = c2;
% e.MarkerFaceColor = 'k';

e = errorbar(kds,kalpha_meso,SEfalpha_meso);
e.Marker = '^';
e.LineWidth = 1.5;
% e.LineStyle = 'none';
e.Color = c2;
e.MarkerEdgeColor = c2;
% e.MarkerFaceColor = 'none';
ylabel('$k_\alpha \tau_0\cdot 10^3$','FontSize',20,'Interpreter','latex')
ylim([0 80])
set(gca, 'YColor',c2); % Set left y-axis color

yyaxis right
e = errorbar(kds,kr_bead,SEkr_bead);
e.Marker = 'o';
e.LineWidth = 1.5;
% e.LineStyle = 'none';
e.Color = c4;
e.MarkerEdgeColor = c4;
% e.MarkerFaceColor = 'k';

e = errorbar(kds,kr_meso,SEkr_meso);
e.Marker = '^';
e.LineWidth = 1.5;
% e.LineStyle = 'none';
e.Color = c4;
e.MarkerEdgeColor = c4;
% e.MarkerFaceColor = 'none';
ylabel('$k_r \tau_0\cdot 10^3$','FontSize',20,'Interpreter','latex')
ylim([-0.05 8])
set(gca, 'YColor',c4); % Set left y-axis color

xlabel('$k_d \tau_0$','FontSize',20,'Interpreter','latex')
set(gcf,'Color','w')
set(gcf,'Position',[1000 100 400 400])
set(gca,'xscale','log')
xticks(flipud(loading_rate))
% xlim([0.02 3.1]*10^7)
xlim([0.007 3.1]*10^7)
xticks(kds)
xticklabels({'0','0.001','0.003','0.01','0.03','0.1'})
pbaspect([2 1 1])
set(gcf,'Position',[1000 100 380 400])


f2_bs_chi = fitlm(log(kds),chi_bead);
disp(['p-value of chi vs kd for bs is ',num2str(f2_bs_chi.anova.pValue(1))]);
f2_ms_chi = fitlm(log(kds),chi_meso);
disp(['p-value of chi vs kd for ms is ',num2str(f2_ms_chi.anova.pValue(1))]);

f2_bs_falpha = fitlm(log(kds),falpha_bead);
disp(['p-value of falpha vs kd for bs is ',num2str(f2_bs_falpha.anova.pValue(1))]);
f2_ms_falpha = fitlm(log(kds),falpha_meso);
disp(['p-value of falpha vs kd for ms is ',num2str(f2_ms_falpha.anova.pValue(1))]);

f2_bs_kalpha = fitlm(log(kds),kalpha_bead);
disp(['p-value of kalpha vs kd for bs is ',num2str(f2_bs_kalpha.anova.pValue(1))]);
f2_ms_kalpha = fitlm(log(kds),kalpha_meso);
disp(['p-value of kalpha vs kd for ms is ',num2str(f2_ms_kalpha.anova.pValue(1))]);

f2_bs_kr = fitlm(log(kds),kr_bead);
disp(['p-value of kr vs kd for bs is ',num2str(f2_bs_kr.anova.pValue(1))]);
f2_ms_kr = fitlm(log(kds),kr_meso);
disp(['p-value of kr vs kd for ms is ',num2str(f2_ms_kr.anova.pValue(1))]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotLargeDeformationResponse(phi,N_Kuhn,kd,loading_rates,folder_name,n_pts,Controls)

global raw_data_folder_name tau0 sig11_lmp_temp

n_fig = 3;
InitiateTheFigures(n_fig,Controls)
xsize = 400;
ysize = 400;

param_rng = loading_rates;%([1 2]);
c1 = [0 0 0];
c2 = [1 0.6471 0];
n_colors = length(loading_rates);
colors = [(linspace(c1(1),c2(1),n_colors))',...
    (linspace(c1(2),c2(2),n_colors))',...
    (linspace(c1(3),c2(3),n_colors))'];

max_stretches = zeros(length(loading_rates),1);

for iv=1:length(loading_rates)
    loading_rate = loading_rates(iv);
    color = colors(iv,:);
    file_tag = ['N_',num2str(N_Kuhn),...
        '.kd_',num2str(kd,'%.2e'),...
        '.W_',num2str(loading_rate,'%.5f'),...
        '.phi_',num2str(phi,'%.3f')];

    % Import the Conslidated Data

    % Time and stretch
    file_name = [raw_data_folder_name,'/','timestretch.',file_tag,'.mat'];
    timestretch = load(file_name,'-mat');
    stretch = timestretch.stretch;
    time_ld = timestretch.time_ld;
    time_eq = timestretch.time_eq;

    max_stretches(iv) = max(stretch);

    % Stress
    file_name = [raw_data_folder_name,'/','stress.',file_tag,'.mat'];
    stress = load(file_name,'-mat');
    sig11 = stress.sig11;
    sig22 = stress.sig22;
    sig12 = stress.sig12;

    % Bond kinetics and topological data
    file_name = [raw_data_folder_name,'/','kinetics.',file_tag,'.mat'];
    kinetics = load(file_name,'-mat');
    ka_out = kinetics.ka;
    kd_out = kinetics.kd;
    fa = kinetics.fa;
    fd = kinetics.fd;

    % Chain alignment
    file_name = [raw_data_folder_name,'/','covariance.',file_tag,'.mat'];
    covariance = load(file_name,'-mat');
    rxr11 = covariance.rxr11_all;
    rxr22 = covariance.rxr22_all;
    rxr12 = covariance.rxr12_all;

    % stress in loading direction from post-processing virial formulation
    figno = 1;
    figure(figno)
    x = time_ld/1e-6;
    y_tmp = 2*sig11(:,:,2)/1e3;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    xlab = '$t$ ($\mu$s)';
    ylab = '$\sigma_{11}$ (kPa)';
    tag = 'sig_11_endtoend';
    ylims = [0 1500];
    xlims = [0 Inf];

    pl(iv) = ...
        PlotPolishedOutputCreep(loading_rate,loading_rates,figno,x,y,err,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,param_rng,xsize);

    pbaspect([2 1 1])
    set(gcf,'position',[1000 100 600 400])


% %     % stress vs strain
% %     figno = figno+1;
% %     figure(figno)
% % %     ld_start_indx = find(stretch>1,1,'first')-1;
% % %     stretch_plot = interp1(time_ld(ld_start_indx:end),stretch(ld_start_indx:end),time_ld);
% %     stretch_plot = (linspace(1,max(stretch),length(time_ld)))';
% %     x = stretch_plot-1;
% %     y_tmp = 2*sig11(:,:,2)/1e3;
% %     y = squeeze(nanmean(y_tmp,2));
% %     err = squeeze(nanstd(y_tmp,1,2));
% % 
% %     xlab = '$\epsilon$';
% %     ylab = '$\sigma_{11}$ (kPa)';
% %     tag = 'sig_11_endtoend_vs_strain';
% % %     ylims = [0 1500];
% % %     xlims = [1 20];
% % %     xticks([1 5 10 15 20])
% %     xlims = [0 1];
% %     ylims = [0 100];
% % 
% %     pl(iv) = ...
% %         PlotPolishedOutputCreep(loading_rate,loading_rates,figno,x,y,err,xlab,ylab,xlims,ylims,color,...
% %         file_tag,folder_name,tag,param_rng,xsize);
% % 
% %     pbaspect([1 1 1])
% %     xticks([0 0.5 1]);
% %     plot([0.05 0.05],[0 100],'k--')

    % stress vs strain
    figno = figno+1;
    figure(figno)
    stretch_plot = (linspace(1,max(stretch),length(time_ld)))';
%     x = stretch_plot-1;
    x = log(stretch_plot);
    y_tmp = 2*sig11(:,:,2)/1e3;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    xlab = '$\varepsilon$';
    ylab = '$\sigma_{11}$ (kPa)';
    tag = 'sig_11_endtoend_vs_strain';
    xlims = [0 2];
    ylims = [0 500];

    pl(iv) = ...
        PlotPolishedOutputCreep(loading_rate,loading_rates,figno,x,y,err,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,param_rng,xsize);

    pbaspect([1 1 1])
    xticks([0 2]);
    yticks([0 500])

    set(gcf,'position',[1000 100 225 225])
    box on
    
    % stress vs stretch
    figno = figno+1;
    figure(figno)
    stretch_plot = (linspace(1,max(stretch),length(time_ld)))';
%     x = stretch_plot-1;
    x = stretch_plot;
    y_tmp = 2*sig11(:,:,2)/1e3;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    xlab = '$\lambda$';
    ylab = '$\sigma_{11}$ (kPa)';
    tag = 'sig_11_endtoend_vs_stretch';
    xlims = [1 3];
    ylims = [0 500];

    pl(iv) = ...
        PlotPolishedOutputCreep(loading_rate,loading_rates,figno,x,y,err,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,param_rng,xsize);

    pbaspect([1 1 1])
    xticks([0 1 3]);
    yticks([0 500])

    set(gcf,'position',[1000 100 225 225])
    box on

    
    % stretch normal to loading direction from post-processing virial formulation
    figno = figno+1;
    figure(figno)
    x = time_ld/1e-6;
    if iv==1
        y_tmp = stretch;
    else
        y_tmp = NaN*ones(size(stretch));
    end
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    xlab = '$t$';
    ylab = 'ln($\lambda$)';
    tag = 'stretch';
    ylims = [0 20];
    xlims = [0 Inf];

    [~] = ...
        PlotPolishedOutputCreep(loading_rate,loading_rates,figno,x,y,err,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,param_rng,xsize);

    pbaspect([4 1 1])
    if iv==1
        plot([0 max(time_ld)/1e-6],[1 1],'k--')
        plot([0 max(time_ld)/1e-6],[20 20],'k--')
        xticks([0 max(time_ld)/1e-6])
        xticklabels({'0','ln(20)/$\dot \varepsilon$'})
        xaxisproperties= get(gca, 'XAxis');
        xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
        yticks([1 20])
        yticklabels({'0','ln(20)'})
        yaxisproperties= get(gca, 'YAxis');
        yaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis        
    end
    set(gca,'yscale','log')
    set(gcf,'position',[1000 100 600 400])
end

disp(max_stretches)
close all
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotLoadingRateSweep(phi,N_Kuhn,kd,loading_rates,folder_name,n_pts,Controls)

global raw_data_folder_name tau0 sig11_lmp_temp

n_fig = 7;
InitiateTheFigures(n_fig,Controls)
xsize = 400;
ysize = 400;

loading_rates = flipud(loading_rates);

param_rng = loading_rates([1 3 5]);
% c1 = [1 0 0];
c1 = [1 0.6471 0];
c2 = [0 0 0];
n_colors = length(loading_rates);
colors = [(linspace(c1(1),c2(1),n_colors))',...
    (linspace(c1(2),c2(2),n_colors))',...
    (linspace(c1(3),c2(3),n_colors))'];

% peak_abs_err_sig11_lmp = zeros(length(loading_rates),1);
peak_abs_err_sig11 = zeros(length(loading_rates),1);
% peak_abs_err_ka = zeros(length(loading_rates),1);
mean_rel_err_sig11 = zeros(length(loading_rates),1);
mean_rel_err_kd = zeros(length(loading_rates),1);

mean_rel_err_g11 = zeros(length(loading_rates),1);
mean_rel_err_g22 = zeros(length(loading_rates),1);
mean_rel_err_g12 = zeros(length(loading_rates),1);

peak_abs_err_g11 = zeros(length(loading_rates),1);
peak_abs_err_g12 = zeros(length(loading_rates),1);
peak_abs_err_g22 = zeros(length(loading_rates),1);
% peak_abs_err_fa = zeros(length(loading_rates),1);
% peak_abs_err_fd = zeros(length(loading_rates),1);
% peak_abs_err_na = zeros(length(loading_rates),1);
% peak_abs_err_nbonds = zeros(length(loading_rates),1);
% peak_abs_err_nself = zeros(length(loading_rates),1);
% peak_abs_err_nother = zeros(length(loading_rates),1);
% peak_abs_err_ca = zeros(length(loading_rates),1);
% peak_abs_err_MSDst = zeros(length(loading_rates),1);
% peak_abs_err_MSDth = zeros(length(loading_rates),1);
f_bead = zeros(length(loading_rates),1);
k_bead = zeros(length(loading_rates),1);
SEf_bead = zeros(length(loading_rates),1);
SEk_bead = zeros(length(loading_rates),1);
f_meso = zeros(length(loading_rates),1);
k_meso = zeros(length(loading_rates),1);
SEf_meso = zeros(length(loading_rates),1);
SEk_meso = zeros(length(loading_rates),1);

for iv=1:length(loading_rates)
    loading_rate = loading_rates(iv);
    color = colors(iv,:);
    file_tag = ['N_',num2str(N_Kuhn),...
        '.kd_',num2str(kd,'%.2e'),...
        '.W_',num2str(loading_rate,'%.3f'),...
        '.phi_',num2str(phi,'%.3f')];

    % Import the Conslidated Data

    % Time and stretch
    file_name = [raw_data_folder_name,'/','timestretch.',file_tag,'.mat'];
    timestretch = load(file_name,'-mat');
    time_ld = timestretch.time_ld;
    time_eq = timestretch.time_eq;
%     stretch = timestretch.stretch;

    % Stress
    file_name = [raw_data_folder_name,'/','stress.',file_tag,'.mat'];
    stress = load(file_name,'-mat');
    sig11 = stress.sig11;
%     sig22 = stress.sig22;
%     sig12 = stress.sig12;
    sig11_lmp = stress.sig11_lmp;
%     sig22_lmp = stress.sig22_lmp;
%     sig33_lmp = stress.sig33_lmp;

    % Bond kinetics and topological data
    file_name = [raw_data_folder_name,'/','kinetics.',file_tag,'.mat'];
    kinetics = load(file_name,'-mat');
%     ka_out = kinetics.ka;
    kd_out = kinetics.kd;
%     fa = kinetics.fa;
%     fd = kinetics.fd;
%     n_attach = kinetics.n_attach;
%     n_detach = kinetics.n_detach;
%     n_bonds = kinetics.n_bonds;
%     n_self = kinetics.n_self;
%     n_others = kinetics.n_others;
%     cluster_coeff = kinetics.cluster_coeff;

    % Chain alignment
    file_name = [raw_data_folder_name,'/','covariance.',file_tag,'.mat'];
    covariance = load(file_name,'-mat');
    rxr11 = covariance.rxr11_all;
    rxr22 = covariance.rxr22_all;
    rxr12 = covariance.rxr12_all;

    % MSD
%     file_name = [raw_data_folder_name,'/','msd.',file_tag,'.mat'];
%     MSD = load(file_name,'-mat');
%     msd_sticker = MSD.msd_st;
%     msd_tether = MSD.msd_th;

    sig11_lmp_temp = (sig11_lmp-sig11_lmp(1,:,:))/1e6;
    sig11_lmp_temp = squeeze(nanmean(sig11_lmp_temp,2));

    % stress in loading direction from LAMMPS virial formulation
%     figno = 1;
%     figure(figno)
%     x = time_ld/1e-6;
%     y_tmp = (sig11_lmp-sig11_lmp(1,:,:))/1e6;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$\sigma_{11}$ (MPa)';
%     tag = 'sig11_lmp';
%     ylims = [0 2.5];
%     xlims = [0 Inf];
%     PlotErr = 1;
%     err_type = 2;
% 
%     [MSE_sig11_lmp(iv),mean_rel_err_sig11_lmp(iv)] = ...
%         PlotPolishedOutputLoadingRates(loading_rate,loading_rates,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr  ,param_rng,xsize,ysize,err_type);
    
    % Stress in loading direction
    figno = 1;
    figure(figno)
    x = time_ld/1e-6;
    y_tmp = (sig11_lmp-sig11_lmp(1,:,:))/1e6;
    y_lmp = squeeze(nanmean(y_tmp,2));
    err_lmp = squeeze(nanstd(y_tmp,1,2));
    y_tmp = 2*sig11/1e6;
    y_ntwrk = squeeze(nanmean(y_tmp,2));
    err_ntwrk = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$\sigma_{11}$ (MPa)';
    tag = 'sig11';
    ylims = [0 2.5];
    xlims = [0 Inf];
    PlotErr = 1;

    if ismember(loading_rate,param_rng)
        PlotPolishedOutputLAMMPS(loading_rate,loading_rates,figno,x,y_lmp(:,1),err_lmp(:,1),...
            y_ntwrk(:,1),err_ntwrk(:,1),xlab,ylab,xlims,ylims,color,...
            file_tag,folder_name,tag,PlotErr);
    end


    % stress in loading direction from post-processing virial formulation
    figno = figno+1;
    figure(figno)
    x = time_ld/1e-6;
    y_tmp = 2*sig11/1e6;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$\sigma_{11}$ (MPa)';
    tag = 'sig_11_endtoend';
    ylims = [0 2.5];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 2;

    [~,mean_rel_err_sig11(iv),~,~,peak_abs_err_sig11(iv)] = ...
        PlotPolishedOutputLoadingRates(loading_rate,loading_rates,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);

    % Compute relaxation timescale
    figno = figno+1;
    figure(figno); hold on
    rlx_indx = find(y(:,2)==max(y(:,2)),1,'first');
    x = time_ld(rlx_indx:end)/(1e3*tau0);
    x = x-x(1);
    sig11_bead = y(rlx_indx:end,1)/y(rlx_indx,1);
%     err_bead = err(rlx_indx:end,1);
    sig11_meso = y(rlx_indx:end,2)/y(rlx_indx,2);
%     err_meso = err(rlx_indx:end,2);
    rng = (ceil(linspace(1,length(x),n_pts)))';

    % Fit the relaxation timescale

%     ft = fittype(['A*exp(-x*B) + (1-A)*exp(-x*',num2str(kd*1e3*tau0),')']);
    ft = fittype(['A*exp(-x^C*B) + (1-A)*exp(-x*',num2str(kd*1e3*tau0),')']);

    if iv==1
        C_max = 0.6;
    else
        C_max = 1;% B_bead(iv-1);
    end
    options = fitoptions('Method', 'NonlinearLeastSquares', ...
                     'Lower',[0 0 0.1], 'Upper',[1 Inf C_max]);

    sig11_ft = sig11_bead(~isnan(sig11_bead));
    x_ft = x(~isnan(sig11_bead));
    func_bead = fit(x_ft,sig11_ft,ft,options);


    sig11_ft = sig11_meso(~isnan(sig11_meso));
    x_ft = x(~isnan(sig11_meso));
    func_meso = fit(x_ft,sig11_ft,ft,options);

    val_bead = coeffvalues(func_bead);
    pm_bead = mean(abs(confint(func_bead)-val_bead),1);
    val_meso = coeffvalues(func_meso);
    pm_meso = mean(abs(confint(func_meso)-val_meso),1);
% 
%     f_bead(iv) = val_bead(1);
%     tau_bead(iv) = val_bead(2);
%     SEf_bead(iv) = pm_bead(1);
%     SEtau_bead(iv) = pm_bead(2);
%     f_meso(iv) = val_meso(1);
%     tau_meso(iv) = val_meso(2);
%     SEf_meso(iv) = pm_meso(1);
%     SEtau_meso(iv) = pm_meso(2);

    f_bead(iv) = val_bead(1);
    k_bead(iv) = val_bead(2)^(1/val_bead(3));
    B_bead(iv) = val_bead(3);
    SEf_bead(iv) = pm_bead(1);
    SEk_bead(iv) = pm_bead(2);
    SEB_bead(iv) = pm_bead(3);
    f_meso(iv) = val_meso(1);
    k_meso(iv) = val_meso(2)^(1/val_meso(3));
    B_meso(iv) = val_meso(3);
    SEf_meso(iv) = pm_meso(1);
    SEk_meso(iv) = pm_meso(2);
    SEB_meso(iv) = pm_meso(3);
    
    if ismember(loading_rate,param_rng)
        [~] = PlotErrorBar(x(rng),sig11_bead(rng),NaN*ones(size(rng)),color,'filled');
        e = PlotErrorBar(x(rng),sig11_meso(rng),NaN*ones(size(rng)),color,'');
        e.Marker = '^';
        p = plot(x,func_bead(x));
        p.LineWidth = 1;
        p.LineStyle = '-';
        p.Color = color;

        p = plot(x,func_meso(x));
        p.LineWidth = 1;
        p.LineStyle = '--';
        p.Color = color;
    end

    if loading_rate==loading_rates(end)
        set(gca,'FontSize',20/1.5)
        set(gcf,'color','w')
        set(gcf,'Position',[1000 100 400 400])
        xlabel('$t/(\tau_0 \cdot 10^3)$','FontSize',20,'Interpreter','latex')
        ylabel('$\sigma_{11}^*$','FontSize',20,'Interpreter','latex')
    end

    % Detachment rate
    figno = figno+1;
    figure(figno)
    x = time_eq/1e-6;
    y_tmp = kd_out*tau0;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$k_d\tau_0$';
    tag = 'kd';
    ylims = [0 Inf];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 2;

    [~,mean_rel_err_kd(iv),~,~,~] = ...
        PlotPolishedOutputLoadingRates(loading_rate,loading_rates,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);


    % Alignment tensor g11
    figno = figno+1;
    figure(figno)
    x = time_ld/1e-6;
    y_tmp = rxr11;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$g_{11}$';
    tag = 'g11';
    ylims = [-0.1 1];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 1;

    plot(x,0.33*ones(size(x)),'k--')

    [~,mean_rel_err_g11(iv),~,~,peak_abs_err_g11(iv)] = ...
        PlotPolishedOutputLoadingRates(loading_rate,loading_rates,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);


    % Alignment tensor g22
    figno = figno+1;
    figure(figno)
    x = time_ld/1e-6;
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$g_{22}$';
    tag = 'g22';
    ylims = [-0.1 1];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 2;

    plot(x,0.33*ones(size(x)),'k--')
    y_tmp = rxr22;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    [~,mean_rel_err_g22(iv),~,~,peak_abs_err_g22(iv)] = ...
        PlotPolishedOutputLoadingRates(loading_rate,loading_rates,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);


    % Alignemnt tensor g12
    figno = figno+1;
    figure(figno)
    x = time_ld/tau0;
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$g_{12}$';
    tag = 'g12';
    ylims = [-0.1 1];
    xlims = [0 Inf];
    PlotErr = 1;
    err_type = 1;

    y_tmp = rxr12;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    [~,mean_rel_err_g12(iv),~,~,peak_abs_err_g12(iv)] = ...
        PlotPolishedOutputLoadingRates(loading_rate,loading_rates,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type);
end

close all
figno = 0;
font_size = 20;
x = loading_rates;
xlab = '$\dot \varepsilon \tau_0$';
xsize = 225; 
ysize = 225;
fig_aspect = [1 1 1];

% outputs = {'sig11','g11','g22','ka','fa','nbonds','nother','nself','ca','MSDst','MSDth'}';
outputs = {'sig11_abs','sig11_rel','g11','g22','kd'}';
slope = zeros(length(outputs),1);
SEm = zeros(length(outputs),1);
intercept = zeros(length(outputs),1);
SEb = zeros(length(outputs),1);
R2 = zeros(length(outputs),1);
p_val = zeros(length(outputs),1);

figno = figno + 1;
y = peak_abs_err_sig11;
xlims = [0 0.11];
ylims = [-Inf Inf];
ylab = '$| Err. |_{max}$';
[~,f2_sig11,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','sig11_error.png'])
saveas(gcf,[folder_name,'/','sig11_error.fig'])
intercept(figno) = f2_sig11.Coefficients.Estimate(1);
slope(figno) = f2_sig11.Coefficients.Estimate(2);
SEm(figno) = f2_sig11.Coefficients.SE(2);
SEb(figno) = f2_sig11.Coefficients.SE(1);
R2(figno) = f2_sig11.Rsquared.Ordinary;
p_val(figno) = f2_sig11.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_sig11;
xlims = [0 0.11];
ylims = [-Inf Inf];
ylab = '$\langle \% Err. \rangle$';
[~,f2_sig11_rel,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','sig11_error.png'])
saveas(gcf,[folder_name,'/','sig11_error.fig'])
intercept(figno) = f2_sig11_rel.Coefficients.Estimate(1);
slope(figno) = f2_sig11_rel.Coefficients.Estimate(2);
SEm(figno) = f2_sig11_rel.Coefficients.SE(1);
SEb(figno) = f2_sig11_rel.Coefficients.SE(2);
R2(figno) = f2_sig11_rel.Rsquared.Ordinary;
p_val(figno) = f2_sig11_rel.anova.pValue(1);

figno = figno + 1;
y = peak_abs_err_g11;
xlims = [0 0.11];
ylims = [-Inf Inf];
ylab = '$| Err. |_{max}$';
[~,f2_g11,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','g11_error.png'])
saveas(gcf,[folder_name,'/','g11_error.fig'])
intercept(figno) = f2_g11.Coefficients.Estimate(1);
slope(figno) = f2_g11.Coefficients.Estimate(2);
SEm(figno) = f2_g11.Coefficients.SE(1);
SEb(figno) = f2_g11.Coefficients.SE(2);
R2(figno) = f2_g11.Rsquared.Ordinary;
p_val(figno) = f2_g11.anova.pValue(1);

figno = figno + 1;
y = peak_abs_err_g22;
xlims = [0 0.11];
ylims = [-Inf Inf];
ylab = '$| Err. |_{max}$';
[~,f2_g22,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','g22_error.png'])
saveas(gcf,[folder_name,'/','g22_error.fig'])
intercept(figno) = f2_g22.Coefficients.Estimate(1);
slope(figno) = f2_g22.Coefficients.Estimate(2);
SEm(figno) = f2_g22.Coefficients.SE(1);
SEb(figno) = f2_g22.Coefficients.SE(2);
R2(figno) = f2_g22.Rsquared.Ordinary;
p_val(figno) = f2_g22.anova.pValue(1);


figno = figno + 1;
y = mean_rel_err_g11;
xlims = [0 0.11];
ylims = [-Inf Inf];
ylab = '$\langle \% Err. \rangle$';
[~,f2_g11,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','g11_error.png'])
saveas(gcf,[folder_name,'/','g11_error.fig'])
% intercept(figno) = f2_g11.Coefficients.Estimate(1);
% slope(figno) = f2_g11.Coefficients.Estimate(2);
% SEm(figno) = f2_g11.Coefficients.SE(1);
% SEb(figno) = f2_g11.Coefficients.SE(2);
% R2(figno) = f2_g11.Rsquared.Ordinary;
% p_val(figno) = f2_g11.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_g22;
xlims = [0 0.11];
ylims = [-Inf Inf];
ylab = '$\langle \% Err. \rangle$';
[~,f2_g22,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','g22_error.png'])
saveas(gcf,[folder_name,'/','g22_error.fig'])
% intercept(figno) = f2_g22.Coefficients.Estimate(1);
% slope(figno) = f2_g22.Coefficients.Estimate(2);
% SEm(figno) = f2_g22.Coefficients.SE(1);
% SEb(figno) = f2_g22.Coefficients.SE(2);
% R2(figno) = f2_g22.Rsquared.Ordinary;
% p_val(figno) = f2_g22.anova.pValue(1);


if nansum(mean_rel_err_kd)~=0
    figno = figno + 1;
    y = mean_rel_err_kd;
    ylab = '$\langle \% Err. \rangle$';
    xlims = [0 0.11];
    ylims = [-Inf Inf];
    [~,f2_kd,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
        xlab,ylab,font_size,fig_aspect);
    saveas(gcf,[folder_name,'/','kd_error.png'])
    saveas(gcf,[folder_name,'/','kd_error.fig'])
    intercept(figno) = f2_kd.Coefficients.Estimate(1);
    slope(figno) = f2_kd.Coefficients.Estimate(2);
    SEm(figno) = f2_kd.Coefficients.SE(1);
    SEb(figno) = f2_kd.Coefficients.SE(2);
    R2(figno) = f2_kd.Rsquared.Ordinary;
    p_val(figno) = f2_kd.anova.pValue(1);
end

tab = table(outputs,slope,SEm,intercept,SEb,R2,p_val);
writetable(tab,[folder_name,'/Error Analysis Table.txt'])

tab = table(loading_rates,f_bead,SEf_bead,k_bead,SEk_bead,...
    f_meso,SEf_meso,k_meso,SEk_meso);
writetable(tab,[folder_name,'/Relaxation Timescale Analysis Table.txt'])

figno = figno + 1; 
c1 = [0.2 0 0.5];
c2 = [0 0.5 0];
c3 = [0.8 0.4 0.1];
figure(figno); clf; hold on
pb_asp = [1 1 1];
% ax1 = axes; 
yyaxis left
set(gca,'FontSize',20/1.5)
e = errorbar(loading_rates,f_bead,SEf_bead);
e.Marker = 'o';
e.LineWidth = 1.5;
% e.LineStyle = 'none';
e.Color = c1;
e.MarkerEdgeColor = c1;
% e.MarkerFaceColor = 'k';

e = errorbar(loading_rates,f_meso,SEf_meso);
e.Marker = '^';
e.LineWidth = 1.5;
% e.LineStyle = 'none';
e.Color = c1;
e.MarkerEdgeColor = c1;
% e.MarkerFaceColor = 'none';
% ylabel('$f_\alpha$','FontSize',20,'Interpreter','latex')
ylim([0 1])
set(gca, 'YColor',c1); % Set left y-axis color

yyaxis right
e = errorbar(loading_rates,k_bead,SEk_bead);
e.Marker = 'o';
e.LineWidth = 1.5;
% e.LineStyle = 'none';
e.Color = c2;
e.MarkerEdgeColor = c2;
% e.MarkerFaceColor = 'k';

e = errorbar(loading_rates,k_meso,SEk_meso);
e.Marker = '^';
e.LineWidth = 1.5;
% e.LineStyle = 'none';
e.Color = c2;
e.MarkerEdgeColor = c2;
% e.MarkerFaceColor = 'none';
ylabel('$k_\alpha\tau_0\cdot 10^3$','FontSize',20,'Interpreter','latex')
ylim([0 100])
set(gca, 'YColor',c2); % Set left y-axis color

xlabel('$\dot \varepsilon \tau_0$','FontSize',20,'Interpreter','latex')
set(gcf,'Color','w')
set(gcf,'Position',[1000 100 400 400])
set(gca,'xscale','log')
xticks(flipud(loading_rates))
xticklabels(round(flipud(loading_rates),2))
xlim([0.009 0.11])
pbaspect(pb_asp)

ax1 = gca;

% create 2nd, transparent axes
ax2 = axes('position', ax1.Position);
hold on

% c3 = [0.5 0.5 0.5];
e = errorbar(loading_rates,B_bead',SEB_bead');
e.Marker = 'o';
e.LineStyle = '-';
e.LineWidth = 1.5;
e.Color = c3;
e.MarkerEdgeColor = c3;

e = errorbar(loading_rates,B_meso',SEB_meso);
e.Marker = '^';
e.LineStyle = '--';
e.LineWidth = 1.5;
e.Color = c3;
e.MarkerEdgeColor = c3;


% plot(ax2,x,y3, 'k')
pause(0.1)                 % see [3]
ax2.Color = 'none'; 
% grid(ax2, 'on')
% Horizontally scale the y axis to alight the grid (again, be careful!)
ax2.XLim = ax1.XLim; 
ax2.XTick = ax1.XTick; 
set(ax2,'ylim',[0 1])
set(ax2,'xscale','log')
ax2.YColor = 'k';
% ax2.YTick = 100;
% horzontally offset y tick labels
ax2.YTick = [0 0.2 0.4 0.6 0.8 1];
% ax2.YTickLabel = strcat(ax2.YTickLabel, {'              '}); 
% ax2.YTickLabel = strcat(ax2.YTickLabel, {'       '}); 
% ax2.YTick = linspace(yl(1), yl(2), length(ytick));      % see [2]

set(gca,'FontSize',20/1.5)
% ylabel('$\chi$','FontSize',20,'Interpreter','latex','Color',c3)
xticks([])
xlim([0.009 0.11])
pbaspect(pb_asp)

figno = figno + 1;
xlims = [0 0.11];
ylims = [-Inf Inf];
ylab = '$| Err. |_{max}$';
[~,f2_g22,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','g22_error.png'])
saveas(gcf,[folder_name,'/','g22_error.fig'])
intercept(figno) = f2_g22.Coefficients.Estimate(1);
slope(figno) = f2_g22.Coefficients.Estimate(2);
SEm(figno) = f2_g22.Coefficients.SE(1);
SEb(figno) = f2_g22.Coefficients.SE(2);
R2(figno) = f2_g22.Rsquared.Ordinary;
p_val(figno) = f2_g22.anova.pValue(1);

% figure(figno+1); clf; hold on
% stretch_indx = find(stretch>stretch(1),1,'first');
% stretch_rng = (stretch_indx:length(time_ld))';
% time_temp = (time_ld(stretch_rng)-time_ld(stretch_indx))/1e-6;
% plot(time_temp,stretch(stretch_rng),'Color','k','LineWidth',1.5)
% ylim([0 4])
% plot([0 time_temp(end)],[1 1],'k--')
% plot([0 time_temp(end)],[3 3],'k--')
% ylabel('$\lambda$','FontSize',20,'Interpreter','latex')
% xlabel('')
% xticks([])
% xticklabels({})
% yticks([1 3])
% yticklabels({'1','3'})
% set(gca,'FontSize',20/1.5)
% set(gcf,'color','w')
% pbaspect([4 1 1])
% set(gcf,'Position',[1000 100 400 400])
% saveas(gcf,[folder_name,'/','stretch_vs_time.png'])
% saveas(gcf,[folder_name,'/','stretch_vs_time.fig'])

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotChainLengthSweep(phi,N_Kuhns,kd,Weissenberg,folder_name,n_pts,Controls)

global raw_data_folder_name tau0 b Nt sig11_lmp_temp

n_fig = 19;
InitiateTheFigures(n_fig,Controls)
xsize = 400;
ysize = 400;

N_Kuhns = [12 18 24 30 36]';
% param_rng = [12 18 36];
param_rng = [12 18 24 30 36];
c1 = [0 0.8 0];
c2 = [0 0 0];
n_colors = length(N_Kuhns);
colors = [(linspace(c1(1),c2(1),n_colors))',...
    (linspace(c1(2),c2(2),n_colors))',...
    (linspace(c1(3),c2(3),n_colors))'];

MSE_sig11_lmp = zeros(length(N_Kuhns),1);
MSE_sig11 = zeros(length(N_Kuhns),1);
mean_rel_err_sig11_lmp = zeros(length(N_Kuhns),1);
mean_rel_err_sig11 = zeros(length(N_Kuhns),1);
mean_rel_err_ka = zeros(length(N_Kuhns),1);
mean_rel_err_kd = zeros(length(N_Kuhns),1);
mean_rel_err_g11 = zeros(length(N_Kuhns),1);
mean_rel_err_g12 = zeros(length(N_Kuhns),1);
mean_rel_err_g22 = zeros(length(N_Kuhns),1);
mean_rel_err_fa = zeros(length(N_Kuhns),1);
mean_rel_err_fd = zeros(length(N_Kuhns),1);
mean_rel_err_na = zeros(length(N_Kuhns),1);
mean_rel_err_nbonds = zeros(length(N_Kuhns),1);
mean_rel_err_nself = zeros(length(N_Kuhns),1);
mean_rel_err_nother = zeros(length(N_Kuhns),1);
mean_rel_err_fother = zeros(length(N_Kuhns),1);
mean_rel_err_ca = zeros(length(N_Kuhns),1);
mean_rel_err_MSDst = zeros(length(N_Kuhns),1);
mean_rel_err_MSDth = zeros(length(N_Kuhns),1);

falpha_bead = zeros(length(N_Kuhns),1);
kalpha_bead = zeros(length(N_Kuhns),1);
SEfalpha_bead = zeros(length(N_Kuhns),1);
SEtau_bead = zeros(length(N_Kuhns),1);
falpha_meso = zeros(length(N_Kuhns),1);
kalpha_meso = zeros(length(N_Kuhns),1);
SEfalpha_meso = zeros(length(N_Kuhns),1);
SEtau_meso = zeros(length(N_Kuhns),1);
sig_peak_bead = zeros(length(N_Kuhns),1);
sig_peak_meso = zeros(length(N_Kuhns),1);

kreconfig_bead = zeros(length(N_Kuhns),1);
SEkreconfig_bead = zeros(length(N_Kuhns),1);
kreconfig_meso = zeros(length(N_Kuhns),1);
SEkreconfig_meso = zeros(length(N_Kuhns),1);

for iv=1:length(N_Kuhns)
    N_Kuhn = N_Kuhns(iv);
    color = colors(iv,:);
    file_tag = ['N_',num2str(N_Kuhn),...
        '.kd_',num2str(kd,'%.2e'),...
        '.W_',num2str(Weissenberg,'%.3f'),...
        '.phi_',num2str(phi,'%.3f')];

    % Import the Conslidated Data

    % Time and stretch
    file_name = [raw_data_folder_name,'/','timestretch.',file_tag,'.mat'];
    timestretch = load(file_name,'-mat');
    time_ld = timestretch.time_ld;
    time_eq = timestretch.time_eq;
    stretch = timestretch.stretch;

    % Stress
    file_name = [raw_data_folder_name,'/','stress.',file_tag,'.mat'];
    stress = load(file_name,'-mat');
    sig11 = stress.sig11;
%     sig22 = stress.sig22;
%     sig12 = stress.sig12;
    sig11_lmp = stress.sig11_lmp;
%     sig22_lmp = stress.sig22_lmp;
%     sig33_lmp = stress.sig33_lmp;

    % Bond kinetics and topological data
    file_name = [raw_data_folder_name,'/','kinetics.',file_tag,'.mat'];
    kinetics = load(file_name,'-mat');
    ka_out = kinetics.ka;
    kd_out = kinetics.kd;
    fa = kinetics.fa;
    fd = kinetics.fd;
    n_attach = kinetics.n_attach;
    n_detach = kinetics.n_detach;
    n_bonds = kinetics.n_bonds;
    n_self = kinetics.n_self;
    n_others = kinetics.n_others;
    cluster_coeff = kinetics.cluster_coeff;

    % Chain alignment
    file_name = [raw_data_folder_name,'/','covariance.',file_tag,'.mat'];
    covariance = load(file_name,'-mat');
    rxr11 = covariance.rxr11_all;
    rxr22 = covariance.rxr22_all;
    rxr12 = covariance.rxr12_all;

    % MSD
    file_name = [raw_data_folder_name,'/','msd.',file_tag,'.mat'];
    MSD = load(file_name,'-mat');
    msd_sticker = MSD.msd_st;
    msd_tether = MSD.msd_th;

    sig11_lmp_temp = (sig11_lmp-sig11_lmp(1,:,:))/1e6;
    sig11_lmp_temp = squeeze(nanmean(sig11_lmp_temp,2));


    % Plot outputs
    % stress in loading direction from LAMMPS virial formulation
%     figno = 1;
%     figure(figno)
%     x = time_ld/1e-6;
%     y_tmp = (sig11_lmp-sig11_lmp(1,:,:))/1e6;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$\sigma_{11}$ (MPa)';
%     tag = 'sig11_lmp';
%     ylims = [0 2.5];
%     xlims = [0 Inf];
%     PlotErr = 1;
% 
%     [MSE_sig11_lmp(iv),mean_rel_err_sig11_lmp(iv)] = ...
%         PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);


    % Stress in loading direction
    figno = 1;
    figure(figno)
    x = time_ld/1e-6;
    y_tmp = (sig11_lmp-sig11_lmp(1,:,:))/1e6;
    y_lmp = squeeze(nanmean(y_tmp,2));
    err_lmp = squeeze(nanstd(y_tmp,1,2));
    y_tmp = 2*sig11/1e6;
    y_ntwrk = squeeze(nanmean(y_tmp,2));
    err_ntwrk = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$\sigma_{11}$ (MPa)';
    tag = 'sig11';
    ylims = [0 2.5];
    xlims = [0 Inf];
    PlotErr = 1;

    if ismember(N_Kuhn,param_rng)
        PlotPolishedOutputLAMMPS(N_Kuhn,N_Kuhns,figno,x,y_lmp(:,1),err_lmp(:,1),...
            y_ntwrk(:,1),err_ntwrk(:,1),xlab,ylab,xlims,ylims,color,...
            file_tag,folder_name,tag,PlotErr);
    end

    % stress in loading direction from post-processing virial formulation
    figno = figno+1;
    figure(figno)
    x = time_ld/1e-6;
    y_tmp = 2*sig11/1e6;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$\sigma_{11}$ (MPa)';
    tag = 'sig_11_endtoend';
    ylims = [0 2.5];
    xlims = [0 Inf];
    PlotErr = 1;

    [MSE_sig11(iv),mean_rel_err_sig11(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);

    % Compute relaxation timescale
    figno = figno+1;
    figure(figno); hold on
    rlx_indx = find(y(:,2)==max(y(:,2)),1,'first');
    x = time_ld(rlx_indx:end)/(1e3*tau0);
    x = x-x(1);
    sig11_bead = y(rlx_indx:end,1)/y(rlx_indx,1);
    sig11_meso = y(rlx_indx:end,2)/y(rlx_indx,2);
    rng = (ceil(linspace(1,length(x),n_pts)))';

    % Fit the relaxation timescale
    %     ft = fittype('S*(f*exp(-x/B) + (1-f)*exp(-x/D))');
    %     ft = fittype('A*exp(-x/B) + (1-A)*exp(-x/C)');
    %     ft = fittype('A*exp(-x/B) + C*exp(-x/D)');
    %     ft = fittype(['A*(B*exp(-x/C) + (1-B)*exp(-x*',num2str(kd*1e-6),'))']);
    %     ft = fittype(['A*exp(-x/',num2str(tau0/1e-6),') + (1-A)*exp(-x/D)']);
    %     ft = fittype('A*exp(-x/B)');
    %     ft = fittype('(1-A)*exp(-x*B) + A*exp(-x/C)');
    %     if iv==1
    %         ft = fittype('A*exp(-x*B)+(1-A)');
    %         ft_meso = ft;
    %         ft_bead = ft;
    %         options = fitoptions('Method', 'NonlinearLeastSquares','Lower',0.1,'Upper',Inf);
    %     else
    % %         ft_meso = fittype(['A*exp(-x*',num2str(kalpha_meso),') + (1-A)*exp(-x*C)']);
    % %         ft_bead = fittype(['A*exp(-x*',num2str(kalpha_bead),') + (1-A)*exp(-x*C)']);
    ft = fittype(['A*exp(-x*B) + (1-A)*exp(-x*',num2str(kd*1e3*tau0),')']);
    options = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[1 Inf]);
    %     end
    %     ft = fittype(['(1-A)*exp(-x*B) + A*exp(-x*',num2str(kd*1e3*tau0),')']);
    sig11_ft = sig11_bead(~isnan(sig11_bead));
    x_ft = x(~isnan(sig11_bead));
    func_bead = fit(x_ft,sig11_ft,ft,options);


    sig11_ft = sig11_meso(~isnan(sig11_meso));
    x_ft = x(~isnan(sig11_meso));
    func_meso = fit(x_ft,sig11_ft,ft,options);

    val_bead = coeffvalues(func_bead);
    pm_bead = mean(abs(confint(func_bead)-val_bead),1);
    val_meso = coeffvalues(func_meso);
    pm_meso = mean(abs(confint(func_meso)-val_meso),1);

    %     if iv==1
    %         falpha_bead(iv) = val_bead(1);
    %         SEfalpha_bead(iv) = pm_bead(1);
    %         falpha_meso(iv) = val_meso(1);
    %         SEfalpha_meso(iv) = pm_meso(1);
    %         kalpha_bead = val_bead(2);
    %         SEkalpha_bead = pm_bead(2);
    %         kalpha_meso = val_meso(2);
    %         SEkalpha_meso = pm_meso(2);
    %         kreconfig_bead(iv) = 0;
    %         SEkreconfig_bead(iv) = 0;
    %         kreconfig_meso(iv) = 0;
    %         SEkreconfig_meso(iv) = 0;
    %     else
    falpha_bead(iv) = val_bead(1);
    SEfalpha_bead(iv) = pm_bead(1);
    falpha_meso(iv) = val_meso(1);
    SEfalpha_meso(iv) = pm_meso(1);
    kreconfig_bead(iv) = val_bead(2);
    SEkreconfig_bead(iv) = pm_bead(2);
    kreconfig_meso(iv) = val_meso(2);
    SEkreconfig_meso(iv) = pm_meso(2);
    sig_peak_bead(iv) = y(rlx_indx,1);
    sig_peak_meso(iv) = y(rlx_indx,2);
    %     end

    if ismember(N_Kuhn,param_rng)
        [~] = PlotErrorBar(x(rng),sig11_bead(rng),NaN*ones(size(rng)),color,'filled');
        e = PlotErrorBar(x(rng),sig11_meso(rng),NaN*ones(size(rng)),color,'');
        e.Marker = '^';
        p = plot(x,func_bead(x));
        p.LineWidth = 1;
        p.LineStyle = '-';
        p.Color = color;

        p = plot(x,func_meso(x));
        p.LineWidth = 1;
        p.LineStyle = '--';
        p.Color = color;
    end

    if N_Kuhn==N_Kuhns(end)
        set(gca,'FontSize',20/1.5)
        set(gcf,'color','w')
        set(gcf,'Position',[1000 100 400 400])
        xlabel('$t/(\tau_0 \cdot 10^3)$','FontSize',20,'Interpreter','latex')
        ylabel('$\sigma_{11}^*$','FontSize',20,'Interpreter','latex')
    end


    % Attachment rate
    figno = figno+1;
    figure(figno)
    x = time_eq/1e-6;
    y_tmp = ka_out*tau0;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$k_a\tau_0$';
    tag = 'ka';
    ylims = [0 Inf];
    xlims = [0 Inf];
    PlotErr = 1;

    [~,mean_rel_err_ka(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);


    % Detachment rate
    figno = figno+1;
    figure(figno)
    x = time_eq/1e-6;
    y_tmp = kd_out*tau0;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$k_d\tau_0$';
    tag = 'kd';
    ylims = [0 Inf];
    xlims = [0 Inf];
    PlotErr = 1;

    [~,mean_rel_err_kd(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);


    % Attached and detached chain fractions
    figno = figno+1;
    figure(figno)
    x = time_eq/tau0;
    y_tmp = fa;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$f_a$, $f_d$';
    tag = 'f';
    ylims = [0 Inf];
    xlims = [0 Inf];
    PlotErr = 1;

    [~,mean_rel_err_fa(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);

    y_tmp = fd;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    [~,mean_rel_err_fd(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);


    % Alignment tensor g11
    figno = figno+1;
    figure(figno)
    x = time_ld/1e-6;
    y_tmp = rxr11;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$g_{11}$';
    tag = 'g11';
    ylims = [-0.1 1];
    xlims = [0 Inf];
    PlotErr = 1;

    plot(x,0.33*ones(size(x)),'k--')

    [~,mean_rel_err_g11(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);

    % Alignment tensor g11 inset
    figno = figno+1;
    figure(figno)
    x = time_ld/1e-6;
    y_tmp = rxr11;
    y = squeeze(nanmean(y_tmp,2));
    err = NaN*ones(size(y));
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$g_{11}$';
    tag = 'g11';
    ylims = [0.5 0.7];
    xlims = [0.25 1.25];
    PlotErr = 1;

    plot(x,0.33*ones(size(x)),'k--')

    [~,mean_rel_err_g11(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);

    set(gcf,'Position',[1000 100 200 200])


    % Alignment tensor g12
    figno = figno+1;
    figure(figno)
    x = time_ld/tau0;
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$g_{22}$';
    tag = 'g22';
    ylims = [-0.1 1];
    xlims = [0 Inf];
    PlotErr = 1;

    plot(x,0.33*ones(size(x)),'k--')
    y_tmp = rxr22;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    [~,mean_rel_err_g22(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);


    % Alignemnt tensor g22
    figno = figno+1;
    figure(figno)
    x = time_ld/tau0;
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$g_{12}$';
    tag = 'g12';
    ylims = [-0.1 1];
    xlims = [0 Inf];
    PlotErr = 1;

    y_tmp = rxr12;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    [~,mean_rel_err_g12(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);

    % MSD of sticker
    figno = figno+1;
    figure(figno)
    x = time_eq/tau0;
    y_tmp = msd_sticker/b^2;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$\langle \Delta x^2 \rangle_{st}/b^2$';
    tag = 'msd_st';
    ylims = [0 Inf];
    xlims = [0 Inf];
    PlotErr = 1;

    [~,mean_rel_err_MSDst(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);

    % MSD of tether
    figno = figno+1;
    figure(figno)
    x = time_eq/tau0;
    y_tmp = msd_tether/b^2;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$\langle \Delta x^2 \rangle_{th}/b^2$';
    tag = 'msd_th';
    ylims = [0 Inf];
    xlims = [0 Inf];
    PlotErr = 1;

    [~,mean_rel_err_MSDth(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);

    % N attached
    figno = figno+1;
    figure(figno)
    x = time_eq/tau0;
    y_tmp = n_attach;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$\langle n_{attach} \rangle$';
    tag = 'n_attach';
    ylims = [0 Inf];
    xlims = [0 Inf];
    PlotErr = 0;

    [~,mean_rel_err_na(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);

    figno = figno+1;
    figure(figno)
    x = time_eq/tau0;
    y_tmp = n_detach;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$\langle n_{detach} \rangle$';
    tag = 'n_detach';
    ylims = [0 Inf];
    xlims = [0 Inf];
    PlotErr = 0;

    PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);

    % N detached
    figno = figno+1;
    figure(figno)
    x = time_eq/tau0;
    y_tmp = n_bonds;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$\langle n_{bonds} \rangle$';
    tag = 'n_bonds';
    ylims = [0 Nt];
    xlims = [0 Inf];
    PlotErr = 1;

    [~,mean_rel_err_nbonds(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);

    % N self-bonded
    figno = figno+1;
    figure(figno)
    x = time_eq/tau0;
    y_tmp = n_self;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$\bar z_{self}$';
%     ylab = '$\langle n_{self} \rangle$';
    tag = 'n_self';
    ylims = [0 Nt];
    xlims = [0 Inf];
    PlotErr = 1;

    [~,mean_rel_err_nself(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);

    % N other-bonded
    figno = figno+1;
    figure(figno)
    x = time_eq/tau0;
    y_tmp = n_others;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
%     ylab = '$\langle n_{others} \rangle$';
    ylab = '$\bar z_{others}$';
    tag = 'n_others';
    ylims = [0 Nt];
    xlims = [0 Inf];
    PlotErr = 1;

    [~,mean_rel_err_nother(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);

    % fraction of other-bonded
    figno = figno+1;
    figure(figno)
    x = time_eq/tau0;
    y_tmp = n_others./n_bonds;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$\bar f_{oth}$';
    tag = 'f_others';
    ylims = [0 Nt];
    xlims = [0 Inf];
    PlotErr = 1;

    [~,mean_rel_err_fother(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);
    
    % Clustering coefficient
    figno = figno+1;
    figure(figno)
    x = time_eq/tau0;
    y_tmp = cluster_coeff;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    y(isnan(y)) = 0;
    err(isnan(err)) = 0;
    rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
    xlab = '$t/\tau_0$';
    ylab = '$\langle c_{\alpha} \rangle$';
    tag = 'cluster_coeff';
    ylims = [0 0.3];
    xlims = [0 Inf];
    PlotErr = 1;

    [~,mean_rel_err_ca(iv)] = ...
        PlotPolishedOutputLength(N_Kuhn,N_Kuhns,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize);
end

close all
figno = 0;
font_size = 20;
x = N_Kuhns;
xlab = '$N$';
xsize = 225; 
ysize = 225;
fig_aspect = [1 1 1];

outputs = {'sig11','g11','g22','ka','fa','nbonds','nother','nself','ca','MSDst','MSDth'}';
slope = zeros(length(outputs),1);
SEm = zeros(length(outputs),1);
intercept = zeros(length(outputs),1);
SEb = zeros(length(outputs),1);
R2 = zeros(length(outputs),1);
p_val = zeros(length(outputs),1);

figno = figno + 1;
y = mean_rel_err_sig11;
ylab = '$\langle \% Err. \rangle$';
xlims = [0 40];
ylims = [-25 5];
[~,f2_sig11,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','sig11_error.png'])
saveas(gcf,[folder_name,'/','sig11_error.fig'])
intercept(figno) = f2_sig11.Coefficients.Estimate(1);
slope(figno) = f2_sig11.Coefficients.Estimate(2);
SEm(figno) = f2_sig11.Coefficients.SE(2);
SEb(figno) = f2_sig11.Coefficients.SE(1);
R2(figno) = f2_sig11.Rsquared.Ordinary;
p_val(figno) = f2_sig11.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_g11;
ylab = '$\langle \% Err. \rangle$';
xlims = [0 40];
ylims = [-10 10];
[~,f2_g11,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','g11_error.png'])
saveas(gcf,[folder_name,'/','g11_error.fig'])
intercept(figno) = f2_g11.Coefficients.Estimate(1);
slope(figno) = f2_g11.Coefficients.Estimate(2);
SEm(figno) = f2_g11.Coefficients.SE(1);
SEb(figno) = f2_g11.Coefficients.SE(2);
R2(figno) = f2_g11.Rsquared.Ordinary;
p_val(figno) = f2_g11.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_g22;
ylab = '$\langle \% Err. \rangle$';
xlims = [0 40];
ylims = [-10 10];
[~,f2_g22,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','g22_error.png'])
saveas(gcf,[folder_name,'/','g22_error.fig'])
intercept(figno) = f2_g22.Coefficients.Estimate(1);
slope(figno) = f2_g22.Coefficients.Estimate(2);
SEm(figno) = f2_g22.Coefficients.SE(1);
SEb(figno) = f2_g22.Coefficients.SE(2);
R2(figno) = f2_g22.Rsquared.Ordinary;
p_val(figno) = f2_g22.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_ka;
ylab = '$\langle \% Err. \rangle$';
xlims = [0 40];
ylims = [-2 2];
[~,f2_ka,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','ka_error.png'])
saveas(gcf,[folder_name,'/','ka_error.fig'])
intercept(figno) = f2_ka.Coefficients.Estimate(1);
slope(figno) = f2_ka.Coefficients.Estimate(2);
SEm(figno) = f2_ka.Coefficients.SE(1);
SEb(figno) = f2_ka.Coefficients.SE(2);
R2(figno) = f2_ka.Rsquared.Ordinary;
p_val(figno) = f2_ka.anova.pValue(1);

if nansum(mean_rel_err_kd)~=0
    figno = figno + 1;
    y = mean_rel_err_kd;
    ylab = '$\langle \% Err. \rangle$';
    xlims = [0 40];
    ylims = [-2 2];
    [~,f2_kd,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
        xlab,ylab,font_size,fig_aspect);
    saveas(gcf,[folder_name,'/','kd_error.png'])
    saveas(gcf,[folder_name,'/','kd_error.fig'])
    intercept(figno) = f2_kd.Coefficients.Estimate(1);
    slope(figno) = f2_kd.Coefficients.Estimate(2);
    SEm(figno) = f2_kd.Coefficients.SE(1);
    SEb(figno) = f2_kd.Coefficients.SE(2);
    R2(figno) = f2_kd.Rsquared.Ordinary;
    p_val(figno) = f2_kd.anova.pValue(1);
end

figno = figno + 1;
y = mean_rel_err_fa;
ylab = '$\langle \% Err. \rangle$';
xlims = [0 40];
ylims = [-50 50];
[~,f2_fa,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','fa_error.png'])
saveas(gcf,[folder_name,'/','fa_error.fig'])
intercept(figno) = f2_fa.Coefficients.Estimate(1);
slope(figno) = f2_fa.Coefficients.Estimate(2);
SEm(figno) = f2_fa.Coefficients.SE(1);
SEb(figno) = f2_fa.Coefficients.SE(2);
R2(figno) = f2_fa.Rsquared.Ordinary;
p_val(figno) = f2_fa.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_nbonds;
ylab = '$\langle \% Err. \rangle$';
xlims = [0 40];
ylims = [-100 100];
[~,f2_nbonds,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','nbonds_error.png'])
saveas(gcf,[folder_name,'/','nbonds_error.fig'])
intercept(figno) = f2_nbonds.Coefficients.Estimate(1);
slope(figno) = f2_nbonds.Coefficients.Estimate(2);
SEm(figno) = f2_nbonds.Coefficients.SE(1);
SEb(figno) = f2_nbonds.Coefficients.SE(2);
R2(figno) = f2_nbonds.Rsquared.Ordinary;
p_val(figno) = f2_nbonds.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_nother;
ylab = '$\langle \% Err. \rangle$';
xlims = [0 40];
ylims = [-25 25];
[~,f2_nother,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','nother_error.png'])
saveas(gcf,[folder_name,'/','nother_error.fig'])
intercept(figno) = f2_nother.Coefficients.Estimate(1);
slope(figno) = f2_nother.Coefficients.Estimate(2);
SEm(figno) = f2_nother.Coefficients.SE(1);
SEb(figno) = f2_nother.Coefficients.SE(2);
R2(figno) = f2_nother.Rsquared.Ordinary;
p_val(figno) = f2_nother.anova.pValue(1);


figno = figno + 1;
y = mean_rel_err_nself;
ylab = '$\langle \% Err. \rangle$';
xlims = [0 40];
ylims = [-100 100];
[~,f2_nself,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','nself_error.png'])
saveas(gcf,[folder_name,'/','nself_error.fig'])
intercept(figno) = f2_nself.Coefficients.Estimate(1);
slope(figno) = f2_nself.Coefficients.Estimate(2);
SEm(figno) = f2_nself.Coefficients.SE(1);
SEb(figno) = f2_nself.Coefficients.SE(2);
R2(figno) = f2_nself.Rsquared.Ordinary;
p_val(figno) = f2_nself.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_ca;
ylab = '$\langle \% Err. \rangle$';
xlims = [0 40];
ylims = [-100 100];
[~,f2_ca,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','ca_error.png'])
saveas(gcf,[folder_name,'/','ca_error.fig'])
intercept(figno) = f2_ca.Coefficients.Estimate(1);
slope(figno) = f2_ca.Coefficients.Estimate(2);
SEm(figno) = f2_ca.Coefficients.SE(1);
SEb(figno) = f2_ca.Coefficients.SE(2);
R2(figno) = f2_ca.Rsquared.Ordinary;
p_val(figno) = f2_ca.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_MSDst;
ylab = '$\langle \% Err. \rangle$';
xlims = [0 40];
ylims = [-50 50];
[~,f2_MSDst,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
saveas(gcf,[folder_name,'/','MSDst_error.png'])
saveas(gcf,[folder_name,'/','MSDth_error.fig'])
intercept(figno) = f2_MSDst.Coefficients.Estimate(1);
slope(figno) = f2_MSDst.Coefficients.Estimate(2);
SEm(figno) = f2_MSDst.Coefficients.SE(1);
SEb(figno) = f2_MSDst.Coefficients.SE(2);
R2(figno) = f2_MSDst.Rsquared.Ordinary;
p_val(figno) = f2_MSDst.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_MSDth;
ylab = '$\langle \% Err. \rangle$';
xlims = [0 40];
ylims = [-50 25];
[~,f2_MSDth,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
intercept(figno) = f2_MSDth.Coefficients.Estimate(1);
slope(figno) = f2_MSDth.Coefficients.Estimate(2);
SEm(figno) = f2_MSDth.Coefficients.SE(1);
SEb(figno) = f2_MSDth.Coefficients.SE(2);
R2(figno) = f2_MSDth.Rsquared.Ordinary;
p_val(figno) = f2_MSDth.anova.pValue(1);

figno = figno + 1;
y = mean_rel_err_fother;
ylab = '$\langle \% Err. \rangle$';
xlims = [0 40];
ylims = [-Inf Inf];
[~,f2_fother,~] = PlotMeanRelativeError(figno,x,y,xlims,ylims,xsize,ysize,...
    xlab,ylab,font_size,fig_aspect);
intercept(figno) = f2_fother.Coefficients.Estimate(1);
slope(figno) = f2_fother.Coefficients.Estimate(2);
SEm(figno) = f2_fother.Coefficients.SE(1);
SEb(figno) = f2_fother.Coefficients.SE(2);
R2(figno) = f2_fother.Rsquared.Ordinary;
p_val(figno) = f2_fother.anova.pValue(1);

tab = table(outputs,slope,SEm,intercept,SEb,R2,p_val);
writetable(tab,[folder_name,'/Error Analysis Table.txt'])

tab = table(N_Kuhns,falpha_bead,SEfalpha_bead,kreconfig_bead,SEkreconfig_bead,...
    falpha_meso,SEfalpha_meso,kreconfig_meso,SEkreconfig_meso);
writetable(tab,[folder_name,'/Relaxation Timescale Analysis Table.txt'])

close all

figno = 0;
figno = figno + 1; 
c1 = [0.2 0 0.5];
c2 = [0 0.5 0];
figure(figno); clf; hold on
yyaxis left
set(gca,'FontSize',20/1.5)
e = errorbar(N_Kuhns,falpha_bead,SEfalpha_bead);
e.Marker = 'o';
e.LineWidth = 1.5;
% e.LineStyle = 'none';
e.Color = c1;
e.MarkerEdgeColor = c1;
% e.MarkerFaceColor = 'k';

e = errorbar(N_Kuhns,falpha_meso,SEfalpha_meso);
e.Marker = '^';
e.LineWidth = 1.5;
% e.LineStyle = 'none';
e.Color = c1;
e.MarkerEdgeColor = c1;
% e.MarkerFaceColor = 'none';
ylabel('$f_\alpha$','FontSize',20,'Interpreter','latex')
ylim([0 1])
set(gca, 'YColor',c1); % Set left y-axis color

yyaxis right
e = errorbar(N_Kuhns,kreconfig_bead,SEkreconfig_bead);
e.Marker = 'o';
e.LineWidth = 1.5;
% e.LineStyle = 'none';
e.Color = c2;
e.MarkerEdgeColor = c2;
% e.MarkerFaceColor = 'k';

e = errorbar(N_Kuhns,kreconfig_meso,SEkreconfig_meso);
e.Marker = '^';
e.LineWidth = 1.5;
% e.LineStyle = 'none';
e.Color = c2;
e.MarkerEdgeColor = c2;
% e.MarkerFaceColor = 'none';
ylabel('$k_\alpha \tau_0\cdot 10^3$','FontSize',20,'Interpreter','latex')
ylim([0 40])
set(gca, 'YColor',c2); % Set left y-axis color

xlabel('$k_d \tau_0$','FontSize',20,'Interpreter','latex')
set(gcf,'Color','w')
set(gcf,'Position',[1000 100 400 400])
% set(gca,'xscale','log')
% xticks(flipud(N_Kuhns))
xlim([-Inf Inf])
xticks(N_Kuhns)
% xticklabels({'0','0.001','0.003','0.01','0.03','0.1'})
pbaspect([1 1 1])


% stored vs dissipated stress and alpha relaxation ate
figno = figno + 1; 
c1 = [0.2 0 0.5];
c2 = [0 0.5 0];
figure(figno); clf; hold on
yyaxis left
set(gca,'FontSize',20/1.5)
e = errorbar(N_Kuhns,falpha_bead.*sig_peak_bead,SEfalpha_bead.*sig_peak_bead);
e.Marker = 'o';
e.LineWidth = 1.5;
e.LineStyle = '--';
e.Color = c1;
e.MarkerEdgeColor = c1;
% e.MarkerFaceColor = 'k';

e = errorbar(N_Kuhns,falpha_meso.*sig_peak_meso,SEfalpha_meso.*sig_peak_meso);
e.Marker = '^';
e.LineWidth = 1.5;
e.LineStyle = '--';
e.Color = c1;
e.MarkerEdgeColor = c1;
% e.MarkerFaceColor = 'none';
ylabel('$f_\alpha$','FontSize',20,'Interpreter','latex')
ylim([0 1])
set(gca, 'YColor',c1); % Set left y-axis color

e = errorbar(N_Kuhns,(1-falpha_bead).*sig_peak_bead,SEfalpha_bead.*sig_peak_bead);
e.Marker = 'o';
e.LineWidth = 1.5;
e.LineStyle = '-';
e.Color = c1;
e.MarkerEdgeColor = c1;
% e.MarkerFaceColor = 'k';

e = errorbar(N_Kuhns,(1-falpha_meso).*sig_peak_meso,SEfalpha_meso.*sig_peak_meso);
e.Marker = '^';
e.LineWidth = 1.5;
e.LineStyle = '-';
e.Color = c1;
e.MarkerEdgeColor = c1;
% e.MarkerFaceColor = 'none';
ylabel('$\sigma_{11}$ (MPa)','FontSize',20,'Interpreter','latex')
ylim([0 1.5])
set(gca, 'YColor',c1); % Set left y-axis color

yyaxis right
e = errorbar(N_Kuhns,kreconfig_bead,SEkreconfig_bead);
e.Marker = 'o';
e.LineWidth = 1.5;
e.LineStyle = '-';
e.Color = c2;
e.MarkerEdgeColor = c2;
% e.MarkerFaceColor = 'k';

e = errorbar(N_Kuhns,kreconfig_meso,SEkreconfig_meso);
e.Marker = '^';
e.LineWidth = 1.5;
e.LineStyle = '-';
e.Color = c2;
e.MarkerEdgeColor = c2;
% e.MarkerFaceColor = 'none';
ylabel('$k_\alpha \tau_0\cdot 10^3$','FontSize',20,'Interpreter','latex')
ylim([0 40])
set(gca, 'YColor',c2); % Set left y-axis color

xlabel('$k_d \tau_0$','FontSize',20,'Interpreter','latex')
set(gcf,'Color','w')
set(gcf,'Position',[1000 100 400 400])
% set(gca,'xscale','log')
% xticks(flipud(N_Kuhns))
xlim([-Inf Inf])
xticks(N_Kuhns)
% xticklabels({'0','0.001','0.003','0.01','0.03','0.1'})
pbaspect([1 1 1])


figure(figno+1); clf; hold on
stretch_indx = find(stretch>stretch(1),1,'first');
stretch_rng = (stretch_indx:length(time_ld))';
time_temp = (time_ld(stretch_rng)-time_ld(stretch_indx))/1e-6;
plot(time_temp,stretch(stretch_rng),'Color','k','LineWidth',1.5)
ylim([0 4])
plot([0 time_temp(end)],[1 1],'k--')
plot([0 time_temp(end)],[3 3],'k--')
ylabel('$\lambda$','FontSize',20,'Interpreter','latex')
xlabel('')
xticks([])
xticklabels({})
yticks([1 3])
yticklabels({'1','3'})
set(gca,'FontSize',20/1.5)
set(gcf,'color','w')
pbaspect([4 1 1])
set(gcf,'Position',[1000 100 400 400])
saveas(gcf,[folder_name,'/','stretch_vs_time.png'])
saveas(gcf,[folder_name,'/','stretch_vs_time.fig'])


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [func,f2,R2] = PlotMeanRelativeError(figno,x,y,xlims,ylims,...
    xsize,ysize,xlab,ylab,font_size,fig_aspect)

figure(figno); clf; hold on 
plot(xlims,[0 0],'k--')
scatter(x,y,'k','filled')
set(gca,'FontSize',font_size/2)
xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
ylabel(ylab,'FontSize',font_size,'Interpreter','latex')
set(gcf,'color','w')
pbaspect(fig_aspect)
xlim(xlims)
ylim(ylims)

f = fittype('poly1');
func = fit(x,y,f);
plot(xlims,func(xlims),'color','r')

f2 = fitlm(x,y);
R2 = f2.Rsquared.Ordinary;

rng = linspace(xlims(1),xlims(2),100)';
p12 = predint(func,rng,0.95,'functional','off');
plot(rng,p12,'k--')
set(gcf,'Position',[1000 100 xsize ysize])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotTheLJStudyOutputs(phis,N_Kuhn,kd,Weissenberg,folder_name,n_pts,Controls)

n_fig = 15;
InitiateTheFigures(n_fig,Controls)

colors = [0      0       0;
          0.9    0.6     0];

for iv=1:length(phis)
    color = colors(iv,:);
    phi = phis(iv);
    file_tag = ['N_',num2str(N_Kuhn),...
        '.kd_',num2str(kd,'%.2e'),...
        '.W_',num2str(Weissenberg,'%.3f'),...
        '.phi_',num2str(phi,'%.3f')];

    % Import the Conslidated Data

    % Time and stretch
    raw_folder_name = 'Consolidated Data';
    [time_ld_withLJ,time_eq_withLJ,stretch_withLJ,...
        sig11_withLJ,sig22_withLJ,sig12_withLJ,...
        ka_out_withLJ,kd_out_withLJ,...
        fa_withLJ,fd_withLJ,...
        n_attach_withLJ,n_detach_withLJ,n_bonds_withLJ,n_self_withLJ,n_others_withLJ,...
        cluster_coeff_withLJ,...
        rxr11_withLJ,rxr22_withLJ,rxr12_withLJ,msd_sticker_withLJ,msd_tether_withLJ] = ...
        ImportThePlottingData(raw_folder_name,file_tag);

    raw_folder_name = 'Consolidated Data/LJ Study';
    [time_ld_woutLJ,time_eq_woutLJ,stretch_woutLJ,...
        sig11_woutLJ,sig22_woutLJ,sig12_woutLJ,...
        ka_out_woutLJ,kd_out_woutLJ,...
        fa_woutLJ,fd_woutLJ,...
        n_attach_woutLJ,n_detach_woutLJ,n_bonds_woutLJ,n_self_woutLJ,n_others_woutLJ,...
        cluster_coeff_woutLJ,...
        rxr11_woutLJ,rxr22_woutLJ,rxr12_woutLJ,msd_sticker_woutLJ,msd_tether_woutLJ] = ...
        ImportThePlottingData(raw_folder_name,file_tag);

    figno = 1;
    figure(figno)
    x_w = time_ld_withLJ/1e-6;
    y_tmp = sig11_withLJ/1e6;
    y_w = squeeze(nanmean(y_tmp,2));
    err_w = squeeze(nanstd(y_tmp,1,2));
    rng_w = (ceil(linspace(1,length(x_w),n_pts)))';

    x_wo = time_ld_woutLJ/1e-6;
    y_tmp = sig11_woutLJ/1e6;
    y_wo = squeeze(nanmean(y_tmp,2));
    err_wo = squeeze(nanstd(y_tmp,1,2));
    rng_wo = (ceil(linspace(1,length(x_wo),n_pts)))';

    xlab = '$t$ ($\mu$s)';
    ylab = '$\sigma_{11}$ (MPa)';
    tag = 'sig11';
    ylims = [0 10];
    xlims = [0 Inf];
    PlotErr = 1;

    PlotPolishedForLJStudy(phi,figno,x_w,y_w,err_w,rng_w,...
        x_wo,y_wo,err_wo,rng_wo,...
        xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr)


    figno = 2;
    figure(figno)
    x_w = time_ld_withLJ/1e-6;
    y_tmp = sig22_withLJ/1e6;
    y_w = squeeze(nanmean(y_tmp,2));
    err_w = squeeze(nanstd(y_tmp,1,2));
    rng_w = (ceil(linspace(1,length(x_w),n_pts)))';

    x_wo = time_ld_woutLJ/1e-6;
    y_tmp = sig22_woutLJ/1e6;
    y_wo = squeeze(nanmean(y_tmp,2));
    err_wo = squeeze(nanstd(y_tmp,1,2));
    rng_wo = (ceil(linspace(1,length(x_wo),n_pts)))';

    xlab = '$t$ ($\mu$s)';
    ylab = '$\sigma_{22}$ (MPa)';
    tag = 'sig22';
    ylims = [-Inf 0];
    xlims = [0 Inf];
    PlotErr = 1;

    PlotPolishedForLJStudy(phi,figno,x_w,y_w,err_w,rng_w,...
        x_wo,y_wo,err_wo,rng_wo,...
        xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr)

%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);

%     figno = 2;
%     figure(figno)
%     x = time_ld/1e-6;
%     y_tmp = sig22/1e6;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$\sigma_{22}$ (MPa)';
%     tag = 'sig22';
%     ylims = [-2.5 0];
%     xlims = [0 Inf];
%     PlotErr = 0;
% 
%     %     y_smth = movmean(y(:,3)-y(1,3),10);
%     %     p = plot(x,y_smth,'--');
%     %     p.Color = color;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
%     figno = 3;
%     figure(figno)
%     x = time_ld/1e-6;
%     y_tmp = sig12/1e6;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$\sigma_{12}$ (MPa)';
%     tag = 'sig12';
%     ylims = [-2.5 2.5];
%     xlims = [0 Inf];
%     PlotErr = 0;
% 
%     %     y_smth = movmean(y(:,3)-y(1,3),10);
%     %     p = plot(x,y_smth,'--');
%     %     p.Color = color;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
% 
%     figno = 4;
%     figure(figno)
%     x = time_eq/1e-6;
%     y_tmp = ka_out*tau0;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$k_a\tau_0$';
%     tag = 'ka';
%     ylims = [0 Inf];
%     xlims = [0 Inf];
%     PlotErr = 1;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
% 
%     figno = 5;
%     figure(figno)
%     x = time_eq/1e-6;
%     y_tmp = kd_out*tau0;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$k_d\tau_0$';
%     tag = 'kd';
%     ylims = [0 Inf];
%     xlims = [0 Inf];
%     PlotErr = 1;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
% 
%     figno = 6;
%     figure(figno)
%     x = time_eq/1e-6;
%     y_tmp = fa;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$f_a$, $f_d$';
%     tag = 'f';
%     ylims = [0 Inf];
%     xlims = [0 Inf];
%     PlotErr = 1;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
%     y_tmp = fd;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
% 
%     figno = 7;
%     figure(figno)
%     x = time_ld/1e-6;
%     y_tmp = rxr11;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$g_{11}$, $g_{22}$, $g_{12}$';
%     tag = 'g';
%     ylims = [-0.1 1];
%     xlims = [0 Inf];
%     PlotErr = 0;
% 
%     plot(x,0.33*ones(size(x)),'k--')
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
%     y_tmp = rxr22;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
%     y_tmp = rxr12;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
% 
%     figno = 8;
%     figure(figno)
%     x = time_eq/1e-6;
%     y_tmp = msd_sticker/b^2;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$\langle \Delta x^2 \rangle_{st}/b^2$';
%     tag = 'msd_st';
%     ylims = [0 Inf];
%     xlims = [0 Inf];
%     PlotErr = 1;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
%     figno = 9;
%     figure(figno)
%     x = time_eq/1e-6;
%     y_tmp = msd_tether/b^2;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$\langle \Delta x^2 \rangle_{th}/b^2$';
%     tag = 'msd_th';
%     ylims = [0 Inf];
%     xlims = [0 Inf];
%     PlotErr = 1;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
%     figno = 10;
%     figure(figno)
%     x = time_eq/1e-6;
%     y_tmp = n_attach;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$\langle n_{attach} \rangle$';
%     tag = 'n_attach';
%     ylims = [0 Inf];
%     xlims = [0 Inf];
%     PlotErr = 0;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
%     figno = 11;
%     figure(figno)
%     x = time_eq/1e-6;
%     y_tmp = n_detach;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$\langle n_{detach} \rangle$';
%     tag = 'n_detach';
%     ylims = [0 Inf];
%     xlims = [0 Inf];
%     PlotErr = 0;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
%     figno = 12;
%     figure(figno)
%     x = time_eq/1e-6;
%     y_tmp = n_bonds;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$\langle n_{bonds} \rangle$';
%     tag = 'n_bonds';
%     ylims = [0 Nt];
%     xlims = [0 Inf];
%     PlotErr = 0;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
%     figno = 13;
%     figure(figno)
%     x = time_eq/1e-6;
%     y_tmp = n_self;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$\langle n_{self} \rangle$';
%     tag = 'n_self';
%     ylims = [0 Nt];
%     xlims = [0 Inf];
%     PlotErr = 0;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
%     figno = 14;
%     figure(figno)
%     x = time_eq/1e-6;
%     y_tmp = n_others;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$\langle n_{others} \rangle$';
%     tag = 'n_others';
%     ylims = [0 Nt];
%     xlims = [0 Inf];
%     PlotErr = 0;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
%     figno = 15;
%     figure(figno)
%     x = time_eq/1e-6;
%     y_tmp = cluster_coeff;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     y(isnan(y)) = 0;
%     err(isnan(err)) = 0;
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$\langle c_{\alpha} \rangle$';
%     tag = 'cluster_coeff';
%     ylims = [0 0.3];
%     xlims = [0 Inf];
%     PlotErr = 1;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time_ld,time_eq,stretch,...
    sig11,sig22,sig12,...
    ka_out,kd_out,...
    fa,fd,...
    n_attach,n_detach,n_bonds,n_self,n_others,...
    cluster_coeff,...
    rxr11,rxr22,rxr12,msd_sticker,msd_tether] = ...
    ImportThePlottingData(folder_name,file_tag)

file_name = [folder_name,'/','timestretch.',file_tag,'.mat'];
timestretch = load(file_name,'-mat');
time_ld = timestretch.time_ld;
time_eq = timestretch.time_eq;
stretch = timestretch.stretch;

% Stress
file_name = [folder_name,'/','stress.',file_tag,'.mat'];
stress = load(file_name,'-mat');
sig11 = stress.sig11;
sig22 = stress.sig22;
sig12 = stress.sig12;

file_name = [folder_name,'/','kinetics.',file_tag,'.mat'];
kinetics = load(file_name,'-mat');
ka_out = kinetics.ka;
kd_out = kinetics.kd;
fa = kinetics.fa;
fd = kinetics.fd;
n_attach = kinetics.n_attach;
n_detach = kinetics.n_detach;
n_bonds = kinetics.n_bonds;
n_self = kinetics.n_self;
n_others = kinetics.n_others;
cluster_coeff = kinetics.cluster_coeff;

file_name = [folder_name,'/','covariance.',file_tag,'.mat'];
covariance = load(file_name,'-mat');
rxr11 = covariance.rxr11_all;
rxr22 = covariance.rxr22_all;
rxr12 = covariance.rxr12_all;

file_name = [folder_name,'/','msd.',file_tag,'.mat'];
MSD = load(file_name,'-mat');
msd_sticker = MSD.msd_st;
msd_tether = MSD.msd_th;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotTheFullParameterSweep(phis,N_Kuhn,kd,Weissenberg,folder_name,n_pts,Controls)

global raw_data_folder_name tau0 b Nt

n_fig = 12;
InitiateTheFigures(n_fig,Controls)

colors = [0      0      0;
         0.5     0      1];

% if Weissenberg~=0.01
%     n_pts = round(n_pts*1.2);
% end

for iv=1:length(phis)
    color = colors(iv,:);
    phi = phis(iv);
    file_tag = ['N_',num2str(N_Kuhn),...
        '.kd_',num2str(kd,'%.2e'),...
        '.W_',num2str(Weissenberg,'%.3f'),...
        '.phi_',num2str(phi,'%.3f')];

    
    % Import the Conslidated Data

    % Time and stretch
    file_name = [raw_data_folder_name,'/','timestretch.',file_tag,'.mat'];
    timestretch = load(file_name,'-mat');
    time_ld = timestretch.time_ld;
    time_all = timestretch.time_all;
    stretch = timestretch.stretch;
    eq_end = find(stretch~=1,1,'first');
    eq_rng = (1:eq_end)';


    % Import stress data
    file_name = [raw_data_folder_name,'/','stress.',file_tag,'.mat'];
    stress = load(file_name,'-mat');
    sig11 = stress.sig11;
    sig22 = stress.sig22;
    sig12 = stress.sig12;
    sig11_lmp = stress.sig11_lmp;
    sig22_lmp = stress.sig22_lmp;
    sig33_lmp = stress.sig33_lmp;

    % Import kinetics/topology data
    file_name = [raw_data_folder_name,'/','kinetics.',file_tag,'.mat'];
    kinetics = load(file_name,'-mat');
    ka_out = kinetics.ka;
    kd_out = kinetics.kd;
    fa = kinetics.fa;
    fd = kinetics.fd;
    n_attach = kinetics.n_attach;
    n_detach = kinetics.n_detach;
    n_bonds = kinetics.n_bonds;
    n_self = kinetics.n_self;
    n_others = kinetics.n_others;
    cluster_coeff = kinetics.cluster_coeff;

    % Import chain alignment data
    file_name = [raw_data_folder_name,'/','covariance.',file_tag,'.mat'];
    covariance = load(file_name,'-mat');
    rxr11 = covariance.rxr11_all;
    rxr22 = covariance.rxr22_all;
    rxr12 = covariance.rxr12_all;

    % Import MSD data
    file_name = [raw_data_folder_name,'/','msd.',file_tag,'.mat'];
    MSD = load(file_name,'-mat');
    msd_sticker = MSD.msd_st;
    msd_tether = MSD.msd_th;

    sig11_lmp_tmp = (sig11_lmp-sig11_lmp(1,:,:))/1e6;
    sig11_lmp_tmp = squeeze(nanmean(sig11_lmp_tmp,2));

    % Stress in loading direction
    figno = 1;
    figure(figno)
    x = time_ld/1e-6;
    y_tmp = (sig11_lmp-sig11_lmp(1,:,:))/1e6;
    y_lmp = squeeze(nanmean(y_tmp,2));
    err_lmp = squeeze(nanstd(y_tmp,1,2));
    y_tmp = 2*sig11/1e6;
    y_ntwrk = squeeze(nanmean(y_tmp,2));
    err_ntwrk = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$\sigma_{11}$ (MPa)';
    tag = 'sig11';
    ylims = [0 10];
    xlims = [0 Inf];
    PlotErr = 1;

    PlotPolishedOutputLAMMPS(phi,phis,figno,x,y_lmp(:,1),err_lmp(:,1),...
        y_ntwrk(:,1),err_ntwrk(:,1),xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr);

    % stress from end-to-end distribution
    figno = figno+1;
    figure(figno)
    x = time_ld/1e-6;
    y_tmp = 2*sig11/1e6;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    if Weissenberg~=0.01
        rng = [rng;find(sig11(:,1)==max(sig11(:,1)))];
        rng = sort(rng); 
    end
    xlab = '$t$ ($\mu$s)';
    ylab = '$\sigma_{11}$ (MPa)';
    tag = 'sig_11_endtoend';
    ylims = [0 2.5];
    xlims = [0 Inf];
    PlotErr = 1;

    PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,sig11_lmp_tmp);

%     figno = figno+1;
%     figure(figno)
%     x = time_ld/1e-6;
%     y_tmp = sig22/1e6;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$\sigma_{22}$ (MPa)';
%     tag = 'sig22';
%     ylims = [-2.5 0];
%     xlims = [0 Inf];
%     PlotErr = 0;
% 
%     %     y_smth = movmean(y(:,3)-y(1,3),10);
%     %     p = plot(x,y_smth,'--');
%     %     p.Color = color;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);
% 
%     figno = figno+1;
%     figure(figno)
%     x = time_ld/1e-6;
%     y_tmp = sig12/1e6;
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$\sigma_{12}$ (MPa)';
%     tag = 'sig12';
%     ylims = [-2.5 2.5];
%     xlims = [0 Inf];
%     PlotErr = 0;
% 
%     %     y_smth = movmean(y(:,3)-y(1,3),10);
%     %     p = plot(x,y_smth,'--');
%     %     p.Color = color;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr);


    figno = figno+1;
    figure(figno)
    x = time_all(eq_rng)/1e-6;
    y_tmp = ka_out(eq_rng,:,:)*tau0;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$k_a\tau_0$';
    tag = 'ka';
    ylims = [0 0.04];
    xlims = [0 Inf];
    PlotErr = 1;

    PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,sig11_lmp_tmp);
    figure(figno)
    set(gcf,'Position',[1000 100 425 300])


    figno = figno+1;
    figure(figno)
    x = time_all(eq_rng)/1e-6;
    y_tmp = kd_out(eq_rng,:,:)*tau0;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$k_d\tau_0$';
    tag = 'kd';
    ylims = [0 Inf];
    xlims = [0 Inf];
    PlotErr = 1;

    PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,sig11_lmp_tmp);
    figure(figno)
    set(gcf,'Position',[1000 100 425 300])


    figno = figno+1;
    figure(figno)
    x = time_all(eq_rng)/1e-6;
    y_tmp = fa(eq_rng,:,:);
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$f_a$, $f_d$';
    tag = 'f';
    ylims = [0 Inf];
    xlims = [0 Inf];
    PlotErr = 1;

    PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,sig11_lmp_tmp);

    y_tmp = fd(eq_rng,:,:);
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,sig11_lmp_tmp);


    figno = figno+1;
    figure(figno)
    x = time_ld/1e-6;
    y_tmp = rxr11;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$g_{11}$, $g_{22}$, $g_{12}$';
    tag = 'g';
    ylims = [-0.1 1];
    xlims = [0 Inf];
    PlotErr = 0;

    plot(x,0.33*ones(size(x)),'k--')

    PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,sig11_lmp_tmp);

    y_tmp = rxr22;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,sig11_lmp_tmp);

    y_tmp = rxr12;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));

    PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,sig11_lmp_tmp);
    pbaspect([1 1 1])


    figno = figno+1;
    figure(figno)
    x = time_all(eq_rng)/1e-6;
    y_tmp = msd_sticker(eq_rng,:,:)/b^2;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$\langle \Delta x^2 \rangle_{st}/b^2$';
    tag = 'msd_st';
    ylims = [0 45];
    xlims = [0 Inf];
    PlotErr = 1;

    PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,sig11_lmp_tmp);

    figno = figno+1;
    figure(figno)
    x = time_all(eq_rng,:,:)/1e-6;
    y_tmp = msd_tether(eq_rng,:,:)/b^2;
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$\langle \Delta x^2 \rangle_{th}/b^2$';
    tag = 'msd_th';
    ylims = [0 20];
    xlims = [0 Inf];
    PlotErr = 1;

    PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,sig11_lmp_tmp);

%     figno = figno+1;
%     figure(figno)
%     x = time_all(eq_rng)/1e-6;
%     y_tmp = n_attach(eq_rng,:,:);
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$\langle n_{attach} \rangle$';
%     tag = 'n_attach';
%     ylims = [0 Inf];
%     xlims = [0 Inf];
%     PlotErr = 0;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr,sig11_lmp_tmp);
% 
%     figno = figno+1;
%     figure(figno)
%     x = time_all(eq_rng)/1e-6;
%     y_tmp = n_detach(eq_rng,:,:);
%     y = squeeze(nanmean(y_tmp,2));
%     err = squeeze(nanstd(y_tmp,1,2));
%     rng = (ceil(linspace(1,length(x),n_pts)))';
%     xlab = '$t$ ($\mu$s)';
%     ylab = '$\langle n_{detach} \rangle$';
%     tag = 'n_detach';
%     ylims = [0 Inf];
%     xlims = [0 Inf];
%     PlotErr = 0;
% 
%     PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
%         file_tag,folder_name,tag,PlotErr,sig11_lmp_tmp);

    figno = figno+1;
    figure(figno)
    x = time_all(eq_rng)/1e-6;
    y_tmp = n_bonds(eq_rng,:,:);
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$\langle n_{bonds} \rangle$';
    tag = 'n_bonds';
    ylims = [0 Nt];
    xlims = [0 Inf];
    PlotErr = 0;

    PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,sig11_lmp_tmp);
    figure(figno)
    pbaspect([4 1 1])
    xlabel('')
    xticklabels({})

    figno = figno+1;
    figure(figno)
    x = time_all(eq_rng)/1e-6;
    y_tmp = n_self(eq_rng,:,:);
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$\langle n_{self} \rangle$';
    tag = 'n_self';
    ylims = [0 Nt];
    xlims = [0 Inf];
    PlotErr = 0;

    PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,sig11_lmp_tmp);
    figure(figno)
    pbaspect([4 1 1])
    xlabel('')
    xticklabels({})

    figno = figno+1;
    figure(figno)
    x = time_all(eq_rng)/1e-6;
    y_tmp = n_others(eq_rng,:,:);
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$\langle n_{others} \rangle$';
    tag = 'n_others';
    ylims = [0 Nt];
    xlims = [0 Inf];
    PlotErr = 0;

    PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,sig11_lmp_tmp);
    figure(figno)
    pbaspect([4 1 1])
    xlabel('')
    xticklabels({})

    figno = figno+1;
    figure(figno)
    x = time_all(eq_rng)/1e-6;
    y_tmp = cluster_coeff(eq_rng,:,:);
    y = squeeze(nanmean(y_tmp,2));
    err = squeeze(nanstd(y_tmp,1,2));
    y(isnan(y)) = 0;
    err(isnan(err)) = 0;
    rng = (ceil(linspace(1,length(x),n_pts)))';
    xlab = '$t$ ($\mu$s)';
    ylab = '$\langle c_{\alpha} \rangle$';
    tag = 'cluster_coeff';
    pbaspect([4 1 1])
    set(gcf,'Position',[400 100 410 300])
    ylims = [0 0.5];
    xlims = [0 Inf];
    PlotErr = 1;

    PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
        file_tag,folder_name,tag,PlotErr,sig11_lmp_tmp);
    figure(figno)
    pbaspect([3 1 1])
    xlabel('')
    xticklabels({})
    set(gcf,'Position',[1000 100 425 300])
    figure(figno+100)
    pbaspect([4 1 1])
    set(gcf,'Position',[1000 100 445 300])
   
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitiateTheFigures(n_fig,Controls)

for fn=1:n_fig
    figure(fn); clf; hold on; 
    if Controls.RunOscillatory~=1 && Controls.RunLargeDeformation~=1
        figure(fn+100); clf; hold on;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotPolishedForLJStudy(phi,figno,x_w,y_w,err_w,rng_w,...
    x_wo,y_wo,err_wo,rng_wo,...
    xlab,ylab,xlims,ylims,color,...
    file_tag,folder_name,tag,PlotErr)

global phis font_size line_width

figure(figno)

[~] = PlotErrorBar(x_w(rng_w),y_w(rng_w,1),err_w(rng_w,1),color,'filled');
[~,~] = PlotCurve(x_w,y_w(:,2),err_w(:,2),color,'-');

[~] = PlotErrorBar(x_wo(rng_wo),y_wo(rng_wo,1),err_wo(rng_wo,1),color,'empty');
% [~,~] = PlotCurve(x_wo,y_wo(:,2),err_wo(:,2),color,'--');
set(gca,'FontSize',font_size/1.5)

if phi==phis(end)
    xlim(xlims)
    ylim(ylims)
    xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
    ylabel(ylab,'FontSize',font_size,'Interpreter','latex')
    pbaspect([1 1 1])
%     set(gcf,'Position',[750 400 400 400])
    set(gcf,'color','w')
    if PlotErr==1
        xlabel('')
        xticklabels({})
    end
    saveas(gcf,[folder_name,'/',tag,'.',file_tag,'.png'])
    saveas(gcf,[folder_name,'/',tag,'.',file_tag,'.fig'])
end

% plot error
if PlotErr==1
    figure(figno+100)
    plot([0 max(x_w)],[0 0],'k--')
    rel_err = (y_w(:,1)-y_w(:,2))./max(y_w(:,1))*100;
    p = plot(x_w,rel_err);
    p.Color = color;
    p.LineWidth = line_width;

    rel_err = (y_wo(:,1)-y_wo(:,2))./max(y_wo(:,1))*100;
    p = plot(x_wo,rel_err);
    p.Color = color;
    p.LineStyle = '--';
    p.LineWidth = line_width;

    set(gca,'FontSize',font_size/1.5)

    % if N_Kuhn==N_Kuhns(end)
    if phi==phis(end)
        xlim(xlims)
        ylim([-100 100])
        xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
        ylabel('$\%$ Err.','FontSize',font_size,'Interpreter','latex')
        pbaspect([4 1 1])
%         set(gcf,'Position',[750 100 400 225])
        set(gcf,'color','w')
        saveas(gcf,[folder_name,'/rel_error.',tag,'.',file_tag,'.png'])
        saveas(gcf,[folder_name,'/rel_error.',tag,'.',file_tag,'.fig'])
    end
else
    if ishandle(figno+100)
        close(figno+100)
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = ...
    PlotPolishedOutputCreep(param_temp,param_temps,figno,...
    x,y,err,xlab,ylab,xlims,ylims,color,...
    file_tag,folder_name,tag,param_rng,xsize)

global font_size 

if ismember(param_temp,param_rng)
    figure(figno)

    err(isnan(y)) = [];
    x(isnan(y)) = [];
    y(isnan(y)) = [];

    if ~isempty(x) && ~isempty(y) && ~isempty(err)
        [p,~] = PlotCurve(x,y,err,color,'-');
    else
        p = [];
    end

    set(gca,'FontSize',font_size/1.5)
end

if param_temp==param_temps(end)
    xlim(xlims)
    ylim(ylims)
    xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
    ylabel(ylab,'FontSize',font_size,'Interpreter','latex')
    pbaspect([1 1 1])
    set(gcf,'Position',[750 400 xsize xsize])
    set(gcf,'color','w')
    saveas(gcf,[folder_name,'/',tag,'.',file_tag,'.png'])
    saveas(gcf,[folder_name,'/',tag,'.',file_tag,'.fig'])
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MSE,mean_rel_err,peak_rel_err,mean_abs_err,peak_abs_err] = ...
    PlotPolishedOutputDetachmentRates(param_temp,param_temps,figno,...
    x,y,err,rng,xlab,ylab,xlims,ylims,color,...
    file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type)

global font_size line_width

if ismember(param_temp,param_rng)
    figure(figno)

    x1 = x;

    y1 = y(:,1);
    err1 = err(:,1);

    err1(isnan(y1)) = [];
    x1(isnan(y1)) = [];
    y1(isnan(y1)) = [];

    rng = rng(rng<=length(x1));
    if contains(tag,'sig11')
        [~,~] = PlotCurve(x1,y1,NaN*ones(size(err1)),color,'-');
    else
        [~] = PlotErrorBar(x1(rng),y1(rng),err1(rng),color,'filled');
    end

    x2 = x;

    y2 = y(:,2);
    err2 = err(:,2);

    err2(isnan(y2)) = [];
    x2(isnan(y2)) = [];
    y2(isnan(y2)) = [];

    if contains(tag,'sig11')
        [~,~] = PlotCurve(x2,y2,err2,color,'--');
    else
        [~,~] = PlotCurve(x2,y2,err2,color,'-');
    end

    % Plot also the virial stress measured directly from LAMMPS
%     if contains(tag,'sig_11_end')
%         p = plot(x2,sig11_lmp_temp(:,1),'LineStyle','-.');
%         p.Color = color;
%         p.LineWidth = 0.5;
%     end

    set(gca,'FontSize',font_size/1.5)
end

if param_temp==param_temps(end)
    xlim(xlims)
    ylim(ylims)
    xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
    ylabel(ylab,'FontSize',font_size,'Interpreter','latex')
    pbaspect([1 1 1])
%     pbaspect([2 1 1])
    set(gcf,'Position',[750 400 xsize xsize])
    set(gcf,'color','w')
    if PlotErr==1
        xlabel('')
        xticklabels({})
    end
    saveas(gcf,[folder_name,'/',tag,'.',file_tag,'.png'])
    saveas(gcf,[folder_name,'/',tag,'.',file_tag,'.fig'])
end

% plot error
if PlotErr==1 && ismember(param_temp,param_rng)
    figure(figno+100)
    set(gca,'FontSize',font_size/1.5)

    if err_type==0 % Bland_Altman
        x = mean(y(:,1:2),2);
        indices = randi(length(x),100,1);
        x = x(indices);
        difference = y(indices,2)-y(indices,1);
        val1 = nanmean(difference);
        val2 = val1+nanstd(difference);
        val3 = val1-nanstd(difference);
        p = plot([0 max(x)],[val1 val1],'-');
        p.Color = color;
        p = plot([0 max(x)],[val2 val2],'--');
        p.Color = color;
        p = plot([0 max(x)],[val3 val3],'--');
        p.Color = color;
        err = difference;
        s = scatter(x,err);
        s.MarkerFaceColor = color;
        s.MarkerEdgeColor = 'k';
        xlabel('$\bar \Delta$','FontSize',font_size,'Interpreter','latex')
        ylabel('$\Delta$','FontSize',font_size,'Interpreter','latex')
        xlim([0 Inf])
        save_err_tag = 'Bland-Altman';
    elseif err_type==1 % Relative error
        plot([0 max(x)],[0 0],'k--')
        rel_err = (y(:,2)-y(:,1))./max(y(:,1))*100;
        err = rel_err;
        p = plot(x,err);
        p.Color = color;
        p.LineWidth = line_width;
        ylim([-100 100])
        xlim(xlims)
        xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
        ylabel('$\%$ Err.','FontSize',font_size,'Interpreter','latex')
        save_err_tag = 'rel_err';
    elseif err_type==2 % Absolute error
        plot([0 max(x)],[0 0],'k--')
        abs_err = y(:,2)-y(:,1);
        err = abs_err;
        p = plot(x,err);
        p.Color = color;
        p.LineWidth = line_width;
        ylim([-Inf Inf])
        xlim(xlims)
        xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
        ylabel('Abs. Err.','FontSize',font_size,'Interpreter','latex')
        save_err_tag = 'rel_err';
    end

    % if N_Kuhn==N_Kuhns(end)
    if param_temp==param_temps(end)
        pbaspect([1 1 1])
%         pbaspect([4 1 1])
        set(gcf,'Position',[750 100 xsize ysize])
        set(gcf,'color','w')
        saveas(gcf,[folder_name,'/',save_err_tag,'.',tag,'.',file_tag,'.png'])
        saveas(gcf,[folder_name,'/',save_err_tag,'.',file_tag,'.fig'])
    end
% else
%     if ishandle(figno+100)
%         close(figno+100)
%     end
end

MSE = nanmean((y(:,1)-y(:,2)).^2);
rel_err = (y(:,2)-y(:,1))./y(:,2)*100;
rel_err(isinf(rel_err)) = 0;
mean_rel_err = nanmean(rel_err);
peak_rel_err = rel_err(abs(rel_err)==max(abs(rel_err)));
abs_err = y(:,2)-y(:,1);
mean_abs_err = nanmean(abs_err);
peak_abs_err = abs_err(abs(abs_err)==max(abs(abs_err)));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MSE,mean_rel_err,peak_rel_err,mean_abs_err,peak_abs_err] = ...
    PlotPolishedOutputLoadingRates(param_temp,param_temps,figno,...
    x,y,err,rng,xlab,ylab,xlims,ylims,color,...
    file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize,err_type)

global font_size line_width

if ismember(param_temp,param_rng)
    figure(figno)

    x1 = x;

    y1 = y(:,1);
    err1 = err(:,1);

    err1(isnan(y1)) = [];
    x1(isnan(y1)) = [];
    y1(isnan(y1)) = [];

    rng = rng(rng<=length(x1));
    if contains(tag,'sig11')
        [~,~] = PlotCurve(x1,y1,NaN*ones(size(err1)),color,'-');
    else
        [~] = PlotErrorBar(x1(rng),y1(rng),err1(rng),color,'filled');
    end

    x2 = x;

    y2 = y(:,2);
    err2 = err(:,2);

    err2(isnan(y2)) = [];
    x2(isnan(y2)) = [];
    y2(isnan(y2)) = [];

    if contains(tag,'sig11')
        [~,~] = PlotCurve(x2,y2,err2,color,'--');
    else
        [~,~] = PlotCurve(x2,y2,err2,color,'-');
    end

    % Plot also the virial stress measured directly from LAMMPS
%     if contains(tag,'sig_11_end')
%         p = plot(x2,sig11_lmp_temp(:,1),'LineStyle','-.');
%         p.Color = color;
%         p.LineWidth = 0.5;
%     end

    set(gca,'FontSize',font_size/1.5)
end

if param_temp==param_temps(end)
    xlim(xlims)
    ylim(ylims)
    xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
    ylabel(ylab,'FontSize',font_size,'Interpreter','latex')
    pbaspect([1 1 1])
%     pbaspect([2 1 1])
    set(gcf,'Position',[750 400 xsize xsize])
    set(gcf,'color','w')
    if PlotErr==1
        xlabel('')
        xticklabels({})
    end
    saveas(gcf,[folder_name,'/',tag,'.',file_tag,'.png'])
    saveas(gcf,[folder_name,'/',tag,'.',file_tag,'.fig'])
end

% plot error
if PlotErr==1 && ismember(param_temp,param_rng)
    figure(figno+100)
    set(gca,'FontSize',font_size/1.5)

    if err_type==0 % Bland_Altman
        x = mean(y(:,1:2),2);
        indices = randi(length(x),100,1);
        x = x(indices);
        difference = y(indices,2)-y(indices,1);
        val1 = nanmean(difference);
        val2 = val1+nanstd(difference);
        val3 = val1-nanstd(difference);
        p = plot([0 max(x)],[val1 val1],'-');
        p.Color = color;
        p = plot([0 max(x)],[val2 val2],'--');
        p.Color = color;
        p = plot([0 max(x)],[val3 val3],'--');
        p.Color = color;
        err = difference;
        s = scatter(x,err);
        s.MarkerFaceColor = color;
        s.MarkerEdgeColor = 'k';
        xlabel('$\bar \Delta$','FontSize',font_size,'Interpreter','latex')
        ylabel('$\Delta$','FontSize',font_size,'Interpreter','latex')
        xlim([0 Inf])
        save_err_tag = 'Bland-Altman';
    elseif err_type==1 % Relative error with respect to peak error
        plot([0 max(x)],[0 0],'k--')
        rel_err = (y(:,2)-y(:,1))./max(y(:,1))*100;
        err = rel_err;
        p = plot(x,err);
        p.Color = color;
        p.LineWidth = line_width;
        ylim([-100 100])
        xlim(xlims)
        xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
        ylabel('$\%$ Err.','FontSize',font_size,'Interpreter','latex')
        save_err_tag = 'rel_err_peak';
    elseif err_type==2 % Absolute error
        plot([0 max(x)],[0 0],'k--')
        abs_err = y(:,2)-y(:,1);
        err = abs_err;
        p = plot(x,err);
        p.Color = color;
        p.LineWidth = line_width;
        ylim([-Inf Inf])
        xlim(xlims)
        xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
        ylabel('Abs. Err.','FontSize',font_size,'Interpreter','latex')
        save_err_tag = 'rel_err';
    elseif err_type==3 % Relative error
        plot([0 max(x)],[0 0],'k--')
        rel_err = (y(:,2)-y(:,1))./y(:,1)*100;
        err = rel_err;
        p = plot(x,err);
        p.Color = color;
        p.LineWidth = line_width;
        ylim([-100 100])
        xlim(xlims)
        xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
        ylabel('$\%$ Err.','FontSize',font_size,'Interpreter','latex')
        save_err_tag = 'rel_err';
    end

    % if N_Kuhn==N_Kuhns(end)
    if param_temp==param_temps(end)
        pbaspect([1 1 1])
%         pbaspect([4 1 1])
        set(gcf,'Position',[750 100 xsize ysize])
        set(gcf,'color','w')
        saveas(gcf,[folder_name,'/',save_err_tag,'.',tag,'.',file_tag,'.png'])
        saveas(gcf,[folder_name,'/',save_err_tag,'.',file_tag,'.fig'])
    end
% else
%     if ishandle(figno+100)
%         close(figno+100)
%     end
end

MSE = nanmean((y(:,1)-y(:,2)).^2);
rel_err = (y(:,2)-y(:,1))./y(:,2)*100;
rel_err(isinf(rel_err)) = 0;
mean_rel_err = nanmean(rel_err);
peak_rel_err = rel_err(abs(rel_err)==max(abs(rel_err)));
abs_err = y(:,2)-y(:,1);
mean_abs_err = nanmean(abs_err);
peak_abs_err = abs_err(abs(abs_err)==max(abs(abs_err)));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MSE,mean_rel_err] = PlotPolishedOutputLength(param_temp,param_temps,figno,...
    x,y,err,rng,xlab,ylab,xlims,ylims,color,...
    file_tag,folder_name,tag,PlotErr,param_rng,xsize,ysize)

global font_size line_width sig11_lmp_temp

if ismember(param_temp,param_rng)
    figure(figno)

    x1 = x;

    y1 = y(:,1);
    err1 = err(:,1);

    err1(isnan(y1)) = [];
    x1(isnan(y1)) = [];
    y1(isnan(y1)) = [];

    rng = rng(rng<=length(x1));
    if contains(tag,'sig11')
        [~,~] = PlotCurve(x1,y1,NaN*ones(size(err1)),color,'-');
    else
        [~] = PlotErrorBar(x1(rng),y1(rng),err1(rng),color,'filled');
    end

    x2 = x;

    y2 = y(:,2);
    err2 = err(:,2);

    err2(isnan(y2)) = [];
    x2(isnan(y2)) = [];
    y2(isnan(y2)) = [];

    if contains(tag,'sig11')
        [~,~] = PlotCurve(x2,y2,err2,color,'--');
    else
        [~,~] = PlotCurve(x2,y2,err2,color,'-');
    end

    % Plot also the virial stress measured directly from LAMMPS
%     if contains(tag,'sig_11_end')
%         p = plot(x2,sig11_lmp_temp(:,1),'LineStyle','-.');
%         p.Color = color;
%         p.LineWidth = 0.5;
%     end

    set(gca,'FontSize',font_size/1.5)
end

if param_temp==param_temps(end)
    xlim(xlims)
    ylim(ylims)
    xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
    ylabel(ylab,'FontSize',font_size,'Interpreter','latex')
    pbaspect([1 1 1])
%     pbaspect([2 1 1])
    set(gcf,'Position',[750 400 xsize xsize])
    set(gcf,'color','w')
    if PlotErr==1
        xlabel('')
        xticklabels({})
    end
    saveas(gcf,[folder_name,'/',tag,'.',file_tag,'.png'])
    saveas(gcf,[folder_name,'/',tag,'.',file_tag,'.fig'])
end

% plot error

if PlotErr==1 && ismember(param_temp,param_rng)
    figure(figno+100)
    set(gca,'FontSize',font_size/1.5)

    err_type = 1;
    if err_type==0 % Bland_Altman
        x = mean(y(:,1:2),2);
        indices = randi(length(x),100,1);
        x = x(indices);
        difference = y(indices,2)-y(indices,1);
        val1 = nanmean(difference);
        val2 = val1+nanstd(difference);
        val3 = val1-nanstd(difference);
%         if ismember(param_temp,param_rng)
            p = plot([0 max(x)],[val1 val1],'-');
            p.Color = color;
            p = plot([0 max(x)],[val2 val2],'--');
            p.Color = color;
            p = plot([0 max(x)],[val3 val3],'--');
            p.Color = color;
            err = difference;
            s = scatter(x,err);
            s.MarkerFaceColor = color;
            s.MarkerEdgeColor = 'k';
            xlabel('$\bar \Delta$','FontSize',font_size,'Interpreter','latex')
            ylabel('$\Delta$','FontSize',font_size,'Interpreter','latex')
            xlim([0 Inf])
%         end
        save_err_tag = 'Bland-Altman';
    elseif err_type==1 % Relative error
        plot([0 max(x)],[0 0],'k--')
        rel_err = (y(:,2)-y(:,1))./max(y(:,1))*100;
        err = rel_err;
%         if ismember(param_temp,param_rng)
            p = plot(x,err);
            p.Color = color;
            p.LineWidth = line_width;
            ylim([-100 100])
            xlim(xlims)
            xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
            ylabel('$\%$ Err.','FontSize',font_size,'Interpreter','latex')
%         end
        save_err_tag = 'rel_err';
    end

    % if N_Kuhn==N_Kuhns(end)
    if param_temp==param_temps(end) 
        pbaspect([1 1 1])
%         pbaspect([4 1 1])
        set(gcf,'Position',[750 100 xsize ysize])
        set(gcf,'color','w')
        saveas(gcf,[folder_name,'/',save_err_tag,'.',tag,'.',file_tag,'.png'])
        saveas(gcf,[folder_name,'/',save_err_tag,'.',file_tag,'.fig'])
    end

    MSE = nanmean((y(:,1)-y(:,2)).^2);
    mean_rel_err = nanmean(rel_err);
else
%     if ishandle(figno+100)
%         close(figno+100)
%     end
    MSE = NaN;
    mean_rel_err = NaN;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotPolishedOutputLAMMPS(var,vars,figno,x,y_lmp,err_lmp,y_ntwk,err_ntwk,...
    xlab,ylab,xlims,ylims,color,...
    file_tag,folder_name,tag,PlotErr)

global font_size line_width


figure(figno)

x1 = x;

y1 = y_lmp;
err1 = err_lmp;

err1(isnan(y1)) = [];
x1(isnan(y1)) = [];
y1(isnan(y1)) = [];

[~,~] = PlotCurve(x1,y1,NaN*ones(size(err1)),color,'-');

x2 = x;

y2 = y_ntwk;
err2 = err_ntwk;

err2(isnan(y2)) = [];
x2(isnan(y2)) = [];
y2(isnan(y2)) = [];

[~,~] = PlotCurve(x2,y2,err2,color,'--');


set(gca,'FontSize',font_size/1.5)

if var==vars(end)
    xlim(xlims)
    ylim(ylims)
    xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
    ylabel(ylab,'FontSize',font_size,'Interpreter','latex')
    pbaspect([1 1 1])
    set(gcf,'Position',[750 400 400 400])
    set(gcf,'color','w')

    saveas(gcf,[folder_name,'/',tag,'.',file_tag,'.png'])
    saveas(gcf,[folder_name,'/',tag,'.',file_tag,'.fig'])
end

% plot error
if PlotErr==1
    figure(figno+100)

    err_type = 1;
    set(gca,'FontSize',font_size/1.5)
    if err_type==0 % Bland_Altman
        x = mean(y(:,1:2),2);
        indices = randi(length(x),100,1);
        x = x(indices);
        difference = y(indices,2)-y(indices,1);
        val1 = nanmean(difference);
        val2 = val1+nanstd(difference);
        val3 = val1-nanstd(difference);
        p = plot([0 max(x)],[val1 val1],'-');
        p.Color = color;
        p = plot([0 max(x)],[val2 val2],'--');
        p.Color = color;
        p = plot([0 max(x)],[val3 val3],'--');
        p.Color = color;
        err = difference;
        s = scatter(x,err);
        s.MarkerFaceColor = color;
        s.MarkerEdgeColor = 'k';
        xlabel('$\bar \Delta$','FontSize',font_size,'Interpreter','latex')
        ylabel('$\Delta$','FontSize',font_size,'Interpreter','latex')
        xlim([0 Inf])
        save_err_tag = 'Bland-Altman';
    elseif err_type==1 % Relative error
        plot([0 max(x)],[0 0],'k--')
        rel_err = (y_ntwk-y_lmp)./max(y_ntwk)*100;
        err = rel_err;
        p = plot(x,err);
        p.Color = color;
        p.LineWidth = line_width;
        ylim([-100 100])
        xlim(xlims)
        xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
        ylabel('$\%$ Err.','FontSize',font_size,'Interpreter','latex')
        save_err_tag = 'rel_err';
    end

    % if N_Kuhn==N_Kuhns(end)
    if var==vars(end)
        pbaspect([3 1 1])
        set(gcf,'Position',[750 100 425 200])
        set(gcf,'color','w')
        saveas(gcf,[folder_name,'/',save_err_tag,'.',tag,'.',file_tag,'.png'])
        saveas(gcf,[folder_name,'/',save_err_tag,'.',tag,'.',file_tag,'.fig'])
    end
else
    if ishandle(figno+100)
        close(figno+100)
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotPolishedOutput(phi,figno,x,y,err,rng,xlab,ylab,xlims,ylims,color,...
    file_tag,folder_name,tag,PlotErr,sig11_lmp)

global phis font_size line_width


figure(figno)

x1 = x;

y1 = y(:,1);
err1 = err(:,1);

err1(isnan(y1)) = [];
x1(isnan(y1)) = [];
y1(isnan(y1)) = [];

rng = rng(rng<=length(x1));
% if contains(tag,'sig11')
%     [~,~] = PlotCurve(x1,y1,NaN*ones(size(err1)),color,'-');
% elseif contains(tag,'sig_11_end')
%     [~,~] = PlotCurve(x1,y1,err1,color,'-');  
% %     [p,~] = PlotCurve(x1,sig11_lmp(:,1),NaN*ones(size(sig11_lmp(:,1))),color,'-.');    
% %     p.LineWidth = 0.5;
% else
    [~] = PlotErrorBar(x1(rng),y1(rng),err1(rng),color,'filled');
% end

x2 = x;

y2 = y(:,2);
err2 = err(:,2);

err2(isnan(y2)) = [];
x2(isnan(y2)) = [];
y2(isnan(y2)) = [];

% if contains(tag,'sig11')
%     [~,~] = PlotCurve(x2,y2,err2,color,'--');
% elseif contains(tag,'sig_11_end')
%     [~,~] = PlotCurve(x2,y2,err2,color,'--');
% else
%     [~,~] = PlotCurve(x2,y2,err2,color,'-');
% end
[~,~] = PlotCurve(x2,y2,err2,color,'-');



set(gca,'FontSize',font_size/1.5)

if phi==phis(end)
    xlim(xlims)
    ylim(ylims)
    xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
    ylabel(ylab,'FontSize',font_size,'Interpreter','latex')
    pbaspect([1.5 1 1])
%     pbaspect([2 1 1])
%     set(gcf,'Position',[750 400 350 350])
    set(gcf,'Position',[750 400 400 400])
    set(gcf,'color','w')
    if PlotErr==1
        xlabel('')
        xticklabels({})
    end
    saveas(gcf,[folder_name,'/',tag,'.',file_tag,'.png'])
    saveas(gcf,[folder_name,'/',tag,'.',file_tag,'.fig'])
end

% plot error
if PlotErr==1
    figure(figno+100)

    err_type = 1;
    set(gca,'FontSize',font_size/1.5)
    if err_type==0 % Bland_Altman
        x = mean(y(:,1:2),2);
        indices = randi(length(x),100,1);
        x = x(indices);
        difference = y(indices,2)-y(indices,1);
        val1 = nanmean(difference);
        val2 = val1+nanstd(difference);
        val3 = val1-nanstd(difference);
        p = plot([0 max(x)],[val1 val1],'-');
        p.Color = color;
        p = plot([0 max(x)],[val2 val2],'--');
        p.Color = color;
        p = plot([0 max(x)],[val3 val3],'--');
        p.Color = color;
        err = difference;
        s = scatter(x,err);
        s.MarkerFaceColor = color;
        s.MarkerEdgeColor = 'k';
        xlabel('$\bar \Delta$','FontSize',font_size,'Interpreter','latex')
        ylabel('$\Delta$','FontSize',font_size,'Interpreter','latex')
        xlim([0 Inf])
        save_err_tag = 'Bland-Altman';
    elseif err_type==1 % Relative error
        plot([0 max(x)],[0 0],'k--')
        rel_err = (y(:,2)-y(:,1))./max(y(:,1))*100;
        err = rel_err;
        p = plot(x,err);
        p.Color = color;
        p.LineWidth = line_width;
        ylim([-100 100])
        xlim(xlims)
        xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
        ylabel('$\%$ Err.','FontSize',font_size,'Interpreter','latex')
        save_err_tag = 'rel_err';
    end

    % if N_Kuhn==N_Kuhns(end)
    if phi==phis(end)
%         pbaspect([2 1 1])
        pbaspect([3 1 1])
%         pbaspect([4 1 1])
%         set(gcf,'Position',[750 100 350 300])
        set(gcf,'Position',[750 100 425 200])
        set(gcf,'color','w')
        saveas(gcf,[folder_name,'/',save_err_tag,'.',tag,'.',file_tag,'.png'])
        saveas(gcf,[folder_name,'/',save_err_tag,'.',file_tag,'.fig'])
    end

    MSE = mean((y(:,1)-y(:,2)).^2);
else
    if ishandle(figno+100)
        close(figno+100)
    end
end

% 
% figure(figno)
% 
% [~] = PlotErrorBar(x(rng),y(rng,1),err(rng,1),color,'filled');
% [~,~] = PlotCurve(x,y(:,2),err(:,2),color,'-');
% 
% set(gca,'FontSize',font_size/1.5)
% 
% if phi==phis(end)
%     xlim(xlims)
%     ylim(ylims)
%     xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
%     ylabel(ylab,'FontSize',font_size,'Interpreter','latex')
%     pbaspect([2 1 1])
%     set(gcf,'Position',[750 400 400 400])
%     set(gcf,'color','w')
%     if PlotErr==1
%         xlabel('')
%         xticklabels({})
%     end
%     saveas(gcf,[folder_name,'/',tag,'.',file_tag,'.png'])
%     saveas(gcf,[folder_name,'/',tag,'.',file_tag,'.fig'])
% end
% 
% % plot error
% if PlotErr==1
%     figure(figno+100)
%     plot([0 max(x)],[0 0],'k--')
%     rel_err = (y(:,1)-y(:,2))./max(y(:,1))*100;
%     p = plot(x,rel_err);
%     p.Color = color;
%     p.LineWidth = line_width;
% 
%     set(gca,'FontSize',font_size/1.5)
% 
%     % if N_Kuhn==N_Kuhns(end)
%     if phi==phis(end)
%         xlim(xlims)
%         ylim([-100 100])
%         xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
%         ylabel('$\%$ Err.','FontSize',font_size,'Interpreter','latex')
%         pbaspect([4 1 1])
%         set(gcf,'Position',[750 100 400 225])
%         set(gcf,'color','w')
%         saveas(gcf,[folder_name,'/rel_error.',tag,'.',file_tag,'.png'])
%         saveas(gcf,[folder_name,'/rel_error.',tag,'.',file_tag,'.fig'])
%     end
% else
%     if ishandle(figno+100)
%         close(figno+100)
%     end
% end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ConsolidateTheData(Override,RunTheLJCase,Controls)

global N_Kuhns phis kds Weissenbergs eq_time_factors...
    N_Kuhn phi kd Weissenberg eq_time_factor...
    sample b ka dt D damp...
    length_conversion damper_conversion...
    force_conversion energy_conversion...
    eaStar edStar...
    BeadSpringOrMeso...
    time stretch...
    sig11 sig22 sig12...
    sig11_vir sig22_vir sig12_vir...
    ka_out kd_out fa fd...
    no_attachments no_detachments...
    rxr11 rxr22 rxr12...
    n_samps n_eqtimes n_pts...
    rx ry rz...
    msd_sticker msd_tether...
    n_bonds n_to_self n_to_others cluster_coeff raw_data_folder_name...
    stress_ldg_filename stress_rlx_filename stress_eq_filename tau0...
    ModelTypes


Override0 = Override;

Directories = DefineFolders(Controls);

% Define tau0 for later
[~,D,tau0,~] = DefineTimeScale(b,length_conversion,damper_conversion,...
    1,0.204,12);

raw_data_folder_name = 'Consolidated Data/Full Sweep';
raw_fig_folder_name = 'Raw Plots/Full Sweep';
ens_avg_folder_name = 'Ensemble Plots/Full Sweep';
end_to_end_folder_name = 'End-to-end Space Comparisons/Full Sweep';
date_str = erase(Directories.beadspring_folder,'_BS/');
if RunTheLJCase==1
    raw_data_folder_name = ['Consolidated Data/LJ Study/',date_str];
    raw_fig_folder_name = ['Raw Plots/LJ Study/',date_str];
    ens_avg_folder_name = ['Ensemble Plots/LJ Study/',date_str];
    end_to_end_folder_name = ['End-to-end Space Comparisons/LJ Study/',date_str];
end
if Controls.RunChainLengthSweep==1
    raw_data_folder_name = ['Consolidated Data/Length Sweep/',date_str];
    raw_fig_folder_name = ['Raw Plots/Length Sweep/',date_str];
    ens_avg_folder_name = ['Ensemble Plots/Length Sweep/',date_str];
    end_to_end_folder_name = ['End-to-end Space Comparisons/Length Sweep/',date_str];
end
if Controls.RunLoadingRateSweep==1
    raw_data_folder_name = ['Consolidated Data/Loading Rate Sweep/',date_str];
    raw_fig_folder_name = ['Raw Plots/Loading Rate Sweep/',date_str];
    ens_avg_folder_name = ['Ensemble Plots/Loading Rate Sweep/',date_str];
    end_to_end_folder_name = ['End-to-end Space Comparisons/Loading Rate Sweep/',date_str];
end
if Controls.RunDetachmentRateSweep==1
    raw_data_folder_name = ['Consolidated Data/Detachment Rate Sweep/',date_str];
    raw_fig_folder_name = ['Raw Plots/Detachment Rate Sweep/',date_str];
    ens_avg_folder_name = ['Ensemble Plots/Detachment Rate Sweep/',date_str];
    end_to_end_folder_name = ['End-to-end Space Comparisons/Detachment Rate Sweep/',date_str];
end
if Controls.RunOscillatory==1
    raw_data_folder_name = ['Consolidated Data/Frequency Sweep/',date_str];
    raw_fig_folder_name = ['Raw Plots/Frequency Sweep/',date_str];
    ens_avg_folder_name = ['Ensemble Plots/Frequency Sweep/',date_str];
    end_to_end_folder_name = ['End-to-end Space Comparisons/Frequency Sweep/',date_str];
end
if Controls.RunLargeDeformation==1
    raw_data_folder_name = ['Consolidated Data/Large Deformation/',date_str];
    raw_fig_folder_name = ['Raw Plots/Large Deformation/',date_str];
    ens_avg_folder_name = ['Ensemble Plots/Large Deformation/',date_str];
    end_to_end_folder_name = ['End-to-end Space Comparisons/Large Deformation/',date_str];
end

% consolidate and print raw data for stress, rxr, ka, kd, fa, fd
perms = length(N_Kuhns)*length(kds)*length(phis)*length(Weissenbergs);
wb = waitbar(0,'Consolidating data...');
ct_wb = 0;
for i=1:length(N_Kuhns)
 for ii=1:length(phis)
  for iii=1:length(kds)
   for iv=1:length(Weissenbergs)
    ct_wb = ct_wb+1;
    waitbar(ct_wb/perms,wb,'Consolidating Data...')

    N_Kuhn = N_Kuhns(i);
    phi = phis(ii);
    kd = kds(iii);
    Weissenberg = Weissenbergs(iv);

    % File tage that will identify the data containing both bead-spring and
    % mesoscale model data
    file_tag = ['N_',num2str(N_Kuhn),...
        '.kd_',num2str(kd,'%.2e'),...
        '.W_',num2str(Weissenberg,'%.3f'),...
        '.phi_',num2str(phi,'%.3f')];
    if Controls.RunLargeDeformation==1
        file_tag = ['N_',num2str(N_Kuhn),...
            '.kd_',num2str(kd,'%.2e'),...
            '.W_',num2str(Weissenberg,'%.5f'),...
            '.phi_',num2str(phi,'%.3f')];
    end
    

    file_name_chk = [raw_data_folder_name,'/','timestretch.',file_tag,'.mat'];

    if ~isfile(file_name_chk) || Override==1
        % Pre-allocate the outputs
        % (e.g., sig11 = zeros(npts,N_samples*N_eq_times)
        sig11_all = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale
        sig22_all = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale
        sig12_all = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale
        kd_all = zeros(n_pts,n_samps*n_eqtimes,2);      % two layers for Bead-spring and mesoscale
        ka_all = zeros(n_pts,n_samps*n_eqtimes,2);      % two layers for Bead-spring and mesoscale
        fd_all = zeros(n_pts,n_samps*n_eqtimes,2);      % two layers for Bead-spring and mesoscale
        fa_all = zeros(n_pts,n_samps*n_eqtimes,2);      % two layers for Bead-spring and mesoscale
        n_attach_all = zeros(n_pts,n_samps*n_eqtimes,2);      % two layers for Bead-spring and mesoscale
        n_detach_all = zeros(n_pts,n_samps*n_eqtimes,2);      % two layers for Bead-spring and mesoscale
        rxr11_all = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale
        rxr22_all = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale
        rxr12_all = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale
        all_stretch = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale
        mean_rx_all = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale
        mean_ry_all = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale
        mean_rz_all = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale
        mean_norms_all = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale
        msd_st_all = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale
        msd_th_all = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale
        avg_n_bonds_all = zeros(n_pts,n_samps*n_eqtimes,2);     % two layers for Bead-spring and mesoscale
        avg_n_to_others_all = zeros(n_pts,n_samps*n_eqtimes,2); % two layers for Bead-spring and mesoscale
        avg_n_to_self_all = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale
        avg_cluster_coeff_all = zeros(n_pts,n_samps*n_eqtimes,2);% two layers for Bead-spring and mesoscale

        sig11_lmp_all = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale
        sig22_lmp_all = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale
        sig33_lmp_all = zeros(n_pts,n_samps*n_eqtimes,2);   % two layers for Bead-spring and mesoscale

        ct = 0;
        % loop overer all the different equilibrium times
        for v=1:length(eq_time_factors)
            eq_time_factor = eq_time_factors(v);
            cols = n_samps*(v-1)+1:v*n_samps;   % Columns where interpolated
            % data will be stored in arrays

            all_rx_bs = [];
            all_ry_bs = [];
            all_rz_bs = [];
            all_rx_ms = [];
            all_ry_ms = [];
            all_rz_ms = [];

            %% Callout Input Script
            for BeadSpringOrMeso=ModelTypes
                if BeadSpringOrMeso==0 % define layer of array
                    lyr = 1;  % store bead-spring data in front layer
                else
                    lyr = 2;  % store mesoscale data in back layer
                end

                [damp,D,tau0,dtFact] = DefineTimeScale(b,length_conversion,damper_conversion,...
                    BeadSpringOrMeso,phi,N_Kuhn);
                dt = tau0/dtFact;

                eaStar = -log(ka*(b*length_conversion)^2/D);
                edStar = -log(kd*(b*length_conversion)^2/D);

                sample = 1;
%                 Directories = DefineFolders(Controls);
                DefineFileNames(Controls);

                % Import data
                ImportComputedOutputs(Controls);

                % Define interpolation range
                % find indx at which loading begins
                load_start_indx = find(stretch>1,1,'first')-1;
                load_end_indx = find(stretch==max(stretch),1,'first')-1;

                % Adjust all data (zero out stresses and times, crop all other
                % data to only the loading and relaxation range)
                dat_rng = (load_start_indx:length(time))';
                time_temp = time(dat_rng)-time(load_start_indx);
                if Controls.CalculateStress==1
                    sig11_temp = sig11(dat_rng,:)-sig11(load_start_indx,:);
                    sig22_temp = sig22(dat_rng,:)-sig22(load_start_indx,:);
                    sig12_temp = sig12(dat_rng,:)-sig12(load_start_indx,:);
                    if BeadSpringOrMeso==0 % Stress measured from virial of Kuhn segments
                        sig11_vir_temp = sig11_vir(dat_rng,:);%-sig11_vir(load_start_indx,:);
                        sig22_vir_temp = sig22_vir(dat_rng,:);%-sig22_vir(load_start_indx,:);
                        sig12_vir_temp = sig12_vir(dat_rng,:);%-sig12_vir(load_start_indx,:);
                    end
                end

                if Controls.CalculateAlignment==1
                    rxr11_temp = rxr11(dat_rng,:);
                    rxr22_temp = rxr22(dat_rng,:);
                    rxr12_temp = rxr12(dat_rng,:);
                end

                if Controls.CalculateMSD==1
                    msd_st_temp = msd_sticker(1:load_start_indx,:);
                    msd_th_temp = msd_tether(1:load_start_indx,:);
                end

                if Controls.CalculateClusteringMetrics==1
                    % Average the bond clustering measures and properties
                    n_bonds_temp = squeeze(nanmean(n_bonds(1:load_start_indx,:,:),2));
                    n_to_self_temp = squeeze(nanmean(n_to_self(1:load_start_indx,:,:),2));
                    n_to_others_temp = squeeze(nanmean(n_to_others(1:load_start_indx,:,:),2));
                    cluster_coeff_temp = cluster_coeff(1:load_start_indx,:,:);
                    cluster_coeff_temp(isinf(cluster_coeff_temp)) = 0;
                    cluster_coeff_temp = squeeze(nanmean(cluster_coeff_temp,2));
                end

                if Controls.CalculateEndtoEnd==1
                    if BeadSpringOrMeso==0
                        all_rx_bs = cat(2,all_rx_bs,rx);
                        all_ry_bs = cat(2,all_ry_bs,ry);
                        all_rz_bs = cat(2,all_rz_bs,rz);
                        mean_rx = squeeze(mean(abs(all_rx_bs),2));
                        mean_ry = squeeze(mean(abs(all_ry_bs),2));
                        mean_rz = squeeze(mean(abs(all_rz_bs),2));
                        norms = (all_rx_bs.^2+all_ry_bs.^2+all_rz_bs.^2).^0.5;
                        mean_norms = squeeze(mean(norms,2));
                    else
                        all_rx_ms = cat(2,all_rx_ms,rx);
                        all_ry_ms = cat(2,all_ry_ms,ry);
                        all_rz_ms = cat(2,all_rz_ms,rz);
                        mean_rx = squeeze(mean(abs(all_rx_ms),2));
                        mean_ry = squeeze(mean(abs(all_ry_ms),2));
                        mean_rz = squeeze(mean(abs(all_rz_ms),2));
                        norms = (all_rx_ms.^2+all_ry_ms.^2+all_rz_ms.^2).^0.5;
                        mean_norms = squeeze(mean(norms,2));
                    end
                end

                % kinetics time and outputs
                interp_time_eq = (linspace(time(1),time(load_start_indx),n_pts))';
                interp_time_all = (linspace(time(1),time(end),n_pts))';

                interp_stretch = interp1(time,stretch,interp_time_all);
                all_stretch(:,cols,lyr) = repmat(interp_stretch,1,n_samps);

                % Interpolate data
                interp_time_ld = linspace(time(load_start_indx),time(end),n_pts);
                interp_time_ld = (interp_time_ld-interp_time_ld(1))';
                for vi=1:n_samps
                    ct = ct+1;
                    sample = vi;
                    Directories = DefineFolders(Controls);
                    DefineFileNames(Controls);

                    if Controls.CalculateStress==1
                        sig11_all(:,cols(vi),lyr) = interp1(time_temp,sig11_temp(:,vi),interp_time_ld);
                        sig22_all(:,cols(vi),lyr) = interp1(time_temp,sig22_temp(:,vi),interp_time_ld);
                        sig12_all(:,cols(vi),lyr) = interp1(time_temp,sig12_temp(:,vi),interp_time_ld);

                        % Remove outliers from bad periodic bounds anomaly
                        window = 25;
                        [temp_vals,indx] = rmoutliers(sig11_all(:,cols(vi),lyr),'movmedian',window);
                        sig11_all(~indx,cols(vi),lyr) = temp_vals;
                        sig11_all(indx,cols(vi),lyr) = NaN;
                        [temp_vals,indx] = rmoutliers(sig22_all(:,cols(vi),lyr),'movmedian',window);
                        sig22_all(~indx,cols(vi),lyr) = temp_vals;
                        sig22_all(indx,cols(vi),lyr) = NaN;
                        [temp_vals,indx] = rmoutliers(sig12_all(:,cols(vi),lyr),'movmedian',window);
                        sig12_all(~indx,cols(vi),lyr) = temp_vals;
                        sig12_all(indx,cols(vi),lyr) = NaN;

                        if BeadSpringOrMeso==0
                            sig11_all(:,cols(vi),3) = interp1(time_temp,sig11_vir_temp(:,vi),interp_time_ld);
                            sig22_all(:,cols(vi),3) = interp1(time_temp,sig22_vir_temp(:,vi),interp_time_ld);
                            sig12_all(:,cols(vi),3) = interp1(time_temp,sig12_vir_temp(:,vi),interp_time_ld);

                            % Stress from LAMMPS
                            eq_stress_filename = [Directories.output_folder_bs,'/',stress_eq_filename];
                            ld_stress_filename = [Directories.output_folder_bs,'/',stress_ldg_filename];
                            rl_stress_filename = [Directories.output_folder_bs,'/',stress_rlx_filename];
                        else
                            % Stress from LAMMPS
                            eq_stress_filename = [Directories.output_folder_ms,'/',stress_eq_filename];
                            ld_stress_filename = [Directories.output_folder_ms,'/',stress_ldg_filename];
                            rl_stress_filename = [Directories.output_folder_ms,'/',stress_rlx_filename];
                        end

                        if isfile(eq_stress_filename)
                            eq_stress = dlmread(eq_stress_filename,'',1,0);
                        else
                            eq_stress = [];
                        end
                        if isfile(ld_stress_filename)
                            ld_stress = dlmread(ld_stress_filename,'',1,0);
                        else
                            ld_stress = [];
                        end
                        if isfile(rl_stress_filename)
                            rl_stress = dlmread(rl_stress_filename,'',1,0);
                        else
                            rl_stress = [];
                        end

                        % Stress from LAMMPS
                        if isfile(ld_stress_filename) && isfile(rl_stress_filename)
                            sig11_lmp_temp = [ld_stress(:,3);rl_stress(:,3)]; % exclude the equilibration stress
                            sig22_lmp_temp = [ld_stress(:,4);rl_stress(:,4)];
                            sig33_lmp_temp = [ld_stress(:,5);rl_stress(:,5)];

                            time_lmp_temp = [ld_stress(:,1);...
                                rl_stress(:,1) + ld_stress(end,1)]-ld_stress(1,1);
                            sig11_lmp_all(:,cols(vi),lyr) = interp1(time_lmp_temp,sig11_lmp_temp,interp_time_ld);
                            sig22_lmp_all(:,cols(vi),lyr) = interp1(time_lmp_temp,sig22_lmp_temp,interp_time_ld);
                            sig33_lmp_all(:,cols(vi),lyr) = interp1(time_lmp_temp,sig33_lmp_temp,interp_time_ld);
                        end
                    end

                    if Controls.CalculateAlignment==1
                        rxr11_all(:,cols(vi),lyr) = interp1(time_temp,rxr11_temp(:,vi),interp_time_ld);
                        rxr22_all(:,cols(vi),lyr) = interp1(time_temp,rxr22_temp(:,vi),interp_time_ld);
                        rxr12_all(:,cols(vi),lyr) = interp1(time_temp,rxr12_temp(:,vi),interp_time_ld);
                    end

                    if Controls.CalculateMSD==1
                        msd_st_all(:,cols(vi),lyr) = interp1(time(1:load_start_indx),...
                            msd_st_temp(:,vi),interp_time_eq);
                        msd_th_all(:,cols(vi),lyr) = interp1(time(1:load_start_indx),...
                            msd_th_temp(:,vi),interp_time_eq);
                    end

                    if Controls.CalculateClusteringMetrics==1
                        avg_n_bonds_all(:,cols(vi),lyr) = interp1(time(1:load_start_indx),...
                            n_bonds_temp(:,vi),interp_time_eq);
                        avg_n_to_others_all(:,cols(vi),lyr) = interp1(time(1:load_start_indx),...
                            n_to_others_temp(:,vi),interp_time_eq);
                        avg_n_to_self_all(:,cols(vi),lyr) = interp1(time(1:load_start_indx),...
                            n_to_self_temp(:,vi),interp_time_eq);
                        avg_cluster_coeff_all(:,cols(vi),lyr) = interp1(time(1:load_start_indx),...
                            cluster_coeff_temp(:,vi),interp_time_eq);
                    end

                    if Controls.CalculateEndtoEnd==1
                        mean_rx_all(:,cols(vi),lyr) = ...
                            interp1(time_temp,mean_rx(load_start_indx:end,vi),...
                            interp_time_ld);
                        mean_ry_all(:,cols(vi),lyr) = ...
                            interp1(time_temp,mean_ry(load_start_indx:end,vi),...
                            interp_time_ld);
                        mean_rz_all(:,cols(vi),lyr) = ...
                            interp1(time_temp,mean_rz(load_start_indx:end,vi),...
                            interp_time_ld);
                        mean_norms_all(:,cols(vi),lyr) = ...
                            interp1(time_temp,mean_norms(load_start_indx:end,vi),...
                            interp_time_ld);
                    end

                    if Controls.CalculateBondKinetics==1
                        ka_all(:,cols(vi),lyr) = interp1(time,ka_out(:,vi),interp_time_all);
                        kd_all(:,cols(vi),lyr) = interp1(time,kd_out(:,vi),interp_time_all);
                        fa_all(:,cols(vi),lyr) = interp1(time,fa(:,vi),interp_time_all);
                        fd_all(:,cols(vi),lyr) = interp1(time,fd(:,vi),interp_time_all);
                        n_attach_all(:,cols(vi),lyr) = interp1(time,no_attachments(:,vi),interp_time_all);
                        n_detach_all(:,cols(vi),lyr) = interp1(time,no_detachments(:,vi),interp_time_all);
                    end
                end
            end
          
            if Controls.CalculateEndtoEnd==1
                PlotEndtoEndSpace(all_rx_bs,all_ry_bs,all_rz_bs,...
                    all_rx_ms,all_ry_ms,all_rz_ms,...
                    load_start_indx,load_end_indx,...
                    file_tag,Override,end_to_end_folder_name);
                Override = Override0;
            end

        end

        if Controls.CalculateStress==1
            if isfile(ld_stress_filename) 
                sig11_lmp_all = -sig11_lmp_all-sig11_lmp_all(1,:,:);
                sig22_lmp_all = -sig22_lmp_all-sig22_lmp_all(1,:,:);
                sig33_lmp_all = -sig33_lmp_all-sig33_lmp_all(1,:,:);

                sig11_lmp_all = movmean(sig11_lmp_all,1,1);
                sig22_lmp_all = movmean(sig22_lmp_all,1,1);
                sig33_lmp_all = movmean(sig33_lmp_all,1,1);

                sig11_lmp_all = sig11_lmp_all*force_conversion/length_conversion^2;
                sig22_lmp_all = sig22_lmp_all*force_conversion/length_conversion^2;
                sig33_lmp_all = sig33_lmp_all*force_conversion/length_conversion^2;
            end
        end

        % Plot raw stress data
        file_name_chk = [raw_fig_folder_name,'/','sig11.',file_tag,'.png'];
        if ~isfile(file_name_chk) || Override==1
            PlotTheRawData(raw_fig_folder_name,file_tag,...
                interp_time_ld,interp_time_all,interp_time_eq,all_stretch,...
                sig11_all,sig22_all,sig12_all,...
                ka_all,kd_all,fa_all,fd_all,...
                n_attach_all,n_detach_all,...
                rxr11_all,rxr22_all,rxr12_all,...
                msd_st_all,msd_th_all,...
                mean_rx_all,mean_ry_all,mean_rz_all,mean_norms_all,...
                avg_n_bonds_all,avg_n_to_self_all,...
                avg_n_to_others_all,avg_cluster_coeff_all,Controls)
            close all
        end

        % Plot a sample of the raw strain and stress data
        if Controls.RunOscillatory==1 && kd==kds(1)
            tau12 = (sig11_all-sig22_all)/2*sin(pi/2);
            eps11 = all_stretch(:,1,2)-1;
            eps22 = -eps11;
            gam12 = (eps11-eps22)/2*sin(pi/2);
            
%             mean_sig = nanmean(sig11_all(:,:,2),2);
%             se_sig = nanstd(sig11_all(:,:,2),1,2)/sqrt(size(sig11_all,2));
            mean_sig = nanmean(tau12(:,:,2),2)/1e3;
            se_sig = nanstd(tau12(:,:,2),1,2)/sqrt(size(tau12,2))/1e3;

            figure(1e3); clf; hold on
  
            subplot(2,1,1); hold on
            plot([0 max(interp_time_ld)/1e-6],[0 0],'k--')
            ld_indx = find(all_stretch(:,1,2)~=1,1,'first');
            time_tmp = interp_time_all(ld_indx:end)-interp_time_all(ld_indx);
%             p = plot(time_tmp/1e-6,eps11(ld_indx:end));
            p = plot(time_tmp/1e-6,gam12(ld_indx:end));
            p.Color = 'k';
            p.LineWidth = 1.5;
            
            set(gca,'FontSize',20/1.5)
            xticks([0 0.5 1 1.5 2])
            xticklabels({})
            ylabel('$\epsilon_{12}$','FontSize',20,'Interpreter','latex')
            pbaspect([3 1 1])

            subplot(2,1,2); hold on
            plot([0 max(interp_time_ld)/1e-6],[0 0],'k--')
%             p = plot(interp_time_ld/1e-6,mean_sig);
%             p.Color = 'k';
            [p,~] = PlotCurve(interp_time_ld/1e-6,mean_sig,se_sig,'k','-');
            p.LineWidth = 1.5;

            set(gca,'FontSize',20/1.5)
            xlabel('$t$ ($\mu$s)','FontSize',20,'Interpreter','latex')
            ylabel('$\sigma_{12}$ (kPa)','FontSize',20,'Interpreter','latex')
            pbaspect([3 1 1])

            set(gcf,'Position',[1000 100 700 400])
            set(gcf,'Color','w')
            file_name = [raw_fig_folder_name,'/','eps12 sig12.',file_tag];
            saveas(gcf,[file_name,'.png'])
            saveas(gcf,[file_name,'.fig'])
        end

        % Plot ensemble average data
        file_name_chk = [ens_avg_folder_name,'/','sig11.',file_tag,'.png'];
        if (~isfile(file_name_chk) || Override==1) && size(sig11_all,2)>1 ...
                && Controls.RunOscillatory~=1 && Controls.RunLargeDeformation~=1    
            PlotTheEnsembleData(ens_avg_folder_name,file_tag,...
                interp_time_ld,interp_time_all,interp_time_eq,...
                sig11_all,sig22_all,sig12_all,...
                sig11_lmp_all,sig22_lmp_all,sig33_lmp_all,...
                ka_all,kd_all,fa_all,fd_all,...
                n_attach_all,n_detach_all,...
                rxr11_all,rxr22_all,rxr12_all,...
                msd_st_all,msd_th_all,...
                mean_rx_all,mean_ry_all,mean_rz_all,mean_norms_all,...
                avg_n_bonds_all,avg_n_to_self_all,...
                avg_n_to_others_all,avg_cluster_coeff_all,Controls,ld_stress_filename)
            close all
        end

        file_name_chk = [raw_data_folder_name,'/','timestretch.',file_tag,'.mat'];
        if ~isfile(file_name_chk) || Override==1
            SaveOrImportConsolidatedData(interp_time_eq,...
                interp_time_ld,interp_time_all,interp_stretch,...
                sig11_all,sig22_all,sig12_all,ka_all,kd_all,...
                sig11_lmp_all,sig22_lmp_all,sig33_lmp_all,...
                n_attach_all,n_detach_all,...
                fa_all,fd_all,rxr11_all,rxr22_all,rxr12_all,...
                msd_st_all,msd_th_all,...
                raw_data_folder_name,file_tag,...
                avg_n_bonds_all,avg_n_to_self_all,...
                avg_n_to_others_all,avg_cluster_coeff_all,Controls);
        end
    end

   end
  end
 end
end
close(wb)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotEndtoEndSpace(all_rx_bs,all_ry_bs,all_rz_bs,...
    all_rx_ms,all_ry_ms,all_rz_ms,...
    load_start_indx,load_end_indx,...
    file_tag,Override,folder_name)

global N_Kuhn b font_size

if ~isfolder(folder_name)
    mkdir(folder_name)
end

lamx_bs = abs(all_rx_bs)/sqrt(N_Kuhn)/b;
lamy_bs = abs(all_ry_bs)/sqrt(N_Kuhn)/b;
lamz_bs = abs(all_rz_bs)/sqrt(N_Kuhn)/b;
lamx_ms = abs(all_rx_ms)/sqrt(N_Kuhn)/b;
lamy_ms = abs(all_ry_ms)/sqrt(N_Kuhn)/b;
lamz_ms = abs(all_rz_ms)/sqrt(N_Kuhn)/b;

for j=1:size(all_rx_bs,1)
    if j==1 || j==load_start_indx || j==load_end_indx || j==size(lamx_bs,1)
        if j==1
            tag = 't0';
        elseif j==load_start_indx
            tag = 't1';
        elseif j==load_end_indx
            tag = 't2';
        elseif j==size(lamx_bs,1)
            tag = 'tf';
        end
        file_name = [folder_name,'/',file_tag,'.',tag];
        if ~isfile([file_name,'.png']) || Override==1

            figure(100); clf; hold on

            % Plot bead-spring end-to-end distribution
            end_pts = [(lamx_bs(j,:))' (lamy_bs(j,:))' (lamz_bs(j,:))'];
            start_pts = zeros(size(end_pts));
            p = plot3([start_pts(:,1) end_pts(:,1)]',...
                [start_pts(:,2) end_pts(:,2)]',...
                [start_pts(:,3) end_pts(:,3)]','Color',[0.5 0 0.5]);

            % Plot mesoscale end-to-end distribution
            end_pts = [(lamx_ms(j,:))' (lamy_ms(j,:))' (lamz_ms(j,:))'];
            start_pts = zeros(size(end_pts));
            p = plot3([start_pts(:,1) end_pts(:,1)]',...
                [start_pts(:,2) end_pts(:,2)]',...
                [start_pts(:,3) end_pts(:,3)]','Color',[0 0.5 0.5]);
            %             p(:).Color = 'c';

            view(135,30)
            pbaspect([1 1 1])
            %                 xlim([0 N_Kuhn*b])
            %                 ylim([0 N_Kuhn*b])
            %                 zlim([0 N_Kuhn*b])
            xlim([0 sqrt(N_Kuhn)])
            ylim([0 sqrt(N_Kuhn)])
            zlim([0 sqrt(N_Kuhn)])

            xticklabels({'0','$\sqrt N$'})
            yticklabels({'0','$\sqrt N$'})
            zticklabels({'0','$\sqrt N$'})

            set(gca,'fontsize',font_size/1.5)
            hAxes = gca;
            hAxes.TickLabelInterpreter = 'latex';

            xlabel('$\lambda_x$','fontsize',font_size,'interpreter','latex')
            ylabel('$\lambda_x$','fontsize',font_size,'interpreter','latex')
            zlabel('$\lambda_x$','fontsize',font_size,'interpreter','latex')

            xticks([0 sqrt(N_Kuhn)])
            yticks([0 sqrt(N_Kuhn)])
            zticks([0 sqrt(N_Kuhn)])

            box on

            saveas(gcf,[file_name,'.png'])
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotTheEnsembleData(folder_name,file_tag,...
    interp_time_ld,interp_time_all,interp_time_eq,...
    sig11_all,sig22_all,sig12_all,...
    sig11_lmp_all,sig22_lmp_all,sig33_lmp_all,...
    ka_all,kd_all,fa_all,fd_all,...
    n_attach_all,n_detach_all,...
    rxr11_all,rxr22_all,rxr12_all,...
    msd_st_all,msd_th_all,...
    mean_rx_all,mean_ry_all,mean_rz_all,mean_norms_all,...
    avg_n_bonds_all,avg_n_to_self_all,...
    avg_n_to_others_all,avg_cluster_coeff_all,Controls,ld_stress_filename)

global kd N_Kuhn b Nt

if ~isfolder(folder_name)
    mkdir(folder_name);
end

if Controls.CalculateStress==1
    % Plot ensemble stresses
    figure(1); clf; hold on
    x = interp_time_ld/1e-9;
    y = sig11_all/1e3;
    xlab = '$t$ (ns)';
    ylab = '$\sigma_{11}$ (kPa)';
    file_name = 'sig11.';
    ylims = [0 1400];
    PlotEnsembleStress(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims,1);

    figure(2); clf; hold on
    x = interp_time_ld/1e-9;
    y = sig22_all/1e3;
    xlab = '$t$ (ns)';
    ylab = '$\sigma_{22}$ (kPa)';
    file_name = 'sig22.';
    ylims = [-Inf Inf];
    plot(x,zeros(size(x)),'k--')
    PlotEnsembleStress(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims,0);

    figure(3); clf; hold on
    x = interp_time_ld/1e-9;
    y = sig12_all/1e3;
    xlab = '$t$ (ns)';
    ylab = '$\sigma_{12}$ (kPa)';
    file_name = 'sig12.';
    ylims = [-Inf Inf];
    plot(x,zeros(size(x)),'k--')
    PlotEnsembleStress(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims,0);
end

if Controls.CalculateAlignment==1
    % Plot rest of ensemble data
    figure(4); clf; hold on
    x = interp_time_ld/1e-9;
    y = rxr11_all;
    xlab = '$t$ (ns)';
    ylab = '$g_{11}$';
    file_name = 'rxr11.';
    ylims = [0 1];
    plot(x,0.33*ones(size(x)),'k--')
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);

    figure(5); clf; hold on
    x = interp_time_ld/1e-9;
    y = rxr22_all;
    xlab = '$t$ (ns)';
    ylab = '$g_{22}$';
    file_name = 'rxr22.';
    ylims = [0 1];
    plot(x,0.33*ones(size(x)),'k--')
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);

    figure(6); clf; hold on
    x = interp_time_ld/1e-9;
    y = rxr12_all;
    xlab = '$t$ (ns)';
    ylab = '$g_{12}$';
    file_name = 'rxr12.';
    ylims = [0 1];
    plot(x,zeros(size(x)),'k--')
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);
end

if Controls.CalculateBondKinetics==1
    figure(7); clf; hold on
    x = interp_time_all/1e-9;
    y = movmean(ka_all,10,1);
    xlab = '$t$ (ns)';
    ylab = '$k_a$ (Hz)';
    file_name = 'ka.';
    ylims = [0 Inf];
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);

    figure(8); clf; hold on
    x = interp_time_all/1e-9;
    if kd==0
        kd_all(kd_all~=0) = 0;
    end
    y = movmean(kd_all,10,1);
    xlab = '$t$ (ns)';
    ylab = '$k_d$ (Hz)';
    file_name = 'kd.';
    ylims = [0 Inf];
    plot(x,kd*ones(size(x)),'k--')
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);

    figure(9); clf; hold on
    x = interp_time_all/1e-9;
    y = fa_all;
    xlab = '$t$ (ns)';
    ylab = '$f_a$, $f_d$';
    file_name = 'f.';
    ylims = [0 1];
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);
    y = fd_all;
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);
end

if Controls.CalculateEndtoEnd==1
    figure(10); clf; hold on
    % x = interp_time_kin/1e-9;
    x = interp_time_ld/1e-9;
    y = mean_rx_all/(sqrt(N_Kuhn)*b);
    xlab = '$t$ (ns)';
    ylab = '$\bar r_x/(\sqrt N b)$';
    file_name = 'rx.';
    ylims = [0 Inf];
    plot(x,1/sqrt(3)*ones(size(x)),'k--')
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);

    figure(11); clf; hold on
    % x = interp_time_kin/1e-9;
    x = interp_time_ld/1e-9;
    y = mean_ry_all/(sqrt(N_Kuhn)*b);
    xlab = '$t$ (ns)';
    ylab = '$\bar r_y/(\sqrt N b)$';
    file_name = 'ry.';
    ylims = [0 Inf];
    plot(x,1/sqrt(3)*ones(size(x)),'k--')
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);

    figure(12); clf; hold on
    % x = interp_time_kin/1e-9;
    x = interp_time_ld/1e-9;
    y = mean_rz_all/(sqrt(N_Kuhn)*b);
    xlab = '$t$ (ns)';
    ylab = '$\bar r_z/(\sqrt N b)$';
    file_name = 'rz.';
    ylims = [0 Inf];
    plot(x,1/sqrt(3)*ones(size(x)),'k--')
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);

    figure(13); clf; hold on
    % x = interp_time_kin/1e-9;
    x = interp_time_ld/1e-9;
    y = mean_norms_all/(sqrt(N_Kuhn)*b);
    xlab = '$t$ (ns)';
    ylab = '$\bar r/(\sqrt N b)$';
    file_name = 'r.';
    ylims = [0 Inf];
    plot(x,ones(size(x)),'k--')
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);
end

if Controls.CalculateMSD==1
    figure(14); clf; hold on
    % x = interp_time_kin/1e-9;
    x = interp_time_ld/1e-9;
    y = msd_st_all;
    xlab = '$t/\tau_0$)';
    ylab = '$\langle x^2 \rangle_{st}/b^2$';
    file_name = 'MSD_st.';
    ylims = [0 Inf];
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);

    figure(15); clf; hold on
    % x = interp_time_kin/1e-9;
    x = interp_time_ld/1e-9;
    y = msd_th_all;
    xlab = '$t$ (ns)';
    ylab = '$\langle x^2 \rangle_{th}/b^2$';
    file_name = 'MSD_th.';
    ylims = [0 Inf];
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);
end

if Controls.CalculateClusteringMetrics==1
    figure(16); clf; hold on
    x = interp_time_eq/1e-9;
    y = avg_n_bonds_all;
    xlab = '$t$ (ns)';
    ylab = '$n_{bonds}$';
    file_name = 'N_bonds.';
    ylims = [0 Nt];
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);

    figure(17); clf; hold on
    x = interp_time_eq/1e-9;
    y = avg_n_to_self_all;
    xlab = '$t$ (ns)';
    ylab = '$n_{self}$';
    file_name = 'N_self.';
    ylims = [0 Nt];
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);

    figure(18); clf; hold on
    x = interp_time_eq/1e-9;
    y = avg_n_to_others_all;
    xlab = '$t$ (ns)';
    ylab = '$n_{other}$';
    file_name = 'other.';
    ylims = [0 Nt];
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);

    figure(19); clf; hold on
    x = interp_time_eq/1e-9;
    y = avg_cluster_coeff_all;
    xlab = '$t$ (ns)';
    ylab = '$C_\alpha$';
    file_name = 'Cluster.';
    ylims = [0 0.5];
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);

    figure(20); clf; hold on
    x = interp_time_all/1e-9;
    y = movmean(n_attach_all,10,1);
    xlab = '$t$ (ns)';
    ylab = '$n_{attach}$';
    file_name = 'n_attach.';
    ylims = [0 Inf];
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);

    figure(22); clf; hold on
    x = interp_time_all/1e-9;
    y = movmean(n_detach_all,10,1);
    xlab = '$t$ (ns)';
    ylab = '$n_{detach}$';
    file_name = 'n_detach.';
    ylims = [0 Inf];
    PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);
end

if Controls.CalculateStress==1 && isfile(ld_stress_filename)
% Plot ensemble stresses from virial formulation in LAMMPS
figure(23); clf; hold on
x = interp_time_ld/1e-9;
y = (sig11_lmp_all-sig11_lmp_all(1,:,:))/1e3;
xlab = '$t$ (ns)';
ylab = '$\sigma_{11}$ (kPa)';
file_name = 'sig11_lmp.';
ylims = [0 Inf];
PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);

figure(24); clf; hold on
x = interp_time_ld/1e-9;
y = (sig22_lmp_all-sig22_lmp_all(1,:,:))/1e3;
xlab = '$t$ (ns)';
ylab = '$\sigma_{22}$ (kPa)';
file_name = 'sig22_lmp.';
ylims = [-Inf Inf];
plot(x,zeros(size(x)),'k--')
PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);

figure(25); clf; hold on
x = interp_time_ld/1e-9;
y = (sig33_lmp_all-sig33_lmp_all(1,:,:))/1e3;
xlab = '$t$ (ns)';
ylab = '$\sigma_{33}$ (kPa)';
file_name = 'sig33_lmp.';
ylims = [-Inf Inf];
plot(x,zeros(size(x)),'k--')
PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotEnsembleData(x,y,xlab,ylab,file_name,folder_name,file_tag,ylims)

global font_size

mean_y = squeeze(nanmean(y,2));
se_y = squeeze(nanstd(y,1,2));

x1 = x;
x2 = x;

y1 = mean_y(:,1);
se_y1 = se_y(:,1);

se_y1(isnan(y1)) = [];
x1(isnan(y1)) = [];
y1(isnan(y1)) = [];

y2 = mean_y(:,2);
se_y2 = se_y(:,2);

se_y2(isnan(y2)) = [];
x2(isnan(y2)) = [];
y2(isnan(y2)) = [];


[~,~] = PlotCurve(x1,y1,se_y1,'k','-');
[~,~] = PlotCurve(x2,y2,se_y2,'k','--');

set(gca,'FontSize',font_size/1.5)

xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
ylabel(ylab,'FontSize',font_size,'Interpreter','latex')

ylim(ylims)

pbaspect([1 1 1])
set(gcf,'color','w')

saveas(gcf,[folder_name,'/',file_name,file_tag,'.png'])
saveas(gcf,[folder_name,'/',file_name,file_tag,'.fig'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotEnsembleStress(x,y,xlab,ylab,file_name,folder_name,file_tag,...
    ylims,PlotVir)

global font_size

mean_y = squeeze(nanmean(y,2));
se_y = squeeze(nanstd(y,1,2));
if PlotVir==1
    y_vir = movmean(mean_y(:,3),10);
    y_vir = y_vir-y_vir(1);
    plot(x,y_vir,'k-.');
end

[~,~] = PlotCurve(x,mean_y(:,2),se_y(:,2),'k','--');
[~,~] = PlotCurve(x,mean_y(:,1),se_y(:,1),'k','-');

set(gca,'FontSize',font_size/1.5)

xlabel(xlab,'FontSize',font_size,'Interpreter','latex')
ylabel(ylab,'FontSize',font_size,'Interpreter','latex')

ylim(ylims)

pbaspect([1 1 1])
set(gcf,'color','w')

saveas(gcf,[folder_name,'/',file_name,file_tag,'.png'])
saveas(gcf,[folder_name,'/',file_name,file_tag,'.fig'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotTheRawData(folder_name,file_tag,...
    interp_time_ld,interp_time_all,interp_time_eq,stretch,...
    sig11_all,sig22_all,sig12_all,...
    ka_all,kd_all,fa_all,fd_all,...
    n_attach_all,n_detach_all,...
    rxr11_all,rxr22_all,rxr12_all,...
    msd_st_all,msd_th_all,...
    mean_rx_all,mean_ry_all,mean_rz_all,mean_norms_all,...
    avg_n_bonds_all,avg_n_to_self_all,...
    avg_n_to_others_all,avg_cluster_coeff_all,Controls)
    
global n_samps n_eqtimes N_Kuhn b tau0 Nt

if ~isfolder(folder_name)
    mkdir(folder_name);
end

figure(1); clf; hold on
figure(2); clf; hold on
figure(3); clf; hold on
figure(4); clf; hold on
figure(5); clf; hold on
figure(6); clf; hold on
figure(7); clf; hold on
figure(8); clf; hold on
figure(9); clf; hold on
figure(10); clf; hold on
figure(11); clf; hold on
figure(12); clf; hold on
figure(13); clf; hold on
figure(14); clf; hold on
figure(15); clf; hold on
figure(16); clf; hold on
figure(17); clf; hold on
figure(18); clf; hold on
figure(19); clf; hold on
figure(20); clf; hold on
figure(21); clf; hold on
figure(22); clf; hold on

cmap1 = cool(n_samps*n_eqtimes);
cmap2 = summer(n_samps*n_eqtimes);
for vi=1:n_samps*n_eqtimes

    if Controls.CalculateStress==1
        figno = 1;
        xlab = '$t$ (ns)';
        ylab = '$\sigma_{11}$ (kPa)';
        ylims = [0 Inf];
        file_name = 'sig11.';
        PlotRawData(figno,interp_time_ld/1e-9,sig11_all/1e3,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)

        figno = 2;
        xlab = '$t$ (ns)';
        ylab = '$\sigma_{22}$ (kPa)';
        ylims = [-Inf 0];
        file_name = 'sig22.';
        PlotRawData(figno,interp_time_ld/1e-9,sig22_all/1e3,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)

        figno = 3;
        xlab = '$t$ (ns)';
        ylab = '$\sigma_{12}$ (kPa)';
        ylims = [-Inf Inf];
        file_name = 'sig12.';
        PlotRawData(figno,interp_time_ld/1e-9,sig12_all/1e3,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)
    end

    if Controls.CalculateAlignment==1
        figno = 4;
        xlab = '$t$ (ns)';
        ylab = '$g_{11}$ (kPa)';
        ylims = [0 1];
        file_name = 'g11.';
        PlotRawData(figno,interp_time_ld/1e-9,rxr11_all,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)

        figno = 5;
        xlab = '$t$ (ns)';
        ylab = '$g_{22}$ (kPa)';
        ylims = [0 1];
        file_name = 'g22.';
        PlotRawData(figno,interp_time_ld/1e-9,rxr22_all,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)

        figno = 6;
        xlab = '$t$ (ns)';
        ylab = '$g_{12}$ (kPa)';
        ylims = [0 1];
        file_name = 'g12.';
        PlotRawData(figno,interp_time_ld/1e-9,rxr12_all,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)
    end

    if Controls.CalculateBondKinetics==1
        figno = 7;
        xlab = '$t$ (ns)';
        ylab = '$k_a$ (Hz)';
        ylims = [0 Inf];
        file_name = 'ka.';
        PlotRawData(figno,interp_time_all/1e-9,movmean(ka_all,100),...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)

        figno = 8;
        xlab = '$t$ (ns)';
        ylab = '$k_d$ (Hz)';
        ylims = [0 Inf];
        file_name = 'kd.';
        PlotRawData(figno,interp_time_all/1e-9,movmean(kd_all,100),...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)

        figno = 9;
        xlab = '$t$ (ns)';
        ylab = '$f$';
        ylims = [0 1];
        file_name = 'attached_frac.';
        PlotRawData(figno,interp_time_all/1e-9,fa_all,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)
        PlotRawData(figno,interp_time_all/1e-9,fd_all,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)
    end

    if Controls.CalculateEndtoEnd==1
        figno = 10;
        xlab = '$t$ (ns)';
        ylab = '$\lambda$';
        ylims = [0 4];
        file_name = 'stretch.';
        PlotRawData(figno,interp_time_all/1e-9,stretch,...
            cmap1*0.5,cmap2*0.5,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)

        figno = 11;
        xlab = '$t$ (ns)';
        ylab = '$\bar r_x/\sqrt Nb$';
        ylims = [0 Inf];
        file_name = 'rx.';
        PlotRawData(figno,interp_time_all/1e-9,mean_rx_all/sqrt(N_Kuhn)/b,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)

        figno = 12;
        xlab = '$t$ (ns)';
        ylab = '$\bar r_y/\sqrt Nb$';
        ylims = [0 Inf];
        file_name = 'ry.';
        PlotRawData(figno,interp_time_all/1e-9,mean_ry_all/sqrt(N_Kuhn)/b,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)

        figno = 13;
        xlab = '$t$ (ns)';
        ylab = '$\bar r_z/\sqrt Nb$';
        ylims = [0 Inf];
        file_name = 'rz.';
        PlotRawData(figno,interp_time_all/1e-9,mean_rz_all/sqrt(N_Kuhn)/b,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)

        figno = 14;
        xlab = '$t$ (ns)';
        ylab = '$\bar r/\sqrt Nb$';
        ylims = [0 Inf];
        file_name = 'r.';
        PlotRawData(figno,interp_time_all/1e-9,mean_norms_all/sqrt(N_Kuhn)/b,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)
    end

    if Controls.CalculateMSD==1
        figno = 15;
        xlab = '$t/\tau_0$';
        ylab = '$\langle \Delta x^2\rangle_{st}/b^2$)';
        ylims = [0 Inf];
        file_name = 'MSD_st.';
        PlotRawData(figno,interp_time_eq/tau0,msd_st_all,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)

        figno = 16;
        xlab = '$t/\tau_0$';
        ylab = '$\langle \Delta x^2\rangle_{th}/b^2$)';
        ylims = [0 Inf];
        file_name = 'MSD_st.';
        PlotRawData(figno,interp_time_eq/tau0,msd_th_all,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)
    end

    if Controls.CalculateClusteringMetrics==1
        figno = 17;
        xlab = '$t/\tau_0$';
        ylab = '$n_{bonds}$';
        ylims = [0 Nt];
        file_name = 'N_bonds.';
        PlotRawData(figno,interp_time_eq/tau0,avg_n_bonds_all,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)

        figno = 18;
        xlab = '$t/\tau_0$';
        ylab = '$n_{self}$';
        ylims = [0 Nt];
        file_name = 'N_self.';
        PlotRawData(figno,interp_time_eq/tau0,avg_n_to_self_all,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)

        figno = 19;
        xlab = '$t/\tau_0$';
        ylab = '$n_{other}$';
        ylims = [0 Nt];
        file_name = 'N_other.';
        PlotRawData(figno,interp_time_eq/tau0,avg_n_to_others_all,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)

        figno = 20;
        xlab = '$t/\tau_0$';
        ylab = '$C_\alpha$)';
        ylims = [0 0.5];
        file_name = 'Cluster_coeff.';
        PlotRawData(figno,interp_time_eq/tau0,avg_cluster_coeff_all,...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)

        figno = 21;
        xlab = '$t$ (ns)';
        ylab = '$n_{attach}$';
        ylims = [0 Inf];
        file_name = 'n_attach.';
        PlotRawData(figno,interp_time_all/1e-9,movmean(n_attach_all,10),...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)

        figno = 22;
        xlab = '$t$ (ns)';
        ylab = '$n_{detach}$';
        ylims = [0 Inf];
        file_name = 'n_detach.';
        PlotRawData(figno,interp_time_all/1e-9,movmean(n_detach_all,10),...
            cmap1,cmap2,vi,xlab,ylab,ylims,file_tag,file_name,folder_name)
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveOrImportConsolidatedData(time_eq,time_ld,time_all,stretch,...
    sig11_all,sig22_all,sig12_all,ka_all,kd_all,...
    sig11_lmp_all,sig22_lmp_all,sig33_lmp_all,...
    n_attach_all,n_detach_all,...
    fa_all,fd_all,rxr11_all,rxr22_all,rxr12_all,...
    msd_st_all,msd_th_all,...
    folder_name,file_tag,...
    avg_n_bonds_all,avg_n_to_self_all,...
    avg_n_to_others_all,avg_cluster_coeff_all,Controls)

if ~isfolder(folder_name)
    mkdir(folder_name);
end

timestretch.time_ld = time_ld;
timestretch.time_eq = time_eq;
timestretch.time_all = time_all;
timestretch.stretch = stretch;

file_name = [folder_name,'/','timestretch.',file_tag,'.mat'];
save(file_name,'-struct','timestretch');

if Controls.CalculateStress==1
    stress.sig11 = sig11_all;
    stress.sig22 = sig22_all;
    stress.sig12 = sig12_all;
    stress.sig11_lmp = sig11_lmp_all;
    stress.sig22_lmp = sig22_lmp_all;
    stress.sig33_lmp = sig33_lmp_all;

    file_name = [folder_name,'/','stress.',file_tag,'.mat'];
    save(file_name,'-struct','stress');
end

if Controls.CalculateBondKinetics==1
    kinetics.ka = ka_all;
    kinetics.kd = kd_all;
    kinetics.fa = fa_all;
    kinetics.fd = fd_all;
    kinetics.n_attach = n_attach_all;
    kinetics.n_detach = n_detach_all;
end
if Controls.CalculateClusteringMetrics==1
    kinetics.n_bonds = avg_n_bonds_all;
    kinetics.n_self = avg_n_to_self_all;
    kinetics.n_others = avg_n_to_others_all;
    kinetics.cluster_coeff = avg_cluster_coeff_all;
end
if Controls.CalculateBondKinetics==1 || Controls.CalculateClusteringMetrics==1
    file_name = [folder_name,'/','kinetics.',file_tag,'.mat'];
    save(file_name,'-struct','kinetics');
end

if Controls.CalculateAlignment==1
    covariance.rxr11_all = rxr11_all;
    covariance.rxr22_all = rxr22_all;
    covariance.rxr12_all = rxr12_all;

    file_name = [folder_name,'/','covariance.',file_tag,'.mat'];
    save(file_name,'-struct','covariance');
end

if Controls.CalculateMSD==1
    MSD.msd_st = msd_st_all;
    MSD.msd_th = msd_th_all;

    file_name = [folder_name,'/','MSD.',file_tag,'.mat'];
    save(file_name,'-struct','MSD');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotRawData(figno,x,y,cmap1,cmap2,vi,xlab,ylab,ylims,...
    file_tag,file_name,folder_name)

global font_size n_samps n_eqtimes

figure(figno)
% Plot bead-spring
p = plot(x,y(:,vi,1));
p.LineStyle = '-';
p.Color = cmap1(vi,:);

% Plot mesoscale
p = plot(x,y(:,vi,2));
p.LineStyle = '-';
p.Color = cmap2(vi,:);
% s = scatter(x,y(:,vi,2)/1e3);
% s.Marker = 'o';
% s.SizeData = 10;
% s.MarkerEdgeColor = cmap2(vi,:);
% s.MarkerFaceColor = 'none';

set(gca,'FontSize',font_size/1.5)

if vi==n_samps*n_eqtimes
    xlabel(xlab,'FontSize',font_size,'interpreter','latex')
    ylabel(ylab,'FontSize',font_size,'interpreter','latex')
    set(gcf,'color','w')
    pbaspect([1 1 1])

    ylim(ylims)

    saveas(gcf,[folder_name,'/',file_name,file_tag,'.png'])
    saveas(gcf,[folder_name,'/',file_name,file_tag,'.fig'])
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImportComputedOutputs(Controls)

global BeadSpringOrMeso time_stretch_data_filename...
    time stretch volume...
    stress_data_filename...
    clustering_data_filename...
    sig11 sig22 sig33 sig12 sig23 sig31...
    sig11_vir sig22_vir sig33_vir sig12_vir sig23_vir sig31_vir...
    bond_kinetics_filename...
    ka_out ka1_out kd_out Na Nd fa fd...
    no_attachments no_detachments...
    attached_bond_lifetimes detached_bond_lifetimes renormalized_bond_lifetimes...
    open_repeat_lifetimes open_exchange_lifetimes...
    total_exchange_events...
    total_repeat_events...
    chain_concentration...
    alignment_data_filename...
    rxr11 rxr22 rxr33 rxr12 rxr23 rxr31...
    rxr11_err rxr22_err rxr33_err rxr12_err rxr23_err rxr31_err...
    endtoend_data_filename  ...
    rx ry rz...
    msd_data_filename...
    msd_sticker msd_st_err...
    msd_st_r0 msd_st_r0_err...
    msd_st_r00 msd_st_r00_err...
    msd_st_t0 msd_st_t0_err...
    msd_st_t00 msd_st_t00_err...
    msd_tether msd_err_tether...
    n_bonds n_to_self n_to_others cluster_coeff

TimeStretchData = load(time_stretch_data_filename,'-mat');
time = TimeStretchData.time;
stretch = TimeStretchData.stretch;
volume = TimeStretchData.V;

if Controls.CalculateStress==1
    Stress = load(stress_data_filename,'-mat');
    sig11 = Stress.sig11;
    sig22 = Stress.sig22;
    sig33 = Stress.sig33;
    sig12 = Stress.sig12;
    sig23 = Stress.sig23;
    sig31 = Stress.sig31;
    if BeadSpringOrMeso==0
        sig11_vir = Stress.sig11_vir;
        sig22_vir = Stress.sig22_vir;
        sig33_vir = Stress.sig33_vir;
        sig12_vir = Stress.sig12_vir;
        sig23_vir = Stress.sig23_vir;
        sig31_vir = Stress.sig31_vir;
    end
end

if Controls.CalculateBondKinetics==1
    BondKinetics = load(bond_kinetics_filename,'-mat');
    ka_out = BondKinetics.ka;
    ka1_out = BondKinetics.ka1;
    kd_out = BondKinetics.kd;
    Na = BondKinetics.Na;
    Nd = BondKinetics.Nd;
    no_attachments = BondKinetics.no_attachments;
    no_detachments = BondKinetics.no_detachments;
    fa = Na./(Na+Nd/2);
    fd = 1-fa;
    total_exchange_events = BondKinetics.ExchangeEvents;
    total_repeat_events = BondKinetics.RepeatEvents;
    chain_concentration = BondKinetics.Concentration;
    attached_bond_lifetimes = BondKinetics.AttachedLifetimes;
    detached_bond_lifetimes = BondKinetics.DetachedLifetimes;
    renormalized_bond_lifetimes = BondKinetics.RenormalizedLifetimes;
    open_repeat_lifetimes = BondKinetics.OpenRepeatLifetimes;
    open_exchange_lifetimes = BondKinetics.OpenExchangeLifetimes;
end

if Controls.CalculateAlignment==1
    RxR = load(alignment_data_filename,'-mat');
    rxr11 = RxR.rxr11;
    rxr22 = RxR.rxr22;
    rxr33 = RxR.rxr33;
    rxr12 = RxR.rxr12;
    rxr23 = RxR.rxr23;
    rxr31 = RxR.rxr31;
    rxr11_err = RxR.rxr11_err;
    rxr22_err = RxR.rxr22_err;
    rxr33_err = RxR.rxr33_err;
    rxr12_err = RxR.rxr12_err;
    rxr23_err = RxR.rxr23_err;
    rxr31_err = RxR.rxr31_err;
end

if Controls.CalculateEndtoEnd==1
    EndtoEnd = load(endtoend_data_filename,'-mat');
    rx = EndtoEnd.rx;
    ry = EndtoEnd.ry;
    rz = EndtoEnd.rz;
end

if Controls.CalculateMSD==1
    MSD = load(msd_data_filename,'-mat');
    msd_sticker = MSD.msd;
    msd_st_err = MSD.msd_err;
    msd_st_r0 = MSD.msd_st_r0;
    msd_st_r0_err = MSD.msd_st_r0_err;
    msd_st_r00 = MSD.msd_st_r00;
    msd_st_r00_err = MSD.msd_st_r00_err;
    msd_st_t0 = MSD.msd_st_t0;
    msd_st_t0_err = MSD.msd_st_t0_err;
    msd_st_t00 = MSD.msd_st_t00;
    msd_st_t00_err = MSD.msd_st_t00_err;
    msd_tether = MSD.msd_tether;
    msd_err_tether = MSD.msd_err_tether;
end

if Controls.CalculateClusteringMetrics==1
    ClusteringData = load(clustering_data_filename,'-mat');
    n_bonds = ClusteringData.n_bonds;
    n_to_self = ClusteringData.n_to_self;
    n_to_others = ClusteringData.n_to_other;
    cluster_coeff = ClusteringData.c_alpha;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotHistograms(rx,ry,rz,bondtypes,time,stretch,MakeMovies,Override)

global FileTagCompiled binMax binMin binWidth font_size ConstantsFileName...
    MovieFolderName

CheckBackBoneVsStickers = 0;
OneDOnly=1;
norms = (rx.^2+ry.^2+rz.^2).^(1/2);

% FileName = ['Matlab Topological Data/Constants.S1',FileTagCompiled,'.txt'];
% 
% Table = readtable(FileName);
Table = readtable(ConstantsFileName);
constants = table2array(Table);
N_Kuhn = constants(4); b = constants(5);

%Histo inputs
L = N_Kuhn*b;
binMax = L;   %nm                   %Histogram inputs
binMin = -binMax;
Nbins = 101;
binWidth = (binMax-binMin)/Nbins;

%Define movie names
MovieName1 = [MovieFolderName,'/End-to-end/Norms',FileTagCompiled];
if OneDOnly~=1
    MovieName2 = [MovieFolderName,'/End-to-end/Norms_xyz',FileTagCompiled];

    MovieName3 = [MovieFolderName,'/End-to-end/xy',FileTagCompiled];
    MovieName4 = [MovieFolderName,'/End-to-end/yz',FileTagCompiled];
    MovieName5 = [MovieFolderName,'/End-to-end/zx',FileTagCompiled];

    if CheckBackBoneVsStickers==1
        MovieName6 = [MovieFolderName,'/End-to-end/Norms_branch',FileTagCompiled];
        MovieName7 = [MovieFolderName,'/End-to-end/xy_branch',FileTagCompiled];
        MovieName8 = [MovieFolderName,'/End-to-end/yz_branch',FileTagCompiled];

        MovieName9 = [MovieFolderName,'/End-to-end/Norms_backbone',FileTagCompiled];
        MovieName10 = [MovieFolderName,'/End-to-end/xy_backbone',FileTagCompiled];
        MovieName11 = [MovieFolderName,'/End-to-end/yz_backbone',FileTagCompiled];
    end
end

if MakeMovies==1 && (~isfile([MovieName1,'.avi']) || Override==1)
    if ~isfolder([MovieFolderName,'/End-to-end'])
        mkdir([MovieFolderName,'/End-to-end'])
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

        figure(201); clf; hold on
        XLim = sqrt(N_Kuhn);
        YLim = 1.5;
        color_t = 'k';
        SymDist = 0;
        PlotIdeal = 1;
        %Exclude dynamic bonds
        norms_temp12 = norms_temp;
        norms_temp12(bondtypes_temp==3) = [];
        Plot1DHistogram(norms_temp12,color_t,0.6,XLim,YLim,N_Kuhn,b,...
            SymDist,PlotIdeal)

        set(gca,'FontSize',font_size/1.5)
        xlabel('$\lambda$','FontSize',font_size,'Interpreter','latex')
        ylabel('$p$','FontSize',font_size,'Interpreter','latex')

        title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
            ', $t$ = ',num2str(time(i),'%.1f'),' s'],'FontSize',font_size/2,...
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
                SymDist,PlotIdeal)
            color_t = 'g';
            Plot1DHistogram(ry_temp,color_t,0.6,XLim,YLim,N_Kuhn,b,...
                SymDist,PlotIdeal)
            color_t = 'm';
            Plot1DHistogram(rz_temp,color_t,0.6,XLim,YLim,N_Kuhn,b,...
                SymDist,PlotIdeal)

            xlabel('$\lambda$','FontSize',font_size,'Interpreter','latex')
            ylabel('$p$','FontSize',font_size,'Interpreter','latex')

            title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
                ', $t$ = ',num2str(time(i),'%.1f'),' s'],'FontSize',font_size/2,...
                'Interpreter','latex')

            l = legend('$|r_x|$','$|r_y|$','$|r_z|$');
            l.Interpreter = 'latex'; l.font_size = font_size/2;
            l.Location = 'NorthEast';

            F2 = getframe(gcf);
            writeVideo(v2,F2);
            mov2(i) = F2;

            figure(203);

            Plot2DHistogram([rx_temp ry_temp],N_Kuhn,b)

            xlabel('$\lambda_x$','FontSize',font_size,'Interpreter','latex')
            ylabel('$\lambda_y$','FontSize',font_size,'Interpreter','latex')

            title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
                ', $t$ = ',num2str(time(i),'%.1f'),' s'],'FontSize',font_size/2,...
                'Interpreter','latex')

            daspect([1 1 1])
            set(gcf,'Color','w')
            box on

            F3 = getframe(gcf);
            writeVideo(v3,F3);
            mov3(i) = F3;

            figure(204);

            Plot2DHistogram([ry_temp rz_temp],N_Kuhn,b)

            xlabel('$\lambda_y$','FontSize',font_size,'Interpreter','latex')
            ylabel('$\lambda_z$','FontSize',font_size,'Interpreter','latex')

            title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
                ', $t$ = ',num2str(time(i),'%.1f'),' s'],'FontSize',font_size/2,...
                'Interpreter','latex')

            daspect([1 1 1])
            set(gcf,'Color','w')
            box on

            F4 = getframe(gcf);
            writeVideo(v4,F4);
            mov4(i) = F4;

            figure(205);

            Plot2DHistogram([rx_temp rz_temp],N_Kuhn,b)

            xlabel('$\lambda_x$','FontSize',font_size,'Interpreter','latex')
            ylabel('$\lambda_z$','FontSize',font_size,'Interpreter','latex')

            title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
                ', $t$ = ',num2str(time(i),'%.1f'),' s'],'FontSize',font_size/2,...
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
                Plot1DHistogram(norms_temp(bondtypes_temp==2),color_t,0.6,XLim,YLim,N_Kuhn,b,...
                    SymDist,PlotIdeal)

                set(gca,'FontSize',font_size/1.5)
                xlabel('$\lambda$','FontSize',font_size,'Interpreter','latex')
                ylabel('$p$','FontSize',font_size,'Interpreter','latex')

                title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
                    ', $t$ = ',num2str(time(i),'%.1f'),' s'],'FontSize',font_size/2,...
                    'Interpreter','latex')

                F6 = getframe(gcf);
                writeVideo(v6,F6);
                mov6(i) = F6;

                figure(207);

                Plot2DHistogram([rx_temp(bondtypes_temp==2) ry_temp(bondtypes_temp==2)],N_Kuhn,b)

                xlabel('$\lambda_x$','FontSize',font_size,'Interpreter','latex')
                ylabel('$\lambda_y$','FontSize',font_size,'Interpreter','latex')

                title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
                    ', $t$ = ',num2str(time(i),'%.1f'),' s'],'FontSize',font_size/2,...
                    'Interpreter','latex')

                daspect([1 1 1])
                set(gcf,'Color','w')
                box on

                F7 = getframe(gcf);
                writeVideo(v7,F7);
                mov7(i) = F7;

                figure(208);

                Plot2DHistogram([ry_temp(bondtypes_temp==2) rz_temp(bondtypes_temp==2)],N_Kuhn,b)

                xlabel('$\lambda_y$','FontSize',font_size,'Interpreter','latex')
                ylabel('$\lambda_z$','FontSize',font_size,'Interpreter','latex')

                title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
                    ', $t$ = ',num2str(time(i),'%.1f'),' s'],'FontSize',font_size/2,...
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
                Plot1DHistogram(norms_temp(bondtypes_temp==1),color_t,0.6,XLim,YLim,N_Kuhn,b,...
                    SymDist,PlotIdeal)

                set(gca,'FontSize',font_size/1.5)
                xlabel('$\lambda$','FontSize',font_size,'Interpreter','latex')
                ylabel('$p$','FontSize',font_size,'Interpreter','latex')

                title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
                    ', $t$ = ',num2str(time(i),'%.1f'),' s'],'FontSize',font_size/2,...
                    'Interpreter','latex')

                F9 = getframe(gcf);
                writeVideo(v9,F9);
                mov9(i) = F9;

                figure(210)
                Plot2DHistogram([rx_temp(bondtypes_temp==1) ry_temp(bondtypes_temp==1)],N_Kuhn,b)

                xlabel('$\lambda_x$','FontSize',font_size,'Interpreter','latex')
                ylabel('$\lambda_y$','FontSize',font_size,'Interpreter','latex')

                title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
                    ', $t$ = ',num2str(time(i),'%.1f'),' s'],'FontSize',font_size/2,...
                    'Interpreter','latex')

                daspect([1 1 1])
                set(gcf,'Color','w')
                box on

                F10 = getframe(gcf);
                writeVideo(v10,F10);
                mov10(i) = F10;

                figure(211);

                Plot2DHistogram([ry_temp(bondtypes_temp==1) rz_temp(bondtypes_temp==1)],N_Kuhn,b)

                xlabel('$\lambda_y$','FontSize',font_size,'Interpreter','latex')
                ylabel('$\lambda_z$','FontSize',font_size,'Interpreter','latex')

                title(['$\lambda$ =',num2str(stretch(i),'%.2f'),...
                    ', $t$ = ',num2str(time(i),'%.1f'),' s'],'FontSize',font_size/2,...
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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plot2DHistogram(HistEndtoEnd,N_Kuhn,b)
% Plot2DHistogram(HistEndtoEnd,v_hat,g_hat)

global binMax binMin binWidth font_size

Edges = {binMin/(sqrt(N_Kuhn)*b):binWidth/(sqrt(N_Kuhn)*b):binMax/(sqrt(N_Kuhn)*b)...
    binMin/(sqrt(N_Kuhn)*b):binWidth/(sqrt(N_Kuhn)*b):binMax/(sqrt(N_Kuhn)*b)};

Lambda = HistEndtoEnd./(ones(size(HistEndtoEnd))*sqrt(N_Kuhn)*b);
hist3(Lambda,'CdataMode','auto',...
    'FaceColor','interp',...
    'Edges',Edges,'line_width',0.25);
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
% TNTFileName = ['Stress Relaxation/TNT Data/Full Weissenberg.',...
%     num2str(Weissenberg,'%.3f'),'.EF.',...
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
set(gca,'FontSize',font_size/1.5)

% line_width = 0.75;
% Color = 'k';
% [~] = circle(0,0,binMax,Color,line_width);
% Color = 'r';
% [~] = circle(0,0,0.95*binMax,Color,line_width);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plot1DHistogram(Norms,Color,Alpha,XLim,YLim,N_Kuhn,b,...
    SymDist,PlotIdeal)

global binMax binMin binWidth font_size

Edges = (binMin:binWidth:binMax)/(sqrt(N_Kuhn)*b);

yshift = 1;

hold on 
% lambda = Norms/NormFact;
% histogram(lambda,'Normalization','Probability','BinEdges',Edges/NormFact);
% plot([0.9 0.9]/NormFact,[0 0.5],'k--')
Lambda = Norms/(sqrt(N_Kuhn)*b);
h = histogram(Lambda,'Normalization','pdf','BinEdges',Edges);
h.FaceColor = Color;
h.FaceAlpha = Alpha;
h.EdgeColor = 'k';
h.EdgeAlpha = 1;

% xlim([0 sqrt(N_Kuhn)])
if SymDist==1
    xlim([-XLim XLim])
else
xlim([0 XLim])
end
ylim([0 YLim])

set(gca,'FontSize',font_size/1.5);
set(gcf,'Color','w')
pbaspect([3 1 1])

% b_temp = 0.005;
% N_temp = L/b_temp;
if PlotIdeal==1
    r = (linspace(0,10*N_Kuhn*b,10000))';%*sqrt(N_Kuhn)*b;
    lam = r/(sqrt(N_Kuhn)*b);
    sigma = sqrt(N_Kuhn/6)*b;
    P = 4*pi*(r.^2)*((sigma*sqrt(2*pi))^(-3)).*exp(-1/2*(r/sigma).^2);
%     P = sqrt(2/pi)*(3/(N_Kuhn*b^2))^3/2 * x.^2 .* exp(-(3*x.^2)/(2*N_Kuhn*b^2));
%     Area = trapz(lam,P);
%     P = P/Area;
%     P_lam = 4*pi*N_Kuhn*(b^2)*((3/(N_Kuhn*(b^2)*pi))^(3/2))*(lam.^2).*exp(-3*lam.^2);
    Area = trapz(lam,P);
    P = P/Area;
    p = plot(lam,P,'k:');
    p.LineWidth = 1.5;

    P = P*yshift;
    p = plot(lam,P,'k-');
    p.LineWidth = 1.5;
    
    yDiscrete = (h.Values)';
    Bins = (h.BinEdges)';
    xDiscrete = (Bins(1:end-1)+Bins(2:end))/2*sqrt(N_Kuhn)*b;
    yIdeal = interp1(r,P,xDiscrete);
    %     yIdeal = sqrt(2/pi)*(3/(N_Kuhn*b^2))^3/2 * xDiscrete.^2 .* exp(-(3*xDiscrete.^2)/(2*N_Kuhn*b^2));
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

    %
    %     x = Edges*sqrt(N_Kuhn)*b;
    % %     P = sqrt(2/pi)*(3/(N_Kuhn*b^2))^3/2 * x.^2 .* exp(-(3*x.^2)/(2*N_Kuhn*b^2));
    %     P = (3/(2*pi*(N_Kuhn*b^2)))^(3/2)*exp(-3*x.^2/(2*N_Kuhn*b^2))*4*pi.*x.^2;
    %     Area = trapz(x,P);
    %     P = P/Area;
    % %     x2 = linspace(0,N_Kuhn*b,1000);
    % %     P2 = interp1(x,P,x2);
    %     p = plot(Edges,P,'b');
    %     p.LineWidth = 1.5;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = PlotErrorBar(X,Y,Err,Color,fill_type)

e = errorbar(X,Y,Err);
e.Marker = 'o';
e.Color = Color;
e.MarkerEdgeColor = Color;
if contains(fill_type,'filled')
    e.MarkerFaceColor = Color;
else
    e.MarkerFaceColor = 'none';
end
% e.LineStyle = 'none';
% e.MarkerSize = 2;
% e.MarkerFaceColor = 'none';
e.LineStyle = 'none';
e.MarkerSize = 4;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,f] = PlotCurve(X,Y,Err,Color,Style)

global line_width

ErrX = [X;flipud(X)];
ErrY = [Y+Err;flipud(Y-Err)];
f = fill(ErrX,ErrY,Color);
f.FaceAlpha = 0.25;
f.LineStyle = 'none';

p = plot(X,Y);
p.Color = Color;
p.LineWidth = line_width;
p.LineStyle = Style;

end
