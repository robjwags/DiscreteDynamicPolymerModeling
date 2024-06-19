function PlotTimingCases(Directories,model_types)

file_folder = Directories.beadspring_folder;
splitString = strsplit(file_folder, '/');
file_folder = splitString{1}; % Extract the first part

wall_cock_file = [file_folder,'/wall_clock.txt'];
wall_clock = readtable(wall_cock_file);
wc_dat = table2array(wall_clock);

cpu_time_file = [file_folder,'/cpu_time.txt'];
cpu_time = readtable(cpu_time_file);
cpu_dat =table2array(cpu_time);

file_sizes_file = [file_folder,'/file_sizes.txt'];
file_sizes = readtable(file_sizes_file);
fs_dat =table2array(file_sizes);

Nps = wc_dat(:,1);

x = Nps;
xlims = [40 2e3];

figno = 1;
line_style = '--';
fill_style = 'filled';
figure(figno); clf; hold on
figure(figno+100); clf; hold on
figure(figno+200); clf; hold on
y_bs = wc_dat(:,2);
y_ms = wc_dat(:,3);
ylab = 'Wall-clock time (s)';
ylims = [1e2 1e4];
[wc_func,R2_wc] = CompareSystems(x,y_bs,y_ms,ylab,ylims,xlims,figno,fill_style,line_style);

figure(figno+100)
xpl = (linspace(x(1),x(end-1),100))';
plot(xpl,100*(1-0.81*xpl.^(-0.44)),'k--')

figno = 2;
figure(figno); clf; hold on
figure(figno+100); clf; hold on
figure(figno+200); clf; hold on
y_bs = cpu_dat(:,2);
y_ms = cpu_dat(:,3);
ylab = 'CPU time (s)';
ylims = [5e3 1e6];
[cpu_func,R2_cpu] = CompareSystems(x,y_bs,y_ms,ylab,ylims,xlims,figno,fill_style,line_style);

figno = 3;
figure(figno); clf; hold on
figure(figno+100); clf; hold on
figure(figno+200); clf; hold on
y_bs = fs_dat(:,3)./fs_dat(:,2)*1e3;
y_ms = fs_dat(:,end-1)./fs_dat(:,end-2)*1e3;
ylab = 'Storage (KB/step)';
xlims = [50 2e3];
ylims = [0 Inf];
[fs_atm_func,R2_atm_cpu] = CompareSystems(x,y_bs,y_ms,ylab,ylims,xlims,figno,fill_style,line_style);

figure(figno+100)
plot(x(1:end-1),ones(size(x(1:end-1)))*nanmean((y_bs-y_ms)./y_bs*100),['k',line_style])

y_bs = fs_dat(:,4)./fs_dat(:,2)*1e3;
y_ms = fs_dat(:,end)./fs_dat(:,end-2)*1e3;
fill_style = 'empty';
line_style = '-.';
ylab = 'Storage (bytes/step)';
xlims = [50 2e3];
ylims = [0 Inf];
[fs_bnd_func,R2_bnd_cpu] = CompareSystems(x,y_bs,y_ms,ylab,ylims,xlims,figno,fill_style,line_style);

figure(figno+100)
plot(x(1:end-1),ones(size(x(1:end-1)))*nanmean((y_bs-y_ms)./y_bs*100),['k',line_style])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [func,R2] = CompareSystems(x,y_bs,y_ms,ylab,ylims,xlims,figno,...
    fill_style,line_style)

% for mt = model_types
%     symb = symbols{mt+1};

ran_indx = find(isnan(y_bs),1,'first')-1;
xplot_ran = (linspace(x(1),x(ran_indx),100))';
xplot_extrap = (linspace(x(ran_indx),x(end),100))';

figure(figno)
s = scatter(x(1:ran_indx),y_ms(1:ran_indx));
s.Marker = '^';
if contains(fill_style,'filled')
    s.MarkerFaceColor = 'k';
end
s.MarkerEdgeColor = 'k';

[func_bs,R2_bs,~] = FitLine(x,y_bs);
plot(xplot_ran,func_bs(1)*xplot_ran + func_bs(2),['k',line_style])

s = scatter(x(1:ran_indx),y_bs(1:ran_indx));
s.Marker = 'o';
if contains(fill_style,'filled')
    s.MarkerFaceColor = 'k';
end
s.MarkerEdgeColor = 'k';

[func_ms,R2_ms,~] = FitLine(x,y_ms);
plot(xplot_ran,func_ms(1)*xplot_ran + func_ms(2),['k',line_style])

set(gca,'xscale','log')
set(gca,'yscale','log')

set(gca,'FontSize',20/1.5)
xlabel('Number of polymers, $n_p$','FontSize',20,'Interpreter','latex')
% xlabel('Number of polymers, $\mathcal{N}$','FontSize',20,'Interpreter','latex')
ylabel(ylab,'FontSize',20,'Interpreter','latex')
pbaspect([1 1 1])
set(gcf,'color','w')
xlim(xlims)
ylim(ylims)
set(gcf,'Position',[1000 100 400 400])


figure(figno+100)
pcnt_red = (y_bs-y_ms)./y_bs*100;
s = scatter(x,pcnt_red);
s.Marker = 'o';
if contains(fill_style,'filled')
    s.MarkerFaceColor = 'k';
end
s.MarkerEdgeColor = 'k';
% s.SizeData = 80;
set(gca,'xscale','log')

[mdl,R2] = FitDecayingExponential(x,pcnt_red);
% [mdl,R2] = FitPwr2(x,pcnt_red);
% plot(xplot_ran,mdl(xplot_ran)*100,'k--')

set(gca,'FontSize',20/1.5)
xlabel('Number of polymers, $n_p$','FontSize',20,'Interpreter','latex')
% xlabel('Number of polymers, $\mathcal{N}$','FontSize',20,'Interpreter','latex')
ylabel('\% Reduction','FontSize',20,'Interpreter','latex')
pbaspect([2 1 1])
set(gcf,'color','w')
ylim([80 100])
xlim(xlims)
set(gcf,'Position',[1000 100 400 400])


func = [func_bs;func_ms];
R2 = [R2_bs;R2_ms];


figure(figno+200)
ratio = y_bs./y_ms;
s = scatter(x,ratio);
s.Marker = 'o';
if contains(fill_style,'filled')
    s.MarkerFaceColor = 'k';
end
s.MarkerEdgeColor = 'k';
% s.SizeData = 80;
set(gca,'xscale','log')
set(gca,'yscale','log')

[mdl,R2] = FitPower(x,ratio);
plot(xplot_ran,mdl(xplot_ran),'k--')

set(gca,'FontSize',20/1.5)
xlabel('Number of polymers, $n_p$','FontSize',20,'Interpreter','latex')
% xlabel('Number of polymers, $\mathcal{N}$','FontSize',20,'Interpreter','latex')
ylabel('ratio','FontSize',20,'Interpreter','latex')
pbaspect([2 1 1])
set(gcf,'color','w')
ylim([0 Inf])
xlim(xlims)
set(gcf,'Position',[1000 100 400 400])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,R2,p_val] = FitLine(x,y)

x(isnan(y)) = [];
y(isnan(y)) = [];

% Fit a line (polynomial of degree 1) to the data
p = polyfit(x, y, 1);

% Extract the coefficients (slope and intercept)
slope = p(1);
intercept = p(2);

% Compute the y-values for the fitted line
y_fit = polyval(p, x);

% Compute the correlation coefficient (R-squared value)
R = corrcoef(y, y_fit);
R2 = R(1, 2)^2;

% Compute the degrees of freedom
n = length(y);
p_val = 2 * (1 - tcdf(abs(R(1,2) * sqrt(n - 2)) / sqrt(1 - R(1, 2)^2), n - 2));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mdl,R2] = FitPower(x,y)

x(isnan(y)) = [];
y(isnan(y)) = [];

% Define the power-law model function
powerLawModel = @(a, b, c, x) a * x.^b + c;

% Set lower and upper limits for the coefficients
startPoints = [1, 1, 80];
lowerLimits = [0, 0, 0]; % Lower limits for [a, b]
upperLimits = [Inf, 2, 100]; % Upper limits for [a, b]

% Fit the power-law model to the data with constraints
mdl = fit(x,y,powerLawModel,'StartPoint',startPoints,'Lower',lowerLimits,'Upper',upperLimits);

% Extract the coefficients
a = mdl.a;
b = mdl.b;
c = mdl.c;

% Compute the R-squared value
y_fit = powerLawModel(a, b, c, x);
R = corrcoef(y, y_fit);
R2 = R(1, 2)^2;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mdl, R2] = FitLogarithmic(x, y)
    % Remove NaN values
    x(isnan(y)) = [];
    y(isnan(y)) = [];

    % Define the logarithmic model function
    logarithmicModel = @(a, b, x) a * log(x) + b;

    % Set start points, lower, and upper limits for the coefficients
    startPoints = [1, 1];
    lowerLimits = [-Inf, -Inf];
    upperLimits = [Inf, Inf];

    % Fit the logarithmic model to the data with constraints
    mdl = fit(x, y, logarithmicModel, 'StartPoint', startPoints, 'Lower', lowerLimits, 'Upper', upperLimits);

    % Extract the coefficients
    a = mdl.a;
    b = mdl.b;

    % Compute the R-squared value
    y_fit = logarithmicModel(a, b, x);
    R = corrcoef(y, y_fit);
    R2 = R(1, 2)^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mdl, R2] = FitPwr2(x, y)
    % Remove NaN values
    x(isnan(y)) = [];
    y(isnan(y)) = [];

    % Define the decaying exponential model function
    decayingExponentialModel = @(a, b, x) (1-a*x.^b);

    % Set start points, lower, and upper limits for the coefficients
    startPoints = [0.81,-0.44]; % starting values for [a, x0]
    lowerLimits = [0.81, -0.43]; % lower limits for [a, x0]
    upperLimits = [0.81, -0.43]; % upper limits for [a, x0]

    % Fit the decaying exponential model to the data with constraints
    mdl = fit(x, y, decayingExponentialModel, 'StartPoint', startPoints, 'Lower', lowerLimits, 'Upper', upperLimits);

    % Extract the coefficients
    a = mdl.a;
    b = mdl.b;

    % Compute the R-squared value
    y_fit = decayingExponentialModel(a, b, x);
    R = corrcoef(y, y_fit);
    R2 = R(1, 2)^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mdl, R2] = FitDecayingExponential(x, y)
    % Remove NaN values
    x(isnan(y)) = [];
    y(isnan(y)) = [];

    % Define the decaying exponential model function
    decayingExponentialModel = @(a, b, c, x) (a*(1-exp(-x/b) + c*exp(-x/b)));

    % Set start points, lower, and upper limits for the coefficients
    startPoints = [96,50,80]; % starting values for [a, x0]
    lowerLimits = [90, 10, 0]; % lower limits for [a, x0]
    upperLimits = [100, Inf, 100]; % upper limits for [a, x0]

    % Fit the decaying exponential model to the data with constraints
    mdl = fit(x, y, decayingExponentialModel, 'StartPoint', startPoints, 'Lower', lowerLimits, 'Upper', upperLimits);

    % Extract the coefficients
    a = mdl.a;
    b = mdl.b;
    c = mdl.c;

    % Compute the R-squared value
    y_fit = decayingExponentialModel(a, b, c, x);
    R = corrcoef(y, y_fit);
    R2 = R(1, 2)^2;
end
