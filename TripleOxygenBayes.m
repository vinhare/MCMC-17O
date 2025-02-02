%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
%%     Triple Oxygen Isotope Parameter estimation by MCMC  %
%                   By Vincent Hare, December 2023                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uses MCMC toolbox (mjlaine.github.io/mcmcstat/) and DERIVESTsuite
% also requires 4 external functions: BWM_Emu_st, TripleOxygenss,
% Body_Water_Model_rh, and SS_over_rh
%
% Parameters are as follows:
% theta = [DeltaOatm, deltaOmw, rh, Leafratio, WEIndex, d18O_atmO2];
%
% Dataset "NAME.mat" has the stucture:
%   NAME.ydata : D'17Obw     
%   NAME.xdata : d'18Obw     
%   NAME.yerr  : uD'17Obw (1 sigma)     
% 
% Note: For non-ostriches, parameters need to be updated in both both the 
%       BWM_Emu_st.m and Body_Water_Model_rh.m files (e.g. mass is 84 kg by
%       default; it needs to be changed in both files for accurate fits.)
%% 
clear all; %data model params options final_D17O fval A B C AllAs AllBs
 
close all;

load TripleO_Central_Quaternary.mat % a .mat file with data.ydata, data.xdata
%data=Central_LateMiocene.data;

tmin = [-0.434 -4.0 0.5 0.75 0.17 24.046]; %initial guess for theta
nsimu = 50000;  % number of simulations
method      = 'dram'; % adaptation method, 'mh', 'dr', 'am', or 'dram'
adaptint    = 50;    % how often to adapt the proposal for am or dram
 
n = length(data.xdata);
p = 3; %number of parameters
 
mse = TripleOxygenss(tmin,data)/(n-p);  % estimate for the error variance

% The following lines define the PRIORS for theta, as well as whether they are targets for the MCMC (target? = 1 means that the posterior distribution will be estimated from the prior. target? = 0 indicates marginalised prior (i.e. the posterior is not estimated).    
% theta(1) = DeltaOatm, theta(2) = deltaOmw, theta(3) = rh, theta(4) =
% Leafratio, theta(5) = WEIndex, theta(6) = d18O_atmO2;
params = {
  %      name,    init,min, max, mu, sig, target?, local?
    {'theta1', tmin(1), -1.5, -0.100, tmin(1), 0.022,   0,      1}
    {'theta2', tmin(2), -30, 15, tmin(2), 3,   1,      1} 
    {'theta3', tmin(3), 0.1, 0.9, tmin(3), Inf, 0, 1}
    {'theta4', tmin(4), 0.0, 1.0, tmin(4), Inf, 1, 1}
    {'theta5', tmin(5), 0.05, 0.6, tmin(5), Inf, 1, 1}
    {'theta6', tmin(6), 20, 27, tmin(6), 0.1, 0, 1}
    };
 
model.ssfun = @TripleOxygenss
model.sigma2 = 1; 
 
options.updatesigma = 1;
 
options.method      = method
options.nsimu       = nsimu; 
%options.qcov        = tcov; %with a simpler model, it would be possible to
%calculate the Jacobian, and propose a covariance matrix, which might improve 
%things somewhat, but as things currently stand the Jacobian is difficult
%to calculate, so I have omitted it. 
 
[res,chain,s2chain] = mcmcrun(model,data,params,options);

%% PLOTS

width = 3 %3.4252;     % Width in inches. 3 is good for multifigures. 3.4252 is golden ratio
height = 2.1168;    % Height in inches
alw = 1.5;    % AxesLineWidth
fsz = 8;      % Fontsize
lw = 0.6;      % LineWidth
msz = 6;       % MarkerSize
set(0,'defaultTextInterpreter','latex')
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
fontname = 'Helvetica';
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);
color1=[215/255 25/255 28/255]; %red
color2=[253/255 174/255 97/255];  %orange
color3=[26/255 150/255 65/255];  %green
color4=[43/255 131/255 186/255];  %blue
 
figure(1); clf
mcmcplot(chain,[],res,'chainpanel');
subplot(2,2,4)
mcmcplot(sqrt(s2chain),[],[],'dens',2);
title('error std');
 
figure(2); clf
mcmcplot(chain,[],res,'denspanel',2);
 
figure(3); clf
subplot(3,3,4);
plot1=dscatter(chain(:,1),chain(:,2));
xLimits1 = get(gca,'XLim');
yLimits1 = get(gca,'YLim');
hold on;
plot([mean(chain(:,1)) mean(chain(:,1))],[yLimits1],'k')
plot([xLimits1],[mean(chain(:,2)) mean(chain(:,2))],'k')
plot([mean(chain(:,1))],[mean(chain(:,2))],'sr')
hold off;
box on;
grid on;
set(gca,'TickDir','out');
ylabel('$r_{l/s}$','Interpreter','latex');
subplot(3,3,7);
plot2=dscatter(chain(:,1),chain(:,3));
hold on;
plot([mean(chain(:,1)) mean(chain(:,1))],[min(chain(:,3)) max(chain(:,3))],'k')
plot([xLimits1],[mean(chain(:,3)) mean(chain(:,3))],'k')
plot([mean(chain(:,1))],[mean(chain(:,3))],'sr')
yLimits2 = get(gca,'YLim');
hold off;
box on;
grid on;
set(gca,'TickDir','out');
ylabel('$\mathrm{WEI}~~(\mathrm{mL/kJ})$','Interpreter','latex');
xlabel('$\delta^{18}\mathrm{O}_{\mathrm{mw}}$','Interpreter','latex');
subplot(3,3,8);
plot3=dscatter(chain(:,2),chain(:,3));
hold on;
xLimits2 = get(gca,'XLim');
xLimits3 = get(gca,'YLim');
plot([mean(chain(:,2)) mean(chain(:,2))],[min(chain(:,3)) max(chain(:,3))],'k')
plot([yLimits1],[mean(chain(:,3)) mean(chain(:,3))],'k')
plot([mean(chain(:,2))],[mean(chain(:,3))],'sr')
xlabel('$r_{l/s}$','Interpreter','latex');
box on;
grid on;
set(gca,'TickDir','out');
subplot(3,3,1);
[f1,xi1] = ksdensity(chain(100:length(chain),1)); 
plot(xi1,f1);
%histogram(chain(100:length(chain),1),'Normalization','probability','DisplayStyle','stairs');
%uncomment the above for a "stair" type histogram, instead of a KDE.
hold on;
plot([mean(chain(:,1)) mean(chain(:,1))],[0 max(f1)],'k')
plot([mean(chain(:,1))],[max(f1)],'sr')
xlim(xLimits1)
set(gca,'YTickLabel',[]);
set(gca,'TickDir','out');
title('$\delta^{18}\mathrm{O}_{\mathrm{mw}}$','Interpreter','latex');
subplot(3,3,5)
[f2,xi2] = ksdensity(chain(100:length(chain),2)); 
plot(xi2,f2);
%histogram(chain(100:length(chain),2),'Normalization','probability','DisplayStyle','stairs');
%uncomment the above for a "stair" type histogram, instead of a KDE.
hold on;
plot([median(chain(:,2)) median(chain(:,2))],[0 max(f2)],'k')
plot([median(chain(:,2))],[max(f2)],'sr')
xlim(xLimits2)
set(gca,'YTickLabel',[])
set(gca,'TickDir','out');
title('$r_{l/s}$','Interpreter','latex');
subplot(3,3,9)
[f3,xi3] = ksdensity(chain(100:length(chain),3)); 
plot(xi3,f3);
%histogram(chain(100:length(chain),3),'Normalization','probability','DisplayStyle','stairs');
%uncomment the above for a "stair" type histogram, instead of a KDE.
hold on;
plot([mean(chain(:,3)) mean(chain(:,3))],[0 max(f3)],'k')
plot([mean(chain(:,3))],[max(f3)],'sr')
xlim(xLimits3);
set(gca,'YTickLabel',[]);
set(gca,'TickDir','out');
title('$\mathrm{WEI}~~(\mathrm{mL/kJ})$','Interpreter','latex');

%output the chain statistics
chainstats(chain,res)

%%Uncomment the following for additional optimisation of D'17O
%initial guess
init_D17O=tmin(1);
options = optimset('PlotFcns',@optimplotfval);
fun = @(x) SS_over_rh(x, chain, data);
[final_D17O, fval] = fminsearch(fun, init_D17O); 

%Bootstrapping to estimate error in final_D17O
numBoot = 100;  % Number of bootstrap iterations (adjust for precision vs. time)
D17O_samples = zeros(numBoot,1);
nData = size(data.xdata,1);
for i = 1:numBoot;
    idx = randi(nData, nData, 1);
    data_boot = data.xdata(idx,:);
    data_boot(:,1) = [1; 2; 3; 4]; % THIS NEEDS To be updated for variable length... 
    datanew.xdata=data_boot; 
    datanew.ydata=data.ydata;
    datanew.yerr=data.yerr;
    datanew.xerr=data.xerr; 
    objectiveFun_boot = @(x) SS_over_rh(x, chain, datanew);
    D17O_samples(i) = fminsearch(objectiveFun_boot, init_D17O);
end;
D17O_error = std(D17O_samples);

figure(4); clf
%Uncomment the following for model prediction (right now, this adds time!)
chain1=chain(1000:50:end,1); %thin the chains a bit 
chain2=chain(1000:50:end,2);
chain3=chain(1000:50:end,3);
for i = 1:length(chain1);    
[A B C] = Body_Water_Model_rh(final_D17O,chain1(i,1), chain2(i,1), chain3(i,1));
AllAs(:,i) = A(1:length(A),1);    %discard burn-in for plots
AllBs(:,i) = B(1:length(B),1);   %discard burn-in for plots
end; 
plot(AllBs,1000.*AllAs,'color',[0 0 1 0.05]);


[A B C] = Body_Water_Model_rh(final_D17O,median(chain1(:,1)), median(chain2(:,1)), median(chain3(:,1)));
hold on
plot(B,1000.*A,'color',[0 0 1],'linewidth',1.5)
errorbar(data.xdata(:,2),1000.*data.ydata,1000.*data.yerr,'ok','MarkerFaceColor','w','MarkerSize',9)
title('$\mathrm{Southern~Namib~-~Late~Miocene}$','Interpreter','latex','FontSize',12);
set(gca,'TickDir','out');
results = analyzeMatrix(AllBs,AllAs, 0.1, 200);
targets = [results.target]';
medians = [results.median]';
confIntLower = arrayfun(@(x) x.confInt(1), results)';
confIntUpper = arrayfun(@(x) x.confInt(2), results)';
%plot(targets,confIntLower);
%plot(targets,confIntUpper);

xlabh=xlabel('$\delta^{18}\mathrm{O_{bw}},$~~~~$\mathrm{vs~VSMOW}$');
ylabh=ylabel('$\Delta''^{17}\mathrm{O_{bw}},$~$\mathrm{per~meg~vs~VSMOW}$');
xy=xlabh.Extent;
text(5, -300, char(8240),'FontName', 'Times','Rotation',0,'HorizontalAlignment','center','Interpreter','none'); 

ylim([-350, 80]);
xlim([-5, 30]);

fprintf('Final estimated atmospheric D17O = %.4f with an error of %.4f (from %d bootstrap iterations)\n', ...
  final_D17O, D17O_error, numBoot);

function indices = findWithinTolerance(M, target, a)
    % Create a logical matrix where the condition is true
    condition = (M > target - a) & (M < target + a);
    
    % Find the indices of elements that satisfy the condition
    [row, col] = find(condition);
    
    % Combine the row and column indices into a single matrix
    indices = [row, col];
end

function results = analyzeMatrix(A, B, a, numSteps)
    % Initialize results structure
    results = struct('target', [], 'median', [], 'confInt', []);
    
    % Determine the range of target values
    minVal = min(A(:));
    maxVal = max(A(:));
    
    % Define the range of target values
    targetValues = linspace(minVal, maxVal, numSteps);
    
    % Loop over each target value
    for i = 1:length(targetValues)
        target = targetValues(i);
        
        % Find indices in A within tolerance of target
        indices = findWithinTolerance(A, target, a);
        
        % Extract corresponding values from B
        if isempty(indices)
            V = [];
        else
            V = arrayfun(@(r, c) B(r, c), indices(:, 1), indices(:, 2));
        end
        
        % Calculate median and 90% confidence intervals
        if isempty(V)
            med = NaN;
            confInt = [NaN NaN];
        else
            med = median(V);
            confInt = prctile(V, [5, 95]);
        end
        
        % Store results
        results(i).target = target;
        results(i).median = med;
        results(i).confInt = confInt;
    end
end
