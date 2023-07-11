% ======================================================================= %
% == INHIBITION EXPERIMENT SERIES ======================================= %
% EXPERIMENT C
% FMRI DATA ANALYSIS - GROUP FIGURES 2 - SPM12
% v0.1
%
% AUTHORS:
% VISUAL PERCEPTION GROUP at CIBIT-ICNAS
% 2018 - 2019
% ======================================================================= %
% ======================================================================= %

clear,clc,close all

% RUN AFTER GROUP ANALYSIS
addpath(fullfile('..','tools','shadedErrorBar'))

% Load data
load(fullfile('..','data','expC-analysis-spm','GroupERA_N20_BilateralMT.mat'))

% Out folder
outputFolder = fullfile('..','data','expC-analysis-spm-images');

%% Figure 1 - PSC
fig1 = figure('Name','Group ERA','Position',[100 100 1100 900]);
xvector = -2:nTrialVols-3;
clrMap = lines;
trialLabels = {'Coherent \rightarrow Coherent';'Incoherent \rightarrow Coherent';'Non-adapting \rightarrow Coherent';'Coherent \rightarrow Incoherent';'Incoherent \rightarrow Incoherent';'Non-adapting \rightarrow Incoherent'};
% ------------------------------------------------------------------------%
s1=subplot(2,1,1);
absMax = 0;
for tt = [1 2 3]
    
    shadedErrorBar(xvector,...
        ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.mean,...
        ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.sem,...
        'lineprops',{'.-','Color',clrMap(tt,:),'LineWidth',1.5,'MarkerSize',10});
    
    if max(ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.mean+ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.sem) > absMax
        absMax = max(ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.mean+ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.sem);
    end
    
    hold on    
end

hold off
xx = [xvector(1)-1 xvector(end)+1];
yy = [-2 ceil(absMax+0.5)];
xlim(xx); ylim(yy);
xlabel('Time (volumes)')
ylabel({'BOLD signal change (%)'})

xticks([1 6 12 13 18 19 21 22 24])

line(xx,zeros(length(xx),1),'LineStyle',':','Color','k')

line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([12.5 12.5],yy,'LineStyle','--','Color','k')
line([18.5 18.5],yy,'LineStyle','--','Color','k')
line([21.5 21.5],yy,'LineStyle','--','Color','k')
line([24.5 24.5],yy,'LineStyle','--','Color','k')

text(6,yy(2)-0.4,'First motion period','HorizontalAlignment','center','FontSize',14)
text(15.5,yy(2)-0.4,{'Test with','Coherent motion'},'HorizontalAlignment','center','FontSize',14)

legend(trialLabels([1 2 3]),'fontSize',14)
s1.FontSize = 14;
% ------------------------------------------------------------------------%
% ------------------------------------------------------------------------%
s2=subplot(2,1,2);
absMax = 0;
auxclrmap = [2 1 3];
for tt = [5 4 6]
    
    shadedErrorBar(xvector,...
        ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.mean,...
        ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.sem,...
        'lineprops',{'.-','Color',clrMap(auxclrmap(tt-3),:),'LineWidth',1.5,'MarkerSize',10});

    if max(ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.mean+ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.sem) > absMax
        absMax = max(ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.mean+ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.sem);
    end
    
    hold on    
end

hold off
xx = [xvector(1)-1 xvector(end)+1];
yy = [-2 ceil(absMax+0.5)];
xlim(xx); ylim(yy);
xlabel('Time (volumes)')
ylabel({'BOLD signal change (%)'})

xticks([1 6 12 13 18 19 21 22 24])

line(xx,zeros(length(xx),1),'LineStyle',':','Color','k')

line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([12.5 12.5],yy,'LineStyle','--','Color','k')
line([18.5 18.5],yy,'LineStyle','--','Color','k')
line([21.5 21.5],yy,'LineStyle','--','Color','k')
line([24.5 24.5],yy,'LineStyle','--','Color','k')

text(6,yy(2)-0.4,'First motion period','HorizontalAlignment','center','FontSize',14)
text(15.5,yy(2)-0.4,{'Test with','Incoherent motion'},'HorizontalAlignment','center','FontSize',14)

legend(trialLabels([5 4 6]),'fontSize',14)
s2.FontSize = 14;
% ------------------------------------------------------------------------%

%% Export fig
print(fig1,fullfile(outputFolder,sprintf('FiguresPaper_1_psc.png')),'-dpng')

%% Figure 11 - psc2static
fig11 = figure('Name','Group ERA','Position',[100 100 1100 900]);
xvector = -2:nTrialVols-3;
clrMap = lines;
trialLabels = {'Coherent \rightarrow Coherent';'Incoherent \rightarrow Coherent';'Non-adapting \rightarrow Coherent';'Coherent \rightarrow Incoherent';'Incoherent \rightarrow Incoherent';'Non-adapting \rightarrow Incoherent'};
% ------------------------------------------------------------------------%
s1=subplot(2,1,1);
absMax = 0;
for tt = [1 2 3]
    
    shadedErrorBar(xvector,...
        ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.mean,...
        ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.sem,...
        'lineprops',{'.-','Color',clrMap(tt,:),'LineWidth',1.5,'MarkerSize',10});
    
    if max(ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.mean+ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.sem) > absMax
        absMax = max(ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.mean+ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.sem);
    end
    
    hold on    
end

hold off
xx = [xvector(1)-1 xvector(end)+1];
yy = [-1 ceil(absMax+0.5)];
xlim(xx); ylim(yy);
xlabel('Time (volumes)')
ylabel({'BOLD signal change (%)'})

xticks([1 6 12 13 18 19 21 22 24])

line(xx,zeros(length(xx),1),'LineStyle',':','Color','k')

line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([12.5 12.5],yy,'LineStyle','--','Color','k')
line([18.5 18.5],yy,'LineStyle','--','Color','k')
line([21.5 21.5],yy,'LineStyle','--','Color','k')
line([24.5 24.5],yy,'LineStyle','--','Color','k')

text(6,yy(2)-0.4,'First motion period','HorizontalAlignment','center','FontSize',14)
text(15.5,yy(2)-0.4,{'Test with','Coherent motion'},'HorizontalAlignment','center','FontSize',14)

legend(trialLabels([1 2 3]),'fontSize',14)
s1.FontSize = 14;
% ------------------------------------------------------------------------%
% ------------------------------------------------------------------------%
s2=subplot(2,1,2);
absMax = 0;
auxclrmap = [2 1 3];
for tt = [5 4 6]
    
    shadedErrorBar(xvector,...
        ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.mean,...
        ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.sem,...
        'lineprops',{'.-','Color',clrMap(auxclrmap(tt-3),:),'LineWidth',1.5,'MarkerSize',10});

    if max(ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.mean+ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.sem) > absMax
        absMax = max(ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.mean+ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.sem);
    end
    
    hold on    
end

hold off
xx = [xvector(1)-1 xvector(end)+1];
yy = [-1 ceil(absMax+0.5)];
xlim(xx); ylim(yy);
xlabel('Time (volumes)')
ylabel({'BOLD signal change (%)'})

xticks([1 6 12 13 18 19 21 22 24])

line(xx,zeros(length(xx),1),'LineStyle',':','Color','k')

line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([12.5 12.5],yy,'LineStyle','--','Color','k')
line([18.5 18.5],yy,'LineStyle','--','Color','k')
line([21.5 21.5],yy,'LineStyle','--','Color','k')
line([24.5 24.5],yy,'LineStyle','--','Color','k')

text(6,yy(2)-0.4,'First motion period','HorizontalAlignment','center','FontSize',14)
text(15.5,yy(2)-0.4,{'Test with','Incoherent motion'},'HorizontalAlignment','center','FontSize',14)

legend(trialLabels([5 4 6]),'fontSize',14)
s2.FontSize = 14;
% ------------------------------------------------------------------------%

%% Export fig
print(fig11,fullfile(outputFolder,sprintf('FiguresPaper_1_psc2static.png')),'-dpng')

%% Figure 2 - psc2static
fig2 = figure('position',[150 150 1200 450]);
trialLabels = {'Coherent \rightarrow Coherent';'Incoherent \rightarrow Coherent';'Non-adapting \rightarrow Coherent';'Coherent \rightarrow Incoherent';'Incoherent \rightarrow Incoherent';'Non-adapting \rightarrow Incoherent'};
clrMap = lines;
x_auc = (15:22) + delay_tc; % 1 before + test block + 1 after
xvector = 0:length(x_auc)-1;
xx = [-1 length(x_auc)];
yy = [-0.5 0.7];
yvector = sort([yy(1):0.2:yy(2) 0]);
lwidth = 1.5;

s1 = subplot(1,2,1);
line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = shadedErrorBar(xvector,Y_AUC_G.Coh_aCoh.psc2static.stats.mean,Y_AUC_G.Coh_aCoh.psc2static.stats.sem,'lineprops',{'Color',clrMap(1,:),'LineWidth',lwidth,'Marker','.','MarkerSize',20});
e2 = shadedErrorBar(xvector,Y_AUC_G.Coh_aInCoh.psc2static.stats.mean,Y_AUC_G.Coh_aInCoh.psc2static.stats.sem,'lineprops',{'Color',clrMap(2,:),'LineWidth',lwidth,'Marker','.','MarkerSize',20});
e3 = shadedErrorBar(xvector,Y_AUC_G.Coh_aNA.psc2static.stats.mean,Y_AUC_G.Coh_aNA.psc2static.stats.sem,'lineprops',{'Color',clrMap(3,:),'LineWidth',lwidth,'Marker','.','MarkerSize',20});

hold on

plot(find(TTestRes(1,:))-1, 0.55*ones(1,length(find(TTestRes(1,:)))),'*','Color','k','LineWidth',1,'MarkerSize',8)

hold off

xlim(xx), ylim(yy)
xticks(xvector)
yticks(yvector)
legend([e1.mainLine e2.mainLine e3.mainLine],trialLabels(1:3),'FontSize',14,'Location','Southwest')
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal change (%)'},'FontSize',12)
s1.FontSize = 15;

%-------------------------------------------------------
s4 = subplot(1,2,2);
line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = shadedErrorBar(xvector,Y_AUC_G.InCoh_aInCoh.psc2static.stats.mean,Y_AUC_G.InCoh_aInCoh.psc2static.stats.sem,'lineprops',{'Color',clrMap(1,:),'LineWidth',lwidth,'Marker','.','MarkerSize',20});
e2 = shadedErrorBar(xvector,Y_AUC_G.InCoh_aCoh.psc2static.stats.mean,Y_AUC_G.InCoh_aCoh.psc2static.stats.sem,'lineprops',{'Color',clrMap(2,:),'LineWidth',lwidth,'Marker','.','MarkerSize',20});
e3 = shadedErrorBar(xvector,Y_AUC_G.InCoh_aNA.psc2static.stats.mean,Y_AUC_G.InCoh_aNA.psc2static.stats.sem,'lineprops',{'Color',clrMap(3,:),'LineWidth',lwidth,'Marker','.','MarkerSize',20});

hold on

plot(find(TTestRes(3,:))-1,0.55*ones(1,length(find(TTestRes(3,:)))),'*','Color','k','LineWidth',1,'MarkerSize',8)

hold off

xlim(xx), ylim(yy)
xticks(xvector)
yticks(yvector)
legend([e1.mainLine e2.mainLine e3.mainLine],trialLabels([5 4 6]),'FontSize',14,'Location','southwest')
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal change (%)'},'FontSize',12)
s4.FontSize = 15;

%% Export fig
print(fig2,fullfile(outputFolder,sprintf('FiguresPaper_2_psc2static.png')),'-dpng')

%% Figure 2 - psc
fig22 = figure('position',[150 150 1200 450]);
trialLabels = {'Coherent \rightarrow Coherent';'Incoherent \rightarrow Coherent';'Non-adapting \rightarrow Coherent';'Coherent \rightarrow Incoherent';'Incoherent \rightarrow Incoherent';'Non-adapting \rightarrow Incoherent'};
clrMap = lines;
x_auc = (15:22) + delay_tc; % 1 before + test block + 1 after
xvector = 0:length(x_auc)-1;
xvectorlabel = 15:22;
xx = [-1 length(x_auc)];
yy = [-0.5 0.7];
yvector = sort([yy(1):0.2:yy(2) 0]);
lwidth = 1.5;

s1 = subplot(1,2,1);
line([0.5 0.5],yy,'LineStyle',':','Color','k')
%line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = shadedErrorBar(xvector,Y_AUC_G.Coh_aCoh.psc.stats.mean,Y_AUC_G.Coh_aCoh.psc.stats.sem,'lineprops',{'Color',clrMap(1,:),'LineWidth',lwidth,'Marker','.','MarkerSize',20});
e2 = shadedErrorBar(xvector,Y_AUC_G.Coh_aInCoh.psc.stats.mean,Y_AUC_G.Coh_aInCoh.psc.stats.sem,'lineprops',{'Color',clrMap(2,:),'LineWidth',lwidth,'Marker','.','MarkerSize',20});
e3 = shadedErrorBar(xvector,Y_AUC_G.Coh_aNA.psc.stats.mean,Y_AUC_G.Coh_aNA.psc.stats.sem,'lineprops',{'Color',clrMap(3,:),'LineWidth',lwidth,'Marker','.','MarkerSize',20});

hold on

plot(find(TTestRes(1,:))-1, 0.55*ones(1,length(find(TTestRes(1,:)))),'*','Color','k','LineWidth',1,'MarkerSize',8)

hold off

xlim(xx), ylim(yy)
xticks(xvector), xticklabels(xvectorlabel)
yticks(yvector)
legend([e1.mainLine e2.mainLine e3.mainLine],trialLabels(1:3),'FontSize',14,'Location','Southwest')
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal change (%)'},'FontSize',12)
s1.FontSize = 15;

%-------------------------------------------------------
s4 = subplot(1,2,2);
line([0.5 0.5],yy,'LineStyle',':','Color','k')
%line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = shadedErrorBar(xvector,Y_AUC_G.InCoh_aInCoh.psc.stats.mean,Y_AUC_G.InCoh_aInCoh.psc.stats.sem,'lineprops',{'Color',clrMap(1,:),'LineWidth',lwidth,'Marker','.','MarkerSize',20});
e2 = shadedErrorBar(xvector,Y_AUC_G.InCoh_aCoh.psc.stats.mean,Y_AUC_G.InCoh_aCoh.psc.stats.sem,'lineprops',{'Color',clrMap(2,:),'LineWidth',lwidth,'Marker','.','MarkerSize',20});
e3 = shadedErrorBar(xvector,Y_AUC_G.InCoh_aNA.psc.stats.mean,Y_AUC_G.InCoh_aNA.psc.stats.sem,'lineprops',{'Color',clrMap(3,:),'LineWidth',lwidth,'Marker','.','MarkerSize',20});

hold on

plot(find(TTestRes(3,:))-1,0.55*ones(1,length(find(TTestRes(3,:)))),'*','Color','k','LineWidth',1,'MarkerSize',8)

hold off

xlim(xx), ylim(yy)
xticks(xvector), xticklabels(xvectorlabel)
yticks(yvector)
legend([e1.mainLine e2.mainLine e3.mainLine],trialLabels([5 4 6]),'FontSize',14,'Location','southwest')
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal change (%)'},'FontSize',12)
s4.FontSize = 15;

%% Export fig
print(fig22,fullfile(outputFolder,sprintf('FiguresPaper_2_psc.png')),'-dpng')
