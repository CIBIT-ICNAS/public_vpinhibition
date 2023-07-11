% ======================================================================= %
% == INHIBITION EXPERIMENT SERIES ======================================= %
% EXPERIMENT C
% FMRI DATA ANALYSIS - GROUP - SPM12
% v0.1
%
% AUTHORS:
% VISUAL PERCEPTION GROUP at CIBIT-ICNAS
% 2018 - 2019
% ======================================================================= %
% ======================================================================= %

clear,clc,close all

%% Initialize stuff

% --- I/O Folders
ioFolder = fullfile('...','ICNAS_VisualPerception','Inhibition','expC-analysis-spm');
outputFolder = fullfile('...','ICNAS_VisualPerception','Inhibition','expC-analysis-spm-images');

delay_plot = 0; % in volumes

N = 20; % number of subjects

ERA_G = struct();

Y_AUC_G = struct();

for ss = 1:N
    
    load(fullfile(ioFolder,sprintf('VPIS%.2i_bold.mat',ss)))
    
    for tt = 1:length(trialTestTps)
    
        ERA_G.BilateralMT.(trialTestTps{tt}).mean(ss,:) = ERA.bilateralMT.(trialTestTps{tt}).stats(1,:);
        
        Y_AUC_G.(trialTestTps{tt}).data(ss,:) = Y_AUC.(trialTestTps{tt}).psc_norm;
                
    end
       
end

%% Calculate stats
for tt = 1:length(trialTestTps)
    
    ERA_G.BilateralMT.(trialTestTps{tt}).stats.mean = mean(ERA_G.BilateralMT.(trialTestTps{tt}).mean);
    ERA_G.BilateralMT.(trialTestTps{tt}).stats.sem = std(ERA_G.BilateralMT.(trialTestTps{tt}).mean) / sqrt(N);
    
    Y_AUC_G.(trialTestTps{tt}).stats.mean = mean(Y_AUC_G.(trialTestTps{tt}).data);
    Y_AUC_G.(trialTestTps{tt}).stats.sem = std(Y_AUC_G.(trialTestTps{tt}).data) / sqrt(N);
    
end

%% Ssafety measures
clear ERA Y_AUC

%% Plot Averages for Bilateral MT - MEAN - Baseline per Run

fig1 = figure('Name','Group ERA','Position',[100 100 1300 1000]);
movegui('center')

xvector = -2:nTrialVols-3;
clrMap = lines;
trialLabels = {'Coh \rightarrow Coh';'InCoh \rightarrow Coh';'NA \rightarrow Coh';'Coh \rightarrow InCoh';'InCoh \rightarrow InCoh';'NA \rightarrow InCoh'};
% ------------------------------------------------------------------------%
s1=subplot(2,1,1);
absMax = 0;
for tt = [1 2 3]
    
    errorbar(xvector,...
        ERA_G.BilateralMT.(trialTestTps{tt}).stats.mean,...
        ERA_G.BilateralMT.(trialTestTps{tt}).stats.sem,...
        '.-','LineWidth',1.5,'MarkerSize',15,'Color',clrMap(tt,:));
    
    if max(ERA_G.BilateralMT.(trialTestTps{tt}).stats.mean+ERA_G.BilateralMT.(trialTestTps{tt}).stats.sem) > absMax
        absMax = max(ERA_G.BilateralMT.(trialTestTps{tt}).stats.mean+ERA_G.BilateralMT.(trialTestTps{tt}).stats.sem);
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

text(6,yy(2)-0.2,'Motion','HorizontalAlignment','center','FontSize',12)
text(15.5,yy(2)-0.2,'Coherent','HorizontalAlignment','center','FontSize',12)
text(20,yy(2)-0.2,'MAE','HorizontalAlignment','center','FontSize',12)
text(23,yy(2)-0.2,'Report','HorizontalAlignment','center','FontSize',12)

legend(trialLabels([1 2 3]),'Location','Northwest')
s1.FontSize = 16;
% ------------------------------------------------------------------------%
% ------------------------------------------------------------------------%
s2=subplot(2,1,2);
absMax = 0;
for tt = [4 5 6]
    
    errorbar(xvector,...
        ERA_G.BilateralMT.(trialTestTps{tt}).stats.mean,...
        ERA_G.BilateralMT.(trialTestTps{tt}).stats.sem,...
        '.-','LineWidth',1.5,'MarkerSize',15,'Color',clrMap(tt,:));
    
    if max(ERA_G.BilateralMT.(trialTestTps{tt}).stats.mean+ERA_G.BilateralMT.(trialTestTps{tt}).stats.sem) > absMax
        absMax = max(ERA_G.BilateralMT.(trialTestTps{tt}).stats.mean+ERA_G.BilateralMT.(trialTestTps{tt}).stats.sem);
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

text(6,yy(2)-0.2,'Motion','HorizontalAlignment','center','FontSize',12)
text(15.5,yy(2)-0.2,'Incoherent','HorizontalAlignment','center','FontSize',12)
text(20,yy(2)-0.2,'MAE','HorizontalAlignment','center','FontSize',12)
text(23,yy(2)-0.2,'Report','HorizontalAlignment','center','FontSize',12)

legend(trialLabels([4 5 6]),'Location','Northwest')
s2.FontSize = 16;
% ------------------------------------------------------------------------%

suptitle(sprintf('Group ERA of Bilateral hMT+ - N = %i - Full Trial',N))

%% Save above figure
print(fig1,fullfile(outputFolder,sprintf('GroupERA_N%i_BilateralMT_D%i_FullTrial',N,delay_plot)),'-dpng')

%% Prepare The ultimate plot
trialLabels = {'Coh \rightarrow Coh';'InCoh \rightarrow Coh';'NA \rightarrow Coh';'Coh \rightarrow InCoh';'InCoh \rightarrow InCoh';'NA \rightarrow InCoh'};
clrMap = lines;
x_auc = 15:22; % 1 before + test block + 1 after
xvector = 0:length(x_auc)-1;
xx = [-1 length(x_auc)];
yy = [-1 1];
lwidth = 1.5;

%% Plot The ultimate plot
fig_ultimate = figure('Name','The Ultimate Figure','units','normalized','outerposition',[0 0 0.6 0.45]);
movegui('center')

% ------------------------------------------------------------------------%
% -- ERA Coherent --------------------------------------------------------%
s1 = subplot(1,2,1);
line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = errorbar(xvector,Y_AUC_G.Coh_aCoh.stats.mean,Y_AUC_G.Coh_aCoh.stats.sem,'Color',clrMap(1,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e2 = errorbar(xvector,Y_AUC_G.Coh_aInCoh.stats.mean,Y_AUC_G.Coh_aInCoh.stats.sem,'Color',clrMap(2,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e3 = errorbar(xvector,Y_AUC_G.Coh_aNA.stats.mean,Y_AUC_G.Coh_aNA.stats.sem,'Color',clrMap(3,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);

hold off

xlim(xx), ylim(yy)
xticks(xvector)
legend([e1 e2 e3],trialLabels(1:3),'FontSize',12)
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal change (%)'},'FontSize',12)
title('ERA')
s1.FontSize = 14;
% ------------------------------------------------------------------------%
% -- ERA InCoherent ------------------------------------------------------%
s4 = subplot(1,2,2);
line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = errorbar(xvector,Y_AUC_G.InCoh_aCoh.stats.mean,Y_AUC_G.InCoh_aCoh.stats.sem,'Color',clrMap(4,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e2 = errorbar(xvector,Y_AUC_G.InCoh_aInCoh.stats.mean,Y_AUC_G.InCoh_aInCoh.stats.sem,'Color',clrMap(5,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e3 = errorbar(xvector,Y_AUC_G.InCoh_aNA.stats.mean,Y_AUC_G.InCoh_aNA.stats.sem,'Color',clrMap(6,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);

hold off

xlim(xx), ylim(yy)
xticks(xvector)
legend([e1 e2 e3],trialLabels(4:6),'FontSize',12)
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal change (%)'},'FontSize',12)
s4.FontSize = 14;
% ------------------------------------------------------------------------%

%% Export The ultimate Figure
print(fig_ultimate,fullfile(outputFolder,sprintf('GroupERA_N%i_BilateralMT_D%i_2MotionBlock',N,delay_plot)),'-dpng')

%% Reorganize to input in GraphPad
Y_AUC_GRAPHPAD = struct();

% Iterate on the data points
for jj = 1:8
    
    % Iterate on the trials
    for cc = 1:6
        Y_AUC_GRAPHPAD.(['T' num2str(jj)])(:,cc) = Y_AUC_G.(trialTestTps{cc}).data(:,jj);
    end
    
end

%% Calculate point to point t-tests
% Bonferroni correction with 3 comparions

TTestRes = zeros(4,8);
alpha = 0.05;
nComp = 3;

for jj = 1:8
    
   [TTestRes(1,jj),TTestRes(2,jj)] = ttest(Y_AUC_G.Coh_aCoh.data(:,jj),Y_AUC_G.Coh_aInCoh.data(:,jj));
   
   [TTestRes(3,jj),TTestRes(4,jj)] = ttest(Y_AUC_G.InCoh_aCoh.data(:,jj),Y_AUC_G.InCoh_aInCoh.data(:,jj));
   
   TTestRes(2,jj) = TTestRes(2,jj) * nComp;
   TTestRes(1,jj) = TTestRes(2,jj) < alpha;
   TTestRes(4,jj) = TTestRes(4,jj) * nComp;
   TTestRes(3,jj) = TTestRes(4,jj) < alpha;
      
end

%% Export workspace
save(fullfile(ioFolder,sprintf('GroupERA_N%i_BilateralMT.mat',N)),'ERA_G','Y_AUC_G','N','TTestRes')
