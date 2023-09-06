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

addpath(fullfile('..','tools','suptitle'))

%% Initialize stuff

% --- I/O Folders
ioFolder = fullfile('..','data','expC-analysis-spm');
outputFolder = fullfile('..','data','expC-analysis-spm-images');

delay_plot = 0; % in volumes

N = 20; % number of subjects

ERA_G = struct();

Y_AUC_G = struct();

for ss = 1:N
    
    load(fullfile(ioFolder,sprintf('VPIS%.2i_bold.mat',ss)))
    
    for tt = 1:length(trialTestTps)
    
        ERA_G.BilateralMT.(trialTestTps{tt}).psc.mean(ss,:) = ERA.bilateralMT.(trialTestTps{tt}).psc.stats(1,:);
        ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.mean(ss,:) = ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(1,:);
        
        Y_AUC_G.(trialTestTps{tt}).psc.data(ss,:) = Y_AUC.(trialTestTps{tt}).psc.psc_norm;
        Y_AUC_G.(trialTestTps{tt}).psc2static.data(ss,:) = Y_AUC.(trialTestTps{tt}).psc2static.psc_norm;
                
    end
       
end

%% Calculate stats
for tt = 1:length(trialTestTps)
    
    ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.mean = mean(ERA_G.BilateralMT.(trialTestTps{tt}).psc.mean);
    ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.sem = std(ERA_G.BilateralMT.(trialTestTps{tt}).psc.mean) / sqrt(N);

    ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.mean = mean(ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.mean);
    ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.sem = std(ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.mean) / sqrt(N);
    
    Y_AUC_G.(trialTestTps{tt}).psc.stats.mean = mean(Y_AUC_G.(trialTestTps{tt}).psc.data);
    Y_AUC_G.(trialTestTps{tt}).psc.stats.sem = std(Y_AUC_G.(trialTestTps{tt}).psc.data) / sqrt(N);

    Y_AUC_G.(trialTestTps{tt}).psc2static.stats.mean = mean(Y_AUC_G.(trialTestTps{tt}).psc2static.data);
    Y_AUC_G.(trialTestTps{tt}).psc2static.stats.sem = std(Y_AUC_G.(trialTestTps{tt}).psc2static.data) / sqrt(N);    
    
end

%% Ssafety measures
clear ERA Y_AUC

%% Plot Averages for Bilateral MT - MEAN - PSC to mean

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
        ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.mean,...
        ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.sem,...
        '.-','LineWidth',1.5,'MarkerSize',15,'Color',clrMap(tt,:));
    
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
ylabel({'BOLD signal change (%) to the mean'})

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
        ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.mean,...
        ERA_G.BilateralMT.(trialTestTps{tt}).psc.stats.sem,...
        '.-','LineWidth',1.5,'MarkerSize',15,'Color',clrMap(tt,:));
    
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
ylabel({'BOLD signal change (%) to the mean'})

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
print(fig1,fullfile(outputFolder,sprintf('GroupERA_N%i_BilateralMT_D%i_FullTrial_psc',N,delay_plot)),'-dpng')

%% Plot Averages for Bilateral MT - MEAN - Baseline per Run

fig2 = figure('Name','Group ERA','Position',[100 100 1300 1000]);
movegui('center')

xvector = -2:nTrialVols-3;
clrMap = lines;
trialLabels = {'Coh \rightarrow Coh';'InCoh \rightarrow Coh';'NA \rightarrow Coh';'Coh \rightarrow InCoh';'InCoh \rightarrow InCoh';'NA \rightarrow InCoh'};
% ------------------------------------------------------------------------%
s1=subplot(2,1,1);
absMax = 0;
for tt = [1 2 3]
    
    errorbar(xvector,...
        ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.mean,...
        ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.sem,...
        '.-','LineWidth',1.5,'MarkerSize',15,'Color',clrMap(tt,:));
    
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
ylabel({'BOLD signal change (%) to the static condition'})

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
        ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.mean,...
        ERA_G.BilateralMT.(trialTestTps{tt}).psc2static.stats.sem,...
        '.-','LineWidth',1.5,'MarkerSize',15,'Color',clrMap(tt,:));
    
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
ylabel({'BOLD signal change (%) to the static condition'})

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
print(fig2,fullfile(outputFolder,sprintf('GroupERA_N%i_BilateralMT_D%i_FullTrial_psc2static',N,delay_plot)),'-dpng')

%% Prepare The ultimate plot
trialLabels = {'Coh \rightarrow Coh';'InCoh \rightarrow Coh';'NA \rightarrow Coh';'Coh \rightarrow InCoh';'InCoh \rightarrow InCoh';'NA \rightarrow InCoh'};
clrMap = lines;
delay_tc = 3;
x_auc = (15:22) + delay_tc; % 1 before + test block + 1 after
xvector = 0:length(x_auc)-1;
xx = [-1 length(x_auc)];
yy = [-1 1];
lwidth = 1.5;

%% Plot The ultimate plot - PSC
fig_ultimate1 = figure('Name','The Ultimate Figure','units','normalized','outerposition',[0 0 0.6 0.45]);
movegui('center')

% ------------------------------------------------------------------------%
% -- ERA Coherent --------------------------------------------------------%
s1 = subplot(1,2,1);
line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = errorbar(xvector,Y_AUC_G.Coh_aCoh.psc.stats.mean,Y_AUC_G.Coh_aCoh.psc.stats.sem,'Color',clrMap(1,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e2 = errorbar(xvector,Y_AUC_G.Coh_aInCoh.psc.stats.mean,Y_AUC_G.Coh_aInCoh.psc.stats.sem,'Color',clrMap(2,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e3 = errorbar(xvector,Y_AUC_G.Coh_aNA.psc.stats.mean,Y_AUC_G.Coh_aNA.psc.stats.sem,'Color',clrMap(3,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);

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

e1 = errorbar(xvector,Y_AUC_G.InCoh_aCoh.psc.stats.mean,Y_AUC_G.InCoh_aCoh.psc.stats.sem,'Color',clrMap(4,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e2 = errorbar(xvector,Y_AUC_G.InCoh_aInCoh.psc.stats.mean,Y_AUC_G.InCoh_aInCoh.psc.stats.sem,'Color',clrMap(5,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e3 = errorbar(xvector,Y_AUC_G.InCoh_aNA.psc.stats.mean,Y_AUC_G.InCoh_aNA.psc.stats.sem,'Color',clrMap(6,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);

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
print(fig_ultimate1,fullfile(outputFolder,sprintf('GroupERA_N%i_BilateralMT_D%i_2MotionBlock',N,delay_plot)),'-dpng')

%% Plot The ultimate plot - psc2static
fig_ultimate2 = figure('Name','The Ultimate Figure','units','normalized','outerposition',[0 0 0.6 0.45]);
movegui('center')

% ------------------------------------------------------------------------%
% -- ERA Coherent --------------------------------------------------------%
s1 = subplot(1,2,1);
line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = errorbar(xvector,Y_AUC_G.Coh_aCoh.psc2static.stats.mean,Y_AUC_G.Coh_aCoh.psc2static.stats.sem,'Color',clrMap(1,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e2 = errorbar(xvector,Y_AUC_G.Coh_aInCoh.psc2static.stats.mean,Y_AUC_G.Coh_aInCoh.psc2static.stats.sem,'Color',clrMap(2,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e3 = errorbar(xvector,Y_AUC_G.Coh_aNA.psc2static.stats.mean,Y_AUC_G.Coh_aNA.psc2static.stats.sem,'Color',clrMap(3,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);

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

e1 = errorbar(xvector,Y_AUC_G.InCoh_aCoh.psc2static.stats.mean,Y_AUC_G.InCoh_aCoh.psc2static.stats.sem,'Color',clrMap(4,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e2 = errorbar(xvector,Y_AUC_G.InCoh_aInCoh.psc2static.stats.mean,Y_AUC_G.InCoh_aInCoh.psc2static.stats.sem,'Color',clrMap(5,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e3 = errorbar(xvector,Y_AUC_G.InCoh_aNA.psc2static.stats.mean,Y_AUC_G.InCoh_aNA.psc2static.stats.sem,'Color',clrMap(6,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);

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
print(fig_ultimate2,fullfile(outputFolder,sprintf('GroupERA_N%i_BilateralMT_D%i_2MotionBlock_psc2static',N,delay_plot)),'-dpng')

%% Reorganize to input in GraphPad
Y_AUC_GRAPHPAD = struct();

% Iterate on the data points
for jj = 1:8
    
    % Iterate on the trials
    for cc = 1:6
        Y_AUC_GRAPHPAD.(['T' num2str(jj)])(:,cc) = Y_AUC_G.(trialTestTps{cc}).psc2static.data(:,jj);
    end
    
end

%% Calculate point to point t-tests in the full curve
% Bonferroni correction with 3 comparions

TTestResAdaptation = zeros(4,30);
alpha = 0.05;
nComp = 3;

for jj = 1:30
    
   [TTestResAdaptation(1,jj),TTestResAdaptation(2,jj)] = ttest(ERA_G.BilateralMT.Coh_aCoh.psc.mean(:,jj),ERA_G.BilateralMT.Coh_aInCoh.psc.mean(:,jj));
   
   [TTestResAdaptation(3,jj),TTestResAdaptation(4,jj)] = ttest(ERA_G.BilateralMT.InCoh_aCoh.psc.mean(:,jj),ERA_G.BilateralMT.InCoh_aInCoh.psc.mean(:,jj));

   [TTestResAdaptation(5,jj),TTestResAdaptation(6,jj)] = ttest(ERA_G.BilateralMT.Coh_aInCoh.psc.mean(:,jj),ERA_G.BilateralMT.Coh_aNA.psc.mean(:,jj));

   [TTestResAdaptation(7,jj),TTestResAdaptation(8,jj)] = ttest(ERA_G.BilateralMT.InCoh_aCoh.psc.mean(:,jj),ERA_G.BilateralMT.InCoh_aNA.psc.mean(:,jj));
   
   TTestResAdaptation(2,jj) = TTestResAdaptation(2,jj) * nComp;
   TTestResAdaptation(1,jj) = TTestResAdaptation(2,jj) < alpha;

   TTestResAdaptation(4,jj) = TTestResAdaptation(4,jj) * nComp;
   TTestResAdaptation(3,jj) = TTestResAdaptation(4,jj) < alpha;

   TTestResAdaptation(6,jj) = TTestResAdaptation(6,jj) * nComp;
   TTestResAdaptation(5,jj) = TTestResAdaptation(6,jj) < alpha;

   TTestResAdaptation(8,jj) = TTestResAdaptation(8,jj) * nComp;
   TTestResAdaptation(7,jj) = TTestResAdaptation(8,jj) < alpha;
      
end

%% Calculate point to point t-tests
% Bonferroni correction with 3 comparions

TTestRes = zeros(4,8);
alpha = 0.05;
nComp = 3;

for jj = 1:8
    
   [TTestRes(1,jj),TTestRes(2,jj)] = ttest(Y_AUC_G.Coh_aCoh.psc.data(:,jj),Y_AUC_G.Coh_aInCoh.psc.data(:,jj));
   
   [TTestRes(3,jj),TTestRes(4,jj)] = ttest(Y_AUC_G.InCoh_aCoh.psc.data(:,jj),Y_AUC_G.InCoh_aInCoh.psc.data(:,jj));
   
   TTestRes(2,jj) = TTestRes(2,jj) * nComp;
   TTestRes(1,jj) = round(TTestRes(2,jj),2) <= alpha;
   TTestRes(4,jj) = TTestRes(4,jj) * nComp;
   TTestRes(3,jj) = round(TTestRes(4,jj),2) <= alpha;
      
end

%% ANOVA + Multiple comparison test (tukey's) - Adaptation period

ANAdaptation = zeros(4,30); % F and p x 2 in rows, 30 points as columns
TukeysAdaptation = zeros(8,30); % h and p values x 4, 8 points as columns

for jj = 1:30

    [p1,tbl1,stats1] = anova1([ERA_G.BilateralMT.Coh_aCoh.psc.mean(:,jj) ERA_G.BilateralMT.Coh_aInCoh.psc.mean(:,jj) ERA_G.BilateralMT.Coh_aNA.psc.mean(:,jj)],{'Coh_aCoh','Coh_aInCoh','Coh_aNA'});

    f1 = tbl1{2,5};

    [c1,~,~,gnames1] = multcompare(stats1);
    % Columns 1 and 2 contain the indices of the two samples being compared.
    % Column 3 contains the lower confidence interval, column 4 contains the estimate, and column 5 contains the upper confidence interval.
    % Column 6 contains the p-value for the hypothesis test that the corresponding mean difference is not equal to 0.

    [p2,tbl2,stats2] = anova1([ERA_G.BilateralMT.InCoh_aCoh.psc.mean(:,jj) ERA_G.BilateralMT.InCoh_aInCoh.psc.mean(:,jj) ERA_G.BilateralMT.InCoh_aNA.psc.mean(:,jj)],{'InCoh_aCoh','InCoh_aInCoh','InCoh_aNA'});

    f2 = tbl2{2,5};

    [c2,~,~,gnames2] = multcompare(stats2);

    TukeysAdaptation(2,jj) = c1(1,6);
    TukeysAdaptation(1,jj) = round(c1(1,6),2) <= 0.05;

    % to do aNA !!!

    TukeysAdaptation(4,jj) = c2(1,6);
    TukeysAdaptation(3,jj) = round(c2(1,6),2) <= 0.05;

    TukeysAdaptation(8,jj) = c2(2,6);
    TukeysAdaptation(7,jj) = round(c2(2,6),2) <= 0.05;

    TukeysAdaptation(6,jj) = c1(2,6);
    TukeysAdaptation(5,jj) = round(c1(2,6),2) <= 0.05;

    ANAdaptation(1,jj) = f1;
    ANAdaptation(2,jj) = p1;
    ANAdaptation(3,jj) = f2;
    ANAdaptation(4,jj) = p2;

end

%% ANOVA + Multiple comparison test (tukey's)

AN = zeros(4,8); % F and p x 2 in rows, 8 points as columns
Tukeys = zeros(8,8); % h and p values x 4, 8 points as columns

for jj = 1:8

    [p1,tbl1,stats1] = anova1([Y_AUC_G.Coh_aCoh.psc.data(:,jj) Y_AUC_G.Coh_aInCoh.psc.data(:,jj) Y_AUC_G.Coh_aNA.psc.data(:,jj)],{'Coh_aCoh','Coh_aInCoh','Coh_aNA'});

    f1 = tbl1{2,5};

    [c1,~,~,gnames1] = multcompare(stats1);
    % Columns 1 and 2 contain the indices of the two samples being compared.
    % Column 3 contains the lower confidence interval, column 4 contains the estimate, and column 5 contains the upper confidence interval.
    % Column 6 contains the p-value for the hypothesis test that the corresponding mean difference is not equal to 0.

    [p2,tbl2,stats2] = anova1([Y_AUC_G.InCoh_aCoh.psc.data(:,jj) Y_AUC_G.InCoh_aInCoh.psc.data(:,jj) Y_AUC_G.InCoh_aNA.psc.data(:,jj)],{'InCoh_aCoh','InCoh_aInCoh','InCoh_aNA'});

    f2 = tbl2{2,5};

    [c2,~,~,gnames2] = multcompare(stats2);

    Tukeys(2,jj) = c1(1,6);
    Tukeys(1,jj) = round(c1(1,6),2) <= 0.05;

    % to do aNA !!!

    Tukeys(4,jj) = c2(1,6);
    Tukeys(3,jj) = round(c2(1,6),2) <= 0.05;

    Tukeys(8,jj) = c2(2,6);
    Tukeys(7,jj) = round(c2(2,6),2) <= 0.05;

    Tukeys(6,jj) = c1(2,6);
    Tukeys(5,jj) = round(c1(2,6),2) <= 0.05;

    AN(1,jj) = f1;
    AN(2,jj) = p1;
    AN(3,jj) = f2;
    AN(4,jj) = p2;

end

%% Slope analysis
conds2fit = {'Coh_aInCoh','Coh_aNA','InCoh_aCoh','InCoh_aNA'};
FITres = zeros(4,2);
x = 0:3;

for jj = 1:length(conds2fit)

    y = Y_AUC_G.(conds2fit{jj}).psc.stats.mean(1:4);

    p = polyfit(x,y,1);
    
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;

    FITres(jj,1) = p(1);
    FITres(jj,2) = rsq;

end

%% Slope Comparison
x1 = x;
x2 = x;
n1 = 20;
n2 = 20;

y1 = Y_AUC_G.Coh_aInCoh.psc.stats.mean(1:4);
y2 = Y_AUC_G.Coh_aNA.psc.stats.mean(1:4);

[t,p] = slopeComparison(x1,y1,n1,x2,y2,n2);

y1 = Y_AUC_G.InCoh_aCoh.psc.stats.mean(1:4);
y2 = Y_AUC_G.InCoh_aNA.psc.stats.mean(1:4);

[t,p] = slopeComparison(x1,y1,n1,x2,y2,n2);

%% Export workspace
save(fullfile(ioFolder,sprintf('GroupERA_N%i_BilateralMT.mat',N)),'ERA_G','Y_AUC_G','N','TTestRes','Tukeys','TTestResAdaptation','TukeysAdaptation','AN','FITres','nTrialVols','trialTestTps','delay_tc')
