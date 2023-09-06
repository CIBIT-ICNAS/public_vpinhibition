%% -- Step03_CalculateERAs.m ---------------------------------------- %%
% ----------------------------------------------------------------------- %
% Script for executing the third step regarding the SPM analysis of fMRI
% data after preprocessing with fmriPrep v21 - reading the time course
% of a number of ROIs and plotting the event-related average curves.
%
% Dataset:
% - Inhibition (Visual Perception)
%
% Warnings:
% - a number of values/steps are custom for this dataset - full code review
% is strongly advised for different datasets
% - this was designed to run on sim01 - a lot of paths must change if run
% at any other computer
% - only compatible with space MNI152NLin2009cAsym
%
% Requirements:
% - SPM12 in path
% - shadedErrorBar
%
% Author: Alexandre Sayal
% CIBIT, University of Coimbra
% February 2022
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

clear,clc

outputFolder = fullfile('..','data','expC-fmriprep-images');

%% Load packages
addpath('/SCRATCH/software/toolboxes/shadedErrorBar/')

%% Load data
% TC is the struct of interest, containing the time-course for all subs,
% runs and ROIs. This time course is demeaned BOLD signal.
% It is TC.ROI(subjects x runs x volumes)
load('TimeSeries_PeakCoordinates.mat')

%% Define hemodynamic delay
% In this case, SPM already compensates for the lag??
hemoDelay = 0; % the TR is 1 second, so I define 4 volumes as an estimate for the hemodynamic delay
hemoDelayBaseline = 3;

%% Fetch baseline indexes
% I know that all runs have the same trial structure, so the 'Static'
% conditions aka baseline have the same volumes for all six runs.

%tsvData = importTSVProtocol(fullfile(bidsFolder,'sub-01','ses-01','func','sub-01_ses-01_task-inhib_run-1_events.tsv'));
tsvData = importTSVProtocol(fullfile(pwd,'tsv','sub-01_ses-01_task-inhib_run-1_events.tsv'));

baselineOnsets = tsvData.onset(tsvData.trial_type == 'Static') + 1; % onsets are zero-based, but we need indexes, which in Matlab are one-based

[t1,t2] = ndgrid(baselineOnsets, 0:5); % I know that the 'Static' duration is 6 volumes

baselineVols = t1 + t2 + hemoDelayBaseline;

%% Calculate Percent Signal Change to Baseline
PSC = struct();

%PSC.leftMT = (TC.leftMT - mean(TC.leftMT(:,:,baselineVols) , 3) )  ./ std(TC.leftMT(:,:,baselineVols), [], 3);
%PSC.rightMT = (TC.rightMT - mean(TC.rightMT(:,:,baselineVols) , 3)  ) ./ std(TC.rightMT(:,:,baselineVols), [], 3);
%PSC.leftSPL = (TC.leftSPL - mean(TC.leftSPL(:,:,baselineVols) , 3) )  ./ std(TC.leftSPL(:,:,baselineVols), [], 3);
%PSC.rightSPL = (TC.rightSPL - mean(TC.rightSPL(:,:,baselineVols) , 3)  ) ./ std(TC.rightSPL(:,:,baselineVols), [], 3);

PSC.leftMT = TC.leftMT;
PSC.rightMT = TC.rightMT;
PSC.leftSPL = TC.leftSPL;
PSC.rightSPL = TC.rightSPL;

%% Fetch indexes for all trial types
trialIndexes = struct();

trialTypes = {'Coh_aCoh','Coh_aInCoh','Coh_aNA','InCoh_aCoh','InCoh_aInCoh','InCoh_aNA'};
nTrialTypes = length(trialTypes);

for rr = 1 : nTasks
    
%     tsvData = importTSVProtocol(fullfile(bidsFolder,'sub-01','ses-01','func',...
%         sprintf('sub-01_ses-01_task-inhib_run-%i_events.tsv',rr)));

    tsvData = importTSVProtocol(fullfile(pwd,'tsv',...
        sprintf('sub-01_ses-01_task-inhib_run-%i_events.tsv',rr)));
    
    for tt = 1:nTrialTypes
        
        if rr == 1 % initialize matrices
            trialIndexes.(trialTypes{tt}) = zeros(nTasks,2,30); % To plot, I want 30 volumes, and each trial occours twice per run
        end
        
        auxOnsets = tsvData.onset(tsvData.trial_type == trialTypes{tt}) + 1; % onsets are zero-based, but we need indexes, which in Matlab are one-based
        
        [t1,t2] = ndgrid(auxOnsets,-15:14); % include volumes before and after the trial type condition, totalling 30 volumes
        
        trialIndexes.(trialTypes{tt})(rr,:,:) = t1 + t2 + hemoDelay;
        
    end
    
end

%% Calculate ERAs

ERA = struct();

for tt = 1:nTrialTypes
    
    ERA.leftMT.(trialTypes{tt}) = zeros(nSubjects,nTasks,30); % initialize - six runs, 30 volumes
    ERA.rightMT.(trialTypes{tt}) = zeros(nSubjects,nTasks,30);
    ERA.bilateralMT.(trialTypes{tt}) = zeros(nSubjects,nTasks,30); % for the average between left and right
    
    ERA.leftSPL.(trialTypes{tt}) = zeros(nSubjects,nTasks,30); % initialize - six runs, 30 volumes
    ERA.rightSPL.(trialTypes{tt}) = zeros(nSubjects,nTasks,30);
    ERA.bilateralSPL.(trialTypes{tt}) = zeros(nSubjects,nTasks,30); % for the average between left and right
    
    for ss = 1:nSubjects
        
        for rr = 1:nTasks
            
            auxIdx = squeeze(trialIndexes.(trialTypes{tt})(rr,:,:));
                     
            ERA.leftMT.(trialTypes{tt})(ss,rr,:) = mean(reshape(PSC.leftMT(ss,rr,auxIdx),2,30),1);
            ERA.rightMT.(trialTypes{tt})(ss,rr,:) = mean(reshape(PSC.rightMT(ss,rr,auxIdx),2,30),1);
            ERA.bilateralMT.(trialTypes{tt})(ss,rr,:) = (ERA.leftMT.(trialTypes{tt})(ss, rr ,:) + ERA.rightMT.(trialTypes{tt})(ss, rr ,:)) ./ 2;
            
            ERA.leftSPL.(trialTypes{tt})(ss,rr,:) = mean(reshape(PSC.leftSPL(ss,rr,auxIdx),2,30),1);
            ERA.rightSPL.(trialTypes{tt})(ss,rr,:) = mean(reshape(PSC.rightSPL(ss,rr,auxIdx),2,30),1);
            ERA.bilateralSPL.(trialTypes{tt})(ss,rr,:) = (ERA.leftSPL.(trialTypes{tt})(ss, rr ,:) + ERA.rightSPL.(trialTypes{tt})(ss, rr ,:)) ./ 2;

        end
        
    end
    
end

%% Figure 1 - Bilateral hMT+
fig1 = figure('Name','Group ERA','Position',[100 100 1100 900]);
xvector = -2:27;
clrMap = lines;
trialLabels = {'Coherent \rightarrow Coherent';'Incoherent \rightarrow Coherent';'Non-adapting \rightarrow Coherent';'Coherent \rightarrow Incoherent';'Incoherent \rightarrow Incoherent';'Non-adapting \rightarrow Incoherent'};
% ------------------------------------------------------------------------%
s1=subplot(2,1,1);
absMax = 0;
for tt = [1 2 3]

    toPlotMean = mean(mean(ERA.bilateralMT.(trialTypes{tt}),2),1,'omitnan');
    toPlotSem = std(mean(ERA.bilateralMT.(trialTypes{tt}),2),1,'omitnan') / sqrt(nSubjects);
    
    shadedErrorBar(xvector,...
        toPlotMean,...
        toPlotSem,...
        'lineprops',{'.-','Color',clrMap(tt,:),'LineWidth',1.5,'MarkerSize',10});
    
    if max(toPlotMean+toPlotSem) > absMax
        absMax = max(toPlotMean+toPlotSem);
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
text(20,  yy(2)-0.4,'MAE','HorizontalAlignment','center','FontSize',14)
text(23,  yy(2)-0.4,'Report','HorizontalAlignment','center','FontSize',14)

legend(trialLabels([1 2 3]),'fontSize',14, 'location','southeast')
s1.FontSize = 14;
% ------------------------------------------------------------------------%
% ------------------------------------------------------------------------%
s2=subplot(2,1,2);
absMax = 0;
auxclrmap = [2 1 3]; % to match the color of 'opposite' trials
for tt = [5 4 6]

    toPlotMean = mean(mean(ERA.bilateralMT.(trialTypes{tt}),2),1,'omitnan');
    toPlotSem = std(mean(ERA.bilateralMT.(trialTypes{tt}),2),1,'omitnan') / sqrt(nSubjects);
    
    shadedErrorBar(xvector,...
        toPlotMean,...
        toPlotSem,...
        'lineprops',{'.-','Color',clrMap(auxclrmap(tt-3),:),'LineWidth',1.5,'MarkerSize',10});
    
    if max(toPlotMean+toPlotSem) > absMax
        absMax = max(toPlotMean+toPlotSem);
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
text(20,  yy(2)-0.4,'MAE','HorizontalAlignment','center','FontSize',14)
text(23,  yy(2)-0.4,'Report','HorizontalAlignment','center','FontSize',14)

legend(trialLabels([1 2 3]),'fontSize',14, 'location','southeast')
s2.FontSize = 14;
% ------------------------------------------------------------------------%
% ------------------------------------------------------------------------%

%% Export fig1
print(fig1,fullfile(outputFolder,sprintf('Figures_fmriprep_1_psc.png')),'-dpng')
print(fig1,fullfile(outputFolder,sprintf('Figures_fmriprep_1_psc.svg')),'-dsvg')

%% Estimate test window and Statistical tests
%toPlotMean = zeros(8,6);
%toPlotSem = zeros(8,6);
toTestMean = struct();
delay_tc = 3;
x_auc = (15:22) + delay_tc; % 1 before + test block + 1 after

for tt = 1:6
    
    auxMean = squeeze(mean(ERA.bilateralMT.(trialTypes{tt}),2));
    
    auxMean = auxMean(~isnan(auxMean(:,1)),:); % remove NaN of one subject
    
    toTestMean.(trialTypes{tt}) = auxMean(:,x_auc) - auxMean(:,x_auc(1));
    
end

nComb = 3;

p12 = zeros(1,8);
h12 = zeros(1,8);
p34 = zeros(1,8);
h34 = zeros(1,8);

for jj = 1:8
    
    [~,p12aux] = ttest(toTestMean.Coh_aInCoh(:,jj),toTestMean.Coh_aCoh(:,jj));

    p12(jj) = p12aux * nComb;
    h12(jj) = p12(jj) <= 0.05;
    
    [~,p34aux] = ttest(toTestMean.InCoh_aCoh(:,jj),toTestMean.InCoh_aInCoh(:,jj));

    p34(jj) = p34aux * nComb;
    h34(jj) = p34(jj) <= 0.05;
    
end

%% Figure 2 - Bilateral hMT+
fig2 = figure('position',[150 150 1200 450]);
trialLabels = {'Coherent \rightarrow Coherent';'Incoherent \rightarrow Coherent';'Non-adapting \rightarrow Coherent';'Coherent \rightarrow Incoherent';'Incoherent \rightarrow Incoherent';'Non-adapting \rightarrow Incoherent'};
clrMap = lines;
auxclrmap = [2 1 3]; % to match the color of 'opposite' trials

xvector = 0:length(x_auc)-1;
xx = [-1 length(x_auc)];
yy = [-0.8 0.8];
yvector = sort([yy(1):0.2:yy(2)]);
lwidth = 1.5;

s1 = subplot(1,2,1);
line([0.5 0.5],yy,'LineStyle',':','Color','k')
% line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

for tt = [1 2 3]
    
    toPlotMean = mean(mean(ERA.bilateralMT.(trialTypes{tt}),2),1,'omitnan');
    toPlotSem = std(mean(ERA.bilateralMT.(trialTypes{tt}),2),1,'omitnan') / sqrt(nSubjects);
    
    toPlotMean = toPlotMean(x_auc) - toPlotMean(x_auc(1));
    toPlotSem = toPlotSem(x_auc);

    e(tt) = shadedErrorBar(xvector,toPlotMean,toPlotSem,...
        'lineprops',{'Color',clrMap(tt,:),'LineWidth',lwidth,'Marker','.','MarkerSize',20});

end

hold on

h12_x = find(h12);
plot(h12_x-1, 0.55*ones(1,length(h12_x)),'*','Color','k','LineWidth',1,'MarkerSize',8)

hold off

xlim(xx), ylim(yy)
xticks(xvector), xticklabels(x_auc-delay_tc)
yticks(yvector)
legend([e(1).mainLine e(2).mainLine e(3).mainLine],trialLabels(1:3),'FontSize',14,'Location','Southwest')
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal change (%)'},'FontSize',12)
s1.FontSize = 15;

%-------------------------------------------------------
s4 = subplot(1,2,2);
line([0.5 0.5],yy,'LineStyle',':','Color','k')
% line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

for tt = [5 4 6]

    toPlotMean = mean(mean(ERA.bilateralMT.(trialTypes{tt}),2),1,'omitnan');
    toPlotSem = std(mean(ERA.bilateralMT.(trialTypes{tt}),2),1,'omitnan') / sqrt(nSubjects);
    
    toPlotMean = toPlotMean(x_auc) - toPlotMean(x_auc(1));
    toPlotSem = toPlotSem(x_auc);

    e(tt) = shadedErrorBar(xvector,toPlotMean,toPlotSem,...
        'lineprops',{'Color',clrMap(auxclrmap(tt-3),:),'LineWidth',lwidth,'Marker','.','MarkerSize',20});

end

hold on

h34_x = find(h34);
plot(h34_x-1, 0.55*ones(1,length(h34_x)),'*','Color','k','LineWidth',1,'MarkerSize',8)

hold off

xlim(xx), ylim(yy)
xticks(xvector), xticklabels(x_auc-delay_tc)
yticks(yvector)
legend([e(5).mainLine e(4).mainLine e(6).mainLine],trialLabels([5 4 6]),'FontSize',14,'Location','southwest')
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal change (%)'},'FontSize',12)
s4.FontSize = 15;

%% Export fig2
print(fig2,fullfile(outputFolder,sprintf('Figures_fmriprep_2_psc.png')),'-dpng')
print(fig2,fullfile(outputFolder,sprintf('Figures_fmriprep_2_psc.svg')),'-dsvg')

%% SAVE
save('ERA-data.mat')
