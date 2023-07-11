% ======================================================================= %
% == INHIBITION EXPERIMENT SERIES ======================================= %
% EXPERIMENT C
% FMRI DATA ANALYSIS - SINGLE SUBJECT - SPM12
% v0.1
%
% AUTHORS:
% VISUAL PERCEPTION GROUP at CIBIT-ICNAS
% 2018 - 2019
% ======================================================================= %
% ======================================================================= %

clear ; clc, close all;
addpath(fullfile('..','tools','suptitle'))

%% Load Stuff and Set Stuff

% --- Configuration data
load('Configs_VP_INHIBITION_SPM.mat')

% --- Subject Name
subjectName = 'VPIS01';
subjectIndex = find(not(cellfun('isempty', strfind(datasetConfigs.subjects, subjectName))));

% --- Select Run Sequence
datasetConfigs.runs = datasetConfigs.runs{subjectIndex};

% --- I/O Folders
ioFolder = fullfile('..','data','expC-analysis-spm');
outputFolder = fullfile('..','data','expC-analysis-spm-images');
prtFilePath = fullfile('..','data','prt');

% --- Runs to process
selectedRuns = datasetConfigs.runs(2:end);
nTrialsPerRun = 2;
nTrials = length(selectedRuns) * nTrialsPerRun;

% --- Structures to keep data
TCP = struct(); % TimeCourse and Protocol data
ERA = struct(); % ERA data

% --- Haemodynamic Delay for the plots
delay_tc = 4; % in volumes
delay_plot = 0; %in volumes

%% Extract experimental Protocol info

for r = selectedRuns
    
    % --- Read Protocol
    prtFileName = ['RunMRI_D12_R' r{:}(end) '.prt'];   
    
    [ cond_names , TCP.(r{:}).intervalsPRT , TCP.(r{:}).intervals ] = ...
        readProtocol( prtFilePath , prtFileName , datasetConfigs.TR );
    
end

%% Load TC data
load(fullfile(ioFolder,'TCs-MT-5mm.mat'))

voiNames = {'leftMT','rightMT','bilateralMT'};

%% Create ERA

% --- Trial Types
trialTestTps = cond_names(5:10);

% --- extra points before
extraB = 3;

% Iterate on the ROIs
for v = 1:nRois
        
    % Iterate on the trial types
    for t = 1:length(trialTestTps)
        
        idx = 1;

        % --- Initialise matrices.
        nTrialVols = 30; % 3 (static) + 12 (adaptation) + 6 (test) + 3 (MAE) + 3 (report) + 3 (static) = 30 volumes
        ERA.(voiNames{v}).(trialTestTps{t}).psc.data = zeros(nTrials,nTrialVols);
        ERA.(voiNames{v}).(trialTestTps{t}).psc2static.data = zeros(nTrials,nTrialVols);
        % For the mean (nTrials+1), sem (nTrials+2), median (nTrials+3) and semedian (nTrials+4) of the above trials.
        ERA.(voiNames{v}).(trialTestTps{t}).psc.stats = zeros(4,nTrialVols);
        ERA.(voiNames{v}).(trialTestTps{t}).psc2stats.stats = zeros(4,nTrialVols);
        
        r_idx = 1;
        
        % Iterate on the runs
        for r = selectedRuns
            
            % --- Volumes of interest
            % I know there are two trials per run, so I spare a for loop :)
            int_aux = [];
            int_aux(1,:) = TCP.(r{:}).intervalsPRT.(trialTestTps{t})(1,1) - 12 - extraB : TCP.(r{:}).intervalsPRT.(trialTestTps{t})(1,2) + 3 + 3 + extraB ;
            int_aux(2,:) = TCP.(r{:}).intervalsPRT.(trialTestTps{t})(2,1) - 12 - extraB : TCP.(r{:}).intervalsPRT.(trialTestTps{t})(2,2) + 3 + 3 + extraB ;

            % --- Compensate Haemodynamic Delay
            int_aux = int_aux + delay_plot;
                        
            % --- Calculate BOLD PSC
            ERA.(voiNames{v}).(trialTestTps{t}).psc.data(idx+0,:) = (tcArray{subjectIndex,r_idx,v}(int_aux(1,:)) - mean(tcArray{subjectIndex,r_idx,v}(int_aux(1,:)))) / mean(tcArray{subjectIndex,r_idx,v}(int_aux(1,:))) * 100;        
            ERA.(voiNames{v}).(trialTestTps{t}).psc.data(idx+1,:) = (tcArray{subjectIndex,r_idx,v}(int_aux(2,:)) - mean(tcArray{subjectIndex,r_idx,v}(int_aux(2,:)))) / mean(tcArray{subjectIndex,r_idx,v}(int_aux(2,:))) * 100;

            % --- Calculate BOLD Signal Var to Baseline per run
            ERA.(voiNames{v}).(trialTestTps{t}).psc2static.data(idx+0,:) = pscArray{subjectIndex,r_idx,v}(int_aux(1,:)) * 100;        
            ERA.(voiNames{v}).(trialTestTps{t}).psc2static.data(idx+1,:) = pscArray{subjectIndex,r_idx,v}(int_aux(2,:)) * 100;
                       
            idx = idx + nTrialsPerRun;
            r_idx = r_idx + 1;
            
        end
        
        % --- Mean of Trials
        ERA.(voiNames{v}).(trialTestTps{t}).psc.stats(1,:) = mean(ERA.(voiNames{v}).(trialTestTps{t}).psc.data(1:idx-1,:));
        ERA.(voiNames{v}).(trialTestTps{t}).psc2static.stats(1,:) = mean(ERA.(voiNames{v}).(trialTestTps{t}).psc2static.data(1:idx-1,:));
        
        % --- Std Error of the Mean of Trials
        ERA.(voiNames{v}).(trialTestTps{t}).psc.stats(2,:) = std(ERA.(voiNames{v}).(trialTestTps{t}).psc.data(1:idx-1,:)) / sqrt(nTrials);
        ERA.(voiNames{v}).(trialTestTps{t}).psc2static.stats(2,:) = std(ERA.(voiNames{v}).(trialTestTps{t}).psc2static.data(1:idx-1,:)) / sqrt(nTrials);
        
        % --- Median of Trials
        ERA.(voiNames{v}).(trialTestTps{t}).psc.stats(3,:) = median(ERA.(voiNames{v}).(trialTestTps{t}).psc.data(1:idx-1,:));
        ERA.(voiNames{v}).(trialTestTps{t}).psc2static.stats(3,:) = median(ERA.(voiNames{v}).(trialTestTps{t}).psc2static.data(1:idx-1,:));
        
        % --- Std Error of the Median of Trials
        ERA.(voiNames{v}).(trialTestTps{t}).psc.stats(4,:) = 1.253*std(ERA.(voiNames{v}).(trialTestTps{t}).psc.data(1:idx-1,:)) / sqrt(nTrials);
        ERA.(voiNames{v}).(trialTestTps{t}).psc2static.stats(4,:) = 1.253*std(ERA.(voiNames{v}).(trialTestTps{t}).psc2static.data(1:idx-1,:)) / sqrt(nTrials);
        
    end
    
end

%% Plot Averages for Bilateral MT - MEAN - PSC to mean
fig1 = figure('Name','ERA','Position',[100 100 1300 1000]);
movegui('center')

xvector = -2:nTrialVols-3;
clrMap = lines;
trialLabels = {'Coh aCoh';'Coh aInCoh';'Coh aNA';'InCoh aCoh';'InCoh aInCoh';'InCoh aNA'};

% ------------------------------------------------------------------------%
subplot(2,1,1)
absMax = 0;
for tt = [1 2 3]
    
    errorbar(xvector,...
        ERA.bilateralMT.(trialTestTps{tt}).psc.stats(1,:),...
        ERA.bilateralMT.(trialTestTps{tt}).psc.stats(2,:),...
        '.-','LineWidth',1.5,'MarkerSize',15,'Color',clrMap(tt,:));
    
    if max(ERA.bilateralMT.(trialTestTps{tt}).psc.stats(1,:)+ERA.bilateralMT.(trialTestTps{tt}).psc.stats(2,:)) > absMax
        absMax = max(ERA.bilateralMT.(trialTestTps{tt}).psc.stats(1,:)+ERA.bilateralMT.(trialTestTps{tt}).psc.stats(2,:));
    end
    
    hold on
    
end
hold off
xx = [xvector(1)-1 xvector(end)+1];
yy = [-2 ceil(absMax+0.5)];
xlim(xx); ylim(yy);
xlabel('Time (volumes)')
ylabel({'BOLD signal change (%)','relative to the mean'})

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

legend(trialLabels([1 2 3]),'FontSize',13)
% ------------------------------------------------------------------------%
% ------------------------------------------------------------------------%
subplot(2,1,2)
absMax = 0;
for tt = [4 5 6]
    
    errorbar(xvector,...
        ERA.bilateralMT.(trialTestTps{tt}).psc.stats(1,:),...
        ERA.bilateralMT.(trialTestTps{tt}).psc.stats(2,:),...
        '.-','LineWidth',1.5,'MarkerSize',15,'Color',clrMap(tt,:));
    
    if max(ERA.bilateralMT.(trialTestTps{tt}).psc.stats(1,:)+ERA.bilateralMT.(trialTestTps{tt}).psc.stats(2,:)) > absMax
        absMax = max(ERA.bilateralMT.(trialTestTps{tt}).psc.stats(1,:)+ERA.bilateralMT.(trialTestTps{tt}).psc.stats(2,:));
    end
    
    hold on    
end

hold off
xx = [xvector(1)-1 xvector(end)+1];
yy = [-2 ceil(absMax+0.5)];
xlim(xx); ylim(yy);
xlabel('Time (volumes)')
ylabel({'BOLD signal change (%)','relative to the mean'})

xticks([1 6 12 13 18 19 21 22 24])

line(xx,zeros(length(xx),1),'LineStyle',':','Color','k')

line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([12.5 12.5],yy,'LineStyle','--','Color','k')
line([18.5 18.5],yy,'LineStyle','--','Color','k')
line([21.5 21.5],yy,'LineStyle','--','Color','k')
line([24.5 24.5],yy,'LineStyle','--','Color','k')

text(6,yy(2)-0.2,'Motion','HorizontalAlignment','center','FontSize',12)
text(15.5,yy(2)-0.2,'InCoherent','HorizontalAlignment','center','FontSize',12)
text(20,yy(2)-0.2,'MAE','HorizontalAlignment','center','FontSize',12)
text(23,yy(2)-0.2,'Report','HorizontalAlignment','center','FontSize',12)

legend(trialLabels([4 5 6]),'FontSize',13)
% ------------------------------------------------------------------------%

suptitle(sprintf('ERA of Bilateral hMT+ - Subject %s - Full trial',subjectName))

%% Save above figure
print(fig1,fullfile(outputFolder,[ subjectName '_BilateralMT_D' num2str(delay_plot) '_FullTrial_psc']),'-dpng')

%% Plot Averages for Bilateral MT - MEAN - PSC to mean
fig2 = figure('Name','ERA','Position',[100 100 1300 1000]);
movegui('center')

xvector = -2:nTrialVols-3;
clrMap = lines;
trialLabels = {'Coh aCoh';'Coh aInCoh';'Coh aNA';'InCoh aCoh';'InCoh aInCoh';'InCoh aNA'};

% ------------------------------------------------------------------------%
subplot(2,1,1)
absMax = 0;
for tt = [1 2 3]
    
    errorbar(xvector,...
        ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(1,:),...
        ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(2,:),...
        '.-','LineWidth',1.5,'MarkerSize',15,'Color',clrMap(tt,:));
    
    if max(ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(1,:)+ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(2,:)) > absMax
        absMax = max(ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(1,:)+ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(2,:));
    end
    
    hold on
    
end
hold off
xx = [xvector(1)-1 xvector(end)+1];
yy = [-1 ceil(absMax+0.5)];
xlim(xx); ylim(yy);
xlabel('Time (volumes)')
ylabel({'BOLD signal change (%)','relative to the static condition'})

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

legend(trialLabels([1 2 3]),'FontSize',13)
% ------------------------------------------------------------------------%
% ------------------------------------------------------------------------%
subplot(2,1,2)
absMax = 0;
for tt = [4 5 6]
    
    errorbar(xvector,...
        ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(1,:),...
        ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(2,:),...
        '.-','LineWidth',1.5,'MarkerSize',15,'Color',clrMap(tt,:));
    
    if max(ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(1,:)+ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(2,:)) > absMax
        absMax = max(ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(1,:)+ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(2,:));
    end
    
    hold on    
end

hold off
xx = [xvector(1)-1 xvector(end)+1];
yy = [-1 ceil(absMax+0.5)];
xlim(xx); ylim(yy);
xlabel('Time (volumes)')
ylabel({'BOLD signal change (%)','relative to the static condition'})

xticks([1 6 12 13 18 19 21 22 24])

line(xx,zeros(length(xx),1),'LineStyle',':','Color','k')

line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([12.5 12.5],yy,'LineStyle','--','Color','k')
line([18.5 18.5],yy,'LineStyle','--','Color','k')
line([21.5 21.5],yy,'LineStyle','--','Color','k')
line([24.5 24.5],yy,'LineStyle','--','Color','k')

text(6,yy(2)-0.2,'Motion','HorizontalAlignment','center','FontSize',12)
text(15.5,yy(2)-0.2,'InCoherent','HorizontalAlignment','center','FontSize',12)
text(20,yy(2)-0.2,'MAE','HorizontalAlignment','center','FontSize',12)
text(23,yy(2)-0.2,'Report','HorizontalAlignment','center','FontSize',12)

legend(trialLabels([4 5 6]),'FontSize',13)
% ------------------------------------------------------------------------%

suptitle(sprintf('ERA of Bilateral hMT+ - Subject %s - Full trial',subjectName))

%% Save above figure
print(fig2,fullfile(outputFolder,[ subjectName '_BilateralMT_D' num2str(delay_plot) '_FullTrial_psc2static']),'-dpng')

%% Calculate AuC of Test blocks
x_auc = (15:22) + delay_tc; % 1 before + test block + 1 after

Y_AUC = struct();
absMin13 = struct();
absMin46 = struct();

absMin13.psc = 99;
absMin13.psc2static = 99;
absMin46.psc = 99;
absMin46.psc2static = 99;

for tt = 1:3
    
    % --- Retrive TC values and SEM values
    Y_AUC.(trialTestTps{tt}).psc.psc = ERA.bilateralMT.(trialTestTps{tt}).psc.stats(1,x_auc);
    Y_AUC.(trialTestTps{tt}).psc.psc_sem = ERA.bilateralMT.(trialTestTps{tt}).psc.stats(2,x_auc);

    % --- Retrive TC values and SEM values (psc2static)
    Y_AUC.(trialTestTps{tt}).psc2static.psc = ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(1,x_auc);
    Y_AUC.(trialTestTps{tt}).psc2static.psc_sem = ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(2,x_auc);
    
    % --- Retrieve TC values normalising with the first value of each trial type
    Y_AUC.(trialTestTps{tt}).psc.psc_norm = Y_AUC.(trialTestTps{tt}).psc.psc - ...
        mean(ERA.bilateralMT.(trialTestTps{tt}).psc.stats(1,x_auc(1)));

    % --- Retrieve TC values normalising with the first value of each trial
    % type (psc2static)
    Y_AUC.(trialTestTps{tt}).psc2static.psc_norm = Y_AUC.(trialTestTps{tt}).psc2static.psc - ...
        mean(ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(1,x_auc(1)));
    
    % --- Find minimum of the three trials
    if min(Y_AUC.(trialTestTps{tt}).psc.psc_norm) < absMin13.psc
        absMin13.psc = min(Y_AUC.(trialTestTps{tt}).psc.psc_norm);
    end

    % --- Find minimum of the three trials (psc2static)
    if min(Y_AUC.(trialTestTps{tt}).psc2static.psc_norm) < absMin13.psc2static
        absMin13.psc2static = min(Y_AUC.(trialTestTps{tt}).psc2static.psc_norm);
    end
    
end

for tt = 4:6
    
    % --- Retrive TC values and SEM values
    Y_AUC.(trialTestTps{tt}).psc.psc = ERA.bilateralMT.(trialTestTps{tt}).psc.stats(1,x_auc);
    Y_AUC.(trialTestTps{tt}).psc.psc_sem = ERA.bilateralMT.(trialTestTps{tt}).psc.stats(2,x_auc);

    % --- Retrive TC values and SEM values (psc2static)
    Y_AUC.(trialTestTps{tt}).psc2static.psc = ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(1,x_auc);
    Y_AUC.(trialTestTps{tt}).psc2static.psc_sem = ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(2,x_auc);
    
    % --- Retrieve TC values normalising with the average of the three last
    % points before test of each trial type
    Y_AUC.(trialTestTps{tt}).psc.psc_norm = Y_AUC.(trialTestTps{tt}).psc.psc - ...
        mean(ERA.bilateralMT.(trialTestTps{tt}).psc.stats(1,x_auc(1)));

    % --- Retrieve TC values normalising with the average of the three last
    % points before test of each trial type (psc2static)
    Y_AUC.(trialTestTps{tt}).psc2static.psc_norm = Y_AUC.(trialTestTps{tt}).psc2static.psc - ...
        mean(ERA.bilateralMT.(trialTestTps{tt}).psc2static.stats(1,x_auc(1)));
    
    % --- Find minimum of the three trials
    if min(Y_AUC.(trialTestTps{tt}).psc.psc_norm) < absMin46.psc
        absMin46.psc = min(Y_AUC.(trialTestTps{tt}).psc.psc_norm);
    end

    % --- Find minimum of the three trials (psc2static)
    if min(Y_AUC.(trialTestTps{tt}).psc2static.psc_norm) < absMin46.psc2static
        absMin46.psc2static = min(Y_AUC.(trialTestTps{tt}).psc2static.psc_norm);
    end    
    
end

%% Remove minimum value of all trials
% This will make the minimum point of the trials = 0 --> push down the
% curves to avoid unnecessary area under the curves :)

for tt = 1:3
    
    Y_AUC.(trialTestTps{tt}).psc.psc_norm0 = Y_AUC.(trialTestTps{tt}).psc.psc_norm - absMin13.psc;  
    Y_AUC.(trialTestTps{tt}).psc2static.psc_norm0 = Y_AUC.(trialTestTps{tt}).psc2static.psc_norm - absMin13.psc2static;  
    
end

for tt = 4:6
    
    Y_AUC.(trialTestTps{tt}).psc.psc_norm0 = Y_AUC.(trialTestTps{tt}).psc.psc_norm - absMin46.psc;  
    Y_AUC.(trialTestTps{tt}).psc2static.psc_norm0 = Y_AUC.(trialTestTps{tt}).psc2static.psc_norm - absMin46.psc2static;  
    
end

%% Prepare The ultimate plot
AUC = struct();
trialLabels = {'Coh aCoh';'Coh aInCoh';'Coh aNA';'InCoh aCoh';'InCoh aInCoh';'InCoh aNA'};
clrMap = lines;
xvector = 0:length(x_auc)-1;
xx = [-1 length(x_auc)];
yy = [-1 1];
lwidth = 1.5;

%% Plot The ultimate plot (psc)
fig_ultimate1 = figure('Name','The Ultimate Figure','units','normalized','outerposition',[0 0 0.5 0.5]);
movegui('center')

% ------------------------------------------------------------------------%
% -- ERA Coherent --------------------------------------------------------%
s1 = subplot(1,2,1);
line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = errorbar(xvector,Y_AUC.Coh_aCoh.psc.psc_norm,Y_AUC.Coh_aCoh.psc.psc_sem,'Color',clrMap(1,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e2 = errorbar(xvector,Y_AUC.Coh_aInCoh.psc.psc_norm,Y_AUC.Coh_aInCoh.psc.psc_sem,'Color',clrMap(2,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e3 = errorbar(xvector,Y_AUC.Coh_aNA.psc.psc_norm,Y_AUC.Coh_aNA.psc.psc_sem,'Color',clrMap(3,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);

hold off

xlim(xx), ylim(yy)
xticks(xvector)
legend([e1 e2 e3],trialLabels(1:3),'FontSize',13,'Location','Southwest')
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal change (%)','relative to time 0'},'FontSize',12)
s1.FontSize = 12;
% ------------------------------------------------------------------------%
% -- ERA InCoherent ------------------------------------------------------%
s4 = subplot(1,2,2);
line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = errorbar(xvector,Y_AUC.InCoh_aCoh.psc.psc_norm,Y_AUC.InCoh_aCoh.psc.psc_sem,'Color',clrMap(4,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e2 = errorbar(xvector,Y_AUC.InCoh_aInCoh.psc.psc_norm,Y_AUC.InCoh_aInCoh.psc.psc_sem,'Color',clrMap(5,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e3 = errorbar(xvector,Y_AUC.InCoh_aNA.psc.psc_norm,Y_AUC.InCoh_aNA.psc.psc_sem,'Color',clrMap(6,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);

hold off

xlim(xx), ylim(yy)
xticks(xvector)
legend([e1 e2 e3],trialLabels(4:6),'FontSize',13,'Location','Southwest')
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal change (%)','relative to time 0'},'FontSize',12)
s4.FontSize = 12;
% ------------------------------------------------------------------------%
suptitle(sprintf('ERA of Bilateral hMT+ - Subject %s - Second motion block',subjectName))

%% Export The ultimate Figure
print(fig_ultimate1,fullfile(outputFolder,[ subjectName '_BilateralMT_D' num2str(delay_tc) '_2MotionBlock_psc']),'-dpng')

%% Plot The ultimate plot (psc)
fig_ultimate2 = figure('Name','The Ultimate Figure','units','normalized','outerposition',[0 0 0.5 0.5]);
movegui('center')

% ------------------------------------------------------------------------%
% -- ERA Coherent --------------------------------------------------------%
s1 = subplot(1,2,1);
line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = errorbar(xvector,Y_AUC.Coh_aCoh.psc2static.psc_norm,Y_AUC.Coh_aCoh.psc2static.psc_sem,'Color',clrMap(1,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e2 = errorbar(xvector,Y_AUC.Coh_aInCoh.psc2static.psc_norm,Y_AUC.Coh_aInCoh.psc2static.psc_sem,'Color',clrMap(2,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e3 = errorbar(xvector,Y_AUC.Coh_aNA.psc2static.psc_norm,Y_AUC.Coh_aNA.psc2static.psc_sem,'Color',clrMap(3,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);

hold off

xlim(xx), ylim(yy)
xticks(xvector)
legend([e1 e2 e3],trialLabels(1:3),'FontSize',13,'Location','Southwest')
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal change (%)','relative to time 0'},'FontSize',12)
s1.FontSize = 12;
% ------------------------------------------------------------------------%
% -- ERA InCoherent ------------------------------------------------------%
s4 = subplot(1,2,2);
line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = errorbar(xvector,Y_AUC.InCoh_aCoh.psc2static.psc_norm,Y_AUC.InCoh_aCoh.psc2static.psc_sem,'Color',clrMap(4,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e2 = errorbar(xvector,Y_AUC.InCoh_aInCoh.psc2static.psc_norm,Y_AUC.InCoh_aInCoh.psc2static.psc_sem,'Color',clrMap(5,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e3 = errorbar(xvector,Y_AUC.InCoh_aNA.psc2static.psc_norm,Y_AUC.InCoh_aNA.psc2static.psc_sem,'Color',clrMap(6,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);

hold off

xlim(xx), ylim(yy)
xticks(xvector)
legend([e1 e2 e3],trialLabels(4:6),'FontSize',13,'Location','Southwest')
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal change (%)','relative to time 0'},'FontSize',12)
s4.FontSize = 12;
% ------------------------------------------------------------------------%
suptitle(sprintf('ERA of Bilateral hMT+ - Subject %s - Second motion block',subjectName))

%% Export The ultimate Figure
print(fig_ultimate2,fullfile(outputFolder,[ subjectName '_BilateralMT_D' num2str(delay_tc) '_2MotionBlock_psc2static']),'-dpng')

%% Save dataset
save(fullfile(ioFolder,[subjectName '_bold.mat']),'ERA','Y_AUC','trialLabels','x_auc','trialLabels','clrMap','trialTestTps','nTrialVols','subjectIndex','subjectName','delay_tc','delay_plot')
    