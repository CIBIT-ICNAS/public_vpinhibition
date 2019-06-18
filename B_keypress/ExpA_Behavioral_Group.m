% ======================================================================= %
% == INHIBITION EXPERIMENT SERIES ======================================= %
% EXPERIMENT A
% BEHAVIORAL DATA ANALYSIS - GROUP ANALYSIS
% v0.1
%
% AUTHORS:
% VISUAL PERCEPTION GROUP at CIBIT-ICNAS
% 2018 - 2019
% ======================================================================= %
% ======================================================================= %

clear,clc,close all

%% Initialize stuff

% Folders
drive_folder = fullfile('...','ICNAS_VisualPerception','Inhibition');

input_prtPath = fullfile(drive_folder,'prt'); % Protocol files
input_output_keyPath = fullfile(drive_folder,'expA-analysis-keypress'); % output folder

inputFolderDir = dir(fullfile(input_output_keyPath,'S*_KeyResults.mat'));
nSubjects = length(inputFolderDir);
nRuns = 3;

Results_D6_G = zeros(nSubjects*nRuns,3);
Results_D12_G = zeros(nSubjects*nRuns,3);

Results_D6_G_PerSub = zeros(nSubjects,3);
Results_D12_G_PerSub = zeros(nSubjects,3);

for ss = 1:nSubjects
    
    load(fullfile(inputFolderDir(ss).folder,inputFolderDir(ss).name));
    
    Results_D6_G(ss*nRuns-2 : ss*nRuns,:) = Results_D6;
    Results_D12_G(ss*nRuns-2 : ss*nRuns,:) = Results_D12;
    
    Results_D6_G_PerSub(ss,:) = mean(Results_D6,1);
    Results_D12_G_PerSub(ss,:) = mean(Results_D12,1);
    
end

%% Plot group results
% Points are layed over a 1.96 SEM (95% confidence interval) in red and a 1 SD in blue.
% Red line is the mean.
figure('units','normalized','position',[0 0 0.2 0.7])
movegui('center')

subplot(2,1,1)

line([0 4],[0 0],'LineStyle','--','color','k'), line([0 4],[1 1],'LineStyle','--','color','k'), hold on
H1 = notBoxPlot(Results_D6_G,[1 2 3]);
ylim([-0.2 1.2]), xlim([0 4])

ylabel('Ratio of Coherent responses')
yticks(0:0.2:1)
xticklabels({'    After\newlineCoherent','    After\newlineIncoherent','    After\newlineNonAdapt'})
title('Ambiguous block of 6 seconds')

set([H1.data],...
    'MarkerFaceColor',[1,1,1]*0.25,...
    'markerEdgeColor',[1,1,1]*0.25,...
    'MarkerSize',4)

set([H1.mu],'color','w')

J = lines(length(H1));
for ii=1:length(H1)
  set(H1(ii).sdPtch,'FaceColor',J(ii,:),...
                   'EdgeColor','none')

  set(H1(ii).semPtch,'FaceColor',J(ii,:)*0.3,...
                   'EdgeColor','none')
  
end

set(gca,'FontSize',16)

box on
% ------------------------------------------------------------------------%
subplot(2,1,2)

line([0 4],[0 0],'LineStyle','--','color','k'), line([0 4],[1 1],'LineStyle','--','color','k'), hold on
H2 = notBoxPlot(Results_D12_G,[1 2 3]);
ylim([-0.2 1.2]), xlim([0 4])
ylabel('Ratio of Coherent responses')
yticks(0:0.2:1)
xticklabels({'    After\newlineCoherent','    After\newlineIncoherent','    After\newlineNonAdapt'})
title('Ambiguous block of 12 seconds')

set([H2.data],...
    'MarkerFaceColor',[1,1,1]*0.25,...
    'markerEdgeColor',[1,1,1]*0.25,...
    'MarkerSize',4)

set([H2.mu],'color','w')

J = lines(length(H2));
for ii=1:length(H2)
  set(H2(ii).sdPtch,'FaceColor',J(ii,:),...
                   'EdgeColor','none')

  set(H2(ii).semPtch,'FaceColor',J(ii,:)*0.3,...
                   'EdgeColor','none')
  
end

set(gca,'FontSize',16)

box on
% ------------------------------------------------------------------------%

%% Export
save(fullfile(input_output_keyPath,...
    sprintf('Group_N%i_ExpA_KeyResults.mat',nSubjects)),...
    'Results_D12_G','Results_D6_G','Results_D12_G_PerSub','Results_D6_G_PerSub')

disp('Done.')
