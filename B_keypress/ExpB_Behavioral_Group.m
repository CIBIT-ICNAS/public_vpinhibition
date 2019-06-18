% ======================================================================= %
% == INHIBITION EXPERIMENT SERIES ======================================= %
% EXPERIMENT B
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
input_output_keyPath = fullfile(drive_folder,'expB-analysis-keypress'); % output folder

inputFolderDir = dir(fullfile(input_output_keyPath,'S*_KeyResults.mat'));
nSubjects = length(inputFolderDir);

nRuns = 3;

Results_D12_G = zeros(nSubjects*nRuns,3);

Results_D12_G_PerSub = zeros(nSubjects,3);

for ss = 1:nSubjects
    
    load(fullfile(inputFolderDir(ss).folder,inputFolderDir(ss).name));
    
    Results_D12_G(ss*nRuns-2 : ss*nRuns,:) = Results_D12;
    
    Results_D12_G_PerSub(ss,:) = mean(Results_D12,1);
    
end

%% Plot group results
% Points are layed over a 1.96 SEM (95% confidence interval) in red and a 1 SD in blue.
% Red line is the mean.
figure('units','normalized','position',[0 0 0.25 0.35])
movegui('center')

line([0 4],[0 0],'LineStyle','--','color','k'), line([0 4],[1 1],'LineStyle','--','color','k'), hold on
H = notBoxPlot(Results_D12_G,[1 2 3]);
ylim([-0.2 1.2]), xlim([0 4])
yticks(0:0.2:1)
ylabel('Ratio of Coherent responses')
xticklabels({'    After\newlineCoherent','    After\newlineIncoherent','    After\newlineNonAdapt'})
% title('D12 - notBoxPlot')

set([H.data],...
    'MarkerFaceColor',[1,1,1]*0.25,...
    'markerEdgeColor',[1,1,1]*0.25,...
    'MarkerSize',4)

set([H.mu],'color','w')

J = lines(length(H));
for ii=1:length(H)
  set(H(ii).sdPtch,'FaceColor',J(ii,:),...
                   'EdgeColor','none')

  set(H(ii).semPtch,'FaceColor',J(ii,:)*0.3,...
                   'EdgeColor','none')
  
end

set(gca,'FontSize',16)

box on

%% Statistics
ttest(Results_D12_G_PerSub(:,1),Results_D12_G_PerSub(:,2))

%% Export
save(fullfile(input_output_keyPath,...
    sprintf('GroupD12_N%i_ExpB_KeyResults.mat',nSubjects)),...
    'Results_D12_G','Results_D12_G_PerSub')

disp('Done.')
