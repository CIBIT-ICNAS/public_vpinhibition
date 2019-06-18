% ======================================================================= %
% == INHIBITION EXPERIMENT SERIES ======================================= %
% EXPERIMENT B
% BEHAVIORAL DATA ANALYSIS - SINGLE SUBJECT
% v0.1
%
% AUTHORS:
% VISUAL PERCEPTION GROUP at CIBIT-ICNAS
% 2018 - 2019
% ======================================================================= %
% ======================================================================= %

clear,clc

%% Initialize stuff

% Folders
drive_folder = fullfile('...','ICNAS_VisualPerception','Inhibition');

input_prtPath = fullfile(drive_folder,'prt'); % Protocol files
input_keyPath = fullfile(drive_folder,'expB-raw-keypress-output'); % raw keypress data (generated by the stimulus scripts)
output_keyPath = fullfile(drive_folder,'expB-analysis-keypress'); % output folder

% Runs to analyze
runs_D12 = {'RunBehav_D12_R1','RunBehav_D12_R2','RunBehav_D12_R3'};
nRuns_D12 = length(runs_D12);

% Ambiguous conditions
cond = {'Ambiguous_aCoh','Ambiguous_aInCoh','Ambiguous_aNA'};
condN_A = [5,6,7];
nCond_A = length(cond);

% Create result matrices
Results_D12 = zeros(nRuns_D12,nCond_A); % nRuns x n Ambiguous Conditions

%% Read keypress data
subject = 'S00';

fileDir = dir(fullfile(input_keyPath,[subject '_OUT_*.mat']));

load(fullfile(fileDir.folder,fileDir.name));

%% Iterate for the D12 runs
% 1) Load protocol data for each run
% 2) Read the keypresses for each ambiguous condition
% 3) Calculate the ratio of coherent responses

for rr = 1:nRuns_D12
    
    load(fullfile(input_prtPath,['Protocols_' runs_D12{rr} '.mat']));
    
    for cc = 1:nCond_A
        
        aux = Output.(runs_D12{rr}).Keys(framesCond == condN_A(cc),1);
        aux(aux == 0) = []; % remove zeros
        
        Results_D12(rr,cc) = sum(aux == Output.(runs_D12{rr}).KeyCodes(1)) / length(aux);
    
    end
     
end

%% Plot results D12
figure('Position',[150 150 600 450])

notBoxPlot(Results_D12,[1 2 3])
ylim([-0.2 1.2]), xlim([0 4])
line([0 4],[0 0],'LineStyle','--','color','k'), line([0 4],[1 1],'LineStyle','--','color','k')
ylabel('Ratio of Coherent responses')
xticklabels({'After Coherent','After Incoherent','After NonAdapt'})
title('D12')

%% Export for group analysis
save(fullfile(output_keyPath,[subject '_ExpB_KeyResults.mat']),'Results_D12')

disp('Done.')
